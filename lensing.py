"""
Differential lensing analysis for the metering-metric hypothesis.

The working observable is a sign-flip morphology around galaxies: a positive
residual-convergence core surrounded by a negative residual ring.

This script can:
- inventory public HST/JWST and Chandra coverage for the target clusters
- score sign-flip morphology from user-supplied residual maps and galaxy tables
- generate a toy dataset for pipeline checks
- run first-pass archive analyses on Abell 370 and RXC J2248 / Abell S1063

Generic `score` inputs:
- galaxies CSV: ra, dec, metallicity, sfr, mass
- residuals CSV: ra, dec, residual
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
from scipy.stats import chi2, spearmanr


MAST_API_URL = "https://mast.stsci.edu/api/v0/invoke"
MAST_ARCHIVE_URL = "https://archive.stsci.edu"
HFF_ARCHIVE_URL = "https://archive.stsci.edu/prepds/frontier/"
CHANDRA_ARCHIVE_URL = "https://cxc.harvard.edu/cda/"
CHANDRA_CSC_SCS_URL = "https://cda.cfa.harvard.edu/csc2scs/coneSearch"

DEFAULT_TARGETS = [
    "Abell 2744",
    "MACS J0416.1-2403",
    "MACS J0717.5+3745",
    "MACS J1149.5+2223",
    "Abell S1063",
    "Abell 370",
]

OFFICIAL_PRODUCTS = {
    "abell2744": {
        "glass_multiband_catalog": "https://archive.stsci.edu/hlsps/glass-jwst/hlsp_glass-jwst_jwst-hst_multi_abell2744_multiband_v1.0_cat.fits",
        "hff_cats_v4.1_kappa": "https://archive.stsci.edu/pub/hlsp/frontier/abell2744/models/cats/v4.1/hlsp_frontier_model_abell2744_cats_v4.1_kappa.fits",
    },
    "abell370": {
        "buffalo_proxy_catalog": "https://archive.stsci.edu/hlsps/buffalo/abell370/catalogs/pagul-v2.0/hlsp_buffalo_hst_ir-weighted_abell370_multi_v2.0_catalog.fits",
        "hff_properties_catalog": "https://archive.stsci.edu/hlsps/hffcatalogs/hlsp_hffcatalogs_multi_multi_abell370_multi_v1_prop-cat.fits",
        "hff_cats_v4_kappa": "https://archive.stsci.edu/pub/hlsp/frontier/abell370/models/cats/v4/hlsp_frontier_model_abell370_cats_v4_kappa.fits",
    },
    "rxcj2248": {
        "buffalo_proxy_catalog": "https://archive.stsci.edu/hlsps/buffalo/abells1063/catalogs/pagul-v2.0/hlsp_buffalo_hst_ir-weighted_abells1063_multi_v2.0_catalog.fits",
        "hff_properties_catalog": "https://archive.stsci.edu/hlsps/hffcatalogs/hlsp_hffcatalogs_multi_multi_rxc-j2248_multi_v1_prop-cat.fits",
        "hff_cats_v4.1_kappa": "https://archive.stsci.edu/pub/hlsp/frontier/abells1063/models/cats/v4.1/hlsp_frontier_model_abells1063_cats_v4.1_kappa.fits",
    },
}

HFF_CLUSTER_CONFIG = {
    "abell370": {
        "archive_cluster_name": "abell370",
        "proxy_key": "buffalo_proxy_catalog",
        "catalog_key": "hff_properties_catalog",
        "kappa_key": "hff_cats_v4_kappa",
        "proxy_filename": "abell370_buffalo_proxy.fits",
        "catalog_filename": "abell370_properties.fits",
        "kappa_filename": "abell370_kappa.fits",
        "label": "Abell 370",
        "default_z_window": (0.35, 0.45),
    },
    "rxcj2248": {
        "archive_cluster_name": "abells1063",
        "proxy_key": "buffalo_proxy_catalog",
        "catalog_key": "hff_properties_catalog",
        "kappa_key": "hff_cats_v4.1_kappa",
        "proxy_filename": "abells1063_buffalo_proxy.fits",
        "catalog_filename": "rxcj2248_properties.fits",
        "kappa_filename": "abells1063_kappa.fits",
        "label": "RXC J2248 / Abell S1063",
        "default_z_window": (0.35, 0.45),
    },
}

HFF_MODEL_VERSIONS = {
    "abell370": {
        "bradac": "v4.1",
        "cats": "v4",
        "diego": "v4.1",
        "glafic": "v4",
        "keeton": "v4",
        "merten": "v1",
        "sharon": "v4",
        "williams": "v4.1",
        "zitrin-ltm-gauss": "v1",
        "zitrin-ltm": "v1",
        "zitrin-nfw": "v1",
    },
    "rxcj2248": {
        "bradac": "v1",
        "cats": "v4.1",
        "diego": "v4.1",
        "glafic": "v4",
        "keeton": "v4",
        "merten": "v1",
        "sharon": "v4",
        "williams": "v4.1",
        "zitrin-ltm-gauss": "v1",
        "zitrin-ltm": "v1",
        "zitrin-nfw": "v1",
    },
}

BUFFALO_PROXY_BANDS = [
    ("FLUX_F275W", "FLUXERR_F275W", 0.275),
    ("FLUX_F336W", "FLUXERR_F336W", 0.336),
    ("FLUX_F435W", "FLUXERR_F435W", 0.435),
    ("FLUX_F475W", "FLUXERR_F475W", 0.475),
    ("FLUX_F606W", "FLUXERR_F606W", 0.606),
    ("FLUX_F625W", "FLUXERR_F625W", 0.625),
    ("FLUX_F814W", "FLUXERR_F814W", 0.814),
    ("FLUX_F105W", "FLUXERR_F105W", 1.05),
    ("FLUX_F110W", "FLUXERR_F110W", 1.10),
    ("FLUX_F125W", "FLUXERR_F125W", 1.25),
    ("FLUX_F140W", "FLUXERR_F140W", 1.40),
    ("FLUX_F160W", "FLUXERR_F160W", 1.54),
]
BUFFALO_PROXY_MATCH_MAX_ARCSEC = 0.3
BUFFALO_PROXY_MIN_BANDS = 6
BUFFALO_PROXY_MIN_SNR = 2.0

CLI_EPILOG = """Examples:
  python lensing.py manifest --out lensing_manifest.json
  python lensing.py simulate --out-dir toy_lensing
  python lensing.py score --galaxies galaxies.csv --residuals residuals.csv --out lensing_results.json
  python lensing.py abell370-firstpass --out-dir abell370_firstpass --residual-mode radial_median
  python lensing.py rxcj2248-firstpass --out-dir rxcj2248_firstpass --residual-mode radial_median
  python lensing.py combined-firstpass --out-dir combined_firstpass --residual-mode radial_median
  python lensing.py full-experiment --out-dir full_experiment
  python lensing.py robustness-sweep --out-dir robustness_sweep
"""


@dataclass
class Target:
    name: str
    ra: float
    dec: float
    radius_deg: float = 0.05


def mast_request(service: str, params: dict, page: int = 1, pagesize: int = 5000) -> dict:
    payload = {
        "service": service,
        "params": params,
        "format": "json",
        "page": page,
        "pagesize": pagesize,
    }
    response = requests.post(
        MAST_API_URL,
        data={"request": json.dumps(payload)},
        timeout=90,
    )
    response.raise_for_status()
    return response.json()


def resolve_target(name: str, radius_deg: float) -> Target:
    payload = mast_request(
        "Mast.Name.Lookup",
        {"input": name, "format": "json"},
        page=1,
        pagesize=20,
    )
    matches = payload.get("resolvedCoordinate", [])
    if not matches:
        raise ValueError(f"Could not resolve target name: {name}")
    best = matches[0]
    return Target(
        name=name,
        ra=float(best["ra"]),
        dec=float(best["decl"]),
        radius_deg=radius_deg,
    )


def query_public_imaging(target: Target, collections: Iterable[str] = ("HST", "JWST")) -> list[dict]:
    payload = mast_request(
        "Mast.Caom.Filtered.Position",
        {
            "columns": "*",
            "position": f"{target.ra}, {target.dec}, {target.radius_deg}",
            "filters": [
                {"paramName": "obs_collection", "values": list(collections)},
                {"paramName": "intentType", "values": ["science"]},
                {"paramName": "dataRights", "values": ["PUBLIC"]},
                {"paramName": "dataproduct_type", "values": ["image"]},
            ],
        },
        page=1,
        pagesize=5000,
    )
    return payload.get("data", [])


def summarize_records(records: list[dict]) -> dict:
    collections = Counter(record.get("obs_collection", "UNKNOWN") for record in records)
    instruments = Counter(record.get("instrument_name", "UNKNOWN") for record in records)
    proposals = Counter(str(record.get("proposal_id", "UNKNOWN")) for record in records)
    return {
        "record_count": len(records),
        "collections": dict(collections),
        "top_instruments": instruments.most_common(8),
        "top_proposals": proposals.most_common(8),
        "sample_obs_ids": [record.get("obs_id") for record in records[:10]],
    }


def query_chandra_sources(target: Target) -> dict:
    response = requests.get(
        CHANDRA_CSC_SCS_URL,
        params={
            "RA": target.ra,
            "DEC": target.dec,
            "SR": target.radius_deg,
            "VERB": 1,
        },
        timeout=90,
    )
    response.raise_for_status()

    root = ET.fromstring(response.text)
    namespace = {"v": "http://www.ivoa.net/xml/VOTable/v1.2"}
    fields = [field.attrib.get("name", "") for field in root.findall(".//v:FIELD", namespace)]
    rows = root.findall(".//v:TABLEDATA/v:TR", namespace)
    names = []
    if rows and fields:
        try:
            name_index = fields.index("name")
        except ValueError:
            name_index = 0
        for row in rows[:10]:
            cells = row.findall("v:TD", namespace)
            if name_index < len(cells):
                names.append(cells[name_index].text)

    return {
        "source_count": len(rows),
        "sample_source_names": names,
        "service_url": CHANDRA_CSC_SCS_URL,
    }


def build_manifest(target_names: Iterable[str], radius_deg: float) -> dict:
    targets = []
    for name in target_names:
        target = resolve_target(name, radius_deg)
        records = query_public_imaging(target)
        chandra = query_chandra_sources(target)
        targets.append(
            {
                **asdict(target),
                **summarize_records(records),
                "chandra": chandra,
                "official_links": {
                    "mast_archive": MAST_ARCHIVE_URL,
                    "hff_archive": HFF_ARCHIVE_URL,
                    "chandra_archive": CHANDRA_ARCHIVE_URL,
                },
            }
        )

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "selection": {
            "targets": list(target_names),
            "radius_deg": radius_deg,
            "collections": ["HST", "JWST"],
            "filters": ["science", "PUBLIC", "image"],
        },
        "targets": targets,
    }


def angular_distance_arcsec(ra1: float, dec1: float, ra2: np.ndarray, dec2: np.ndarray) -> np.ndarray:
    ra1_rad = math.radians(ra1)
    dec1_rad = math.radians(dec1)
    ra2_rad = np.radians(ra2)
    dec2_rad = np.radians(dec2)

    cos_sep = (
        math.sin(dec1_rad) * np.sin(dec2_rad)
        + math.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra2_rad - ra1_rad)
    )
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    return np.degrees(np.arccos(cos_sep)) * 3600.0


def angular_distance_arcsec_scalar(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    return float(angular_distance_arcsec(ra1, dec1, np.asarray([ra2]), np.asarray([dec2]))[0])


def read_csv_rows(path: Path) -> list[dict]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def download_file(url: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        return path
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return path


def frontier_model_kappa_url(cluster_key: str, model_family: str) -> str:
    config = HFF_CLUSTER_CONFIG[cluster_key]
    archive_cluster_name = config["archive_cluster_name"]
    version = HFF_MODEL_VERSIONS[cluster_key][model_family]
    return (
        f"{MAST_ARCHIVE_URL}/pub/hlsp/frontier/{archive_cluster_name}/models/{model_family}/{version}/"
        f"hlsp_frontier_model_{archive_cluster_name}_{model_family}_{version}_kappa.fits"
    )


def finite_positive_float(value: object) -> float | None:
    numeric = float(value)
    if np.isfinite(numeric) and numeric > 0.0:
        return numeric
    return None


def flux_microjy_to_ab_magnitude(flux_microjy: float) -> float:
    return 23.9 - 2.5 * math.log10(flux_microjy)


def buffalo_continuum_proxy(buffalo_row: object, redshift: float) -> tuple[float | None, int]:
    safe_redshift = redshift if np.isfinite(redshift) and redshift >= 0.0 else 0.0
    xs = []
    ys = []
    weights = []
    for flux_key, fluxerr_key, wavelength_microns in BUFFALO_PROXY_BANDS:
        flux = finite_positive_float(buffalo_row[flux_key])
        fluxerr = finite_positive_float(buffalo_row[fluxerr_key])
        if flux is None or fluxerr is None:
            continue
        snr = flux / fluxerr
        if snr < BUFFALO_PROXY_MIN_SNR:
            continue
        rest_wavelength_microns = wavelength_microns / (1.0 + safe_redshift)
        magnitude = flux_microjy_to_ab_magnitude(flux)
        sigma_mag = 2.5 / math.log(10.0) * (fluxerr / flux)
        xs.append(math.log10(rest_wavelength_microns))
        ys.append(magnitude)
        weights.append(1.0 / max(sigma_mag * sigma_mag, 1e-6))
    if len(xs) < BUFFALO_PROXY_MIN_BANDS:
        return None, len(xs)
    slope = float(
        np.polyfit(
            np.asarray(xs, dtype=float),
            np.asarray(ys, dtype=float),
            1,
            w=np.sqrt(np.asarray(weights, dtype=float)),
        )[0]
    )
    return slope, len(xs)


def load_photometric_proxy_lookup(catalog_path: Path, proxy_catalog_path: Path) -> dict[int, dict]:
    lookup = {}
    with fits.open(catalog_path) as hdul:
        hff_rows = hdul[1].data
        hff_coords = SkyCoord(
            ra=np.asarray(hff_rows["alpha_j2000"], dtype=float) * u.deg,
            dec=np.asarray(hff_rows["delta_j2000"], dtype=float) * u.deg,
        )
    with fits.open(proxy_catalog_path) as hdul:
        proxy_rows = hdul[1].data
        proxy_coords = SkyCoord(
            ra=np.asarray(proxy_rows["ALPHA_J2000_STACK"], dtype=float) * u.deg,
            dec=np.asarray(proxy_rows["DELTA_J2000_STACK"], dtype=float) * u.deg,
        )

    nearest_proxy_idx, nearest_sep, _ = hff_coords.match_to_catalog_sky(proxy_coords)
    reverse_hff_idx, _, _ = proxy_coords.match_to_catalog_sky(hff_coords)

    for row_index, record in enumerate(hff_rows):
        proxy_index = int(nearest_proxy_idx[row_index])
        match_sep_arcsec = float(nearest_sep[row_index].arcsec)
        if reverse_hff_idx[proxy_index] != row_index or match_sep_arcsec > BUFFALO_PROXY_MATCH_MAX_ARCSEC:
            continue
        proxy_row = proxy_rows[proxy_index]
        proxy, band_count = buffalo_continuum_proxy(proxy_row, redshift=float(record["ZBEST"]))
        lookup[int(record["IDcat"])] = {
            "photometric_proxy": proxy,
            "photometric_proxy_band_count": band_count,
            "photometric_proxy_match_sep_arcsec": match_sep_arcsec,
            "photometric_proxy_catalog_id": int(proxy_row["ID"]),
            "photometric_proxy_ebv": (
                float(proxy_row["E_BV"])
                if np.isfinite(float(proxy_row["E_BV"]))
                else float("nan")
            ),
            "photometric_proxy_nb_used": int(proxy_row["NB_USED"]),
            "photometric_proxy_catalog_name": "BUFFALO v2.0",
        }
    return lookup


def score_sign_flip(
    galaxies_path: Path,
    residuals_path: Path,
    out_path: Path,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
) -> dict:
    galaxies = read_csv_rows(galaxies_path)
    residuals = read_csv_rows(residuals_path)
    residual_ra = np.array([float(row["ra"]) for row in residuals], dtype=float)
    residual_dec = np.array([float(row["dec"]) for row in residuals], dtype=float)
    residual_value = np.array([float(row["residual"]) for row in residuals], dtype=float)

    scored_rows = []
    design = []
    response = []

    for galaxy in galaxies:
        ra = float(galaxy["ra"])
        dec = float(galaxy["dec"])
        metallicity = float(galaxy["metallicity"])
        sfr = float(galaxy["sfr"])
        mass = float(galaxy["mass"])

        separations = angular_distance_arcsec(ra, dec, residual_ra, residual_dec)
        inner_mask = separations <= inner_radius_arcsec
        ring_mask = (separations >= ring_inner_arcsec) & (separations <= ring_outer_arcsec)

        if not inner_mask.any() or not ring_mask.any():
            continue

        inner_mean = float(np.mean(residual_value[inner_mask]))
        ring_mean = float(np.mean(residual_value[ring_mask]))
        sign_flip_score = max(inner_mean, 0.0) + max(-ring_mean, 0.0)
        contrast = inner_mean - ring_mean
        sign_flip_pass = inner_mean > 0.0 and ring_mean < 0.0

        scored = {
            **galaxy,
            "inner_mean": inner_mean,
            "ring_mean": ring_mean,
            "sign_flip_score": sign_flip_score,
            "contrast": contrast,
            "sign_flip_pass": int(sign_flip_pass),
        }
        scored_rows.append(scored)
        design.append([metallicity, math.log10(max(sfr, 1e-6)), math.log10(max(mass, 1e-6))])
        response.append(sign_flip_score)

    if not scored_rows:
        raise ValueError("No galaxies retained after aperture scoring.")

    X = np.asarray(design, dtype=float)
    y = np.asarray(response, dtype=float)
    means = X.mean(axis=0)
    stds = X.std(axis=0)
    stds[stds == 0.0] = 1.0
    Xz = (X - means) / stds
    Xreg = np.column_stack([np.ones(len(Xz)), Xz])
    beta, *_ = np.linalg.lstsq(Xreg, y, rcond=None)
    yhat = Xreg @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "inner_radius_arcsec": inner_radius_arcsec,
            "ring_inner_arcsec": ring_inner_arcsec,
            "ring_outer_arcsec": ring_outer_arcsec,
        },
        "n_scored": len(scored_rows),
        "regression": {
            "response": "sign_flip_score",
            "predictors": ["metallicity", "log10_sfr", "log10_mass"],
            "coefficients": {
                "intercept": float(beta[0]),
                "metallicity_z": float(beta[1]),
                "log10_sfr_z": float(beta[2]),
                "log10_mass_z": float(beta[3]),
            },
            "predictor_means": {
                "metallicity": float(means[0]),
                "log10_sfr": float(means[1]),
                "log10_mass": float(means[2]),
            },
            "predictor_stds": {
                "metallicity": float(stds[0]),
                "log10_sfr": float(stds[1]),
                "log10_mass": float(stds[2]),
            },
            "r_squared": r_squared,
        },
    }

    write_csv_rows(
        out_path.with_suffix(".csv"),
        scored_rows,
        list(scored_rows[0].keys()),
    )
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def radial_kernel(separation_arcsec: np.ndarray, core_arcsec: float, boundary_arcsec: float) -> np.ndarray:
    core = np.exp(-0.5 * (separation_arcsec / core_arcsec) ** 2)
    ring = np.exp(-0.5 * ((separation_arcsec - boundary_arcsec) / core_arcsec) ** 2)
    return core - 0.7 * ring


def build_radius_map_arcsec(
    shape: tuple[int, int],
    center_xpix: float,
    center_ypix: float,
    pixel_scale_arcsec: float,
) -> np.ndarray:
    yy, xx = np.indices(shape, dtype=float)
    return np.sqrt((xx - center_xpix) ** 2 + (yy - center_ypix) ** 2) * pixel_scale_arcsec


def build_residual_map(
    kappa: np.ndarray,
    footprint: np.ndarray,
    center_xpix: float,
    center_ypix: float,
    pixel_scale_arcsec: float,
    residual_mode: str,
    smooth_sigma_px: float,
    radial_bin_arcsec: float,
) -> tuple[np.ndarray, dict]:
    if residual_mode == "gaussian_highpass":
        smooth = gaussian_filter(kappa, sigma=smooth_sigma_px)
        residual = kappa - smooth
        return residual, {
            "residual_mode": residual_mode,
            "smooth_sigma_px": smooth_sigma_px,
        }

    if residual_mode == "radial_median":
        radius_map_arcsec = build_radius_map_arcsec(
            shape=kappa.shape,
            center_xpix=center_xpix,
            center_ypix=center_ypix,
            pixel_scale_arcsec=pixel_scale_arcsec,
        )
        finite_mask = footprint & np.isfinite(kappa) & np.isfinite(radius_map_arcsec)
        if not np.any(finite_mask):
            raise ValueError("No finite pixels available for radial-median subtraction.")

        max_radius_arcsec = float(np.max(radius_map_arcsec[finite_mask]))
        bin_edges = np.arange(0.0, max_radius_arcsec + radial_bin_arcsec, radial_bin_arcsec)
        if bin_edges.size < 2:
            bin_edges = np.array([0.0, max_radius_arcsec + radial_bin_arcsec], dtype=float)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        radial_medians = np.full(bin_centers.shape, np.nan, dtype=float)

        for index in range(bin_centers.size):
            in_bin = finite_mask & (radius_map_arcsec >= bin_edges[index]) & (radius_map_arcsec < bin_edges[index + 1])
            if np.any(in_bin):
                radial_medians[index] = float(np.median(kappa[in_bin]))

        valid = np.isfinite(radial_medians)
        if np.count_nonzero(valid) < 2:
            raise ValueError("Radial-median subtraction did not find enough populated bins.")

        baseline = np.interp(
            radius_map_arcsec,
            bin_centers[valid],
            radial_medians[valid],
            left=float(radial_medians[valid][0]),
            right=float(radial_medians[valid][-1]),
        )
        residual = kappa - baseline
        residual[~footprint] = np.nan
        return residual, {
            "residual_mode": residual_mode,
            "radial_bin_arcsec": radial_bin_arcsec,
            "radial_bin_count": int(bin_centers.size),
            "populated_radial_bin_count": int(np.count_nonzero(valid)),
        }

    raise ValueError(f"Unsupported residual mode: {residual_mode}")


def fit_linear_model(design: np.ndarray, response: np.ndarray, predictor_names: list[str]) -> dict:
    means = design.mean(axis=0)
    stds = design.std(axis=0)
    stds[stds == 0.0] = 1.0
    Xz = (design - means) / stds
    Xreg = np.column_stack([np.ones(len(Xz)), Xz])
    beta, *_ = np.linalg.lstsq(Xreg, response, rcond=None)
    fitted = Xreg @ beta
    ss_res = float(np.sum((response - fitted) ** 2))
    ss_tot = float(np.sum((response - np.mean(response)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    coefficients = {"intercept": float(beta[0])}
    predictor_means = {}
    predictor_stds = {}
    for index, name in enumerate(predictor_names):
        coefficients[f"{name}_z"] = float(beta[index + 1])
        predictor_means[name] = float(means[index])
        predictor_stds[name] = float(stds[index])

    return {
        "predictors": predictor_names,
        "coefficients": coefficients,
        "predictor_means": predictor_means,
        "predictor_stds": predictor_stds,
        "r_squared": r_squared,
    }


def sample_image_value(image: np.ndarray, xpix: float, ypix: float) -> float | None:
    if not np.isfinite(xpix) or not np.isfinite(ypix):
        return None
    height, width = image.shape
    if height == 0 or width == 0:
        return None

    x = float(np.clip(xpix, 0.0, width - 1))
    y = float(np.clip(ypix, 0.0, height - 1))
    x0 = int(math.floor(x))
    x1 = min(x0 + 1, width - 1)
    y0 = int(math.floor(y))
    y1 = min(y0 + 1, height - 1)
    dx = x - x0
    dy = y - y0

    weighted_values = []
    weighted_weights = []
    neighbors = [
        ((1.0 - dx) * (1.0 - dy), image[y0, x0]),
        (dx * (1.0 - dy), image[y0, x1]),
        ((1.0 - dx) * dy, image[y1, x0]),
        (dx * dy, image[y1, x1]),
    ]
    for weight, value in neighbors:
        if weight > 0.0 and np.isfinite(value):
            weighted_values.append(float(value))
            weighted_weights.append(float(weight))

    if weighted_weights:
        return float(np.average(weighted_values, weights=weighted_weights))

    xnearest = int(round(x))
    ynearest = int(round(y))
    nearest = float(image[ynearest, xnearest])
    if np.isfinite(nearest):
        return nearest
    return None


def score_cutout(
    residual: np.ndarray,
    footprint: np.ndarray,
    xpix: float,
    ypix: float,
    pixel_scale_arcsec: float,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
) -> dict | None:
    outer_radius_px = ring_outer_arcsec / pixel_scale_arcsec
    xmin = max(int(math.floor(xpix - outer_radius_px)) - 1, 0)
    xmax = min(int(math.ceil(xpix + outer_radius_px)) + 2, residual.shape[1])
    ymin = max(int(math.floor(ypix - outer_radius_px)) - 1, 0)
    ymax = min(int(math.ceil(ypix + outer_radius_px)) + 2, residual.shape[0])
    if xmin >= xmax or ymin >= ymax:
        return None

    cutout = residual[ymin:ymax, xmin:xmax]
    cutout_footprint = footprint[ymin:ymax, xmin:xmax]
    yy, xx = np.indices(cutout.shape)
    xx = xx + xmin
    yy = yy + ymin
    radius_px = np.sqrt((xx - xpix) ** 2 + (yy - ypix) ** 2)

    inner_mask = radius_px <= inner_radius_arcsec / pixel_scale_arcsec
    ring_mask = (
        (radius_px >= ring_inner_arcsec / pixel_scale_arcsec)
        & (radius_px <= outer_radius_px)
    )
    if not inner_mask.any() or not ring_mask.any():
        return None

    inner_values = cutout[inner_mask & cutout_footprint]
    ring_values = cutout[ring_mask & cutout_footprint]
    if inner_values.size == 0 or ring_values.size == 0:
        return None

    inner_mean = float(np.mean(inner_values))
    ring_mean = float(np.mean(ring_values))
    return {
        "inner_mean": inner_mean,
        "ring_mean": ring_mean,
        "sign_flip_score": max(inner_mean, 0.0) + max(-ring_mean, 0.0),
        "contrast": inner_mean - ring_mean,
        "sign_flip_pass": int(inner_mean > 0.0 and ring_mean < 0.0),
    }


def matched_group_difference(
    predictor: np.ndarray,
    response: np.ndarray,
    controls: np.ndarray,
    seed: int,
    permutations: int,
    high_quantile: float = 0.75,
    low_quantile: float = 0.25,
) -> dict:
    result = {
        "high_threshold": float("nan"),
        "low_threshold": float("nan"),
        "n_high_pool": 0,
        "n_low_pool": 0,
        "n_pairs": 0,
        "mean_high_score": float("nan"),
        "mean_matched_low_score": float("nan"),
        "mean_paired_difference": float("nan"),
        "median_paired_difference": float("nan"),
        "mean_control_distance": float("nan"),
        "permutation_p_value": float("nan"),
        "permutations": permutations,
    }

    finite_mask = np.isfinite(predictor) & np.isfinite(response)
    if controls.ndim == 1:
        controls = controls.reshape(-1, 1)
    finite_mask &= np.all(np.isfinite(controls), axis=1)
    predictor = predictor[finite_mask]
    response = response[finite_mask]
    controls = controls[finite_mask]
    if len(predictor) < 8:
        return result

    high_threshold = float(np.quantile(predictor, high_quantile))
    low_threshold = float(np.quantile(predictor, low_quantile))
    result["high_threshold"] = high_threshold
    result["low_threshold"] = low_threshold
    if not high_threshold > low_threshold:
        return result

    high_indices = np.flatnonzero(predictor >= high_threshold)
    low_indices = np.flatnonzero(predictor <= low_threshold)
    low_indices = np.asarray([index for index in low_indices if predictor[index] < high_threshold], dtype=int)
    high_indices = np.asarray([index for index in high_indices if predictor[index] > low_threshold], dtype=int)
    result["n_high_pool"] = int(len(high_indices))
    result["n_low_pool"] = int(len(low_indices))
    if len(high_indices) == 0 or len(low_indices) == 0:
        return result

    control_means = controls.mean(axis=0)
    control_stds = controls.std(axis=0)
    control_stds[control_stds == 0.0] = 1.0
    controls_z = (controls - control_means) / control_stds

    rng = np.random.default_rng(seed)
    high_order = rng.permutation(high_indices)
    available_low = list(low_indices.tolist())
    pairs = []
    pair_distances = []
    for high_index in high_order:
        if not available_low:
            break
        low_array = np.asarray(available_low, dtype=int)
        distances = np.linalg.norm(controls_z[low_array] - controls_z[high_index], axis=1)
        best_position = int(np.argmin(distances))
        low_index = int(low_array[best_position])
        pairs.append((int(high_index), low_index))
        pair_distances.append(float(distances[best_position]))
        del available_low[best_position]

    if not pairs:
        return result

    high_pair_indices = np.asarray([high_index for high_index, _ in pairs], dtype=int)
    low_pair_indices = np.asarray([low_index for _, low_index in pairs], dtype=int)
    paired_differences = response[high_pair_indices] - response[low_pair_indices]
    observed = float(np.mean(paired_differences))
    sign_flips = rng.choice(np.asarray([-1.0, 1.0]), size=(permutations, len(paired_differences)))
    permuted_means = np.mean(sign_flips * paired_differences, axis=1)
    p_value = float((1 + np.sum(np.abs(permuted_means) >= abs(observed))) / (permutations + 1))

    result.update(
        {
            "n_pairs": int(len(pairs)),
            "mean_high_score": float(np.mean(response[high_pair_indices])),
            "mean_matched_low_score": float(np.mean(response[low_pair_indices])),
            "mean_paired_difference": observed,
            "median_paired_difference": float(np.median(paired_differences)),
            "mean_control_distance": float(np.mean(pair_distances)),
            "permutation_p_value": p_value,
        }
    )
    return result


def partial_correlation_and_permutation_pvalue(
    predictor: np.ndarray,
    response: np.ndarray,
    controls: np.ndarray,
    seed: int,
    permutations: int,
) -> dict:
    rng = np.random.default_rng(seed)
    Xc = np.column_stack([np.ones(len(response)), controls])
    predictor_resid = predictor - Xc @ np.linalg.lstsq(Xc, predictor, rcond=None)[0]
    response_resid = response - Xc @ np.linalg.lstsq(Xc, response, rcond=None)[0]
    if np.std(predictor_resid) == 0.0 or np.std(response_resid) == 0.0:
        return {
            "partial_correlation": float("nan"),
            "permutation_p_value": 1.0,
            "permutations": permutations,
        }
    observed = float(np.corrcoef(predictor_resid, response_resid)[0, 1])

    permuted = np.empty(permutations, dtype=float)
    for i in range(permutations):
        shuffled = rng.permutation(predictor)
        shuffled_resid = shuffled - Xc @ np.linalg.lstsq(Xc, shuffled, rcond=None)[0]
        if np.std(shuffled_resid) == 0.0 or np.std(response_resid) == 0.0:
            permuted[i] = 0.0
        else:
            permuted[i] = np.corrcoef(shuffled_resid, response_resid)[0, 1]

    p_value = float((1 + np.sum(np.abs(permuted) >= abs(observed))) / (permutations + 1))
    return {
        "partial_correlation": observed,
        "permutation_p_value": p_value,
        "permutations": permutations,
    }


def random_position_null(
    residual: np.ndarray,
    footprint: np.ndarray,
    wcs: WCS,
    n_positions: int,
    repeats: int,
    pixel_scale_arcsec: float,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)
    ys, xs = np.where(footprint)
    valid_scores = []
    candidate_attempts = 0
    target_pool_size = max(n_positions * 2, 200)
    max_candidate_attempts = max(target_pool_size * 20, 5000)
    while len(valid_scores) < target_pool_size and candidate_attempts < max_candidate_attempts:
        candidate_attempts += 1
        index = int(rng.integers(0, len(xs)))
        xpix = xs[index]
        ypix = ys[index]
        score = score_cutout(
            residual=residual,
            footprint=footprint,
            xpix=float(xpix),
            ypix=float(ypix),
            pixel_scale_arcsec=pixel_scale_arcsec,
            inner_radius_arcsec=inner_radius_arcsec,
            ring_inner_arcsec=ring_inner_arcsec,
            ring_outer_arcsec=ring_outer_arcsec,
        )
        if score is not None:
            valid_scores.append(score["sign_flip_score"])

    if not valid_scores:
        raise ValueError("Random-position null failed to generate valid draws.")

    effective_n_positions = min(n_positions, len(valid_scores))
    draws = []
    for _ in range(repeats):
        sample = rng.choice(valid_scores, size=effective_n_positions, replace=True)
        draws.append(float(np.mean(sample)))

    return {
        "requested_position_count": n_positions,
        "used_position_count": effective_n_positions,
        "valid_position_pool_size": len(valid_scores),
        "draw_count": len(draws),
        "mean_of_draw_means": float(np.mean(draws)),
        "std_of_draw_means": float(np.std(draws)),
        "draw_means": draws,
    }


def combine_pvalues_fisher(pvalues: Iterable[float]) -> dict:
    cleaned = [min(max(float(p), 1e-300), 1.0) for p in pvalues if np.isfinite(p)]
    if not cleaned:
        raise ValueError("No finite p-values supplied for Fisher combination.")
    statistic = float(-2.0 * np.sum(np.log(cleaned)))
    return {
        "k": len(cleaned),
        "statistic": statistic,
        "p_value": float(chi2.sf(statistic, 2 * len(cleaned))),
    }


def fisher_meta_correlation(correlations: Iterable[float], sample_sizes: Iterable[int]) -> dict:
    weighted_terms = []
    total_weight = 0.0
    for correlation, sample_size in zip(correlations, sample_sizes):
        if not np.isfinite(correlation) or sample_size <= 3:
            continue
        clipped = float(np.clip(correlation, -0.999999, 0.999999))
        weight = float(sample_size - 3)
        weighted_terms.append(weight * np.arctanh(clipped))
        total_weight += weight
    if total_weight == 0.0:
        raise ValueError("No valid correlations supplied for meta-analysis.")
    combined_z = float(np.sum(weighted_terms) / total_weight)
    return {
        "combined_correlation": float(np.tanh(combined_z)),
        "weight_sum": total_weight,
        "k": len(weighted_terms),
    }


def simulate_dataset(out_dir: Path, seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_galaxies = 48
    galaxy_ra = 3.55 + 0.06 * rng.random(n_galaxies)
    galaxy_dec = -30.42 + 0.06 * rng.random(n_galaxies)
    metallicity = 8.2 + 0.9 * rng.random(n_galaxies)
    sfr = 10 ** (-0.4 + 2.2 * rng.random(n_galaxies))
    mass = 10 ** (9.0 + 2.0 * rng.random(n_galaxies))

    proxy_strength = (
        0.6 * (metallicity - np.mean(metallicity)) / np.std(metallicity)
        + 0.8 * (np.log10(sfr) - np.mean(np.log10(sfr))) / np.std(np.log10(sfr))
        + 0.2 * (np.log10(mass) - np.mean(np.log10(mass))) / np.std(np.log10(mass))
    )
    amplitudes = 0.12 + 0.08 * proxy_strength

    grid_ra = np.linspace(3.54, 3.62, 72)
    grid_dec = np.linspace(-30.43, -30.35, 72)

    residual_rows = []
    for ra in grid_ra:
        for dec in grid_dec:
            separation_stack = angular_distance_arcsec(ra, dec, galaxy_ra, galaxy_dec)
            field = np.sum(amplitudes * radial_kernel(separation_stack, core_arcsec=2.8, boundary_arcsec=7.5))
            field += rng.normal(0.0, 0.01)
            residual_rows.append({"ra": ra, "dec": dec, "residual": float(field)})

    galaxy_rows = [
        {
            "galaxy_id": f"G{i:03d}",
            "ra": float(galaxy_ra[i]),
            "dec": float(galaxy_dec[i]),
            "metallicity": float(metallicity[i]),
            "sfr": float(sfr[i]),
            "mass": float(mass[i]),
        }
        for i in range(n_galaxies)
    ]

    residual_path = out_dir / "toy_residuals.csv"
    galaxies_path = out_dir / "toy_galaxies.csv"
    write_csv_rows(residual_path, residual_rows, ["ra", "dec", "residual"])
    write_csv_rows(
        galaxies_path,
        galaxy_rows,
        ["galaxy_id", "ra", "dec", "metallicity", "sfr", "mass"],
    )

    result = score_sign_flip(
        galaxies_path=galaxies_path,
        residuals_path=residual_path,
        out_path=out_dir / "toy_results.json",
        inner_radius_arcsec=2.5,
        ring_inner_arcsec=5.0,
        ring_outer_arcsec=9.0,
    )
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "paths": {
            "galaxies": str(galaxies_path),
            "residuals": str(residual_path),
            "results": str(out_dir / "toy_results.json"),
            "scored_rows": str(out_dir / "toy_results.csv"),
        },
        "result": result,
    }


def analyze_hff_cluster_map(
    cluster_key: str,
    catalog_path: Path,
    photometric_proxy_lookup: dict[int, dict],
    kappa_path: Path,
    out_dir: Path,
    result_prefix: str,
    products_metadata: dict,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    smooth_sigma_px: float,
    residual_mode: str,
    radial_bin_arcsec: float,
    cluster_z_min: float,
    cluster_z_max: float,
    permutations: int,
    random_repeats: int,
    seed: int,
) -> dict:
    config = HFF_CLUSTER_CONFIG[cluster_key]

    with fits.open(kappa_path) as hdul:
        kappa = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header
        wcs = WCS(header)

    pixel_scale_deg = header.get("CD1_1", header.get("CDELT1"))
    if pixel_scale_deg is None:
        raise ValueError("Could not determine pixel scale from FITS header.")
    pixel_scale_arcsec = abs(float(pixel_scale_deg)) * 3600.0
    footprint = np.isfinite(kappa)
    cluster_center_ra = float(header["CRVAL1"])
    cluster_center_dec = float(header["CRVAL2"])
    center_xpix, center_ypix = wcs.world_to_pixel_values(cluster_center_ra, cluster_center_dec)
    residual, residual_model = build_residual_map(
        kappa=kappa,
        footprint=footprint,
        center_xpix=float(center_xpix),
        center_ypix=float(center_ypix),
        pixel_scale_arcsec=pixel_scale_arcsec,
        residual_mode=residual_mode,
        smooth_sigma_px=smooth_sigma_px,
        radial_bin_arcsec=radial_bin_arcsec,
    )

    with fits.open(catalog_path) as hdul:
        table = hdul[1].data
        rows = []
        design = []
        response = []
        footprint = np.isfinite(residual)
        for record in table:
            mass = float(record["MASS_NEB"])
            sfr = float(record["SFR_NEB"])
            zbest = float(record["ZBEST"])
            magnif = float(record["MAGNIF"])
            ra = float(record["alpha_j2000"])
            dec = float(record["delta_j2000"])

            if not np.isfinite(mass) or not np.isfinite(sfr) or mass <= 0.0 or sfr <= 0.0:
                continue

            xpix, ypix = wcs.world_to_pixel_values(ra, dec)
            if not np.isfinite(xpix) or not np.isfinite(ypix):
                continue
            if xpix < 0 or ypix < 0 or xpix >= residual.shape[1] or ypix >= residual.shape[0]:
                continue

            score = score_cutout(
                residual=residual,
                footprint=footprint,
                xpix=xpix,
                ypix=ypix,
                pixel_scale_arcsec=pixel_scale_arcsec,
                inner_radius_arcsec=inner_radius_arcsec,
                ring_inner_arcsec=ring_inner_arcsec,
                ring_outer_arcsec=ring_outer_arcsec,
            )
            if score is None:
                continue
            local_kappa = sample_image_value(kappa, xpix, ypix)
            if local_kappa is None:
                continue
            inner_mean = score["inner_mean"]
            ring_mean = score["ring_mean"]
            sign_flip_score = score["sign_flip_score"]
            contrast = score["contrast"]
            sign_flip_pass = score["sign_flip_pass"]
            cluster_radius_arcsec = angular_distance_arcsec_scalar(
                cluster_center_ra,
                cluster_center_dec,
                ra,
                dec,
            )
            is_cluster_member = int(cluster_z_min <= zbest <= cluster_z_max)
            photometric_info = photometric_proxy_lookup.get(int(record["IDcat"]), {})

            row = {
                "IDcat": int(record["IDcat"]),
                "ra": ra,
                "dec": dec,
                "zbest": zbest,
                "magnif": magnif,
                "mass_neb": mass,
                "sfr_neb": sfr,
                "cluster_radius_arcsec": cluster_radius_arcsec,
                "local_kappa": local_kappa,
                "is_cluster_member": is_cluster_member,
                "photometric_proxy": (
                    float(photometric_info["photometric_proxy"])
                    if photometric_info.get("photometric_proxy") is not None
                    else float("nan")
                ),
                "photometric_proxy_band_count": photometric_info.get("photometric_proxy_band_count", 0),
                "photometric_proxy_match_sep_arcsec": photometric_info.get("photometric_proxy_match_sep_arcsec", float("nan")),
                "photometric_proxy_ebv": photometric_info.get("photometric_proxy_ebv", float("nan")),
                "photometric_proxy_nb_used": photometric_info.get("photometric_proxy_nb_used", 0),
                "inner_mean": inner_mean,
                "ring_mean": ring_mean,
                "sign_flip_score": sign_flip_score,
                "contrast": contrast,
                "sign_flip_pass": sign_flip_pass,
            }
            rows.append(row)
            design.append([mass, sfr, zbest, math.log10(max(magnif, 1e-6)), cluster_radius_arcsec, local_kappa])
            response.append(sign_flip_score)

        if not rows:
            raise ValueError(f"No {config['label']} sources survived the first-pass filters.")

    regression = fit_linear_model(
        np.asarray(design, dtype=float),
        np.asarray(response, dtype=float),
        ["mass_neb", "sfr_neb", "zbest", "log10_magnif", "cluster_radius_arcsec", "local_kappa"],
    )

    cluster_rows = [row for row in rows if row["is_cluster_member"] == 1]
    if not cluster_rows:
        raise ValueError("No cluster members survived the first-pass filters.")

    cluster_score = np.asarray([row["sign_flip_score"] for row in cluster_rows], dtype=float)
    cluster_mass = np.asarray([row["mass_neb"] for row in cluster_rows], dtype=float)
    cluster_sfr = np.asarray([row["sfr_neb"] for row in cluster_rows], dtype=float)
    cluster_zbest = np.asarray([row["zbest"] for row in cluster_rows], dtype=float)
    cluster_logmagnif = np.log10(np.maximum(np.asarray([row["magnif"] for row in cluster_rows], dtype=float), 1e-6))
    cluster_radius = np.asarray([row["cluster_radius_arcsec"] for row in cluster_rows], dtype=float)
    cluster_local_kappa = np.asarray([row["local_kappa"] for row in cluster_rows], dtype=float)

    cluster_regression = fit_linear_model(
        np.column_stack([cluster_mass, cluster_sfr, cluster_zbest, cluster_logmagnif, cluster_radius, cluster_local_kappa]),
        cluster_score,
        ["mass_neb", "sfr_neb", "zbest", "log10_magnif", "cluster_radius_arcsec", "local_kappa"],
    )

    controls_for_mass = np.column_stack([cluster_sfr, cluster_zbest, cluster_logmagnif, cluster_radius, cluster_local_kappa])
    controls_for_sfr = np.column_stack([cluster_mass, cluster_zbest, cluster_logmagnif, cluster_radius, cluster_local_kappa])
    controls_for_local_kappa = np.column_stack([cluster_mass, cluster_sfr, cluster_zbest, cluster_logmagnif, cluster_radius])
    mass_partial = partial_correlation_and_permutation_pvalue(
        predictor=cluster_mass,
        response=cluster_score,
        controls=controls_for_mass,
        seed=seed,
        permutations=permutations,
    )
    sfr_partial = partial_correlation_and_permutation_pvalue(
        predictor=cluster_sfr,
        response=cluster_score,
        controls=controls_for_sfr,
        seed=seed + 1,
        permutations=permutations,
    )
    local_kappa_partial = partial_correlation_and_permutation_pvalue(
        predictor=cluster_local_kappa,
        response=cluster_score,
        controls=controls_for_local_kappa,
        seed=seed + 2,
        permutations=permutations,
    )
    matched_group_effects = {
        "mass_high_vs_low": matched_group_difference(
            predictor=cluster_mass,
            response=cluster_score,
            controls=controls_for_mass,
            seed=seed + 20,
            permutations=permutations,
        ),
        "sfr_high_vs_low": matched_group_difference(
            predictor=cluster_sfr,
            response=cluster_score,
            controls=controls_for_sfr,
            seed=seed + 21,
            permutations=permutations,
        ),
        "local_kappa_high_vs_low": matched_group_difference(
            predictor=cluster_local_kappa,
            response=cluster_score,
            controls=controls_for_local_kappa,
            seed=seed + 22,
            permutations=permutations,
        ),
    }

    mass_threshold = float(np.quantile(cluster_mass, 0.75))
    sfr_threshold = float(np.quantile(cluster_sfr, 0.75))
    high_mass_mask = cluster_mass >= mass_threshold
    high_sfr_mask = cluster_sfr >= sfr_threshold

    proxy_cluster_rows = [row for row in cluster_rows if np.isfinite(float(row["photometric_proxy"]))]
    if proxy_cluster_rows:
        proxy_cluster_score = np.asarray([row["sign_flip_score"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_mass = np.asarray([row["mass_neb"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_sfr = np.asarray([row["sfr_neb"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster = np.asarray([row["photometric_proxy"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_zbest = np.asarray([row["zbest"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_logmagnif = np.log10(
            np.maximum(np.asarray([row["magnif"] for row in proxy_cluster_rows], dtype=float), 1e-6)
        )
        proxy_cluster_radius = np.asarray([row["cluster_radius_arcsec"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_local_kappa = np.asarray([row["local_kappa"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_match_sep = np.asarray([row["photometric_proxy_match_sep_arcsec"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_ebv = np.asarray([row["photometric_proxy_ebv"] for row in proxy_cluster_rows], dtype=float)
        proxy_cluster_band_count = np.asarray([row["photometric_proxy_band_count"] for row in proxy_cluster_rows], dtype=float)

        proxy_regression = fit_linear_model(
            np.column_stack(
                [
                    proxy_cluster_mass,
                    proxy_cluster_sfr,
                    proxy_cluster,
                    proxy_cluster_zbest,
                    proxy_cluster_logmagnif,
                    proxy_cluster_radius,
                    proxy_cluster_local_kappa,
                ]
            ),
            proxy_cluster_score,
            [
                "mass_neb",
                "sfr_neb",
                "photometric_proxy",
                "zbest",
                "log10_magnif",
                "cluster_radius_arcsec",
                "local_kappa",
            ],
        )

        proxy_controls = np.column_stack(
            [
                proxy_cluster_mass,
                proxy_cluster_sfr,
                proxy_cluster_zbest,
                proxy_cluster_logmagnif,
                proxy_cluster_radius,
                proxy_cluster_local_kappa,
            ]
        )
        proxy_partial = partial_correlation_and_permutation_pvalue(
            predictor=proxy_cluster,
            response=proxy_cluster_score,
            controls=proxy_controls,
            seed=seed + 3,
            permutations=permutations,
        )
        proxy_matched = matched_group_difference(
            predictor=proxy_cluster,
            response=proxy_cluster_score,
            controls=proxy_controls,
            seed=seed + 23,
            permutations=permutations,
        )
        proxy_threshold = float(np.quantile(proxy_cluster, 0.75))
        high_proxy_mask = proxy_cluster >= proxy_threshold
        photometric_proxy_analysis = {
            "n_cluster_members_with_proxy": len(proxy_cluster_rows),
            "mean_proxy": float(np.mean(proxy_cluster)),
            "std_proxy": float(np.std(proxy_cluster)),
            "median_proxy_match_sep_arcsec": float(np.nanmedian(proxy_cluster_match_sep)),
            "median_proxy_band_count": float(np.nanmedian(proxy_cluster_band_count)),
            "mean_proxy_ebv": float(np.nanmean(proxy_cluster_ebv)),
            "high_proxy_threshold": proxy_threshold,
            "high_proxy_mean_score": float(np.mean(proxy_cluster_score[high_proxy_mask])),
            "high_proxy_pass_rate": float(
                np.mean(np.asarray([row["sign_flip_pass"] for row in proxy_cluster_rows], dtype=float)[high_proxy_mask])
            ),
            "regression": proxy_regression,
            "proxy_partial": proxy_partial,
            "matched_high_vs_low": proxy_matched,
            "spearman": {
                "proxy_vs_score": float(spearmanr(proxy_cluster, proxy_cluster_score).statistic),
            },
        }
    else:
        photometric_proxy_analysis = {
            "n_cluster_members_with_proxy": 0,
            "mean_proxy": float("nan"),
            "std_proxy": float("nan"),
            "median_proxy_match_sep_arcsec": float("nan"),
            "median_proxy_band_count": float("nan"),
            "mean_proxy_ebv": float("nan"),
            "high_proxy_threshold": float("nan"),
            "high_proxy_mean_score": float("nan"),
            "high_proxy_pass_rate": float("nan"),
            "regression": {},
            "proxy_partial": {"partial_correlation": float("nan"), "permutation_p_value": float("nan"), "permutations": permutations},
            "matched_high_vs_low": {
                "high_threshold": float("nan"),
                "low_threshold": float("nan"),
                "n_high_pool": 0,
                "n_low_pool": 0,
                "n_pairs": 0,
                "mean_high_score": float("nan"),
                "mean_matched_low_score": float("nan"),
                "mean_paired_difference": float("nan"),
                "median_paired_difference": float("nan"),
                "mean_control_distance": float("nan"),
                "permutation_p_value": float("nan"),
                "permutations": permutations,
            },
            "spearman": {"proxy_vs_score": float("nan")},
        }

    cluster_mean_score = float(np.mean(cluster_score))
    try:
        random_null = random_position_null(
            residual=residual,
            footprint=footprint,
            wcs=wcs,
            n_positions=len(cluster_rows),
            repeats=random_repeats,
            pixel_scale_arcsec=pixel_scale_arcsec,
            inner_radius_arcsec=inner_radius_arcsec,
            ring_inner_arcsec=ring_inner_arcsec,
            ring_outer_arcsec=ring_outer_arcsec,
            seed=seed,
        )
        random_draws = np.asarray(random_null["draw_means"], dtype=float)
        random_p_value = float((1 + np.sum(random_draws >= cluster_mean_score)) / (len(random_draws) + 1))
    except ValueError:
        random_null = {
            "requested_position_count": len(cluster_rows),
            "used_position_count": 0,
            "valid_position_pool_size": 0,
            "mean_of_draw_means": float("nan"),
            "std_of_draw_means": float("nan"),
            "draw_count": 0,
        }
        random_p_value = float("nan")

    write_csv_rows(
        out_dir / f"{result_prefix}_rows.csv",
        rows,
        list(rows[0].keys()),
    )
    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cluster": config["label"],
        "cluster_key": cluster_key,
        "products": products_metadata,
        "parameters": {
            "inner_radius_arcsec": inner_radius_arcsec,
            "ring_inner_arcsec": ring_inner_arcsec,
            "ring_outer_arcsec": ring_outer_arcsec,
            "smooth_sigma_px": smooth_sigma_px,
            "residual_mode": residual_mode,
            "radial_bin_arcsec": radial_bin_arcsec,
            "pixel_scale_arcsec": pixel_scale_arcsec,
            "cluster_member_z_window": [cluster_z_min, cluster_z_max],
        },
        "n_scored": len(rows),
        "sign_flip_pass_rate": float(np.mean([row["sign_flip_pass"] for row in rows])),
        "residual_summary": {
            "kappa_mean": float(np.mean(kappa)),
            "residual_mean": float(np.mean(residual)),
            "residual_std": float(np.std(residual)),
            **residual_model,
        },
        "regression": regression,
        "cluster_member_analysis": {
            "n_cluster_members": len(cluster_rows),
            "mean_score": cluster_mean_score,
            "pass_rate": float(np.mean([row["sign_flip_pass"] for row in cluster_rows])),
            "high_mass_threshold": mass_threshold,
            "high_mass_mean_score": float(np.mean(cluster_score[high_mass_mask])),
            "high_mass_pass_rate": float(np.mean(np.asarray([row["sign_flip_pass"] for row in cluster_rows], dtype=float)[high_mass_mask])),
            "high_sfr_threshold": sfr_threshold,
            "high_sfr_mean_score": float(np.mean(cluster_score[high_sfr_mask])),
            "high_sfr_pass_rate": float(np.mean(np.asarray([row["sign_flip_pass"] for row in cluster_rows], dtype=float)[high_sfr_mask])),
            "regression": cluster_regression,
            "mass_partial": mass_partial,
            "sfr_partial": sfr_partial,
            "local_kappa_partial": local_kappa_partial,
            "matched_group_effects": matched_group_effects,
            "photometric_proxy_analysis": photometric_proxy_analysis,
            "random_position_null": {
                "mean_of_draw_means": random_null["mean_of_draw_means"],
                "std_of_draw_means": random_null["std_of_draw_means"],
                "p_value_cluster_mean_gt_random": random_p_value,
                "draw_count": random_null["draw_count"],
            },
            "spearman": {
                "mass_vs_score": float(spearmanr(cluster_mass, cluster_score).statistic),
                "sfr_vs_score": float(spearmanr(cluster_sfr, cluster_score).statistic),
                "radius_vs_score": float(spearmanr(cluster_radius, cluster_score).statistic),
                "local_kappa_vs_score": float(spearmanr(cluster_local_kappa, cluster_score).statistic),
            },
        },
    }
    (out_dir / f"{result_prefix}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def analyze_hff_cluster_firstpass(
    cluster_key: str,
    out_dir: Path,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    smooth_sigma_px: float,
    residual_mode: str,
    radial_bin_arcsec: float,
    cluster_z_min: float,
    cluster_z_max: float,
    permutations: int,
    random_repeats: int,
    seed: int,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = out_dir / "inputs"
    config = HFF_CLUSTER_CONFIG[cluster_key]
    products = OFFICIAL_PRODUCTS[cluster_key]
    proxy_catalog_path = download_file(
        products[config["proxy_key"]],
        inputs_dir / config["proxy_filename"],
    )
    catalog_path = download_file(
        products[config["catalog_key"]],
        inputs_dir / config["catalog_filename"],
    )
    kappa_path = download_file(
        products[config["kappa_key"]],
        inputs_dir / config["kappa_filename"],
    )
    photometric_proxy_lookup = load_photometric_proxy_lookup(catalog_path, proxy_catalog_path)
    return analyze_hff_cluster_map(
        cluster_key=cluster_key,
        catalog_path=catalog_path,
        photometric_proxy_lookup=photometric_proxy_lookup,
        kappa_path=kappa_path,
        out_dir=out_dir,
        result_prefix=f"{cluster_key}_firstpass",
        products_metadata={
            "proxy_catalog_url": products[config["proxy_key"]],
            "proxy_catalog_filename": proxy_catalog_path.name,
            "proxy_catalog_name": "BUFFALO v2.0",
            "properties_catalog_url": products[config["catalog_key"]],
            "properties_catalog_filename": catalog_path.name,
            "kappa_map_url": products[config["kappa_key"]],
            "kappa_map_filename": kappa_path.name,
            "model_family": "cats",
        },
        inner_radius_arcsec=inner_radius_arcsec,
        ring_inner_arcsec=ring_inner_arcsec,
        ring_outer_arcsec=ring_outer_arcsec,
        smooth_sigma_px=smooth_sigma_px,
        residual_mode=residual_mode,
        radial_bin_arcsec=radial_bin_arcsec,
        cluster_z_min=cluster_z_min,
        cluster_z_max=cluster_z_max,
        permutations=permutations,
        random_repeats=random_repeats,
        seed=seed,
    )


def summarize_cluster_result(result: dict) -> dict:
    cluster = result["cluster_member_analysis"]
    proxy = cluster.get("photometric_proxy_analysis", {})
    matched = cluster.get("matched_group_effects", {})
    return {
        "cluster": result["cluster"],
        "cluster_key": result["cluster_key"],
        "model_family": result["products"].get("model_family", "unknown"),
        "residual_mode": result["parameters"]["residual_mode"],
        "inner_radius_arcsec": result["parameters"]["inner_radius_arcsec"],
        "ring_inner_arcsec": result["parameters"]["ring_inner_arcsec"],
        "ring_outer_arcsec": result["parameters"]["ring_outer_arcsec"],
        "smooth_sigma_px": result["parameters"]["smooth_sigma_px"],
        "radial_bin_arcsec": result["parameters"]["radial_bin_arcsec"],
        "n_scored": result["n_scored"],
        "n_cluster_members": cluster["n_cluster_members"],
        "mean_score": cluster["mean_score"],
        "pass_rate": cluster["pass_rate"],
        "high_mass_mean_score": cluster["high_mass_mean_score"],
        "regression_r_squared": cluster["regression"]["r_squared"],
        "mass_partial_correlation": cluster["mass_partial"]["partial_correlation"],
        "mass_partial_p_value": cluster["mass_partial"]["permutation_p_value"],
        "sfr_partial_correlation": cluster["sfr_partial"]["partial_correlation"],
        "sfr_partial_p_value": cluster["sfr_partial"]["permutation_p_value"],
        "local_kappa_partial_correlation": cluster["local_kappa_partial"]["partial_correlation"],
        "local_kappa_partial_p_value": cluster["local_kappa_partial"]["permutation_p_value"],
        "mass_matched_mean_difference": matched.get("mass_high_vs_low", {}).get("mean_paired_difference", float("nan")),
        "mass_matched_p_value": matched.get("mass_high_vs_low", {}).get("permutation_p_value", float("nan")),
        "sfr_matched_mean_difference": matched.get("sfr_high_vs_low", {}).get("mean_paired_difference", float("nan")),
        "sfr_matched_p_value": matched.get("sfr_high_vs_low", {}).get("permutation_p_value", float("nan")),
        "n_cluster_members_with_proxy": proxy.get("n_cluster_members_with_proxy", 0),
        "median_proxy_match_sep_arcsec": proxy.get("median_proxy_match_sep_arcsec", float("nan")),
        "median_proxy_band_count": proxy.get("median_proxy_band_count", float("nan")),
        "mean_proxy_ebv": proxy.get("mean_proxy_ebv", float("nan")),
        "photometric_proxy_partial_correlation": proxy.get("proxy_partial", {}).get("partial_correlation", float("nan")),
        "photometric_proxy_partial_p_value": proxy.get("proxy_partial", {}).get("permutation_p_value", float("nan")),
        "photometric_proxy_matched_mean_difference": proxy.get("matched_high_vs_low", {}).get("mean_paired_difference", float("nan")),
        "photometric_proxy_matched_p_value": proxy.get("matched_high_vs_low", {}).get("permutation_p_value", float("nan")),
        "random_null_p_value": cluster["random_position_null"]["p_value_cluster_mean_gt_random"],
    }


def analyze_hff_cluster_model_ensemble(
    cluster_key: str,
    out_dir: Path,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    smooth_sigma_px: float,
    residual_mode: str,
    radial_bin_arcsec: float,
    cluster_z_min: float,
    cluster_z_max: float,
    permutations: int,
    random_repeats: int,
    seed: int,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = out_dir / "inputs"
    config = HFF_CLUSTER_CONFIG[cluster_key]
    products = OFFICIAL_PRODUCTS[cluster_key]
    proxy_catalog_path = download_file(
        products[config["proxy_key"]],
        inputs_dir / config["proxy_filename"],
    )
    catalog_path = download_file(
        products[config["catalog_key"]],
        inputs_dir / config["catalog_filename"],
    )
    photometric_proxy_lookup = load_photometric_proxy_lookup(catalog_path, proxy_catalog_path)

    model_summaries = []
    run_counter = 0
    for model_family in HFF_MODEL_VERSIONS[cluster_key]:
        run_counter += 1
        result_prefix = f"{cluster_key}_{model_family}"
        result_path = out_dir / f"{result_prefix}.json"
        if result_path.exists():
            model_summaries.append(summarize_cluster_result(json.loads(result_path.read_text(encoding="utf-8"))))
            continue
        kappa_url = frontier_model_kappa_url(cluster_key, model_family)
        kappa_path = download_file(
            kappa_url,
            inputs_dir / f"{cluster_key}_{model_family}_kappa.fits",
        )
        result = analyze_hff_cluster_map(
            cluster_key=cluster_key,
            catalog_path=catalog_path,
            photometric_proxy_lookup=photometric_proxy_lookup,
            kappa_path=kappa_path,
            out_dir=out_dir,
            result_prefix=result_prefix,
            products_metadata={
                "proxy_catalog_url": products[config["proxy_key"]],
                "proxy_catalog_filename": proxy_catalog_path.name,
                "proxy_catalog_name": "BUFFALO v2.0",
                "properties_catalog_url": products[config["catalog_key"]],
                "properties_catalog_filename": catalog_path.name,
                "kappa_map_url": kappa_url,
                "kappa_map_filename": kappa_path.name,
                "model_family": model_family,
            },
            inner_radius_arcsec=inner_radius_arcsec,
            ring_inner_arcsec=ring_inner_arcsec,
            ring_outer_arcsec=ring_outer_arcsec,
            smooth_sigma_px=smooth_sigma_px,
            residual_mode=residual_mode,
            radial_bin_arcsec=radial_bin_arcsec,
            cluster_z_min=cluster_z_min,
            cluster_z_max=cluster_z_max,
            permutations=permutations,
            random_repeats=random_repeats,
            seed=seed + run_counter,
        )
        model_summaries.append(summarize_cluster_result(result))

    write_csv_rows(out_dir / "model_ensemble_summary.csv", model_summaries, list(model_summaries[0].keys()))

    mass_corr = np.asarray([row["mass_partial_correlation"] for row in model_summaries], dtype=float)
    mass_p = np.asarray([row["mass_partial_p_value"] for row in model_summaries], dtype=float)
    sfr_corr = np.asarray([row["sfr_partial_correlation"] for row in model_summaries], dtype=float)
    sfr_p = np.asarray([row["sfr_partial_p_value"] for row in model_summaries], dtype=float)
    local_kappa_corr = np.asarray([row["local_kappa_partial_correlation"] for row in model_summaries], dtype=float)
    local_kappa_p = np.asarray([row["local_kappa_partial_p_value"] for row in model_summaries], dtype=float)
    proxy_corr = np.asarray([row["photometric_proxy_partial_correlation"] for row in model_summaries], dtype=float)
    proxy_p = np.asarray([row["photometric_proxy_partial_p_value"] for row in model_summaries], dtype=float)
    random_p = np.asarray([row["random_null_p_value"] for row in model_summaries], dtype=float)
    member_counts = np.asarray([row["n_cluster_members"] for row in model_summaries], dtype=float)
    mean_scores = np.asarray([row["mean_score"] for row in model_summaries], dtype=float)
    pass_rates = np.asarray([row["pass_rate"] for row in model_summaries], dtype=float)

    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cluster": config["label"],
        "cluster_key": cluster_key,
        "parameters": {
            "inner_radius_arcsec": inner_radius_arcsec,
            "ring_inner_arcsec": ring_inner_arcsec,
            "ring_outer_arcsec": ring_outer_arcsec,
            "smooth_sigma_px": smooth_sigma_px,
            "residual_mode": residual_mode,
            "radial_bin_arcsec": radial_bin_arcsec,
            "cluster_member_z_window": [cluster_z_min, cluster_z_max],
            "permutations": permutations,
            "random_repeats": random_repeats,
        },
        "models": model_summaries,
        "ensemble_summary": {
            "model_count": len(model_summaries),
            "weighted_mean_score": float(np.average(mean_scores, weights=member_counts)),
            "weighted_pass_rate": float(np.average(pass_rates, weights=member_counts)),
            "mass_positive_fraction": float(np.mean(mass_corr > 0.0)),
            "mass_significant_fraction": float(np.mean((mass_corr > 0.0) & (mass_p < 0.05))),
            "mass_median_partial_correlation": float(np.median(mass_corr)),
            "mass_min_partial_correlation": float(np.min(mass_corr)),
            "mass_max_partial_correlation": float(np.max(mass_corr)),
            "sfr_positive_fraction": float(np.mean(sfr_corr > 0.0)),
            "sfr_significant_fraction": float(np.mean((sfr_corr > 0.0) & (sfr_p < 0.05))),
            "sfr_median_partial_correlation": float(np.median(sfr_corr)),
            "sfr_min_partial_correlation": float(np.min(sfr_corr)),
            "sfr_max_partial_correlation": float(np.max(sfr_corr)),
            "local_kappa_positive_fraction": float(np.mean(local_kappa_corr > 0.0)),
            "local_kappa_significant_fraction": float(np.mean((local_kappa_corr > 0.0) & (local_kappa_p < 0.05))),
            "local_kappa_median_partial_correlation": float(np.median(local_kappa_corr)),
            "local_kappa_min_partial_correlation": float(np.min(local_kappa_corr)),
            "local_kappa_max_partial_correlation": float(np.max(local_kappa_corr)),
            "photometric_proxy_positive_fraction": float(np.mean(proxy_corr > 0.0)),
            "photometric_proxy_significant_fraction": float(np.mean((proxy_corr > 0.0) & (proxy_p < 0.05))),
            "photometric_proxy_median_partial_correlation": float(np.nanmedian(proxy_corr)),
            "photometric_proxy_min_partial_correlation": float(np.nanmin(proxy_corr)),
            "photometric_proxy_max_partial_correlation": float(np.nanmax(proxy_corr)),
            "random_null_significant_fraction": float(np.mean(random_p < 0.05)),
        },
    }
    (out_dir / "model_ensemble.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def analyze_combined_hff_firstpass(
    cluster_keys: list[str],
    out_dir: Path,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    smooth_sigma_px: float,
    residual_mode: str,
    radial_bin_arcsec: float,
    cluster_z_min: float,
    cluster_z_max: float,
    permutations: int,
    random_repeats: int,
    seed: int,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    cluster_results = []
    cluster_summaries = []
    for offset, cluster_key in enumerate(cluster_keys):
        result = analyze_hff_cluster_firstpass(
            cluster_key=cluster_key,
            out_dir=out_dir / cluster_key,
            inner_radius_arcsec=inner_radius_arcsec,
            ring_inner_arcsec=ring_inner_arcsec,
            ring_outer_arcsec=ring_outer_arcsec,
            smooth_sigma_px=smooth_sigma_px,
            residual_mode=residual_mode,
            radial_bin_arcsec=radial_bin_arcsec,
            cluster_z_min=cluster_z_min,
            cluster_z_max=cluster_z_max,
            permutations=permutations,
            random_repeats=random_repeats,
            seed=seed + offset * 1000,
        )
        cluster_results.append(result)
        cluster_summaries.append(summarize_cluster_result(result))

    member_counts = np.asarray([item["n_cluster_members"] for item in cluster_summaries], dtype=float)
    mean_scores = np.asarray([item["mean_score"] for item in cluster_summaries], dtype=float)
    pass_rates = np.asarray([item["pass_rate"] for item in cluster_summaries], dtype=float)
    mass_corrs = [item["mass_partial_correlation"] for item in cluster_summaries]
    mass_ps = [item["mass_partial_p_value"] for item in cluster_summaries]
    sfr_corrs = [item["sfr_partial_correlation"] for item in cluster_summaries]
    sfr_ps = [item["sfr_partial_p_value"] for item in cluster_summaries]
    local_kappa_corrs = [item["local_kappa_partial_correlation"] for item in cluster_summaries]
    local_kappa_ps = [item["local_kappa_partial_p_value"] for item in cluster_summaries]
    proxy_corrs = [item["photometric_proxy_partial_correlation"] for item in cluster_summaries]
    proxy_ps = [item["photometric_proxy_partial_p_value"] for item in cluster_summaries]
    random_ps = [item["random_null_p_value"] for item in cluster_summaries]

    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "clusters": cluster_summaries,
        "parameters": {
            "cluster_keys": cluster_keys,
            "inner_radius_arcsec": inner_radius_arcsec,
            "ring_inner_arcsec": ring_inner_arcsec,
            "ring_outer_arcsec": ring_outer_arcsec,
            "smooth_sigma_px": smooth_sigma_px,
            "residual_mode": residual_mode,
            "radial_bin_arcsec": radial_bin_arcsec,
            "cluster_member_z_window": [cluster_z_min, cluster_z_max],
            "permutations": permutations,
            "random_repeats": random_repeats,
        },
        "combined_summary": {
            "cluster_count": len(cluster_summaries),
            "combined_cluster_members": int(np.sum(member_counts)),
            "weighted_mean_score": float(np.average(mean_scores, weights=member_counts)),
            "weighted_pass_rate": float(np.average(pass_rates, weights=member_counts)),
            "mass_partial_meta": fisher_meta_correlation(mass_corrs, member_counts.astype(int)),
            "mass_partial_fisher": combine_pvalues_fisher(mass_ps),
            "sfr_partial_meta": fisher_meta_correlation(sfr_corrs, member_counts.astype(int)),
            "sfr_partial_fisher": combine_pvalues_fisher(sfr_ps),
            "local_kappa_partial_meta": fisher_meta_correlation(local_kappa_corrs, member_counts.astype(int)),
            "local_kappa_partial_fisher": combine_pvalues_fisher(local_kappa_ps),
            "photometric_proxy_meta": fisher_meta_correlation(proxy_corrs, member_counts.astype(int)),
            "photometric_proxy_fisher": combine_pvalues_fisher(proxy_ps),
            "random_null_fisher": combine_pvalues_fisher(random_ps),
            "replication": {
                "mass_positive_cluster_count": int(sum(value > 0.0 for value in mass_corrs)),
                "mass_significant_cluster_count": int(
                    sum(corr > 0.0 and p_value < 0.05 for corr, p_value in zip(mass_corrs, mass_ps))
                ),
                "sfr_positive_cluster_count": int(sum(value > 0.0 for value in sfr_corrs)),
                "sfr_significant_cluster_count": int(
                    sum(corr > 0.0 and p_value < 0.05 for corr, p_value in zip(sfr_corrs, sfr_ps))
                ),
                "local_kappa_positive_cluster_count": int(sum(value > 0.0 for value in local_kappa_corrs)),
                "local_kappa_significant_cluster_count": int(
                    sum(corr > 0.0 and p_value < 0.05 for corr, p_value in zip(local_kappa_corrs, local_kappa_ps))
                ),
                "photometric_proxy_positive_cluster_count": int(sum(value > 0.0 for value in proxy_corrs)),
                "photometric_proxy_significant_cluster_count": int(
                    sum(corr > 0.0 and p_value < 0.05 for corr, p_value in zip(proxy_corrs, proxy_ps))
                ),
                "random_null_significant_cluster_count": int(sum(p_value < 0.05 for p_value in random_ps)),
            },
        },
    }
    (out_dir / "combined_firstpass.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def run_full_experiment(
    cluster_keys: list[str],
    out_dir: Path,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    smooth_sigma_px: float,
    residual_mode: str,
    radial_bin_arcsec: float,
    cluster_z_min: float,
    cluster_z_max: float,
    firstpass_permutations: int,
    firstpass_random_repeats: int,
    ensemble_permutations: int,
    ensemble_random_repeats: int,
    seed: int,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)

    firstpass = analyze_combined_hff_firstpass(
        cluster_keys=cluster_keys,
        out_dir=out_dir / "primary_firstpass",
        inner_radius_arcsec=inner_radius_arcsec,
        ring_inner_arcsec=ring_inner_arcsec,
        ring_outer_arcsec=ring_outer_arcsec,
        smooth_sigma_px=smooth_sigma_px,
        residual_mode=residual_mode,
        radial_bin_arcsec=radial_bin_arcsec,
        cluster_z_min=cluster_z_min,
        cluster_z_max=cluster_z_max,
        permutations=firstpass_permutations,
        random_repeats=firstpass_random_repeats,
        seed=seed,
    )

    ensemble_results = []
    for offset, cluster_key in enumerate(cluster_keys):
        ensemble_results.append(
            analyze_hff_cluster_model_ensemble(
                cluster_key=cluster_key,
                out_dir=out_dir / f"{cluster_key}_model_ensemble",
                inner_radius_arcsec=inner_radius_arcsec,
                ring_inner_arcsec=ring_inner_arcsec,
                ring_outer_arcsec=ring_outer_arcsec,
                smooth_sigma_px=smooth_sigma_px,
                residual_mode=residual_mode,
                radial_bin_arcsec=radial_bin_arcsec,
                cluster_z_min=cluster_z_min,
                cluster_z_max=cluster_z_max,
                permutations=ensemble_permutations,
                random_repeats=ensemble_random_repeats,
                seed=seed + 10000 + offset * 1000,
            )
        )

    ensemble_summaries = []
    all_model_rows = []
    for result in ensemble_results:
        summary = result["ensemble_summary"]
        ensemble_summaries.append(
            {
                "cluster": result["cluster"],
                "cluster_key": result["cluster_key"],
                **summary,
            }
        )
        all_model_rows.extend(result["models"])

    write_csv_rows(out_dir / "full_experiment_model_runs.csv", all_model_rows, list(all_model_rows[0].keys()))
    write_csv_rows(out_dir / "full_experiment_cluster_consensus.csv", ensemble_summaries, list(ensemble_summaries[0].keys()))

    mass_corr = np.asarray([row["mass_partial_correlation"] for row in all_model_rows], dtype=float)
    mass_p = np.asarray([row["mass_partial_p_value"] for row in all_model_rows], dtype=float)
    sfr_corr = np.asarray([row["sfr_partial_correlation"] for row in all_model_rows], dtype=float)
    sfr_p = np.asarray([row["sfr_partial_p_value"] for row in all_model_rows], dtype=float)
    local_kappa_corr = np.asarray([row["local_kappa_partial_correlation"] for row in all_model_rows], dtype=float)
    local_kappa_p = np.asarray([row["local_kappa_partial_p_value"] for row in all_model_rows], dtype=float)
    proxy_corr = np.asarray([row["photometric_proxy_partial_correlation"] for row in all_model_rows], dtype=float)
    proxy_p = np.asarray([row["photometric_proxy_partial_p_value"] for row in all_model_rows], dtype=float)
    random_p = np.asarray([row["random_null_p_value"] for row in all_model_rows], dtype=float)

    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "cluster_keys": cluster_keys,
            "inner_radius_arcsec": inner_radius_arcsec,
            "ring_inner_arcsec": ring_inner_arcsec,
            "ring_outer_arcsec": ring_outer_arcsec,
            "smooth_sigma_px": smooth_sigma_px,
            "residual_mode": residual_mode,
            "radial_bin_arcsec": radial_bin_arcsec,
            "cluster_member_z_window": [cluster_z_min, cluster_z_max],
            "firstpass_permutations": firstpass_permutations,
            "firstpass_random_repeats": firstpass_random_repeats,
            "ensemble_permutations": ensemble_permutations,
            "ensemble_random_repeats": ensemble_random_repeats,
        },
        "primary_firstpass": firstpass,
        "cluster_model_consensus": ensemble_summaries,
        "full_model_run_summary": {
            "model_run_count": len(all_model_rows),
            "mass_positive_run_fraction": float(np.mean(mass_corr > 0.0)),
            "mass_significant_run_fraction": float(np.mean((mass_corr > 0.0) & (mass_p < 0.05))),
            "mass_median_partial_correlation": float(np.median(mass_corr)),
            "sfr_positive_run_fraction": float(np.mean(sfr_corr > 0.0)),
            "sfr_significant_run_fraction": float(np.mean((sfr_corr > 0.0) & (sfr_p < 0.05))),
            "sfr_median_partial_correlation": float(np.median(sfr_corr)),
            "local_kappa_positive_run_fraction": float(np.mean(local_kappa_corr > 0.0)),
            "local_kappa_significant_run_fraction": float(np.mean((local_kappa_corr > 0.0) & (local_kappa_p < 0.05))),
            "local_kappa_median_partial_correlation": float(np.median(local_kappa_corr)),
            "photometric_proxy_positive_run_fraction": float(np.mean(proxy_corr > 0.0)),
            "photometric_proxy_significant_run_fraction": float(np.mean((proxy_corr > 0.0) & (proxy_p < 0.05))),
            "photometric_proxy_median_partial_correlation": float(np.nanmedian(proxy_corr)),
            "random_null_significant_run_fraction": float(np.mean(random_p < 0.05)),
        },
    }
    (out_dir / "full_experiment.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def run_robustness_sweep(
    cluster_keys: list[str],
    out_dir: Path,
    inner_radii_arcsec: list[float],
    ring_inner_radii_arcsec: list[float],
    ring_outer_radii_arcsec: list[float],
    smooth_sigmas_px: list[float],
    residual_modes: list[str],
    radial_bin_arcsec: float,
    cluster_z_min: float,
    cluster_z_max: float,
    permutations: int,
    random_repeats: int,
    seed: int,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    run_index = 0
    for cluster_key in cluster_keys:
        for residual_mode, inner_radius_arcsec, ring_inner_arcsec, ring_outer_arcsec, smooth_sigma_px in product(
            residual_modes,
            inner_radii_arcsec,
            ring_inner_radii_arcsec,
            ring_outer_radii_arcsec,
            smooth_sigmas_px,
        ):
            if not (inner_radius_arcsec < ring_inner_arcsec < ring_outer_arcsec):
                continue
            run_index += 1
            run_label = (
                f"{cluster_key}_{residual_mode}_"
                f"i{str(inner_radius_arcsec).replace('.', 'p')}_"
                f"ri{str(ring_inner_arcsec).replace('.', 'p')}_"
                f"ro{str(ring_outer_arcsec).replace('.', 'p')}_"
                f"s{str(smooth_sigma_px).replace('.', 'p')}"
            )
            result = analyze_hff_cluster_firstpass(
                cluster_key=cluster_key,
                out_dir=out_dir / run_label,
                inner_radius_arcsec=inner_radius_arcsec,
                ring_inner_arcsec=ring_inner_arcsec,
                ring_outer_arcsec=ring_outer_arcsec,
                smooth_sigma_px=smooth_sigma_px,
                residual_mode=residual_mode,
                radial_bin_arcsec=radial_bin_arcsec,
                cluster_z_min=cluster_z_min,
                cluster_z_max=cluster_z_max,
                permutations=permutations,
                random_repeats=random_repeats,
                seed=seed + run_index,
            )
            rows.append(summarize_cluster_result(result))

    if not rows:
        raise ValueError("Robustness sweep produced no valid runs.")

    write_csv_rows(out_dir / "robustness_sweep.csv", rows, list(rows[0].keys()))

    summary = {}
    for cluster_key in cluster_keys:
        for residual_mode in residual_modes:
            subset = [row for row in rows if row["cluster_key"] == cluster_key and row["residual_mode"] == residual_mode]
            if not subset:
                continue
            mass_corr = np.asarray([row["mass_partial_correlation"] for row in subset], dtype=float)
            mass_p = np.asarray([row["mass_partial_p_value"] for row in subset], dtype=float)
            sfr_corr = np.asarray([row["sfr_partial_correlation"] for row in subset], dtype=float)
            sfr_p = np.asarray([row["sfr_partial_p_value"] for row in subset], dtype=float)
            local_kappa_corr = np.asarray([row["local_kappa_partial_correlation"] for row in subset], dtype=float)
            local_kappa_p = np.asarray([row["local_kappa_partial_p_value"] for row in subset], dtype=float)
            proxy_corr = np.asarray([row["photometric_proxy_partial_correlation"] for row in subset], dtype=float)
            proxy_p = np.asarray([row["photometric_proxy_partial_p_value"] for row in subset], dtype=float)
            random_p = np.asarray([row["random_null_p_value"] for row in subset], dtype=float)
            summary[f"{cluster_key}:{residual_mode}"] = {
                "run_count": len(subset),
                "mass_positive_fraction": float(np.mean(mass_corr > 0.0)),
                "mass_significant_fraction": float(np.mean((mass_corr > 0.0) & (mass_p < 0.05))),
                "mass_median_partial_correlation": float(np.median(mass_corr)),
                "mass_median_p_value": float(np.median(mass_p)),
                "sfr_positive_fraction": float(np.mean(sfr_corr > 0.0)),
                "sfr_significant_fraction": float(np.mean((sfr_corr > 0.0) & (sfr_p < 0.05))),
                "sfr_median_partial_correlation": float(np.median(sfr_corr)),
                "sfr_median_p_value": float(np.median(sfr_p)),
                "local_kappa_positive_fraction": float(np.mean(local_kappa_corr > 0.0)),
                "local_kappa_significant_fraction": float(np.mean((local_kappa_corr > 0.0) & (local_kappa_p < 0.05))),
                "local_kappa_median_partial_correlation": float(np.median(local_kappa_corr)),
                "local_kappa_median_p_value": float(np.median(local_kappa_p)),
                "photometric_proxy_positive_fraction": float(np.mean(proxy_corr > 0.0)),
                "photometric_proxy_significant_fraction": float(np.mean((proxy_corr > 0.0) & (proxy_p < 0.05))),
                "photometric_proxy_median_partial_correlation": float(np.nanmedian(proxy_corr)),
                "photometric_proxy_median_p_value": float(np.nanmedian(proxy_p)),
                "random_null_significant_fraction": float(np.mean(random_p < 0.05)),
            }

    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "parameter_grid": {
            "cluster_keys": cluster_keys,
            "inner_radii_arcsec": inner_radii_arcsec,
            "ring_inner_radii_arcsec": ring_inner_radii_arcsec,
            "ring_outer_radii_arcsec": ring_outer_radii_arcsec,
            "smooth_sigmas_px": smooth_sigmas_px,
            "residual_modes": residual_modes,
            "radial_bin_arcsec": radial_bin_arcsec,
            "cluster_member_z_window": [cluster_z_min, cluster_z_max],
            "permutations": permutations,
            "random_repeats": random_repeats,
        },
        "n_runs": len(rows),
        "summary": summary,
    }
    (out_dir / "robustness_sweep.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def command_manifest(args: argparse.Namespace) -> None:
    manifest = build_manifest(args.targets, args.radius_deg)
    Path(args.out).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote archive manifest to {args.out}")
    for target in manifest["targets"]:
        summary = target["collections"]
        print(f"{target['name']}: {target['record_count']} public HST/JWST image products {summary}")


def command_score(args: argparse.Namespace) -> None:
    result = score_sign_flip(
        galaxies_path=Path(args.galaxies),
        residuals_path=Path(args.residuals),
        out_path=Path(args.out),
        inner_radius_arcsec=args.inner_radius_arcsec,
        ring_inner_arcsec=args.ring_inner_arcsec,
        ring_outer_arcsec=args.ring_outer_arcsec,
    )
    print(json.dumps(result["regression"], indent=2))


def command_simulate(args: argparse.Namespace) -> None:
    result = simulate_dataset(Path(args.out_dir), seed=args.seed)
    print(json.dumps(result, indent=2))


def command_abell370(args: argparse.Namespace) -> None:
    result = analyze_hff_cluster_firstpass(
        cluster_key="abell370",
        out_dir=Path(args.out_dir),
        inner_radius_arcsec=args.inner_radius_arcsec,
        ring_inner_arcsec=args.ring_inner_arcsec,
        ring_outer_arcsec=args.ring_outer_arcsec,
        smooth_sigma_px=args.smooth_sigma_px,
        residual_mode=args.residual_mode,
        radial_bin_arcsec=args.radial_bin_arcsec,
        cluster_z_min=args.cluster_z_min,
        cluster_z_max=args.cluster_z_max,
        permutations=args.permutations,
        random_repeats=args.random_repeats,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2))


def command_rxcj2248(args: argparse.Namespace) -> None:
    result = analyze_hff_cluster_firstpass(
        cluster_key="rxcj2248",
        out_dir=Path(args.out_dir),
        inner_radius_arcsec=args.inner_radius_arcsec,
        ring_inner_arcsec=args.ring_inner_arcsec,
        ring_outer_arcsec=args.ring_outer_arcsec,
        smooth_sigma_px=args.smooth_sigma_px,
        residual_mode=args.residual_mode,
        radial_bin_arcsec=args.radial_bin_arcsec,
        cluster_z_min=args.cluster_z_min,
        cluster_z_max=args.cluster_z_max,
        permutations=args.permutations,
        random_repeats=args.random_repeats,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2))


def command_combined_firstpass(args: argparse.Namespace) -> None:
    result = analyze_combined_hff_firstpass(
        cluster_keys=args.clusters,
        out_dir=Path(args.out_dir),
        inner_radius_arcsec=args.inner_radius_arcsec,
        ring_inner_arcsec=args.ring_inner_arcsec,
        ring_outer_arcsec=args.ring_outer_arcsec,
        smooth_sigma_px=args.smooth_sigma_px,
        residual_mode=args.residual_mode,
        radial_bin_arcsec=args.radial_bin_arcsec,
        cluster_z_min=args.cluster_z_min,
        cluster_z_max=args.cluster_z_max,
        permutations=args.permutations,
        random_repeats=args.random_repeats,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2))


def command_robustness_sweep(args: argparse.Namespace) -> None:
    result = run_robustness_sweep(
        cluster_keys=args.clusters,
        out_dir=Path(args.out_dir),
        inner_radii_arcsec=args.inner_radii_arcsec,
        ring_inner_radii_arcsec=args.ring_inner_radii_arcsec,
        ring_outer_radii_arcsec=args.ring_outer_radii_arcsec,
        smooth_sigmas_px=args.smooth_sigmas_px,
        residual_modes=args.residual_modes,
        radial_bin_arcsec=args.radial_bin_arcsec,
        cluster_z_min=args.cluster_z_min,
        cluster_z_max=args.cluster_z_max,
        permutations=args.permutations,
        random_repeats=args.random_repeats,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2))


def command_full_experiment(args: argparse.Namespace) -> None:
    result = run_full_experiment(
        cluster_keys=args.clusters,
        out_dir=Path(args.out_dir),
        inner_radius_arcsec=args.inner_radius_arcsec,
        ring_inner_arcsec=args.ring_inner_arcsec,
        ring_outer_arcsec=args.ring_outer_arcsec,
        smooth_sigma_px=args.smooth_sigma_px,
        residual_mode=args.residual_mode,
        radial_bin_arcsec=args.radial_bin_arcsec,
        cluster_z_min=args.cluster_z_min,
        cluster_z_max=args.cluster_z_max,
        firstpass_permutations=args.firstpass_permutations,
        firstpass_random_repeats=args.firstpass_random_repeats,
        ensemble_permutations=args.ensemble_permutations,
        ensemble_random_repeats=args.ensemble_random_repeats,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Archive manifesting and sign-flip analysis for the metering-metric lensing proposal.",
        epilog=CLI_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    manifest = subparsers.add_parser("manifest", help="Query public HST/JWST imaging near the target clusters.")
    manifest.add_argument("--out", default="lensing_manifest.json")
    manifest.add_argument("--radius-deg", type=float, default=0.05)
    manifest.add_argument("--targets", nargs="+", default=DEFAULT_TARGETS)
    manifest.set_defaults(func=command_manifest)

    score = subparsers.add_parser("score", help="Score convergence-divergence sign flips and fit the proxy regression.")
    score.add_argument("--galaxies", required=True, help="CSV with ra, dec, metallicity, sfr, mass.")
    score.add_argument("--residuals", required=True, help="CSV with ra, dec, residual.")
    score.add_argument("--out", default="lensing_results.json")
    score.add_argument("--inner-radius-arcsec", type=float, default=2.5)
    score.add_argument("--ring-inner-arcsec", type=float, default=5.0)
    score.add_argument("--ring-outer-arcsec", type=float, default=9.0)
    score.set_defaults(func=command_score)

    simulate = subparsers.add_parser("simulate", help="Generate a toy sign-flip dataset and run the pipeline.")
    simulate.add_argument("--out-dir", default="toy_lensing")
    simulate.add_argument("--seed", type=int, default=7)
    simulate.set_defaults(func=command_simulate)

    abell370 = subparsers.add_parser(
        "abell370-firstpass",
        help="Run a first-pass sign-flip analysis using the official Abell 370 HFF catalog and kappa map.",
    )
    abell370.add_argument("--out-dir", default="abell370_firstpass")
    abell370.add_argument("--inner-radius-arcsec", type=float, default=2.5)
    abell370.add_argument("--ring-inner-arcsec", type=float, default=5.0)
    abell370.add_argument("--ring-outer-arcsec", type=float, default=9.0)
    abell370.add_argument("--smooth-sigma-px", type=float, default=18.0)
    abell370.add_argument("--residual-mode", choices=["gaussian_highpass", "radial_median"], default="radial_median")
    abell370.add_argument("--radial-bin-arcsec", type=float, default=8.0)
    abell370.add_argument("--cluster-z-min", type=float, default=0.35)
    abell370.add_argument("--cluster-z-max", type=float, default=0.45)
    abell370.add_argument("--permutations", type=int, default=500)
    abell370.add_argument("--random-repeats", type=int, default=200)
    abell370.add_argument("--seed", type=int, default=7)
    abell370.set_defaults(func=command_abell370)

    rxcj2248 = subparsers.add_parser(
        "rxcj2248-firstpass",
        help="Run the same first-pass sign-flip analysis on RXC J2248 / Abell S1063.",
    )
    rxcj2248.add_argument("--out-dir", default="rxcj2248_firstpass")
    rxcj2248.add_argument("--inner-radius-arcsec", type=float, default=2.5)
    rxcj2248.add_argument("--ring-inner-arcsec", type=float, default=5.0)
    rxcj2248.add_argument("--ring-outer-arcsec", type=float, default=9.0)
    rxcj2248.add_argument("--smooth-sigma-px", type=float, default=18.0)
    rxcj2248.add_argument("--residual-mode", choices=["gaussian_highpass", "radial_median"], default="radial_median")
    rxcj2248.add_argument("--radial-bin-arcsec", type=float, default=8.0)
    rxcj2248.add_argument("--cluster-z-min", type=float, default=0.35)
    rxcj2248.add_argument("--cluster-z-max", type=float, default=0.45)
    rxcj2248.add_argument("--permutations", type=int, default=500)
    rxcj2248.add_argument("--random-repeats", type=int, default=200)
    rxcj2248.add_argument("--seed", type=int, default=7)
    rxcj2248.set_defaults(func=command_rxcj2248)

    combined = subparsers.add_parser(
        "combined-firstpass",
        help="Run the first-pass HFF analysis on multiple clusters and combine the cluster-level evidence.",
    )
    combined.add_argument("--out-dir", default="combined_firstpass")
    combined.add_argument("--clusters", nargs="+", default=["abell370", "rxcj2248"])
    combined.add_argument("--inner-radius-arcsec", type=float, default=2.5)
    combined.add_argument("--ring-inner-arcsec", type=float, default=5.0)
    combined.add_argument("--ring-outer-arcsec", type=float, default=9.0)
    combined.add_argument("--smooth-sigma-px", type=float, default=18.0)
    combined.add_argument("--residual-mode", choices=["gaussian_highpass", "radial_median"], default="radial_median")
    combined.add_argument("--radial-bin-arcsec", type=float, default=8.0)
    combined.add_argument("--cluster-z-min", type=float, default=0.35)
    combined.add_argument("--cluster-z-max", type=float, default=0.45)
    combined.add_argument("--permutations", type=int, default=500)
    combined.add_argument("--random-repeats", type=int, default=200)
    combined.add_argument("--seed", type=int, default=7)
    combined.set_defaults(func=command_combined_firstpass)

    robustness = subparsers.add_parser(
        "robustness-sweep",
        help="Sweep aperture and residual-model settings across the HFF first-pass analysis.",
    )
    robustness.add_argument("--out-dir", default="robustness_sweep")
    robustness.add_argument("--clusters", nargs="+", default=["abell370", "rxcj2248"])
    robustness.add_argument("--inner-radii-arcsec", nargs="+", type=float, default=[2.0, 2.5])
    robustness.add_argument("--ring-inner-radii-arcsec", nargs="+", type=float, default=[4.5, 5.0])
    robustness.add_argument("--ring-outer-radii-arcsec", nargs="+", type=float, default=[8.5])
    robustness.add_argument("--smooth-sigmas-px", nargs="+", type=float, default=[14.0, 18.0])
    robustness.add_argument(
        "--residual-modes",
        nargs="+",
        choices=["gaussian_highpass", "radial_median"],
        default=["gaussian_highpass", "radial_median"],
    )
    robustness.add_argument("--radial-bin-arcsec", type=float, default=8.0)
    robustness.add_argument("--cluster-z-min", type=float, default=0.35)
    robustness.add_argument("--cluster-z-max", type=float, default=0.45)
    robustness.add_argument("--permutations", type=int, default=100)
    robustness.add_argument("--random-repeats", type=int, default=50)
    robustness.add_argument("--seed", type=int, default=7)
    robustness.set_defaults(func=command_robustness_sweep)

    full_experiment = subparsers.add_parser(
        "full-experiment",
        help="Run the full archival experiment: primary two-cluster analysis plus all available independent HFF model families.",
    )
    full_experiment.add_argument("--out-dir", default="full_experiment")
    full_experiment.add_argument("--clusters", nargs="+", default=["abell370", "rxcj2248"])
    full_experiment.add_argument("--inner-radius-arcsec", type=float, default=2.5)
    full_experiment.add_argument("--ring-inner-arcsec", type=float, default=5.0)
    full_experiment.add_argument("--ring-outer-arcsec", type=float, default=9.0)
    full_experiment.add_argument("--smooth-sigma-px", type=float, default=18.0)
    full_experiment.add_argument("--residual-mode", choices=["gaussian_highpass", "radial_median"], default="radial_median")
    full_experiment.add_argument("--radial-bin-arcsec", type=float, default=8.0)
    full_experiment.add_argument("--cluster-z-min", type=float, default=0.35)
    full_experiment.add_argument("--cluster-z-max", type=float, default=0.45)
    full_experiment.add_argument("--firstpass-permutations", type=int, default=500)
    full_experiment.add_argument("--firstpass-random-repeats", type=int, default=200)
    full_experiment.add_argument("--ensemble-permutations", type=int, default=100)
    full_experiment.add_argument("--ensemble-random-repeats", type=int, default=50)
    full_experiment.add_argument("--seed", type=int, default=7)
    full_experiment.set_defaults(func=command_full_experiment)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

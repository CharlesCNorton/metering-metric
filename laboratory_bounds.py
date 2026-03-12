"""
Laboratory-side bound workbench for the metering-metric project.

This module turns the repository's Casimir-scale benchmark claim into an
executed, measurement-anchored null bound on the physical coupling alpha.

The benchmark is explicit:

- anomaly anchor: at alpha = 1, the metering-induced pressure correction is
  taken to be 16 times the standard Casimir pressure
- weak-coupling scaling: the anomalous pressure scales as alpha^q with q = 2
  by default
- measurement anchors: use published precision Casimir-pressure uncertainties

This is a benchmark-model laboratory bound, not a final first-principles
derivation of the physical alpha.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


HBAR = 1.054_571_817e-34
C = 299_792_458.0
PI = math.pi
NM = 1.0e-9
MICRON = 1.0e-6
MPA = 1.0e-3


@dataclass(frozen=True)
class AbsolutePressureMeasurement:
    label: str
    absolute_error_pa: float
    separation_min_m: float
    separation_max_m: float
    source: str
    note: str


@dataclass(frozen=True)
class RelativePressureMeasurement:
    label: str
    relative_error: float
    separation_m: float
    source: str
    note: str


DECCA_2003_PRESSURE = AbsolutePressureMeasurement(
    label="Decca et al. dynamic parallel-plate Casimir pressure",
    absolute_error_pa=0.6 * MPA,
    separation_min_m=0.2 * MICRON,
    separation_max_m=1.2 * MICRON,
    source="Phys. Rev. D 68 (2003) 116003 / arXiv:hep-ph/0310157",
    note="Published absolute pressure uncertainty of about 0.6 mPa over 0.2-1.2 um.",
)

DECCA_2007_RELATIVE = RelativePressureMeasurement(
    label="Decca et al. gold-plate Casimir pressure minimum relative error",
    relative_error=0.0019,
    separation_m=160.0 * NM,
    source="Phys. Rev. D 75 (2007) 077101 / arXiv:hep-ph/0703290",
    note="Published minimum relative experimental error of about 0.19% in the 160-750 nm range.",
)

DEFAULT_ALPHA_TARGETS = (1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4)


def perfect_conductor_casimir_pressure(separation_m: float) -> float:
    if separation_m <= 0.0:
        raise ValueError("separation_m must be positive.")
    return (PI * PI * HBAR * C) / (240.0 * separation_m**4)


def benchmark_anomalous_pressure(
    alpha: float,
    separation_m: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> float:
    if alpha < 0.0:
        raise ValueError("alpha must be nonnegative.")
    if anomaly_factor_at_alpha_one <= 0.0:
        raise ValueError("anomaly_factor_at_alpha_one must be positive.")
    if scaling_power <= 0.0:
        raise ValueError("scaling_power must be positive.")
    casimir_pressure = perfect_conductor_casimir_pressure(separation_m)
    return anomaly_factor_at_alpha_one * casimir_pressure * alpha**scaling_power


def alpha_upper_bound_from_null(
    absolute_error_pa: float,
    separation_m: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> float:
    if absolute_error_pa <= 0.0:
        raise ValueError("absolute_error_pa must be positive.")
    reference_pressure = anomaly_factor_at_alpha_one * perfect_conductor_casimir_pressure(separation_m)
    return (absolute_error_pa / reference_pressure) ** (1.0 / scaling_power)


def alpha_upper_bound_from_relative_null(
    relative_error: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> float:
    if relative_error <= 0.0:
        raise ValueError("relative_error must be positive.")
    if anomaly_factor_at_alpha_one <= 0.0:
        raise ValueError("anomaly_factor_at_alpha_one must be positive.")
    if scaling_power <= 0.0:
        raise ValueError("scaling_power must be positive.")
    return (relative_error / anomaly_factor_at_alpha_one) ** (1.0 / scaling_power)


def build_absolute_scan_report(
    measurement: AbsolutePressureMeasurement,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    samples: int,
) -> dict:
    if samples < 2:
        raise ValueError("samples must be at least 2.")
    separations = np.geomspace(measurement.separation_min_m, measurement.separation_max_m, samples)
    rows = []
    for separation_m in separations:
        casimir_pressure = perfect_conductor_casimir_pressure(float(separation_m))
        alpha_upper = alpha_upper_bound_from_null(
            absolute_error_pa=measurement.absolute_error_pa,
            separation_m=float(separation_m),
            anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
            scaling_power=scaling_power,
        )
        rows.append(
            {
                "separation_m": float(separation_m),
                "separation_um": float(separation_m / MICRON),
                "standard_casimir_pressure_pa": float(casimir_pressure),
                "absolute_measurement_error_pa": float(measurement.absolute_error_pa),
                "relative_error_vs_standard_casimir": float(measurement.absolute_error_pa / casimir_pressure),
                "alpha_upper_bound": float(alpha_upper),
            }
        )

    best_row = min(rows, key=lambda row: row["alpha_upper_bound"])
    weakest_row = max(rows, key=lambda row: row["alpha_upper_bound"])
    reference_alpha = 1.0
    reference_rows = [
        {
            **row,
            "benchmark_anomalous_pressure_at_alpha_one_pa": float(
                benchmark_anomalous_pressure(
                    alpha=reference_alpha,
                    separation_m=row["separation_m"],
                    anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
                    scaling_power=scaling_power,
                )
            ),
        }
        for row in rows
    ]
    return {
        "type": "absolute_pressure_scan",
        "measurement": {
            "label": measurement.label,
            "source": measurement.source,
            "note": measurement.note,
            "absolute_error_pa": float(measurement.absolute_error_pa),
            "separation_min_um": float(measurement.separation_min_m / MICRON),
            "separation_max_um": float(measurement.separation_max_m / MICRON),
        },
        "benchmark_model": {
            "anomaly_factor_at_alpha_one": float(anomaly_factor_at_alpha_one),
            "scaling_power": float(scaling_power),
            "benchmark_statement": (
                "At alpha = 1, the anomalous metering-induced pressure is taken to be "
                f"{anomaly_factor_at_alpha_one:.6g} times the standard Casimir pressure."
            ),
        },
        "best_bound": best_row,
        "weakest_bound": weakest_row,
        "rows": reference_rows,
    }


def build_relative_point_report(
    measurement: RelativePressureMeasurement,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> dict:
    standard_pressure = perfect_conductor_casimir_pressure(measurement.separation_m)
    alpha_upper = alpha_upper_bound_from_relative_null(
        relative_error=measurement.relative_error,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
    )
    return {
        "type": "relative_pressure_point",
        "measurement": {
            "label": measurement.label,
            "source": measurement.source,
            "note": measurement.note,
            "relative_error": float(measurement.relative_error),
            "separation_m": float(measurement.separation_m),
            "separation_um": float(measurement.separation_m / MICRON),
        },
        "derived": {
            "standard_casimir_pressure_pa": float(standard_pressure),
            "absolute_error_pa": float(measurement.relative_error * standard_pressure),
            "alpha_upper_bound": float(alpha_upper),
            "benchmark_anomalous_pressure_at_alpha_one_pa": float(
                benchmark_anomalous_pressure(
                    alpha=1.0,
                    separation_m=measurement.separation_m,
                    anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
                    scaling_power=scaling_power,
                )
            ),
        },
    }


def build_alpha_target_table(
    reference_separation_m: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    alpha_targets: Iterable[float],
    comparison_absolute_error_pa: float,
) -> list[dict]:
    standard_pressure = perfect_conductor_casimir_pressure(reference_separation_m)
    rows = []
    for alpha_target in alpha_targets:
        required_relative = anomaly_factor_at_alpha_one * alpha_target**scaling_power
        required_absolute = required_relative * standard_pressure
        rows.append(
            {
                "alpha_target": float(alpha_target),
                "reference_separation_m": float(reference_separation_m),
                "reference_separation_um": float(reference_separation_m / MICRON),
                "standard_casimir_pressure_pa": float(standard_pressure),
                "required_relative_anomaly": float(required_relative),
                "required_absolute_anomaly_pa": float(required_absolute),
                "required_absolute_anomaly_mpa": float(required_absolute / MPA),
                "decca_2003_error_ratio": float(required_absolute / comparison_absolute_error_pa),
            }
        )
    return rows


def build_casimir_benchmark_report(
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    samples: int,
    reference_separation_m: float,
    alpha_targets: Iterable[float],
) -> dict:
    absolute_scan = build_absolute_scan_report(
        measurement=DECCA_2003_PRESSURE,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
        samples=samples,
    )
    relative_point = build_relative_point_report(
        measurement=DECCA_2007_RELATIVE,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
    )
    strongest_anchor = min(
        [
            {
                "label": absolute_scan["measurement"]["label"],
                "source": absolute_scan["measurement"]["source"],
                "alpha_upper_bound": absolute_scan["best_bound"]["alpha_upper_bound"],
                "anchor_type": absolute_scan["type"],
                "reference_separation_um": absolute_scan["best_bound"]["separation_um"],
            },
            {
                "label": relative_point["measurement"]["label"],
                "source": relative_point["measurement"]["source"],
                "alpha_upper_bound": relative_point["derived"]["alpha_upper_bound"],
                "anchor_type": relative_point["type"],
                "reference_separation_um": relative_point["measurement"]["separation_um"],
            },
        ],
        key=lambda row: row["alpha_upper_bound"],
    )
    alpha_target_table = build_alpha_target_table(
        reference_separation_m=reference_separation_m,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
        alpha_targets=alpha_targets,
        comparison_absolute_error_pa=DECCA_2003_PRESSURE.absolute_error_pa,
    )
    return {
        "benchmark_model": {
            "anomaly_factor_at_alpha_one": float(anomaly_factor_at_alpha_one),
            "scaling_power": float(scaling_power),
            "statement": (
                "At alpha = 1, the anomalous metering-induced pressure is benchmarked as "
                f"{anomaly_factor_at_alpha_one:.6g} times the standard Casimir pressure, "
                f"with weak-coupling scaling proportional to alpha^{scaling_power:.6g}."
            ),
        },
        "measurement_anchors": {
            "decca_2003_absolute_scan": absolute_scan,
            "decca_2007_relative_point": relative_point,
        },
        "strongest_bound": strongest_anchor,
        "reference_target_scan": {
            "reference_separation_m": float(reference_separation_m),
            "reference_separation_um": float(reference_separation_m / MICRON),
            "alpha_target_table": alpha_target_table,
        },
        "executive_summary": {
            "benchmark_lab_bound_alpha_upper": float(strongest_anchor["alpha_upper_bound"]),
            "recommended_reading": (
                "This is a benchmark-model null bound. It is stronger than the older internal "
                "WEC consistency cutoff, but it still rests on the repository's present "
                "Casimir-scale anomaly model rather than a completed first-principles lab derivation."
            ),
        },
    }


def write_json_report(report: dict, out_path: str | None) -> None:
    rendered = json.dumps(report, indent=2)
    if out_path is None:
        print(rendered)
        return
    destination = Path(out_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


def command_casimir_benchmark_report(args: argparse.Namespace) -> None:
    report = build_casimir_benchmark_report(
        anomaly_factor_at_alpha_one=args.anomaly_factor_at_alpha_one,
        scaling_power=args.scaling_power,
        samples=args.samples,
        reference_separation_m=args.reference_separation_um * MICRON,
        alpha_targets=args.alpha_targets,
    )
    write_json_report(report, args.out)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Laboratory-side benchmark bounds for the metering-metric coupling.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    decca = subparsers.add_parser(
        "casimir-benchmark-report",
        help="Build a measurement-anchored Casimir benchmark report for alpha.",
    )
    decca.add_argument("--anomaly-factor-at-alpha-one", type=float, default=16.0)
    decca.add_argument("--scaling-power", type=float, default=2.0)
    decca.add_argument("--samples", type=int, default=64)
    decca.add_argument("--reference-separation-um", type=float, default=0.2)
    decca.add_argument("--alpha-targets", type=float, nargs="+", default=list(DEFAULT_ALPHA_TARGETS))
    decca.add_argument("--out", type=str, default=None)
    decca.set_defaults(func=command_casimir_benchmark_report)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

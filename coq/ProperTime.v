(** * ProperTime: Theorem 3 -- Duration requires a meter *)
(**
   Theorem 3a: For a static observer at position x where mu(x) = 0,
   the proper time accumulated over any coordinate time interval is zero.

   Theorem 3b: For a worldline passing through a region where mu = 0,
   the proper time contribution from that region is zero.

   These are the formal statements of the core claim: regions without
   metering subsystems accumulate zero proper time.
*)

Require Import Reals.
Require Import Lra.
Require Import Coquelicot.Coquelicot.
Require Import MeteringMetric.Definitions.
Open Scope R_scope.

(** ** Theorem 3a: Static observer in meterless region *)

Theorem proper_time_vanishes_static :
  forall (f : R -> R) (x T : R),
  f 0 = 0 ->
  forall (mu : R -> R), mu x = 0 ->
  proper_time_static f mu x T = 0.
Proof.
  intros f x T Hf0 mu Hmu0.
  unfold proper_time_static, lapse.
  rewrite Hmu0.
  rewrite Hf0.
  ring.
Qed.

(** ** Theorem 3b: Worldline through meterless region *)
(**
   If a worldline gamma passes through a region where mu(gamma(t)) = 0
   for all t in [t1, t2], then the proper time accumulated over that
   interval is zero.

   The proof uses the fact that the integrand f(mu(gamma(t))) is
   identically zero on the interval, so the Riemann integral is zero.
*)

Theorem proper_time_zero_interval :
  forall (f : R -> R) (mu : R -> R) (gamma : R -> R) (t1 t2 : R),
  f 0 = 0 ->
  t1 <= t2 ->
  (forall t, t1 <= t <= t2 -> mu (gamma t) = 0) ->
  RInt (fun t => lapse f mu (gamma t)) t1 t2 = 0.
Proof.
  intros f mu gamma t1 t2 Hf0 Hle Hzero.
  (* Replace integrand with the zero function *)
  transitivity (RInt (fun _ : R => 0) t1 t2).
  - apply RInt_ext.
    intros x Hx.
    unfold lapse.
    (* x is in [Rmin t1 t2, Rmax t1 t2] = [t1, t2] since t1 <= t2 *)
    assert (Hmin : Rmin t1 t2 = t1) by (apply Rmin_left; exact Hle).
    assert (Hmax : Rmax t1 t2 = t2) by (apply Rmax_right; exact Hle).
    rewrite Hmin, Hmax in Hx.
    rewrite Hzero.
    + exact Hf0.
    + split; left; apply Hx.
  - (* Integral of zero is zero *)
    rewrite RInt_const.
    unfold scal; simpl.
    unfold mult; simpl.
    lra.
Qed.

(** ** Corollary: proper time is bounded by coordinate time *)
(**
   Since N(x) <= 1 for all x, proper time never exceeds coordinate time.
   This is weaker than the vanishing theorem but useful as a sanity check.
*)

Lemma proper_time_bounded :
  forall (f : R -> R) (mu : R -> R) (x T : R),
  ActivationProps f -> MeteringProps mu ->
  0 <= T ->
  proper_time_static f mu x T <= T.
Proof.
  intros f mu x T Hf Hmu HT.
  unfold proper_time_static, lapse.
  assert (HN : f (mu x) <= 1) by (apply act_le_one; exact Hf).
  assert (HNge : 0 <= f (mu x)).
  { apply act_nonneg.
    - exact Hf.
    - left. apply met_pos. exact Hmu. }
  nra.
Qed.

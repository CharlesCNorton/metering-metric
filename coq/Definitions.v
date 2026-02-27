(** * Definitions: Core concepts of the metering metric framework *)
(**
   This file defines the mathematical objects underlying the conditional
   temporal metric: metering density, activation functions, the lapse
   function, tortoise coordinate, effective potential, and proper time.

   All definitions are parametric in the activation function f and
   metering density mu. The theorems in subsequent files prove structural
   properties that hold for ANY monotone activation with f(0)=0 and
   ANY strictly positive decaying metering density.
*)

Require Import Reals.
Require Import Lra.
Require Import Coquelicot.Coquelicot.
Open Scope R_scope.

(** ** Activation function hypotheses *)
(**
   An activation function f : R -> R satisfies:
   - f(0) = 0
   - f is nonneg on nonneg inputs
   - f is bounded above by 1
   - f is monotone increasing on nonneg inputs
   - f is continuous everywhere
   - f is strictly positive on strictly positive inputs

   Examples: tanh, sigmoid shifted to f(0)=0, arctan/pi, erf.
*)

Record ActivationProps (f : R -> R) : Prop := mkActivation {
  act_zero : f 0 = 0;
  act_nonneg : forall z, 0 <= z -> 0 <= f z;
  act_le_one : forall z, f z <= 1;
  act_mono : forall z1 z2, 0 <= z1 -> z1 <= z2 -> f z1 <= f z2;
  act_cont : forall z, continuous f z;
  act_pos : forall z, 0 < z -> 0 < f z;
  (** Epsilon-delta continuity at 0 for nonneg inputs.
      This follows from [act_cont] and [act_zero] and [act_nonneg],
      but is provided as an explicit hypothesis to avoid entanglement
      with the Coquelicot filter API in downstream proofs. *)
  act_cont_0 : forall eps, eps > 0 ->
    exists delta, delta > 0 /\ forall z, 0 <= z < delta -> f z < eps
}.

(** ** Metering density hypotheses *)
(**
   A metering density mu : R -> R satisfies:
   - mu(x) > 0 for all x (strictly positive, "rapidly decreasing" case)
   - mu is continuous
   - mu decays to 0 as x -> +infinity

   The strict positivity ensures N(x) = f(mu(x)) > 0 everywhere,
   so 1/N(x) is well-defined and the tortoise integral is a proper
   Riemann integral on any bounded interval. The "horizon" is that
   this integral diverges, not that N literally reaches zero at any
   finite point.
*)

Record MeteringProps (mu : R -> R) : Prop := mkMetering {
  met_pos : forall x, 0 < mu x;
  met_cont : forall x, continuous mu x;
  met_decay : forall eps, eps > 0 ->
    exists X, forall x, X < x -> mu x < eps
}.

(** ** Lapse function *)

Definition lapse (f : R -> R) (mu : R -> R) (x : R) : R :=
  f (mu x).

(** ** Proper time for a static observer *)

Definition proper_time_static (f : R -> R) (mu : R -> R)
  (x : R) (T : R) : R :=
  lapse f mu x * T.

(** ** Tortoise coordinate *)
(**
   The tortoise coordinate is the integral of 1/N(x) from 0 to L.
   Theorem 1 (Horizon) shows this diverges as L -> infinity.
*)

Definition tortoise (N : R -> R) (L : R) : R :=
  RInt (fun x => / N x) 0 L.

(** ** Effective potential for temporal modes *)

Definition V_eff (N : R -> R) (x : R) : R :=
  / (N x * N x).

(** ** Derived properties *)

Lemma lapse_pos : forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 < lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  unfold lapse.
  apply act_pos.
  - exact Hf.
  - apply met_pos. exact Hmu.
Qed.

Lemma lapse_le_one : forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, lapse f mu x <= 1.
Proof.
  intros f mu Hf Hmu x.
  unfold lapse.
  apply act_le_one. exact Hf.
Qed.

Lemma lapse_nonneg : forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 <= lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  left. apply lapse_pos; assumption.
Qed.

Lemma inv_lapse_pos : forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 < / lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  apply Rinv_0_lt_compat.
  apply lapse_pos; assumption.
Qed.

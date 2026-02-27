(** * Auxiliary: Supporting lemmas for the Horizon and Confinement theorems *)

Require Import Reals.
Require Import Lra.
Require Import Coquelicot.Coquelicot.
Require Import MeteringMetric.Definitions.
Open Scope R_scope.

(** ** Lapse decay: if mu decays, then N = f(mu) decays *)

Lemma lapse_decays :
  forall (f : R -> R) (mu : R -> R),
  ActivationProps f -> MeteringProps mu ->
  forall eps, eps > 0 ->
  exists X, forall x, X < x -> lapse f mu x < eps.
Proof.
  intros f mu Hf Hmu eps Heps.
  unfold lapse.
  (* From act_cont_0: exists delta > 0 such that 0 <= z < delta -> f z < eps *)
  destruct (act_cont_0 f Hf eps Heps) as [delta [Hdelta Hfd]].
  (* From met_decay: exists X such that x > X -> mu(x) < delta *)
  destruct (met_decay mu Hmu delta Hdelta) as [X HX].
  exists X.
  intros x Hx.
  apply Hfd.
  split.
  - left. apply met_pos. exact Hmu.
  - apply HX. exact Hx.
Qed.

(** ** Inverse lapse grows without bound *)

Lemma inv_lapse_grows :
  forall (f : R -> R) (mu : R -> R),
  ActivationProps f -> MeteringProps mu ->
  forall M, M > 0 ->
  exists X, forall x, X < x -> M < / lapse f mu x.
Proof.
  intros f mu Hf Hmu M HM.
  assert (Heps : / M > 0) by (apply Rinv_0_lt_compat; exact HM).
  destruct (lapse_decays f mu Hf Hmu (/ M) Heps) as [X HX].
  exists X.
  intros x Hx.
  assert (HNx : lapse f mu x < / M) by (apply HX; exact Hx).
  assert (HNpos : 0 < lapse f mu x) by (apply lapse_pos; assumption).
  rewrite <- (Rinv_inv M).
  apply Rinv_lt_contravar.
  - apply Rmult_lt_0_compat; [exact HNpos | exact Heps].
  - exact HNx.
Qed.

(** ** Continuity of the inverse lapse *)

Lemma continuous_inv_lapse :
  forall (f : R -> R) (mu : R -> R),
  ActivationProps f -> MeteringProps mu ->
  forall x, continuous (fun y => / lapse f mu y) x.
Proof.
  intros f mu Hf Hmu x.
  unfold lapse.
  apply continuous_comp.
  - apply continuous_comp.
    + apply met_cont. exact Hmu.
    + apply act_cont. exact Hf.
  - apply continuous_Rinv.
    assert (H : 0 < f (mu x)).
    { apply act_pos; [exact Hf | apply met_pos; exact Hmu]. }
    lra.
Qed.

(** ** Inverse lapse is nonneg *)

Lemma inv_lapse_nonneg :
  forall (f : R -> R) (mu : R -> R),
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 <= / lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  left. apply inv_lapse_pos; assumption.
Qed.

(** ** Integrability of the inverse lapse on bounded intervals *)

Lemma ex_RInt_inv_lapse :
  forall (f : R -> R) (mu : R -> R),
  ActivationProps f -> MeteringProps mu ->
  forall a b, ex_RInt (fun x => / lapse f mu x) a b.
Proof.
  intros f mu Hf Hmu a b.
  apply (@ex_RInt_continuous R_CompleteNormedModule).
  intros x _.
  apply continuous_inv_lapse; assumption.
Qed.

(** ** Lower bound on integral of a function bounded below by a constant *)

Lemma RInt_ge_const :
  forall (g : R -> R) (a b c : R),
  a <= b ->
  (forall x, a <= x <= b -> continuous g x) ->
  (forall x, a <= x <= b -> c <= g x) ->
  c * (b - a) <= RInt g a b.
Proof.
  intros g a b c Hab Hcont Hge.
  (* Rmin/Rmax simplification since a <= b *)
  assert (Hmin : Rmin a b = a) by (apply Rmin_left; exact Hab).
  assert (Hmax : Rmax a b = b) by (apply Rmax_right; exact Hab).
  (* Step 1: RInt (fun _ => c) a b <= RInt g a b *)
  assert (Hle : RInt (fun _ : R => c) a b <= RInt g a b).
  { apply RInt_le.
    - exact Hab.
    - apply ex_RInt_const.
    - apply (@ex_RInt_continuous R_CompleteNormedModule).
      intros z Hz. rewrite Hmin, Hmax in Hz.
      apply Hcont. exact Hz.
    - intros x Hx.
      apply Hge. lra. }
  (* Step 2: RInt (fun _ => c) a b = c * (b - a) *)
  assert (Hconst : RInt (fun _ : R => c) a b = scal (b - a) c).
  { apply (@RInt_const R_CompleteNormedModule). }
  (* Step 3: scal (b-a) c = c * (b-a) for R *)
  unfold scal in Hconst; simpl in Hconst.
  unfold mult in Hconst; simpl in Hconst.
  lra.
Qed.

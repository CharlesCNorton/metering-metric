(** * meteringmetric *)
(**
   Abstract theorem package for the metering metric.

   Hypotheses:
   - [f : R -> R] is an activation function with [f 0 = 0], continuity,
     monotonicity on nonnegative inputs, boundedness by [1], and strict
     positivity on positive inputs.
   - [mu : R -> R] is a strictly positive continuous metering density
     that decays to [0] as [x -> +infinity].

   Main results:
   1. Proper time vanishes on meterless regions.
   2. The tortoise coordinate is unbounded.
   3. The effective potential diverges, and the corresponding WKB barrier
      integral can be made arbitrarily large.
*)

Require Import Reals.
Require Import Lra.
Require Import Coquelicot.Coquelicot.
Open Scope R_scope.

(** ** Activation function hypotheses *)

Record ActivationProps (f : R -> R) : Prop := mkActivation {
  act_zero : f 0 = 0;
  act_nonneg : forall z, 0 <= z -> 0 <= f z;
  act_le_one : forall z, f z <= 1;
  act_mono : forall z1 z2, 0 <= z1 -> z1 <= z2 -> f z1 <= f z2;
  act_cont : forall z, continuous f z;
  act_pos : forall z, 0 < z -> 0 < f z;
  act_cont_0 : forall eps, eps > 0 ->
    exists delta, delta > 0 /\ forall z, 0 <= z < delta -> f z < eps
}.

(** ** Metering density hypotheses *)

Record MeteringProps (mu : R -> R) : Prop := mkMetering {
  met_pos : forall x, 0 < mu x;
  met_cont : forall x, continuous mu x;
  met_decay : forall eps, eps > 0 ->
    exists X, forall x, X < x -> mu x < eps
}.

(** ** Core definitions *)

Definition lapse (f : R -> R) (mu : R -> R) (x : R) : R :=
  f (mu x).

Definition proper_time_static (f : R -> R) (mu : R -> R)
  (x : R) (T : R) : R :=
  lapse f mu x * T.

Definition proper_time_density_path (f : R -> R) (mu : R -> R)
  (gamma v : R -> R) (t : R) : R :=
  sqrt ((lapse f mu (gamma t))^2 - (v t)^2).

Definition tortoise (N : R -> R) (L : R) : R :=
  RInt (fun x => / N x) 0 L.

Definition V_eff (N : R -> R) (x : R) : R :=
  / (N x * N x).

(** ** Basic positivity and continuity facts *)

Lemma lapse_pos :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 < lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  unfold lapse.
  apply act_pos.
  - exact Hf.
  - apply met_pos. exact Hmu.
Qed.

Lemma lapse_le_one :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, lapse f mu x <= 1.
Proof.
  intros f mu Hf Hmu x.
  unfold lapse.
  apply act_le_one. exact Hf.
Qed.

Lemma lapse_nonneg :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 <= lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  left. apply lapse_pos; assumption.
Qed.

Lemma inv_lapse_pos :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 < / lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  apply Rinv_0_lt_compat.
  apply lapse_pos; assumption.
Qed.

Lemma continuous_lapse :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, continuous (lapse f mu) x.
Proof.
  intros f mu Hf Hmu x.
  unfold lapse.
  eapply continuous_comp.
  - apply met_cont. exact Hmu.
  - apply act_cont. exact Hf.
Qed.

(** ** Decay and inverse-lapse control *)

Lemma lapse_decays :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall eps, eps > 0 ->
  exists X, forall x, X < x -> lapse f mu x < eps.
Proof.
  intros f mu Hf Hmu eps Heps.
  unfold lapse.
  destruct (act_cont_0 f Hf eps Heps) as [delta [Hdelta Hfd]].
  destruct (met_decay mu Hmu delta Hdelta) as [X HX].
  exists X.
  intros x Hx.
  apply Hfd.
  split.
  - left. apply met_pos. exact Hmu.
  - apply HX. exact Hx.
Qed.

Lemma inv_lapse_grows :
  forall f mu,
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

Lemma continuous_inv_lapse :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, continuous (fun y => / lapse f mu y) x.
Proof.
  intros f mu Hf Hmu x.
  eapply continuous_comp.
  - apply continuous_lapse; assumption.
  - apply continuous_Rinv.
    assert (H : 0 < lapse f mu x) by (apply lapse_pos; assumption).
    lra.
Qed.

Lemma inv_lapse_nonneg :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, 0 <= / lapse f mu x.
Proof.
  intros f mu Hf Hmu x.
  left. apply inv_lapse_pos; assumption.
Qed.

Lemma ex_RInt_inv_lapse :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall a b, ex_RInt (fun x => / lapse f mu x) a b.
Proof.
  intros f mu Hf Hmu a b.
  apply (@ex_RInt_continuous R_CompleteNormedModule).
  intros x _.
  apply continuous_inv_lapse; assumption.
Qed.

Lemma RInt_ge_const :
  forall (g : R -> R) (a b c : R),
  a <= b ->
  (forall x, a <= x <= b -> continuous g x) ->
  (forall x, a <= x <= b -> c <= g x) ->
  c * (b - a) <= RInt g a b.
Proof.
  intros g a b c Hab Hcont Hge.
  assert (Hmin : Rmin a b = a) by (apply Rmin_left; exact Hab).
  assert (Hmax : Rmax a b = b) by (apply Rmax_right; exact Hab).
  assert (Hle : RInt (fun _ : R => c) a b <= RInt g a b).
  { apply RInt_le.
    - exact Hab.
    - apply ex_RInt_const.
    - apply (@ex_RInt_continuous R_CompleteNormedModule).
      intros z Hz. rewrite Hmin, Hmax in Hz.
      apply Hcont. exact Hz.
    - intros x Hx.
      apply Hge. lra. }
  assert (Hconst : RInt (fun _ : R => c) a b = scal (b - a) c).
  { apply (@RInt_const R_CompleteNormedModule). }
  unfold scal in Hconst; simpl in Hconst.
  unfold mult in Hconst; simpl in Hconst.
  lra.
Qed.

(** ** Proper-time theorem family *)

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

Theorem proper_time_zero_interval :
  forall (f : R -> R) (mu : R -> R) (gamma : R -> R) (t1 t2 : R),
  f 0 = 0 ->
  t1 <= t2 ->
  (forall t, t1 <= t <= t2 -> mu (gamma t) = 0) ->
  RInt (fun t => lapse f mu (gamma t)) t1 t2 = 0.
Proof.
  intros f mu gamma t1 t2 Hf0 Hle Hzero.
  transitivity (RInt (fun _ : R => 0) t1 t2).
  - apply RInt_ext.
    intros x Hx.
    unfold lapse.
    assert (Hmin : Rmin t1 t2 = t1) by (apply Rmin_left; exact Hle).
    assert (Hmax : Rmax t1 t2 = t2) by (apply Rmax_right; exact Hle).
    rewrite Hmin, Hmax in Hx.
    rewrite Hzero.
    + exact Hf0.
    + lra.
  - rewrite RInt_const.
    unfold scal; simpl.
    unfold mult; simpl.
    lra.
Qed.

Corollary proper_time_zero_interval_prescribed_path :
  forall (f : R -> R) (mu : R -> R) (gamma : R -> R) (t1 t2 : R),
  f 0 = 0 ->
  t1 <= t2 ->
  (forall t, t1 <= t <= t2 -> mu (gamma t) = 0) ->
  RInt (fun t => lapse f mu (gamma t)) t1 t2 = 0.
Proof.
  apply proper_time_zero_interval.
Qed.

Theorem proper_time_zero_interval_admissible_path :
  forall (f : R -> R) (mu gamma v : R -> R) (t1 t2 : R),
  f 0 = 0 ->
  t1 <= t2 ->
  (forall t, t1 <= t <= t2 -> mu (gamma t) = 0) ->
  (forall t, t1 <= t <= t2 -> 0 <= v t <= lapse f mu (gamma t)) ->
  RInt (fun t => proper_time_density_path f mu gamma v t) t1 t2 = 0.
Proof.
  intros f mu gamma v t1 t2 Hf0 Hle Hzero Hadm.
  transitivity (RInt (fun _ : R => 0) t1 t2).
  - apply RInt_ext.
    intros x Hx.
    assert (Hmin : Rmin t1 t2 = t1) by (apply Rmin_left; exact Hle).
    assert (Hmax : Rmax t1 t2 = t2) by (apply Rmax_right; exact Hle).
    rewrite Hmin, Hmax in Hx.
    assert (Hmu0 : mu (gamma x) = 0).
    { apply Hzero. lra. }
    assert (Hlap0 : lapse f mu (gamma x) = 0).
    { unfold lapse. rewrite Hmu0, Hf0. reflexivity. }
    assert (Hv0 : v x = 0).
    {
      assert (Hvx : 0 <= v x <= lapse f mu (gamma x)).
      { apply Hadm. lra. }
      rewrite Hlap0 in Hvx.
      lra.
    }
    unfold proper_time_density_path.
    rewrite Hlap0.
    rewrite Hv0.
    replace (0 ^ 2 - 0 ^ 2) with 0 by ring.
    rewrite sqrt_0.
    reflexivity.
  - rewrite RInt_const.
    unfold scal; simpl.
    unfold mult; simpl.
    lra.
Qed.

(**
  Scope note: the proper-time theorem family here treats static observers and
  admissible prescribed paths supported entirely in meterless regions. It
  does not yet formalize dynamically solved general worldlines.
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

(** ** Horizon theorem family *)

Theorem horizon_existence :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall B, B > 0 ->
  exists L, 0 < L /\ B < RInt (fun x => / lapse f mu x) 0 L.
Proof.
  intros f mu Hf Hmu B HB.
  assert (HBp : B + 1 > 0) by lra.
  destruct (inv_lapse_grows f mu Hf Hmu (B + 1) HBp) as [X HX].
  set (L := Rmax (X + 2) 1).
  exists L.
  split.
  - unfold L. apply Rlt_le_trans with 1; [lra | apply Rmax_r].
  - assert (Hsub_raw : (B + 1) * ((X + 2) - (X + 1)) <=
      RInt (fun x => / lapse f mu x) (X + 1) (X + 2)).
    { apply RInt_ge_const.
      - lra.
      - intros x _. apply continuous_inv_lapse; assumption.
      - intros x Hx.
        assert (Hxgt : X < x) by lra.
        left. apply HX. exact Hxgt. }
    assert (Hsub : B + 1 <= RInt (fun x => / lapse f mu x) (X + 1) (X + 2))
      by lra.
    assert (HX2_le_L : X + 2 <= L).
    { unfold L. apply Rmax_l. }
    assert (Hpos : forall x, 0 <= / lapse f mu x).
    { intro x. left. apply inv_lapse_pos; assumption. }
    assert (HintXX : ex_RInt (fun x => / lapse f mu x) (X + 1) (X + 2)).
    { apply ex_RInt_inv_lapse; assumption. }
    assert (HintX2L : ex_RInt (fun x => / lapse f mu x) (X + 2) L).
    { apply ex_RInt_inv_lapse; assumption. }
    assert (HintX1L : ex_RInt (fun x => / lapse f mu x) (X + 1) L).
    { apply ex_RInt_inv_lapse; assumption. }
    assert (Htail : (B + 1) * 1 <= RInt (fun x => / lapse f mu x) (X + 1) L).
    { rewrite <- (RInt_Chasles _ (X + 1) (X + 2) L HintXX HintX2L).
      unfold plus; simpl.
      assert (H0part : 0 <= RInt (fun x => / lapse f mu x) (X + 2) L).
      { apply RInt_ge_0.
        - exact HX2_le_L.
        - exact HintX2L.
        - intros x _. apply Hpos. }
      lra. }
    destruct (Rle_dec 0 (X + 1)) as [H0X1 | H0X1].
    + assert (Hint0X1 : ex_RInt (fun x => / lapse f mu x) 0 (X + 1)).
      { apply ex_RInt_inv_lapse; assumption. }
      rewrite <- (RInt_Chasles _ 0 (X + 1) L Hint0X1 HintX1L).
      unfold plus; simpl.
      assert (H0first : 0 <= RInt (fun x => / lapse f mu x) 0 (X + 1)).
      { apply RInt_ge_0.
        - exact H0X1.
        - exact Hint0X1.
        - intros x _. apply Hpos. }
      lra.
    + assert (Hint01 : ex_RInt (fun x => / lapse f mu x) 0 1).
      { apply ex_RInt_inv_lapse; assumption. }
      assert (Hint1L : ex_RInt (fun x => / lapse f mu x) 1 L).
      { apply ex_RInt_inv_lapse; assumption. }
      assert (H1L : 1 <= L).
      { unfold L. apply Rmax_r. }
      rewrite <- (RInt_Chasles _ 0 1 L Hint01 Hint1L).
      unfold plus; simpl.
      assert (Hge1L : 0 <= RInt (fun x => / lapse f mu x) 1 L).
      { apply RInt_ge_0.
        - exact H1L.
        - exact Hint1L.
        - intros x _. apply Hpos. }
      assert (Hge01_raw : (B + 1) * (1 - 0) <= RInt (fun x => / lapse f mu x) 0 1).
      { apply RInt_ge_const.
        - lra.
        - intros x _. apply continuous_inv_lapse; assumption.
        - intros x Hx.
          apply Rlt_le.
          apply HX.
          lra. }
      lra.
Qed.

Corollary tortoise_unbounded :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall B, B > 0 ->
  exists L, B < tortoise (lapse f mu) L.
Proof.
  intros f mu Hf Hmu B HB.
  destruct (horizon_existence f mu Hf Hmu B HB) as [L [_ HL]].
  exists L.
  unfold tortoise.
  exact HL.
Qed.

(** ** Confinement theorem family *)

Lemma continuous_V_eff :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall x, continuous (fun y => V_eff (lapse f mu) y) x.
Proof.
  intros f mu Hf Hmu x.
  unfold V_eff.
  eapply continuous_comp.
  - apply (continuous_mult (fun y => lapse f mu y) (fun y => lapse f mu y) x).
    + apply continuous_lapse; assumption.
    + apply continuous_lapse; assumption.
  - apply continuous_Rinv.
    assert (Hpos : 0 < lapse f mu x) by (apply lapse_pos; assumption).
    nra.
Qed.

Lemma wkb_integrand_continuous :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall E x,
  continuous (fun y => sqrt (V_eff (lapse f mu) y - E)) x.
Proof.
  intros f mu Hf Hmu E x.
  eapply continuous_comp.
  - apply (continuous_minus
      (fun y => V_eff (lapse f mu) y)
      (fun _ => E) x).
    + apply continuous_V_eff; assumption.
    + apply continuous_const.
  - apply continuous_sqrt.
Qed.

Lemma ex_RInt_wkb_integrand :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall E a b, ex_RInt (fun x => sqrt (V_eff (lapse f mu) x - E)) a b.
Proof.
  intros f mu Hf Hmu E a b.
  apply (@ex_RInt_continuous R_CompleteNormedModule).
  intros z _.
  apply wkb_integrand_continuous; assumption.
Qed.

Theorem V_eff_diverges :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall M, M > 0 ->
  exists X, forall x, X < x -> M < V_eff (lapse f mu) x.
Proof.
  intros f mu Hf Hmu M HM.
  assert (HsM : sqrt M + 1 > 0).
  { assert (Hs : 0 <= sqrt M) by apply sqrt_pos.
    lra. }
  destruct (inv_lapse_grows f mu Hf Hmu (sqrt M + 1) HsM) as [X HX].
  exists X.
  intros x Hx.
  unfold V_eff.
  assert (HiN : sqrt M + 1 < / lapse f mu x) by (apply HX; exact Hx).
  assert (HsqM : sqrt M < / lapse f mu x) by lra.
  assert (Hsqrtpos : 0 <= sqrt M) by apply sqrt_pos.
  assert (Hinvpos : 0 < / lapse f mu x).
  { apply inv_lapse_pos; assumption. }
  assert (Hsq : sqrt M * sqrt M < / lapse f mu x * / lapse f mu x).
  { nra. }
  replace (sqrt M * sqrt M) with M in Hsq.
  - rewrite <- Rinv_mult in Hsq.
    exact Hsq.
  - rewrite sqrt_sqrt; lra.
Qed.

Theorem forbidden_region_unbounded :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall E, E > 0 ->
  exists X_E, forall x, X_E < x -> E < V_eff (lapse f mu) x.
Proof.
  intros f mu Hf Hmu E HE.
  exact (V_eff_diverges f mu Hf Hmu E HE).
Qed.

Theorem wkb_integral_diverges :
  forall f mu,
  ActivationProps f -> MeteringProps mu ->
  forall E, E > 0 ->
  forall B, B > 0 ->
  exists L, 0 < L /\
    B < RInt (fun x => sqrt (V_eff (lapse f mu) x - E)) 0 L.
Proof.
  intros f mu Hf Hmu E HE B HB.
  assert (HEp : E + (B + 1) * (B + 1) > 0) by nra.
  destruct (V_eff_diverges f mu Hf Hmu
    (E + (B + 1) * (B + 1)) HEp) as [X HX].
  set (g := fun x => sqrt (V_eff (lapse f mu) x - E)).
  set (L := Rmax (X + 2) 1).
  exists L.
  split.
  - unfold L. apply Rlt_le_trans with 1; [lra | apply Rmax_r].
  - assert (Hg_cont : forall z, continuous g z).
    { intro z.
      unfold g.
      apply wkb_integrand_continuous; assumption. }
    assert (Hg_ex : forall a b, ex_RInt g a b).
    { intros a b.
      apply (@ex_RInt_continuous R_CompleteNormedModule).
      intros z _.
      apply Hg_cont. }
    assert (Hsub_raw : (B + 1) * ((X + 2) - (X + 1)) <=
      RInt g (X + 1) (X + 2)).
    { apply RInt_ge_const.
      - lra.
      - intros x _. apply Hg_cont.
      - intros x Hx.
        assert (HV : E + (B + 1) * (B + 1) < V_eff (lapse f mu) x).
        { apply HX. lra. }
        assert (HBp1 : 0 <= B + 1) by lra.
        replace (B + 1) with (sqrt ((B + 1) * (B + 1))).
        apply sqrt_le_1.
        + nra.
        + nra.
        + nra.
        + rewrite sqrt_square; [reflexivity | exact HBp1]. }
    assert (Hsub : B + 1 <= RInt g (X + 1) (X + 2)) by lra.
    destruct (Rle_dec 0 (X + 1)) as [H0X1 | H0X1].
    + assert (Hex1 : ex_RInt g 0 (X + 1)) by apply Hg_ex.
      assert (Hex2 : ex_RInt g (X + 1) (X + 2)) by apply Hg_ex.
      assert (Hex3 : ex_RInt g (X + 2) L) by apply Hg_ex.
      assert (Hex12 : ex_RInt g 0 (X + 2)) by apply Hg_ex.
      assert (HX2L : X + 2 <= L) by (unfold L; apply Rmax_l).
      rewrite <- (RInt_Chasles g 0 (X + 2) L Hex12 Hex3).
      rewrite <- (RInt_Chasles g 0 (X + 1) (X + 2) Hex1 Hex2).
      unfold plus; simpl.
      assert (Hge0_1 : 0 <= RInt g 0 (X + 1)).
      { apply RInt_ge_0.
        - exact H0X1.
        - exact Hex1.
        - intros x _. unfold g. apply sqrt_pos. }
      assert (Hge0_3 : 0 <= RInt g (X + 2) L).
      { apply RInt_ge_0.
        - exact HX2L.
        - exact Hex3.
        - intros x _. unfold g. apply sqrt_pos. }
      lra.
    + assert (Hex01 : ex_RInt g 0 1) by apply Hg_ex.
      assert (Hex1L : ex_RInt g 1 L) by apply Hg_ex.
      assert (H1L : 1 <= L) by (unfold L; apply Rmax_r).
      rewrite <- (RInt_Chasles g 0 1 L Hex01 Hex1L).
      unfold plus; simpl.
      assert (Hge1L : 0 <= RInt g 1 L).
      { apply RInt_ge_0.
        - exact H1L.
        - exact Hex1L.
        - intros x _. unfold g. apply sqrt_pos. }
      assert (Hge01_raw : (B + 1) * (1 - 0) <= RInt g 0 1).
      { apply RInt_ge_const.
        - lra.
        - intros x _. apply Hg_cont.
        - intros x Hx.
          assert (HV : E + (B + 1) * (B + 1) < V_eff (lapse f mu) x).
          { apply HX. lra. }
          assert (HBp1 : 0 <= B + 1) by lra.
          replace (B + 1) with (sqrt ((B + 1) * (B + 1))).
          apply sqrt_le_1.
          + nra.
          + nra.
          + nra.
          + rewrite sqrt_square; [reflexivity | exact HBp1]. }
      lra.
Qed.

(** ** Informal consequence *)
(**
   The abstract WKB theorem above is the formal confinement statement
   carried by this development: for every finite energy [E], the WKB
   barrier integral can be made arbitrarily large. That is the precise
   sense in which the barrier to escape from the metered region becomes
   infinite in this model.

   Scope note: this proves an unbounded continuum barrier under the stated
   positivity and decay hypotheses. It does not, by itself, establish any
   finite-dimensional continuum mode count.
*)

(** * Horizon: Theorem 1 -- The metering boundary is a horizon *)
(**
   The tortoise coordinate integral diverges: for any bound B,
   there exists L such that the integral of 1/N from 0 to L exceeds B.

   This is profile-independent and activation-function-independent.
   It holds for any monotone activation f with f(0) = 0 and any
   strictly positive decaying metering density mu.

   Physical interpretation: signals from the metered region take
   infinite coordinate time to reach the unmetered region. The
   metering boundary is a true causal horizon.
*)

Require Import Reals.
Require Import Lra.
Require Import Coquelicot.Coquelicot.
Require Import MeteringMetric.Definitions.
Require Import MeteringMetric.Auxiliary.
Open Scope R_scope.

(** ** Main theorem: tortoise coordinate diverges *)

Theorem horizon_existence :
  forall (f : R -> R) (mu : R -> R),
  ActivationProps f -> MeteringProps mu ->
  forall B, B > 0 ->
  exists L, 0 < L /\ B < RInt (fun x => / lapse f mu x) 0 L.
Proof.
  intros f mu Hf Hmu B HB.
  (* Strategy:
     1. By inv_lapse_grows, there exists X such that for x > X,
        1/N(x) > B + 1.
     2. Then the integral from X to X+1 is at least (B+1)*1 = B+1 > B.
     3. Since 1/N >= 0, the integral from 0 to X+1 >= integral from X to X+1.
     4. So L = max(X+1, 1) works. *)

  (* Step 1: find X where 1/N > B + 1 *)
  assert (HBp : B + 1 > 0) by lra.
  destruct (inv_lapse_grows f mu Hf Hmu (B + 1) HBp) as [X HX].

  (* We'll use L = X + 2, which is > 0 if X >= -1, but we need L > 0.
     Use max to ensure positivity. *)
  set (L := Rmax (X + 2) 1).

  exists L.
  split.
  - (* L > 0 *)
    unfold L. apply Rlt_le_trans with 1; [lra | apply Rmax_r].
  - (* B < integral from 0 to L *)
    (* Split the integral: [0, X+1] + [X+1, L] using Chasles *)
    (* First just use [X+1, X+2] which is contained in [0, L] *)

    (* Key bound: integral from X+1 to X+2 >= B+1 *)
    assert (Hsub_raw : (B + 1) * ((X + 2) - (X + 1)) <=
      RInt (fun x => / lapse f mu x) (X + 1) (X + 2)).
    { apply RInt_ge_const.
      - lra.
      - intros x _. apply continuous_inv_lapse; assumption.
      - intros x Hx.
        assert (Hxgt : X < x) by lra.
        left. apply HX. exact Hxgt. }
    assert (Hsub : B + 1 <= RInt (fun x => / lapse f mu x) (X + 1) (X + 2))
      by lra.

    (* integral from 0 to L >= integral from 0 to X+2 >= integral from X+1 to X+2 *)
    assert (HX1_le_X2 : X + 1 <= X + 2) by lra.
    assert (H0_le_X1 : 0 <= X + 1 \/ X + 1 < 0) by lra.

    (* We need: integral from 0 to L >= B + 1 > B *)
    (* Use: integral [0, L] >= integral [X+1, X+2] for appropriate arrangement *)

    (* Since 1/N > 0, integral is monotone in the integration domain *)
    (* integral [0, L] = integral [0, X+1] + integral [X+1, X+2] + integral [X+2, L] *)
    (* All three parts are >= 0, so integral [0, L] >= integral [X+1, X+2] >= B+1 *)

    assert (HX2_le_L : X + 2 <= L).
    { unfold L. apply Rmax_l. }

    assert (HintXX : ex_RInt (fun x => / lapse f mu x) (X + 1) (X + 2)).
    { apply ex_RInt_inv_lapse; assumption. }

    (* We need 0 <= X+1 for the Chasles decomposition to work simply.
       If X+1 < 0, the integral from 0 to X+2 still contains [X+1, X+2]
       but the decomposition is different. Let's handle both cases. *)

    (* Use the fact that for any a <= b <= c with positive integrand,
       integral [a,c] >= integral [b,c] *)

    assert (Hpos : forall x, 0 <= / lapse f mu x).
    { intro. left. apply inv_lapse_pos; assumption. }

    (* integral from X+1 to L >= integral from X+1 to X+2 *)
    assert (HintX1L : (B + 1) * 1 <= RInt (fun x => / lapse f mu x) (X + 1) L).
    { assert (HintX2L : ex_RInt (fun x => / lapse f mu x) (X + 2) L).
      { apply ex_RInt_inv_lapse; assumption. }
      assert (HintX1L_ex : ex_RInt (fun x => / lapse f mu x) (X + 1) L).
      { apply ex_RInt_inv_lapse; assumption. }
      rewrite <- (RInt_Chasles _ (X + 1) (X + 2) L HintXX HintX2L).
      unfold plus; simpl.
      assert (H0part : 0 <= RInt (fun x => / lapse f mu x) (X + 2) L).
      { apply RInt_ge_0.
        - exact HX2_le_L.
        - apply HintX2L.
        - intros x _. apply Hpos. }
      lra. }

    (* Now: integral from 0 to L >= integral from X+1 to L *)
    (* Case split on whether 0 <= X+1 *)
    destruct (Rle_dec 0 (X + 1)) as [H0X1 | H0X1].
    + (* 0 <= X+1: integral [0, L] = integral [0, X+1] + integral [X+1, L] *)
      assert (Hint0X1 : ex_RInt (fun x => / lapse f mu x) 0 (X + 1)).
      { apply ex_RInt_inv_lapse; assumption. }
      assert (HintX1L_ex : ex_RInt (fun x => / lapse f mu x) (X + 1) L).
      { apply ex_RInt_inv_lapse; assumption. }
      rewrite <- (RInt_Chasles _ 0 (X + 1) L Hint0X1 HintX1L_ex).
      unfold plus; simpl.
      assert (H0first : 0 <= RInt (fun x => / lapse f mu x) 0 (X + 1)).
      { apply RInt_ge_0.
        - exact H0X1.
        - exact Hint0X1.
        - intros x _. apply Hpos. }
      lra.
    + (* X+1 < 0: X < -1, so all x >= 0 satisfy x > X.
         1/N(x) > B+1 for all x in [0, 1]. Integral [0, L] >= integral [0, 1] >= B+1. *)
      assert (Hint01 : ex_RInt (fun x => / lapse f mu x) 0 1).
      { apply ex_RInt_inv_lapse; assumption. }
      assert (Hint1L : ex_RInt (fun x => / lapse f mu x) 1 L).
      { apply ex_RInt_inv_lapse; assumption. }
      assert (H1L : 1 <= L).
      { unfold L. apply Rmax_r. }
      rewrite <- (RInt_Chasles _ 0 1 L Hint01 Hint1L).
      unfold plus; simpl.
      assert (Hge1L : 0 <= RInt (fun x => / lapse f mu x) 1 L).
      { apply RInt_ge_0; [exact H1L | exact Hint1L | intros x _; apply Hpos]. }
      assert (Hge01_raw : (B + 1) * (1 - 0) <= RInt (fun x => / lapse f mu x) 0 1).
      { apply RInt_ge_const.
        - lra.
        - intros x _. apply continuous_inv_lapse; assumption.
        - intros x Hx.
          left. apply HX. lra. }
      lra.
Qed.

(** ** Corollary: tortoise coordinate is unbounded *)

Corollary tortoise_unbounded :
  forall (f : R -> R) (mu : R -> R),
  ActivationProps f -> MeteringProps mu ->
  forall B, B > 0 ->
  exists L, B < tortoise (lapse f mu) L.
Proof.
  intros f mu Hf Hmu B HB.
  destruct (horizon_existence f mu Hf Hmu B HB) as [L [HL Hint]].
  exists L.
  unfold tortoise. exact Hint.
Qed.

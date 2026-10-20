import Riemann.Deep.Hasse.Main
import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.SetTheory.Cardinal.NatCard

/-!
# Hasse's theorem: main module

* hasse_weierstrass: every Weierstrass curve over a finite field of odd characteristic has
  ((q + 1) − #E)² ≤ 4q, where #E = #{affine solutions} + 1 (no nonsingularity needed).
* hasse_elliptic: the same for Nat.card of the Mathlib point group of an elliptic curve.
* Rung 1 (CM curve y² = x³ − x over ZMod p): bound, and #E = p + 1 for p ≡ 3 (mod 4).
* Rung 2 (families y² = x³ + D, y² = x³ − Dx): bounds, and the supersingular counts
  #E = q + 1 for q ≡ 2 (mod 3), resp. q ≡ 3 (mod 4).

The proof (Riemann.Deep.Hasse.Main) is Manin's twist argument in upper-bound form.
-/

open Polynomial

namespace Hasse

open scoped Classical

variable {F : Type*} [Field F] [Fintype F]

omit [Fintype F] in
lemma eval_cub (c₂ c₁ c₀ x : F) : (cub c₂ c₁ c₀).eval x = x ^ 3 + c₂ * x ^ 2 + c₁ * x + c₀ := by
  simp [cub]

/-- Completing the square: the affine points of W correspond to those of y² = f(x). -/
lemma card_weierstrass (hF : ringChar F ≠ 2) (W : WeierstrassCurve F) :
    (Finset.univ.filter (fun p : F × F ↦ W.toAffine.Equation p.1 p.2)).card =
      (Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = (cub (W.a₂ + W.a₁ ^ 2 * (2 : F)⁻¹ ^ 2)
        (W.a₄ + 2 * W.a₁ * W.a₃ * (2 : F)⁻¹ ^ 2) (W.a₆ + W.a₃ ^ 2 * (2 : F)⁻¹ ^ 2)).eval p.1)).card := by
  have hi : (2 : F) * (2 : F)⁻¹ = 1 := mul_inv_cancel₀ (Ring.two_ne_zero hF)
  set g : F → F := fun x ↦ (W.a₁ * x + W.a₃) * (2 : F)⁻¹ with hg
  let e : F × F ≃ F × F :=
    { toFun := fun p ↦ (p.1, p.2 + g p.1)
      invFun := fun p ↦ (p.1, p.2 - g p.1)
      left_inv := fun p ↦ by simp
      right_inv := fun p ↦ by simp }
  refine Finset.card_equiv e (fun p ↦ ?_)
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, e, Equiv.coe_fn_mk, eval_cub,
    WeierstrassCurve.Affine.equation_iff, hg]
  constructor
  · intro h
    linear_combination h + ((W.a₁ * p.1 + W.a₃) * p.2) * hi
  · intro h
    linear_combination h - ((W.a₁ * p.1 + W.a₃) * p.2) * hi

/-- **Hasse's theorem** for an arbitrary Weierstrass curve over a finite field of odd
characteristic (affine count + point at infinity). -/
theorem hasse_weierstrass (hF : ringChar F ≠ 2) (W : WeierstrassCurve F) :
    ((Fintype.card F : ℤ) + 1 -
      ((Finset.univ.filter (fun p : F × F ↦ W.toAffine.Equation p.1 p.2)).card + 1)) ^ 2 ≤
      4 * Fintype.card F := by
  rw [card_weierstrass hF W]
  have hc := card_affine hF (W.a₂ + W.a₁ ^ 2 * (2 : F)⁻¹ ^ 2)
    (W.a₄ + 2 * W.a₁ * W.a₃ * (2 : F)⁻¹ ^ 2) (W.a₆ + W.a₃ ^ 2 * (2 : F)⁻¹ ^ 2)
  have hb := sum_quadraticChar_sq_le hF (W.a₂ + W.a₁ ^ 2 * (2 : F)⁻¹ ^ 2)
    (W.a₄ + 2 * W.a₁ * W.a₃ * (2 : F)⁻¹ ^ 2) (W.a₆ + W.a₃ ^ 2 * (2 : F)⁻¹ ^ 2)
  rw [hc]
  nlinarith [hb]

/-- **Hasse's theorem** for elliptic curves, stated with Mathlib's group of points:
((q + 1) − #E(F))² ≤ 4q. -/
theorem hasse_elliptic (hF : ringChar F ≠ 2) (W : WeierstrassCurve F) [W.IsElliptic] :
    ((Fintype.card F : ℤ) + 1 - Nat.card W.toAffine.Point) ^ 2 ≤ 4 * Fintype.card F := by
  have hcard : Nat.card W.toAffine.Point =
      (Finset.univ.filter (fun p : F × F ↦ W.toAffine.Equation p.1 p.2)).card + 1 := by
    rw [Nat.card_congr (WeierstrassCurve.Affine.nonsingularPointEquiv W.toAffine)]
    have e2 : {xy : F × F // W.toAffine.Nonsingular xy.1 xy.2} ≃
        {xy : F × F // W.toAffine.Equation xy.1 xy.2} :=
      Equiv.subtypeEquivRight (fun xy ↦ WeierstrassCurve.Affine.equation_iff_nonsingular.symm)
    rw [show Nat.card (WithZero {xy : F × F // W.toAffine.Nonsingular xy.1 xy.2}) =
      Nat.card (Option {xy : F × F // W.toAffine.Nonsingular xy.1 xy.2}) from rfl,
      Finite.card_option, Nat.card_congr e2, Nat.card_eq_fintype_card, Fintype.card_subtype]
  rw [hcard]
  push_cast
  exact hasse_weierstrass hF W

/-! ## Rung 2: the CM families y² = x³ + D and y² = x³ − Dx -/

/-- Bound for y² = x³ + D. -/
theorem hasse_x3_add (hF : ringChar F ≠ 2) (D : F) :
    |((Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 + D)).card + 1 : ℝ) -
      (Fintype.card F + 1)| ≤ 2 * Real.sqrt (Fintype.card F) := by
  have h := hasse hF 0 0 D
  simpa using h

/-- Bound for y² = x³ − Dx. -/
theorem hasse_x3_sub (hF : ringChar F ≠ 2) (D : F) :
    |((Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 - D * p.1)).card + 1 : ℝ) -
      (Fintype.card F + 1)| ≤ 2 * Real.sqrt (Fintype.card F) := by
  have h := hasse hF 0 (-D) 0
  have hfilt : Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 + 0 * p.1 ^ 2 + -D * p.1 + 0) =
      Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 - D * p.1) := by
    refine Finset.filter_congr (fun p _ ↦ ?_)
    constructor <;> intro h' <;> linear_combination h'
  rwa [hfilt] at h

/-- Supersingular count: if q ≡ 3 (mod 4) then y² = x³ − Dx has exactly q affine points. -/
theorem card_x3_sub_of_mod_four (hF : ringChar F ≠ 2) (hq : Fintype.card F % 4 = 3) (D : F) :
    (Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 - D * p.1)).card = Fintype.card F := by
  have hc := card_affine hF 0 (-D) 0
  have hS : ∑ a : F, quadraticChar F ((cub 0 (-D) 0).eval a) = 0 := by
    have hneg : quadraticChar F (-1) = -1 :=
      quadraticChar_neg_one_iff_not_isSquare.mpr (by rw [FiniteField.isSquare_neg_one_iff]; omega)
    have hsym : ∑ a : F, quadraticChar F ((cub 0 (-D) 0).eval a) =
        - ∑ a : F, quadraticChar F ((cub 0 (-D) 0).eval a) := by
      calc ∑ a : F, quadraticChar F ((cub 0 (-D) 0).eval a)
          = ∑ a : F, quadraticChar F ((cub 0 (-D) 0).eval ((Equiv.neg F) a)) :=
            (Equiv.sum_comp (Equiv.neg F) (fun a ↦ quadraticChar F ((cub 0 (-D) 0).eval a))).symm
        _ = ∑ a : F, -quadraticChar F ((cub 0 (-D) 0).eval a) := by
            refine Finset.sum_congr rfl (fun a _ ↦ ?_)
            have : (cub 0 (-D) 0).eval ((Equiv.neg F) a) = (-1) * (cub 0 (-D) 0).eval a := by
              simp only [Equiv.neg_apply, eval_cub]; ring
            rw [this, map_mul, hneg, neg_one_mul]
        _ = - ∑ a : F, quadraticChar F ((cub 0 (-D) 0).eval a) := Finset.sum_neg_distrib _
    linarith
  have hfilt : Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = (cub 0 (-D) 0).eval p.1) =
      Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 - D * p.1) := by
    refine Finset.filter_congr (fun p _ ↦ ?_)
    rw [eval_cub]
    constructor <;> intro h' <;> linear_combination h'
  rw [hS, add_zero, hfilt] at hc
  exact_mod_cast hc

/-- Supersingular count: if q ≡ 2 (mod 3) then y² = x³ + D has exactly q affine points. -/
theorem card_x3_add_of_mod_three (hF : ringChar F ≠ 2) (hq : Fintype.card F % 3 = 2) (D : F) :
    (Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 + D)).card = Fintype.card F := by
  have hc := card_affine hF 0 0 D
  set q := Fintype.card F with hqdef
  have hq2 : 2 ≤ q := Fintype.one_lt_card
  have hk : 3 * ((2 * q - 1) / 3) = 2 * q - 1 := by omega
  have hinv : ∀ x : F, (x ^ 3) ^ ((2 * q - 1) / 3) = x := by
    intro x
    rw [← pow_mul, hk]
    by_cases hx : x = 0
    · rw [hx, zero_pow (by omega)]
    · have h1 : x ^ (q - 1) = 1 := FiniteField.pow_card_sub_one_eq_one x hx
      have : 2 * q - 1 = 2 * (q - 1) + 1 := by omega
      rw [this, pow_succ, pow_mul', h1, one_pow, one_mul]
  have hbij : Function.Bijective (fun x : F ↦ x ^ 3) := by
    refine Finite.injective_iff_bijective.mp (fun x y hxy ↦ ?_)
    rw [← hinv x, ← hinv y]
    exact congrArg (· ^ ((2 * q - 1) / 3)) hxy
  have hS : ∑ a : F, quadraticChar F ((cub 0 0 D).eval a) = 0 := by
    calc ∑ a : F, quadraticChar F ((cub 0 0 D).eval a)
        = ∑ a : F, quadraticChar F ((Equiv.ofBijective _ hbij) a + D) := by
          refine Finset.sum_congr rfl (fun a _ ↦ ?_)
          simp [eval_cub]
      _ = ∑ u : F, quadraticChar F (u + D) :=
          Equiv.sum_comp (Equiv.ofBijective _ hbij) (fun u ↦ quadraticChar F (u + D))
      _ = ∑ u : F, quadraticChar F u :=
          Equiv.sum_comp (Equiv.addRight D) (fun u ↦ quadraticChar F u)
      _ = 0 := quadraticChar_sum_zero hF
  have hfilt : Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = (cub 0 0 D).eval p.1) =
      Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 + D) := by
    refine Finset.filter_congr (fun p _ ↦ ?_)
    simp [eval_cub]
  rw [hS, add_zero, hfilt] at hc
  exact_mod_cast hc

/-! ## Rung 1: the CM curve y² = x³ − x over ZMod p -/

/-- Hasse bound for y² = x³ − x over F_p, p an odd prime. -/
theorem rung1_bound (p : ℕ) [Fact p.Prime] (hp : p ≠ 2) :
    |((Finset.univ.filter (fun P : ZMod p × ZMod p ↦ P.2 ^ 2 = P.1 ^ 3 - P.1)).card + 1 : ℝ) -
      (p + 1)| ≤ 2 * Real.sqrt p := by
  have hF : ringChar (ZMod p) ≠ 2 := by rw [ZMod.ringChar_zmod_n]; exact hp
  have h := hasse_x3_sub hF (1 : ZMod p)
  simp only [one_mul, ZMod.card] at h
  exact h

/-- Exact count: for p ≡ 3 (mod 4), y² = x³ − x has p + 1 points over F_p (a_p = 0). -/
theorem rung1_exact (p : ℕ) [Fact p.Prime] (hp : p % 4 = 3) :
    (Finset.univ.filter (fun P : ZMod p × ZMod p ↦ P.2 ^ 2 = P.1 ^ 3 - P.1)).card = p := by
  have hF : ringChar (ZMod p) ≠ 2 := by rw [ZMod.ringChar_zmod_n]; omega
  have h := card_x3_sub_of_mod_four hF (by rw [ZMod.card]; exact hp) (1 : ZMod p)
  simp only [one_mul, ZMod.card] at h
  exact h

end Hasse

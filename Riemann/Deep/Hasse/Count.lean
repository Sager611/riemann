import Riemann.Deep.Hasse.Twist
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic

/-!
# The Frobenius point P₁ and the heights of P₁ ± P₀
-/

open Polynomial

namespace Hasse

open scoped Classical
open Function

variable {F : Type*} [Field F] [Fintype F]

local notation "ι" => algebraMap F[X] (RatFunc F)

variable {c₂ c₁ c₀ : F}

lemma card_odd (hF : ringChar F ≠ 2) : Fintype.card F = 2 * (Fintype.card F / 2) + 1 := by
  have := FiniteField.odd_card_of_char_ne_two hF; omega

lemma expand_cub : expand F (Fintype.card F) (cub c₂ c₁ c₀) =
    (X ^ Fintype.card F) ^ 3 + C c₂ * (X ^ Fintype.card F) ^ 2 + C c₁ * X ^ Fintype.card F +
      C c₀ := by
  simp [cub]

/-- Frobenius: f(t)^q = f(t^q). -/
lemma Dk_pow_card : Dk c₂ c₁ c₀ ^ Fintype.card F =
    (Tk ^ Fintype.card F) ^ 3 + ι (C c₂) * (Tk ^ Fintype.card F) ^ 2 +
      ι (C c₁) * Tk ^ Fintype.card F + ι (C c₀) := by
  rw [Dk, ← map_pow, ← FiniteField.Polynomial.expand_card, expand_cub]
  simp [Tk]

lemma P1_equation (hF : ringChar F ≠ 2) :
    (twist c₂ c₁ c₀).toAffine.Equation (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F)
      (Dk c₂ c₁ c₀ ^ 2 * Dk c₂ c₁ c₀ ^ (Fintype.card F / 2)) := by
  rw [twist_equation_iff]
  have hm := card_odd hF
  have hG : Dk c₂ c₁ c₀ * (Dk c₂ c₁ c₀ ^ (Fintype.card F / 2)) ^ 2 =
      Dk c₂ c₁ c₀ ^ Fintype.card F := by
    conv_rhs => rw [hm]
    ring
  linear_combination (Dk c₂ c₁ c₀ ^ 3) * hG +
    (Dk c₂ c₁ c₀ ^ 3) * Dk_pow_card (c₂ := c₂) (c₁ := c₁) (c₀ := c₀)

variable (c₂ c₁ c₀) in
/-- The Frobenius point P₁ = (t^q, f(t)^((q-1)/2)) of the twist. -/
noncomputable def P1 (hF : ringChar F ≠ 2) : (twist c₂ c₁ c₀).toAffine.Point :=
  .some _ _ (twist_nonsingular (Ring.two_ne_zero hF)
    (mul_ne_zero (pow_ne_zero 2 (Dk_ne_zero c₂ c₁ c₀)) (pow_ne_zero _ (Dk_ne_zero c₂ c₁ c₀)))
    (P1_equation hF))

lemma ht_P1 (hF : ringChar F ≠ 2) : ht c₂ c₁ c₀ (P1 c₂ c₁ c₀ hF) = Fintype.card F := by
  rw [P1, ht_some, mul_div_cancel_left₀ _ (Dk_ne_zero c₂ c₁ c₀), Tk, ← map_pow, htK_algebraMap,
    natDegree_X_pow]

lemma Tk_pow_card_ne : (Tk : RatFunc F) ^ Fintype.card F ≠ Tk := by
  intro h
  have h1 : htK ((Tk : RatFunc F) ^ Fintype.card F) = htK (Tk : RatFunc F) := by rw [h]
  rw [Tk, ← map_pow, htK_algebraMap, htK_algebraMap, natDegree_X_pow, natDegree_X] at h1
  have := Fintype.one_lt_card (α := F)
  omega

lemma X1_ne_X0 : Dk c₂ c₁ c₀ * Tk ^ Fintype.card F ≠ Dk c₂ c₁ c₀ * Tk := fun h ↦
  Tk_pow_card_ne (mul_left_cancel₀ (Dk_ne_zero c₂ c₁ c₀) h)

lemma E_ne_zero : (X ^ Fintype.card F - X : F[X]) ≠ 0 :=
  FiniteField.X_pow_card_sub_X_ne_zero F Fintype.one_lt_card

/-- The x-coordinate of P₁ + Q for Q = (t, s), as an explicit fraction. -/
lemma x_formula (s : F) {ℓ : RatFunc F}
    (hℓ : ℓ * (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk) =
      Dk c₂ c₁ c₀ ^ 2 * Dk c₂ c₁ c₀ ^ (Fintype.card F / 2) - ι (C s) * Dk c₂ c₁ c₀ ^ 2) :
    (ℓ ^ 2 - ι (C c₂) * Dk c₂ c₁ c₀ - Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk) /
        Dk c₂ c₁ c₀ =
      ι (cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ (Fintype.card F / 2) - C s) ^ 2 -
          (C c₂ + X ^ Fintype.card F + X) * (X ^ Fintype.card F - X) ^ 2) /
        ι ((X ^ Fintype.card F - X) ^ 2) := by
  have hD := Dk_ne_zero c₂ c₁ c₀
  have hℓ' : ℓ * (Tk ^ Fintype.card F - Tk) =
      Dk c₂ c₁ c₀ * (Dk c₂ c₁ c₀ ^ (Fintype.card F / 2) - ι (C s)) := by
    apply mul_left_cancel₀ hD
    linear_combination hℓ
  rw [div_eq_div_iff hD (RatFunc.algebraMap_ne_zero (pow_ne_zero 2 E_ne_zero))]
  simp only [map_sub, map_mul, map_add, map_pow]
  rw [← Dk, ← Tk]
  linear_combination (ℓ * (Tk ^ Fintype.card F - Tk) +
    Dk c₂ c₁ c₀ * (Dk c₂ c₁ c₀ ^ (Fintype.card F / 2) - ι (C s))) * hℓ'

/-! ## Divisibility and degrees -/

/-- Exponent of (X - a) that divides the numerator of x(P₁ + σ P₀). -/
noncomputable def ex (c₂ c₁ c₀ : F) (σ : ℤ) (a : F) : ℕ :=
  (1 + σ * quadraticChar F ((cub c₂ c₁ c₀).eval a)).toNat

lemma ex_cast (σ : ℤ) (hσ : σ = 1 ∨ σ = -1) (a : F) :
    (ex c₂ c₁ c₀ σ a : ℤ) = 1 + σ * quadraticChar F ((cub c₂ c₁ c₀).eval a) := by
  unfold ex
  by_cases h0 : (cub c₂ c₁ c₀).eval a = 0
  · rw [h0, quadraticChar_zero]; simp
  · rcases quadraticChar_dichotomy h0 with h1 | h1 <;> rcases hσ with rfl | rfl <;> rw [h1] <;>
      decide

lemma sum_ex (σ : ℤ) (hσ : σ = 1 ∨ σ = -1) :
    ((∑ a : F, ex c₂ c₁ c₀ σ a : ℕ) : ℤ) =
      Fintype.card F + σ * ∑ a : F, quadraticChar F ((cub c₂ c₁ c₀).eval a) := by
  push_cast
  simp_rw [ex_cast σ hσ]
  rw [Finset.sum_add_distrib, Finset.mul_sum]
  simp

lemma pow_half_eq (hF : ringChar F ≠ 2) (b : F) :
    b ^ (Fintype.card F / 2) = (quadraticChar F b : F) :=
  (quadraticChar_eq_pow_of_char_ne_two' hF b).symm

/-- The common divisor G = ∏ (X - a)^e(a). -/
noncomputable def Gpoly (c₂ c₁ c₀ : F) (σ : ℤ) : F[X] := ∏ a : F, (X - C a) ^ ex c₂ c₁ c₀ σ a

lemma Gpoly_natDegree (σ : ℤ) : (Gpoly c₂ c₁ c₀ σ).natDegree = ∑ a : F, ex c₂ c₁ c₀ σ a := by
  rw [Gpoly, natDegree_prod_of_monic _ _ (fun a _ ↦ (monic_X_sub_C a).pow _)]
  simp [natDegree_pow]

lemma Gpoly_dvd (hF : ringChar F ≠ 2) (σ : ℤ) (hσ : σ = 1 ∨ σ = -1) :
    Gpoly c₂ c₁ c₀ σ ∣ cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ (Fintype.card F / 2) - C (σ : F)) ^ 2 -
        (C c₂ + X ^ Fintype.card F + X) * (X ^ Fintype.card F - X) ^ 2 ∧
      Gpoly c₂ c₁ c₀ σ ∣ (X ^ Fintype.card F - X) ^ 2 := by
  have hcop : ((Finset.univ : Finset F) : Set F).Pairwise
      (IsCoprime on fun a ↦ (X - C a) ^ ex c₂ c₁ c₀ σ a) := by
    intro a _ b _ hab
    exact ((pairwise_coprime_X_sub_C (s := id) Function.injective_id) hab).pow
  have hE : ∀ a : F, X - C a ∣ (X ^ Fintype.card F - X : F[X]) := by
    intro a
    rw [dvd_iff_isRoot]
    simp [FiniteField.pow_card]
  have key : ∀ a : F, (X - C a) ^ ex c₂ c₁ c₀ σ a ∣
      cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ (Fintype.card F / 2) - C (σ : F)) ^ 2 -
        (C c₂ + X ^ Fintype.card F + X) * (X ^ Fintype.card F - X) ^ 2 ∧
      (X - C a) ^ ex c₂ c₁ c₀ σ a ∣ (X ^ Fintype.card F - X) ^ 2 := by
    intro a
    have hcast := ex_cast (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) σ hσ a
    have hE2 : (X - C a) ^ 2 ∣ (X ^ Fintype.card F - X : F[X]) ^ 2 := pow_dvd_pow_of_dvd (hE a) 2
    by_cases h0 : (cub c₂ c₁ c₀).eval a = 0
    · have he : ex c₂ c₁ c₀ σ a = 1 := by
        rw [h0, quadraticChar_zero, mul_zero, add_zero] at hcast; omega
      rw [he, pow_one]
      have hf : X - C a ∣ cub c₂ c₁ c₀ := dvd_iff_isRoot.mpr h0
      exact ⟨dvd_sub (dvd_mul_of_dvd_left hf _)
        (dvd_mul_of_dvd_right (dvd_pow (hE a) two_ne_zero) _), dvd_pow (hE a) two_ne_zero⟩
    · have hsq : ∀ τ : ℤ, quadraticChar F ((cub c₂ c₁ c₀).eval a) = τ →
          X - C a ∣ cub c₂ c₁ c₀ ^ (Fintype.card F / 2) - C (τ : F) := by
        intro τ hτ
        rw [dvd_iff_isRoot]
        simp [eval_pow, pow_half_eq hF, hτ]
      rcases quadraticChar_dichotomy h0 with h1 | h1 <;> rcases hσ with rfl | rfl
      · have he : ex c₂ c₁ c₀ 1 a = 2 := by rw [h1] at hcast; omega
        rw [he]
        have hg := hsq 1 h1
        exact ⟨dvd_sub (dvd_mul_of_dvd_right (pow_dvd_pow_of_dvd hg 2) _)
          (dvd_mul_of_dvd_right hE2 _), hE2⟩
      · have he : ex c₂ c₁ c₀ (-1) a = 0 := by rw [h1] at hcast; omega
        rw [he, pow_zero]; exact ⟨one_dvd _, one_dvd _⟩
      · have he : ex c₂ c₁ c₀ 1 a = 0 := by rw [h1] at hcast; omega
        rw [he, pow_zero]; exact ⟨one_dvd _, one_dvd _⟩
      · have he : ex c₂ c₁ c₀ (-1) a = 2 := by rw [h1] at hcast; omega
        rw [he]
        have hg := hsq (-1) h1
        exact ⟨dvd_sub (dvd_mul_of_dvd_right (pow_dvd_pow_of_dvd hg 2) _)
          (dvd_mul_of_dvd_right hE2 _), hE2⟩
  exact ⟨Finset.prod_dvd_of_coprime hcop (fun a _ ↦ (key a).1),
    Finset.prod_dvd_of_coprime hcop (fun a _ ↦ (key a).2)⟩

lemma num_natDegree (hF : ringChar F ≠ 2) (σ : ℤ) (hσ : σ = 1 ∨ σ = -1) :
    (cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ (Fintype.card F / 2) - C (σ : F)) ^ 2 -
        (C c₂ + X ^ Fintype.card F + X) * (X ^ Fintype.card F - X) ^ 2).natDegree ≤
      2 * Fintype.card F + 1 := by
  have hm : Fintype.card F = 2 * (Fintype.card F / 2) + 1 := card_odd hF
  set q := Fintype.card F with hqdef
  set m := q / 2 with hmdef
  have hq1 : 1 ≤ q := by omega
  have hfq : cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ m) ^ 2 =
      (X ^ q) ^ 3 + C c₂ * (X ^ q) ^ 2 + C c₁ * X ^ q + C c₀ := by
    have h1 : cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ m) ^ 2 = cub c₂ c₁ c₀ ^ q := by
      rw [hm]; ring
    rw [h1, hqdef, ← FiniteField.Polynomial.expand_card, expand_cub]
  have hσ2 : C (σ : F) ^ 2 = 1 := by
    rcases hσ with rfl | rfl <;> simp
  have hcub : cub c₂ c₁ c₀ = X ^ 3 + C c₂ * X ^ 2 + C c₁ * X + C c₀ := rfl
  have hN : cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ m - C (σ : F)) ^ 2 -
      (C c₂ + X ^ q + X) * (X ^ q - X) ^ 2 =
      ((X ^ q) ^ 2 * X + X ^ q * X ^ 2 + 2 * C c₂ * X ^ q * X + C c₁ * (X ^ q + X) +
        2 * C c₀) - 2 * C (σ : F) * cub c₂ c₁ c₀ ^ (m + 1) := by
    linear_combination hfq + (cub c₂ c₁ c₀) * hσ2 + hcub
  rw [hN]
  have n2 : (2 : F[X]).natDegree ≤ 0 := by simp
  have hXq : (X ^ q : F[X]).natDegree ≤ q := natDegree_X_pow_le q
  have hX1 : (X : F[X]).natDegree ≤ 1 := natDegree_X_le
  refine nd_sub ?_ ?_
  · refine nd_add (nd_add (nd_add (nd_add ?_ ?_) ?_) ?_) ?_
    · exact nd_mono (nd_mul (natDegree_pow_le_of_le 2 hXq) hX1) (by omega)
    · exact nd_mono (nd_mul hXq (natDegree_pow_le_of_le 2 hX1)) (by omega)
    · exact nd_mono (nd_mul (nd_mul (nd_mul n2 (nd_C c₂ (n := 0))) hXq) hX1) (by omega)
    · exact nd_mono (nd_mul (nd_C c₁ (n := 0)) (nd_add hXq (nd_mono hX1 hq1))) (by omega)
    · exact nd_mono (nd_mul n2 (nd_C c₀ (n := 0))) (by omega)
  · have h3 : (cub c₂ c₁ c₀ ^ (m + 1)).natDegree ≤ (m + 1) * 3 :=
      natDegree_pow_le_of_le _ (cub_natDegree c₂ c₁ c₀).le
    exact nd_mono (nd_mul (nd_mul n2 (nd_C _ (n := 0))) h3) (by omega)

lemma E2_natDegree :
    ((X ^ Fintype.card F - X : F[X]) ^ 2).natDegree ≤ 2 * Fintype.card F + 1 := by
  have h1 : (X ^ Fintype.card F - X : F[X]).natDegree ≤ Fintype.card F :=
    nd_sub (natDegree_X_pow_le _) (natDegree_X_le.trans Fintype.card_pos)
  have := natDegree_pow_le_of_le 2 h1
  omega

/-- Height bound for x(P₁ + σP₀) with the exact number of cancelled poles. -/
lemma htK_formula_le (hF : ringChar F ≠ 2) (σ : ℤ) (hσ : σ = 1 ∨ σ = -1) :
    (htK (ι (cub c₂ c₁ c₀ * (cub c₂ c₁ c₀ ^ (Fintype.card F / 2) - C (σ : F)) ^ 2 -
        (C c₂ + X ^ Fintype.card F + X) * (X ^ Fintype.card F - X) ^ 2) /
      ι ((X ^ Fintype.card F - X) ^ 2)) : ℤ) ≤
      Fintype.card F + 1 - σ * ∑ a : F, quadraticChar F ((cub c₂ c₁ c₀).eval a) := by
  obtain ⟨d1, d2⟩ := Gpoly_dvd (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) hF σ hσ
  have h := htK_add_le_of_dvd (pow_ne_zero 2 E_ne_zero) rfl d1 d2
  have hb := max_le (num_natDegree (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) hF σ hσ) E2_natDegree
  rw [Gpoly_natDegree] at h
  have hs := sum_ex (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) σ hσ
  omega

/-- ht(P₁ + P₀) ≤ q + 1 - S. -/
lemma ht_P1_add (hF : ringChar F ≠ 2) :
    (ht c₂ c₁ c₀ (P1 c₂ c₁ c₀ hF + P0 c₂ c₁ c₀ (Ring.two_ne_zero hF)) : ℤ) ≤
      Fintype.card F + 1 - ∑ a : F, quadraticChar F ((cub c₂ c₁ c₀).eval a) := by
  rw [P1, P0, WeierstrassCurve.Affine.Point.add_of_X_ne X1_ne_X0, ht_some,
    WeierstrassCurve.Affine.slope_of_X_ne X1_ne_X0]
  have hX : (twist c₂ c₁ c₀).toAffine.addX (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F)
      (Dk c₂ c₁ c₀ * Tk)
      ((Dk c₂ c₁ c₀ ^ 2 * Dk c₂ c₁ c₀ ^ (Fintype.card F / 2) - Dk c₂ c₁ c₀ ^ 2) /
        (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk)) =
      ((Dk c₂ c₁ c₀ ^ 2 * Dk c₂ c₁ c₀ ^ (Fintype.card F / 2) - Dk c₂ c₁ c₀ ^ 2) /
        (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk)) ^ 2 - ι (C c₂) * Dk c₂ c₁ c₀ -
        Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk := by
    simp [twist]
  rw [hX, x_formula ((1 : ℤ) : F)]
  · have := htK_formula_le (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) hF 1 (Or.inl rfl)
    simpa using this
  · rw [div_mul_cancel₀ _ (sub_ne_zero.mpr X1_ne_X0)]
    simp

/-- ht(P₁ - P₀) ≤ q + 1 + S. -/
lemma ht_P1_sub (hF : ringChar F ≠ 2) :
    (ht c₂ c₁ c₀ (P1 c₂ c₁ c₀ hF - P0 c₂ c₁ c₀ (Ring.two_ne_zero hF)) : ℤ) ≤
      Fintype.card F + 1 + ∑ a : F, quadraticChar F ((cub c₂ c₁ c₀).eval a) := by
  rw [sub_eq_add_neg, P1, P0, WeierstrassCurve.Affine.Point.neg_some,
    WeierstrassCurve.Affine.Point.add_of_X_ne X1_ne_X0, ht_some,
    WeierstrassCurve.Affine.slope_of_X_ne X1_ne_X0, twist_negY]
  have hX : (twist c₂ c₁ c₀).toAffine.addX (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F)
      (Dk c₂ c₁ c₀ * Tk)
      ((Dk c₂ c₁ c₀ ^ 2 * Dk c₂ c₁ c₀ ^ (Fintype.card F / 2) - -Dk c₂ c₁ c₀ ^ 2) /
        (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk)) =
      ((Dk c₂ c₁ c₀ ^ 2 * Dk c₂ c₁ c₀ ^ (Fintype.card F / 2) - -Dk c₂ c₁ c₀ ^ 2) /
        (Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk)) ^ 2 - ι (C c₂) * Dk c₂ c₁ c₀ -
        Dk c₂ c₁ c₀ * Tk ^ Fintype.card F - Dk c₂ c₁ c₀ * Tk := by
    simp [twist]
  rw [hX, x_formula ((-1 : ℤ) : F)]
  · have := htK_formula_le (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) hF (-1) (Or.inr rfl)
    simpa using this
  · rw [div_mul_cancel₀ _ (sub_ne_zero.mpr X1_ne_X0)]
    simp

end Hasse

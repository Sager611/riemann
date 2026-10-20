import Mathlib.FieldTheory.RatFunc.AsPolynomial
import Mathlib.RingTheory.Coprime.Lemmas
import Mathlib.Tactic.LinearCombination

/-!
# Heights of rational functions (projective degree) and the Vieta height bound
-/

open Polynomial

namespace Hasse

variable {F : Type*} [Field F]

/-- Projective degree of a rational function: the larger of the degrees of its reduced
numerator and denominator. -/
noncomputable def htK (x : RatFunc F) : ℕ := max x.num.natDegree x.denom.natDegree

local notation "ι" => algebraMap F[X] (RatFunc F)

lemma htK_algebraMap (p : F[X]) : htK (ι p) = p.natDegree := by
  simp [htK]

lemma htK_le_of_eq_div {x : RatFunc F} {p q : F[X]} (hq : q ≠ 0) (hx : x = ι p / ι q) :
    htK x ≤ max p.natDegree q.natDegree := by
  by_cases hp : p = 0
  · subst hp
    simp only [map_zero, zero_div] at hx
    subst hx
    simp [htK]
  · have h1 : x.num ∣ p := hx ▸ RatFunc.num_div_dvd p hq
    have h2 : x.denom ∣ q := hx ▸ RatFunc.denom_div_dvd p q
    exact max_le_max (natDegree_le_of_dvd h1 hp) (natDegree_le_of_dvd h2 hq)

/-- Cancelling a common factor `G` of `p` and `q` lowers the height bound by `deg G`. -/
lemma htK_add_le_of_dvd {x : RatFunc F} {p q G : F[X]} (hq : q ≠ 0) (hx : x = ι p / ι q)
    (hGp : G ∣ p) (hGq : G ∣ q) :
    htK x + G.natDegree ≤ max p.natDegree q.natDegree := by
  obtain ⟨p', rfl⟩ := hGp
  obtain ⟨q', rfl⟩ := hGq
  have hG : G ≠ 0 := left_ne_zero_of_mul hq
  have hq' : q' ≠ 0 := right_ne_zero_of_mul hq
  have hx' : x = ι p' / ι q' := by
    rw [hx, map_mul, map_mul, mul_div_mul_left _ _ (RatFunc.algebraMap_ne_zero hG)]
  have h1 := htK_le_of_eq_div hq' hx'
  rw [natDegree_mul hG hq']
  by_cases hp' : p' = 0
  · subst hp'
    have : htK x = 0 := by
      rw [hx']; simp [htK]
    rw [this]; omega
  · rw [natDegree_mul hG hp']
    omega

/-- Explicit Bezout identity for the coefficients of `(v₃ X - u₃)(v₄ X - u₄)`. -/
lemma bezout_three {R : Type*} [CommRing R] {u₃ v₃ u₄ v₄ : R} (h₃ : IsCoprime u₃ v₃)
    (h₄ : IsCoprime u₄ v₄) :
    ∃ a₀ a₁ a₂ : R, a₀ * (v₃ * v₄) + a₁ * (u₃ * v₄ + u₄ * v₃) + a₂ * (u₃ * u₄) = 1 := by
  obtain ⟨α, β, hab⟩ := h₃
  obtain ⟨γ, δ, hcd⟩ := h₄
  set L₁ := α * u₃ + β * v₃
  set L₂ := γ * u₄ + δ * v₄
  set e := α * δ * u₃ * v₄ + β * γ * v₃ * u₄
  refine ⟨β * δ * (L₁ * L₂ + e) + (2 * α * β * γ * δ - α ^ 2 * δ ^ 2 - β ^ 2 * γ ^ 2) * (u₃ * u₄),
    α ^ 2 * δ ^ 2 * u₃ * v₄ + β ^ 2 * γ ^ 2 * v₃ * u₄, α * γ * (L₁ * L₂ + e), ?_⟩
  have h1 : L₁ = 1 := hab
  have h2 : L₂ = 1 := hcd
  linear_combination (L₂ * (L₁ * L₂ + 1)) * h1 + (L₁ * L₂ + 1) * h2

lemma coeff_ne_zero_of_htK (x : RatFunc F) :
    x.num.coeff (htK x) ≠ 0 ∨ x.denom.coeff (htK x) ≠ 0 := by
  unfold htK
  by_cases h : x.num.natDegree ≤ x.denom.natDegree
  · right
    rw [max_eq_right h]
    exact (leadingCoeff_ne_zero.mpr x.denom_ne_zero)
  · left
    push Not at h
    rw [max_eq_left h.le]
    have hne : x.num ≠ 0 := by
      intro h0; rw [h0, natDegree_zero] at h; omega
    exact leadingCoeff_ne_zero.mpr hne

/-- **Vieta height bound.** If `A (x₃ + x₄) = B` and `A x₃ x₄ = C` with `A ≠ 0`, then the sum of
the heights of `x₃` and `x₄` is at most the largest degree among `A, B, C`. -/
theorem htK_add_htK_le {x₃ x₄ : RatFunc F} {A B C : F[X]} (hA : A ≠ 0)
    (hs : ι A * (x₃ + x₄) = ι B) (hp : ι A * (x₃ * x₄) = ι C) :
    htK x₃ + htK x₄ ≤ max A.natDegree (max B.natDegree C.natDegree) := by
  set u₃ := x₃.num with hu₃
  set v₃ := x₃.denom with hv₃
  set u₄ := x₄.num with hu₄
  set v₄ := x₄.denom with hv₄
  have hv₃0 : v₃ ≠ 0 := x₃.denom_ne_zero
  have hv₄0 : v₄ ≠ 0 := x₄.denom_ne_zero
  have e₃ : ι u₃ = x₃ * ι v₃ := by
    have := RatFunc.num_div_denom x₃
    rwa [div_eq_iff (RatFunc.algebraMap_ne_zero x₃.denom_ne_zero)] at this
  have e₄ : ι u₄ = x₄ * ι v₄ := by
    have := RatFunc.num_div_denom x₄
    rwa [div_eq_iff (RatFunc.algebraMap_ne_zero x₄.denom_ne_zero)] at this
  set T₀ := v₃ * v₄ with hT₀
  set T₁ := u₃ * v₄ + u₄ * v₃ with hT₁
  set T₂ := u₃ * u₄ with hT₂
  have hT₀0 : T₀ ≠ 0 := mul_ne_zero hv₃0 hv₄0
  have r1 : A * T₁ = B * T₀ := by
    apply RatFunc.algebraMap_injective F
    simp only [T₁, T₀, map_mul, map_add]
    rw [e₃, e₄]
    linear_combination (ι v₃ * ι v₄) * hs
  have r2 : A * T₂ = C * T₀ := by
    apply RatFunc.algebraMap_injective F
    simp only [T₂, T₀, map_mul]
    rw [e₃, e₄]
    linear_combination (ι v₃ * ι v₄) * hp
  obtain ⟨a₀, a₁, a₂, hbez⟩ := bezout_three (x₃.isCoprime_num_denom) (x₄.isCoprime_num_denom)
  set c := a₀ * A + a₁ * B + a₂ * C with hc
  have hc0 : c * T₀ = A := by
    linear_combination (-a₁) * r1 + (-a₂) * r2 + A * hbez
  have hc1 : c * T₁ = B := by
    apply mul_right_cancel₀ hT₀0
    linear_combination T₁ * hc0 + r1
  have hc2 : c * T₂ = C := by
    apply mul_right_cancel₀ hT₀0
    linear_combination T₂ * hc0 + r2
  have hcne : c ≠ 0 := by
    rintro h0; rw [h0, zero_mul] at hc0; exact hA hc0.symm
  have dle : ∀ T : F[X], T.natDegree ≤ (c * T).natDegree := by
    intro T
    by_cases hT : T = 0
    · simp [hT]
    · rw [natDegree_mul hcne hT]; omega
  -- leading forms
  set n₃ := htK x₃
  set n₄ := htK x₄
  have hu₃d : u₃.natDegree ≤ n₃ := le_max_left _ _
  have hv₃d : v₃.natDegree ≤ n₃ := le_max_right _ _
  have hu₄d : u₄.natDegree ≤ n₄ := le_max_left _ _
  have hv₄d : v₄.natDegree ≤ n₄ := le_max_right _ _
  have k0 : T₀.coeff (n₃ + n₄) = v₃.coeff n₃ * v₄.coeff n₄ :=
    coeff_mul_add_eq_of_natDegree_le hv₃d hv₄d
  have k2 : T₂.coeff (n₃ + n₄) = u₃.coeff n₃ * u₄.coeff n₄ :=
    coeff_mul_add_eq_of_natDegree_le hu₃d hu₄d
  have k1 : T₁.coeff (n₃ + n₄) = u₃.coeff n₃ * v₄.coeff n₄ + u₄.coeff n₄ * v₃.coeff n₃ := by
    rw [coeff_add, coeff_mul_add_eq_of_natDegree_le hu₃d hv₄d, add_comm n₃ n₄,
      coeff_mul_add_eq_of_natDegree_le hu₄d hv₃d]
  have key : T₀.coeff (n₃ + n₄) ≠ 0 ∨ T₁.coeff (n₃ + n₄) ≠ 0 ∨ T₂.coeff (n₃ + n₄) ≠ 0 := by
    rw [k0, k1, k2]
    have h3 := coeff_ne_zero_of_htK x₃
    have h4 := coeff_ne_zero_of_htK x₄
    by_contra hcon
    push Not at hcon
    obtain ⟨z0, z1, z2⟩ := hcon
    rcases h3 with h3 | h3 <;> rcases h4 with h4 | h4
    · exact (mul_ne_zero h3 h4) z2
    · rcases mul_eq_zero.mp z0 with h | h
      · rw [h, mul_zero, add_zero] at z1; exact (mul_ne_zero h3 h4) z1
      · exact h4 h
    · rcases mul_eq_zero.mp z0 with h | h
      · exact h3 h
      · rw [h, mul_zero, zero_add] at z1; exact (mul_ne_zero h4 h3) z1
    · exact (mul_ne_zero h3 h4) z0
  rcases key with h | h | h
  · have := le_natDegree_of_ne_zero h
    have := dle T₀; rw [hc0] at this; omega
  · have := le_natDegree_of_ne_zero h
    have := dle T₁; rw [hc1] at this; omega
  · have := le_natDegree_of_ne_zero h
    have := dle T₂; rw [hc2] at this; omega

end Hasse

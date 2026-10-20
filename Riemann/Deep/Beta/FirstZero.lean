import Riemann.Deep.Beta.Cert

/-!
# The first zero of the Dirichlet beta function on the critical line, certified in Lean

**Theorem.** There is t ∈ (6.01, 6.03) with L(1/2 + it, χ₋₄) = 0, χ₋₄ = Mathlib's ZMod.χ₄ in ℂ
(the first zero is at t ≈ 6.020948904697597).

Proof. F(t) = e^{i t log 4 / 2} Λ(1/2 + it, χ₋₄) satisfies conj F = ν F (Reality.lean: conjugation
symmetry from the first-order Euler–Maclaurin formula, plus the odd functional equation). At
t₁ = 6.012, t₂ = 6.028 the order-8 Euler–Maclaurin enclosures of ζ(s, 1/4) - ζ(s, 3/4) (Cert.lean)
give L(s₂) conj L(s₁) ≈ -1.08·10⁻⁴ ± 2.8·10⁻⁶, and e^{iδ log 4} Γℝ(s₂+1)/Γℝ(s₁+1) is within 0.148
of 1 (GammaRatio.lean). So Re F(t₂) conj F(t₁) < 0 and zero_of_conj_sign gives a zero of F, hence
of Λ, hence of L, in between.
-/

open Complex Set

namespace Fayed.Beta

open Fayed.FirstZero

lemma chiM4_ne_one : chiM4 ≠ 1 := by
  intro h
  have h3 := chiM4_vals.2.2.2
  have hu : IsUnit (3 : ZMod 4) := ⟨⟨3, 3, by decide, by decide⟩, rfl⟩
  rw [h, MulChar.one_apply hu] at h3
  norm_num at h3

lemma Fb_continuous : Continuous Fb := by
  have hd := DirichletCharacter.differentiable_completedLFunction chiM4_ne_one
  unfold Fb
  exact (Continuous.cexp (by fun_prop)).mul (hd.continuous.comp (by fun_prop))

lemma re_pos_add_one (t : ℝ) : 0 < ((1 / 2 : ℂ) + (t : ℂ) * I + 1).re := by
  simp; norm_num

/-- **The first zero of the Dirichlet beta function L(s, χ₋₄) on the critical line lies in
(6.01, 6.03).** -/
theorem first_zero : ∃ t : ℝ, 6.01 < t ∧ t < 6.03 ∧
    DirichletCharacter.LFunction χ₋₄ (1 / 2 + t * I) = 0 := by
  have hlt : (1503 / 250 : ℝ) < 1507 / 250 := by norm_num
  have hsign : (Fb (1507 / 250) * starRingEnd ℂ (Fb (1503 / 250))).re < 0 := by
    have hL1 := Cert.L_a
    have hL2 := Cert.L_b
    have hρ := GammaR1_ratio (t1 := 1503 / 250) (t2 := 1507 / 250) (δ := 1 / 125) (by norm_num)
      (by norm_num) (by norm_num) (by norm_num) (S := 126167 / 1000000) (by norm_num) (by norm_num)
      (by norm_num) (c := 147131 / 1000000) (by norm_num)
    have key := re_mul_conj_neg_gen (by norm_num) hL1 hL2 hρ
      (R := 2760143973 / 1000000000000000) (by norm_num) (by norm_num)
    set s1 : ℂ := (1 / 2 : ℂ) + (((1503 / 250 : ℝ)) : ℂ) * I with hs1
    set s2 : ℂ := (1 / 2 : ℂ) + (((1507 / 250 : ℝ)) : ℂ) * I with hs2
    set G1 := Gammaℝ (s1 + 1)
    set G2 := Gammaℝ (s2 + 1)
    have hG1 : G1 ≠ 0 := Gammaℝ_ne_zero_of_re_pos (re_pos_add_one _)
    have hexp : cexp (I * (((1507 / 250 : ℝ) * Real.log 4 / 2 : ℝ) : ℂ)) *
        starRingEnd ℂ (cexp (I * (((1503 / 250 : ℝ) * Real.log 4 / 2 : ℝ) : ℂ))) =
        cexp (I * (((1 / 125 : ℝ) * Real.log 4 : ℝ) : ℂ)) := by
      rw [← exp_conj, ← Complex.exp_add]
      congr 1
      rw [map_mul, conj_I, conj_ofReal]
      push_cast
      ring
    have e : Fb (1507 / 250) * starRingEnd ℂ (Fb (1503 / 250)) =
        (DirichletCharacter.LFunction chiM4 s2 * starRingEnd ℂ (DirichletCharacter.LFunction chiM4 s1)
          * (cexp (I * (((1 / 125 : ℝ) * Real.log 4 : ℝ) : ℂ)) * (G2 / G1))) *
          ((normSq G1 : ℝ) : ℂ) := by
      simp only [Fb]
      rw [← hs1, ← hs2, completed_eq_mul (re_pos_add_one _), completed_eq_mul (re_pos_add_one _),
        ← mul_conj, ← hexp]
      simp only [map_mul]
      have hc : G2 / G1 * G1 = G2 := div_mul_cancel₀ G2 hG1
      linear_combination (-(DirichletCharacter.LFunction chiM4 s2 *
        starRingEnd ℂ (DirichletCharacter.LFunction chiM4 s1) *
        cexp (I * (((1507 / 250 : ℝ) * Real.log 4 / 2 : ℝ) : ℂ)) *
        starRingEnd ℂ (cexp (I * (((1503 / 250 : ℝ) * Real.log 4 / 2 : ℝ) : ℂ))) *
        starRingEnd ℂ G1)) * hc
    rw [e, re_mul_ofReal]
    exact mul_neg_of_neg_of_pos key (normSq_pos.mpr hG1)
  obtain ⟨t, ht, hF⟩ := zero_of_conj_sign hlt Fb_continuous.continuousOn Fb_conj hsign
  refine ⟨t, ?_, ?_, ?_⟩
  · have := ht.1; norm_num at this ⊢; linarith
  · have := ht.2; norm_num at this ⊢; linarith
  · have hΛ : DirichletCharacter.completedLFunction chiM4 ((1 / 2 : ℂ) + (t : ℂ) * I) = 0 := by
      simp only [Fb] at hF
      rcases mul_eq_zero.mp hF with h | h
      · exact absurd h (Complex.exp_ne_zero _)
      · exact h
    rw [completed_eq_mul (re_pos_add_one t)] at hΛ
    rcases mul_eq_zero.mp hΛ with h | h
    · exact h
    · exact absurd h (Gammaℝ_ne_zero_of_re_pos (re_pos_add_one t))

end Fayed.Beta

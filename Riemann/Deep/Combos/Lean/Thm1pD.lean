module

public import Riemann.Deep.Combos.Lean.Thm1pLg
public import Riemann.Deep.Combos.Lean.Thm1pSm
public import Riemann.Deep.Combos.Lean.Thm1pK

@[expose] public section

/-! # Theorem 1′(a): the bound for one D (both 𝔔 cases combined), term mode

Arithmetic is in `Thm1pK` (no `Real.rpow` import); here only lemma applications. -/

open Real Finset

namespace RingClass

/-- the constant of Theorem 1′(a): depends only on ε, c₀, C_I, C_II, 𝔔₀ -/
noncomputable def C1p (ε c₀ : ℝ) (CI CII Q0 : ℝ → ℝ) : ℝ :=
  81 * (c₀ * (2 / (ε / 10)) * (CI (ε / 10) * (8 / (ε / 10) ^ 2) + CII (ε / 10) + 1 +
    |3 * (Q0 (ε / 10) / 9) ^ (1 / 2 : ℝ) * Real.log (3 * (Q0 (ε / 10) / 9)) ^ 2 + 1|))

lemma per_D {X : Type*} [Fintype X] [DecidableEq X] (d : RCZeros X) {c₀ c₁ : ℝ}
    {CI CII Q0 : ℝ → ℝ} (h : d.Hyp c₀ c₁ CI CII Q0) (ε : ℝ) (hε : 0 < ε) (hε1 : ε ≤ 1)
    (σ t : ℝ) (h0 : 0.6 ≤ σ) (h1 : σ ≤ 1) :
    (d.N σ t : ℝ) ≤ C1p ε c₀ CI CII Q0 * (d.absD * (|t| + 3)) ^ (ee σ + ε) := by
  obtain ⟨he, he1⟩ := eps_facts ε hε hε1
  have hT := T_ge t
  have hD := one_le_of_three_le d.absD d.absD_ge
  have hP : (9 * d.absD * (|t| + 3)) ^ (ee σ + 5 * (ε / 10))
      ≤ 81 * (d.absD * (|t| + 3)) ^ (ee σ + ε) := by
    rw [eps5 (ee σ) ε]
    exact final_exp d.absD (|t| + 3) (ee σ) ε hD hT (ee_nonneg σ h1 h0) (ee_le σ h1 h0) hε hε1
  unfold C1p
  exact case_combine _ _ _ _ _ _ _ (A_nonneg c₀ (ε / 10) h.c₀_nonneg he)
    (Real.rpow_pos_of_pos (W_pos d.absD (|t| + 3) hD hT) _).le
    (KI_nonneg (CI (ε / 10)) (ε / 10) (h.CI_nonneg _)) (h.CII_nonneg _)
    ((le_or_gt (Q0 (ε / 10)) (9 * d.absD)).imp
      (fun hQ => large_case d h σ t (ε / 10) h0 h1 he he1 hQ)
      (fun hQ => small_case d h σ t (ε / 10) h0 h1 he hQ)) hP

end RingClass

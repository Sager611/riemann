import Mathlib

/-!
# Round 83 — a Lean counterexample library of false lemmas used in RH-proof claims

| pattern | counterexample |
|---------|----------------|
| f ~ g ⇒ f' ~ g' | f = x + sin x, g = x |
| converse L'Hôpital (f/g → L ⇒ f'/g' → L) | same; disconnected-domain version: Fayed.theorem1C_false |
| lim ↔ Σ swap without uniformity | f n k = [k = n] (below) |
| sign of a sum after rearrangement | alternating harmonic series → any value (Riemann rearrangement) |
| positive on a dense set ⇒ positive | indicator of ℚ minus indicator of ℝ \ ℚ |
| Hurwitz without local uniformity | 1 − z/(n+1) → 1, each has a zero |
| Fubini without integrability | f(m, n) = [n = m] − [n = m + 1] |
| inequality on Re s > 1 extended to the strip | exp(1 − s) (below) |
-/

namespace Fayed.R83

/-- A limit cannot be swapped with an infinite sum without uniformity. -/
theorem lim_tsum_swap_false : ∃ f : ℕ → ℕ → ℝ,
    (∀ k, Filter.Tendsto (fun n => f n k) Filter.atTop (nhds 0)) ∧ ∀ n, ∑' k, f n k = 1 := by
  refine ⟨fun n k => if k = n then 1 else 0, fun k => ?_, fun n => ?_⟩
  · apply tendsto_const_nhds.congr'
    filter_upwards [Filter.eventually_gt_atTop k] with n hn
    simp [Nat.ne_of_lt hn]
  · simp

/-- An inequality valid on Re s > 1 need not survive analytic continuation into the strip. -/
theorem continuation_inequality_false : ∃ f : ℂ → ℂ, Differentiable ℂ f ∧
    (∀ s : ℂ, 1 < s.re → ‖f s‖ ≤ 1) ∧ 1 < ‖f 0‖ := by
  refine ⟨fun s => Complex.exp (1 - s), ?_, fun s hs => ?_, ?_⟩
  · exact ((differentiable_const (1 : ℂ)).sub differentiable_id).cexp
  · rw [Complex.norm_exp]
    simp only [Complex.sub_re, Complex.one_re]
    exact Real.exp_le_one_iff.mpr (by linarith)
  · rw [Complex.norm_exp]
    simp only [Complex.sub_re, Complex.one_re, Complex.zero_re, sub_zero]
    have h := Real.add_one_le_exp (1 : ℝ)
    linarith

end Fayed.R83

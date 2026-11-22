module

public import Riemann.Deep.LeanLi.Thm2
public import Riemann.Deep.LeanLi.Corollaries

/-!
# LeanLi: Claim part 3 and the unconditional upper bound of part 1 (from Theorem 2 and Lemma 0)
-/

@[expose] public section

open Real Finset Filter

namespace LeanLi

/-- **Claim part 1, second clause** (Theorem 2). -/
theorem part1_upper (L : LitCore) (hRH : RiemannHypothesis) :
    ∃ C : ℝ, ∀ᶠ N : ℕ in atTop,
      ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤ C * ((N : ℝ) ^ (3 / 2 : ℝ) * Real.log N ^ (3 / 2 : ℝ)) := by
  obtain ⟨C, N₀, h⟩ := thm2_upper L hRH
  refine ⟨C, ?_⟩
  filter_upwards [eventually_ge_atTop (max N₀ 1)] with N hN
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast le_trans (le_max_right _ _) hN
  rw [rpow_three_halves _ (by linarith), rpow_three_halves _ (Real.log_nonneg hN1)]
  have := h N (le_trans (le_max_left _ _) hN)
  calc ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤ C * (N * √(N : ℝ)) * (Real.log N * √(Real.log N)) := this
    _ = C * (N * √(N : ℝ) * (Real.log N * √(Real.log N))) := by ring

/-- **Claim part 3.** -/
theorem part3 (L : LitCore) : Part3 := by
  refine ⟨?_, ?_, part3_half L⟩
  · simp only [summable_keiper_iff]
    exact LR.rh_iff_exists R RiemannHypothesis (thm2_upper L) (lemma0_R L)
  · simp only [summable_keiper_iff]
    exact LR.rh_iff_forall R RiemannHypothesis (thm2_upper L) (lemma0_R L)

end LeanLi


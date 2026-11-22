import Riemann.Deep.LeanNBRate.Bridge

/-!
# Assembly of the main claim from the upper bound
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The upper bound for Nyman's constrained distance (to be proved by U1 + U2 + closure). -/
def UpperBoundNyman : Prop :=
  RiemannHypothesis → ∀ ε > 0, ∀ᶠ lam in 𝓝[>] (0 : ℝ), Real.log (1 / lam) * nymanDist lam ^ 2 ≤ zeroSum + ε

theorem mainClaim_of_upper (hB : Burnol2002_Thm) (hc : Titchmarsh_Thm9_2) (hD : Davenport_Ch12)
    (hU : UpperBoundNyman) : MainClaim := by
  intro hRH
  have hev0 : ∀ᶠ lam in 𝓝[>] (0 : ℝ), 0 < lam := self_mem_nhdsWithin
  have hev1 : ∀ᶠ lam in 𝓝[>] (0 : ℝ), lam < 1 := by
    filter_upwards [Ioo_mem_nhdsGT (show (0 : ℝ) < 1 by norm_num)] with lam hlam using hlam.2
  have hmono : ∀ᶠ lam in 𝓝[>] (0 : ℝ),
      Real.log (1 / lam) * nbDist lam ^ 2 ≤ Real.log (1 / lam) * nymanDist lam ^ 2 := by
    filter_upwards [hev0, hev1] with lam h0 h1
    have hL : 0 ≤ Real.log (1 / lam) := Real.log_nonneg (by rw [le_div_iff₀ h0]; linarith)
    have h := nbDist_le_nymanDist lam h0
    have hn : 0 ≤ nbDist lam := Metric.infDist_nonneg
    exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hn h 2) hL
  refine ⟨summable_zeroSum hc hRH, ?_, fun hs => zeroSum_value hD hRH hs, ?_⟩
  · rw [tendsto_order]
    refine ⟨fun a ha => ?_, fun a ha => ?_⟩
    · have hε : 0 < (zeroSum - a) / 2 := by linarith
      filter_upwards [lower_bound hB hc hRH _ hε] with lam h
      linarith
    · have hε : 0 < (a - zeroSum) / 2 := by linarith
      filter_upwards [hU hRH _ hε, hmono] with lam h1 h2
      linarith
  · rw [tendsto_order]
    refine ⟨fun a ha => ?_, fun a ha => ?_⟩
    · have hε : 0 < (zeroSum - a) / 2 := by linarith
      filter_upwards [lower_bound hB hc hRH _ hε, hmono] with lam h1 h2
      linarith
    · have hε : 0 < (a - zeroSum) / 2 := by linarith
      filter_upwards [hU hRH _ hε] with lam h1
      linarith

end Riemann.Deep.LeanNBRate


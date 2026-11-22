import Riemann.Deep.LeanNBRate.Final4

/-!
# UpperBoundNyman: L·D_N(e^{-L})² ≤ S + ε eventually, then λ = e^{-L} → 0
-/

open MeasureTheory Filter Topology Set Complex
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- **The upper bound for Nyman's constrained distance.** -/
theorem upper_bound (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) (h97 : Titchmarsh_Thm9_7)
    (hT : Nazarov1993_Turan) (hGG : GoldstonGonek2007_Thm1) (hc : Titchmarsh_Thm9_2) : UpperBoundNyman := by
  intro hRH ε hε
  have hS := zeroSum_nonneg
  obtain ⟨θ, hθdef⟩ : ∃ θ : ℝ, θ = min 1 (ε / (2 * (4 + 3 * zeroSum))) := ⟨_, rfl⟩
  have hθ0 : 0 < θ := by rw [hθdef]; exact lt_min one_pos (by positivity)
  have hθ1 : θ ≤ 1 := by rw [hθdef]; exact min_le_left _ _
  have hθε : θ ≤ ε / (2 * (4 + 3 * zeroSum)) := by rw [hθdef]; exact min_le_right _ _
  obtain ⟨ℓ₀, -, hU2⟩ := u2 hT hGG hc hRH θ hθ0
  obtain ⟨A, C, T₀, -, h97'⟩ := h97
  have hev : ∀ᶠ L : ℝ in atTop, L * nymanDist (Real.exp (-L)) ^ 2 ≤ zeroSum + ε := by
    have e1 := eventually_ge_atTop (max (max 256 (max T₀ 1)) (max (2 * ℓ₀) (max (4 * Real.exp 5) (pE A : ℝ))))
    have e2 := (tendsto_order.1 tendsto_L_div_ell).2 (1 + θ) (by linarith)
    have e3 := (tendsto_order.1 (tendsto_err C A θ)).2 (ε / 2) (by linarith)
    filter_upwards [e1, e2, e3] with L hL1 hL2 hL3
    have hL256 : 256 ≤ L := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hL1
    have hLT : max T₀ 1 ≤ L := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hL1
    have hLℓ : 2 * ℓ₀ ≤ L := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hL1
    have hLe : 4 * Real.exp 5 ≤ L :=
      le_trans (le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) (le_max_right _ _)) hL1
    have hLp : (pE A : ℝ) ≤ L :=
      le_trans (le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) (le_max_right _ _)) hL1
    have hL1' : 1 ≤ L := by linarith
    have hsL : Real.sqrt L ≤ L / 2 := by
      have h1 : L ≤ (L / 2) ^ 2 := by nlinarith
      calc Real.sqrt L ≤ Real.sqrt ((L / 2) ^ 2) := Real.sqrt_le_sqrt h1
        _ = L / 2 := Real.sqrt_sq (by linarith)
    have hLT₀ : T₀ ≤ L ^ 6 := (le_max_left T₀ 1).trans (hLT.trans (le_self_pow₀ hL1' (by norm_num)))
    have hLℓ₀ : ℓ₀ ≤ L - Real.sqrt L := by linarith
    have hL7 : L ^ 6 + 1 ≤ (L - Real.sqrt L) ^ 7 := by
      have h1 : L / 2 ≤ L - Real.sqrt L := by linarith
      have h2 : (L / 2) ^ 7 ≤ (L - Real.sqrt L) ^ 7 := pow_le_pow_left₀ (by linarith) h1 7
      have h3 : L ^ 6 + 1 ≤ (L / 2) ^ 7 := by
        have e : (L / 2) ^ 7 = L ^ 6 * (L / 128) := by ring
        rw [e]
        have h4 : 2 ≤ L / 128 := by linarith
        have h5 : 1 ≤ L ^ 6 := one_le_pow₀ hL1'
        nlinarith
      linarith
    have hL5 : 4 * Real.exp 5 ≤ L ^ 3 := hLe.trans (le_self_pow₀ hL1' (by norm_num))
    have hb := bound_at hP h215 hRH hθ0 hθ1 hU2 h97' (by linarith) hLT₀ hLℓ₀ hL7 hL5 hLp
    have h1 : (1 + θ) * (zeroSum + θ) * (L / (L - Real.sqrt L)) ≤ zeroSum + ε / 2 := by
      have h4 : (1 + θ) * (zeroSum + θ) * (L / (L - Real.sqrt L)) ≤ (1 + θ) * (zeroSum + θ) * (1 + θ) :=
        mul_le_mul_of_nonneg_left hL2.le (by positivity)
      have hθ2 : θ ^ 2 ≤ θ := by nlinarith
      have hθ2' : θ ^ 2 ≤ 1 := by nlinarith
      have hθ3 : θ ^ 3 ≤ θ := by nlinarith
      have hθS : θ ^ 2 * zeroSum ≤ θ * zeroSum := mul_le_mul_of_nonneg_right hθ2 hS
      have h5 : (1 + θ) * (zeroSum + θ) * (1 + θ) ≤ zeroSum + θ * (4 + 3 * zeroSum) := by
        nlinarith [hθ2, hθ3, hθS]
      have h6 : θ * (4 + 3 * zeroSum) ≤ ε / 2 := by
        have h7 := (le_div_iff₀ (by positivity)).mp hθε
        linarith
      linarith
    linarith
  have hlog : Tendsto (fun lam : ℝ => Real.log (1 / lam)) (𝓝[>] 0) atTop := by
    have h2 := tendsto_neg_atBot_atTop.comp Real.tendsto_log_nhdsGT_zero
    refine h2.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with lam hlam
    simp [Function.comp, one_div, Real.log_inv]
  filter_upwards [hlog.eventually hev, self_mem_nhdsWithin] with lam h hlam
  have hlam' : (0 : ℝ) < lam := hlam
  have e : Real.exp (-Real.log (1 / lam)) = lam := by
    rw [one_div, Real.log_inv, neg_neg, Real.exp_log hlam']
  rw [e] at h
  exact h

end Riemann.Deep.LeanNBRate


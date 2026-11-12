module

public import Riemann.Deep.Combos.Lean.Thm1pW

@[expose] public section

/-! # Theorem 1′(a), small 𝔔 (𝔔 < 𝔔₀(ε₀)): RCApps2 l.33 ("we may assume 𝔔 ≥ 𝔔₀(ε)") -/

open Real Finset

namespace RingClass

/-- small 𝔔: N(σ;t) ≤ c₀(2/ε₀)|K_S| W^{e(σ)+5ε₀}, K_S = 3 q^{1/2} log²(3q) + 1, q = 𝔔₀(ε₀)/9 -/
lemma small_case {X : Type*} [Fintype X] [DecidableEq X] (d : RCZeros X) {c₀ c₁ : ℝ}
    {CI CII Q0 : ℝ → ℝ} (h : d.Hyp c₀ c₁ CI CII Q0) (σ t ε₀ : ℝ) (h0 : 0.6 ≤ σ) (h1 : σ ≤ 1)
    (he : 0 < ε₀) (hQ : 9 * d.absD < Q0 ε₀) :
    (d.N σ t : ℝ) ≤ c₀ * (2 / ε₀) *
      |3 * (Q0 ε₀ / 9) ^ (1 / 2 : ℝ) * Real.log (3 * (Q0 ε₀ / 9)) ^ 2 + 1| *
      (9 * d.absD * (|t| + 3)) ^ (ee σ + 5 * ε₀) := by
  set T := |t| + 3 with hTdef
  have hT : 1 ≤ T := by have := abs_nonneg t; linarith only [this, hTdef]
  have hD3 := d.absD_ge
  set W := 9 * d.absD * T with hWdef
  have hW : 1 ≤ W := one_le_mul_of_one_le_of_one_le (by linarith only [hD3]) hT
  have hW0 : 0 < W := lt_of_lt_of_le one_pos hW
  have hE : 0 ≤ W ^ ε₀ := (Real.rpow_pos_of_pos hW0 _).le
  have hc₀ := h.c₀_nonneg
  have hA0 : 0 ≤ c₀ * (2 / ε₀) := mul_nonneg hc₀ (div_nonneg zero_le_two he.le)
  have hN := d.N_le σ t (c₀ * (2 / ε₀) * W ^ ε₀) (fun χ => d.window h ε₀ he χ t)
  have hB0 : 0 ≤ c₀ * (2 / ε₀) * W ^ ε₀ := mul_nonneg hA0 hE
  have hee0 := ee_nonneg σ h1 h0
  have f3 : W ^ ε₀ ≤ W ^ (ee σ + 5 * ε₀) :=
    Real.rpow_le_rpow_of_exponent_le hW (by linarith only [hee0, he])
  have hq : d.absD ≤ Q0 ε₀ / 9 := by linarith only [hQ]
  have hx : ((d.Xs σ t).card : ℝ) ≤ Fintype.card X := by exact_mod_cast Finset.card_le_univ _
  have hF3 := h.F3
  have hs1 : d.absD ^ (1 / 2 : ℝ) ≤ (Q0 ε₀ / 9) ^ (1 / 2 : ℝ) :=
    Real.rpow_le_rpow (by linarith only [hD3]) hq (by norm_num)
  have hl0 : 0 ≤ Real.log (3 * d.absD) := Real.log_nonneg (by linarith only [hD3])
  have hl1 : Real.log (3 * d.absD) ≤ Real.log (3 * (Q0 ε₀ / 9)) :=
    Real.log_le_log (by linarith only [hD3]) (by linarith only [hq])
  have hl2 : Real.log (3 * d.absD) ^ 2 ≤ Real.log (3 * (Q0 ε₀ / 9)) ^ 2 :=
    pow_le_pow_left₀ hl0 hl1 2
  have hs0 : 0 ≤ d.absD ^ (1 / 2 : ℝ) := Real.rpow_nonneg (by linarith only [hD3]) _
  have hm0 := mul_le_mul hs1 hl2 (sq_nonneg _) (le_trans hs0 hs1)
  have hKS1 : (Fintype.card X : ℝ) + 1 ≤
      3 * (Q0 ε₀ / 9) ^ (1 / 2 : ℝ) * Real.log (3 * (Q0 ε₀ / 9)) ^ 2 + 1 := by
    have e1 : 3 * d.absD ^ (1 / 2 : ℝ) * Real.log (3 * d.absD) ^ 2
        = 3 * (d.absD ^ (1 / 2 : ℝ) * Real.log (3 * d.absD) ^ 2) := by ring
    have e2 : 3 * (Q0 ε₀ / 9) ^ (1 / 2 : ℝ) * Real.log (3 * (Q0 ε₀ / 9)) ^ 2
        = 3 * ((Q0 ε₀ / 9) ^ (1 / 2 : ℝ) * Real.log (3 * (Q0 ε₀ / 9)) ^ 2) := by ring
    linarith only [hF3, hm0, e1, e2]
  exact small_arith _ _ _ _ c₀ ε₀ _ (W ^ ε₀) (W ^ (ee σ + 5 * ε₀)) he hc₀ hE hN hx hKS1
    hB0 le_rfl f3

end RingClass


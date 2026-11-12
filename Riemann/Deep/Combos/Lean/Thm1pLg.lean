import Riemann.Deep.Combos.Lean.Thm1pW

/-! # Theorem 1′(a), large 𝔔 (𝔔 ≥ 𝔔₀(ε₀)): RCApps2 l.35–61 -/

open Real Finset

namespace RingClass

/-- large 𝔔: N(σ;t) ≤ c₀(2/ε₀)(K_I + K_II + 1) W^{e(σ)+5ε₀}, W = 9|D|T -/
lemma large_case {X : Type*} [Fintype X] [DecidableEq X] (d : RCZeros X) {c₀ c₁ : ℝ}
    {CI CII Q0 : ℝ → ℝ} (h : d.Hyp c₀ c₁ CI CII Q0) (σ t ε₀ : ℝ) (h0 : 0.6 ≤ σ) (h1 : σ ≤ 1)
    (he : 0 < ε₀) (he1 : ε₀ ≤ 1 / 10) (hQ : Q0 ε₀ ≤ 9 * d.absD) :
    (d.N σ t : ℝ) ≤ c₀ * (2 / ε₀) * (CI ε₀ * (8 / ε₀ ^ 2) + CII ε₀ + 1) *
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
  have hKI0 : 0 ≤ CI ε₀ * (8 / ε₀ ^ 2) :=
    mul_nonneg (h.CI_nonneg ε₀) (div_nonneg (by norm_num) (sq_nonneg _))
  have hQ4 : 4 ≤ 9 * d.absD := by linarith only [hD3]
  obtain ⟨hV2, hVY⟩ := V_le_Y (9 * d.absD) T σ hQ4 hT h0 h1
  have hx := d.Xs_card_le σ t
    (fun χ ρ => 1 / 4 ≤ d.detI ((9 * d.absD) ^ (1 / 2 : ℝ)) (W ^ (1 / (3 - 2 * σ))) χ ρ)
    (fun χ ρ => c₁ * (W ^ (1 / (3 - 2 * σ))) ^ (σ - 1 / 2) ≤
      d.detII ((9 * d.absD) ^ (1 / 2 : ℝ)) χ ρ.im)
    (fun χ hχ ρ hρ hσ => h.L5 σ h0 h1 _ _ hV2 hVY χ hχ t ρ hρ hσ)
    _ _
    (fun S ρ hS => h.classI ε₀ he he1 hQ σ h0 h1 t _ hVY S ρ hS)
    (fun S ρ hS => h.classII ε₀ he he1 σ h0 h1 t _ hVY S ρ hS)
  have e2 := classII_exp (9 * d.absD) T σ ε₀ (by linarith only [hD3]) hT h1
  have hBII : CII ε₀ * (W ^ (1 / (3 - 2 * σ))) ^ (1 / 2 - σ) * (9 * d.absD) ^ (1 / 2 + ε₀) *
      T ^ (1 / 2 + ε₀) = CII ε₀ * W ^ (ee σ + ε₀) := by
    rw [← e2]; ring
  have c1 := classI_exp W σ ε₀ hW h0 h1 he he1
  have m := mul_le_mul_of_nonneg_left c1 (mul_nonneg (h.CI_nonneg ε₀) hE)
  have e4 : W ^ ε₀ * W ^ (ee σ + 3 * ε₀) = W ^ (ee σ + 4 * ε₀) := by
    rw [← Real.rpow_add hW0]; congr 1; ring
  have e6 : CI ε₀ * W ^ ε₀ * (8 / ε₀ ^ 2 * W ^ (ee σ + 3 * ε₀))
      = CI ε₀ * (8 / ε₀ ^ 2) * (W ^ ε₀ * W ^ (ee σ + 3 * ε₀)) := by ring
  rw [e6, e4] at m
  have f1 : W ^ ε₀ * W ^ (ee σ + 4 * ε₀) = W ^ (ee σ + 5 * ε₀) := by
    rw [← Real.rpow_add hW0]; congr 1; ring
  have f2 : W ^ ε₀ * W ^ (ee σ + ε₀) ≤ W ^ (ee σ + 5 * ε₀) := by
    rw [← Real.rpow_add hW0]
    exact Real.rpow_le_rpow_of_exponent_le hW (by linarith only [he])
  have f3 : W ^ ε₀ ≤ W ^ (ee σ + 5 * ε₀) :=
    Real.rpow_le_rpow_of_exponent_le hW (by linarith only [hee0, he])
  exact large_arith _ _ _ _ _ c₀ ε₀ (CI ε₀ * (8 / ε₀ ^ 2)) (CII ε₀) (W ^ ε₀)
    (W ^ (ee σ + 4 * ε₀)) (W ^ (ee σ + ε₀)) (W ^ (ee σ + 5 * ε₀)) he hc₀ hKI0 (h.CII_nonneg ε₀)
    (Real.rpow_pos_of_pos hW0 _).le (Real.rpow_pos_of_pos hW0 _).le f1 f2 f3 hN hx hB0
    le_rfl m (le_of_eq hBII)

end RingClass


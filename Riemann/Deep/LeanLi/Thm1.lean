module

public import Riemann.Deep.LeanLi.Tiny

/-!
# LeanLi: **Theorem 1** (RH): `Σ_{n≤N} R(n)² ≥ (c₁ - ε) N^{3/2} (log N)^{1/2}`
-/

@[expose] public section

open Real Finset Filter

namespace LeanLi

lemma trendT_eq_liTrend : trendT = liTrend := by funext n; rfl

lemma rvm_ZK (L : LitCore) : LR.RvM ZK := by
  obtain ⟨C, T₀, _, h⟩ := count_asymp L
  exact ⟨C, T₀, fun T hT => by rw [card_ZK]; exact h T hT⟩

/-- The flat-top parameter for quality `η`. -/
noncomputable def aη (η : ℝ) : ℝ := 1 / 2 - min η 1 / 4

lemma aη_pos {η : ℝ} (hη : 0 < η) : 0 < aη η := by
  unfold aη; have := min_le_right η 1; linarith

lemma aη_lt {η : ℝ} (hη : 0 < η) : aη η < 1 / 2 := by
  unfold aη; have := lt_min hη one_pos; linarith

/-- The weights of quality `η`. -/
noncomputable def wη (η : ℝ) (hη : 0 < η) : ℕ → ℕ → ℝ :=
  fun N n => wts (aη η) (aη_pos hη) (aη_lt hη) N n

lemma good_wη (η : ℝ) (hη : 0 < η) : LR.GoodWeights (wη η hη) η := by
  refine ⟨fun N n => wts_nonneg _ _ _ N n, fun N n => wts_le_one _ _ _ N n,
    fun N x => phi_nonneg _ _ _ N x, ?_⟩
  set m := min η 1 with hm
  have hm0 : 0 < m := lt_min hη one_pos
  have hmη : m ≤ η := min_le_left _ _
  filter_upwards [eventually_ge_atTop (max 1 ⌈2 / m⌉₊)] with N hN
  have hN1 : 1 ≤ N := le_trans (le_max_left _ _) hN
  have hNm : 2 / m ≤ N := le_trans (Nat.le_ceil _) (by exact_mod_cast le_trans (le_max_right _ _) hN)
  have h := phi_zero_ge (aη η) (aη_pos hη) (aη_lt hη) hN1
  show (1 - η) * N ≤ LR.Phi N (wts (aη η) (aη_pos hη) (aη_lt hη) N) 0
  have h2 : 1 ≤ m / 2 * N := by
    rw [div_le_iff₀ hm0] at hNm; linarith
  have h3 : m / 2 * N ≤ (η - m / 2) * N := by
    apply mul_le_mul_of_nonneg_right (by linarith) (Nat.cast_nonneg _)
  have ha : 2 * aη η * N = (1 - m / 2) * N := by unfold aη; rw [← hm]; ring
  linarith

lemma phiGood_wη (η : ℝ) (hη : 0 < η) : ∃ C : ℝ, ∀ N : ℕ, 1 ≤ N → PhiGood N (wη η hη N) C := by
  obtain ⟨C, _, hC⟩ := phi_decay (aη η) (aη_pos hη) (aη_lt hη) 15
  refine ⟨C, fun N hN => ⟨wts_zero _ _ _ hN, fun n => wts_nonneg _ _ _ N n,
    fun n => wts_le_one _ _ _ N n, fun x => phi_nonneg _ _ _ N x, fun x hx => ?_⟩⟩
  have := hC N hN x hx
  show LR.Phi N (wts (aη η) (aη_pos hη) (aη_lt hη) N) x * (1 + N * |x|) ^ 30 ≤ C * N
  simpa using this

lemma trendPt_wη (η : ℝ) (hη : 0 < η) : TrendPt' (wη η hη) thK liTrend ZK := by
  obtain ⟨C, hC⟩ := phiGood_wη η hη
  obtain ⟨c, hc0, _, hc⟩ := exists_im_lower
  set δ₀ := π - liTheta c with hδ₀
  have hδ0 : 0 < δ₀ := by have := liTheta_lt_pi hc0; linarith
  refine ⟨Kc C δ₀, ?_⟩
  filter_upwards [eventually_ge_atTop (max 2 ⌈(1 / δ₀) ^ 2⌉₊)] with N hN
  intro T _ hT2 j hj
  have hN2 : 2 ≤ N := le_trans (le_max_left _ _) hN
  have hN1 : 1 ≤ N := by omega
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  have hNδ : (1 / δ₀) ^ 2 ≤ N := le_trans (Nat.le_ceil _) (by exact_mod_cast le_trans (le_max_right _ _) hN)
  have hs : 0 < √(N : ℝ) := Real.sqrt_pos.2 hNpos
  have hsδ : 1 / δ₀ ≤ √(N : ℝ) := by
    rw [Real.le_sqrt (by positivity) hNpos.le]; exact hNδ
  have hδ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2 := by
    rw [div_le_iff₀ (by positivity)]
    rw [div_le_iff₀ hδ0] at hsδ
    nlinarith
  have hγ0 : 0 < ordK j := j.1.2.2
  have hγT : ordK j ≤ T := mem_ZK.1 hj
  have ha : 1 / √(N : ℝ) ≤ thK j := liTheta_ge_inv_sqrt hN1 hγ0 (le_trans hγT hT2)
  have ha' : thK j ≤ π - δ₀ := by
    rw [hδ₀, sub_sub_cancel]
    exact liTheta_anti (hc _ j.1.2.1 j.1.2.2)
  have h := trendT_lower (hC N hN1) hN2 ha ha' hδ
  rw [trendT_eq_liTrend] at h
  exact h

open Classical in
/-- **Theorem 1** (RH; proofs.md Theorem 1): for every `ε > 0`, eventually
`Σ_{n≤N} R(n)² ≥ (c₁ - ε) N^{3/2} (log N)^{1/2}`. -/
theorem thm1_li (L : LitCore) (hRH : RiemannHypothesis) : LR.Thm1Bound R := by
  have h := thm1_of_lemmas' thK lam liTrend (hasSum_lam L hRH) (by simp [liTrend]) ZK (rvm_ZK L)
    (tauT L) (fun η hη => by
      obtain ⟨C, hC⟩ := phiGood_wη η hη
      exact ⟨wη η hη, fun N => ZK ((N : ℝ) ^ 8), good_wη η hη,
        tiny_terms L hRH (wη η hη) C hC η hη, trendPt_wη η hη⟩)
  exact h

end LeanLi

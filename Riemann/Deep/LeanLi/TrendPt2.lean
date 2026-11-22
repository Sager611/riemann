module

public import Riemann.Deep.LeanLi.LR.Trend

/-!
# LeanLi: (H4) with error `C √N (1 + log N)` suffices for Theorem 1

`TrendPt'` weakens LiRate's `TrendPt` (error `C log N` per zero) to `C √N (1 + log N)` per zero; the
total `N(T) · C √N log N = O(T² log³ T)` is still `o(T³ L²)`.
-/

@[expose] public section

open Finset Real Filter

namespace LeanLi

/-- (H4') proofs.md step 4 with error `C √N (1 + log N)`. -/
def TrendPt' {κ : Type*} (w : ℕ → ℕ → ℝ) (θ : κ → ℝ) (t : ℕ → ℝ) (Z : ℝ → Finset κ) : Prop :=
  ∃ C : ℝ, ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 → ∀ j ∈ Z T,
    Real.log (2 * π * θ j) / θ j ^ 2 - C * (√(N : ℝ) * (1 + Real.log N)) ≤ LR.trend N (w N) t (θ j)

lemma log_le_two_sqrt {T : ℝ} (hT : 0 < T) : Real.log T ≤ 2 * √T := by
  have hs : 0 < √T := Real.sqrt_pos.2 hT
  have h := Real.log_le_sub_one_of_pos hs
  have e : Real.log T = 2 * Real.log √T := by
    rw [← Real.log_rpow hs, show ((2 : ℝ)) = ((2 : ℕ) : ℝ) by norm_num, Real.rpow_natCast,
      Real.sq_sqrt hT.le]
  linarith

/-- The trend half of (H_τ) from (H4'), (H7), RvM. -/
theorem htrend_of' {κ : Type*} (w : ℕ → ℕ → ℝ) (θ : κ → ℝ) (t : ℕ → ℝ) (Z : ℝ → Finset κ)
    (hN : LR.RvM Z) (h4 : TrendPt' w θ t Z) (h7 : LR.TauT θ Z) (η : ℝ) (hη : 0 < η) :
    ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 →
      -∑ j ∈ Z T, LR.trend N (w N) t (θ j) ≤
        (1 + η / 2) * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π)) := by
  obtain ⟨C₄, h4⟩ := h4
  obtain ⟨C, T₀, hC⟩ := hN
  obtain ⟨T₇, hT₇⟩ := eventually_atTop.1 (h7 (η / 4) (by linarith))
  set K : ℝ := 3456 * π * (|C₄| + 1) / η with hK
  have hKpos : 0 < K := by positivity
  set T₁ : ℝ := max (max (max T₀ ((2 * π) ^ 2)) (2 * |C| + 1)) (max T₇ (K ^ 2 + 1)) with hT₁
  have hT₁s : K ^ 2 + 1 ≤ T₁ := le_trans (le_max_right _ _) (le_max_right _ _)
  have hT₁pos : 0 < T₁ := lt_of_lt_of_le (by positivity) hT₁s
  filter_upwards [h4, LR.eventually_basic 1 (1 / T₁) (2 / 3) (one_div_pos.2 hT₁pos) (by norm_num)]
    with N hN4 hNb
  obtain ⟨hl1, hll, hsq⟩ := hNb
  intro T hT1 hT2
  have hℓ0 : 0 < Real.log N := by linarith
  have hx0 : (0 : ℝ) < N := by
    by_contra h; push Not at h
    have : (N : ℝ) = 0 := le_antisymm h (Nat.cast_nonneg N)
    rw [this, Real.log_zero] at hl1; linarith
  have hsx : 0 < √(N : ℝ) := Real.sqrt_pos.2 hx0
  have hTbig : T₁ ≤ T := by
    have h : T₁ * Real.log N ≤ √N := by
      have := mul_le_mul_of_nonneg_left hsq hT₁pos.le
      have e : T₁ * (1 / T₁ * √(N : ℝ)) = √N := by field_simp
      linarith
    have : T₁ ≤ √N / Real.log N := by rw [le_div_iff₀ hℓ0]; linarith
    linarith
  have hTpos : 0 < T := lt_of_lt_of_le hT₁pos hTbig
  have hlogN : Real.log N ≤ 3 * Real.log T := by
    have h : Real.log (√N / Real.log N) ≤ Real.log T := Real.log_le_log (div_pos hsx hℓ0) hT1
    rw [Real.log_div hsx.ne' hℓ0.ne', Real.log_sqrt hx0.le] at h
    linarith
  obtain ⟨hcount, hL⟩ := LR.count_le Z hC (le_trans (le_max_left _ _) hTbig)
  have hT7 : T₇ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hTbig
  have hTK : K ^ 2 + 1 ≤ T := le_trans hT₁s hTbig
  have hT1' : 1 ≤ T := by nlinarith
  have hlogT : 0 ≤ Real.log T := Real.log_nonneg hT1'
  have hsum : -∑ j ∈ Z T, LR.trend N (w N) t (θ j) ≤
      -∑ j ∈ Z T, Real.log (2 * π * θ j) / θ j ^ 2 +
        (Z T).card * (C₄ * (√(N : ℝ) * (1 + Real.log N))) := by
    have h := sum_le_sum fun j hj => hN4 T hT1 hT2 j hj
    rw [sum_sub_distrib, sum_const, nsmul_eq_mul] at h
    linarith
  have h7' := hT₇ T hT7
  set L := Real.log (T / (2 * π))
  have hL0 : 0 ≤ L := le_trans (div_nonneg hlogT (by norm_num)) hL
  have herr : (Z T).card * (C₄ * (√(N : ℝ) * (1 + Real.log N))) ≤ η / 4 * (T ^ 3 * L ^ 2 / (6 * π)) := by
    have hc0 : (0 : ℝ) ≤ (Z T).card := Nat.cast_nonneg _
    have hsN : √(N : ℝ) ≤ T * Real.log N := by
      have := mul_le_mul_of_nonneg_right hT1 hℓ0.le
      rw [div_mul_cancel₀ _ hℓ0.ne'] at this; exact this
    have h1l : 1 + Real.log N ≤ 2 * Real.log N := by linarith
    have hE : C₄ * (√(N : ℝ) * (1 + Real.log N)) ≤ 18 * |C₄| * (T * Real.log T ^ 2) := by
      calc C₄ * (√(N : ℝ) * (1 + Real.log N)) ≤ |C₄| * (√(N : ℝ) * (1 + Real.log N)) :=
            mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
        _ ≤ |C₄| * ((T * Real.log N) * (2 * Real.log N)) := by
            apply mul_le_mul_of_nonneg_left _ (abs_nonneg _)
            exact mul_le_mul hsN h1l (by positivity) (by positivity)
        _ ≤ |C₄| * ((T * (3 * Real.log T)) * (2 * (3 * Real.log T))) := by
            apply mul_le_mul_of_nonneg_left _ (abs_nonneg _)
            apply mul_le_mul (mul_le_mul_of_nonneg_left hlogN hTpos.le) (by linarith) (by positivity)
              (by positivity)
        _ = 18 * |C₄| * (T * Real.log T ^ 2) := by ring
    have e2 : (Z T).card * (C₄ * (√(N : ℝ) * (1 + Real.log N))) ≤
        T * Real.log T * (18 * |C₄| * (T * Real.log T ^ 2)) := by
      calc (Z T).card * (C₄ * (√(N : ℝ) * (1 + Real.log N))) ≤
            (Z T).card * (18 * |C₄| * (T * Real.log T ^ 2)) := mul_le_mul_of_nonneg_left hE hc0
        _ ≤ T * Real.log T * (18 * |C₄| * (T * Real.log T ^ 2)) :=
            mul_le_mul_of_nonneg_right hcount (by positivity)
    have hLsq : Real.log T ^ 2 / 4 ≤ L ^ 2 := by
      have := pow_le_pow_left₀ (div_nonneg hlogT (by norm_num)) hL 2
      linarith [show (Real.log T / 2) ^ 2 = Real.log T ^ 2 / 4 by ring]
    -- 1728 π |C₄| log T ≤ η T
    have hlog2 := log_le_two_sqrt hTpos
    have hsT : K ≤ √T := by
      rw [Real.le_sqrt hKpos.le hTpos.le]; linarith
    have hkey : 1728 * π * |C₄| * Real.log T ≤ η * T := by
      have h1 : 1728 * π * |C₄| * Real.log T ≤ 3456 * π * (|C₄| + 1) * √T := by
        have : 1728 * π * |C₄| * Real.log T ≤ 1728 * π * |C₄| * (2 * √T) :=
          mul_le_mul_of_nonneg_left hlog2 (by positivity)
        have h0 : 0 ≤ 3456 * π * √T := by positivity
        have : 1728 * π * |C₄| * (2 * √T) + 3456 * π * √T = 3456 * π * (|C₄| + 1) * √T := by ring
        linarith
      have h2 : 3456 * π * (|C₄| + 1) * √T = η * K * √T := by
        rw [hK]; field_simp
      have h3 : η * K * √T ≤ η * (√T * √T) := by
        have := mul_le_mul_of_nonneg_left hsT (mul_nonneg hη.le (Real.sqrt_nonneg T))
        have e : η * √T * K = η * K * √T := by ring
        have e' : η * √T * √T = η * (√T * √T) := by ring
        linarith
      rw [Real.mul_self_sqrt hTpos.le] at h3
      linarith
    have e3 : T * Real.log T * (18 * |C₄| * (T * Real.log T ^ 2)) ≤ η / 4 * (T ^ 3 * L ^ 2 / (6 * π)) := by
      have hTl : 0 ≤ T ^ 2 * Real.log T ^ 2 := by positivity
      have key : (1728 * π * |C₄| * Real.log T) * (T ^ 2 * Real.log T ^ 2) ≤
          (η * T) * (T ^ 2 * Real.log T ^ 2) := mul_le_mul_of_nonneg_right hkey hTl
      have hL2 : η * T ^ 3 * (Real.log T ^ 2 / 4) ≤ η * T ^ 3 * L ^ 2 :=
        mul_le_mul_of_nonneg_left hLsq (mul_pos hη (pow_pos hTpos 3)).le
      rw [show η / 4 * (T ^ 3 * L ^ 2 / (6 * π)) = η * T ^ 3 * L ^ 2 / (24 * π) by ring,
        le_div_iff₀ (by positivity)]
      have e4 : T * Real.log T * (18 * |C₄| * (T * Real.log T ^ 2)) * (24 * π) =
          (1728 * π * |C₄| * Real.log T) * (T ^ 2 * Real.log T ^ 2) / 4 := by ring
      have e5 : (η * T) * (T ^ 2 * Real.log T ^ 2) = 4 * (η * T ^ 3 * (Real.log T ^ 2 / 4)) := by ring
      rw [e4]
      linarith [key, hL2]
    linarith
  linarith

/-- **Theorem 1** from the lemmas, with (H4') in place of (H4). -/
theorem thm1_of_lemmas' {κ : Type*} [DecidableEq κ] (θ : κ → ℝ) (lam t : ℕ → ℝ)
    (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)) (ht0 : t 0 = 0)
    (Z : ℝ → Finset κ) (hN : LR.RvM Z) (h7 : LR.TauT θ Z)
    (hW : ∀ η > 0, ∃ w : ℕ → ℕ → ℝ, ∃ Xs : ℕ → Finset κ,
      LR.GoodWeights w η ∧ LR.TinyTerms w θ lam Z Xs η ∧ TrendPt' w θ t Z) :
    LR.Thm1Bound (fun n => lam n - t n) :=
  LR.thm1 θ lam t hlam ht0 Z hN fun η hη => by
    obtain ⟨w, Xs, hg, ⟨hsub, htiny⟩, h4⟩ := hW η hη
    exact ⟨w, hg, LR.trendTail_of w θ lam t Z η hg.2.2.1 Xs hsub htiny
      (htrend_of' w θ t Z hN h4 h7 η hη)⟩

end LeanLi

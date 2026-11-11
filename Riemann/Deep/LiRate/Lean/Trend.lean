import Riemann.Deep.LiRate.Lean.Reduce

/-!
# LiRate: the trend half of (H_τ) from step 4 and step 7 (proofs.md lines 42–45, 51, 57)

(H4) `TrendPt`: `τ(θ_γ) ≥ log(2πθ_γ)/θ_γ² - C log N` for `γ ≤ T` (step 4: `τ(θ) = log(2πθ)/θ² +
O(log N)` on `[N^{-1/2}, 1]`, which contains `θ_γ` for `γ ≤ T ≤ √N/2`).
(H7) `TauT`: `τ_T = -Σ_{γ≤T} log(2πθ_γ)/θ_γ² ≤ (1+η) T³L²/(6π)` for large `T` (line 51, step 7:
`τ_T = T³L²/(6π)(1 + O(1/L))`, partial summation with `N(T)`).
With RvM, the error `N(T) · C log N` is `o(T³L²)` uniformly for `√N/log N ≤ T ≤ √N/2`.
-/

open Finset Real Filter

namespace LiRate

/-- (H4) proofs.md step 4. -/
def TrendPt {κ : Type*} (w : ℕ → ℕ → ℝ) (θ : κ → ℝ) (t : ℕ → ℝ) (Z : ℝ → Finset κ) : Prop :=
  ∃ C : ℝ, ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 → ∀ j ∈ Z T,
    Real.log (2 * π * θ j) / θ j ^ 2 - C * Real.log N ≤ trend N (w N) t (θ j)

/-- (H7) proofs.md line 51 and step 7: `τ_T ≤ (1 + η) T³L²/(6π)` eventually in `T`. -/
def TauT {κ : Type*} (θ : κ → ℝ) (Z : ℝ → Finset κ) : Prop :=
  ∀ η > 0, ∀ᶠ T : ℝ in atTop,
    -∑ j ∈ Z T, Real.log (2 * π * θ j) / θ j ^ 2 ≤
      (1 + η) * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π))

/-- Facts at a large `T`: `#Z(T) ≤ T log T` (RvM) and `log(T/2π) ≥ ½ log T`. -/
lemma count_le {κ : Type*} (Z : ℝ → Finset κ) {C T₀ T : ℝ}
    (hC : ∀ T ≥ T₀, |((Z T).card : ℝ) - T / (2 * π) * Real.log (T / (2 * π * exp 1))| ≤
      C * Real.log T)
    (hT : max (max T₀ ((2 * π) ^ 2)) (2 * |C| + 1) ≤ T) :
    ((Z T).card : ℝ) ≤ T * Real.log T ∧ Real.log T / 2 ≤ Real.log (T / (2 * π)) := by
  have hT0 : T₀ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hT
  have hT2π : (2 * π) ^ 2 ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hT
  have hTC : 2 * |C| + 1 ≤ T := le_trans (le_max_right _ _) hT
  have hπ : 3 < π := Real.pi_gt_three
  have h2π : 0 < 2 * π := by positivity
  have hT1 : 1 ≤ T := by linarith [abs_nonneg C]
  have hTpos : 0 < T := by linarith
  have hlogT : 0 ≤ Real.log T := Real.log_nonneg hT1
  constructor
  · have h := (abs_le.1 (hC T hT0)).2
    have hle : Real.log (T / (2 * π * exp 1)) ≤ Real.log T :=
      Real.log_le_log (div_pos hTpos (by positivity)) (div_le_self hTpos.le
        (by nlinarith [Real.add_one_le_exp 1]))
    have h1 : T / (2 * π) * Real.log (T / (2 * π * exp 1)) ≤ T / 2 * Real.log T := by
      have hd : T / (2 * π) ≤ T / 2 := div_le_div_of_nonneg_left hTpos.le (by norm_num) (by linarith)
      rcases le_total 0 (Real.log (T / (2 * π * exp 1))) with h0 | h0
      · calc T / (2 * π) * Real.log (T / (2 * π * exp 1)) ≤ T / 2 * Real.log (T / (2 * π * exp 1)) :=
              mul_le_mul_of_nonneg_right hd h0
          _ ≤ T / 2 * Real.log T := mul_le_mul_of_nonneg_left hle (by linarith)
      · have : T / (2 * π) * Real.log (T / (2 * π * exp 1)) ≤ 0 :=
          mul_nonpos_of_nonneg_of_nonpos (div_pos hTpos h2π).le h0
        have : 0 ≤ T / 2 * Real.log T := mul_nonneg (by linarith) hlogT
        linarith
    have h2 : C * Real.log T ≤ T / 2 * Real.log T :=
      mul_le_mul_of_nonneg_right (by linarith [le_abs_self C]) hlogT
    linarith
  · rw [Real.log_div hTpos.ne' h2π.ne']
    have h := Real.log_le_log (by positivity) hT2π
    rw [Real.log_pow] at h
    push_cast at h
    linarith

/-- The trend half of (H_τ): `-Σ_{γ≤T} τ(θ_γ) ≤ (1 + η/2) T³L²/(6π)`, from (H4), (H7), RvM. -/
theorem htrend_of {κ : Type*} (w : ℕ → ℕ → ℝ) (θ : κ → ℝ) (t : ℕ → ℝ) (Z : ℝ → Finset κ)
    (hN : RvM Z) (h4 : TrendPt w θ t Z) (h7 : TauT θ Z) (η : ℝ) (hη : 0 < η) :
    ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 →
      -∑ j ∈ Z T, trend N (w N) t (θ j) ≤
        (1 + η / 2) * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π)) := by
  obtain ⟨C₄, h4⟩ := h4
  obtain ⟨C, T₀, hC⟩ := hN
  obtain ⟨T₇, hT₇⟩ := eventually_atTop.1 (h7 (η / 4) (by linarith))
  set T₁ : ℝ := max (max (max T₀ ((2 * π) ^ 2)) (2 * |C| + 1))
    (max T₇ (√(288 * π * (|C₄| + 1) / η) + 1)) with hT₁
  have hT₁s : √(288 * π * (|C₄| + 1) / η) + 1 ≤ T₁ := le_trans (le_max_right _ _) (le_max_right _ _)
  have hT₁pos : 0 < T₁ := lt_of_lt_of_le (by positivity) hT₁s
  have hT₁C : 2 * |C| + 1 ≤ T₁ := le_trans (le_max_right _ _) (le_max_left _ _)
  filter_upwards [h4, eventually_basic 1 (1 / T₁) (2 / 3) (one_div_pos.2 hT₁pos) (by norm_num)]
    with N hN4 hNb
  obtain ⟨hl1, hll, hsq⟩ := hNb
  intro T hT1 hT2
  have hℓ0 : 0 < Real.log N := by linarith
  have hx0 : (0 : ℝ) < N := by
    by_contra h; push Not at h
    have : (N : ℝ) = 0 := le_antisymm h (Nat.cast_nonneg N)
    rw [this, Real.log_zero] at hl1; linarith
  have hsx : 0 < √(N : ℝ) := Real.sqrt_pos.2 hx0
  -- T ≥ T₁
  have hTbig : T₁ ≤ T := by
    have h : T₁ * Real.log N ≤ √N := by
      have := mul_le_mul_of_nonneg_left hsq hT₁pos.le
      have e : T₁ * (1 / T₁ * √(N : ℝ)) = √N := by field_simp
      linarith
    have : T₁ ≤ √N / Real.log N := by rw [le_div_iff₀ hℓ0]; linarith
    linarith
  have hTpos : 0 < T := lt_of_lt_of_le hT₁pos hTbig
  -- log N ≤ 3 log T
  have hlogN : Real.log N ≤ 3 * Real.log T := by
    have h : Real.log (√N / Real.log N) ≤ Real.log T := Real.log_le_log (div_pos hsx hℓ0) hT1
    rw [Real.log_div hsx.ne' hℓ0.ne', Real.log_sqrt hx0.le] at h
    linarith
  obtain ⟨hcount, hL⟩ := count_le Z hC (le_trans (le_max_left _ _) hTbig)
  have hT7 : T₇ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hTbig
  have hsq' : √(288 * π * (|C₄| + 1) / η) + 1 ≤ T := le_trans hT₁s hTbig
  have hT2' : 288 * π * (|C₄| + 1) / η ≤ T ^ 2 := by
    have h0 : 0 ≤ √(288 * π * (|C₄| + 1) / η) := Real.sqrt_nonneg _
    have h1 : √(288 * π * (|C₄| + 1) / η) ≤ T := by linarith
    calc 288 * π * (|C₄| + 1) / η = √(288 * π * (|C₄| + 1) / η) ^ 2 :=
          (Real.sq_sqrt (by positivity)).symm
      _ ≤ T ^ 2 := pow_le_pow_left₀ h0 h1 2
  have hT2'' : 288 * π * (|C₄| + 1) ≤ η * T ^ 2 := by
    rw [div_le_iff₀ hη] at hT2'; linarith
  have hT1' : 1 ≤ T := by linarith [abs_nonneg C]
  have hlogT : 0 ≤ Real.log T := Real.log_nonneg hT1'
  -- the pointwise bound summed
  have hsum : -∑ j ∈ Z T, trend N (w N) t (θ j) ≤
      -∑ j ∈ Z T, Real.log (2 * π * θ j) / θ j ^ 2 + (Z T).card * (C₄ * Real.log N) := by
    have h := sum_le_sum fun j hj => hN4 T hT1 hT2 j hj
    rw [sum_sub_distrib, sum_const, nsmul_eq_mul] at h
    linarith
  have h7' := hT₇ T hT7
  -- the error term
  set L := Real.log (T / (2 * π))
  have hL0 : 0 ≤ L := le_trans (div_nonneg hlogT (by norm_num)) hL
  have herr : (Z T).card * (C₄ * Real.log N) ≤ η / 4 * (T ^ 3 * L ^ 2 / (6 * π)) := by
    have hc0 : (0 : ℝ) ≤ (Z T).card := Nat.cast_nonneg _
    have e1 : (Z T).card * (C₄ * Real.log N) ≤ (Z T).card * (|C₄| * (3 * Real.log T)) :=
      mul_le_mul_of_nonneg_left (mul_le_mul (le_abs_self C₄) hlogN hℓ0.le (abs_nonneg _)) hc0
    have e2 : (Z T).card * (|C₄| * (3 * Real.log T)) ≤ T * Real.log T * (|C₄| * (3 * Real.log T)) :=
      mul_le_mul_of_nonneg_right hcount (mul_nonneg (abs_nonneg _) (by linarith))
    have hLsq : Real.log T ^ 2 / 4 ≤ L ^ 2 := by
      have := pow_le_pow_left₀ (div_nonneg hlogT (by norm_num)) hL 2
      linarith [show (Real.log T / 2) ^ 2 = Real.log T ^ 2 / 4 by ring]
    have e3 : T * Real.log T * (|C₄| * (3 * Real.log T)) ≤ η / 4 * (T ^ 3 * L ^ 2 / (6 * π)) := by
      have hTl : 0 ≤ T * Real.log T ^ 2 := mul_nonneg hTpos.le (sq_nonneg _)
      have key : 288 * π * |C₄| * (T * Real.log T ^ 2) ≤ η * T ^ 2 * (T * Real.log T ^ 2) :=
        mul_le_mul_of_nonneg_right (by linarith [mul_pos Real.pi_pos (by norm_num : (0:ℝ) < 288)]) hTl
      have hL2 : η * T ^ 3 * (Real.log T ^ 2 / 4) ≤ η * T ^ 3 * L ^ 2 :=
        mul_le_mul_of_nonneg_left hLsq (mul_pos hη (pow_pos hTpos 3)).le
      rw [show η / 4 * (T ^ 3 * L ^ 2 / (6 * π)) = η * T ^ 3 * L ^ 2 / (24 * π) by ring,
        le_div_iff₀ (by positivity)]
      linarith [key, hL2]
    linarith
  linarith

/-- (H_tiny) proofs.md step 5 (lines 48–49, 52): `Z T ⊆ Xs N` (`Xs N` = the zeros with
`γ ≤ X = N³`), and the terms that the rapid decay of `Φ` off `0` makes tiny are
`≤ (η/2) T³L²/(6π)`: `2N(T)Σ_{γ≤T}Φ(θ_γ)`, `Σ_{T<γ'≤X}[2Φ(θ_γ) - Φ(θ_γ + θ_{γ'})]` and the
far tail `Σ_{γ'>X} B(θ_γ, θ_{γ'})`. -/
def TinyTerms {κ : Type*} [DecidableEq κ] (w : ℕ → ℕ → ℝ) (θ : κ → ℝ) (lam : ℕ → ℝ)
    (Z : ℝ → Finset κ) (Xs : ℕ → Finset κ) (η : ℝ) : Prop :=
  (∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 → Z T ⊆ Xs N) ∧
    ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 →
      2 * (Z T).card * ∑ j ∈ Z T, Phi N (w N) (θ j) +
        ∑ j ∈ Z T, (∑ k ∈ Xs N \ Z T, (2 * Phi N (w N) (θ j) - Phi N (w N) (θ j + θ k)) +
          tail N (w N) θ lam (Xs N) (θ j)) ≤
        η / 2 * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π))

/-- **Theorem 1** from the lemmas of proofs.md: (H_Li) Li's formula, (H_N) RvM, (H7) `τ_T`
asymptotics, and for every `η > 0` weights with (H_w) `GoodWeights`, (H_tiny) `TinyTerms`,
(H4) `TrendPt`. -/
theorem thm1_of_lemmas {κ : Type*} [DecidableEq κ] (θ : κ → ℝ) (lam t : ℕ → ℝ)
    (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)) (ht0 : t 0 = 0)
    (Z : ℝ → Finset κ) (hN : RvM Z) (h7 : TauT θ Z)
    (hW : ∀ η > 0, ∃ w : ℕ → ℕ → ℝ, ∃ Xs : ℕ → Finset κ,
      GoodWeights w η ∧ TinyTerms w θ lam Z Xs η ∧ TrendPt w θ t Z) :
    Thm1Bound (fun n => lam n - t n) :=
  thm1 θ lam t hlam ht0 Z hN fun η hη => by
    obtain ⟨w, Xs, hg, ⟨hsub, htiny⟩, h4⟩ := hW η hη
    exact ⟨w, hg, trendTail_of w θ lam t Z η hg.2.2.1 Xs hsub htiny
      (htrend_of w θ t Z hN h4 h7 η hη)⟩

end LiRate

module

public import Riemann.Deep.LeanLi.Tails
public import Riemann.Deep.LeanLi.TrendFinal
public import Riemann.Deep.LeanLi.TrendPt2

/-!
# LeanLi: (H_tiny) for the zeros of ζ (proofs.md Theorem 1, step 5)
-/

@[expose] public section

open Real Finset Filter

namespace LeanLi

open Classical in
/-- **(H_tiny)** with `X = N⁸`. -/
theorem tiny_terms (L : LitCore) (hRH : RiemannHypothesis) (w : ℕ → ℕ → ℝ) (C : ℝ)
    (hgood : ∀ N : ℕ, 1 ≤ N → PhiGood N (w N) C) (η : ℝ) (hη : 0 < η) :
    LR.TinyTerms w thK lam ZK (fun N => ZK ((N : ℝ) ^ 8)) η := by
  have hπ := pi_gt_three
  have hπ0 := pi_pos
  constructor
  · filter_upwards [eventually_ge_atTop 1] with N hN
    intro T _ hT2 k hk
    rw [mem_ZK] at hk ⊢
    have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
    have hsN : √(N : ℝ) ≤ N := by
      rw [Real.sqrt_le_left (by positivity)]; nlinarith
    have : (N : ℝ) ≤ (N : ℝ) ^ 8 := le_self_pow₀ hN1 (by norm_num)
    have : 0 ≤ √(N : ℝ) := Real.sqrt_nonneg _
    linarith
  obtain ⟨Kt, hKt0, hKt⟩ := tsum_tail_le L
  obtain ⟨T₁c, hT₁c4, hcnt⟩ := count_le L
  have hC0 : 0 ≤ C := (hgood 1 le_rfl).C_nonneg le_rfl
  set Kti : ℝ := 18 * (C * 2 ^ 30) + 2 * Kt + 40 with hKti
  have hKti0 : 0 ≤ Kti := by positivity
  set T₁ : ℝ := max (max ((2 * π) ^ 2) (exp 2)) (max T₁c (12 * π * Kti / η + 1)) with hT₁
  have hT₁pos : 0 < T₁ := lt_of_lt_of_le (exp_pos 2) (le_trans (le_max_right _ _) (le_max_left _ _))
  filter_upwards [LR.eventually_basic 1 (1 / T₁) (2 / 3) (one_div_pos.2 hT₁pos) (by norm_num),
    eventually_ge_atTop (max 2 ⌈T₁c⌉₊)] with N hNb hN0
  obtain ⟨hl1, -, hsq⟩ := hNb
  have hN2 : 2 ≤ N := le_trans (le_max_left _ _) hN0
  have hN1 : 1 ≤ N := by omega
  have hNc : T₁c ≤ N := le_trans (Nat.le_ceil _) (by exact_mod_cast le_trans (le_max_right _ _) hN0)
  intro T hT1 hT2
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  have hN1' : (1 : ℝ) ≤ N := by exact_mod_cast hN1
  have hℓ0 : 0 < Real.log N := by linarith
  have hsx : 0 < √(N : ℝ) := Real.sqrt_pos.2 hNpos
  have hTbig : T₁ ≤ T := by
    have h : T₁ * Real.log N ≤ √N := by
      have := mul_le_mul_of_nonneg_left hsq hT₁pos.le
      have e : T₁ * (1 / T₁ * √(N : ℝ)) = √N := by field_simp
      linarith
    have : T₁ ≤ √N / Real.log N := by rw [le_div_iff₀ hℓ0]; linarith
    linarith
  have hT2π : (2 * π) ^ 2 ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hTbig
  have hTe : exp 2 ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hTbig
  have hTc : T₁c ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hTbig
  have hTK : 12 * π * Kti / η + 1 ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hTbig
  have hTpos : 0 < T := lt_of_lt_of_le (exp_pos 2) hTe
  have hg := hgood N hN1
  set τ := tinyB C N with hτdef
  have hτ0 : 0 ≤ τ := tinyB_nonneg hg hN1
  -- facts on the zeros up to T
  have hord : ∀ j ∈ ZK T, 0 < ordK j ∧ ordK j ≤ T := fun j hj => ⟨j.1.2.2, mem_ZK.1 hj⟩
  have hsN : √(N : ℝ) ≤ N := by rw [Real.sqrt_le_left (by positivity)]; nlinarith
  have hθ : ∀ j ∈ ZK T, 1 / √(N : ℝ) ≤ thK j ∧ thK j < π := by
    intro j hj
    obtain ⟨h0, h1⟩ := hord j hj
    exact ⟨liTheta_ge_inv_sqrt hN1 h0 (le_trans h1 hT2), liTheta_lt_pi h0⟩
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN1
  have hΦj : ∀ j ∈ ZK T, LR.Phi N (w N) (thK j) ≤ τ := by
    intro j hj
    obtain ⟨h1, h2⟩ := hθ j hj
    have hp : 0 < thK j := lt_of_lt_of_le (by positivity) h1
    exact hg.le_tiny hN1 (by rw [abs_of_pos hp]; linarith) (by rw [abs_of_pos hp]; linarith)
  -- counts
  have hlogT : Real.log T ≤ Real.log N := by
    apply Real.log_le_log hTpos
    have : 0 ≤ √(N : ℝ) := Real.sqrt_nonneg _
    linarith
  have hcardZ : ((ZK T).card : ℝ) ≤ √(N : ℝ) * Real.log N := by
    rw [card_ZK]
    have h := hcnt T hTc
    have hlT : 0 ≤ Real.log T := Real.log_nonneg (by linarith [add_one_le_exp (2 : ℝ)])
    calc (zeroCount T : ℝ) ≤ T * Real.log T := h
      _ ≤ √(N : ℝ) * Real.log N := by
          apply mul_le_mul (by linarith [Real.sqrt_nonneg (N : ℝ)]) hlogT hlT (Real.sqrt_nonneg _)
  set X : ℝ := (N : ℝ) ^ 8 with hX
  have hXpos : 0 < X := by positivity
  have hXc : T₁c ≤ X := le_trans hNc (le_self_pow₀ hN1' (by norm_num))
  have hcardX : ((ZK X).card : ℝ) ≤ 8 * (N : ℝ) ^ 8 * Real.log N := by
    rw [card_ZK]
    have h := hcnt X hXc
    rw [hX, Real.log_pow] at h; push_cast at h
    linarith
  have hlN : Real.log N ≤ N := log_le_self' hNpos
  -- S1
  have hS1 : 2 * (ZK T).card * ∑ j ∈ ZK T, LR.Phi N (w N) (thK j) ≤ 2 * (C * 2 ^ 30) := by
    have h1 : ∑ j ∈ ZK T, LR.Phi N (w N) (thK j) ≤ (ZK T).card * τ := by
      have := Finset.sum_le_sum hΦj
      rwa [Finset.sum_const, nsmul_eq_mul] at this
    have hc0 : (0 : ℝ) ≤ (ZK T).card := Nat.cast_nonneg _
    have h2 : ((ZK T).card : ℝ) ^ 2 ≤ (N : ℝ) ^ 4 := by
      have : ((ZK T).card : ℝ) ≤ (N : ℝ) * N := by
        calc ((ZK T).card : ℝ) ≤ √(N : ℝ) * Real.log N := hcardZ
          _ ≤ N * N := mul_le_mul hsN hlN hℓ0.le hNpos.le
      have := pow_le_pow_left₀ hc0 this 2
      calc ((ZK T).card : ℝ) ^ 2 ≤ ((N : ℝ) * N) ^ 2 := this
        _ = (N : ℝ) ^ 4 := by ring
    have h3 := tinyB_mul_pow_le hC0 hN1 (k := 4) (by norm_num)
    calc 2 * (ZK T).card * ∑ j ∈ ZK T, LR.Phi N (w N) (thK j) ≤ 2 * (ZK T).card * ((ZK T).card * τ) :=
          mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = 2 * (((ZK T).card : ℝ) ^ 2 * τ) := by ring
      _ ≤ 2 * ((N : ℝ) ^ 4 * τ) := by gcongr
      _ = 2 * (τ * (N : ℝ) ^ 4) := by ring
      _ ≤ 2 * (C * 2 ^ 30) := by gcongr
  -- S2
  have hS2 : ∑ j ∈ ZK T, ∑ k ∈ ZK X \ ZK T,
      (2 * LR.Phi N (w N) (thK j) - LR.Phi N (w N) (thK j + thK k)) ≤ 16 * (C * 2 ^ 30) := by
    have h1 : ∀ j ∈ ZK T, ∑ k ∈ ZK X \ ZK T,
        (2 * LR.Phi N (w N) (thK j) - LR.Phi N (w N) (thK j + thK k)) ≤ (ZK X).card * (2 * τ) := by
      intro j hj
      have : ∀ k ∈ ZK X \ ZK T, 2 * LR.Phi N (w N) (thK j) - LR.Phi N (w N) (thK j + thK k) ≤ 2 * τ :=
        fun k _ => by linarith [hg.phi_nonneg (thK j + thK k), hΦj j hj]
      calc ∑ k ∈ ZK X \ ZK T, (2 * LR.Phi N (w N) (thK j) - LR.Phi N (w N) (thK j + thK k))
          ≤ ∑ k ∈ ZK X \ ZK T, 2 * τ := Finset.sum_le_sum this
        _ = (ZK X \ ZK T).card * (2 * τ) := by rw [Finset.sum_const, nsmul_eq_mul]
        _ ≤ (ZK X).card * (2 * τ) := by
            gcongr; exact Finset.sdiff_subset
    have h2 := Finset.sum_le_sum h1
    rw [Finset.sum_const, nsmul_eq_mul] at h2
    refine h2.trans ?_
    have hc0 : (0 : ℝ) ≤ (ZK T).card := Nat.cast_nonneg _
    have h11 := tinyB_mul_pow_le hC0 hN1 (k := 11) (by norm_num)
    have hb1 : ((ZK T).card : ℝ) ≤ (N : ℝ) * N := by
      calc ((ZK T).card : ℝ) ≤ √(N : ℝ) * Real.log N := hcardZ
        _ ≤ N * N := mul_le_mul hsN hlN hℓ0.le hNpos.le
    have hb2 : ((ZK X).card : ℝ) ≤ 8 * (N : ℝ) ^ 9 := by
      calc ((ZK X).card : ℝ) ≤ 8 * (N : ℝ) ^ 8 * Real.log N := hcardX
        _ ≤ 8 * (N : ℝ) ^ 8 * N := by gcongr
        _ = 8 * (N : ℝ) ^ 9 := by ring
    calc ((ZK T).card : ℝ) * ((ZK X).card * (2 * τ)) ≤ ((N : ℝ) * N) * ((8 * (N : ℝ) ^ 9) * (2 * τ)) := by
          gcongr
      _ = 16 * (τ * (N : ℝ) ^ 11) := by ring
      _ ≤ 16 * (C * 2 ^ 30) := by gcongr
  -- S3
  have hlam := hasSum_lam L hRH
  have hS3 : ∑ j ∈ ZK T, LR.tail N (w N) thK lam (ZK X) (thK j) ≤ 2 * Kt + 40 := by
    have htail : ∀ j ∈ ZK T, LR.tail N (w N) thK lam (ZK X) (thK j) ≤
        2 * (N : ℝ) ^ 3 * (Kt / X ^ 2 + 10 / √X) := by
      intro j _
      rw [LR.tail_eq_tsum_compl N (w N) thK lam hlam]
      have hsB : Summable (fun k : ↑((ZK X : Set ZeroIdx)ᶜ) => LR.B N (w N) (thK j) (thK k)) :=
        (LR.hasSum_B N (w N) thK lam hlam (thK j)).summable.subtype _
      have hsI : Summable (fun k : ↑((ZK X : Set ZeroIdx)ᶜ) => 2 * (N : ℝ) ^ 3 * (1 / ordK k ^ 2)) :=
        ((summable_inv_ordK_sq L).subtype _).mul_left _
      calc ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), LR.B N (w N) (thK j) (thK k)
          ≤ ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), 2 * (N : ℝ) ^ 3 * (1 / ordK k ^ 2) := by
            refine hsB.tsum_le_tsum (fun k => ?_) hsI
            have hB := abs_B_le_quad N (w N) hg.nonneg hg.le_one (thK j) (thK k)
            have hk0 : 0 < ordK (k : ZeroIdx) := (k : ZeroIdx).1.2.2
            have hθk := liTheta_le_inv hk0
            have hθk0 := liTheta_pos hk0
            have : thK (k : ZeroIdx) ^ 2 ≤ 1 / ordK (k : ZeroIdx) ^ 2 := by
              have := pow_le_pow_left₀ hθk0.le hθk 2
              rwa [div_pow, one_pow] at this
            have h3 : 2 * (N : ℝ) ^ 3 * thK (k : ZeroIdx) ^ 2 ≤ 2 * (N : ℝ) ^ 3 * (1 / ordK (k : ZeroIdx) ^ 2) :=
              mul_le_mul_of_nonneg_left this (by positivity)
            linarith [le_abs_self (LR.B N (w N) (thK j) (thK k))]
        _ = 2 * (N : ℝ) ^ 3 * ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), 1 / ordK k ^ 2 := tsum_mul_left
        _ ≤ 2 * (N : ℝ) ^ 3 * (Kt / X ^ 2 + 10 / √X) :=
            mul_le_mul_of_nonneg_left (hKt X hXpos) (by positivity)
    have h2 := Finset.sum_le_sum htail
    rw [Finset.sum_const, nsmul_eq_mul] at h2
    refine h2.trans ?_
    have hsX : √X = (N : ℝ) ^ 4 := by
      rw [hX, show (8 : ℕ) = 4 * 2 by norm_num, pow_mul, Real.sqrt_sq (by positivity)]
    rw [hsX]
    have hc0 : (0 : ℝ) ≤ (ZK T).card := Nat.cast_nonneg _
    have hl2 : Real.log N ≤ 2 * √(N : ℝ) := log_le_two_sqrt hNpos
    have hss : √(N : ℝ) * √(N : ℝ) = N := Real.mul_self_sqrt hNpos.le
    have e : 2 * (N : ℝ) ^ 3 * (Kt / X ^ 2 + 10 / (N : ℝ) ^ 4) = 2 * Kt / (N : ℝ) ^ 13 + 20 / N := by
      rw [hX]; field_simp; ring
    rw [e]
    have t1 : ((ZK T).card : ℝ) * (2 * Kt / (N : ℝ) ^ 13) ≤ 2 * Kt := by
      have hb : ((ZK T).card : ℝ) ≤ (N : ℝ) ^ 13 := by
        calc ((ZK T).card : ℝ) ≤ √(N : ℝ) * Real.log N := hcardZ
          _ ≤ N * N := mul_le_mul hsN hlN hℓ0.le hNpos.le
          _ = (N : ℝ) ^ 2 := by ring
          _ ≤ (N : ℝ) ^ 13 := pow_le_pow_right₀ hN1' (by norm_num)
      calc ((ZK T).card : ℝ) * (2 * Kt / (N : ℝ) ^ 13) ≤ (N : ℝ) ^ 13 * (2 * Kt / (N : ℝ) ^ 13) :=
            mul_le_mul_of_nonneg_right hb (by positivity)
        _ = 2 * Kt := by field_simp
    have t2 : ((ZK T).card : ℝ) * (20 / N) ≤ 40 := by
      have hb : ((ZK T).card : ℝ) ≤ 2 * N := by
        calc ((ZK T).card : ℝ) ≤ √(N : ℝ) * Real.log N := hcardZ
          _ ≤ √(N : ℝ) * (2 * √(N : ℝ)) := mul_le_mul_of_nonneg_left hl2 (Real.sqrt_nonneg _)
          _ = 2 * (√(N : ℝ) * √(N : ℝ)) := by ring
          _ = 2 * N := by rw [hss]
      calc ((ZK T).card : ℝ) * (20 / N) ≤ (2 * N) * (20 / N) :=
            mul_le_mul_of_nonneg_right hb (by positivity)
        _ = 40 := by field_simp; ring
    rw [mul_add]
    linarith
  -- RHS
  have hL : 1 ≤ Real.log (T / (2 * π)) := by
    rw [Real.log_div hTpos.ne' (by positivity)]
    have h1 : 2 ≤ Real.log T := by rw [← Real.log_exp 2]; exact Real.log_le_log (exp_pos 2) hTe
    have h2 := Real.log_le_log (by positivity) hT2π
    rw [Real.log_pow] at h2; push_cast at h2
    linarith
  have hRHS : Kti ≤ η / 2 * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π)) := by
    have hT3 : T ≤ T ^ 3 * Real.log (T / (2 * π)) ^ 2 := by
      have h1 : 1 ≤ Real.log (T / (2 * π)) ^ 2 := one_le_pow₀ hL
      have hT1' : 1 ≤ T := by linarith [add_one_le_exp (2 : ℝ)]
      have h2 : T ≤ T ^ 3 := by
        have := le_mul_of_one_le_right hTpos.le (one_le_pow₀ (n := 2) hT1')
        calc T ≤ T * T ^ 2 := this
          _ = T ^ 3 := by ring
      exact h2.trans (le_mul_of_one_le_right (by positivity) h1)
    have hK : 12 * π * Kti ≤ η * T := by
      have : 12 * π * Kti / η ≤ T := by linarith
      rw [div_le_iff₀ hη] at this; linarith
    rw [show η / 2 * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π)) =
      η * (T ^ 3 * Real.log (T / (2 * π)) ^ 2) / (12 * π) by ring, le_div_iff₀ (by positivity)]
    have := mul_le_mul_of_nonneg_left hT3 hη.le
    linarith
  have hsplit : ∑ j ∈ ZK T, (∑ k ∈ ZK X \ ZK T,
      (2 * LR.Phi N (w N) (thK j) - LR.Phi N (w N) (thK j + thK k)) +
        LR.tail N (w N) thK lam (ZK X) (thK j)) =
      (∑ j ∈ ZK T, ∑ k ∈ ZK X \ ZK T,
        (2 * LR.Phi N (w N) (thK j) - LR.Phi N (w N) (thK j + thK k))) +
      ∑ j ∈ ZK T, LR.tail N (w N) thK lam (ZK X) (thK j) := Finset.sum_add_distrib
  show 2 * (ZK T).card * ∑ j ∈ ZK T, LR.Phi N (w N) (thK j) +
      ∑ j ∈ ZK T, (∑ k ∈ ZK X \ ZK T,
        (2 * LR.Phi N (w N) (thK j) - LR.Phi N (w N) (thK j + thK k)) +
        LR.tail N (w N) thK lam (ZK X) (thK j)) ≤
      η / 2 * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π))
  rw [hsplit]
  linarith

end LeanLi

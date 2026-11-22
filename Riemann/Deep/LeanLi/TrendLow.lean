module

public import Riemann.Deep.LeanLi.TrendI

/-!
# LeanLi: lower bound for the trend transform of `I(n)` (proofs.md Theorem 1, step 4)

For weights with `Φ ≥ 0`, `∫Φ = 2π` and `Φ(x)(1 + N|x|)³⁰ ≤ CN`:
`τ_I(a) ≥ log(2πa)/a² - O(√N log N)` for `N^{-1/2} ≤ a ≤ π - δ₀`.
-/

@[expose] public section

open Real Finset MeasureTheory Set

namespace LeanLi

/-- Hypotheses on the weights used for the trend estimate. -/
structure PhiGood (N : ℕ) (w : ℕ → ℝ) (C : ℝ) : Prop where
  w0 : w 0 = 1
  nonneg : ∀ n, 0 ≤ w n
  le_one : ∀ n, w n ≤ 1
  phi_nonneg : ∀ x, 0 ≤ LR.Phi N w x
  decay : ∀ x, |x| ≤ π → LR.Phi N w x * (1 + N * |x|) ^ 30 ≤ C * N

/-- The tiny bound `C 2³⁰/N¹⁴`. -/
noncomputable def tinyB (C : ℝ) (N : ℕ) : ℝ := C * 2 ^ 30 / (N : ℝ) ^ 14

section

variable {N : ℕ} {w : ℕ → ℝ} {C : ℝ}

lemma PhiGood.C_nonneg (hg : PhiGood N w C) (hN : 1 ≤ N) : 0 ≤ C := by
  have h := hg.decay 0 (by simp [pi_pos.le])
  have h0 := hg.phi_nonneg 0
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have : 0 ≤ C * N := le_trans (mul_nonneg h0 (by positivity)) h
  exact nonneg_of_mul_nonneg_left this hNpos

lemma tinyB_nonneg (hg : PhiGood N w C) (hN : 1 ≤ N) : 0 ≤ tinyB C N := by
  unfold tinyB; have := hg.C_nonneg hN; positivity

lemma PhiGood.le_tiny (hg : PhiGood N w C) (hN : 1 ≤ N) {x : ℝ}
    (hx : 1 / (2 * √(N : ℝ)) ≤ |x|) (hx' : |x| ≤ π) : LR.Phi N w x ≤ tinyB C N :=
  phi_far hN hg.decay hg.phi_nonneg hx hx'

/-- `Φ ≤ τ` on `[ε, 2π - ε]`, `ε = 1/(2√N)`. -/
lemma PhiGood.le_tiny_pos (hg : PhiGood N w C) (hN : 1 ≤ N) {x : ℝ}
    (hx : 1 / (2 * √(N : ℝ)) ≤ x) (hx' : x ≤ 2 * π - 1 / (2 * √(N : ℝ))) :
    LR.Phi N w x ≤ tinyB C N := by
  have hε : 0 < 1 / (2 * √(N : ℝ)) := by
    have : (0 : ℝ) < N := by exact_mod_cast hN
    positivity
  rcases le_or_gt x π with h | h
  · exact hg.le_tiny hN (by rw [abs_of_pos (by linarith)]; exact hx) (by rw [abs_of_pos (by linarith)]; exact h)
  · rw [← phi_two_pi_sub]
    exact hg.le_tiny hN (by rw [abs_of_pos (by linarith)]; linarith) (by rw [abs_of_pos (by linarith)]; linarith)

/-- `Φ ≤ τ` on `[-π, -ε]`. -/
lemma PhiGood.le_tiny_neg (hg : PhiGood N w C) (hN : 1 ≤ N) {x : ℝ}
    (hx : x ≤ -(1 / (2 * √(N : ℝ)))) (hx' : -π ≤ x) : LR.Phi N w x ≤ tinyB C N := by
  have hε : 0 < 1 / (2 * √(N : ℝ)) := by
    have : (0 : ℝ) < N := by exact_mod_cast hN
    positivity
  rw [← phi_neg]
  have h1 : 0 < -x := by linarith
  exact hg.le_tiny hN (by rw [abs_of_pos h1]; linarith) (by rw [abs_of_pos h1]; linarith)

end

/-- `|B(a, v)| ≤ 2N³v²` for `0 ≤ w ≤ 1`. -/
lemma abs_B_le_quad (N : ℕ) (w : ℕ → ℝ) (hw0 : ∀ n, 0 ≤ w n) (hw1 : ∀ n, w n ≤ 1) (a v : ℝ) :
    |LR.B N w a v| ≤ 2 * (N : ℝ) ^ 3 * v ^ 2 := by
  rw [B_eq_sum, abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 2)]
  have hterm : ∀ n ∈ Finset.Icc 1 N, |w n * cos (n * a) * (2 * (1 - cos (n * v)))| ≤ (N : ℝ) ^ 2 * v ^ 2 := by
    intro n hn
    have hnN : (n : ℝ) ≤ N := by exact_mod_cast (Finset.mem_Icc.1 hn).2
    have hc := one_sub_sq_div_two_le_cos (x := n * v)
    have hc1 := cos_le_one (n * v)
    rw [abs_mul, abs_mul, abs_of_nonneg (hw0 n), abs_of_nonneg (by linarith : (0 : ℝ) ≤ 2 * (1 - cos (n * v)))]
    have h1 : 2 * (1 - cos (n * v)) ≤ (n * v) ^ 2 := by linarith
    have h2 : (n * v) ^ 2 ≤ (N : ℝ) ^ 2 * v ^ 2 := by
      rw [mul_pow]; exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (Nat.cast_nonneg _) hnN 2) (sq_nonneg _)
    calc w n * |cos (n * a)| * (2 * (1 - cos (n * v))) ≤ 1 * 1 * (2 * (1 - cos (n * v))) := by
          gcongr
          · exact hw1 n
          · exact abs_cos_le_one _
      _ ≤ (N : ℝ) ^ 2 * v ^ 2 := by linarith
  calc 2 * |∑ n ∈ Finset.Icc 1 N, w n * cos (n * a) * (2 * (1 - cos (n * v)))|
      ≤ 2 * ∑ n ∈ Finset.Icc 1 N, |w n * cos (n * a) * (2 * (1 - cos (n * v)))| := by
        gcongr; exact Finset.abs_sum_le_sum_abs _ _
    _ ≤ 2 * ∑ n ∈ Finset.Icc 1 N, (N : ℝ) ^ 2 * v ^ 2 := by gcongr with n hn; exact hterm n hn
    _ = 2 * (N : ℝ) ^ 3 * v ^ 2 := by
        rw [Finset.sum_const, Nat.card_Icc, nsmul_eq_mul]; push_cast; ring

/-- `L₁ = ∫₀^π |ℓ|`. -/
noncomputable def L1 : ℝ := ∫ v in (0 : ℝ)..π, |ell v|

lemma L1_nonneg : 0 ≤ L1 := intervalIntegral.integral_nonneg pi_pos.le fun _ _ => abs_nonneg _

lemma integral_abs_ell_le_L1 {c d : ℝ} (hc : 0 ≤ c) (hcd : c ≤ d) (hd : d ≤ π) :
    ∫ v in c..d, |ell v| ≤ L1 :=
  intervalIntegral.integral_mono_interval hc hcd hd
    (Filter.Eventually.of_forall fun _ => abs_nonneg _) (intervalIntegrable_ell 0 π).abs

/-- `∫₀^δ |ℓ| ≤ δ(2 + log π + log 4π - log δ)` for `0 < δ ≤ 1`. -/
lemma integral_abs_ell_small {δ : ℝ} (hδ : 0 < δ) (hδ1 : δ ≤ 1) :
    ∫ v in (0 : ℝ)..δ, |ell v| ≤ δ * (2 + log π + log (4 * π) - log δ) := by
  have hπ := pi_gt_three
  have hle : ∫ v in (0 : ℝ)..δ, |ell v| ≤ ∫ v in (0 : ℝ)..δ, (1 + log π + log (4 * π) - log v) := by
    apply intervalIntegral.integral_mono_on_of_le_Ioo hδ.le (intervalIntegrable_ell 0 δ).abs
    · exact (intervalIntegrable_const.sub intervalIntegral.intervalIntegrable_log')
    · intro v hv
      have h := abs_ell_le hv.1 (by linarith [hv.2])
      have h1 : log (π / v) = log π - log v := log_div (by positivity) hv.1.ne'
      have h2 : log (π / (π - v)) ≤ 1 := by
        have hpv : 0 < π - v := by linarith [hv.2]
        rw [log_div (by positivity) hpv.ne']
        have : log (π - 1) ≤ log (π - v) := log_le_log (by linarith) (by linarith [hv.2])
        have : log π - log (π - 1) ≤ 1 := by
          rw [← log_div (by positivity) (by linarith)]
          have := log_le_sub_one_of_pos (show 0 < π / (π - 1) by apply div_pos <;> linarith)
          have hp1 : π - 1 ≠ 0 := by linarith
          have : π / (π - 1) - 1 = 1 / (π - 1) := by field_simp; ring
          have : 1 / (π - 1) ≤ 1 := by rw [div_le_one (by linarith)]; linarith
          linarith
        linarith
      linarith
  refine hle.trans (le_of_eq ?_)
  rw [intervalIntegral.integral_sub intervalIntegrable_const intervalIntegral.intervalIntegrable_log',
    intervalIntegral.integral_const, integral_log]
  simp; ring

lemma eps_facts {N : ℕ} (hN : 1 ≤ N) :
    0 < 1 / (2 * √(N : ℝ)) ∧ 1 / (2 * √(N : ℝ)) ≤ 1 / 2 ∧ 1 / (2 * √(N : ℝ)) ≤ 1 / √(N : ℝ) ∧
      1 / (2 * √(N : ℝ)) = 1 / √(N : ℝ) / 2 := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hs1 : 1 ≤ √(N : ℝ) := by rw [Real.one_le_sqrt]; exact_mod_cast hN
  have hs : 0 < √(N : ℝ) := by linarith
  refine ⟨by positivity, ?_, ?_, by rw [div_div, mul_comm]⟩
  · rw [div_le_div_iff₀ (by positivity) (by norm_num)]; linarith
  · rw [div_le_div_iff₀ (by positivity) hs]; linarith

/-- `|B(a, v)| ≤ 2τ` for `0 ≤ v ≤ a/2`, `N^{-1/2} ≤ a ≤ π`. -/
lemma abs_B_le_tiny {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 1 ≤ N) {a v : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π) (hv : 0 ≤ v) (hv' : v ≤ a / 2) :
    |LR.B N w a v| ≤ 2 * tinyB C N := by
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN
  have hπ := pi_gt_three
  have h1 : LR.Phi N w a ≤ tinyB C N := hg.le_tiny_pos hN (by linarith) (by linarith)
  have h2 : LR.Phi N w (a - v) ≤ tinyB C N := hg.le_tiny_pos hN (by linarith) (by linarith)
  have h3 : LR.Phi N w (a + v) ≤ tinyB C N := hg.le_tiny_pos hN (by linarith) (by linarith)
  have p1 := hg.phi_nonneg a
  have p2 := hg.phi_nonneg (a - v)
  have p3 := hg.phi_nonneg (a + v)
  unfold LR.B
  rw [abs_le]; constructor <;> linarith

/-- **Near `v = 0`**: `∫₀^{a/2} B(a,v) a₀(v) dv ≥ -(π/4)(2 + log π + log 4π + 3 log N) - (π/4)τN⁶L₁`. -/
lemma near_zero_bound {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 2 ≤ N) {a : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π) :
    -((π / 4) * (2 + log π + log (4 * π) + 3 * log N) + π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1) ≤
      ∫ v in (0 : ℝ)..(a / 2), LR.B N w a v * a0 v := by
  have hN1 : 1 ≤ N := by omega
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN1
  have hπ := pi_gt_three
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  have hN2 : (2 : ℝ) ≤ N := by exact_mod_cast hN
  have ha0 : 0 < a := lt_of_lt_of_le (by positivity) ha
  set δ : ℝ := 1 / (N : ℝ) ^ 3 with hδ
  have hδpos : 0 < δ := by positivity
  have hδ1 : δ ≤ 1 := by
    rw [hδ, div_le_one (by positivity)]; exact one_le_pow₀ (by linarith)
  have hδa : δ ≤ a / 2 := by
    have hsN : √(N : ℝ) ≤ N := by
      rw [Real.sqrt_le_left]; nlinarith
      exact hNpos.le
    have : δ ≤ 1 / (2 * √(N : ℝ)) := by
      rw [hδ, div_le_div_iff₀ (by positivity) (by positivity)]
      nlinarith [Real.sqrt_nonneg (N : ℝ), mul_nonneg hNpos.le (by nlinarith : (0:ℝ) ≤ (N:ℝ)^2 - 2)]
    linarith
  have hτ := tinyB_nonneg hg hN1
  have hI := intervalIntegrable_B_a0 N w a
  have hI1 : IntervalIntegrable (fun v => LR.B N w a v * a0 v) volume 0 δ :=
    hI.mono_set (by rw [uIcc_of_le hδpos.le, uIcc_of_le pi_pos.le]; exact Icc_subset_Icc le_rfl (by linarith))
  have hI2 : IntervalIntegrable (fun v => LR.B N w a v * a0 v) volume δ (a / 2) :=
    hI.mono_set (by rw [uIcc_of_le hδa, uIcc_of_le pi_pos.le]; exact Icc_subset_Icc hδpos.le (by linarith))
  rw [← intervalIntegral.integral_add_adjacent_intervals hI1 hI2]
  -- piece 1
  have hp1 : -(π / 4 * (2 + log π + log (4 * π) + 3 * log N)) ≤
      ∫ v in (0 : ℝ)..δ, LR.B N w a v * a0 v := by
    have hle : ∫ v in (0 : ℝ)..δ, -(π * (N : ℝ) ^ 3 / 4 * |ell v|) ≤
        ∫ v in (0 : ℝ)..δ, LR.B N w a v * a0 v := by
      apply intervalIntegral.integral_mono_on_of_le_Ioo hδpos.le
        (((intervalIntegrable_ell 0 δ).abs.const_mul _).neg) hI1
      intro v hv
      show -(π * (N : ℝ) ^ 3 / 4 * |ell v|) ≤ LR.B N w a v * a0 v
      have hB := abs_B_le_quad N w hg.nonneg hg.le_one a v
      have hA := abs_a0_le hv.1 (by linarith [hv.2])
      have hv0 := hv.1
      have : |LR.B N w a v * a0 v| ≤ π * (N : ℝ) ^ 3 / 4 * |ell v| := by
        rw [abs_mul]
        calc |LR.B N w a v| * |a0 v| ≤ (2 * (N : ℝ) ^ 3 * v ^ 2) * (π * |ell v| / (8 * v ^ 2)) :=
              mul_le_mul hB hA (abs_nonneg _) (by positivity)
          _ = π * (N : ℝ) ^ 3 / 4 * |ell v| := by field_simp; ring
      linarith [neg_abs_le (LR.B N w a v * a0 v)]
    have hsm := integral_abs_ell_small hδpos hδ1
    have hlogδ : log δ = -(3 * log N) := by
      rw [hδ, one_div, log_inv, log_pow]; push_cast; ring
    rw [intervalIntegral.integral_neg, intervalIntegral.integral_const_mul] at hle
    have hN3 : (N : ℝ) ^ 3 * δ = 1 := by rw [hδ]; field_simp
    have : π * (N : ℝ) ^ 3 / 4 * ∫ v in (0 : ℝ)..δ, |ell v| ≤
        π / 4 * (2 + log π + log (4 * π) + 3 * log N) := by
      calc π * (N : ℝ) ^ 3 / 4 * ∫ v in (0 : ℝ)..δ, |ell v| ≤
            π * (N : ℝ) ^ 3 / 4 * (δ * (2 + log π + log (4 * π) - log δ)) :=
            mul_le_mul_of_nonneg_left hsm (by positivity)
        _ = π / 4 * ((N : ℝ) ^ 3 * δ) * (2 + log π + log (4 * π) - log δ) := by ring
        _ = π / 4 * (2 + log π + log (4 * π) + 3 * log N) := by rw [hN3, hlogδ]; ring
    linarith
  -- piece 2
  have hp2 : -(π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1) ≤
      ∫ v in δ..(a / 2), LR.B N w a v * a0 v := by
    have hle : ∫ v in δ..(a / 2), -(π * tinyB C N * (N : ℝ) ^ 6 / 4 * |ell v|) ≤
        ∫ v in δ..(a / 2), LR.B N w a v * a0 v := by
      apply intervalIntegral.integral_mono_on_of_le_Ioo hδa
        (((intervalIntegrable_ell δ (a / 2)).abs.const_mul _).neg) hI2
      intro v hv
      show -(π * tinyB C N * (N : ℝ) ^ 6 / 4 * |ell v|) ≤ LR.B N w a v * a0 v
      have hv0 : 0 < v := lt_trans hδpos hv.1
      have hB := abs_B_le_tiny hg hN1 ha ha' hv0.le hv.2.le
      have hA := abs_a0_le hv0 (by linarith [hv.2])
      have hvδ : δ ^ 2 ≤ v ^ 2 := pow_le_pow_left₀ hδpos.le hv.1.le 2
      have hinv : 1 / v ^ 2 ≤ (N : ℝ) ^ 6 := by
        rw [div_le_iff₀ (by positivity)]
        have : (N : ℝ) ^ 6 * δ ^ 2 = 1 := by rw [hδ]; field_simp
        have h6 : 0 ≤ (N : ℝ) ^ 6 := by positivity
        nlinarith [mul_le_mul_of_nonneg_left hvδ h6]
      have : |LR.B N w a v * a0 v| ≤ π * tinyB C N * (N : ℝ) ^ 6 / 4 * |ell v| := by
        rw [abs_mul]
        calc |LR.B N w a v| * |a0 v| ≤ (2 * tinyB C N) * (π * |ell v| / (8 * v ^ 2)) :=
              mul_le_mul hB hA (abs_nonneg _) (by positivity)
          _ = π * tinyB C N / 4 * |ell v| * (1 / v ^ 2) := by field_simp; ring
          _ ≤ π * tinyB C N / 4 * |ell v| * (N : ℝ) ^ 6 :=
              mul_le_mul_of_nonneg_left hinv (by positivity)
          _ = π * tinyB C N * (N : ℝ) ^ 6 / 4 * |ell v| := by ring
      linarith [neg_abs_le (LR.B N w a v * a0 v)]
    rw [intervalIntegral.integral_neg, intervalIntegral.integral_const_mul] at hle
    have hL := integral_abs_ell_le_L1 hδpos.le hδa (by linarith)
    have : π * tinyB C N * (N : ℝ) ^ 6 / 4 * ∫ v in δ..(a / 2), |ell v| ≤
        π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1 := mul_le_mul_of_nonneg_left hL (by positivity)
    linarith
  linarith

lemma intervalIntegrable_a0 {c d : ℝ} (hc : 0 < c) (hcd : c ≤ d) (hd : d ≤ π) :
    IntervalIntegrable a0 volume c d := by
  have hcont : ContinuousOn (fun v : ℝ => 1 / (8 * π * sin (v / 2) ^ 2)) (uIcc c d) := by
    rw [uIcc_of_le hcd]
    refine ContinuousOn.div continuousOn_const (Continuous.continuousOn (by fun_prop)) ?_
    intro v hv
    have := sin_half_pos (lt_of_lt_of_le hc hv.1) (le_trans hv.2 hd)
    positivity
  have := (intervalIntegrable_ell c d).mul_continuousOn hcont
  refine this.congr_uIoo fun v _ => ?_
  simp only [a0]; ring

/-- `|a₀(v)| ≤ (πN/2)|ℓ(v)|` for `v ≥ 1/(2√N)`. -/
lemma abs_a0_le_N {N : ℕ} (hN : 1 ≤ N) {v : ℝ} (hv : 1 / (2 * √(N : ℝ)) ≤ v) (hv' : v ≤ π) :
    |a0 v| ≤ π * N / 2 * |ell v| := by
  obtain ⟨he0, -, -, -⟩ := eps_facts hN
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hv0 : 0 < v := lt_of_lt_of_le he0 hv
  have h := abs_a0_le hv0 hv'
  have hsq : 1 / (4 * N) ≤ v ^ 2 := by
    have := pow_le_pow_left₀ he0.le hv 2
    rw [div_pow, mul_pow, Real.sq_sqrt hNpos.le] at this
    norm_num at this ⊢; linarith
  have hinv : 1 / v ^ 2 ≤ 4 * N := by
    rw [div_le_iff₀ (by positivity)]
    rw [div_le_iff₀ (by positivity)] at hsq
    linarith
  calc |a0 v| ≤ π * |ell v| / (8 * v ^ 2) := h
    _ = π * |ell v| / 8 * (1 / v ^ 2) := by field_simp
    _ ≤ π * |ell v| / 8 * (4 * N) := mul_le_mul_of_nonneg_left hinv (by positivity)
    _ = π * N / 2 * |ell v| := by ring

lemma integral_abs_a0_le {N : ℕ} (hN : 1 ≤ N) {c : ℝ} (hc : 1 / (2 * √(N : ℝ)) ≤ c) (hcπ : c ≤ π) :
    ∫ v in c..π, |a0 v| ≤ π * N / 2 * L1 := by
  obtain ⟨he0, -, -, -⟩ := eps_facts hN
  have hc0 : 0 < c := lt_of_lt_of_le he0 hc
  have h1 : ∫ v in c..π, |a0 v| ≤ ∫ v in c..π, π * N / 2 * |ell v| := by
    apply intervalIntegral.integral_mono_on_of_le_Ioo hcπ (intervalIntegrable_a0 hc0 hcπ le_rfl).abs
      ((intervalIntegrable_ell c π).abs.const_mul _)
    intro v hv
    exact abs_a0_le_N hN (by linarith [hv.1]) hv.2.le
  rw [intervalIntegral.integral_const_mul] at h1
  have h2 := integral_abs_ell_le_L1 hc0.le hcπ le_rfl
  have : π * N / 2 * ∫ v in c..π, |ell v| ≤ π * N / 2 * L1 :=
    mul_le_mul_of_nonneg_left h2 (by positivity)
  linarith

/-- **Middle terms**: `|∫_{a/2}^π (2Φ(a) - Φ(a+v)) a₀(v) dv| ≤ (3π/2) τ N L₁`. -/
lemma mid_bound {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 1 ≤ N) {a δ₀ : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π - δ₀) (hδ₀ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2) :
    |∫ v in (a / 2)..π, (2 * LR.Phi N w a - LR.Phi N w (a + v)) * a0 v| ≤
      3 * π / 2 * tinyB C N * N * L1 := by
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN
  have hπ := pi_gt_three
  have hτ := tinyB_nonneg hg hN
  have ha0 : 0 < a := lt_of_lt_of_le (by linarith) ha
  have hi : IntervalIntegrable (fun v => (2 * LR.Phi N w a - LR.Phi N w (a + v)) * a0 v) volume (a / 2) π :=
    (intervalIntegrable_a0 (by linarith) (by linarith) le_rfl).continuousOn_mul
      (Continuous.continuousOn (by unfold LR.Phi; fun_prop))
  refine (intervalIntegral.abs_integral_le_integral_abs (by linarith)).trans ?_
  have h1 : ∫ v in (a / 2)..π, |(2 * LR.Phi N w a - LR.Phi N w (a + v)) * a0 v| ≤
      ∫ v in (a / 2)..π, 3 * tinyB C N * |a0 v| := by
    apply intervalIntegral.integral_mono_on_of_le_Ioo (by linarith) hi.abs
      ((intervalIntegrable_a0 (by linarith) (by linarith) le_rfl).abs.const_mul _)
    intro v hv
    rw [abs_mul]
    apply mul_le_mul_of_nonneg_right _ (abs_nonneg _)
    have p1 := hg.phi_nonneg a
    have p2 := hg.phi_nonneg (a + v)
    have t1 : LR.Phi N w a ≤ tinyB C N := hg.le_tiny_pos hN (by linarith) (by linarith)
    have t2 : LR.Phi N w (a + v) ≤ tinyB C N := hg.le_tiny_pos hN (by linarith [hv.1]) (by linarith [hv.2])
    rw [abs_le]; constructor <;> linarith
  rw [intervalIntegral.integral_const_mul] at h1
  have h2 := integral_abs_a0_le hN (c := a / 2) (by linarith) (by linarith)
  have : 3 * tinyB C N * ∫ v in (a / 2)..π, |a0 v| ≤ 3 * tinyB C N * (π * N / 2 * L1) :=
    mul_le_mul_of_nonneg_left h2 (by positivity)
  calc _ ≤ 3 * tinyB C N * ∫ v in (a / 2)..π, |a0 v| := h1
    _ ≤ 3 * tinyB C N * (π * N / 2 * L1) := this
    _ = 3 * π / 2 * tinyB C N * N * L1 := by ring

/-- **Mass**: `2π - 2πτ ≤ ∫_{a/2}^π Φ(a - v) dv ≤ 2π`. -/
lemma mass_bounds {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 1 ≤ N) {a δ₀ : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π - δ₀) (hδ₀ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2) :
    2 * π - 2 * π * tinyB C N ≤ ∫ v in (a / 2)..π, LR.Phi N w (a - v) ∧
      ∫ v in (a / 2)..π, LR.Phi N w (a - v) ≤ 2 * π := by
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN
  have hπ := pi_gt_three
  have hτ := tinyB_nonneg hg hN
  have ha0 : 0 < a := lt_of_lt_of_le (by linarith) ha
  have hsub : ∫ v in (a / 2)..π, LR.Phi N w (a - v) = ∫ x in (a - π)..(a / 2), LR.Phi N w x := by
    rw [intervalIntegral.integral_comp_sub_left (fun x => LR.Phi N w x) a]
    congr 1; ring
  have hII := intervalIntegrable_phi N w
  have htot : (∫ x in (-π)..(a - π), LR.Phi N w x) + (∫ x in (a - π)..(a / 2), LR.Phi N w x) +
      (∫ x in (a / 2)..π, LR.Phi N w x) = 2 * π := by
    rw [intervalIntegral.integral_add_adjacent_intervals (hII _ _) (hII _ _),
      intervalIntegral.integral_add_adjacent_intervals (hII _ _) (hII _ _), integral_phi, hg.w0, mul_one]
  have hA0 : 0 ≤ ∫ x in (-π)..(a - π), LR.Phi N w x :=
    intervalIntegral.integral_nonneg (by linarith) fun x _ => hg.phi_nonneg x
  have hB0 : 0 ≤ ∫ x in (a / 2)..π, LR.Phi N w x :=
    intervalIntegral.integral_nonneg (by linarith) fun x _ => hg.phi_nonneg x
  have hA1 : ∫ x in (-π)..(a - π), LR.Phi N w x ≤ ∫ x in (-π)..(a - π), tinyB C N := by
    apply intervalIntegral.integral_mono_on (by linarith) (hII _ _) intervalIntegrable_const
    intro x hx
    exact hg.le_tiny_neg hN (by linarith [hx.2]) hx.1
  have hB1 : ∫ x in (a / 2)..π, LR.Phi N w x ≤ ∫ x in (a / 2)..π, tinyB C N := by
    apply intervalIntegral.integral_mono_on (by linarith) (hII _ _) intervalIntegrable_const
    intro x hx
    exact hg.le_tiny_pos hN (by linarith [hx.1]) (by linarith [hx.2])
  rw [intervalIntegral.integral_const, smul_eq_mul] at hA1 hB1
  rw [hsub]
  constructor
  · have : (a - π - -π) * tinyB C N + (π - a / 2) * tinyB C N ≤ 2 * π * tinyB C N := by
      nlinarith
    linarith
  · linarith

/-- `Φ(x)|x| ≤ C/(1 + N|x|)²` for `|x| ≤ π`. -/
lemma phi_mul_abs_le {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 1 ≤ N) {x : ℝ}
    (hx : |x| ≤ π) : LR.Phi N w x * |x| ≤ C / (1 + N * |x|) ^ 2 := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hp := hg.phi_nonneg x
  have hd := hg.decay x hx
  have hu : 1 ≤ 1 + N * |x| := by have := abs_nonneg x; nlinarith
  rw [le_div_iff₀ (by positivity)]
  have h1 : (N : ℝ) * |x| ≤ 1 + N * |x| := by linarith
  have h2 : (1 + N * |x|) ^ 3 ≤ (1 + N * |x|) ^ 30 := pow_le_pow_right₀ hu (by norm_num)
  have h3 : N * (LR.Phi N w x * |x| * (1 + N * |x|) ^ 2) ≤ LR.Phi N w x * (1 + N * |x|) ^ 3 := by
    have key : LR.Phi N w x * (N * |x|) * (1 + N * |x|) ^ 2 ≤
        LR.Phi N w x * (1 + N * |x|) * (1 + N * |x|) ^ 2 :=
      mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left h1 hp) (by positivity)
    calc N * (LR.Phi N w x * |x| * (1 + N * |x|) ^ 2) = LR.Phi N w x * (N * |x|) * (1 + N * |x|) ^ 2 := by
          ring
      _ ≤ LR.Phi N w x * (1 + N * |x|) * (1 + N * |x|) ^ 2 := key
      _ = LR.Phi N w x * (1 + N * |x|) ^ 3 := by ring
  have h4 : LR.Phi N w x * (1 + N * |x|) ^ 3 ≤ LR.Phi N w x * (1 + N * |x|) ^ 30 :=
    mul_le_mul_of_nonneg_left h2 hp
  have : (N : ℝ) * (LR.Phi N w x * |x| * (1 + N * |x|) ^ 2) ≤ N * C := by linarith
  exact le_of_mul_le_mul_left this hNpos

/-- `∫_c^d C/(1 + N(v - c))² dv ≤ C/N`. -/
lemma integral_inv_sq_le {N : ℕ} (hN : 1 ≤ N) {C c d : ℝ} (hC : 0 ≤ C) (hcd : c ≤ d) :
    ∫ v in c..d, C / (1 + N * (v - c)) ^ 2 ≤ C / N := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hpos : ∀ v ∈ Set.uIcc c d, 0 < 1 + N * (v - c) := by
    intro v hv
    rw [Set.uIcc_of_le hcd] at hv
    have : 0 ≤ v - c := by linarith [hv.1]
    positivity
  have hderiv : ∀ v ∈ Set.uIcc c d, HasDerivAt (fun v : ℝ => -(C / (N * (1 + N * (v - c)))))
      (C / (1 + N * (v - c)) ^ 2) v := by
    intro v hv
    have hu := hpos v hv
    have h1 : HasDerivAt (fun v : ℝ => N * (1 + N * (v - c))) (N * (N * 1)) v := by
      have := (((hasDerivAt_id v).sub_const c).const_mul (N : ℝ)).const_add 1 |>.const_mul (N : ℝ)
      simpa using this
    have h2 := ((h1.inv (by positivity)).const_mul C).neg
    convert h2 using 1
    · funext v; simp [div_eq_mul_inv]
    · field_simp
  have hint : IntervalIntegrable (fun v => C / (1 + N * (v - c)) ^ 2) volume c d := by
    apply ContinuousOn.intervalIntegrable
    apply ContinuousOn.div continuousOn_const (Continuous.continuousOn (by fun_prop))
    intro v hv
    have := hpos v hv
    positivity
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv hint]
  simp only [sub_self, mul_zero, add_zero, mul_one]
  have : 0 ≤ C / (N * (1 + N * (d - c))) := by
    have : 0 < 1 + N * (d - c) := hpos d (by rw [Set.uIcc_of_le hcd]; exact ⟨hcd, le_rfl⟩)
    positivity
  linarith

/-- **First moment**: `∫_{a/2}^π Φ(a - v)|a - v| dv ≤ 2C/N`. -/
lemma first_moment {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 1 ≤ N) {a : ℝ}
    (ha0 : 0 < a) (ha' : a ≤ π) :
    ∫ v in (a / 2)..π, LR.Phi N w (a - v) * |a - v| ≤ 2 * C / N := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hC := hg.C_nonneg hN
  have hcont : ∀ c d, IntervalIntegrable (fun v => LR.Phi N w (a - v) * |a - v|) volume c d :=
    fun c d => Continuous.intervalIntegrable (by unfold LR.Phi; fun_prop) c d
  rw [← intervalIntegral.integral_add_adjacent_intervals (hcont (a / 2) a) (hcont a π)]
  -- left piece, reflected: v ↦ a/2 + a - v
  have hL : ∫ v in (a / 2)..a, LR.Phi N w (a - v) * |a - v| ≤ C / N := by
    have hrefl : ∫ v in (a / 2)..a, LR.Phi N w (a - v) * |a - v| =
        ∫ u in (a / 2)..a, LR.Phi N w (u - a / 2) * |u - a / 2| := by
      have := intervalIntegral.integral_comp_sub_left
        (fun u => LR.Phi N w (u - a / 2) * |u - a / 2|) (a := a / 2) (b := a) (a / 2 + a)
      rw [show a / 2 + a - a = a / 2 by ring, show a / 2 + a - a / 2 = a by ring] at this
      rw [← this]
      refine intervalIntegral.integral_congr fun v _ => ?_
      try dsimp only
      rw [show a / 2 + a - v - a / 2 = a - v by ring]
    rw [hrefl]
    have hle : ∫ u in (a / 2)..a, LR.Phi N w (u - a / 2) * |u - a / 2| ≤
        ∫ u in (a / 2)..a, C / (1 + N * (u - a / 2)) ^ 2 := by
      apply intervalIntegral.integral_mono_on (by linarith)
        (Continuous.intervalIntegrable (by unfold LR.Phi; fun_prop) _ _)
      · apply ContinuousOn.intervalIntegrable
        apply ContinuousOn.div continuousOn_const (Continuous.continuousOn (by fun_prop))
        intro u hu
        rw [Set.uIcc_of_le (by linarith)] at hu
        have : 0 ≤ u - a / 2 := by linarith [hu.1]
        positivity
      intro u hu
      have hav : 0 ≤ u - a / 2 := by linarith [hu.1]
      have h := phi_mul_abs_le hg hN (x := u - a / 2) (by rw [abs_of_nonneg hav]; linarith [hu.2])
      rwa [abs_of_nonneg hav] at h ⊢
    exact hle.trans (integral_inv_sq_le hN hC (by linarith))
  -- right piece
  have hR : ∫ v in a..π, LR.Phi N w (a - v) * |a - v| ≤ C / N := by
    have hle : ∫ v in a..π, LR.Phi N w (a - v) * |a - v| ≤
        ∫ v in a..π, C / (1 + N * (v - a)) ^ 2 := by
      apply intervalIntegral.integral_mono_on ha' (hcont _ _)
      · apply ContinuousOn.intervalIntegrable
        apply ContinuousOn.div continuousOn_const (Continuous.continuousOn (by fun_prop))
        intro v hv
        rw [Set.uIcc_of_le ha'] at hv
        have : 0 ≤ v - a := by linarith [hv.1]
        positivity
      intro v hv
      have hav : 0 ≤ v - a := by linarith [hv.1]
      have h := phi_mul_abs_le hg hN (x := a - v)
        (by rw [abs_sub_comm, abs_of_nonneg hav]; linarith [hv.2])
      rwa [abs_sub_comm, abs_of_nonneg hav] at h ⊢
    exact hle.trans (integral_inv_sq_le hN hC ha')
  have : C / N + C / N = 2 * C / N := by ring
  linarith

/-- Lipschitz constant of `a₀` near `a`. -/
noncomputable def M1 (a δ₀ : ℝ) : ℝ :=
  (π / δ₀ + 1 / 2 + (log (2 * π / a) + log (2 * π / δ₀) + log (4 * π))) / (8 * π * (a / (2 * π)) ^ 3)

lemma a0_lip {a δ₀ v : ℝ} (ha0 : 0 < a) (hδ₀ : 0 < δ₀) (ha' : a ≤ π - δ₀)
    (hv : |v - a| ≤ min (a / 2) (δ₀ / 2)) : |a0 v - a0 a| ≤ M1 a δ₀ * |v - a| := by
  have hπ := pi_pos
  set h := min (a / 2) (δ₀ / 2) with hh
  have hh1 : h ≤ a / 2 := min_le_left _ _
  have hh2 : h ≤ δ₀ / 2 := min_le_right _ _
  have hh0 : 0 ≤ h := le_min (by linarith) (by linarith)
  set S := Icc (a - h) (a + h)
  have hmem : ∀ x ∈ S, a / 2 ≤ x ∧ x ≤ π - δ₀ / 2 := fun x hx => ⟨by linarith [hx.1], by linarith [hx.2]⟩
  have hf : ∀ x ∈ S, HasDerivWithinAt a0 (a0' x) S x := fun x hx =>
    (hasDerivAt_a0 (by linarith [(hmem x hx).1]) (by linarith [(hmem x hx).2])).hasDerivWithinAt
  have hbound : ∀ x ∈ S, ‖a0' x‖ ≤ M1 a δ₀ := by
    intro x hx
    obtain ⟨hx1, hx2⟩ := hmem x hx
    have hx0 : 0 < x := by linarith
    have hs : a / (2 * π) ≤ sin (x / 2) := by
      have := sin_half_ge hx0.le (by linarith)
      have : a / (2 * π) ≤ x / π := by rw [div_le_div_iff₀ (by positivity) hπ]; nlinarith
      linarith
    have hc : δ₀ / (2 * π) ≤ cos (x / 2) := by
      have := cos_half_ge hx0.le (by linarith)
      have : δ₀ / (2 * π) ≤ (π - x) / π := by rw [div_le_div_iff₀ (by positivity) hπ]; nlinarith
      linarith
    have h1 := abs_a0'_le (by positivity) (by positivity) hs hc
    have hell := abs_ell_le hx0 (by linarith)
    have l1 : log (π / x) ≤ log (2 * π / a) := log_le_log (by positivity) (by
      rw [div_le_div_iff₀ hx0 ha0]; nlinarith)
    have l2 : log (π / (π - x)) ≤ log (2 * π / δ₀) := log_le_log (div_pos hπ (by linarith)) (by
      rw [div_le_div_iff₀ (by linarith) hδ₀]; nlinarith)
    rw [Real.norm_eq_abs]
    refine h1.trans ?_
    unfold M1
    apply div_le_div_of_nonneg_right _ (by positivity)
    have : 1 / (2 * (δ₀ / (2 * π))) = π / δ₀ := by field_simp
    rw [this]
    linarith
  have := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le hf hbound (convex_Icc _ _)
    (⟨by linarith, by linarith⟩ : a ∈ S)
    (⟨by linarith [(abs_le.1 hv).1], by linarith [(abs_le.1 hv).2]⟩ : v ∈ S)
  simpa [Real.norm_eq_abs] using this

lemma M1_nonneg {a δ₀ : ℝ} (ha0 : 0 < a) (hδ₀ : 0 < δ₀) (ha' : a ≤ π - δ₀) : 0 ≤ M1 a δ₀ := by
  have hπ := pi_pos
  unfold M1
  apply div_nonneg _ (by positivity)
  have l1 : 0 ≤ log (2 * π / a) := log_nonneg (by rw [le_div_iff₀ ha0]; linarith)
  have l2 : 0 ≤ log (2 * π / δ₀) := log_nonneg (by rw [le_div_iff₀ hδ₀]; linarith)
  have l3 : 0 ≤ log (4 * π) := log_nonneg (by linarith [pi_gt_three])
  positivity

/-- Pointwise bound behind the main term. -/
lemma main_pointwise {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 1 ≤ N) {a δ₀ v : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π - δ₀) (hδ₀ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2)
    (hv : a / 2 ≤ v) (hv' : v ≤ π) :
    LR.Phi N w (a - v) * (a0 v - a0 a) ≤
      M1 a δ₀ * (LR.Phi N w (a - v) * |a - v|) + tinyB C N * (|a0 v| + |a0 a|) := by
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN
  have ha0 : 0 < a := lt_of_lt_of_le (by linarith) ha
  have hδ0 : 0 < δ₀ := by linarith
  have hp := hg.phi_nonneg (a - v)
  have hτ := tinyB_nonneg hg hN
  have hM := M1_nonneg ha0 hδ0 ha'
  have t1 : LR.Phi N w (a - v) * (a0 v - a0 a) ≤ LR.Phi N w (a - v) * |a0 v - a0 a| :=
    mul_le_mul_of_nonneg_left (le_abs_self _) hp
  have hr : 0 ≤ tinyB C N * (|a0 v| + |a0 a|) := by positivity
  have hl : 0 ≤ M1 a δ₀ * (LR.Phi N w (a - v) * |a - v|) := by positivity
  by_cases hc : |v - a| ≤ min (a / 2) (δ₀ / 2)
  · have hlip := a0_lip ha0 hδ0 ha' hc
    have : LR.Phi N w (a - v) * |a0 v - a0 a| ≤ M1 a δ₀ * (LR.Phi N w (a - v) * |a - v|) := by
      rw [abs_sub_comm a v]
      calc LR.Phi N w (a - v) * |a0 v - a0 a| ≤ LR.Phi N w (a - v) * (M1 a δ₀ * |v - a|) :=
            mul_le_mul_of_nonneg_left hlip hp
        _ = M1 a δ₀ * (LR.Phi N w (a - v) * |v - a|) := by ring
    linarith
  · push Not at hc
    have hfar : 1 / (2 * √(N : ℝ)) ≤ |a - v| := by
      rw [abs_sub_comm]
      rcases le_total (a / 2) (δ₀ / 2) with h | h
      · rw [min_eq_left h] at hc; linarith
      · rw [min_eq_right h] at hc; linarith
    have hle : |a - v| ≤ π := by rw [abs_le]; constructor <;> linarith
    have hΦ := hg.le_tiny hN hfar hle
    have : LR.Phi N w (a - v) * |a0 v - a0 a| ≤ tinyB C N * (|a0 v| + |a0 a|) := by
      calc LR.Phi N w (a - v) * |a0 v - a0 a| ≤ tinyB C N * |a0 v - a0 a| :=
            mul_le_mul_of_nonneg_right hΦ (abs_nonneg _)
        _ ≤ tinyB C N * (|a0 v| + |a0 a|) := mul_le_mul_of_nonneg_left (abs_sub _ _) hτ
    linarith

/-- **Main term**: `∫_{a/2}^π Φ(a-v) a₀(v) dv ≤ 2π a₀(a) + E₃`. -/
lemma main_bound {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 1 ≤ N) {a δ₀ : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π - δ₀) (hδ₀ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2) :
    ∫ v in (a / 2)..π, LR.Phi N w (a - v) * a0 v ≤
      2 * π * a0 a + 2 * π * tinyB C N * |a0 a| + M1 a δ₀ * (2 * C / N) +
        tinyB C N * (π * N / 2 * L1 + π * |a0 a|) := by
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN
  have hπ := pi_gt_three
  have ha0 : 0 < a := lt_of_lt_of_le (by linarith) ha
  have hδ0 : 0 < δ₀ := by linarith
  have hτ := tinyB_nonneg hg hN
  have hM := M1_nonneg ha0 hδ0 ha'
  have hA := intervalIntegrable_a0 (c := a / 2) (d := π) (by linarith) (by linarith) le_rfl
  have hcΦ : Continuous fun v => LR.Phi N w (a - v) := by unfold LR.Phi; fun_prop
  have hi1 : IntervalIntegrable (fun v => LR.Phi N w (a - v) * (a0 v - a0 a)) volume (a / 2) π :=
    (hA.sub intervalIntegrable_const).continuousOn_mul hcΦ.continuousOn
  have hi2 : IntervalIntegrable (fun v => LR.Phi N w (a - v)) volume (a / 2) π :=
    hcΦ.intervalIntegrable _ _
  have hsplit : ∫ v in (a / 2)..π, LR.Phi N w (a - v) * a0 v =
      (∫ v in (a / 2)..π, LR.Phi N w (a - v) * (a0 v - a0 a)) +
        a0 a * ∫ v in (a / 2)..π, LR.Phi N w (a - v) := by
    rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_add hi1 (hi2.const_mul _)]
    refine intervalIntegral.integral_congr fun v _ => ?_
    ring
  rw [hsplit]
  -- first part
  have hR1 : IntervalIntegrable (fun v => M1 a δ₀ * (LR.Phi N w (a - v) * |a - v|)) volume (a / 2) π :=
    (Continuous.intervalIntegrable (by unfold LR.Phi; fun_prop) _ _)
  have hR2 : IntervalIntegrable (fun v => tinyB C N * (|a0 v| + |a0 a|)) volume (a / 2) π :=
    (hA.abs.add intervalIntegrable_const).const_mul _
  have hP : ∫ v in (a / 2)..π, LR.Phi N w (a - v) * (a0 v - a0 a) ≤
      ∫ v in (a / 2)..π, (M1 a δ₀ * (LR.Phi N w (a - v) * |a - v|) + tinyB C N * (|a0 v| + |a0 a|)) := by
    apply intervalIntegral.integral_mono_on_of_le_Ioo (by linarith) hi1 (hR1.add hR2)
    intro v hv
    exact main_pointwise hg hN ha ha' hδ₀ hv.1.le hv.2.le
  rw [intervalIntegral.integral_add hR1 hR2, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_add hA.abs intervalIntegrable_const,
    intervalIntegral.integral_const, smul_eq_mul] at hP
  have hfm := first_moment hg hN ha0 (by linarith)
  have hab := integral_abs_a0_le hN (c := a / 2) (by linarith) (by linarith)
  have hP' : ∫ v in (a / 2)..π, LR.Phi N w (a - v) * (a0 v - a0 a) ≤
      M1 a δ₀ * (2 * C / N) + tinyB C N * (π * N / 2 * L1 + π * |a0 a|) := by
    have e1 : M1 a δ₀ * ∫ v in (a / 2)..π, LR.Phi N w (a - v) * |a - v| ≤ M1 a δ₀ * (2 * C / N) :=
      mul_le_mul_of_nonneg_left hfm hM
    have e2 : tinyB C N * ((∫ v in (a / 2)..π, |a0 v|) + (π - a / 2) * |a0 a|) ≤
        tinyB C N * (π * N / 2 * L1 + π * |a0 a|) := by
      apply mul_le_mul_of_nonneg_left _ hτ
      have : (π - a / 2) * |a0 a| ≤ π * |a0 a| := by
        apply mul_le_mul_of_nonneg_right (by linarith) (abs_nonneg _)
      linarith
    linarith
  -- second part
  obtain ⟨hJ1, hJ2⟩ := mass_bounds hg hN ha ha' hδ₀
  set J := ∫ v in (a / 2)..π, LR.Phi N w (a - v)
  have hJ : a0 a * J ≤ 2 * π * a0 a + 2 * π * tinyB C N * |a0 a| := by
    have h1 : |J - 2 * π| ≤ 2 * π * tinyB C N := by rw [abs_le]; constructor <;> linarith
    have h2 : a0 a * (J - 2 * π) ≤ |a0 a| * (2 * π * tinyB C N) := by
      calc a0 a * (J - 2 * π) ≤ |a0 a * (J - 2 * π)| := le_abs_self _
        _ = |a0 a| * |J - 2 * π| := abs_mul _ _
        _ ≤ |a0 a| * (2 * π * tinyB C N) := mul_le_mul_of_nonneg_left h1 (abs_nonneg _)
    nlinarith
  linarith

/-- **Trend lower bound with explicit error**. -/
theorem trendI_lower {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 2 ≤ N) {a δ₀ : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π - δ₀) (hδ₀ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2) :
    log (2 * π * a) / a ^ 2 - (3 + |log (2 * π * a)| +
      ((π / 4) * (2 + log π + log (4 * π) + 3 * log N) + π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1) +
      3 * π / 2 * tinyB C N * N * L1 +
      (2 * π * tinyB C N * |a0 a| + M1 a δ₀ * (2 * C / N) +
        tinyB C N * (π * N / 2 * L1 + π * |a0 a|))) ≤ LR.trend N w Iclosed a := by
  have hN1 : 1 ≤ N := by omega
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN1
  have hπ := pi_gt_three
  have ha0 : 0 < a := lt_of_lt_of_le (by linarith) ha
  have hδ0 : 0 < δ₀ := by linarith
  rw [trend_Iclosed_eq]
  have hI := intervalIntegrable_B_a0 N w a
  have hs1 : Set.uIcc 0 (a / 2) ⊆ Set.uIcc 0 π := by
    rw [Set.uIcc_of_le (by linarith), Set.uIcc_of_le pi_pos.le]
    exact Icc_subset_Icc le_rfl (by linarith)
  have hs2 : Set.uIcc (a / 2) π ⊆ Set.uIcc 0 π := by
    rw [Set.uIcc_of_le (by linarith), Set.uIcc_of_le pi_pos.le]
    exact Icc_subset_Icc (by linarith) le_rfl
  have hI1 : IntervalIntegrable (fun v => LR.B N w a v * a0 v) volume 0 (a / 2) := hI.mono_set hs1
  have hI2 : IntervalIntegrable (fun v => LR.B N w a v * a0 v) volume (a / 2) π := hI.mono_set hs2
  rw [← intervalIntegral.integral_add_adjacent_intervals hI1 hI2]
  have h1 := near_zero_bound hg hN ha (by linarith)
  have hA := intervalIntegrable_a0 (c := a / 2) (d := π) (by linarith) (by linarith) le_rfl
  have hc1 : IntervalIntegrable (fun v => (2 * LR.Phi N w a - LR.Phi N w (a + v)) * a0 v) volume (a / 2) π :=
    hA.continuousOn_mul (Continuous.continuousOn (by unfold LR.Phi; fun_prop))
  have hc2 : IntervalIntegrable (fun v => LR.Phi N w (a - v) * a0 v) volume (a / 2) π :=
    hA.continuousOn_mul (Continuous.continuousOn (by unfold LR.Phi; fun_prop))
  have hsplit : ∫ v in (a / 2)..π, LR.B N w a v * a0 v =
      (∫ v in (a / 2)..π, (2 * LR.Phi N w a - LR.Phi N w (a + v)) * a0 v) -
        ∫ v in (a / 2)..π, LR.Phi N w (a - v) * a0 v := by
    rw [← intervalIntegral.integral_sub hc1 hc2]
    refine intervalIntegral.integral_congr fun v _ => ?_
    simp only [LR.B]; ring
  rw [hsplit]
  have h2 := mid_bound hg hN1 ha ha' hδ₀
  have h3 := main_bound hg hN1 ha ha' hδ₀
  have h4 := ell_div_le ha0 (by linarith)
  have e : 2 * π * a0 a = ell a / (4 * sin (a / 2) ^ 2) := by
    unfold a0; field_simp; ring
  rw [← e, neg_div] at h4
  have h2' := (abs_le.1 h2).1
  linarith

end LeanLi

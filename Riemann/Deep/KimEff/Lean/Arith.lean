import Riemann.Deep.KimEff.Lean.Core
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# Explicit Kim–Lee wedge, part 3: the arithmetic of DERIVATION §10 ((ii) ⇒ (iii))

Elementary inequalities, with exp/log/π bounds from Mathlib (exp 1 to 9 digits, π > 3.14,
Taylor lower bounds of exp) and exact rational arithmetic (norm_num): H2, H3 and the
power form with c = 1.
-/

namespace KimEff

noncomputable section

open Real

/-- Z(u) = u (π e^{2u} - 9/4) (DERIVATION §0). -/
def Zfun (u : ℝ) : ℝ := u * (π * exp (2 * u) - 9 / 4)

/-- R(n, u) = (9/10)(n - 1/2)^2/(u + 1)^2; R_n = R(n, u(n)) (DERIVATION §0, §7). -/
def Rfun (n u : ℝ) : ℝ := 9 / 10 * (n - 1 / 2) ^ 2 / (u + 1) ^ 2

/-- n_* = ceil Z(12) = 998,616,377,771 (certify.py). -/
def nStar : ℕ := 998616377771

/-- a = 0.37338196661320752: Bebel's Lean-checked half-plane Re ρ ≤ 1/2 + a. -/
def a0 : ℝ := 0.37338196661320752

lemma a0_pos : 0 < a0 := by norm_num [a0]

/-- Z is strictly increasing on [0, ∞), so u(n) is unique. -/
lemma Zfun_strictMonoOn : StrictMonoOn Zfun (Set.Ici 0) := by
  intro u hu v hv huv
  simp only [Set.mem_Ici] at hu hv
  unfold Zfun
  have h1 : exp (2 * u) ≤ exp (2 * v) := exp_le_exp.mpr (by linarith)
  have h2 : 1 ≤ exp (2 * v) := one_le_exp (by linarith)
  have hpi : 3 < π := pi_gt_three
  have hA : 0 ≤ π * (u * (exp (2 * v) - exp (2 * u))) :=
    mul_nonneg pi_pos.le (mul_nonneg hu (sub_nonneg.mpr h1))
  have hB : 0 ≤ π * ((v - u) * (exp (2 * v) - 1)) :=
    mul_nonneg pi_pos.le (mul_nonneg (sub_nonneg.mpr huv.le) (sub_nonneg.mpr h2))
  have hC : 0 < (π - 9 / 4) * (v - u) := mul_pos (by linarith) (sub_pos.mpr huv)
  nlinarith

lemma nStar_le {n : ℝ} (hn : (nStar : ℝ) ≤ n) : (998616377771 : ℝ) ≤ n := by
  simpa [nStar] using hn

/-- log n ≥ 27 for n ≥ n_*. -/
lemma log_ge_27 {n : ℝ} (hn : (nStar : ℝ) ≤ n) : 27 ≤ log n := by
  have hn' := nStar_le hn
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) hn'
  rw [le_log_iff_exp_le hn0]
  have h1 : exp 27 = exp 1 ^ 27 := by rw [← exp_nat_mul]; norm_num
  have h2 : exp 1 ^ 27 < 2.7182818286 ^ 27 :=
    pow_lt_pow_left₀ exp_one_lt_d9 (exp_pos 1).le (by norm_num)
  have h3 : (2.7182818286 : ℝ) ^ 27 ≤ 998616377771 := by norm_num
  linarith

/-- H2 (DERIVATION §10): for n ≥ n_* there is u > 0 with Z(u) = n and u + 1 ≤ (log n)/2. -/
theorem exists_u {n : ℝ} (hn : (nStar : ℝ) ≤ n) :
    ∃ u : ℝ, 0 < u ∧ Zfun u = n ∧ u + 1 ≤ log n / 2 := by
  have hn' := nStar_le hn
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) hn'
  have hlog := log_ge_27 hn
  have hL0 : 0 ≤ log n / 2 - 1 := by linarith
  have hZL : n ≤ Zfun (log n / 2 - 1) := by
    unfold Zfun
    have he : exp (2 * (log n / 2 - 1)) = n / exp 2 := by
      rw [show 2 * (log n / 2 - 1) = log n - 2 by ring, exp_sub, exp_log hn0]
    rw [he]
    have he2 : exp 2 < 7.39 := by
      have h : exp 2 = exp 1 ^ 2 := by rw [← exp_nat_mul]; norm_num
      rw [h]
      nlinarith [exp_one_lt_d9, exp_pos 1]
    have hpi : 3.14 < π := pi_gt_d2
    have hq : 0.42 * n ≤ π * (n / exp 2) := by
      rw [mul_div_assoc', le_div_iff₀ (exp_pos 2)]
      nlinarith [mul_pos (sub_pos.mpr hpi) hn0, mul_pos (sub_pos.mpr he2) hn0]
    nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ log n / 2 - 1 - 12.5)
      (by linarith : (0:ℝ) ≤ π * (n / exp 2) - 9 / 4)]
  have hcont : ContinuousOn Zfun (Set.Icc 0 (log n / 2 - 1)) := by
    unfold Zfun
    fun_prop
  have hZ0 : Zfun 0 = 0 := by simp [Zfun]
  obtain ⟨u, ⟨hu0, huL⟩, hu⟩ :=
    intermediate_value_Icc hL0 hcont ⟨by rw [hZ0]; exact hn0.le, hZL⟩
  refine ⟨u, ?_, hu, by linarith⟩
  rcases hu0.lt_or_eq with h | h
  · exact h
  · rw [← h, hZ0] at hu
    linarith

/-- D(R) ≥ (R - a^2)/(4a^2) (DERIVATION §10). -/
lemma sub_div_le_Dfun {a R : ℝ} (ha : 0 < a) (hR : 0 < R) :
    (R - a ^ 2) / (4 * a ^ 2) ≤ Dfun a R := by
  unfold Dfun
  have h4 : 0 < 4 * a ^ 2 := by have := pow_pos ha 2; linarith
  have hX : 0 < R + a ^ 2 := add_pos hR (pow_pos ha 2)
  rw [div_le_div_iff₀ h4 (mul_pos h4 hX)]
  nlinarith [pow_pos ha 6]

/-- log^2 n ≤ 2n. -/
lemma log_sq_le {n : ℝ} (hn0 : 1 ≤ n) : log n ^ 2 ≤ 2 * n := by
  have hl : 0 ≤ log n := log_nonneg hn0
  have h := quadratic_le_exp_of_nonneg hl
  rw [exp_log (by linarith)] at h
  nlinarith

/-- H3 (DERIVATION §10): 6.45 n^2/log^2 n ≤ D(R(n, u)) when u > 0, u + 1 ≤ (log n)/2, n ≥ n_*. -/
theorem le_Dfun_of_le_log {n u : ℝ} (hn : (nStar : ℝ) ≤ n) (hu : 0 < u)
    (hub : u + 1 ≤ log n / 2) : 6.45 * n ^ 2 / log n ^ 2 ≤ Dfun a0 (Rfun n u) := by
  have hn' := nStar_le hn
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) hn'
  have hlog := log_ge_27 hn
  have hl : 0 < log n := by linarith
  have hl2 : 0 < log n ^ 2 := pow_pos hl 2
  have hnh : 0 < n - 1 / 2 := by linarith
  have hR : 0 < Rfun n u := by
    unfold Rfun
    exact div_pos (mul_pos (by norm_num) (pow_pos hnh 2)) (pow_pos (by linarith) 2)
  have hRlow : 3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 ≤ Rfun n u := by
    unfold Rfun
    have h1 : (u + 1) ^ 2 ≤ (log n / 2) ^ 2 := pow_le_pow_left₀ (by linarith) hub 2
    calc 3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 = 9 / 10 * (n - 1 / 2) ^ 2 / (log n / 2) ^ 2 := by
          field_simp
          ring
      _ ≤ 9 / 10 * (n - 1 / 2) ^ 2 / (u + 1) ^ 2 :=
          div_le_div_of_nonneg_left (by positivity) (pow_pos (by linarith) 2) h1
  have h4 : 0 < 4 * a0 ^ 2 := by have := pow_pos a0_pos 2; linarith
  have hmain : 6.45 * n ^ 2 / log n ^ 2
      ≤ (3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 - a0 ^ 2) / (4 * a0 ^ 2) := by
    rw [div_le_div_iff₀ hl2 h4, sub_mul, div_mul_cancel₀ _ hl2.ne']
    have hls := log_sq_le (by linarith : (1:ℝ) ≤ n)
    have hnn : 998616377771 * n ≤ n * n := mul_le_mul_of_nonneg_right hn' hn0.le
    norm_num [a0]
    nlinarith
  calc 6.45 * n ^ 2 / log n ^ 2
      ≤ (3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 - a0 ^ 2) / (4 * a0 ^ 2) := hmain
    _ ≤ (Rfun n u - a0 ^ 2) / (4 * a0 ^ 2) :=
        div_le_div_of_nonneg_right (by linarith) h4.le
    _ ≤ Dfun a0 (Rfun n u) := sub_div_le_Dfun a0_pos hR

/-- log n_* ≤ 27.6305 (true value 27.62963...). -/
lemma log_nStar_le : log (998616377771 : ℝ) ≤ 27.6305 := by
  rw [log_le_iff_le_exp (by norm_num)]
  have h1 : exp 27.6305 = exp 1 ^ 27 * exp 0.6305 := by
    rw [← exp_nat_mul, ← exp_add]; norm_num
  have h2 : (2.7182818283 : ℝ) ^ 27 ≤ exp 1 ^ 27 :=
    pow_le_pow_left₀ (by norm_num) exp_one_gt_d9.le 27
  have h3 : ∑ i ∈ Finset.range 10, (0.6305 : ℝ) ^ i / (i.factorial : ℝ) ≤ exp 0.6305 :=
    sum_le_exp_of_nonneg (by norm_num) 10
  have h4 : (998616377771 : ℝ) ≤
      2.7182818283 ^ 27 * ∑ i ∈ Finset.range 10, (0.6305 : ℝ) ^ i / (i.factorial : ℝ) := by
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    norm_num
  rw [h1]
  calc (998616377771 : ℝ)
      ≤ 2.7182818283 ^ 27 * ∑ i ∈ Finset.range 10, (0.6305 : ℝ) ^ i / (i.factorial : ℝ) := h4
    _ ≤ exp 1 ^ 27 * exp 0.6305 :=
        mul_le_mul h2 h3 (Finset.sum_nonneg fun i _ => by positivity) (by positivity)

/-- log^2 n ≤ 6.45 · 1.1853e-10 · n for n ≥ n_* (power form c = 1, DERIVATION §10). -/
lemma log_sq_le_power_one {n : ℝ} (hn : (nStar : ℝ) ≤ n) :
    log n ^ 2 ≤ 6.45 * 1.1853e-10 * n := by
  have hn' := nStar_le hn
  set N : ℝ := 998616377771 with hN
  have hN0 : 0 < N := by norm_num [hN]
  set t := n / N with ht
  have ht1 : 1 ≤ t := by rw [ht, le_div_iff₀ hN0]; linarith
  have hnt : n = N * t := by rw [ht]; field_simp
  have hst : 1 ≤ Real.sqrt t := Real.one_le_sqrt.mpr ht1
  have hlogt : log t ≤ 2 * (Real.sqrt t - 1) := by
    have h := log_le_sub_one_of_pos (Real.sqrt_pos.mpr (by linarith : (0:ℝ) < t))
    rw [log_sqrt (by linarith)] at h
    linarith
  have hlogN := log_nStar_le
  have hlogn : log n = log N + log t := by
    rw [hnt, log_mul hN0.ne' (by linarith)]
  have hl0 : 0 ≤ log n := log_nonneg (by linarith)
  have hle : log n ≤ 27.6305 * Real.sqrt t := by
    rw [hlogn]
    nlinarith
  have hsq : Real.sqrt t ^ 2 = t := Real.sq_sqrt (by linarith)
  have h2 : log n ^ 2 ≤ 27.6305 ^ 2 * t := by
    calc log n ^ 2 ≤ (27.6305 * Real.sqrt t) ^ 2 := pow_le_pow_left₀ hl0 hle 2
      _ = 27.6305 ^ 2 * t := by rw [mul_pow, hsq]
  have h3 : (27.6305 : ℝ) ^ 2 * t ≤ 6.45 * 1.1853e-10 * n := by
    rw [hnt]
    have : (27.6305 : ℝ) ^ 2 ≤ 6.45 * 1.1853e-10 * N := by norm_num [hN]
    nlinarith
  linarith

/-- Power form c = 1: n ≥ n_* and n ≥ 1.1853e-10 d imply d ≤ 6.45 n^2/log^2 n. -/
theorem le_of_power_one {n d : ℝ} (hn : (nStar : ℝ) ≤ n) (hd : 1.1853e-10 * d ≤ n) :
    d ≤ 6.45 * n ^ 2 / log n ^ 2 := by
  have hn' := nStar_le hn
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) hn'
  have hlog := log_ge_27 hn
  have hl2 : 0 < log n ^ 2 := pow_pos (by linarith) 2
  have h := log_sq_le_power_one hn
  rw [le_div_iff₀ hl2]
  nlinarith [mul_le_mul_of_nonneg_left h (by linarith : (0:ℝ) ≤ 6.45 * n), mul_le_mul_of_nonneg_left hd (by linarith : (0:ℝ) ≤ log n ^ 2)]

end

end KimEff


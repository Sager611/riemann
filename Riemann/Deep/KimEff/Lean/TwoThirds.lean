import Riemann.Deep.KimEff.Lean.Power

/-!
# Explicit Kim–Lee wedge, part 10: the two-thirds wedge is contained (DERIVATION §10)

Self-contained comparison: for d > d0 = 1.04710e8, every (d, n) of the previous record
n^2 log(n+2) ≥ 24 d^3 has n > n_* and d ≤ 6.45 n^2/log^2 n, so it lies in the new region.
-/

namespace KimEff

noncomputable section

open Real

/-- log(n_* + 2) ≤ 27.62964. -/
lemma log_nStar_add_two_le : log (998616377773 : ℝ) ≤ 27.62964 := by
  rw [log_le_iff_le_exp (by norm_num)]
  have h1 : exp 27.62964 = exp 1 ^ 27 * exp 0.62964 := by
    rw [← exp_nat_mul, ← exp_add]; norm_num
  have h2 : (2.7182818283 : ℝ) ^ 27 ≤ exp 1 ^ 27 :=
    pow_le_pow_left₀ (by norm_num) exp_one_gt_d9.le 27
  have h3 : ∑ i ∈ Finset.range 10, (0.62964 : ℝ) ^ i / (i.factorial : ℝ) ≤ exp 0.62964 :=
    sum_le_exp_of_nonneg (by norm_num) 10
  have h4 : (998616377773 : ℝ) ≤
      2.7182818283 ^ 27 * ∑ i ∈ Finset.range 10, (0.62964 : ℝ) ^ i / (i.factorial : ℝ) := by
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    norm_num
  rw [h1]
  calc (998616377773 : ℝ)
      ≤ 2.7182818283 ^ 27 * ∑ i ∈ Finset.range 10, (0.62964 : ℝ) ^ i / (i.factorial : ℝ) := h4
    _ ≤ exp 1 ^ 27 * exp 0.62964 :=
        mul_le_mul h2 h3 (Finset.sum_nonneg fun i _ => by positivity) (by positivity)

/-- For d > 1.04710e8: n^2 log(n+2) ≥ 24 d^3 implies n > n_* and d ≤ 6.45 n^2/log^2 n. -/
theorem twoThirds_sub {n d : ℝ} (hn0 : 0 ≤ n) (hd : 1.04710e8 < d)
    (hw : 24 * d ^ 3 ≤ n ^ 2 * log (n + 2)) :
    (998616377771 : ℝ) < n ∧ d ≤ 6.45 * n ^ 2 / log n ^ 2 := by
  have hd0 : 0 < d := by linarith
  have hN : (998616377771 : ℝ) < n := by
    by_contra h
    have h' := not_lt.mp h
    have h1 : log (n + 2) ≤ log 998616377773 := log_le_log (by linarith) (by linarith)
    have h2 : 0 ≤ log (n + 2) := log_nonneg (by linarith)
    have h3 : n ^ 2 ≤ 998616377771 ^ 2 := pow_le_pow_left₀ hn0 h' 2
    have h4 : n ^ 2 * log (n + 2) ≤ 998616377771 ^ 2 * 27.62964 :=
      mul_le_mul h3 (h1.trans log_nStar_add_two_le) h2 (by norm_num)
    have h5 : (998616377771 : ℝ) ^ 2 * 27.62964 ≤ 24 * 1.04710e8 ^ 3 := by norm_num
    have h6 : (1.04710e8 : ℝ) ^ 3 < d ^ 3 := pow_lt_pow_left₀ hd (by norm_num) (by norm_num)
    linarith
  refine ⟨hN, ?_⟩
  have hnS : (nStar : ℝ) ≤ n := by
    have : (nStar : ℝ) = 998616377771 := by norm_num [nStar]
    rw [this]
    exact hN.le
  have hl := log_ge_27 hnS
  have hn1 : 1 ≤ n := by linarith
  have hl0 : 0 < log n := by linarith
  have hls : log n ^ 2 ≤ 2 * n := log_sq_le hn1
  have hl6 : log n ^ 6 ≤ 8 * n ^ 3 := by
    have e : log n ^ 6 = (log n ^ 2) ^ 3 := by ring
    rw [e]
    calc (log n ^ 2) ^ 3 ≤ (2 * n) ^ 3 := pow_le_pow_left₀ (sq_nonneg _) hls 3
      _ = 8 * n ^ 3 := by ring
  have hlog2 : log (n + 2) ≤ 2 * n := by
    have := log_le_sub_one_of_pos (by linarith : (0:ℝ) < n + 2)
    linarith
  have hM0 : 0 ≤ 6.45 * n ^ 2 / log n ^ 2 := by positivity
  have hcube : d ^ 3 ≤ (6.45 * n ^ 2 / log n ^ 2) ^ 3 := by
    rw [div_pow, mul_pow, ← pow_mul, ← pow_mul, le_div_iff₀ (pow_pos hl0 _)]
    show d ^ 3 * log n ^ 6 ≤ 6.45 ^ 3 * n ^ 6
    have hA : 24 * d ^ 3 * log n ^ 6 ≤ n ^ 2 * log (n + 2) * log n ^ 6 :=
      mul_le_mul_of_nonneg_right hw (by positivity)
    have hB : n ^ 2 * log (n + 2) * log n ^ 6 ≤ n ^ 2 * (2 * n) * (8 * n ^ 3) :=
      mul_le_mul (mul_le_mul_of_nonneg_left hlog2 (sq_nonneg n)) hl6 (by positivity)
        (mul_nonneg (sq_nonneg n) (by linarith))
    have hn6 : 0 ≤ n ^ 6 := by positivity
    nlinarith [hA, hB, hn6]
  exact (pow_le_pow_iff_left₀ hd0.le hM0 (by norm_num)).mp hcube

/-- Given the inputs, the two-thirds wedge with d > 1.04710e8 lies in the hyperbolic region. -/
theorem jensen_hyperbolic_twoThirds {γ : ℕ → ℝ} {F : ℕ → ℂ → ℂ} (H : Inputs γ F a0) {n d : ℕ}
    (hd : 1.04710e8 < (d : ℝ)) (hw : 24 * (d : ℝ) ^ 3 ≤ (n : ℝ) ^ 2 * log ((n : ℝ) + 2)) :
    Hyperbolic (jensen γ d n) := by
  obtain ⟨hN, hle⟩ := twoThirds_sub (Nat.cast_nonneg n) hd hw
  have hn : nStar ≤ n := by
    have h1 : ((998616377771 : ℕ) : ℝ) < n := by exact_mod_cast hN
    have h2 : 998616377771 < n := by exact_mod_cast h1
    show 998616377771 ≤ n
    omega
  exact jensen_hyperbolic_iii H hn hle

end

end KimEff


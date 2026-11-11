import Riemann.Deep.KimEff.Lean.Theorem

/-!
# Explicit Kim–Lee wedge, part 6: the power forms (c, C) of DERIVATION §10 (iii)

n ≥ max(n_*, C d^c) implies d ≤ 6.45 n^2/log^2 n, hence (by jensen_hyperbolic_iii) J^{d,n} is
hyperbolic, for (c, C) = (0.51, 7.6890), (0.55, 0.87160), (0.6, 0.069829) (c = 1 is in Arith).
Writing c = p/q, the claim reduces to log^{2p} n ≤ 6.45^p C^q n^{2p-q}, proved with
x ≤ e^{x-1}, exp 1 > 2.7182818283, log n_* ≤ 27.62964 and exact rational arithmetic.
-/

namespace KimEff

noncomputable section

open Real

/-- log n_* ≤ 27.62964 (true value 27.6296365...). -/
lemma log_nStar_le_fine : log (998616377771 : ℝ) ≤ 27.62964 := by
  rw [log_le_iff_le_exp (by norm_num)]
  have h1 : exp 27.62964 = exp 1 ^ 27 * exp 0.62964 := by
    rw [← exp_nat_mul, ← exp_add]; norm_num
  have h2 : (2.7182818283 : ℝ) ^ 27 ≤ exp 1 ^ 27 :=
    pow_le_pow_left₀ (by norm_num) exp_one_gt_d9.le 27
  have h3 : ∑ i ∈ Finset.range 10, (0.62964 : ℝ) ^ i / (i.factorial : ℝ) ≤ exp 0.62964 :=
    sum_le_exp_of_nonneg (by norm_num) 10
  have h4 : (998616377771 : ℝ) ≤
      2.7182818283 ^ 27 * ∑ i ∈ Finset.range 10, (0.62964 : ℝ) ^ i / (i.factorial : ℝ) := by
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    norm_num
  rw [h1]
  calc (998616377771 : ℝ)
      ≤ 2.7182818283 ^ 27 * ∑ i ∈ Finset.range 10, (0.62964 : ℝ) ^ i / (i.factorial : ℝ) := h4
    _ ≤ exp 1 ^ 27 * exp 0.62964 :=
        mul_le_mul h2 h3 (Finset.sum_nonneg fun i _ => by positivity) (by positivity)

/-- If M < d were true, d^p > M^p ≥ (n/C)^q ≥ d^p. -/
lemma le_of_pow_le {n d C ℓ : ℝ} {p q : ℕ} (hp : p ≠ 0)
    (hM : 0 ≤ 6.45 * n ^ 2 / ℓ ^ 2) (hdq : d ^ p ≤ (n / C) ^ q)
    (hkey : (n / C) ^ q ≤ (6.45 * n ^ 2 / ℓ ^ 2) ^ p) : d ≤ 6.45 * n ^ 2 / ℓ ^ 2 := by
  by_contra h
  have h' := not_le.mp h
  have := pow_lt_pow_left₀ h' hM hp
  linarith

lemma key_of_bound {n C ℓ : ℝ} {p q : ℕ} (hn : 0 < n) (hC : 0 < C) (hℓ : 0 < ℓ)
    (hqp : q ≤ 2 * p) (hb : ℓ ^ (2 * p) ≤ 6.45 ^ p * C ^ q * n ^ (2 * p - q)) :
    (n / C) ^ q ≤ (6.45 * n ^ 2 / ℓ ^ 2) ^ p := by
  rw [div_pow, div_pow, mul_pow, ← pow_mul, ← pow_mul,
    div_le_div_iff₀ (pow_pos hC q) (pow_pos hℓ (2 * p))]
  have hsplit : n ^ (2 * p) = n ^ q * n ^ (2 * p - q) := by
    rw [← pow_add, Nat.add_sub_cancel' hqp]
  rw [hsplit]
  have h := mul_le_mul_of_nonneg_left hb (pow_pos hn q).le
  nlinarith [h]

/-- From C d^{p/q} ≤ n: d^p ≤ (n/C)^q. -/
lemma pow_le_of_rpow {n d C : ℝ} {p q : ℕ} (hq : q ≠ 0) (hd0 : 0 ≤ d) (hC : 0 < C)
    (hd : C * d ^ ((p : ℝ) / q) ≤ n) : d ^ p ≤ (n / C) ^ q := by
  have h1 : d ^ ((p : ℝ) / q) ≤ n / C := by
    rw [le_div_iff₀ hC]
    linarith
  have h2 := pow_le_pow_left₀ (Real.rpow_nonneg hd0 _) h1 q
  have h3 : (d ^ ((p : ℝ) / q)) ^ q = d ^ p := by
    rw [← Real.rpow_natCast (d ^ ((p : ℝ) / q)) q, ← Real.rpow_mul hd0, ← Real.rpow_natCast d p]
    congr 1
    field_simp
  rw [h3] at h2
  exact h2

/-- For n ≥ n_* and k ≤ 27 m: log^k n ≤ 27.62964^k (n/n_*)^m (log^k n / n^m decreasing). -/
lemma log_pow_le {n : ℝ} (hn : (nStar : ℝ) ≤ n) (k m : ℕ) (hkm : (k : ℝ) ≤ 27 * m) :
    log n ^ k ≤ 27.62964 ^ k * (n / 998616377771) ^ m := by
  have hn' := nStar_le hn
  have hN0 : (0 : ℝ) < 998616377771 := by norm_num
  have hn0 : 0 < n := lt_of_lt_of_le hN0 hn'
  set L0 := log (998616377771 : ℝ) with hL0def
  have hL0 : 27 ≤ L0 := log_ge_27 (by norm_num [nStar])
  have hL0U : L0 ≤ 27.62964 := log_nStar_le_fine
  have hl : L0 ≤ log n := log_le_log hN0 hn'
  have hL0pos : 0 < L0 := by linarith
  have h1 : log n / L0 ≤ exp (log n / L0 - 1) := by
    have := add_one_le_exp (log n / L0 - 1)
    linarith
  have h2 : (log n / L0) ^ k ≤ exp (log n / L0 - 1) ^ k :=
    pow_le_pow_left₀ (div_nonneg (by linarith) hL0pos.le) h1 k
  have h3 : exp (log n / L0 - 1) ^ k ≤ (n / 998616377771) ^ m := by
    rw [← exp_nat_mul]
    have hq : n / 998616377771 = exp (log n - L0) := by
      rw [exp_sub, exp_log hn0, hL0def, exp_log hN0]
    rw [hq, ← exp_nat_mul]
    apply exp_le_exp.mpr
    have e1 : log n / L0 - 1 = (log n - L0) / L0 := by field_simp
    rw [e1]
    have hd : 0 ≤ log n - L0 := by linarith
    have e2 : (k : ℝ) * ((log n - L0) / L0) ≤ (m : ℝ) * (log n - L0) := by
      rw [mul_div_assoc', div_le_iff₀ hL0pos]
      have hm0 : (0 : ℝ) ≤ m := Nat.cast_nonneg m
      nlinarith [mul_le_mul_of_nonneg_right hkm hd, mul_le_mul_of_nonneg_left hL0 (mul_nonneg hm0 hd)]
    exact e2
  have h4 : (log n / L0) ^ k ≤ (n / 998616377771) ^ m := h2.trans h3
  rw [div_pow, div_le_iff₀ (pow_pos hL0pos k)] at h4
  calc log n ^ k ≤ (n / 998616377771) ^ m * L0 ^ k := h4
    _ ≤ (n / 998616377771) ^ m * 27.62964 ^ k :=
        mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hL0pos.le hL0U k)
          (pow_nonneg (div_nonneg hn0.le hN0.le) m)
    _ = 27.62964 ^ k * (n / 998616377771) ^ m := mul_comm _ _

/-- c = 0.51: log^{102} n ≤ 6.45^51 · 7.6890^100 · n^2 for every n > 0
(the maximum of log^{102} n / n^2 is at log n = 51). -/
lemma hb_051 {n : ℝ} (hn0 : 0 < n) (hl : 0 ≤ log n) :
    log n ^ (2 * 51) ≤ 6.45 ^ 51 * 7.6890 ^ 100 * n ^ (2 * 51 - 100) := by
  have h1 : log n / 51 ≤ exp (log n / 51 - 1) := by
    have := add_one_le_exp (log n / 51 - 1)
    linarith
  have h2 : (log n / 51) ^ 102 ≤ exp (log n / 51 - 1) ^ 102 :=
    pow_le_pow_left₀ (by positivity) h1 102
  have h3 : exp (log n / 51 - 1) ^ 102 = n ^ 2 / exp 1 ^ 102 := by
    have e1 : exp (log n / 51 - 1) ^ 102 = exp ((102 : ℕ) * (log n / 51 - 1)) :=
      (exp_nat_mul _ 102).symm
    have e2 : ((102 : ℕ) : ℝ) * (log n / 51 - 1) = (2 : ℕ) * log n - (102 : ℕ) * 1 := by
      push_cast
      ring
    rw [e1, e2, exp_sub, exp_nat_mul, exp_nat_mul, exp_log hn0]
  have he : (2.7182818283 : ℝ) ^ 102 ≤ exp 1 ^ 102 :=
    pow_le_pow_left₀ (by norm_num) exp_one_gt_d9.le 102
  have hc : (51 : ℝ) ^ 102 / 2.7182818283 ^ 102 ≤ 6.45 ^ 51 * 7.6890 ^ 100 := by norm_num
  have h4 : log n ^ 102 = 51 ^ 102 * (log n / 51) ^ 102 := by
    rw [div_pow]
    field_simp
  have hn2 : 0 ≤ n ^ 2 := sq_nonneg n
  show log n ^ 102 ≤ 6.45 ^ 51 * 7.6890 ^ 100 * n ^ 2
  calc log n ^ 102 = 51 ^ 102 * (log n / 51) ^ 102 := h4
    _ ≤ 51 ^ 102 * (n ^ 2 / exp 1 ^ 102) := by
        rw [← h3]
        exact mul_le_mul_of_nonneg_left h2 (by positivity)
    _ ≤ 51 ^ 102 * (n ^ 2 / 2.7182818283 ^ 102) := by
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        exact div_le_div_of_nonneg_left hn2 (by positivity) he
    _ = (51 ^ 102 / 2.7182818283 ^ 102) * n ^ 2 := by ring
    _ ≤ (6.45 ^ 51 * 7.6890 ^ 100) * n ^ 2 := mul_le_mul_of_nonneg_right hc hn2

/-- c = 0.55: log^{22} n ≤ 6.45^11 · 0.87160^20 · n^2 for n ≥ n_*. -/
lemma hb_055 {n : ℝ} (hn : (nStar : ℝ) ≤ n) :
    log n ^ (2 * 11) ≤ 6.45 ^ 11 * 0.87160 ^ 20 * n ^ (2 * 11 - 20) := by
  have h := log_pow_le hn 22 2 (by norm_num)
  have hc : (27.62964 : ℝ) ^ 22 / 998616377771 ^ 2 ≤ 6.45 ^ 11 * 0.87160 ^ 20 := by norm_num
  have hn2 : 0 ≤ n ^ 2 := sq_nonneg n
  show log n ^ 22 ≤ 6.45 ^ 11 * 0.87160 ^ 20 * n ^ 2
  calc log n ^ 22 ≤ 27.62964 ^ 22 * (n / 998616377771) ^ 2 := h
    _ = (27.62964 ^ 22 / 998616377771 ^ 2) * n ^ 2 := by ring
    _ ≤ (6.45 ^ 11 * 0.87160 ^ 20) * n ^ 2 := mul_le_mul_of_nonneg_right hc hn2

/-- c = 0.6: log^6 n ≤ 6.45^3 · 0.069829^5 · n for n ≥ n_*. -/
lemma hb_06 {n : ℝ} (hn : (nStar : ℝ) ≤ n) :
    log n ^ (2 * 3) ≤ 6.45 ^ 3 * 0.069829 ^ 5 * n ^ (2 * 3 - 5) := by
  have h := log_pow_le hn 6 1 (by norm_num)
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) (nStar_le hn)
  have hc : (27.62964 : ℝ) ^ 6 / 998616377771 ≤ 6.45 ^ 3 * 0.069829 ^ 5 := by norm_num
  have hn1 : 0 ≤ n ^ 1 := (pow_pos hn0 1).le
  show log n ^ 6 ≤ 6.45 ^ 3 * 0.069829 ^ 5 * n ^ 1
  calc log n ^ 6 ≤ 27.62964 ^ 6 * (n / 998616377771) ^ 1 := h
    _ = (27.62964 ^ 6 / 998616377771) * n ^ 1 := by ring
    _ ≤ (6.45 ^ 3 * 0.069829 ^ 5) * n ^ 1 := mul_le_mul_of_nonneg_right hc hn1

/-- Power form for c = p/q: C d^c ≤ n, n ≥ n_* and log^{2p} n ≤ 6.45^p C^q n^{2p-q} give
d ≤ 6.45 n^2/log^2 n. -/
theorem le_of_power_rat {n d C : ℝ} {p q : ℕ} (hp : p ≠ 0) (hq : q ≠ 0) (hqp : q ≤ 2 * p)
    (hn : (nStar : ℝ) ≤ n) (hd0 : 0 ≤ d) (hC : 0 < C)
    (hb : log n ^ (2 * p) ≤ 6.45 ^ p * C ^ q * n ^ (2 * p - q))
    (hd : C * d ^ ((p : ℝ) / q) ≤ n) : d ≤ 6.45 * n ^ 2 / log n ^ 2 := by
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) (nStar_le hn)
  have hl : 0 < log n := by linarith [log_ge_27 hn]
  exact le_of_pow_le hp (by positivity) (pow_le_of_rpow hq hd0 hC hd)
    (key_of_bound hn0 hC hl hqp hb)

theorem le_of_power_051 {n d : ℝ} (hn : (nStar : ℝ) ≤ n) (hd0 : 0 ≤ d)
    (hd : 7.6890 * d ^ (0.51 : ℝ) ≤ n) : d ≤ 6.45 * n ^ 2 / log n ^ 2 := by
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) (nStar_le hn)
  have hl : 0 ≤ log n := by linarith [log_ge_27 hn]
  refine le_of_power_rat (p := 51) (q := 100) (by norm_num) (by norm_num) (by norm_num) hn hd0
    (by norm_num) (hb_051 hn0 hl) ?_
  have e : ((51 : ℕ) : ℝ) / ((100 : ℕ) : ℝ) = 0.51 := by norm_num
  rw [e]
  exact hd

theorem le_of_power_055 {n d : ℝ} (hn : (nStar : ℝ) ≤ n) (hd0 : 0 ≤ d)
    (hd : 0.87160 * d ^ (0.55 : ℝ) ≤ n) : d ≤ 6.45 * n ^ 2 / log n ^ 2 := by
  refine le_of_power_rat (p := 11) (q := 20) (by norm_num) (by norm_num) (by norm_num) hn hd0
    (by norm_num) (hb_055 hn) ?_
  have e : ((11 : ℕ) : ℝ) / ((20 : ℕ) : ℝ) = 0.55 := by norm_num
  rw [e]
  exact hd

theorem le_of_power_06 {n d : ℝ} (hn : (nStar : ℝ) ≤ n) (hd0 : 0 ≤ d)
    (hd : 0.069829 * d ^ (0.6 : ℝ) ≤ n) : d ≤ 6.45 * n ^ 2 / log n ^ 2 := by
  refine le_of_power_rat (p := 3) (q := 5) (by norm_num) (by norm_num) (by norm_num) hn hd0
    (by norm_num) (hb_06 hn) ?_
  have e : ((3 : ℕ) : ℝ) / ((5 : ℕ) : ℝ) = 0.6 := by norm_num
  rw [e]
  exact hd

variable {γ : ℕ → ℝ} {F : ℕ → ℂ → ℂ}

/-- (iii), power form (0.51, 7.6890). -/
theorem jensen_hyperbolic_power_051 (H : Inputs γ F a0) {n d : ℕ} (hn : nStar ≤ n)
    (hd : 7.6890 * (d : ℝ) ^ (0.51 : ℝ) ≤ n) : Hyperbolic (jensen γ d n) :=
  jensen_hyperbolic_iii H hn (le_of_power_051 (by exact_mod_cast hn) (Nat.cast_nonneg d) hd)

/-- (iii), power form (0.55, 0.87160). -/
theorem jensen_hyperbolic_power_055 (H : Inputs γ F a0) {n d : ℕ} (hn : nStar ≤ n)
    (hd : 0.87160 * (d : ℝ) ^ (0.55 : ℝ) ≤ n) : Hyperbolic (jensen γ d n) :=
  jensen_hyperbolic_iii H hn (le_of_power_055 (by exact_mod_cast hn) (Nat.cast_nonneg d) hd)

/-- (iii), power form (0.6, 0.069829). -/
theorem jensen_hyperbolic_power_06 (H : Inputs γ F a0) {n d : ℕ} (hn : nStar ≤ n)
    (hd : 0.069829 * (d : ℝ) ^ (0.6 : ℝ) ≤ n) : Hyperbolic (jensen γ d n) :=
  jensen_hyperbolic_iii H hn (le_of_power_06 (by exact_mod_cast hn) (Nat.cast_nonneg d) hd)

end

end KimEff


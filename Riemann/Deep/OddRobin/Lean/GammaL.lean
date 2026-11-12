import Riemann.Deep.OddRobin.Lean.LogLike
import Mathlib.NumberTheory.Harmonic.Defs

/-!
# A sharp bracket for Euler's constant (light file)

g plays Mathlib's Real.eulerMascheroniConstant; the two hypotheses are Mathlib's
Real.eulerMascheroniSeq_lt_eulerMascheroniConstant (H_n − log(n+1) < gamma) and
Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' (gamma < H_n − log n, n ≥ 1).
d_n = H_n − L(n + 1/2) decreases and d_n − 1/(6(2n+1)^2) increases; with n = 10 and M = 10^7
this gives 0.5772145 ≤ g ≤ 0.5776.
-/

namespace OddRobin

variable {L : ℝ → ℝ}

/-- d_n = H_n − L(n + 1/2). -/
noncomputable def dSeq (L : ℝ → ℝ) (n : ℕ) : ℝ := (harmonic n : ℝ) - L ((n : ℝ) + 1 / 2)

/-- The correction 1/(24 (n + 1/2)^2) = 1/(6 (2n+1)^2). -/
noncomputable def cSeq (n : ℕ) : ℝ := 1 / (6 * (2 * (n : ℝ) + 1) ^ 2)

lemma tail_ineq (m : ℝ) (hm : 1 < m) :
    2 * (1 / m) ^ 3 / (3 * (1 - (1 / m) ^ 2)) ≤ 1 / (6 * (m - 1) ^ 2) - 1 / (6 * (m + 1) ^ 2) := by
  have h0 : 0 < m := by linarith
  have h1 : 0 < m - 1 := by linarith
  have h2 : 0 < m + 1 := by linarith
  have h3 : 0 < m ^ 2 - 1 := by nlinarith
  have hm0 : m ≠ 0 := h0.ne'
  have hm1 : m - 1 ≠ 0 := h1.ne'
  have hm2 : m + 1 ≠ 0 := h2.ne'
  have hm3 : m ^ 2 - 1 ≠ 0 := h3.ne'
  have e1 : 2 * (1 / m) ^ 3 / (3 * (1 - (1 / m) ^ 2)) = 2 / (3 * m * (m ^ 2 - 1)) := by
    have : (1 : ℝ) - (1 / m) ^ 2 = (m ^ 2 - 1) / m ^ 2 := by field_simp
    rw [this]
    field_simp
  have e2 : 1 / (6 * (m - 1) ^ 2) - 1 / (6 * (m + 1) ^ 2) = 2 * m / (3 * (m ^ 2 - 1) ^ 2) := by
    have h4 : m ^ 2 - 1 = (m - 1) * (m + 1) := by ring
    rw [h4]
    field_simp
    ring
  rw [e1, e2, div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith

lemma step_q (n : ℕ) : ((n : ℝ) + 3 / 2) / ((n : ℝ) + 1 / 2) = 1 + 1 / ((n : ℝ) + 1 / 2) := by
  have h : (n : ℝ) + 1 / 2 ≠ 0 := by positivity
  field_simp
  ring

lemma step_y (n : ℕ) :
    1 / ((n : ℝ) + 1 / 2) / (1 / ((n : ℝ) + 1 / 2) + 2) = 1 / (2 * (n : ℝ) + 2) := by
  have h : (n : ℝ) + 1 / 2 ≠ 0 := by positivity
  have h' : 2 * (n : ℝ) + 2 ≠ 0 := by positivity
  field_simp
  ring

lemma two_div_eq (n : ℕ) : 2 * (1 / (2 * (n : ℝ) + 2)) = 1 / ((n : ℝ) + 1) := by
  have h : (n : ℝ) + 1 ≠ 0 := by positivity
  have h' : 2 * (n : ℝ) + 2 ≠ 0 := by positivity
  field_simp

lemma log_step_ge (hL : LogLike L) (n : ℕ) :
    1 / ((n : ℝ) + 1) ≤ L (((n : ℝ) + 3 / 2) / ((n : ℝ) + 1 / 2)) := by
  have h := hL.two_y_le (a := 1 / ((n : ℝ) + 1 / 2)) (by positivity)
  rw [step_y, two_div_eq] at h
  rw [step_q]
  exact h

lemma cSeq_succ (n : ℕ) : cSeq (n + 1) = 1 / (6 * (2 * (n : ℝ) + 3) ^ 2) := by
  unfold cSeq
  push_cast
  rw [show 2 * ((n : ℝ) + 1) + 1 = 2 * (n : ℝ) + 3 by ring]

lemma log_step_le (hL : LogLike L) (n : ℕ) :
    L (((n : ℝ) + 3 / 2) / ((n : ℝ) + 1 / 2)) ≤ 1 / ((n : ℝ) + 1) + (cSeq n - cSeq (n + 1)) := by
  have h := hL.le_two_y_add (a := 1 / ((n : ℝ) + 1 / 2)) (by positivity)
  rw [step_y, two_div_eq] at h
  rw [step_q]
  have hm : (1 : ℝ) < 2 * (n : ℝ) + 2 := by
    have := n.cast_nonneg (α := ℝ)
    linarith
  have ht := tail_ineq (2 * (n : ℝ) + 2) hm
  rw [show 2 * (n : ℝ) + 2 - 1 = 2 * (n : ℝ) + 1 by ring,
    show 2 * (n : ℝ) + 2 + 1 = 2 * (n : ℝ) + 3 by ring] at ht
  rw [cSeq_succ]
  unfold cSeq
  linarith

lemma dSeq_succ (hL : LogLike L) (n : ℕ) :
    dSeq L (n + 1) = dSeq L n + 1 / ((n : ℝ) + 1) - L (((n : ℝ) + 3 / 2) / ((n : ℝ) + 1 / 2)) := by
  unfold dSeq
  rw [harmonic_succ, hL.map_div (x := (n : ℝ) + 3 / 2) (y := (n : ℝ) + 1 / 2) (by positivity)
    (by positivity)]
  push_cast
  rw [show (n : ℝ) + 1 + 1 / 2 = (n : ℝ) + 3 / 2 by ring]
  ring

lemma dSeq_anti (hL : LogLike L) : Antitone (dSeq L) := by
  refine antitone_nat_of_succ_le (fun n => ?_)
  rw [dSeq_succ hL]
  linarith [log_step_ge hL n]

lemma uSeq_mono (hL : LogLike L) : Monotone (fun n => dSeq L n - cSeq n) := by
  refine monotone_nat_of_le_succ (fun n => ?_)
  show dSeq L n - cSeq n ≤ dSeq L (n + 1) - cSeq (n + 1)
  rw [dSeq_succ hL]
  linarith [log_step_le hL n]

/-- Upper bound: g ≤ d_n + 1/(2M) for n ≤ M, M ≥ 1. -/
lemma g_le (hL : LogLike L) {g : ℝ} (hup : ∀ n : ℕ, 0 < n → g < (harmonic n : ℝ) - L n)
    {n M : ℕ} (hnM : n ≤ M) (hM : 0 < M) : g ≤ dSeq L n + 1 / (2 * (M : ℝ)) := by
  have h1 := hup M hM
  have hM' : (0 : ℝ) < M := by exact_mod_cast hM
  have h2 := hL.le_sub_one (x := ((M : ℝ) + 1 / 2) / M) (by positivity)
  rw [hL.map_div (by positivity) hM'] at h2
  have h3 : ((M : ℝ) + 1 / 2) / M - 1 = 1 / (2 * (M : ℝ)) := by field_simp; ring
  have h4 := dSeq_anti hL hnM
  unfold dSeq at h4 ⊢
  linarith

/-- Lower bound: d_n − c_n − 1/(2M+1) ≤ g for n ≤ M. -/
lemma le_g (hL : LogLike L) {g : ℝ} (hlo : ∀ n : ℕ, (harmonic n : ℝ) - L ((n : ℝ) + 1) < g)
    {n M : ℕ} (hnM : n ≤ M) : dSeq L n - cSeq n - 1 / (2 * (M : ℝ) + 1) ≤ g := by
  have h1 := hlo M
  have h2 := hL.le_sub_one (x := ((M : ℝ) + 1) / ((M : ℝ) + 1 / 2)) (by positivity)
  rw [hL.map_div (by positivity) (by positivity)] at h2
  have h3 : ((M : ℝ) + 1) / ((M : ℝ) + 1 / 2) - 1 = 1 / (2 * (M : ℝ) + 1) := by field_simp; ring
  have h4 : dSeq L n - cSeq n ≤ dSeq L M - cSeq M := uSeq_mono hL hnM
  have h5 : 0 ≤ cSeq M := by unfold cSeq; positivity
  unfold dSeq at h4 ⊢
  linarith

lemma harmonic_ten : harmonic 10 = 7381 / 2520 := by
  norm_num [harmonic, Finset.sum_range_succ]

lemma L_21_2_eq (hL : LogLike L) : L (21 / 2 : ℝ) = 3 * L 2 + L (21 / 16 : ℝ) := by
  have h := hL.two_pow_mul 3 (r := (21 / 16 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 3 * (21 / 16 : ℝ) = 21 / 2 by norm_num] at h
  exact_mod_cast h

theorem gL2_ge (hL : LogLike L) : (6931471805599 / 10000000000000 : ℝ) ≤ L 2 := hL.ge_of (q := (2 : ℝ)) (a := (1 : ℝ)) (by norm_num) (by norm_num) 14 (by norm_num [lsum, Finset.sum_range_succ])

theorem gL2_le (hL : LogLike L) : L 2 ≤ (8664339757 / 12500000000 : ℝ) := hL.le_of (q := (2 : ℝ)) (a := (1 : ℝ)) (by norm_num) (by norm_num) 14 (by norm_num [lsum, ltail, Finset.sum_range_succ])

lemma L_21_16_ge (hL : LogLike L) : (13596685774181 / 50000000000000 : ℝ) ≤ L (21 / 16 : ℝ) := hL.ge_of (q := (21 / 16 : ℝ)) (a := (5 / 16 : ℝ)) (by norm_num) (by norm_num) 7 (by norm_num [lsum, Finset.sum_range_succ])

lemma L_21_16_le (hL : LogLike L) : L (21 / 16 : ℝ) ≤ (5438674309673 / 20000000000000 : ℝ) := hL.le_of (q := (21 / 16 : ℝ)) (a := (5 / 16 : ℝ)) (by norm_num) (by norm_num) 7 (by norm_num [lsum, ltail, Finset.sum_range_succ])

lemma dSeq_ten (hL : LogLike L) : dSeq L 10 = 7381 / 2520 - (3 * L 2 + L (21 / 16 : ℝ)) := by
  rw [← L_21_2_eq hL]
  unfold dSeq
  rw [harmonic_ten]
  norm_num

/-- 0.5772145 ≤ g ≤ 0.5776 from Mathlib's two bounds for gamma (at M = 10^7). -/
theorem gamma_bounds (hL : LogLike L) {g : ℝ}
    (hlo : ∀ n : ℕ, (harmonic n : ℝ) - L ((n : ℝ) + 1) < g)
    (hup : ∀ n : ℕ, 0 < n → g < (harmonic n : ℝ) - L n) :
    (5772145 / 10000000 : ℝ) ≤ g ∧ g ≤ (5776 / 10000 : ℝ) := by
  have h1 := le_g hL hlo (n := 10) (M := 10000000) (by norm_num)
  have h2 := g_le hL hup (n := 10) (M := 10000000) (by norm_num) (by norm_num)
  have hc : cSeq 10 = 1 / 2646 := by unfold cSeq; norm_num
  rw [dSeq_ten hL, hc] at h1
  rw [dSeq_ten hL] at h2
  have a1 := gL2_ge hL
  have a2 := gL2_le hL
  have a3 := L_21_16_ge hL
  have a4 := L_21_16_le hL
  push_cast at h1 h2
  constructor <;> linarith

end OddRobin


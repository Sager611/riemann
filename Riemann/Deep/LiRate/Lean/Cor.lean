import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.SpecificLimits.Normed

/-!
# LiRate corollaries (proofs.md Corollaries 1–4, from Theorem 1, Theorem 2's upper bound, Lemma 0)

`R : ℕ → ℝ` is arbitrary; for the Li remainder `R(n) = λ_n - (n/2)(log n - 1 - log 2π + γ)` the three
inputs are `Thm1Bound R` (Theorem 1, RH), `Thm2Upper R` (Theorem 2's upper bound, RH) and
`Lemma0 R` (¬RH). `y_m = R(m)/m`. `N^{3/2}` is written `N √N`, `(n log n)^{1/4}` is `√√(n log n)`.
-/

open Finset Real Filter

namespace LiRate

/-- `c₁ = (3/5)^{1/2} (4/5)² / (2π√2) = 0.05579…` -/
noncomputable def c₁ : ℝ := √(3 / 5) * (4 / 5) ^ 2 / (2 * π * √2)

/-- Theorem 1 (RH): `∀ ε > 0`, eventually `Σ_{n≤N} R(n)² ≥ (c₁ - ε) N^{3/2} (log N)^{1/2}`. -/
def Thm1Bound (R : ℕ → ℝ) : Prop :=
  ∀ ε > 0, ∀ᶠ N : ℕ in atTop, (c₁ - ε) * (N * √N * √(Real.log N)) ≤ ∑ n ∈ Icc 1 N, R n ^ 2

/-- Lemma 0 (¬RH): `|R(n)| ≥ e^{cn}` for infinitely many `n`, some `c > 0`. -/
def Lemma0 (R : ℕ → ℝ) : Prop := ∃ c > 0, ∃ᶠ n : ℕ in atTop, exp (c * n) ≤ |R n|

/-- Theorem 2's upper bound (RH): `Σ_{N<n≤2N} R(n)² ≪ N^{3/2} (log N)^{3/2}`. -/
def Thm2Upper (R : ℕ → ℝ) : Prop :=
  ∃ C : ℝ, ∃ N₀ : ℕ, ∀ N ≥ N₀,
    ∑ n ∈ Ioc N (2 * N), R n ^ 2 ≤ C * (N * √N) * (Real.log N * √(Real.log N))

lemma c₁_pos : 0 < c₁ := by unfold c₁; positivity

lemma c₁_gt : (0.05569 : ℝ) < c₁ := by
  have ha : (0.7745 : ℝ) < √(3 / 5) := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hb : √2 < (1.4143 : ℝ) := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hb0 : 0 < √(2 : ℝ) := by positivity
  have hp : π * √2 < 3.1416 * 1.4143 := mul_lt_mul'' Real.pi_lt_d4 hb Real.pi_pos.le hb0.le
  unfold c₁
  rw [lt_div_iff₀ (by positivity)]
  nlinarith

/-! ## Corollary 1 -/

lemma step_ineq (x : ℝ) (hx : 0 ≤ x) :
    2 / 3 * (x * √x) + √x ≤ 2 / 3 * ((x + 1) * √(x + 1)) := by
  have ha2 : √x ^ 2 = x := Real.sq_sqrt hx
  have hb2 : √(x + 1) ^ 2 = x + 1 := Real.sq_sqrt (by linarith)
  have hab : √x ≤ √(x + 1) := Real.sqrt_le_sqrt (by linarith)
  set a := √x
  set b := √(x + 1)
  have hq : 0 ≤ 2 - a * b + a ^ 2 := by nlinarith [sq_nonneg (a - b)]
  have key := mul_nonneg (sub_nonneg.2 hab) hq
  have e3 : (b ^ 2 - a ^ 2 - 1) * (2 * b + a) = 0 := by rw [hb2, ha2]; ring
  have e1 : x * a = a ^ 3 := by rw [← ha2]; ring
  have e2 : (x + 1) * b = b ^ 3 := by rw [← hb2]; ring
  rw [e1, e2]
  nlinarith [key, e3]

/-- `Σ_{n=1}^N √n ≤ (2/3) N √N + √N`. -/
lemma sum_sqrt_le (N : ℕ) : ∑ n ∈ Icc 1 N, √(n : ℝ) ≤ 2 / 3 * (N * √N) + √N := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [sum_Icc_succ_top (by omega)]
    have := step_ineq N (Nat.cast_nonneg N)
    push_cast
    linarith

/-- **Corollary 1** (from Theorem 1): `limsup |R(n)|/(n log n)^{1/4} ≥ (3c₁/2)^{1/2}`. -/
theorem cor1 (R : ℕ → ℝ) (h : Thm1Bound R) {c : ℝ} (hc : c < √(3 * c₁ / 2)) :
    ∃ᶠ n : ℕ in atTop, c * √(√(n * Real.log n)) ≤ |R n| := by
  by_cases hc0 : c ≤ 0
  · exact Frequently.of_forall fun n =>
      (mul_nonpos_of_nonpos_of_nonneg hc0 (Real.sqrt_nonneg _)).trans (abs_nonneg _)
  push Not at hc0
  have hc2 : c ^ 2 < 3 * c₁ / 2 := (Real.lt_sqrt hc0.le).1 hc
  have hD : 0 < c₁ - 2 / 3 * c ^ 2 := by linarith
  intro hev
  obtain ⟨n₀, hn₀⟩ := eventually_atTop.1 hev
  set C₀ : ℝ := ∑ n ∈ range n₀, R n ^ 2 with hC₀
  have hC₀0 : 0 ≤ C₀ := sum_nonneg fun n _ => sq_nonneg _
  set K : ℝ := 3 * (C₀ + c ^ 2) / (2 * (c₁ - 2 / 3 * c ^ 2)) with hK
  obtain ⟨N, hN, hNbig⟩ := ((h ((c₁ - 2 / 3 * c ^ 2) / 3) (by linarith)).and
    (eventually_ge_atTop (max (max n₀ 3) (⌈K⌉₊ + 1)))).exists
  have hN3 : 3 ≤ N := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hNbig
  have hNK : K < N := by
    have h1 : (⌈K⌉₊ + 1 : ℕ) ≤ N := le_trans (le_max_right _ _) hNbig
    have h2 : K ≤ ⌈K⌉₊ := Nat.le_ceil K
    have h3 : ((⌈K⌉₊ + 1 : ℕ) : ℝ) ≤ N := by exact_mod_cast h1
    push_cast at h3
    linarith
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
  have hL1 : 1 ≤ Real.log N := by
    rw [Real.le_log_iff_exp_le hNpos]
    have : (3 : ℝ) ≤ N := by exact_mod_cast hN3
    linarith [Real.exp_one_lt_d9]
  have hterm : ∀ n ∈ Icc 1 N, R n ^ 2 ≤
      (if n < n₀ then R n ^ 2 else 0) + c ^ 2 * (√(n : ℝ) * √(Real.log N)) := by
    intro n hn
    rw [mem_Icc] at hn
    split_ifs with hlt
    · have : 0 ≤ c ^ 2 * (√(n : ℝ) * √(Real.log N)) := by positivity
      linarith
    · push Not at hlt
      have h1 := hn₀ n hlt
      push Not at h1
      have hn0 : (0 : ℝ) < n := by exact_mod_cast (by omega : 0 < n)
      have hlogle : Real.log n ≤ Real.log N := Real.log_le_log hn0 (by exact_mod_cast hn.2)
      have hsq : R n ^ 2 ≤ c ^ 2 * √(n * Real.log n) := by
        have hR : |R n| ^ 2 ≤ (c * √(√(n * Real.log n))) ^ 2 := by
          nlinarith [h1, abs_nonneg (R n)]
        rw [sq_abs, mul_pow, Real.sq_sqrt (Real.sqrt_nonneg _)] at hR
        exact hR
      have hs : √(n * Real.log n) ≤ √(n : ℝ) * √(Real.log N) := by
        rw [Real.sqrt_mul hn0.le]
        exact mul_le_mul_of_nonneg_left (Real.sqrt_le_sqrt hlogle) (Real.sqrt_nonneg _)
      have := mul_le_mul_of_nonneg_left hs (sq_nonneg c)
      linarith
  have hsum1 : ∑ n ∈ Icc 1 N, R n ^ 2 ≤
      C₀ + c ^ 2 * √(Real.log N) * ∑ n ∈ Icc 1 N, √(n : ℝ) := by
    have h1 := sum_le_sum hterm
    rw [sum_add_distrib] at h1
    have h2 : ∑ n ∈ Icc 1 N, (if n < n₀ then R n ^ 2 else 0) ≤ C₀ := by
      rw [← sum_filter]
      apply sum_le_sum_of_subset_of_nonneg
      · intro n hn; simp only [mem_filter] at hn; exact mem_range.2 hn.2
      · intro n _ _; exact sq_nonneg _
    have h3 : ∑ n ∈ Icc 1 N, c ^ 2 * (√(n : ℝ) * √(Real.log N)) =
        c ^ 2 * √(Real.log N) * ∑ n ∈ Icc 1 N, √(n : ℝ) := by
      rw [mul_sum]; exact sum_congr rfl fun n _ => by ring
    linarith
  have hS := sum_sqrt_le N
  have hsN : 1 ≤ √(N : ℝ) := Real.one_le_sqrt.2 (by exact_mod_cast (by omega : 1 ≤ N))
  have hsL : 1 ≤ √(Real.log N) := Real.one_le_sqrt.2 hL1
  have hA : (c₁ - (c₁ - 2 / 3 * c ^ 2) / 3) * (N * √N * √(Real.log N)) ≤
      C₀ + c ^ 2 * √(Real.log N) * (2 / 3 * (N * √N) + √N) := by
    have := mul_le_mul_of_nonneg_left hS (by positivity : 0 ≤ c ^ 2 * √(Real.log N))
    linarith
  have hs1 : 1 ≤ √(N : ℝ) * √(Real.log N) := one_le_mul_of_one_le_of_one_le hsN hsL
  have hB : C₀ ≤ C₀ * (√(N : ℝ) * √(Real.log N)) := by
    nlinarith [mul_nonneg hC₀0 (sub_nonneg.2 hs1)]
  have hC : (2 / 3 * (c₁ - 2 / 3 * c ^ 2) * N) * (√(N : ℝ) * √(Real.log N)) ≤
      (C₀ + c ^ 2) * (√(N : ℝ) * √(Real.log N)) := by nlinarith [hA, hB]
  have hE := le_of_mul_le_mul_right hC (by positivity)
  rw [hK, div_lt_iff₀ (by positivity)] at hNK
  nlinarith [hE, hNK]

/-- Corollary 1, numerically: `limsup |R(n)|/(n log n)^{1/4} ≥ 0.289` (RH). -/
theorem cor1_num (R : ℕ → ℝ) (h : Thm1Bound R) :
    ∃ᶠ n : ℕ in atTop, (0.289 : ℝ) * √(√(n * Real.log n)) ≤ |R n| := by
  apply cor1 R h
  rw [Real.lt_sqrt (by norm_num)]
  have := c₁_gt
  norm_num
  linarith

/-- Corollary 1, unconditional (Theorem 1 under RH, Lemma 0 otherwise):
`R(n) = Ω((n log n)^{1/4})`. -/
theorem omega_uncond (R : ℕ → ℝ) (P : Prop) (h1 : P → Thm1Bound R) (h0 : ¬P → Lemma0 R) :
    ∃ c > 0, ∃ᶠ n : ℕ in atTop, c * √(√(n * Real.log n)) ≤ |R n| := by
  by_cases hP : P
  · exact ⟨0.289, by norm_num, cor1_num R (h1 hP)⟩
  · obtain ⟨c, hc, hf⟩ := h0 hP
    refine ⟨1, one_pos, ?_⟩
    have hev : ∀ᶠ n : ℕ in atTop, √(√(n * Real.log n)) ≤ exp (c * n) := by
      filter_upwards [eventually_ge_atTop (⌈1 / c ^ 2⌉₊ + 1)] with n hn
      have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast (by omega : 1 ≤ n)
      have hn0 : (0 : ℝ) < n := by linarith
      have hbig : 1 / c ^ 2 ≤ n := by
        have := Nat.le_ceil (1 / c ^ 2)
        have h2 : ((⌈1 / c ^ 2⌉₊ + 1 : ℕ) : ℝ) ≤ n := by exact_mod_cast hn
        push_cast at h2; linarith
      have h1 : √(√(n * Real.log n)) ≤ √(n : ℝ) := by
        apply Real.sqrt_le_sqrt
        have : n * Real.log n ≤ n * n :=
          mul_le_mul_of_nonneg_left (by linarith [Real.log_le_sub_one_of_pos hn0]) hn0.le
        calc √(n * Real.log n) ≤ √(n * n) := Real.sqrt_le_sqrt this
          _ = n := Real.sqrt_mul_self hn0.le
      have h2 : √(n : ℝ) ≤ c * n := by
        have hs : 1 / c ≤ √(n : ℝ) := by
          rw [Real.le_sqrt (by positivity) hn0.le, div_pow, one_pow]; exact hbig
        have hsn : √(n : ℝ) * √(n : ℝ) = n := Real.mul_self_sqrt hn0.le
        have h3 : 1 ≤ c * √(n : ℝ) := by
          rw [div_le_iff₀ hc] at hs; linarith
        nlinarith [mul_nonneg (sub_nonneg.2 h3) (Real.sqrt_nonneg (n : ℝ))]
      have h3 : c * n ≤ exp (c * n) := by linarith [Real.add_one_le_exp (c * n)]
      linarith
    exact (hf.and_eventually hev).mono fun n hn => by rw [one_mul]; exact hn.2.trans hn.1

/-! ## Corollary 4 (Merlini–Sala–Sala) -/

lemma rpow_quarter {x : ℝ} (hx : 0 ≤ x) : x ^ (1 / 4 : ℝ) = √(√x) := by
  rw [Real.sqrt_eq_rpow, Real.sqrt_eq_rpow, ← Real.rpow_mul hx]; norm_num

lemma quarter_le (n : ℕ) (hn : 3 ≤ n) : √(√(n : ℝ)) ≤ √(√(n * Real.log n)) := by
  have hn3 : (3 : ℝ) ≤ n := by exact_mod_cast hn
  have hL : 1 ≤ Real.log n := by
    rw [Real.le_log_iff_exp_le (by linarith)]; linarith [Real.exp_one_lt_d9]
  exact Real.sqrt_le_sqrt (Real.sqrt_le_sqrt (by nlinarith))

/-- `a log n + C = o((n log n)^{1/4})`. -/
lemma eventually_log_lt (a C c : ℝ) (hc : 0 < c) :
    ∀ᶠ n : ℕ in atTop, a * Real.log n + C < c * √(√(n * Real.log n)) := by
  have hlo := (isLittleO_log_rpow_atTop (by norm_num : (0 : ℝ) < 1 / 4)).bound
    (show 0 < c / 2 / (|a| + 1) by positivity)
  have hbig := (tendsto_rpow_atTop (by norm_num : (0 : ℝ) < 1 / 4)).eventually_gt_atTop
    (2 * |C| / c)
  have hR : ∀ᶠ x : ℝ in atTop, a * Real.log x + C < c * x ^ (1 / 4 : ℝ) := by
    filter_upwards [hlo, hbig, eventually_gt_atTop 0] with x h1 h2 hx
    have hy : 0 < x ^ (1 / 4 : ℝ) := Real.rpow_pos_of_pos hx _
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos hy] at h1
    have h3 : a * Real.log x ≤ |a| * |Real.log x| := by
      rw [← abs_mul]; exact le_abs_self _
    have h4 : |a| * |Real.log x| ≤ |a| * (c / 2 / (|a| + 1) * x ^ (1 / 4 : ℝ)) :=
      mul_le_mul_of_nonneg_left h1 (abs_nonneg a)
    have h5 : |a| * (c / 2 / (|a| + 1) * x ^ (1 / 4 : ℝ)) ≤ c / 2 * x ^ (1 / 4 : ℝ) := by
      have e : |a| * (c / 2 / (|a| + 1) * x ^ (1 / 4 : ℝ)) =
          c / 2 * (|a| / (|a| + 1)) * x ^ (1 / 4 : ℝ) := by ring
      have hq : |a| / (|a| + 1) ≤ 1 := (div_le_one (by positivity)).2 (by linarith)
      rw [e]
      have : c / 2 * (|a| / (|a| + 1)) ≤ c / 2 := by nlinarith
      exact mul_le_mul_of_nonneg_right this hy.le
    have h6 : C < c / 2 * x ^ (1 / 4 : ℝ) := by
      rw [div_lt_iff₀ hc] at h2; linarith [le_abs_self C]
    linarith
  filter_upwards [tendsto_natCast_atTop_atTop.eventually hR, eventually_ge_atTop 3] with n h1 h2
  rw [rpow_quarter (Nat.cast_nonneg n)] at h1
  nlinarith [mul_le_mul_of_nonneg_left (quarter_le n h2) hc.le]

/-- `C' n^ε + C = o((n log n)^{1/4})` for `ε < 1/4`. -/
lemma eventually_rpow_lt (C' C c ε : ℝ) (hc : 0 < c) (hε : ε < 1 / 4) :
    ∀ᶠ n : ℕ in atTop, C' * (n : ℝ) ^ ε + C < c * √(√(n * Real.log n)) := by
  have hz := (tendsto_rpow_neg_atTop (by linarith : (0 : ℝ) < 1 / 4 - ε)).eventually
    (gt_mem_nhds (show (0 : ℝ) < c / 2 / (|C'| + 1) by positivity))
  have hbig := (tendsto_rpow_atTop (by norm_num : (0 : ℝ) < 1 / 4)).eventually_gt_atTop
    (2 * |C| / c)
  have hR : ∀ᶠ x : ℝ in atTop, C' * x ^ ε + C < c * x ^ (1 / 4 : ℝ) := by
    filter_upwards [hz, hbig, eventually_gt_atTop 0] with x h1 h2 hx
    have hy : 0 < x ^ (1 / 4 : ℝ) := Real.rpow_pos_of_pos hx _
    have hsplit : x ^ ε = x ^ (-(1 / 4 - ε)) * x ^ (1 / 4 : ℝ) := by
      rw [← Real.rpow_add hx]; congr 1; ring
    have h3 : C' * x ^ ε ≤ |C'| * x ^ ε :=
      mul_le_mul_of_nonneg_right (le_abs_self C') (Real.rpow_pos_of_pos hx _).le
    have h4 : |C'| * x ^ ε ≤ c / 2 * x ^ (1 / 4 : ℝ) := by
      have hq : |C'| / (|C'| + 1) ≤ 1 := (div_le_one (by positivity)).2 (by linarith)
      have h7 : |C'| * x ^ (-(1 / 4 - ε)) ≤ c / 2 :=
        calc |C'| * x ^ (-(1 / 4 - ε)) ≤ |C'| * (c / 2 / (|C'| + 1)) :=
              mul_le_mul_of_nonneg_left h1.le (abs_nonneg _)
          _ = c / 2 * (|C'| / (|C'| + 1)) := by ring
          _ ≤ c / 2 := by nlinarith
      rw [hsplit]
      calc |C'| * (x ^ (-(1 / 4 - ε)) * x ^ (1 / 4 : ℝ)) =
            (|C'| * x ^ (-(1 / 4 - ε))) * x ^ (1 / 4 : ℝ) := by ring
        _ ≤ c / 2 * x ^ (1 / 4 : ℝ) := mul_le_mul_of_nonneg_right h7 hy.le
    have h6 : C < c / 2 * x ^ (1 / 4 : ℝ) := by
      rw [div_lt_iff₀ hc] at h2; linarith [le_abs_self C]
    linarith
  filter_upwards [tendsto_natCast_atTop_atTop.eventually hR, eventually_ge_atTop 3] with n h1 h2
  rw [rpow_quarter (Nat.cast_nonneg n)] at h1
  nlinarith [mul_le_mul_of_nonneg_left (quarter_le n h2) hc.le]

/-- **Corollary 4** (refutation of Merlini–Sala–Sala, arXiv 1903.06233 §2): if
`λ̃_n - R(n)` is bounded (`λ̃_n = R(n) - 3/4 + O(1/n)`), then `|λ̃_n| ≤ a log n` fails
infinitely often, for every `a`. Unconditional (RH: Theorem 1; ¬RH: Lemma 0). -/
theorem cor4 (R lt : ℕ → ℝ) (P : Prop) (h1 : P → Thm1Bound R) (h0 : ¬P → Lemma0 R)
    (hlt : ∃ C, ∀ n, |lt n - R n| ≤ C) (a : ℝ) :
    ¬ ∀ᶠ n : ℕ in atTop, |lt n| ≤ a * Real.log n := by
  intro hev
  obtain ⟨c, hc, hf⟩ := omega_uncond R P h1 h0
  obtain ⟨C, hC⟩ := hlt
  obtain ⟨n, hn1, hn2, hn3⟩ := (hf.and_eventually (hev.and (eventually_log_lt a C c hc))).exists
  have h4 : |R n| ≤ |lt n| + C := by
    have := abs_sub_abs_le_abs_sub (R n) (lt n)
    rw [abs_sub_comm] at this
    linarith [hC n]
  linarith

/-- Corollary 4, second form: `λ̃_n = O(n^ε)` fails for every `ε < 1/4`. -/
theorem cor4_rpow (R lt : ℕ → ℝ) (P : Prop) (h1 : P → Thm1Bound R) (h0 : ¬P → Lemma0 R)
    (hlt : ∃ C, ∀ n, |lt n - R n| ≤ C) {ε : ℝ} (hε : ε < 1 / 4) (C' : ℝ) :
    ¬ ∀ᶠ n : ℕ in atTop, |lt n| ≤ C' * (n : ℝ) ^ ε := by
  intro hev
  obtain ⟨c, hc, hf⟩ := omega_uncond R P h1 h0
  obtain ⟨C, hC⟩ := hlt
  obtain ⟨n, hn1, hn2, hn3⟩ :=
    (hf.and_eventually (hev.and (eventually_rpow_lt C' C c ε hc hε))).exists
  have h4 : |R n| ≤ |lt n| + C := by
    have := abs_sub_abs_le_abs_sub (R n) (lt n)
    rw [abs_sub_comm] at this
    linarith [hC n]
  linarith

end LiRate

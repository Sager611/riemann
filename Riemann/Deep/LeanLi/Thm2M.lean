module

public import Riemann.Deep.LeanLi.PairCorr

/-!
# LeanLi: Theorem 2 with Montgomery (RH): `Σ_{N<n≤2N} R(n)² ≪ N^{3/2} (log N)^{1/2} log log N`

Smooth majorant `ψ(n/N)` of `1_{(N, 2N]}` (bump on `(1/2, 5/2)`, `= 1` on `[1, 2]`):
`Σ ψ W_T² ≤ Σ_{γ,γ'} |Σ_n ψ(n/N) e^{in(θ_γ - θ_γ')}| ≪ N Σ_{γ,γ'} (1 + log T |γ - γ'|)^{-2} ≪ N T log T` (Montgomery),
`Σ (2na)² ≪ N² log log T/T` (Bessel, Selberg), `Rest = O(log N)`.
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory Topology
open scoped ContDiff

namespace LeanLi

/-- The bump `ψ`: `1` on `[1, 2]`, supported in `(1/2, 5/2)`. -/
noncomputable def psiB : ContDiffBump ((3 : ℝ) / 2) := ⟨1 / 2, 1, by norm_num, by norm_num⟩

noncomputable def psi (y : ℝ) : ℝ := psiB y

lemma psi_nonneg (y : ℝ) : 0 ≤ psi y := psiB.nonneg
lemma psi_le_one (y : ℝ) : psi y ≤ 1 := psiB.le_one

lemma psi_eq_one {y : ℝ} (h1 : 1 ≤ y) (h2 : y ≤ 2) : psi y = 1 := by
  apply psiB.one_of_mem_closedBall
  rw [Metric.mem_closedBall, Real.dist_eq, abs_le]
  show -(1 / 2 : ℝ) ≤ y - 3 / 2 ∧ y - 3 / 2 ≤ 1 / 2
  constructor <;> linarith

lemma psi_supp {y : ℝ} (h : psi y ≠ 0) : 1 / 2 < y ∧ y < 5 / 2 := by
  have : y ∈ Function.support psiB := h
  rw [psiB.support_eq, Metric.mem_ball, Real.dist_eq, abs_lt] at this
  have h1 := this.1; have h2 := this.2
  change -(1 : ℝ) < y - 3 / 2 at h1
  change y - 3 / 2 < (1 : ℝ) at h2
  constructor <;> linarith

lemma testFn_psi : TestFn psi 3 where
  smooth := psiB.contDiff
  compact := psiB.hasCompactSupport
  pos := fun y h => by linarith [(psi_supp h).1]
  le := fun y h => by have := (psi_supp h).2; push_cast; linarith

/-- `log N ≤ 8 T` at `T = √(N/log N)`. -/
lemma log_le_eight_Tf {N : ℕ} (hN : 3 ≤ N) : Real.log N ≤ 8 * Tf N := by
  have hN0 : (0 : ℝ) < N := by
    have h3 : (3 : ℝ) ≤ N := by exact_mod_cast hN
    linarith
  have hl := log_nat_ge_one hN
  set q := √(√(N : ℝ)) with hq
  have hs0 : 0 < √(N : ℝ) := Real.sqrt_pos.2 hN0
  have hq0 : 0 < q := Real.sqrt_pos.2 hs0
  have hq2 : q ^ 2 = √(N : ℝ) := Real.sq_sqrt hs0.le
  have hq4 : q ^ 4 = N := by
    rw [show q ^ 4 = (q ^ 2) ^ 2 by ring, hq2, Real.sq_sqrt hN0.le]
  have hq1 : 1 ≤ q := by
    by_contra h; push Not at h
    have : q ^ 4 < 1 := by
      calc q ^ 4 < 1 ^ 4 := pow_lt_pow_left₀ h hq0.le (by norm_num)
        _ = 1 := one_pow _
    rw [hq4] at this
    have : (3 : ℝ) ≤ N := by exact_mod_cast hN
    linarith
  have h1 : Real.log N ≤ 4 * q := by
    have e : Real.log N = 2 * Real.log √(N : ℝ) := by rw [Real.log_sqrt hN0.le]; ring
    have := log_le_two_sqrt hs0
    rw [e]; linarith
  have h3 : Real.log N ^ 3 ≤ 64 * N := by
    have := pow_le_pow_left₀ (by linarith) h1 3
    have hq34 : q ^ 3 ≤ q ^ 4 := pow_le_pow_right₀ hq1 (by norm_num)
    calc Real.log N ^ 3 ≤ (4 * q) ^ 3 := this
      _ = 64 * q ^ 3 := by ring
      _ ≤ 64 * q ^ 4 := by linarith
      _ = 64 * N := by rw [hq4]
  unfold Tf
  have : Real.log N / 8 ≤ √((N : ℝ) / Real.log N) := by
    rw [Real.le_sqrt (by positivity) (by positivity), le_div_iff₀ (by linarith)]
    nlinarith
  linarith

lemma log_log_ge_one {N : ℕ} (hN : 16 ≤ N) : 1 ≤ Real.log (Real.log N) := by
  have h16 : (16 : ℝ) ≤ N := by exact_mod_cast hN
  have hl : Real.log 16 ≤ Real.log N := Real.log_le_log (by norm_num) h16
  have h2 : Real.log 16 = 4 * Real.log 2 := by
    rw [show (16 : ℝ) = 2 ^ 4 by norm_num, Real.log_pow]; push_cast; ring
  have hl2 := Real.log_two_gt_d9
  have he := Real.exp_one_lt_d9
  rw [Real.le_log_iff_exp_le (by linarith)]
  linarith

section
variable (L : LiteratureInputs)
include L

/-- Smoothed mean square of `W_T` under Montgomery: `Σ ψ(n/N) W_T(n)² ≤ C N T log N` at `T = √(N/log N)`. -/
lemma sumW_psi (hRH : RiemannHypothesis) : ∃ C : ℝ, ∃ N₀ : ℕ, 3 ≤ N₀ ∧ 0 ≤ C ∧ ∀ N : ℕ, N₀ ≤ N →
    ∑ n ∈ range (3 * N + 1), psi (n / N) * Wv (Tf N) n ^ 2 ≤ C * (N * Tf N * Real.log N) := by
  obtain ⟨CE, hCE0, hCE⟩ := Ej_decay testFn_psi 0 2
  obtain ⟨KP, TP, hKP0, hTP3, hKP⟩ := pair_sum_bound L hRH
  obtain ⟨N₁, hN₁⟩ : ∃ N₁ : ℕ, 4 * TP ^ 4 ≤ N₁ := exists_nat_ge _
  refine ⟨CE * KP / 2, max N₁ 3, le_max_right _ _, by positivity, fun N hN => ?_⟩
  have hN3 : 3 ≤ N := le_trans (le_max_right _ _) hN
  have hN1 : 1 ≤ N := by omega
  have hN1r : (1 : ℝ) ≤ N := by exact_mod_cast hN1
  have hNN₁ : (N₁ : ℝ) ≤ N := by exact_mod_cast le_trans (le_max_left _ _) hN
  have hℓN := log_nat_ge_one hN3
  set T := Tf N with hT
  have hTP : TP ≤ T := Tf_ge (by linarith) hN3 (by linarith)
  have hT3 : 3 ≤ T := by linarith
  have hT2 : T ^ 2 = N / Real.log N := Tf_sq hN3
  have hlogT := log_Tf_le hN3
  have hlT0 : 0 < Real.log T := Real.log_pos (by linarith)
  -- (1) W² ≤ Σ cos, swap sums
  have h1 : ∑ n ∈ range (3 * N + 1), psi (n / N) * Wv T n ^ 2 ≤
      ∑ j ∈ ZK T, ∑ k ∈ ZK T, ‖Ej psi 3 N 0 (thK j - thK k)‖ := by
    calc ∑ n ∈ range (3 * N + 1), psi (n / N) * Wv T n ^ 2 ≤
          ∑ n ∈ range (3 * N + 1), psi (n / N) * ∑ j ∈ ZK T, ∑ k ∈ ZK T, Real.cos (n * (thK j - thK k)) :=
          Finset.sum_le_sum fun n _ => mul_le_mul_of_nonneg_left (Wv_sq_le T n) (psi_nonneg _)
      _ = ∑ j ∈ ZK T, ∑ k ∈ ZK T, (Ej psi 3 N 0 (thK j - thK k)).re := by
          simp only [Ej_re, floor_cast_mul, pow_zero, mul_one, Finset.mul_sum]
          rw [Finset.sum_comm]
          refine Finset.sum_congr rfl fun j _ => ?_
          rw [Finset.sum_comm]
      _ ≤ _ := Finset.sum_le_sum fun j _ => Finset.sum_le_sum fun k _ => Complex.re_le_norm _
  -- (2) decay and separation
  have hsepN : N / (T ^ 2 + 1 / 4) ≥ Real.log T := by
    rw [ge_iff_le, le_div_iff₀ (by positivity)]
    have h4 : T ^ 2 + 1 / 4 ≤ 2 * T ^ 2 := by nlinarith
    have h5 : Real.log T * (2 * T ^ 2) ≤ N := by
      rw [hT2]; field_simp; nlinarith
    calc Real.log T * (T ^ 2 + 1 / 4) ≤ Real.log T * (2 * T ^ 2) := mul_le_mul_of_nonneg_left h4 hlT0.le
      _ ≤ N := h5
  have h2 : ∀ j ∈ ZK T, ∀ k ∈ ZK T, ‖Ej psi 3 N 0 (thK j - thK k)‖ ≤
      CE * N * (1 / (1 + Real.log T * |ordK j - ordK k|) ^ 2) := by
    intro j hj k hk
    have hj0 : 0 < ordK j := j.1.2.2
    have hk0 : 0 < ordK k := k.1.2.2
    have hjT : ordK j ≤ T := mem_ZK.1 hj
    have hkT : ordK k ≤ T := mem_ZK.1 hk
    have hxπ : |thK j - thK k| ≤ Real.pi := by
      have := liTheta_pos hj0; have := liTheta_lt_pi hj0
      have := liTheta_pos hk0; have := liTheta_lt_pi hk0
      unfold thK; rw [abs_le]; constructor <;> linarith
    have hd := hCE N hN1r (thK j - thK k) hxπ
    have hsep := liTheta_sep hj0.le hjT hk0.le hkT
    have hNsep : Real.log T * |ordK j - ordK k| ≤ N * |thK j - thK k| := by
      have hw : 0 < T ^ 2 + 1 / 4 := by positivity
      calc Real.log T * |ordK j - ordK k| ≤ N / (T ^ 2 + 1 / 4) * |ordK j - ordK k| :=
            mul_le_mul_of_nonneg_right hsepN (abs_nonneg _)
        _ = N * (|ordK j - ordK k| / (T ^ 2 + 1 / 4)) := by ring
        _ ≤ N * |thK j - thK k| := mul_le_mul_of_nonneg_left hsep (by linarith)
    have hpos : 0 < (1 + Real.log T * |ordK j - ordK k|) ^ 2 := by positivity
    have hle : (1 + Real.log T * |ordK j - ordK k|) ^ 2 ≤ (1 + N * |thK j - thK k|) ^ 2 :=
      pow_le_pow_left₀ (by positivity) (by linarith) 2
    rw [zero_add, pow_one] at hd
    have h3 : ‖Ej psi 3 N 0 (thK j - thK k)‖ * (1 + Real.log T * |ordK j - ordK k|) ^ 2 ≤ CE * N :=
      le_trans (mul_le_mul_of_nonneg_left hle (norm_nonneg _)) hd
    rw [mul_one_div, le_div_iff₀ hpos]
    exact h3
  have h3 := hKP T hTP
  calc ∑ n ∈ range (3 * N + 1), psi (n / N) * Wv T n ^ 2 ≤
        ∑ j ∈ ZK T, ∑ k ∈ ZK T, ‖Ej psi 3 N 0 (thK j - thK k)‖ := h1
    _ ≤ ∑ j ∈ ZK T, ∑ k ∈ ZK T, CE * N * (1 / (1 + Real.log T * |ordK j - ordK k|) ^ 2) :=
        Finset.sum_le_sum fun j hj => Finset.sum_le_sum fun k hk => h2 j hj k hk
    _ = CE * N * ∑ j ∈ ZK T, ∑ k ∈ ZK T, 1 / (1 + Real.log T * |ordK j - ordK k|) ^ 2 := by
        simp only [Finset.mul_sum]
    _ ≤ CE * N * (KP * (T * Real.log T)) := mul_le_mul_of_nonneg_left h3 (by positivity)
    _ ≤ CE * N * (KP * (T * (Real.log N / 2))) := by
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        apply mul_le_mul_of_nonneg_left _ hKP0
        exact mul_le_mul_of_nonneg_left hlogT (by linarith)
    _ = CE * KP / 2 * (N * T * Real.log N) := by ring

omit L in
lemma aCoef_zero (T X : ℝ) : aCoef T X 0 = 0 := by
  unfold aCoef; simp

set_option maxHeartbeats 4000000 in
/-- **Theorem 2 with Montgomery** (RH): `Σ_{N<n≤2N} R(n)² ≤ C N^{3/2} (log N)^{1/2} log log N`. -/
theorem thm2_mont (hRH : RiemannHypothesis) : ∃ C : ℝ, ∀ᶠ N : ℕ in atTop,
    ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤ C * ((N : ℝ) * Tf N * Real.log N * Real.log (Real.log N)) := by
  obtain ⟨CW, NW, hNW3, hCW0, hW⟩ := sumW_psi L hRH
  obtain ⟨CS, TS, hCS0, hTSe, hS⟩ := selberg_tail L
  obtain ⟨K, NR, hNR3, hR⟩ := rest_bound L
  obtain ⟨N₁, hN₁⟩ : ∃ N₁ : ℕ, 4 * TS ^ 4 ≤ N₁ := exists_nat_ge _
  have hπ := Real.pi_gt_three
  have hπ0 := Real.pi_pos
  refine ⟨12 * CW + 162 * Real.pi * CS + 96 * K ^ 2, ?_⟩
  filter_upwards [eventually_ge_atTop (max (max NW NR) (max N₁ 16))] with N hN
  have hNW : NW ≤ N := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hN
  have hNR : NR ≤ N := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hN
  have hN16 : 16 ≤ N := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hN
  have hNN₁ : (N₁ : ℝ) ≤ N := by exact_mod_cast le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hN
  have hN3 : 3 ≤ N := by omega
  have hN1 : 1 ≤ N := by omega
  have hN3r : (3 : ℝ) ≤ N := by exact_mod_cast hN3
  have hℓ := log_nat_ge_one hN3
  have hll := log_log_ge_one hN16
  have hlogT := log_Tf_le hN3
  have hT2 := Tf_sq hN3
  have hτT : TS ≤ Tf N := Tf_ge (le_trans (Real.exp_pos 1).le hTSe) hN3 (by linarith)
  have hTsq := Tf_le_sqrt hN3
  have h8T := log_le_eight_Tf hN3
  set ℓ := Real.log N with hℓdef
  set T := Tf N with hTdef
  set X := (N : ℝ) ^ 8 with hXdef
  have hT1 : 1 ≤ T := by
    have := Real.exp_one_gt_d9
    have : Real.exp 1 ≤ T := le_trans hTSe hτT
    linarith
  have hT0 : 0 < T := by linarith
  have hT2ℓ : T ^ 2 * ℓ = N := by rw [hT2]; field_simp
  have hTX : T ≤ X := by
    have : √(N : ℝ) ≤ N := by
      rw [Real.sqrt_le_left (by linarith)]; nlinarith
    have : (N : ℝ) ≤ X := by rw [hXdef]; exact le_self_pow₀ (by linarith) (by norm_num)
    linarith
  -- majorize by the smooth sum
  have hsub : Finset.Ioc N (2 * N) ⊆ range (3 * N + 1) := by
    intro n hn; rw [Finset.mem_Ioc] at hn; rw [Finset.mem_range]; omega
  have hmaj : ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤ ∑ n ∈ range (3 * N + 1), psi (n / N) * R n ^ 2 := by
    have e : ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 = ∑ n ∈ Finset.Ioc N (2 * N), psi (n / N) * R n ^ 2 := by
      refine Finset.sum_congr rfl fun n hn => ?_
      rw [Finset.mem_Ioc] at hn
      have hN0 : (0 : ℝ) < N := by linarith
      have h1 : (1 : ℝ) ≤ n / N := by
        rw [le_div_iff₀ hN0, one_mul]; exact_mod_cast hn.1.le
      have h2 : (n : ℝ) / N ≤ 2 := by
        rw [div_le_iff₀ hN0]; exact_mod_cast hn.2
      rw [psi_eq_one h1 h2, one_mul]
    rw [e]
    exact Finset.sum_le_sum_of_subset_of_nonneg hsub fun n _ _ => mul_nonneg (psi_nonneg _) (sq_nonneg _)
  -- pointwise decomposition
  have hpt : ∀ n ∈ range (3 * N + 1), psi (n / N) * R n ^ 2 ≤
      12 * (psi (n / N) * Wv T n ^ 2) + 108 * (N : ℝ) ^ 2 * aCoef T X n ^ 2 + 3 * (K * ℓ) ^ 2 := by
    intro n hn
    have hnN : (n : ℝ) ≤ 3 * N := by
      have := Finset.mem_range.1 hn
      have : n ≤ 3 * N := by omega
      exact_mod_cast this
    by_cases hψ : psi (n / N) = 0
    · rw [hψ]; simp only [zero_mul, mul_zero, zero_add]; positivity
    have hsupp := psi_supp hψ
    have hN0 : (0 : ℝ) < N := by linarith
    have hn1 : (N : ℝ) ≤ 2 * n := by
      have := hsupp.1; rw [lt_div_iff₀ hN0] at this; linarith
    have hn2 : (n : ℝ) ≤ 3 * N := hnN
    have hd := decomp L hRH n hT0 hTX
    have hr := hR N hNR n hn1 hn2
    have hr2 : Rest n T X ^ 2 ≤ (K * ℓ) ^ 2 := by
      calc Rest n T X ^ 2 = |Rest n T X| ^ 2 := (sq_abs _).symm
        _ ≤ (K * ℓ) ^ 2 := pow_le_pow_left₀ (abs_nonneg _) hr 2
    have hR2 : R n ^ 2 ≤ 3 * (4 * Wv T n ^ 2 + (2 * n * aCoef T X n) ^ 2 + Rest n T X ^ 2) := by
      rw [hd]
      nlinarith [sq_nonneg (-2 * Wv T n - 2 * n * aCoef T X n), sq_nonneg (-2 * Wv T n - Rest n T X),
        sq_nonneg (2 * n * aCoef T X n - Rest n T X)]
    have ha2 : (2 * n * aCoef T X n) ^ 2 ≤ 36 * (N : ℝ) ^ 2 * aCoef T X n ^ 2 := by
      have h4 : (n : ℝ) ^ 2 ≤ 9 * N ^ 2 := by nlinarith [Nat.cast_nonneg (α := ℝ) n]
      have : (2 * n * aCoef T X n) ^ 2 = 4 * (n : ℝ) ^ 2 * aCoef T X n ^ 2 := by ring
      rw [this]; nlinarith [sq_nonneg (aCoef T X n)]
    have hψ0 := psi_nonneg (n / N)
    have hψ1 := psi_le_one (n / N)
    have hb : psi (n / N) * R n ^ 2 ≤ psi (n / N) * (3 * (4 * Wv T n ^ 2 + 36 * (N : ℝ) ^ 2 * aCoef T X n ^ 2 +
        (K * ℓ) ^ 2)) := mul_le_mul_of_nonneg_left (by linarith) hψ0
    have hA : 0 ≤ 36 * (N : ℝ) ^ 2 * aCoef T X n ^ 2 := by positivity
    have hB : 0 ≤ (K * ℓ) ^ 2 := sq_nonneg _
    have e1 : psi (n / N) * (36 * (N : ℝ) ^ 2 * aCoef T X n ^ 2) ≤ 36 * (N : ℝ) ^ 2 * aCoef T X n ^ 2 :=
      mul_le_of_le_one_left hA hψ1
    have e2 : psi (n / N) * (K * ℓ) ^ 2 ≤ (K * ℓ) ^ 2 := mul_le_of_le_one_left hB hψ1
    nlinarith
  have hsum := Finset.sum_le_sum hpt
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, Finset.sum_const,
    Finset.card_range, nsmul_eq_mul] at hsum
  -- the three pieces
  have hWb := hW N hNW
  have hfilt : ∑ n ∈ range (3 * N + 1), aCoef T X n ^ 2 =
      ∑ n ∈ (range (3 * N + 1)).filter (fun n => 1 ≤ n), aCoef T X n ^ 2 := by
    refine (Finset.sum_filter_of_ne fun n _ hne => ?_).symm
    by_contra h; push Not at h
    have : n = 0 := by omega
    rw [this, aCoef_zero] at hne; simp at hne
  have hBa := bessel_a hT0.le hTX ((range (3 * N + 1)).filter (fun n => 1 ≤ n)) (fun n hn => (Finset.mem_filter.1 hn).2)
  have hSb := hS T X hτT hTX
  have hlT : Real.log T ≤ ℓ := by linarith
  have hllT : Real.log (Real.log T) ≤ Real.log ℓ := by
    have h1 : 0 < Real.log T := by
      have : 1 ≤ Real.log T := by rw [← Real.log_exp 1]; exact Real.log_le_log (Real.exp_pos 1) (le_trans hTSe hτT)
      linarith
    exact Real.log_le_log h1 hlT
  have hSfin : ∫ t in T..X, zetaS t ^ 2 / (t ^ 2 + 1 / 4) ≤ CS * (3 * Real.log ℓ) / T :=
    hSb.trans (div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left (by linarith) hCS0) hT0.le)
  have hNT2 : (N : ℝ) ^ 2 / T = N * T * ℓ := by
    rw [div_eq_iff hT0.ne']; linear_combination (-(N : ℝ)) * hT2ℓ
  have hafin : 108 * (N : ℝ) ^ 2 * ∑ n ∈ range (3 * N + 1), aCoef T X n ^ 2 ≤
      162 * Real.pi * CS * (N * T * ℓ * Real.log ℓ) := by
    rw [hfilt]
    have h1 : 108 * (N : ℝ) ^ 2 * ∑ n ∈ (range (3 * N + 1)).filter (fun n => 1 ≤ n), aCoef T X n ^ 2 ≤
        108 * (N : ℝ) ^ 2 * (Real.pi / 2 * (CS * (3 * Real.log ℓ) / T)) :=
      mul_le_mul_of_nonneg_left (hBa.trans (mul_le_mul_of_nonneg_left hSfin (by positivity))) (by positivity)
    refine h1.trans (le_of_eq ?_)
    have : 108 * (N : ℝ) ^ 2 * (Real.pi / 2 * (CS * (3 * Real.log ℓ) / T)) =
        162 * Real.pi * CS * Real.log ℓ * ((N : ℝ) ^ 2 / T) := by
      field_simp; ring
    rw [this, hNT2]; ring
  have hRfin : ((3 * N + 1 : ℕ) : ℝ) * (3 * (K * ℓ) ^ 2) ≤ 96 * K ^ 2 * (N * T * ℓ * Real.log ℓ) := by
    push_cast
    have h1 : (3 * (N : ℝ) + 1) ≤ 4 * N := by linarith
    have h2 : (3 * (N : ℝ) + 1) * (3 * (K * ℓ) ^ 2) ≤ 4 * N * (3 * (K * ℓ) ^ 2) :=
      mul_le_mul_of_nonneg_right h1 (by positivity)
    have h3 : 4 * (N : ℝ) * (3 * (K * ℓ) ^ 2) = 12 * K ^ 2 * N * ℓ * ℓ := by ring
    have h4 : 12 * K ^ 2 * N * ℓ * ℓ ≤ 12 * K ^ 2 * N * ℓ * (8 * T) :=
      mul_le_mul_of_nonneg_left h8T (by positivity)
    have h5 : 12 * K ^ 2 * N * ℓ * (8 * T) = 96 * K ^ 2 * (N * T * ℓ) := by ring
    have h6 : 96 * K ^ 2 * (N * T * ℓ) ≤ 96 * K ^ 2 * (N * T * ℓ * Real.log ℓ) := by
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      exact le_mul_of_one_le_right (by positivity) hll
    linarith
  have hWfin : 12 * ∑ n ∈ range (3 * N + 1), psi (n / N) * Wv T n ^ 2 ≤ 12 * CW * (N * T * ℓ * Real.log ℓ) := by
    have h1 : CW * (N * T * ℓ) ≤ CW * (N * T * ℓ * Real.log ℓ) := by
      apply mul_le_mul_of_nonneg_left _ hCW0
      exact le_mul_of_one_le_right (by positivity) hll
    linarith
  calc ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤ ∑ n ∈ range (3 * N + 1), psi (n / N) * R n ^ 2 := hmaj
    _ ≤ _ := hsum
    _ ≤ 12 * CW * (N * T * ℓ * Real.log ℓ) + 162 * Real.pi * CS * (N * T * ℓ * Real.log ℓ) +
        96 * K ^ 2 * (N * T * ℓ * Real.log ℓ) := by linarith
    _ = (12 * CW + 162 * Real.pi * CS + 96 * K ^ 2) * ((N : ℝ) * T * ℓ * Real.log ℓ) := by ring

end

end LeanLi

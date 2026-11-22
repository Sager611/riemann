module

public import Riemann.Deep.LeanLi.MeanW

/-!
# LeanLi: Theorem 2 (RH), unconditional form: `Σ_{N<n≤2N} R(n)² ≪ N^{3/2} (log N)^{3/2}`

Lemma A at height `T = √(N/log N)`, `X = N⁸`: `R = -2W_T + 2n a + Rest` with
`Σ W_T² ≪ N·N(T) log T` (local counts), `Σ (2na)² ≪ N² ∫_T^∞ S²/t² ≪ N² log log T/T` (Bessel, Selberg),
`Rest = O(log N)`.
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory

namespace LeanLi

/-- The height `T = √(N/log N)`. -/
noncomputable def Tf (N : ℕ) : ℝ := √((N : ℝ) / Real.log N)

lemma log_nat_ge_one {N : ℕ} (hN : 3 ≤ N) : 1 ≤ Real.log N := by
  have h3 : (3 : ℝ) ≤ N := by exact_mod_cast hN
  rw [Real.le_log_iff_exp_le (by linarith)]
  have := Real.exp_one_lt_d9; linarith

lemma Tf_pos {N : ℕ} (hN : 3 ≤ N) : 0 < Tf N := by
  have := log_nat_ge_one hN
  have h3 : (3 : ℝ) ≤ N := by exact_mod_cast hN
  exact Real.sqrt_pos.2 (by positivity)

lemma Tf_sq {N : ℕ} (hN : 3 ≤ N) : Tf N ^ 2 = N / Real.log N := by
  have := log_nat_ge_one hN
  exact Real.sq_sqrt (by positivity)

lemma Tf_le_sqrt {N : ℕ} (hN : 3 ≤ N) : Tf N ≤ √(N : ℝ) :=
  Real.sqrt_le_sqrt (div_le_self (Nat.cast_nonneg _) (log_nat_ge_one hN))

lemma log_Tf_le {N : ℕ} (hN : 3 ≤ N) : Real.log (Tf N) ≤ Real.log N / 2 := by
  have h := Real.log_le_log (Tf_pos hN) (Tf_le_sqrt hN)
  rwa [Real.log_sqrt (Nat.cast_nonneg _)] at h

lemma Tf_ge {τ : ℝ} (hτ : 0 ≤ τ) {N : ℕ} (hN3 : 3 ≤ N) (hN : 4 * τ ^ 4 ≤ N) : τ ≤ Tf N := by
  have hN0 : (0 : ℝ) < N := by
    have h3 : (3 : ℝ) ≤ N := by exact_mod_cast hN3
    linarith
  have hl := log_le_two_sqrt hN0
  have hl1 := log_nat_ge_one hN3
  unfold Tf
  rw [Real.le_sqrt hτ (by positivity), le_div_iff₀ (by linarith)]
  have hs : 2 * τ ^ 2 ≤ √(N : ℝ) := by
    rw [Real.le_sqrt (by positivity) hN0.le]; nlinarith
  have hsN : √(N : ℝ) * √(N : ℝ) = N := Real.mul_self_sqrt hN0.le
  have h1 : τ ^ 2 * Real.log N ≤ τ ^ 2 * (2 * √(N : ℝ)) := mul_le_mul_of_nonneg_left hl (by positivity)
  have h2 : 2 * τ ^ 2 * √(N : ℝ) ≤ √(N : ℝ) * √(N : ℝ) := mul_le_mul_of_nonneg_right hs (Real.sqrt_nonneg _)
  nlinarith

/-- `N^{3/2} (log N)^{3/2} = N T (log N)²` at `T = √(N/log N)`. -/
lemma Tf_id {N : ℕ} (hN : 3 ≤ N) :
    (N : ℝ) * √(N : ℝ) * (Real.log N * √(Real.log N)) = N * Tf N * Real.log N ^ 2 := by
  have hl := log_nat_ge_one hN
  have hq : 0 < √(Real.log N) := Real.sqrt_pos.2 (by linarith)
  obtain ⟨q, hq0, hℓq⟩ : ∃ q : ℝ, 0 < q ∧ Real.log N = q ^ 2 := ⟨√(Real.log N), hq, (Real.sq_sqrt (by linarith)).symm⟩
  have hsq : √(Real.log N) = q := by rw [hℓq, Real.sqrt_sq hq0.le]
  unfold Tf
  rw [Real.sqrt_div' _ (by linarith), hsq, hℓq]
  field_simp

/-- Bessel for `a_{T,X}(n)`. -/
lemma bessel_a {T X : ℝ} (hT : 0 ≤ T) (hTX : T ≤ X) (s : Finset ℕ) (hs : ∀ n ∈ s, 1 ≤ n) :
    ∑ n ∈ s, aCoef T X n ^ 2 ≤ Real.pi / 2 * ∫ t in T..X, zetaS t ^ 2 / (t ^ 2 + 1 / 4) := by
  have hc : Continuous fun t : ℝ => 1 / (t ^ 2 + 1 / 4) := continuous_inv_sq_add
  have hh : IntervalIntegrable (fun t => zetaS t / (t ^ 2 + 1 / 4)) volume T X := by
    have := (intervalIntegrable_zetaS T X).mul_continuousOn hc.continuousOn
    refine this.congr_uIoo fun t _ => ?_
    simp only [div_eq_mul_one_div (zetaS t)]
  have e : ∀ t, (zetaS t / (t ^ 2 + 1 / 4)) ^ 2 * (t ^ 2 + 1 / 4) = zetaS t ^ 2 / (t ^ 2 + 1 / 4) := by
    intro t
    have : 0 < t ^ 2 + 1 / 4 := by positivity
    field_simp
  have hh2 : IntervalIntegrable (fun t => (zetaS t / (t ^ 2 + 1 / 4)) ^ 2 * (t ^ 2 + 1 / 4)) volume T X := by
    have := (intervalIntegrable_zetaS_sq T X).mul_continuousOn hc.continuousOn
    refine this.congr_uIoo fun t _ => ?_
    rw [e]; simp only [div_eq_mul_one_div (zetaS t ^ 2)]
  have hb := bessel_t hT hTX (fun t => zetaS t / (t ^ 2 + 1 / 4)) hh hh2 s hs
  simp only [e] at hb
  exact hb

section
variable (L : LitCore)
include L

set_option maxHeartbeats 2000000 in
/-- `Rest = O(log N)` for `N/2 ≤ n ≤ 3N` at `T = √(N/log N)`, `X = N⁸`. -/
lemma rest_bound : ∃ K : ℝ, ∃ N₀ : ℕ, 3 ≤ N₀ ∧ ∀ N : ℕ, N₀ ≤ N → ∀ n : ℕ, (N : ℝ) ≤ 2 * n →
    (n : ℝ) ≤ 3 * N → |Rest n (Tf N) ((N : ℝ) ^ 8)| ≤ K * Real.log N := by
  obtain ⟨C₁, T₁, hT₁, hNA⟩ := NA_bound L
  obtain ⟨Kr, hKr⟩ := rfun_L1 L
  obtain ⟨c₀, hc₀⟩ := cos_log_bound
  obtain ⟨Kt, hKt0, hKt⟩ := li_tail_bound L
  obtain ⟨N₁, hN₁⟩ : ∃ N₁ : ℕ, 4 * T₁ ^ 4 ≤ N₁ := exists_nat_ge _
  have hπ := Real.pi_gt_three
  have hπ0 := Real.pi_pos
  refine ⟨34 * |C₁| + 16 / Real.pi + (2 + 6 * |Kr| + 3 + 18 / Real.pi + (|c₀| + 48) / Real.pi + 9 * Kt + 90),
    max N₁ 3, le_max_right _ _, fun N hN n hn1 hn2 => ?_⟩
  have hN3 : 3 ≤ N := le_trans (le_max_right _ _) hN
  have hNN₁ : (N₁ : ℝ) ≤ N := by exact_mod_cast le_trans (le_max_left _ _) hN
  have hN3r : (3 : ℝ) ≤ N := by exact_mod_cast hN3
  have hℓ := log_nat_ge_one hN3
  set ℓ := Real.log N with hℓdef
  set T := Tf N with hTdef
  set X := (N : ℝ) ^ 8 with hXdef
  have hT₁T : T₁ ≤ T := Tf_ge (by linarith) hN3 (by linarith)
  have hT1 : 1 ≤ T := by linarith
  have hT0 : 0 < T := by linarith
  have hlogT0 : 0 ≤ Real.log T := Real.log_nonneg hT1
  have hlogT : Real.log T ≤ ℓ := by have := log_Tf_le hN3; linarith
  have hT2 : T ^ 2 = N / ℓ := Tf_sq hN3
  have hT2N : T ^ 2 ≤ N := by rw [hT2]; exact div_le_self (by linarith) hℓ
  have hNX : (N : ℝ) ≤ X := by
    rw [hXdef]; exact le_self_pow₀ (by linarith) (by norm_num)
  have hT4 : T₁ ≤ T₁ ^ 4 := le_self_pow₀ (by linarith) (by norm_num)
  have hT₁X : T₁ ≤ X := by linarith
  have hTX : T ≤ X := by
    have : T ≤ √(N : ℝ) := Tf_le_sqrt hN3
    have : √(N : ℝ) ≤ N := by
      rw [Real.sqrt_le_left (by linarith)]; nlinarith
    linarith
  have hX0 : 0 < X := by linarith
  have hX9 : 9 ≤ X := by
    have : (N : ℝ) ^ 2 ≤ X := by rw [hXdef]; exact pow_le_pow_right₀ (by linarith) (by norm_num)
    nlinarith
  have hπ4 := Real.pi_lt_d2
  have hlogX : Real.log X = 8 * ℓ := by rw [hXdef, Real.log_pow]; push_cast; ring
  have hsX : √X = (N : ℝ) ^ 4 := by
    rw [hXdef, show (N : ℝ) ^ 8 = ((N : ℝ) ^ 4) ^ 2 by ring, Real.sqrt_sq (by positivity)]
  have hn0 : (0 : ℝ) < n := by linarith
  have hn1' : 1 ≤ n := by exact_mod_cast (show (1 : ℝ) ≤ n by linarith)
  have hnsq : (n : ℝ) ^ 2 ≤ 9 * N ^ 2 := by nlinarith
  have hN2 : (N : ℝ) ^ 2 ≤ (N : ℝ) ^ 4 := pow_le_pow_right₀ (by linarith) (by norm_num)
  have hN2' : (N : ℝ) ^ 4 ≤ X ^ 2 := by
    rw [hXdef]; rw [← pow_mul]; exact pow_le_pow_right₀ (by linarith) (by norm_num)
  -- the eight terms
  have t1 : |2 * Real.cos (n * liTheta T) * ((zeroCount T : ℝ) - Afun T)| ≤ 2 * (|C₁| * ℓ) := by
    rw [abs_mul, abs_mul, abs_two]
    have h2 := (hNA T hT₁T).trans ((mul_le_mul_of_nonneg_right (le_abs_self C₁) hlogT0).trans
      (mul_le_mul_of_nonneg_left hlogT (abs_nonneg _)))
    calc 2 * |Real.cos (n * liTheta T)| * |(zeroCount T : ℝ) - Afun T| ≤ 2 * 1 * (|C₁| * ℓ) :=
          mul_le_mul (mul_le_mul_of_nonneg_left (abs_cos_le_one _) (by norm_num)) h2 (abs_nonneg _) (by norm_num)
      _ = 2 * (|C₁| * ℓ) := by ring
  have hr0 : |∫ t in (0 : ℝ)..T, rfun t| ≤ |Kr| :=
    (intervalIntegral.abs_integral_le_integral_abs hT0.le).trans ((hKr T hT0.le).trans (le_abs_self _))
  have t3 : |∫ t in T..X, gfun n t * rfun t| ≤ 4 * |Kr| := by
    have hgc : Continuous (gfun n) := by unfold gfun; have := continuous_liTheta; fun_prop
    have hi := (intervalIntegrable_rfun hT0.le hX0.le).continuousOn_mul hgc.continuousOn
    refine (intervalIntegral.abs_integral_le_integral_abs hTX).trans ?_
    have hm : ∫ t in T..X, |gfun n t * rfun t| ≤ ∫ t in T..X, 4 * |rfun t| := by
      apply intervalIntegral.integral_mono_on hTX hi.abs ((intervalIntegrable_rfun hT0.le hX0.le).abs.const_mul 4)
      intro t _
      rw [abs_mul]
      have : |gfun n t| ≤ 4 := by
        unfold gfun
        have := Real.cos_le_one (n * liTheta t); have := Real.neg_one_le_cos (n * liTheta t)
        rw [abs_le]; constructor <;> linarith
      exact mul_le_mul_of_nonneg_right this (abs_nonneg _)
    refine hm.trans ?_
    rw [intervalIntegral.integral_const_mul]
    have h1 : ∫ t in T..X, |rfun t| ≤ ∫ t in (0 : ℝ)..X, |rfun t| :=
      intervalIntegral.integral_mono_interval hT0.le hTX le_rfl
        (Filter.Eventually.of_forall fun t => abs_nonneg (rfun t)) (intervalIntegrable_rfun le_rfl hX0.le).abs
    have h2 := hKr X hX0.le
    have : Kr ≤ |Kr| := le_abs_self _
    linarith
  have t4 : |Iclosed n - liTrend n| ≤ 3 := Iclosed_sub_trend_le n
  have t5 : |∫ t in Ioi X, hI n t| ≤ 18 / Real.pi := by
    have h := I_tail_bound n (show 2 * Real.pi ≤ X by linarith)
    refine h.trans ?_
    rw [hsX, div_le_div_iff₀ (by positivity) (by positivity)]
    have := mul_le_mul_of_nonneg_left hnsq hπ0.le
    have := mul_le_mul_of_nonneg_left hN2 hπ0.le
    nlinarith
  have t6 : |1 / Real.pi * (∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi)))| ≤
      (|c₀| + 16 * ℓ + 48) / Real.pi := by
    rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 1 / Real.pi)]
    have h := hc₀ n hn1' T hT1
    have h2 : 4 * (T ^ 2 + 1) * (Real.log T + 3) / n ≤ 16 * ℓ + 48 := by
      rw [div_le_iff₀ hn0]
      have hA : T ^ 2 + 1 ≤ 2 * N := by linarith
      have hB : Real.log T + 3 ≤ ℓ + 3 := by linarith
      have hAB : (T ^ 2 + 1) * (Real.log T + 3) ≤ (2 * N) * (ℓ + 3) :=
        mul_le_mul hA hB (by linarith) (by linarith)
      nlinarith
    have h3 : |∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))| ≤ |c₀| + 16 * ℓ + 48 := by
      have := le_abs_self c₀; linarith
    calc 1 / Real.pi * |∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))| ≤
          1 / Real.pi * (|c₀| + 16 * ℓ + 48) := mul_le_mul_of_nonneg_left h3 (by positivity)
      _ = (|c₀| + 16 * ℓ + 48) / Real.pi := by ring
  have t7 : |gfun n X * ((zeroCount X : ℝ) - Afun X)| ≤ 32 * (|C₁| * ℓ) := by
    rw [abs_mul]
    have hg : |gfun n X| ≤ 4 := by
      unfold gfun
      have := Real.cos_le_one (n * liTheta X); have := Real.neg_one_le_cos (n * liTheta X)
      rw [abs_le]; constructor <;> linarith
    have h2 := (hNA X hT₁X).trans (mul_le_mul_of_nonneg_right (le_abs_self C₁) (by rw [hlogX]; linarith))
    rw [hlogX] at h2
    calc |gfun n X| * |(zeroCount X : ℝ) - Afun X| ≤ 4 * (|C₁| * (8 * ℓ)) :=
          mul_le_mul hg h2 (abs_nonneg _) (by norm_num)
      _ = 32 * (|C₁| * ℓ) := by ring
  have t8 : |∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k)| ≤ 9 * Kt + 90 := by
    have hnn : 0 ≤ ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k) := tsum_nonneg fun k => by
      unfold gfun; have := Real.cos_le_one (n * liTheta (ordK k)); linarith
    rw [abs_of_nonneg hnn]
    have h : ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k) ≤ (n : ℝ) ^ 2 * (Kt / X ^ 2 + 10 / √X) :=
      hKt n X hX0
    refine h.trans ?_
    rw [hsX, mul_add]
    have e1 : (n : ℝ) ^ 2 * (Kt / X ^ 2) ≤ 9 * Kt := by
      rw [mul_div_assoc', div_le_iff₀ (pow_pos hX0 2)]
      have h1 : (n : ℝ) ^ 2 ≤ 9 * X ^ 2 := by linarith
      have := mul_le_mul_of_nonneg_right h1 hKt0
      linarith
    have e2 : (n : ℝ) ^ 2 * (10 / (N : ℝ) ^ 4) ≤ 90 := by
      rw [mul_div_assoc', div_le_iff₀ (pow_pos (by linarith) 4)]
      linarith
    linarith
  -- assemble
  unfold Rest
  set A1 := 2 * Real.cos (n * liTheta T) * ((zeroCount T : ℝ) - Afun T)
  set A2 := ∫ t in (0 : ℝ)..T, rfun t
  set A3 := ∫ t in T..X, gfun n t * rfun t
  set A4 := Iclosed n - liTrend n
  set A5 := ∫ t in Ioi X, hI n t
  set A6 := 1 / Real.pi * (∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi)))
  set A7 := gfun n X * ((zeroCount X : ℝ) - Afun X)
  set A8 := ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k)
  have hC0 : 0 ≤ |C₁| * ℓ := by positivity
  have hp : 0 < 1 / Real.pi := by positivity
  have e16 : 16 / Real.pi * ℓ = 16 * ℓ / Real.pi := by ring
  have hconst : 0 ≤ 2 + 6 * |Kr| + 3 + 18 / Real.pi + (|c₀| + 48) / Real.pi + 9 * Kt + 90 := by positivity
  have hK : (34 * |C₁| + 16 / Real.pi + (2 + 6 * |Kr| + 3 + 18 / Real.pi + (|c₀| + 48) / Real.pi + 9 * Kt + 90)) * ℓ ≥
      34 * (|C₁| * ℓ) + 16 * ℓ / Real.pi + (2 + 6 * |Kr| + 3 + 18 / Real.pi + (|c₀| + 48) / Real.pi + 9 * Kt + 90) := by
    have : (2 + 6 * |Kr| + 3 + 18 / Real.pi + (|c₀| + 48) / Real.pi + 9 * Kt + 90) ≤
        (2 + 6 * |Kr| + 3 + 18 / Real.pi + (|c₀| + 48) / Real.pi + 9 * Kt + 90) * ℓ := le_mul_of_one_le_right hconst hℓ
    linarith
  have e6 : (|c₀| + 16 * ℓ + 48) / Real.pi = 16 * ℓ / Real.pi + (|c₀| + 48) / Real.pi := by ring
  rw [abs_le] at t1 hr0 t3 t4 t5 t6 t7 t8 ⊢
  constructor <;> linarith [t1.1, t1.2, hr0.1, hr0.2, t3.1, t3.2, t4.1, t4.2, t5.1, t5.2, t6.1, t6.2, t7.1, t7.2, t8.1, t8.2]

set_option maxHeartbeats 4000000 in
/-- **Theorem 2** (RH): `Σ_{N<n≤2N} R(n)² ≪ N^{3/2} (log N)^{3/2}`. -/
theorem thm2_upper (hRH : RiemannHypothesis) : LR.Thm2Upper R := by
  obtain ⟨CW, TW, hTW, hCW0, hW⟩ := sumW_bound L
  obtain ⟨CS, TS, hCS0, hTSe, hS⟩ := selberg_tail L
  obtain ⟨Tc, hTc, hc⟩ := count_le L
  obtain ⟨K, NR, hNR3, hR⟩ := rest_bound L
  set τ := max TW (max TS Tc) with hτ
  have hτ4 : 4 ≤ τ := le_trans hTW (le_max_left _ _)
  obtain ⟨N₁, hN₁⟩ : ∃ N₁ : ℕ, 4 * τ ^ 4 ≤ N₁ := exists_nat_ge _
  have hπ := Real.pi_gt_three
  have hπ0 := Real.pi_pos
  refine ⟨12 * (CW / 2 + 6 * Real.pi * CW) + 24 * Real.pi * (3 * CS) + 3 * K ^ 2, max N₁ NR, fun N hN => ?_⟩
  have hNR : NR ≤ N := le_trans (le_max_right _ _) hN
  have hN3 : 3 ≤ N := le_trans hNR3 hNR
  have hN3r : (3 : ℝ) ≤ N := by exact_mod_cast hN3
  have hNN₁ : (N₁ : ℝ) ≤ N := by exact_mod_cast le_trans (le_max_left _ _) hN
  have hℓ := log_nat_ge_one hN3
  have hP := Tf_id hN3
  have hlogT := log_Tf_le hN3
  have hT2 := Tf_sq hN3
  have hτT : τ ≤ Tf N := Tf_ge (by linarith) hN3 (by linarith)
  have hTsq := Tf_le_sqrt hN3
  set ℓ := Real.log N with hℓdef
  set T := Tf N with hTdef
  set X := (N : ℝ) ^ 8 with hXdef
  have hTWT : TW ≤ T := le_trans (le_max_left _ _) hτT
  have hTST : TS ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hτT
  have hTcT : Tc ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hτT
  have hT1 : 1 ≤ T := by linarith
  have hT0 : 0 < T := by linarith
  have hlogT0 : 0 ≤ Real.log T := Real.log_nonneg hT1
  have hT2ℓ : T ^ 2 * ℓ = N := by rw [hT2]; field_simp
  have hTX : T ≤ X := by
    have : √(N : ℝ) ≤ N := by
      rw [Real.sqrt_le_left (by linarith)]; nlinarith
    have : (N : ℝ) ≤ X := by rw [hXdef]; exact le_self_pow₀ (by linarith) (by norm_num)
    linarith
  -- pointwise decomposition
  have hdec : ∀ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤
      3 * (4 * Wv T n ^ 2 + 16 * (N : ℝ) ^ 2 * aCoef T X n ^ 2 + (K * ℓ) ^ 2) := by
    intro n hn
    rw [Finset.mem_Ioc] at hn
    have hn1 : (N : ℝ) < n := by exact_mod_cast hn.1
    have hn2 : (n : ℝ) ≤ 2 * N := by exact_mod_cast hn.2
    have hd := decomp L hRH n hT0 hTX
    have hr := hR N hNR n (by linarith) (by linarith)
    rw [hd]
    have hr2 : Rest n T X ^ 2 ≤ (K * ℓ) ^ 2 := by
      have h0 := abs_nonneg (Rest n T X)
      calc Rest n T X ^ 2 = |Rest n T X| ^ 2 := (sq_abs _).symm
        _ ≤ (K * ℓ) ^ 2 := pow_le_pow_left₀ h0 hr 2
    have ha2 : (2 * n * aCoef T X n) ^ 2 ≤ 16 * (N : ℝ) ^ 2 * aCoef T X n ^ 2 := by
      have h4 : (n : ℝ) ^ 2 ≤ 4 * N ^ 2 := by nlinarith
      have h0 := sq_nonneg (aCoef T X n)
      have : (2 * n * aCoef T X n) ^ 2 = 4 * (n : ℝ) ^ 2 * aCoef T X n ^ 2 := by ring
      rw [this]
      nlinarith
    nlinarith [sq_nonneg (-2 * Wv T n - 2 * n * aCoef T X n), sq_nonneg (-2 * Wv T n - Rest n T X),
      sq_nonneg (2 * n * aCoef T X n - Rest n T X)]
  have hcard : (Finset.Ioc N (2 * N)).card = N := by rw [Nat.card_Ioc]; omega
  have hsum : ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤
      3 * (4 * ∑ n ∈ Finset.Ioc N (2 * N), Wv T n ^ 2 + 16 * (N : ℝ) ^ 2 * ∑ n ∈ Finset.Ioc N (2 * N), aCoef T X n ^ 2 +
        N * (K * ℓ) ^ 2) := by
    refine (Finset.sum_le_sum hdec).trans (le_of_eq ?_)
    rw [← Finset.mul_sum, Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
      Finset.sum_const, hcard, nsmul_eq_mul]
  -- the three pieces
  have hWb := hW T hTWT N
  have hBa := bessel_a hT0.le hTX (Finset.Ioc N (2 * N)) (fun n hn => by
    have := (Finset.mem_Ioc.1 hn).1; omega)
  have hSb := hS T X hTST hTX
  have hNT : (zeroCount T : ℝ) ≤ T * (ℓ / 2) := (hc T hTcT).trans (mul_le_mul_of_nonneg_left hlogT hT0.le)
  have hNT0 : 0 ≤ (zeroCount T : ℝ) := Nat.cast_nonneg _
  have hB : CW * Real.log T ≤ CW * (ℓ / 2) := mul_le_mul_of_nonneg_left hlogT hCW0
  have hB0 : 0 ≤ CW * Real.log T := by positivity
  have hQ : Real.pi * (T ^ 2 + 1 / 4) ≤ 2 * Real.pi * T ^ 2 := by
    have h1 : (1 : ℝ) ≤ T ^ 2 := by nlinarith
    have := mul_le_mul_of_nonneg_left (show T ^ 2 + 1 / 4 ≤ 2 * T ^ 2 by linarith) hπ0.le
    linarith
  have hH0 : 0 ≤ 1 + Real.log (T + 1) := by have := Real.log_nonneg (show (1 : ℝ) ≤ T + 1 by linarith); linarith
  have hH : 1 + Real.log (T + 1) ≤ 3 * ℓ := by
    have h1 : Real.log (T + 1) ≤ Real.log (2 * T) := Real.log_le_log (by linarith) (by linarith)
    rw [Real.log_mul (by norm_num) hT0.ne'] at h1
    have h2 : Real.log 2 ≤ 1 := by
      have := Real.log_le_sub_one_of_pos (show (0 : ℝ) < 2 by norm_num); linarith
    linarith
  have hin : 2 * (CW * Real.log T) * N + 2 * (Real.pi * (T ^ 2 + 1 / 4)) *
      (2 * (CW * Real.log T) * (1 + Real.log (T + 1))) ≤ CW * ℓ * N + 12 * Real.pi * CW * T ^ 2 * ℓ ^ 2 := by
    have h1 : 2 * (CW * Real.log T) * N ≤ CW * ℓ * N := by
      have := mul_le_mul_of_nonneg_right hB (show (0 : ℝ) ≤ N by linarith); linarith
    have h2 : (CW * Real.log T) * (1 + Real.log (T + 1)) ≤ (CW * (ℓ / 2)) * (3 * ℓ) :=
      mul_le_mul hB hH hH0 (by positivity)
    have h3 : Real.pi * (T ^ 2 + 1 / 4) * ((CW * Real.log T) * (1 + Real.log (T + 1))) ≤
        (2 * Real.pi * T ^ 2) * ((CW * (ℓ / 2)) * (3 * ℓ)) :=
      mul_le_mul hQ h2 (by positivity) (by positivity)
    linarith
  have hin0 : 0 ≤ 2 * (CW * Real.log T) * N + 2 * (Real.pi * (T ^ 2 + 1 / 4)) *
      (2 * (CW * Real.log T) * (1 + Real.log (T + 1))) := by positivity
  have hWfin : ∑ n ∈ Finset.Ioc N (2 * N), Wv T n ^ 2 ≤ (CW / 2 + 6 * Real.pi * CW) * (N * T * ℓ ^ 2) := by
    refine hWb.trans ((mul_le_mul hNT hin hin0 (by positivity)).trans (le_of_eq ?_))
    linear_combination (6 * Real.pi * CW * T * ℓ ^ 2) * hT2ℓ
  have hll : Real.log (Real.log T) + 2 ≤ 3 * ℓ := by
    have hl1 : 1 ≤ Real.log T := by
      rw [← Real.log_exp 1]; exact Real.log_le_log (Real.exp_pos 1) (le_trans hTSe hTST)
    have := Real.log_le_sub_one_of_pos (show 0 < Real.log T by linarith)
    linarith
  have hSfin : ∫ t in T..X, zetaS t ^ 2 / (t ^ 2 + 1 / 4) ≤ CS * (3 * ℓ) / T :=
    hSb.trans (div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left hll hCS0) hT0.le)
  have hNT2 : (N : ℝ) ^ 2 / T = N * T * ℓ := by
    rw [div_eq_iff hT0.ne']; linear_combination (-(N : ℝ)) * hT2ℓ
  have hafin : 16 * (N : ℝ) ^ 2 * ∑ n ∈ Finset.Ioc N (2 * N), aCoef T X n ^ 2 ≤ 8 * Real.pi * (3 * CS) * (N * T * ℓ ^ 2) := by
    have h1 : 16 * (N : ℝ) ^ 2 * ∑ n ∈ Finset.Ioc N (2 * N), aCoef T X n ^ 2 ≤
        16 * (N : ℝ) ^ 2 * (Real.pi / 2 * (CS * (3 * ℓ) / T)) :=
      mul_le_mul_of_nonneg_left (hBa.trans (mul_le_mul_of_nonneg_left hSfin (by positivity))) (by positivity)
    refine h1.trans (le_of_eq ?_)
    have : 16 * (N : ℝ) ^ 2 * (Real.pi / 2 * (CS * (3 * ℓ) / T)) = 24 * Real.pi * CS * ℓ * ((N : ℝ) ^ 2 / T) := by
      field_simp; ring
    rw [this, hNT2]; ring
  have hRfin : (N : ℝ) * (K * ℓ) ^ 2 ≤ K ^ 2 * (N * T * ℓ ^ 2) := by
    have : (N : ℝ) * (K * ℓ) ^ 2 = K ^ 2 * (N * 1 * ℓ ^ 2) := by ring
    rw [this]
    apply mul_le_mul_of_nonneg_left _ (sq_nonneg K)
    apply mul_le_mul_of_nonneg_right _ (sq_nonneg ℓ)
    exact mul_le_mul_of_nonneg_left hT1 (by linarith)
  rw [mul_assoc _ ((N : ℝ) * √(N : ℝ)), hP]
  calc ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤ _ := hsum
    _ ≤ 3 * (4 * ((CW / 2 + 6 * Real.pi * CW) * (N * T * ℓ ^ 2)) + 8 * Real.pi * (3 * CS) * (N * T * ℓ ^ 2) +
        K ^ 2 * (N * T * ℓ ^ 2)) := by linarith
    _ = (12 * (CW / 2 + 6 * Real.pi * CW) + 24 * Real.pi * (3 * CS) + 3 * K ^ 2) * (N * T * ℓ ^ 2) := by ring

end

end LeanLi

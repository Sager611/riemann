module

public import Riemann.Deep.LeanLi.Thm4

/-!
# LeanLi: Theorem 4 (RH), Claim part 4: the smoothed mean of `R(n)` is `7/4`
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory Topology
open scoped ContDiff

namespace LeanLi

section
variable (L : LiteratureInputs)
include L

/-- The zero term `∫_{T₀}^X S G = O(log X)`. -/
lemma zero_term {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) : ∃ K₁ K₂ X₀ : ℝ, 1 ≤ X₀ ∧ 0 ≤ K₁ ∧ 0 ≤ K₂ ∧
    ∀ N : ℕ, 1 ≤ N → ∀ T₀ X : ℝ, 0 ≤ T₀ → T₀ ≤ X → X₀ ≤ X →
      |∫ t in T₀..X, zetaS t * Gf φ B N t| ≤ K₁ + K₂ * Real.log X := by
  obtain ⟨S₁, S₂, X₁, hX₁, hS₁0, hS₂0, hS⟩ := S1_sup L
  obtain ⟨CG, hCG0, hG⟩ := Gf_bound h
  obtain ⟨Cd, hCd0, hd⟩ := Gd_bound h
  refine ⟨(2 * CG + Cd / 2) * S₁, (2 * CG + Cd / 2) * S₂, X₁, hX₁, by positivity, by positivity,
    fun N hN T₀ X hT₀ hTX hX => ?_⟩
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  rw [ibp_S hT₀ hTX (Gf φ B N) (Gd φ B N) (hasDerivAt_Gf φ B N) (continuous_Gd φ B N)]
  set M := S₁ + S₂ * Real.log X with hM
  have hM0 : 0 ≤ M := by have := Real.log_nonneg (show (1 : ℝ) ≤ X by linarith); positivity
  have hSX := hS X hX X (by linarith) le_rfl
  have hST := hS X hX T₀ hT₀ hTX
  have hGX := hG N hN X (by linarith)
  have hGT := hG N hN T₀ hT₀
  have h1 : |S1 X * Gf φ B N X| ≤ M * CG := by rw [abs_mul]; exact mul_le_mul hSX hGX (abs_nonneg _) hM0
  have h2 : |S1 T₀ * Gf φ B N T₀| ≤ M * CG := by rw [abs_mul]; exact mul_le_mul hST hGT (abs_nonneg _) hM0
  have hi : IntervalIntegrable (fun t => S1 t * Gd φ B N t) volume T₀ X :=
    (continuous_S1.mul (continuous_Gd φ B N)).intervalIntegrable _ _
  have hk0 : IntervalIntegrable (fun t => (N : ℝ) ^ 2 / (t + N) ^ 3) volume T₀ X := by
    apply ContinuousOn.intervalIntegrable
    intro t ht
    rw [Set.uIcc_of_le hTX] at ht
    have hp : 0 < t + N := by linarith [ht.1]
    exact (continuousAt_const.div ((continuousAt_id.add continuousAt_const).pow 3)
      (pow_ne_zero 3 hp.ne')).continuousWithinAt
  have hk : IntervalIntegrable (fun t => M * (Cd * ((N : ℝ) ^ 2 / (t + N) ^ 3))) volume T₀ X :=
    (hk0.const_mul Cd).const_mul M
  have h3 : |∫ t in T₀..X, S1 t * Gd φ B N t| ≤ M * (Cd / 2) := by
    refine (intervalIntegral.abs_integral_le_integral_abs hTX).trans ?_
    have hm : ∫ t in T₀..X, |S1 t * Gd φ B N t| ≤ ∫ t in T₀..X, M * (Cd * ((N : ℝ) ^ 2 / (t + N) ^ 3)) := by
      apply intervalIntegral.integral_mono_on hTX hi.abs hk
      intro t ht
      rw [abs_mul]
      have ht0 : 0 ≤ t := le_trans hT₀ ht.1
      have := hd N hN t ht0
      rw [mul_div_assoc] at this
      exact mul_le_mul (hS X hX t ht0 ht.2) this (abs_nonneg _) hM0
    refine hm.trans ?_
    rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul]
    have h4 := integral_Nsq_div_cube_le hT₀ hTX hN1
    have h5 : Cd * ∫ t in T₀..X, (N : ℝ) ^ 2 / (t + N) ^ 3 ≤ Cd * (1 / 2) := mul_le_mul_of_nonneg_left h4 hCd0
    have h6 := mul_le_mul_of_nonneg_left h5 hM0
    linarith
  calc |S1 X * Gf φ B N X - S1 T₀ * Gf φ B N T₀ - ∫ t in T₀..X, S1 t * Gd φ B N t| ≤
        |S1 X * Gf φ B N X| + |S1 T₀ * Gf φ B N T₀| + |∫ t in T₀..X, S1 t * Gd φ B N t| := by
        have := abs_sub (S1 X * Gf φ B N X - S1 T₀ * Gf φ B N T₀) (∫ t in T₀..X, S1 t * Gd φ B N t)
        have := abs_sub (S1 X * Gf φ B N X) (S1 T₀ * Gf φ B N T₀)
        linarith
    _ ≤ M * CG + M * CG + M * (Cd / 2) := by linarith
    _ = (2 * CG + Cd / 2) * S₁ + (2 * CG + Cd / 2) * S₂ * Real.log X := by rw [hM]; ring

/-- The `r` term `∫_{T₀}^X F₀ r = O(log X)`. -/
lemma r_term {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) : ∃ K₁ K₂ X₀ : ℝ, 1 ≤ X₀ ∧ 0 ≤ K₁ ∧ 0 ≤ K₂ ∧
    ∀ N : ℕ, 1 ≤ N → ∀ T₀ X : ℝ, 0 ≤ T₀ → T₀ ≤ 1 → X₀ ≤ X →
      |∫ t in T₀..X, F0 φ B N t * rfun t| ≤ K₁ + K₂ * Real.log X := by
  obtain ⟨C₀, hC₀0, hF⟩ := F0_bound h
  obtain ⟨K₁, K₂, X₁, hX₁, hW⟩ := int_weighted_r L
  refine ⟨C₀ * |K₁|, C₀ * |K₂|, X₁, hX₁, by positivity, by positivity, fun N hN T₀ X hT₀ hT₀1 hX => ?_⟩
  have hTX : T₀ ≤ X := by linarith
  have hX0 : 0 ≤ X := by linarith
  have hlogX : 0 ≤ Real.log X := Real.log_nonneg (by linarith)
  have hi : IntervalIntegrable (fun t => F0 φ B N t * rfun t) volume T₀ X :=
    (intervalIntegrable_rfun hT₀ hX0).continuousOn_mul (continuous_F0 φ B N).continuousOn
  have hc : Continuous fun t : ℝ => 2 * t + 1 := by fun_prop
  have hk : IntervalIntegrable (fun t => C₀ * ((2 * t + 1) * |rfun t|)) volume T₀ X :=
    ((intervalIntegrable_rfun hT₀ hX0).abs.continuousOn_mul hc.continuousOn).const_mul _
  refine (intervalIntegral.abs_integral_le_integral_abs hTX).trans ?_
  have hm : ∫ t in T₀..X, |F0 φ B N t * rfun t| ≤ ∫ t in T₀..X, C₀ * ((2 * t + 1) * |rfun t|) := by
    apply intervalIntegral.integral_mono_on hTX hi.abs hk
    intro t ht
    rw [abs_mul, ← mul_assoc]
    exact mul_le_mul_of_nonneg_right (hF N hN t (le_trans hT₀ ht.1)) (abs_nonneg _)
  refine hm.trans ?_
  rw [intervalIntegral.integral_const_mul]
  have h1 := hW T₀ X hT₀ hT₀1 hX
  have h2 : K₁ + K₂ * Real.log X ≤ |K₁| + |K₂| * Real.log X := by
    have := mul_le_mul_of_nonneg_right (le_abs_self K₂) hlogX
    linarith [le_abs_self K₁]
  have h3 := mul_le_mul_of_nonneg_left (h1.trans h2) hC₀0
  linarith

end

/-- The `log` term on `[0, T₀]`. -/
lemma log_term {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) : ∃ K : ℝ, ∀ N : ℕ, 1 ≤ N → ∀ T₀ : ℝ, 0 ≤ T₀ → T₀ ≤ 1 →
    |1 / Real.pi * ∫ t in (0 : ℝ)..T₀, F0 φ B N t * Real.log (t / (2 * Real.pi))| ≤ K := by
  obtain ⟨C₀, hC₀0, hF⟩ := F0_bound h
  set J := ∫ t in (0 : ℝ)..1, |Real.log (t / (2 * Real.pi))| with hJ
  refine ⟨3 * C₀ * J, fun N hN T₀ hT₀ hT₀1 => ?_⟩
  have hπ := Real.pi_gt_three
  have hli := intervalIntegrable_log_div_two_pi (a := 0) (b := T₀) le_rfl hT₀
  have hi : IntervalIntegrable (fun t => F0 φ B N t * Real.log (t / (2 * Real.pi))) volume 0 T₀ :=
    hli.continuousOn_mul (continuous_F0 φ B N).continuousOn
  rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 1 / Real.pi)]
  have h1 : |∫ t in (0 : ℝ)..T₀, F0 φ B N t * Real.log (t / (2 * Real.pi))| ≤ 3 * C₀ * J := by
    refine (intervalIntegral.abs_integral_le_integral_abs hT₀).trans ?_
    have hm : ∫ t in (0 : ℝ)..T₀, |F0 φ B N t * Real.log (t / (2 * Real.pi))| ≤
        ∫ t in (0 : ℝ)..T₀, 3 * C₀ * |Real.log (t / (2 * Real.pi))| := by
      apply intervalIntegral.integral_mono_on hT₀ hi.abs (hli.abs.const_mul _)
      intro t ht
      rw [abs_mul]
      apply mul_le_mul_of_nonneg_right _ (abs_nonneg _)
      have h2 := hF N hN t ht.1
      have h3 : C₀ * (2 * t + 1) ≤ C₀ * 3 :=
        mul_le_mul_of_nonneg_left (by linarith [ht.2]) hC₀0
      linarith
    refine hm.trans ?_
    rw [intervalIntegral.integral_const_mul]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    exact intervalIntegral.integral_mono_interval le_rfl hT₀ hT₀1
      (Filter.Eventually.of_forall fun t => abs_nonneg _) (intervalIntegrable_log_div_two_pi le_rfl zero_le_one).abs
  have h2 : 1 / Real.pi ≤ 1 := by rw [div_le_one (by positivity)]; linarith
  calc 1 / Real.pi * |∫ t in (0 : ℝ)..T₀, F0 φ B N t * Real.log (t / (2 * Real.pi))| ≤ 1 * (3 * C₀ * J) :=
        mul_le_mul h2 h1 (abs_nonneg _) (by norm_num)
    _ = 3 * C₀ * J := one_mul _

section
variable {φ : ℝ → ℝ} {B : ℕ} (hB : 1 ≤ B) {M : ℝ} (hM0 : 0 ≤ M) (hM : ∀ y, |φ y| ≤ M)
  {N : ℕ} (hN : 1 ≤ N)
include hB hM0 hM hN

lemma card_le_two : ((B * N + 1 : ℕ) : ℝ) ≤ 2 * B * N := by
  have h1 : (1 : ℝ) ≤ B := by exact_mod_cast hB
  have h2 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  push_cast; nlinarith

lemma mem_range_le {n : ℕ} (hn : n ∈ range (B * N + 1)) : (n : ℝ) ≤ B * N := by
  have := Finset.mem_range.1 hn
  have : n ≤ B * N := by omega
  exact_mod_cast this

lemma const_term (c : ℝ) : |∑ n ∈ range (B * N + 1), φ (n / N) * c| ≤ 2 * B * N * M * |c| := by
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  have hpt : ∀ n ∈ range (B * N + 1), |φ (n / N) * c| ≤ M * |c| := fun n _ => by
    rw [abs_mul]; exact mul_le_mul_of_nonneg_right (hM _) (abs_nonneg _)
  refine (Finset.sum_le_sum hpt).trans ?_
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have := card_le_two hB hM0 hM hN
  have h0 : 0 ≤ M * |c| := by positivity
  calc ((B * N + 1 : ℕ) : ℝ) * (M * |c|) ≤ 2 * B * N * (M * |c|) := mul_le_mul_of_nonneg_right this h0
    _ = 2 * B * N * M * |c| := by ring

lemma itail_term {X : ℝ} (hX : 2 * Real.pi ≤ X) (hXN : (N : ℝ) ^ 8 ≤ X) :
    |∑ n ∈ range (B * N + 1), φ (n / N) * (-∫ t in Ioi X, hI n t)| ≤ 4 * M * (B : ℝ) ^ 3 := by
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hB1 : (1 : ℝ) ≤ B := by exact_mod_cast hB
  have hπ := Real.pi_gt_three
  have hX0 : 0 < X := by linarith
  have hsX : (N : ℝ) ^ 4 ≤ √X := by
    rw [Real.le_sqrt (by positivity) hX0.le]
    calc ((N : ℝ) ^ 4) ^ 2 = (N : ℝ) ^ 8 := by ring
      _ ≤ X := hXN
  have hpt : ∀ n ∈ range (B * N + 1), |φ (n / N) * (-∫ t in Ioi X, hI n t)| ≤
      M * (2 * ((B : ℝ) * N) ^ 2 / (Real.pi * (N : ℝ) ^ 4)) := by
    intro n hn
    have hnB := mem_range_le hB hM0 hM hN hn
    rw [abs_mul, abs_neg]
    apply mul_le_mul (hM _) _ (abs_nonneg _) hM0
    have h1 : |∫ t in Ioi X, hI n t| ≤ 2 * (n : ℝ) ^ 2 / (Real.pi * √X) := I_tail_bound n hX
    refine h1.trans ?_
    apply div_le_div₀ (by positivity) _ (by positivity) _
    · exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (Nat.cast_nonneg _) hnB 2) (by norm_num)
    · exact mul_le_mul_of_nonneg_left hsX Real.pi_pos.le
  refine (Finset.abs_sum_le_sum_abs _ _).trans ((Finset.sum_le_sum hpt).trans ?_)
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have hc := card_le_two hB hM0 hM hN
  have e : 2 * ((B : ℝ) * N) ^ 2 / (Real.pi * (N : ℝ) ^ 4) = 2 * (B : ℝ) ^ 2 / (Real.pi * (N : ℝ) ^ 2) := by
    field_simp
  rw [e]
  have h2 : 2 * (B : ℝ) ^ 2 / (Real.pi * (N : ℝ) ^ 2) ≤ 2 * (B : ℝ) ^ 2 / N := by
    apply div_le_div_of_nonneg_left (by positivity) (by positivity)
    nlinarith
  calc ((B * N + 1 : ℕ) : ℝ) * (M * (2 * (B : ℝ) ^ 2 / (Real.pi * (N : ℝ) ^ 2))) ≤
        (2 * B * N) * (M * (2 * (B : ℝ) ^ 2 / N)) :=
        mul_le_mul hc (mul_le_mul_of_nonneg_left h2 hM0) (by positivity) (by positivity)
    _ = 4 * M * (B : ℝ) ^ 3 := by field_simp; ring

lemma gfun_abs_le (n : ℕ) {t : ℝ} (ht : 0 < t) : |gfun n t| ≤ (n : ℝ) ^ 2 / t ^ 2 := by
  have hθ := liTheta_le_inv ht
  have hθ0 := liTheta_pos ht
  have hc := Real.one_sub_sq_div_two_le_cos (x := n * liTheta t)
  have hc1 := Real.cos_le_one (n * liTheta t)
  unfold gfun
  rw [abs_of_nonneg (by linarith)]
  have : (n * liTheta t) ^ 2 ≤ (n : ℝ) ^ 2 / t ^ 2 := by
    rw [mul_pow, div_eq_mul_one_div ((n : ℝ) ^ 2)]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    have := pow_le_pow_left₀ hθ0.le hθ 2
    rwa [div_pow, one_pow] at this
  linarith

lemma gX_term {X : ℝ} (hX1 : 1 ≤ X) (hXN : (N : ℝ) ^ 8 ≤ X) {D K : ℝ} (hK : 0 ≤ K) (hD : |D| ≤ K * X) :
    |∑ n ∈ range (B * N + 1), φ (n / N) * (gfun n X * D)| ≤ 2 * M * K * (B : ℝ) ^ 3 := by
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hB1 : (1 : ℝ) ≤ B := by exact_mod_cast hB
  have hX0 : 0 < X := by linarith
  have hN3 : (N : ℝ) ^ 3 ≤ X := le_trans (pow_le_pow_right₀ hN1 (by norm_num)) hXN
  have hpt : ∀ n ∈ range (B * N + 1), |φ (n / N) * (gfun n X * D)| ≤ M * (((B : ℝ) * N) ^ 2 * K / X) := by
    intro n hn
    have hnB := mem_range_le hB hM0 hM hN hn
    rw [abs_mul, abs_mul]
    apply mul_le_mul (hM _) _ (by positivity) hM0
    have h1 := gfun_abs_le hB hM0 hM hN n hX0
    have h2 : (n : ℝ) ^ 2 ≤ ((B : ℝ) * N) ^ 2 := pow_le_pow_left₀ (Nat.cast_nonneg _) hnB 2
    calc |gfun n X| * |D| ≤ ((n : ℝ) ^ 2 / X ^ 2) * (K * X) := mul_le_mul h1 hD (abs_nonneg _) (by positivity)
      _ = (n : ℝ) ^ 2 * K / X := by field_simp
      _ ≤ ((B : ℝ) * N) ^ 2 * K / X := by
          apply div_le_div_of_nonneg_right _ hX0.le
          exact mul_le_mul_of_nonneg_right h2 hK
  refine (Finset.abs_sum_le_sum_abs _ _).trans ((Finset.sum_le_sum hpt).trans ?_)
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have hc := card_le_two hB hM0 hM hN
  have h3 : ((B : ℝ) * N) ^ 2 * K / X ≤ (B : ℝ) ^ 2 * K / N := by
    rw [div_le_div_iff₀ hX0 (by positivity)]
    have : (B : ℝ) ^ 2 * K * X ≥ (B : ℝ) ^ 2 * K * (N : ℝ) ^ 3 := mul_le_mul_of_nonneg_left hN3 (by positivity)
    nlinarith
  calc ((B * N + 1 : ℕ) : ℝ) * (M * (((B : ℝ) * N) ^ 2 * K / X)) ≤ (2 * B * N) * (M * ((B : ℝ) ^ 2 * K / N)) :=
        mul_le_mul hc (mul_le_mul_of_nonneg_left h3 hM0) (by positivity) (by positivity)
    _ = 2 * M * K * (B : ℝ) ^ 3 := by field_simp

lemma ztail_term {X : ℝ} (hX1 : 1 ≤ X) (hXN : (N : ℝ) ^ 8 ≤ X) {Kt : ℝ} (hKt0 : 0 ≤ Kt)
    (hKt : ∀ n : ℕ, ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k) ≤ (n : ℝ) ^ 2 * (Kt / X ^ 2 + 10 / √X)) :
    |∑ n ∈ range (B * N + 1), φ (n / N) * ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k)| ≤
      2 * M * (B : ℝ) ^ 3 * (Kt + 10) := by
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hB1 : (1 : ℝ) ≤ B := by exact_mod_cast hB
  have hX0 : 0 < X := by linarith
  have hsX : (N : ℝ) ^ 4 ≤ √X := by
    rw [Real.le_sqrt (by positivity) hX0.le]
    calc ((N : ℝ) ^ 4) ^ 2 = (N : ℝ) ^ 8 := by ring
      _ ≤ X := hXN
  have hX2 : (N : ℝ) ^ 4 ≤ X ^ 2 := by
    have h1 : (N : ℝ) ^ 4 ≤ (N : ℝ) ^ 8 := pow_le_pow_right₀ hN1 (by norm_num)
    have h2 : X ≤ X ^ 2 := by nlinarith
    linarith
  have hfac : Kt / X ^ 2 + 10 / √X ≤ (Kt + 10) / (N : ℝ) ^ 4 := by
    have h1 : Kt / X ^ 2 ≤ Kt / (N : ℝ) ^ 4 := div_le_div_of_nonneg_left hKt0 (by positivity) hX2
    have h2 : 10 / √X ≤ 10 / (N : ℝ) ^ 4 := div_le_div_of_nonneg_left (by norm_num) (by positivity) hsX
    have e : (Kt + 10) / (N : ℝ) ^ 4 = Kt / (N : ℝ) ^ 4 + 10 / (N : ℝ) ^ 4 := by ring
    linarith
  have hpt : ∀ n ∈ range (B * N + 1), |φ (n / N) * ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k)| ≤
      M * ((B : ℝ) ^ 2 * (Kt + 10) / (N : ℝ) ^ 2) := by
    intro n hn
    have hnB := mem_range_le hB hM0 hM hN hn
    have hnn : 0 ≤ ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k) := tsum_nonneg fun k => by
      unfold gfun; have := Real.cos_le_one (n * liTheta (ordK k)); linarith
    rw [abs_mul, abs_of_nonneg hnn]
    apply mul_le_mul (hM _) _ hnn hM0
    refine (hKt n).trans ?_
    have h2 : (n : ℝ) ^ 2 ≤ ((B : ℝ) * N) ^ 2 := pow_le_pow_left₀ (Nat.cast_nonneg _) hnB 2
    calc (n : ℝ) ^ 2 * (Kt / X ^ 2 + 10 / √X) ≤ ((B : ℝ) * N) ^ 2 * ((Kt + 10) / (N : ℝ) ^ 4) :=
          mul_le_mul h2 hfac (by positivity) (by positivity)
      _ = (B : ℝ) ^ 2 * (Kt + 10) / (N : ℝ) ^ 2 := by field_simp
  refine (Finset.abs_sum_le_sum_abs _ _).trans ((Finset.sum_le_sum hpt).trans ?_)
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have hc := card_le_two hB hM0 hM hN
  calc ((B * N + 1 : ℕ) : ℝ) * (M * ((B : ℝ) ^ 2 * (Kt + 10) / (N : ℝ) ^ 2)) ≤
        (2 * B * N) * (M * ((B : ℝ) ^ 2 * (Kt + 10) / (N : ℝ) ^ 2)) :=
        mul_le_mul_of_nonneg_right hc (by positivity)
    _ = 2 * M * (B : ℝ) ^ 3 * (Kt + 10) / N := by field_simp
    _ ≤ 2 * M * (B : ℝ) ^ 3 * (Kt + 10) := div_le_self (by positivity) hN1

end

lemma trend_term {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) {M : ℝ} (hM0 : 0 ≤ M) (hM : ∀ y, |φ y| ≤ M) (N : ℕ) :
    |∑ n ∈ range (B * N + 1), φ (n / N) * (Iclosed n - liTrend n)| ≤ 4 * M * (1 + Real.log (B * N)) := by
  have hpt : ∀ n ∈ range (B * N + 1), |φ (n / N) * (Iclosed n - liTrend n)| ≤ 4 * M * (1 / (n : ℝ)) := by
    intro n _
    rcases Nat.eq_zero_or_pos n with h0 | h0
    · subst h0
      rw [h.zero N, zero_mul, abs_zero]; simp
    · rw [abs_mul]
      have h1 : |Iclosed n - liTrend n| ≤ 4 / n := Iclosed_sub_trend_le_div n h0
      calc |φ (n / N)| * |Iclosed n - liTrend n| ≤ M * (4 / n) := mul_le_mul (hM _) h1 (abs_nonneg _) hM0
        _ = 4 * M * (1 / (n : ℝ)) := by ring
  refine (Finset.abs_sum_le_sum_abs _ _).trans ((Finset.sum_le_sum hpt).trans ?_)
  rw [← Finset.mul_sum]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  rw [Finset.sum_range_succ']
  simp only [Nat.cast_zero, div_zero, add_zero]
  have e : ∑ i ∈ range (B * N), 1 / ((i + 1 : ℕ) : ℝ) = (harmonic (B * N) : ℝ) := by
    simp only [harmonic]; push_cast; simp only [one_div]
  rw [e]
  have := harmonic_le_one_add_log (B * N)
  push_cast at this ⊢
  exact this

set_option maxHeartbeats 8000000 in
/-- **Claim part 4** (Theorem 4, RH): `Σ φ(n/N) R(n) = (7/4) Σ φ(n/N) + O(log N)` for `φ ∈ C_c^∞(0, ∞)`. -/
theorem part4 (L : LiteratureInputs) : Part4 := by
  intro hRH φ hφ hc hs
  obtain ⟨B, hB1, h⟩ := testFn_of hφ hc hs
  obtain ⟨M, hM0, hM⟩ := h.bound
  obtain ⟨T₀, hT₀0, hT₀1, hz⟩ := exists_T0
  obtain ⟨Za, Zb, XZ, hXZ, hZa, hZb, hZ⟩ := zero_term L h
  obtain ⟨Wa, Wb, XW, hXW, hWa, hWb, hW⟩ := r_term L h
  obtain ⟨K7, hK7⟩ := log_term h
  obtain ⟨C₀, hC₀0, hF0⟩ := F0_bound h
  obtain ⟨Cr, Xr, hXr, hCr⟩ := rfun_tail L
  obtain ⟨C₁, T₁, hT₁, hNA⟩ := NA_bound L
  obtain ⟨Kt, hKt0, hKt⟩ := li_tail_bound L
  set X₀ := max (max XZ XW) (max Xr (max T₁ 7)) with hX₀
  have hX₀7 : 7 ≤ X₀ := le_trans (le_max_right _ _) (le_trans (le_max_right _ _) (le_max_right _ _))
  have hX₀1 : 1 ≤ X₀ := by linarith
  have hπ := Real.pi_gt_three
  have hπ4 := Real.pi_lt_d2
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hB1r : (1 : ℝ) ≤ B := by exact_mod_cast hB1
  set Q₀ := Za + Zb * Real.log X₀ + 2 * |Afun T₀| * (C₀ * (2 * T₀ + 1)) + 4 * B * M * |Cr| +
    2 * (Wa + Wb * Real.log X₀) + 4 * M * (1 + Real.log B) + 4 * M * (B : ℝ) ^ 3 + |K7| +
    2 * M * |C₁| * (B : ℝ) ^ 3 + 2 * M * (B : ℝ) ^ 3 * (Kt + 10) with hQ₀
  set Q₁ := 8 * Zb + 16 * Wb + 4 * M with hQ₁
  have hlB : 0 ≤ Real.log B := Real.log_nonneg hB1r
  have hlX₀ : 0 ≤ Real.log X₀ := Real.log_nonneg hX₀1
  have hQ₀0 : 0 ≤ Q₀ := by positivity
  have hQ₁0 : 0 ≤ Q₁ := by positivity
  refine ⟨Q₀ / Real.log 2 + Q₁, fun N hN => ?_⟩
  have hN1 : 1 ≤ N := by omega
  have hNr : (2 : ℝ) ≤ N := by exact_mod_cast hN
  have hN1r : (1 : ℝ) ≤ N := by linarith
  have hlN : Real.log 2 ≤ Real.log N := Real.log_le_log (by norm_num) hNr
  set X := X₀ * (N : ℝ) ^ 8 with hXdef
  have hN8 : (1 : ℝ) ≤ (N : ℝ) ^ 8 := one_le_pow₀ hN1r
  have hXN : (N : ℝ) ^ 8 ≤ X := le_mul_of_one_le_left (by positivity) hX₀1
  have hX₀X : X₀ ≤ X := le_mul_of_one_le_right (by linarith) hN8
  have hX1 : 1 ≤ X := by linarith
  have hX0 : 0 < X := by linarith
  have hT₀X : T₀ ≤ X := by linarith
  have hNX : (N : ℝ) ≤ X := le_trans (le_self_pow₀ hN1r (by norm_num)) hXN
  have hlogX : Real.log X = Real.log X₀ + 8 * Real.log N := by
    rw [hXdef, Real.log_mul (by positivity) (by positivity), Real.log_pow]; push_cast; ring
  have hXZX : XZ ≤ X := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hX₀X
  have hXWX : XW ≤ X := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hX₀X
  have hXrX : Xr ≤ X := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hX₀X
  have hT₁X : T₁ ≤ X := le_trans (le_trans (le_max_left _ _) (le_trans (le_max_right _ _) (le_max_right _ _))) hX₀X
  -- finite sums
  have hts1 : ∑' n : ℕ, φ (n / N) * R n = ∑ n ∈ range (B * N + 1), φ (n / N) * R n := by
    apply tsum_eq_sum; intro n hn
    rw [Finset.mem_range, not_lt] at hn
    rw [h.vanish hN1 (by omega), zero_mul]
  have hts2 : ∑' n : ℕ, φ (n / N) = ∑ n ∈ range (B * N + 1), φ (n / N) := by
    apply tsum_eq_sum; intro n hn
    rw [Finset.mem_range, not_lt] at hn
    exact h.vanish hN1 (by omega)
  set S := range (B * N + 1) with hS
  have key : ∑ n ∈ S, φ (n / N) * R n - 7 / 4 * ∑ n ∈ S, φ (n / N) =
      (∑ n ∈ S, φ (n / N) * (2 * n * aCoef T₀ X n)) +
      (∑ n ∈ S, φ (n / N) * (-2 * Afun T₀ * Real.cos (n * liTheta T₀))) +
      (∑ n ∈ S, φ (n / N) * (2 * ((∫ t in (0 : ℝ)..X, rfun t) + 1 / 8))) +
      (∑ n ∈ S, φ (n / N) * (-2 * ∫ t in T₀..X, Real.cos (n * liTheta t) * rfun t)) +
      (∑ n ∈ S, φ (n / N) * (Iclosed n - liTrend n)) +
      (∑ n ∈ S, φ (n / N) * (-∫ t in Ioi X, hI n t)) +
      (∑ n ∈ S, φ (n / N) * (1 / Real.pi *
        (∫ t in (0 : ℝ)..T₀, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))))) +
      (∑ n ∈ S, φ (n / N) * (gfun n X * ((zeroCount X : ℝ) - Afun X))) +
      (∑ n ∈ S, φ (n / N) * ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k)) := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
    simp only [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun n _ => ?_
    have := R_sub_eq L hRH hT₀0 hT₀X hz n
    linear_combination φ (n / N) * this
  rw [hts1, hts2, key]
  -- the nine bounds
  have b1 : |∑ n ∈ S, φ (n / N) * (2 * n * aCoef T₀ X n)| ≤ Za + Zb * Real.log X := by
    rw [sum_mul_aCoef]; exact hZ N hN1 T₀ X hT₀0.le hT₀X hXZX
  have b2 : |∑ n ∈ S, φ (n / N) * (-2 * Afun T₀ * Real.cos (n * liTheta T₀))| ≤
      2 * |Afun T₀| * (C₀ * (2 * T₀ + 1)) := by
    have e : ∑ n ∈ S, φ (n / N) * (-2 * Afun T₀ * Real.cos (n * liTheta T₀)) = -2 * Afun T₀ * F0 φ B N T₀ := by
      unfold F0; rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun n _ => by ring
    rw [e, abs_mul, abs_mul, abs_neg, abs_two]
    exact mul_le_mul_of_nonneg_left (hF0 N hN1 T₀ hT₀0.le) (by positivity)
  have b3 : |∑ n ∈ S, φ (n / N) * (2 * ((∫ t in (0 : ℝ)..X, rfun t) + 1 / 8))| ≤ 4 * B * M * |Cr| := by
    refine (const_term hB1 hM0 hM hN1 _).trans ?_
    have h1 := hCr X hXrX
    rw [abs_mul, abs_two]
    have h2 : |(∫ t in (0 : ℝ)..X, rfun t) + 1 / 8| ≤ |Cr| / X :=
      h1.trans (div_le_div_of_nonneg_right (le_abs_self _) hX0.le)
    calc 2 * B * N * M * (2 * |(∫ t in (0 : ℝ)..X, rfun t) + 1 / 8|) ≤ 2 * B * N * M * (2 * (|Cr| / X)) :=
          mul_le_mul_of_nonneg_left (by linarith) (by positivity)
      _ = 4 * B * M * |Cr| * (N / X) := by ring
      _ ≤ 4 * B * M * |Cr| * 1 :=
          mul_le_mul_of_nonneg_left ((div_le_one hX0).2 hNX) (by positivity)
      _ = 4 * B * M * |Cr| := mul_one _
  have b4 : |∑ n ∈ S, φ (n / N) * (-2 * ∫ t in T₀..X, Real.cos (n * liTheta t) * rfun t)| ≤
      2 * (Wa + Wb * Real.log X) := by
    have e : ∑ n ∈ S, φ (n / N) * (-2 * ∫ t in T₀..X, Real.cos (n * liTheta t) * rfun t) =
        -2 * ∫ t in T₀..X, F0 φ B N t * rfun t := by
      rw [← sum_mul_integral_cos_r φ B N hT₀0.le hX0.le, Finset.mul_sum]
      exact Finset.sum_congr rfl fun n _ => by ring
    rw [e, abs_mul, abs_neg, abs_two]
    exact mul_le_mul_of_nonneg_left (hW N hN1 T₀ X hT₀0.le hT₀1 hXWX) (by norm_num)
  have b5 := trend_term h hM0 hM N
  have b6 := itail_term hB1 hM0 hM hN1 (X := X) (by linarith) hXN
  have b7 : |∑ n ∈ S, φ (n / N) * (1 / Real.pi *
      (∫ t in (0 : ℝ)..T₀, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))))| ≤ |K7| := by
    rw [sum_mul_integral_cos_log φ B N le_rfl hT₀0.le]
    exact (hK7 N hN1 T₀ hT₀0.le hT₀1).trans (le_abs_self _)
  have hD : |(zeroCount X : ℝ) - Afun X| ≤ |C₁| * X := by
    have h1 := hNA X hT₁X
    have h2 : Real.log X ≤ X := (Real.log_le_sub_one_of_pos hX0).trans (by linarith)
    have h3 : 0 ≤ Real.log X := Real.log_nonneg hX1
    calc |(zeroCount X : ℝ) - Afun X| ≤ C₁ * Real.log X := h1
      _ ≤ |C₁| * Real.log X := mul_le_mul_of_nonneg_right (le_abs_self _) h3
      _ ≤ |C₁| * X := mul_le_mul_of_nonneg_left h2 (abs_nonneg _)
  have b8 := gX_term hB1 hM0 hM hN1 hX1 hXN (abs_nonneg C₁) hD
  have b9 := ztail_term hB1 hM0 hM hN1 hX1 hXN hKt0 (fun n => hKt n X hX0)
  -- combine
  have hlBN : Real.log ((B : ℝ) * N) = Real.log B + Real.log N := Real.log_mul (by positivity) (by positivity)
  rw [hlBN] at b5
  rw [hlogX] at b1 b4
  have htot : Q₀ + Q₁ * Real.log N ≤ (Q₀ / Real.log 2 + Q₁) * Real.log N := by
    have : Q₀ ≤ Q₀ / Real.log 2 * Real.log N := by
      rw [div_mul_eq_mul_div, le_div_iff₀ hl2]
      exact mul_le_mul_of_nonneg_left hlN hQ₀0
    linarith
  refine le_trans ?_ htot
  rw [abs_le] at b1 b2 b3 b4 b5 b6 b7 b8 b9 ⊢
  constructor <;> linarith [b1.1, b1.2, b2.1, b2.2, b3.1, b3.2, b4.1, b4.2, b5.1, b5.2, b6.1, b6.2,
    b7.1, b7.2, b8.1, b8.2, b9.1, b9.2]

end LeanLi

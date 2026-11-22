module

public import Riemann.Deep.LeanLi.LiTheta
public import Riemann.Deep.LeanLi.ZeroSums
public import Riemann.Deep.LeanLi.LR.Trend

/-!
# LeanLi: (H7) `τ_T = -Σ_{γ≤T} log(2πθ_γ)/θ_γ² ≤ (1 + η) T³L²/(6π)` (proofs.md step 7)
-/

@[expose] public section

open Real Finset Filter Set

namespace LeanLi

/-- A sum over `Z(T)` (with multiplicity) as a sum over the zeros of `zset T`. -/
lemma sum_ZK (T : ℝ) (g : ℝ → ℝ) :
    ∑ j ∈ ZK T, g (ordK j) = ∑ ρ ∈ zset T, (mult ρ : ℝ) * g ρ.im := by
  rw [← upperUpTo_map, Finset.sum_map]
  unfold ZK
  rw [Finset.sum_sigma]
  refine Finset.sum_congr rfl fun ρ _ => ?_
  simp [ordK, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]

/-- `F(t) = t² log(t/2π)`. -/
noncomputable def Fq (t : ℝ) : ℝ := t ^ 2 * log (t / (2 * π))

lemma hasDerivAt_Fq {t : ℝ} (ht : 0 < t) :
    HasDerivAt Fq (2 * t * log (t / (2 * π)) + t) t := by
  have h1 : HasDerivAt (fun t : ℝ => t / (2 * π)) (1 / (2 * π)) t := by
    simpa using (hasDerivAt_id t).div_const (2 * π)
  have h2 := h1.log (by positivity)
  have h3 := (hasDerivAt_pow 2 t).mul h2
  unfold Fq
  convert h3 using 1
  field_simp
  ring

lemma hasDerivAt_mainN {t : ℝ} (ht : 0 < t) :
    HasDerivAt mainN (1 / (2 * π) * log (t / (2 * π))) t := by
  have h1 : HasDerivAt (fun t : ℝ => t / (2 * π * exp 1)) (1 / (2 * π * exp 1)) t := by
    simpa using (hasDerivAt_id t).div_const (2 * π * exp 1)
  have h2 := h1.log (by positivity)
  have h3 : HasDerivAt (fun y : ℝ => y / (2 * π) * log (y / (2 * π * exp 1))) _ t :=
    ((hasDerivAt_id' t).div_const (2 * π)).mul h2
  unfold mainN
  convert h3 using 1
  have e : log (t / (2 * π * exp 1)) = log (t / (2 * π)) - 1 := by
    rw [show t / (2 * π * exp 1) = t / (2 * π) / exp 1 by ring,
      log_div (by positivity) (exp_pos 1).ne', log_exp]
  rw [e]
  field_simp
  ring

/-- **Abel summation for `F`**: `Σ_{T₀<γ≤T} m F(γ) ≤ T³L²/(6π) + F(T₀) mainN(T₀) + 2|C| T² log(T/2π) log T`. -/
lemma sum_Fq_le {C T₀ : ℝ} (hT₀ : 2 * π ≤ T₀)
    (hC : ∀ t, T₀ ≤ t → |(zeroCount t : ℝ) - mainN t| ≤ C * log t) {T : ℝ} (hT : T₀ ≤ T) :
    ∑ ρ ∈ (zset T).filter (fun ρ => T₀ < ρ.im), (mult ρ : ℝ) * Fq ρ.im ≤
      T ^ 3 * log (T / (2 * π)) ^ 2 / (6 * π) + Fq T₀ * mainN T₀ +
        2 * |C| * (T ^ 2 * log (T / (2 * π)) * log T) := by
  have hπ := pi_pos
  have hπ3 := pi_gt_three
  have hT₀pos : 0 < T₀ := by linarith
  have hTpos : 0 < T := by linarith
  have hmem : ∀ t ∈ Set.uIcc T₀ T, T₀ ≤ t ∧ t ≤ T := fun t ht => by
    rw [Set.uIcc_of_le hT] at ht; exact ht
  set F' : ℝ → ℝ := fun t => 2 * t * log (t / (2 * π)) + t with hF'
  have hFd : ∀ t ∈ Icc T₀ T, HasDerivAt Fq (F' t) t := fun t ht =>
    hasDerivAt_Fq (by linarith [ht.1])
  have hFd' : ∀ t ∈ Set.uIcc T₀ T, HasDerivAt Fq (F' t) t := fun t ht =>
    hasDerivAt_Fq (by linarith [(hmem t ht).1])
  have hcF' : ContinuousOn F' (Set.uIcc T₀ T) := by
    intro t ht
    have h0 : 0 < t := by linarith [(hmem t ht).1]
    have hne : t / (2 * π) ≠ 0 := (div_pos h0 (by positivity)).ne'
    apply ContinuousAt.continuousWithinAt
    simp only [hF']
    fun_prop (disch := assumption)
  have hiF' : IntervalIntegrable F' MeasureTheory.volume T₀ T := hcF'.intervalIntegrable
  have hab := zsum_abel hT hFd hiF'
  rw [hab]
  -- notation
  have hF'nn : ∀ t ∈ Icc T₀ T, 0 ≤ F' t := by
    intro t ht
    have ht0 : 0 ≤ t := by linarith [ht.1]
    have hl : 0 ≤ log (t / (2 * π)) := log_nonneg (by rw [le_div_iff₀ (by positivity)]; linarith [ht.1])
    simp only [hF']
    have := mul_nonneg (mul_nonneg (by norm_num : (0:ℝ) ≤ 2) ht0) hl
    linarith
  have hCabs : ∀ t, T₀ ≤ t → |(zeroCount t : ℝ) - mainN t| ≤ |C| * log t := fun t ht =>
    (hC t ht).trans (mul_le_mul_of_nonneg_right (le_abs_self C) (log_nonneg (by linarith)))
  -- lower bound for the integral
  have hNi : IntervalIntegrable (fun t => (zeroCount t : ℝ)) MeasureTheory.volume T₀ T :=
    zeroCount_mono.intervalIntegrable
  have hcm : ContinuousOn mainN (Set.uIcc T₀ T) := fun t ht =>
    (hasDerivAt_mainN (by linarith [(hmem t ht).1])).continuousAt.continuousWithinAt
  have hcl : ContinuousOn log (Set.uIcc T₀ T) := fun t ht =>
    (continuousAt_log (by linarith [(hmem t ht).1])).continuousWithinAt
  have hint1 : IntervalIntegrable (fun t => F' t * ((zeroCount t : ℝ) - zeroCount T₀)) MeasureTheory.volume T₀ T :=
    (hNi.sub intervalIntegrable_const).continuousOn_mul hcF'
  have hint2 : IntervalIntegrable (fun t => F' t * (mainN t - |C| * log t - zeroCount T₀)) MeasureTheory.volume T₀ T :=
    ((hcm.sub (hcl.const_smul |C|)).sub continuousOn_const).intervalIntegrable.continuousOn_mul hcF'
  have hlow : ∫ t in T₀..T, F' t * (mainN t - |C| * log t - zeroCount T₀) ≤
      ∫ t in T₀..T, F' t * ((zeroCount t : ℝ) - zeroCount T₀) := by
    apply intervalIntegral.integral_mono_on hT hint2 hint1
    intro t ht
    apply mul_le_mul_of_nonneg_left _ (hF'nn t ht)
    have := (abs_le.1 (hCabs t ht.1)).1
    linarith
  -- evaluate the lower integral
  have hiFm : IntervalIntegrable (fun t => F' t * mainN t) MeasureTheory.volume T₀ T :=
    (hcF'.mul hcm).intervalIntegrable
  have hiFl : IntervalIntegrable (fun t => F' t * log t) MeasureTheory.volume T₀ T :=
    (hcF'.mul hcl).intervalIntegrable
  have hsplit : ∫ t in T₀..T, F' t * (mainN t - |C| * log t - zeroCount T₀) =
      (∫ t in T₀..T, F' t * mainN t) - |C| * (∫ t in T₀..T, F' t * log t) -
        (zeroCount T₀ : ℝ) * ∫ t in T₀..T, F' t := by
    rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_const_mul,
      ← intervalIntegral.integral_sub hiFm (hiFl.const_mul _),
      ← intervalIntegral.integral_sub (hiFm.sub (hiFl.const_mul _)) (hiF'.const_mul _)]
    refine intervalIntegral.integral_congr fun t _ => ?_
    ring
  have hftc : ∫ t in T₀..T, F' t = Fq T - Fq T₀ :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt hFd' hiF'
  have hFl : ∫ t in T₀..T, F' t * log t ≤ log T * (Fq T - Fq T₀) := by
    rw [← hftc, ← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_mono_on hT hiFl (hiF'.const_mul _)
    intro t ht
    rw [mul_comm (log T)]
    exact mul_le_mul_of_nonneg_left (log_le_log (by linarith [ht.1]) ht.2) (hF'nn t ht)
  -- by parts: ∫ F' mainN = F(T)mainN(T) - F(T₀)mainN(T₀) - ∫ F mainN'
  have hmd : ∀ t ∈ Set.uIcc T₀ T, HasDerivAt mainN (1 / (2 * π) * log (t / (2 * π))) t := fun t ht =>
    hasDerivAt_mainN (by linarith [(hmem t ht).1])
  have hcmd : ContinuousOn (fun t => 1 / (2 * π) * log (t / (2 * π))) (Set.uIcc T₀ T) := by
    intro t ht
    have h0 : 0 < t := by linarith [(hmem t ht).1]
    have hne : t / (2 * π) ≠ 0 := (div_pos h0 (by positivity)).ne'
    apply ContinuousAt.continuousWithinAt
    fun_prop (disch := assumption)
  have hparts := intervalIntegral.integral_mul_deriv_eq_deriv_mul hmd hFd' hcmd.intervalIntegrable hiF'
  -- hparts : ∫ mainN * F' = mainN T * Fq T - mainN T₀ * Fq T₀ - ∫ (1/2π log) * Fq
  have hFm : ∫ t in T₀..T, F' t * mainN t =
      mainN T * Fq T - mainN T₀ * Fq T₀ - ∫ t in T₀..T, 1 / (2 * π) * log (t / (2 * π)) * Fq t := by
    rw [← hparts]
    refine intervalIntegral.integral_congr fun t _ => ?_
    ring
  -- ∫ F mainN' ≤ T³L²/(6π)
  have hcFq : ContinuousOn Fq (Set.uIcc T₀ T) := fun t ht => (hFd' t ht).continuousAt.continuousWithinAt
  have hmain : ∫ t in T₀..T, 1 / (2 * π) * log (t / (2 * π)) * Fq t ≤
      T ^ 3 * log (T / (2 * π)) ^ 2 / (6 * π) := by
    have hLT : 0 ≤ log (T / (2 * π)) := log_nonneg (by rw [le_div_iff₀ (by positivity)]; linarith)
    have h1 : ∫ t in T₀..T, 1 / (2 * π) * log (t / (2 * π)) * Fq t ≤
        ∫ t in T₀..T, 1 / (2 * π) * log (T / (2 * π)) ^ 2 * t ^ 2 := by
      apply intervalIntegral.integral_mono_on hT (hcmd.mul hcFq).intervalIntegrable
        (Continuous.intervalIntegrable (by fun_prop) _ _)
      intro t ht
      show 1 / (2 * π) * log (t / (2 * π)) * Fq t ≤ 1 / (2 * π) * log (T / (2 * π)) ^ 2 * t ^ 2
      have ht0 : 0 < t := by linarith [ht.1]
      have hl0 : 0 ≤ log (t / (2 * π)) := log_nonneg (by
        rw [le_div_iff₀ (by positivity)]; linarith [ht.1])
      have hl1 : log (t / (2 * π)) ≤ log (T / (2 * π)) :=
        log_le_log (by positivity) (div_le_div_of_nonneg_right ht.2 (by positivity))
      unfold Fq
      have : log (t / (2 * π)) * (t ^ 2 * log (t / (2 * π))) ≤ log (T / (2 * π)) ^ 2 * t ^ 2 := by
        have := mul_le_mul hl1 hl1 hl0 hLT
        nlinarith [sq_nonneg t]
      have h2π : 0 ≤ 1 / (2 * π) := by positivity
      nlinarith [mul_le_mul_of_nonneg_left this h2π]
    have h2 : ∫ t in T₀..T, 1 / (2 * π) * log (T / (2 * π)) ^ 2 * t ^ 2 =
        1 / (2 * π) * log (T / (2 * π)) ^ 2 * ((T ^ 3 - T₀ ^ 3) / 3) := by
      rw [intervalIntegral.integral_const_mul, integral_pow]; ring
    rw [h2] at h1
    have h3 : 1 / (2 * π) * log (T / (2 * π)) ^ 2 * ((T ^ 3 - T₀ ^ 3) / 3) ≤
        T ^ 3 * log (T / (2 * π)) ^ 2 / (6 * π) := by
      have : 0 ≤ T₀ ^ 3 := by positivity
      rw [show T ^ 3 * log (T / (2 * π)) ^ 2 / (6 * π) =
        1 / (2 * π) * log (T / (2 * π)) ^ 2 * (T ^ 3 / 3) by field_simp; ring]
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      linarith
    linarith
  -- combine
  have hFT₀ : 0 ≤ Fq T₀ := by
    unfold Fq; have := log_nonneg (show 1 ≤ T₀ / (2 * π) by rw [le_div_iff₀ (by positivity)]; linarith)
    positivity
  have hNT := hCabs T hT
  have hN0 : (0 : ℝ) ≤ zeroCount T₀ := Nat.cast_nonneg _
  have hFT : 0 ≤ Fq T := by
    unfold Fq; have := log_nonneg (show 1 ≤ T / (2 * π) by rw [le_div_iff₀ (by positivity)]; linarith)
    positivity
  have hlogT : 0 ≤ log T := log_nonneg (by linarith)
  have hNTu : (zeroCount T : ℝ) ≤ mainN T + |C| * log T := by linarith [(abs_le.1 hNT).2]
  have hFTeq : Fq T * log T = T ^ 2 * log (T / (2 * π)) * log T := by unfold Fq; ring
  have hk1 : Fq T * ((zeroCount T : ℝ) - zeroCount T₀) ≤ Fq T * (mainN T + |C| * log T - zeroCount T₀) :=
    mul_le_mul_of_nonneg_left (by linarith) hFT
  have hk2 : |C| * (log T * (Fq T - Fq T₀)) ≤ |C| * (Fq T * log T) := by
    apply mul_le_mul_of_nonneg_left _ (abs_nonneg _)
    nlinarith
  have hk3 : |C| * (∫ t in T₀..T, F' t * log t) ≤ |C| * (log T * (Fq T - Fq T₀)) :=
    mul_le_mul_of_nonneg_left hFl (abs_nonneg _)
  rw [hsplit, hFm, hftc] at hlow
  have hFTeq' : |C| * (Fq T * log T) = |C| * (T ^ 2 * log (T / (2 * π)) * log T) := by rw [hFTeq]
  have hnn := mul_nonneg hN0 hFT₀
  linarith [hk1, hk2, hk3, hlow, hmain, hnn, hFTeq']

lemma zset_filter_le {T₀ T : ℝ} (h : T₀ ≤ T) :
    (zset T).filter (fun ρ => ρ.im ≤ T₀) = zset T₀ := by
  ext ρ
  simp only [Finset.mem_filter, mem_zset, zerosUpTo, Set.mem_setOf_eq]
  constructor
  · rintro ⟨⟨h1, h2, _⟩, h4⟩; exact ⟨h1, h2, h4⟩
  · rintro ⟨h1, h2, h3⟩; exact ⟨⟨h1, h2, h3.trans h⟩, h3⟩

set_option maxHeartbeats 2000000 in
/-- **(H7)** `τ_T ≤ (1 + η) T³L²/(6π)` eventually. -/
theorem tauT (L : LiteratureInputs) : LR.TauT thK ZK := by
  intro η hη
  have hπ := pi_gt_three
  have hπ0 := pi_pos
  obtain ⟨C, T₀c, hT₀c4, hC⟩ := count_asymp L
  obtain ⟨T₁c, hT₁c, hcount⟩ := count_le L
  set T₀ := max (2 * π) T₀c with hT₀
  have hT₀2π : 2 * π ≤ T₀ := le_max_left _ _
  have hC' : ∀ t, T₀ ≤ t → |(zeroCount t : ℝ) - mainN t| ≤ C * log t :=
    fun t ht => hC t (le_trans (le_max_right _ _) ht)
  set f : ℝ → ℝ := fun γ => -log (2 * π * liTheta γ) / liTheta γ ^ 2 with hf
  set Ka : ℝ := ∑ ρ ∈ zset T₀, (mult ρ : ℝ) * f ρ.im
  set Kb : ℝ := Fq T₀ * mainN T₀
  set T₂ : ℝ := max (max (max T₀ T₁c) (max ((2 * π) ^ 2) (exp 1)))
    (max (48 * π * (2 * |C| + 2) / η) (48 * π * (|Ka| + |Kb|) / η + 1)) with hT₂
  filter_upwards [eventually_ge_atTop T₂] with T hT
  have hTT₀ : T₀ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) (le_trans (le_max_left _ _) hT)
  have hTT₁ : T₁c ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) (le_trans (le_max_left _ _) hT)
  have hT2π2 : (2 * π) ^ 2 ≤ T :=
    le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) (le_trans (le_max_left _ _) hT)
  have hTe : exp 1 ≤ T :=
    le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) (le_trans (le_max_left _ _) hT)
  have hTC : 48 * π * (2 * |C| + 2) / η ≤ T :=
    le_trans (le_max_left _ _) (le_trans (le_max_right _ _) hT)
  have hTK : 48 * π * (|Ka| + |Kb|) / η + 1 ≤ T :=
    le_trans (le_max_right _ _) (le_trans (le_max_right _ _) hT)
  have hTpos : 0 < T := lt_of_lt_of_le (exp_pos 1) hTe
  have hlogT1 : 1 ≤ log T := by rw [← log_exp 1]; exact log_le_log (exp_pos 1) hTe
  have hT1 : 1 ≤ T := by linarith [add_one_le_exp (1 : ℝ)]
  -- rewrite the sum
  have hsum : -∑ j ∈ ZK T, log (2 * π * thK j) / thK j ^ 2 = ∑ ρ ∈ zset T, (mult ρ : ℝ) * f ρ.im := by
    rw [← Finset.sum_neg_distrib]
    rw [← sum_ZK T f]
    refine Finset.sum_congr rfl fun j _ => ?_
    simp only [hf, thK, neg_div]
  rw [hsum, ← Finset.sum_filter_add_sum_filter_not (zset T) (fun ρ => ρ.im ≤ T₀),
    zset_filter_le hTT₀]
  -- the tail
  have hfilt : ∀ ρ ∈ (zset T).filter (fun ρ => ¬ ρ.im ≤ T₀), T₀ < ρ.im ∧ ρ.im ≤ T := by
    intro ρ hρ
    rw [Finset.mem_filter, mem_zset] at hρ
    exact ⟨not_le.1 hρ.2, hρ.1.2.2⟩
  have hfilt' : (zset T).filter (fun ρ => ¬ ρ.im ≤ T₀) = (zset T).filter (fun ρ => T₀ < ρ.im) := by
    ext ρ; simp only [Finset.mem_filter, not_le]
  have hpt : ∑ ρ ∈ (zset T).filter (fun ρ => ¬ ρ.im ≤ T₀), (mult ρ : ℝ) * f ρ.im ≤
      ∑ ρ ∈ (zset T).filter (fun ρ => T₀ < ρ.im), (mult ρ : ℝ) * Fq ρ.im +
        (log T + 1) * ∑ ρ ∈ zset T, (mult ρ : ℝ) := by
    rw [hfilt', Finset.mul_sum]
    have h1 : ∑ ρ ∈ (zset T).filter (fun ρ => T₀ < ρ.im), (mult ρ : ℝ) * f ρ.im ≤
        ∑ ρ ∈ (zset T).filter (fun ρ => T₀ < ρ.im), ((mult ρ : ℝ) * Fq ρ.im + (log T + 1) * (mult ρ : ℝ)) := by
      apply Finset.sum_le_sum
      intro ρ hρ
      have hρ' := hfilt ρ (by rw [hfilt']; exact hρ)
      have hγ : 2 * π ≤ ρ.im := by linarith [hρ'.1]
      have hb := neg_log_div_sq_le hγ
      have hlg : log ρ.im ≤ log T := log_le_log (by linarith) hρ'.2
      have hm : (0 : ℝ) ≤ mult ρ := Nat.cast_nonneg _
      have : f ρ.im ≤ Fq ρ.im + (log T + 1) := by
        simp only [hf]; unfold Fq; linarith
      nlinarith
    rw [Finset.sum_add_distrib] at h1
    have h2 : ∑ ρ ∈ (zset T).filter (fun ρ => T₀ < ρ.im), (log T + 1) * (mult ρ : ℝ) ≤
        ∑ ρ ∈ zset T, (log T + 1) * (mult ρ : ℝ) :=
      Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) fun ρ _ _ => by
        have : (0 : ℝ) ≤ mult ρ := Nat.cast_nonneg _
        positivity
    linarith
  have hF := sum_Fq_le hT₀2π hC' hTT₀
  have hN : ∑ ρ ∈ zset T, (mult ρ : ℝ) ≤ T * log T := by
    rw [← zeroCount_eq]; exact hcount T hTT₁
  -- error terms
  set Lg := log (T / (2 * π)) with hLg
  have hL : log T / 2 ≤ Lg := by
    rw [hLg, log_div hTpos.ne' (by positivity)]
    have h := log_le_log (by positivity) hT2π2
    rw [log_pow] at h; push_cast at h
    linarith
  have hL0 : 0 ≤ Lg := by linarith
  have hLle : Lg ≤ log T := by
    rw [hLg]; exact log_le_log (by positivity) (div_le_self hTpos.le (by linarith))
  have e1 : (log T + 1) * ∑ ρ ∈ zset T, (mult ρ : ℝ) ≤ 2 * (T ^ 2 * log T ^ 2) := by
    calc (log T + 1) * ∑ ρ ∈ zset T, (mult ρ : ℝ) ≤ (log T + 1) * (T * log T) :=
          mul_le_mul_of_nonneg_left hN (by linarith)
      _ ≤ (2 * log T) * (T * log T) := by
          apply mul_le_mul_of_nonneg_right (by linarith) (by positivity)
      _ = 2 * (T * log T ^ 2) := by ring
      _ ≤ 2 * (T ^ 2 * log T ^ 2) := by
          have : T ≤ T ^ 2 := by nlinarith
          nlinarith [sq_nonneg (log T)]
  have e2 : 2 * |C| * (T ^ 2 * Lg * log T) ≤ 2 * |C| * (T ^ 2 * log T ^ 2) := by
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    have := mul_le_mul_of_nonneg_left hLle (by positivity : (0 : ℝ) ≤ T ^ 2 * log T)
    nlinarith
  have hkey : |Ka| + |Kb| + (2 * |C| + 2) * (T ^ 2 * log T ^ 2) ≤ η * (T ^ 3 * Lg ^ 2 / (6 * π)) := by
    have hlow : η * (T ^ 3 * log T ^ 2) / (24 * π) ≤ η * (T ^ 3 * Lg ^ 2 / (6 * π)) := by
      have : log T ^ 2 / 4 ≤ Lg ^ 2 := by
        have := pow_le_pow_left₀ (by positivity) hL 2
        nlinarith
      have h' : η * (T ^ 3 * log T ^ 2) / (24 * π) = η * (T ^ 3 * (log T ^ 2 / 4) / (6 * π)) := by ring
      rw [h']
      apply mul_le_mul_of_nonneg_left _ hη.le
      apply div_le_div_of_nonneg_right _ (by positivity)
      exact mul_le_mul_of_nonneg_left this (by positivity)
    have hA : (2 * |C| + 2) * (T ^ 2 * log T ^ 2) ≤ η * (T ^ 3 * log T ^ 2) / (48 * π) := by
      have h1 : 48 * π * (2 * |C| + 2) ≤ η * T := by
        rw [div_le_iff₀ hη] at hTC; linarith
      rw [le_div_iff₀ (by positivity)]
      have h2 : 0 ≤ T ^ 2 * log T ^ 2 := by positivity
      have h3 := mul_le_mul_of_nonneg_right h1 h2
      have e1 : (2 * |C| + 2) * (T ^ 2 * log T ^ 2) * (48 * π) =
          48 * π * (2 * |C| + 2) * (T ^ 2 * log T ^ 2) := by ring
      have e2 : η * T * (T ^ 2 * log T ^ 2) = η * (T ^ 3 * log T ^ 2) := by ring
      rw [e1]; linarith
    have hB : |Ka| + |Kb| ≤ η * (T ^ 3 * log T ^ 2) / (48 * π) := by
      have h1 : 48 * π * (|Ka| + |Kb|) ≤ η * T := by
        have : 48 * π * (|Ka| + |Kb|) / η ≤ T := by linarith
        rw [div_le_iff₀ hη] at this; linarith
      have h3 : T ≤ T ^ 3 * log T ^ 2 := by
        have h4 : 1 ≤ T ^ 2 * log T ^ 2 :=
          one_le_mul_of_one_le_of_one_le (one_le_pow₀ hT1) (one_le_pow₀ hlogT1)
        have := le_mul_of_one_le_right hTpos.le h4
        have e : T * (T ^ 2 * log T ^ 2) = T ^ 3 * log T ^ 2 := by ring
        linarith
      rw [le_div_iff₀ (by positivity)]
      have h5 := mul_le_mul_of_nonneg_left h3 hη.le
      have e : (|Ka| + |Kb|) * (48 * π) = 48 * π * (|Ka| + |Kb|) := by ring
      rw [e]; linarith
    have : η * (T ^ 3 * log T ^ 2) / (48 * π) + η * (T ^ 3 * log T ^ 2) / (48 * π) =
        η * (T ^ 3 * log T ^ 2) / (24 * π) := by ring
    linarith
  have hKa : Ka ≤ |Ka| := le_abs_self Ka
  have hKb : Fq T₀ * mainN T₀ ≤ |Kb| := le_abs_self Kb
  have e3 : (1 + η) * (T ^ 3 * Lg ^ 2 / (6 * π)) = T ^ 3 * Lg ^ 2 / (6 * π) + η * (T ^ 3 * Lg ^ 2 / (6 * π)) := by
    ring
  have e4 : (2 * |C| + 2) * (T ^ 2 * log T ^ 2) = 2 * |C| * (T ^ 2 * log T ^ 2) + 2 * (T ^ 2 * log T ^ 2) := by
    ring
  rw [e3]
  linarith [hpt, hF, e1, e2, hkey]

end LeanLi

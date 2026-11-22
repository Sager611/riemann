module

public import Riemann.Deep.LeanLi.SelbergT
public import Riemann.Deep.LeanLi.Bessel
public import Riemann.Deep.LeanLi.Tails
public import Riemann.Deep.LeanLi.TrendPt2

/-!
# LeanLi: the auxiliary bounds of Lemma A (proofs.md) used by Theorem 2
-/

@[expose] public section

open Real Filter Set MeasureTheory

namespace LeanLi

lemma intervalIntegrable_log_div_two_pi {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    IntervalIntegrable (fun t => Real.log (t / (2 * Real.pi))) volume a b := by
  have h := (intervalIntegral.intervalIntegrable_log' (a := a) (b := b)).sub
    (intervalIntegrable_const (c := Real.log (2 * Real.pi)))
  refine h.congr_uIoo fun t ht => ?_
  have ht0 : 0 < t := by
    rcases ht with ⟨h1, _⟩
    have : 0 ≤ min a b := le_min ha hb
    linarith
  simp only
  rw [Real.log_div ht0.ne' (by positivity)]

lemma intervalIntegrable_rfun {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    IntervalIntegrable rfun volume a b :=
  ((continuous_Ad.intervalIntegrable a b).sub
    ((intervalIntegrable_log_div_two_pi ha hb).const_mul _)).congr_uIoo fun t _ => rfl

lemma integral_inv_sq_le' {a b C : ℝ} (ha : 0 < a) (hab : a ≤ b) (hC : 0 ≤ C) :
    ∫ t in a..b, C / t ^ 2 ≤ C / a := by
  have hd : ∀ t ∈ Set.uIcc a b, HasDerivAt (fun t : ℝ => -(C / t)) (C / t ^ 2) t := by
    intro t ht
    rw [Set.uIcc_of_le hab] at ht
    have ht0 : t ≠ 0 := by linarith [ht.1]
    have := ((hasDerivAt_inv ht0).const_mul C).neg
    convert this using 1
    · funext x; simp [div_eq_mul_inv]
    · field_simp
  have hi : IntervalIntegrable (fun t : ℝ => C / t ^ 2) volume a b := by
    apply ContinuousOn.intervalIntegrable
    apply ContinuousOn.div continuousOn_const (continuousOn_pow 2)
    intro t ht; rw [Set.uIcc_of_le hab] at ht
    have : 0 < t := by linarith [ht.1]
    positivity
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd hi]
  have : 0 ≤ C / b := div_nonneg hC (by linarith)
  linarith

section
variable (L : LiteratureInputs)
include L

/-- `|N(T) - A(T)| ≤ C log T`. -/
lemma NA_bound : ∃ C T₁ : ℝ, 4 ≤ T₁ ∧ ∀ T, T₁ ≤ T → |(zeroCount T : ℝ) - Afun T| ≤ C * Real.log T := by
  obtain ⟨C₁, T₁, hS⟩ := L.rvm
  obtain ⟨C₂, T₂, hT₂, hJ⟩ := jump_le L
  refine ⟨|C₁| + |C₂|, max T₂ (max T₁ 4), le_trans (le_max_right _ _) (le_max_right _ _), fun T hT => ?_⟩
  have hT1 : T₁ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hT
  have hT2 : T₂ ≤ T := le_trans (le_max_left _ _) hT
  have hT4 : 4 ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hT
  have hlog : 0 ≤ Real.log T := Real.log_nonneg (by linarith)
  have hS' := hS T hT1
  have hJ' := hJ T hT2
  have e : (zeroCount T : ℝ) - Afun T = ((zeroCount T : ℝ) - zeroCountAvg T) + zetaS T := by
    rw [zetaS_eq]; ring
  have hd0 : 0 ≤ (zeroCount T : ℝ) - zeroCountAvg T := by
    have := avg_le_zeroCount T; linarith
  have hd1 : (zeroCount T : ℝ) - zeroCountAvg T ≤ C₂ * Real.log T := by
    unfold zeroCountAvg
    rw [zeroCountBelow_eq]
    have h0 : 0 ≤ ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) := Finset.sum_nonneg fun _ _ => Nat.cast_nonneg _
    have : (zeroCount T : ℝ) - (∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) + zeroCount T) / 2 =
        ((zeroCount T : ℝ) - ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ)) / 2 := by ring
    rw [this]
    have hj0 : 0 ≤ (zeroCount T : ℝ) - ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) := by
      have := avg_le_zeroCount T
      unfold zeroCountAvg at this; rw [zeroCountBelow_eq] at this; linarith
    linarith
  rw [e]
  calc |(zeroCount T : ℝ) - zeroCountAvg T + zetaS T| ≤ |(zeroCount T : ℝ) - zeroCountAvg T| + |zetaS T| :=
        abs_add_le _ _
    _ ≤ |C₂| * Real.log T + |C₁| * Real.log T := by
        rw [abs_of_nonneg hd0]
        have h1 : C₂ * Real.log T ≤ |C₂| * Real.log T := mul_le_mul_of_nonneg_right (le_abs_self _) hlog
        have h2 : C₁ * Real.log T ≤ |C₁| * Real.log T := mul_le_mul_of_nonneg_right (le_abs_self _) hlog
        linarith
    _ = (|C₁| + |C₂|) * Real.log T := by ring

/-- `∫₀^Y |r| ≤ K_r` uniformly. -/
lemma rfun_L1 : ∃ Kr : ℝ, ∀ Y : ℝ, 0 ≤ Y → ∫ t in (0 : ℝ)..Y, |rfun t| ≤ Kr := by
  obtain ⟨C, t₀, ht₀, hC⟩ := rfun_bound L
  have ht0 : 0 < t₀ := by linarith
  refine ⟨(∫ t in (0 : ℝ)..t₀, |rfun t|) + |C| / t₀, fun Y hY => ?_⟩
  have hpos : 0 ≤ |C| / t₀ := by positivity
  rcases le_or_gt Y t₀ with hYt | hYt
  · have := intervalIntegral.integral_mono_interval (c := 0) (d := t₀) (a := 0) (b := Y) le_rfl hY hYt
      (Filter.Eventually.of_forall fun t => abs_nonneg (rfun t)) (intervalIntegrable_rfun le_rfl ht0.le).abs
    linarith
  · have hi1 := (intervalIntegrable_rfun le_rfl ht0.le).abs
    have hi2 := (intervalIntegrable_rfun ht0.le hY).abs
    rw [← intervalIntegral.integral_add_adjacent_intervals hi1 hi2]
    have h2 : ∫ t in t₀..Y, |rfun t| ≤ ∫ t in t₀..Y, |C| / t ^ 2 := by
      apply intervalIntegral.integral_mono_on hYt.le hi2
      · apply ContinuousOn.intervalIntegrable
        apply ContinuousOn.div continuousOn_const (continuousOn_pow 2)
        intro t ht; rw [Set.uIcc_of_le hYt.le] at ht
        have : 0 < t := by linarith [ht.1]
        positivity
      intro t ht
      have := hC t ht.1
      have ht0' : 0 < t := by linarith [ht.1]
      exact this.trans (div_le_div_of_nonneg_right (le_abs_self C) (by positivity))
    have h3 := integral_inv_sq_le' ht0 hYt.le (abs_nonneg C)
    linarith

omit L in
/-- **The `cos nθ · log` integral** (one integration by parts):
`|∫₀^T cos(nθ(t)) log(t/2π) dt| ≤ c₀ + 4(T² + 1)(log T + 3)/n`. -/
lemma cos_log_bound : ∃ c₀ : ℝ, ∀ n : ℕ, 1 ≤ n → ∀ T : ℝ, 1 ≤ T →
    |∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))| ≤
      c₀ + 4 * (T ^ 2 + 1) * (Real.log T + 3) / n := by
  have hπ := Real.pi_gt_three
  have hπ0 := Real.pi_pos
  set Lf : ℝ → ℝ := fun t => Real.log (t / (2 * Real.pi)) with hLf
  refine ⟨∫ t in (0 : ℝ)..1, |Lf t|, fun n hn T hT => ?_⟩
  have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
  have hθc : Continuous liTheta := continuous_iff_continuousAt.2 fun t => (hasDerivAt_liTheta t).continuousAt
  have hLi : ∀ a b, 0 ≤ a → 0 ≤ b → IntervalIntegrable Lf volume a b :=
    fun a b ha hb => intervalIntegrable_log_div_two_pi ha hb
  have hgi : ∀ a b, 0 ≤ a → 0 ≤ b →
      IntervalIntegrable (fun t => Real.cos (n * liTheta t) * Lf t) volume a b := fun a b ha hb =>
    (hLi a b ha hb).continuousOn_mul (Continuous.continuousOn (by fun_prop))
  rw [← intervalIntegral.integral_add_adjacent_intervals (hgi 0 1 le_rfl zero_le_one)
    (hgi 1 T zero_le_one (by linarith))]
  -- piece on [0,1]
  have h01 : |∫ t in (0 : ℝ)..1, Real.cos (n * liTheta t) * Lf t| ≤ ∫ t in (0 : ℝ)..1, |Lf t| := by
    refine (intervalIntegral.abs_integral_le_integral_abs zero_le_one).trans ?_
    apply intervalIntegral.integral_mono_on zero_le_one (hgi 0 1 le_rfl zero_le_one).abs
      (hLi 0 1 le_rfl zero_le_one).abs
    intro t _
    rw [abs_mul]
    exact mul_le_of_le_one_left (abs_nonneg _) (abs_cos_le_one _)
  -- by parts on [1, T]
  set u : ℝ → ℝ := fun t => (t ^ 2 + 1 / 4) * Lf t with hu
  set u' : ℝ → ℝ := fun t => 2 * t * Lf t + (t ^ 2 + 1 / 4) * (1 / t) with hu'
  set v : ℝ → ℝ := fun t => Real.sin (n * liTheta t) with hv
  set v' : ℝ → ℝ := fun t => Real.cos (n * liTheta t) * (n * -(1 / (t ^ 2 + 1 / 4))) with hv'
  have hmem : ∀ t ∈ Set.uIcc 1 T, 1 ≤ t ∧ t ≤ T := fun t ht => by
    rw [Set.uIcc_of_le hT] at ht; exact ht
  have hud : ∀ t ∈ Set.uIcc 1 T, HasDerivAt u (u' t) t := by
    intro t ht
    have ht0 : 0 < t := by linarith [(hmem t ht).1]
    have h1 : HasDerivAt (fun t : ℝ => t / (2 * Real.pi)) (1 / (2 * Real.pi)) t := by
      simpa using (hasDerivAt_id t).div_const (2 * Real.pi)
    have h2 := h1.log (by positivity)
    have h3 := ((hasDerivAt_pow 2 t).add_const (1 / 4)).mul h2
    simp only [hu, hu', hLf]
    convert h3 using 1
    field_simp
    ring
  have hvd : ∀ t ∈ Set.uIcc 1 T, HasDerivAt v (v' t) t := by
    intro t _
    have := ((hasDerivAt_liTheta t).const_mul (n : ℝ)).sin
    simp only [hv, hv']
    convert this using 1
  have hcu' : ContinuousOn u' (Set.uIcc 1 T) := by
    intro t ht
    have ht0 : 0 < t := by linarith [(hmem t ht).1]
    have hne : t / (2 * Real.pi) ≠ 0 := by positivity
    have hne2 : t ≠ 0 := ht0.ne'
    apply ContinuousAt.continuousWithinAt
    simp only [hu', hLf]
    fun_prop (disch := assumption)
  have hcv' : Continuous v' := by
    simp only [hv']
    exact (Real.continuous_cos.comp (continuous_const.mul hθc)).mul
      (continuous_const.mul continuous_inv_sq_add.neg)
  have hparts := intervalIntegral.integral_mul_deriv_eq_deriv_mul hud hvd hcu'.intervalIntegrable
    (hcv'.intervalIntegrable _ _)
  have e : ∫ t in (1 : ℝ)..T, Real.cos (n * liTheta t) * Lf t = -(1 / n) * ∫ t in (1 : ℝ)..T, u t * v' t := by
    rw [← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_congr fun t ht => ?_
    simp only [hu, hv']
    have hw : 0 < t ^ 2 + 1 / 4 := by positivity
    field_simp
  rw [e, hparts]
  -- bounds on [1,T]
  have hlT : 0 ≤ Real.log T := Real.log_nonneg hT
  have hl2π : Real.log (2 * Real.pi) ≤ 2 := by
    rw [Real.log_le_iff_le_exp (by positivity)]
    have he := Real.exp_one_gt_d9
    have hp := Real.pi_lt_d2
    have : Real.exp 2 = Real.exp 1 * Real.exp 1 := by rw [← Real.exp_add]; norm_num
    rw [this]; nlinarith
  have hl2π0 : 0 ≤ Real.log (2 * Real.pi) := Real.log_nonneg (by linarith)
  have hLb : ∀ t, 1 ≤ t → t ≤ T → |Lf t| ≤ Real.log T + 2 := by
    intro t h1 h2
    simp only [hLf]
    rw [Real.log_div (by linarith) (by positivity)]
    have : 0 ≤ Real.log t := Real.log_nonneg h1
    have : Real.log t ≤ Real.log T := Real.log_le_log (by linarith) h2
    rw [abs_le]; constructor <;> linarith
  have hT2 : T ≤ T ^ 2 := by nlinarith
  have huT : |u T * v T| ≤ (T ^ 2 + 1) * (Real.log T + 3) := by
    rw [abs_mul]
    have h1 : |v T| ≤ 1 := abs_sin_le_one _
    have h2 : |u T| ≤ (T ^ 2 + 1 / 4) * (Real.log T + 2) := by
      simp only [hu]; rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < T ^ 2 + 1 / 4)]
      exact mul_le_mul_of_nonneg_left (hLb T hT le_rfl) (by positivity)
    calc |u T| * |v T| ≤ (T ^ 2 + 1 / 4) * (Real.log T + 2) * 1 :=
          mul_le_mul h2 h1 (abs_nonneg _) (by positivity)
      _ ≤ (T ^ 2 + 1) * (Real.log T + 3) := by nlinarith
  have hu1 : |u 1 * v 1| ≤ (T ^ 2 + 1) * (Real.log T + 3) := by
    rw [abs_mul]
    have h1 : |v 1| ≤ 1 := abs_sin_le_one _
    have h2 : |u 1| ≤ (1 ^ 2 + 1 / 4) * (Real.log T + 2) := by
      simp only [hu]; rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 1 ^ 2 + 1 / 4)]
      exact mul_le_mul_of_nonneg_left (hLb 1 le_rfl hT) (by positivity)
    calc |u 1| * |v 1| ≤ (1 ^ 2 + 1 / 4) * (Real.log T + 2) * 1 :=
          mul_le_mul h2 h1 (abs_nonneg _) (by positivity)
      _ ≤ (T ^ 2 + 1) * (Real.log T + 3) := by nlinarith
  have hint : |∫ t in (1 : ℝ)..T, u' t * v t| ≤ 2 * (T ^ 2 + 1) * (Real.log T + 3) := by
    have hcv : Continuous v := by simp only [hv]; fun_prop
    have hi : IntervalIntegrable (fun t => u' t * v t) volume 1 T :=
      (hcu'.mul hcv.continuousOn).intervalIntegrable
    refine (intervalIntegral.abs_integral_le_integral_abs hT).trans ?_
    have hpt : ∀ t ∈ Set.Icc 1 T, |u' t * v t| ≤ 2 * T * (Real.log T + 2) + T + 1 := by
      intro t ht
      have h1 : |v t| ≤ 1 := abs_sin_le_one _
      have ht0 : 0 < t := by linarith [ht.1]
      have h2 : |u' t| ≤ 2 * T * (Real.log T + 2) + T + 1 := by
        simp only [hu']
        have hl := hLb t ht.1 ht.2
        have e1 : |2 * t * Lf t| ≤ 2 * T * (Real.log T + 2) := by
          rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 2 * t)]
          exact mul_le_mul (by linarith [ht.2]) hl (abs_nonneg _) (by positivity)
        have e2 : |(t ^ 2 + 1 / 4) * (1 / t)| ≤ T + 1 := by
          rw [abs_of_pos (by positivity)]
          have : (t ^ 2 + 1 / 4) * (1 / t) = t + 1 / (4 * t) := by field_simp
          rw [this]
          have : 1 / (4 * t) ≤ 1 := by rw [div_le_one (by positivity)]; linarith [ht.1]
          linarith [ht.2]
        calc |2 * t * Lf t + (t ^ 2 + 1 / 4) * (1 / t)| ≤ |2 * t * Lf t| + |(t ^ 2 + 1 / 4) * (1 / t)| :=
              abs_add_le _ _
          _ ≤ 2 * T * (Real.log T + 2) + (T + 1) := add_le_add e1 e2
          _ = 2 * T * (Real.log T + 2) + T + 1 := by ring
      rw [abs_mul]
      calc |u' t| * |v t| ≤ (2 * T * (Real.log T + 2) + T + 1) * 1 :=
            mul_le_mul h2 h1 (abs_nonneg _) (by positivity)
        _ = 2 * T * (Real.log T + 2) + T + 1 := by ring
    have := intervalIntegral.integral_mono_on hT hi.abs intervalIntegrable_const hpt
    rw [intervalIntegral.integral_const, smul_eq_mul] at this
    refine this.trans ?_
    have h0 : 0 ≤ T - 1 := by linarith
    have : (T - 1) * (2 * T * (Real.log T + 2) + T + 1) ≤ T * (2 * T * (Real.log T + 2) + T + 1) := by
      apply mul_le_mul_of_nonneg_right (by linarith) (by positivity)
    nlinarith
  have hfin : |-(1 / (n : ℝ)) * (u T * v T - u 1 * v 1 - ∫ t in (1 : ℝ)..T, u' t * v t)| ≤
      4 * (T ^ 2 + 1) * (Real.log T + 3) / n := by
    rw [abs_mul, abs_neg, abs_of_pos (by positivity : (0 : ℝ) < 1 / n)]
    have h3 : |u T * v T - u 1 * v 1 - ∫ t in (1 : ℝ)..T, u' t * v t| ≤ 4 * (T ^ 2 + 1) * (Real.log T + 3) := by
      calc |u T * v T - u 1 * v 1 - ∫ t in (1 : ℝ)..T, u' t * v t| ≤
            |u T * v T| + |u 1 * v 1| + |∫ t in (1 : ℝ)..T, u' t * v t| := by
            have := abs_sub (u T * v T - u 1 * v 1) (∫ t in (1 : ℝ)..T, u' t * v t)
            have := abs_sub (u T * v T) (u 1 * v 1)
            linarith
        _ ≤ (T ^ 2 + 1) * (Real.log T + 3) + (T ^ 2 + 1) * (Real.log T + 3) +
            2 * (T ^ 2 + 1) * (Real.log T + 3) := by linarith
        _ = 4 * (T ^ 2 + 1) * (Real.log T + 3) := by ring
    calc 1 / (n : ℝ) * |u T * v T - u 1 * v 1 - ∫ t in (1 : ℝ)..T, u' t * v t| ≤
          1 / n * (4 * (T ^ 2 + 1) * (Real.log T + 3)) := mul_le_mul_of_nonneg_left h3 (by positivity)
      _ = 4 * (T ^ 2 + 1) * (Real.log T + 3) / n := by ring
  calc |(∫ t in (0 : ℝ)..1, Real.cos (n * liTheta t) * Lf t) +
        -(1 / (n : ℝ)) * (u T * v T - u 1 * v 1 - ∫ t in (1 : ℝ)..T, u' t * v t)| ≤
      |∫ t in (0 : ℝ)..1, Real.cos (n * liTheta t) * Lf t| +
        |-(1 / (n : ℝ)) * (u T * v T - u 1 * v 1 - ∫ t in (1 : ℝ)..T, u' t * v t)| := abs_add_le _ _
    _ ≤ (∫ t in (0 : ℝ)..1, |Lf t|) + 4 * (T ^ 2 + 1) * (Real.log T + 3) / n := add_le_add h01 hfin

omit L in
/-- The integrand of `I(n)` in `t`. -/
lemma integrableOn_hI (n : ℕ) :
    IntegrableOn (fun t => 1 / Real.pi * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * Real.pi))))
      (Ioi 0) := by
  have hd : ∀ t ∈ Ioi (0 : ℝ), HasDerivWithinAt liTheta (-(1 / (t ^ 2 + 1 / 4))) (Ioi 0) t :=
    fun t _ => (hasDerivAt_liTheta t).hasDerivWithinAt
  have hinj : InjOn liTheta (Ioi 0) := liTheta_strictAnti.injective.injOn
  have h := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hd hinj
    (fun v => 2 * (1 - Real.cos (n * v)) * a0 v)).1 (by
      rw [liTheta_image]
      have := (intervalIntegrable_Imod_integrand n).1
      exact this.mono_set Ioo_subset_Ioc_self)
  refine h.congr_fun (fun t ht => ?_) measurableSet_Ioi
  have ht0 : 0 < t := ht
  simp only [smul_eq_mul]
  rw [abs_neg, abs_of_pos (by positivity : (0 : ℝ) < 1 / (t ^ 2 + 1 / 4))]
  have := a0_liTheta ht0
  have hπ := Real.pi_pos
  calc 1 / (t ^ 2 + 1 / 4) * (2 * (1 - Real.cos (n * liTheta t)) * a0 (liTheta t))
      = 2 * (1 - Real.cos (n * liTheta t)) * (1 / (t ^ 2 + 1 / 4) * a0 (liTheta t)) := by ring
    _ = 2 * (1 - Real.cos (n * liTheta t)) * (1 / (2 * Real.pi) * Real.log (t / (2 * Real.pi))) := by rw [this]
    _ = 1 / Real.pi * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * Real.pi))) := by field_simp

omit L in
/-- **Tail of `I(n)`**: `|∫_{t>X} π⁻¹(1 - cos nθ) log(t/2π)| ≤ 2n²/(π√X)` for `X ≥ 2π`. -/
lemma I_tail_bound (n : ℕ) {X : ℝ} (hX : 2 * Real.pi ≤ X) :
    |∫ t in Ioi X, 1 / Real.pi * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * Real.pi)))| ≤
      2 * (n : ℝ) ^ 2 / (Real.pi * √X) := by
  have hπ := Real.pi_pos
  have hπ3 := Real.pi_gt_three
  have hX0 : 0 < X := lt_of_lt_of_le (by positivity) hX
  set h := fun t => 1 / Real.pi * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * Real.pi))) with hh
  have hint : IntegrableOn h (Ioi X) := (integrableOn_hI n).mono_set (Ioi_subset_Ioi hX0.le)
  have hrp : IntegrableOn (fun t : ℝ => (n : ℝ) ^ 2 / Real.pi * t ^ (-(3 / 2 : ℝ))) (Ioi X) :=
    (integrableOn_Ioi_rpow_of_lt (by norm_num) hX0).const_mul _
  refine (abs_integral_le_integral_abs).trans ?_
  have hmono : ∫ t in Ioi X, |h t| ≤ ∫ t in Ioi X, (n : ℝ) ^ 2 / Real.pi * t ^ (-(3 / 2 : ℝ)) := by
    apply setIntegral_mono_on hint.abs hrp measurableSet_Ioi
    intro t ht
    have ht0 : 0 < t := lt_trans hX0 ht
    have htX : X < t := ht
    simp only [hh]
    have hθ := liTheta_le_inv ht0
    have hθ0 := liTheta_pos ht0
    have hc := Real.one_sub_sq_div_two_le_cos (x := n * liTheta t)
    have hc1 := Real.cos_le_one (n * liTheta t)
    have hl0 : 0 ≤ Real.log (t / (2 * Real.pi)) := Real.log_nonneg (by rw [le_div_iff₀ (by positivity)]; linarith)
    have hl1 : Real.log (t / (2 * Real.pi)) ≤ 2 * √t := by
      have := log_le_two_sqrt (show 0 < t / (2 * Real.pi) by positivity)
      have : √(t / (2 * Real.pi)) ≤ √t := Real.sqrt_le_sqrt (div_le_self ht0.le (by linarith))
      linarith
    have h1 : 1 - Real.cos (n * liTheta t) ≤ (n : ℝ) ^ 2 / (2 * t ^ 2) := by
      have : (n * liTheta t) ^ 2 ≤ (n : ℝ) ^ 2 / t ^ 2 := by
        rw [mul_pow, div_eq_mul_one_div ((n : ℝ) ^ 2)]
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        have := pow_le_pow_left₀ hθ0.le hθ 2
        rwa [div_pow, one_pow] at this
      have e : (n : ℝ) ^ 2 / (2 * t ^ 2) = (n : ℝ) ^ 2 / t ^ 2 / 2 := by ring
      rw [e]; linarith
    rw [abs_mul, abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 1 / Real.pi),
      abs_of_nonneg (by linarith : (0 : ℝ) ≤ 1 - Real.cos (n * liTheta t)), abs_of_nonneg hl0]
    have hrp' : t ^ (-(3 / 2 : ℝ)) = 1 / (t * √t) := by
      rw [Real.rpow_neg ht0.le, show (3 / 2 : ℝ) = 1 + 1 / 2 by norm_num, Real.rpow_add ht0,
        Real.rpow_one, ← Real.sqrt_eq_rpow, one_div]
    rw [hrp']
    have hs : 0 < √t := Real.sqrt_pos.2 ht0
    have hss : √t * √t = t := Real.mul_self_sqrt ht0.le
    calc 1 / Real.pi * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * Real.pi)))
        ≤ 1 / Real.pi * ((n : ℝ) ^ 2 / (2 * t ^ 2) * (2 * √t)) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          exact mul_le_mul h1 hl1 hl0 (by positivity)
      _ = (n : ℝ) ^ 2 / Real.pi * (1 / (t * √t)) := by
          field_simp
          rw [Real.sq_sqrt ht0.le]
  refine hmono.trans (le_of_eq ?_)
  rw [integral_const_mul, integral_Ioi_rpow_of_lt (by norm_num) hX0]
  have hsX : X ^ (-(3 / 2 : ℝ) + 1) = 1 / √X := by
    rw [show (-(3 / 2 : ℝ) + 1) = -(1 / 2) by norm_num, Real.rpow_neg hX0.le, ← Real.sqrt_eq_rpow, one_div]
  rw [hsX]
  have hs : 0 < √X := Real.sqrt_pos.2 hX0
  field_simp
  ring

/-- **Tail of Li's sum**: `Σ'_{γ > X} 2(1 - cos nθ_k) ≤ n²(K/X² + 10/√X)`. -/
lemma li_tail_bound : ∃ K : ℝ, 0 ≤ K ∧ ∀ (n : ℕ) (X : ℝ), 0 < X →
    ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), 2 * (1 - Real.cos (n * thK k)) ≤ (n : ℝ) ^ 2 * (K / X ^ 2 + 10 / √X) := by
  obtain ⟨K, hK0, hK⟩ := tsum_tail_le L
  refine ⟨K, hK0, fun n X hX => ?_⟩
  have hs := (summable_inv_ordK_sq L).subtype ((ZK X : Set ZeroIdx)ᶜ)
  have hpt : ∀ k : ↑((ZK X : Set ZeroIdx)ᶜ), 2 * (1 - Real.cos (n * thK k)) ≤ (n : ℝ) ^ 2 * (1 / ordK k ^ 2) := by
    intro k
    have hk0 : 0 < ordK (k : ZeroIdx) := (k : ZeroIdx).1.2.2
    have hθ := liTheta_le_inv hk0
    have hθ0 := liTheta_pos hk0
    have hc := Real.one_sub_sq_div_two_le_cos (x := n * thK k)
    have : (n * thK (k : ZeroIdx)) ^ 2 ≤ (n : ℝ) ^ 2 * (1 / ordK (k : ZeroIdx) ^ 2) := by
      rw [mul_pow]
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      have := pow_le_pow_left₀ hθ0.le hθ 2
      rwa [div_pow, one_pow] at this
    linarith
  have hnn : ∀ k : ↑((ZK X : Set ZeroIdx)ᶜ), 0 ≤ 2 * (1 - Real.cos (n * thK k)) := fun k => by
    have := Real.cos_le_one (n * thK k); linarith
  have hsum1 : Summable fun k : ↑((ZK X : Set ZeroIdx)ᶜ) => 2 * (1 - Real.cos (n * thK k)) :=
    Summable.of_nonneg_of_le hnn hpt (hs.mul_left _)
  calc ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), 2 * (1 - Real.cos (n * thK k))
      ≤ ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), (n : ℝ) ^ 2 * (1 / ordK k ^ 2) :=
        hsum1.tsum_le_tsum hpt (hs.mul_left _)
    _ = (n : ℝ) ^ 2 * ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), 1 / ordK k ^ 2 := tsum_mul_left
    _ ≤ (n : ℝ) ^ 2 * (K / X ^ 2 + 10 / √X) := mul_le_mul_of_nonneg_left (hK X hX) (by positivity)

end

end LeanLi

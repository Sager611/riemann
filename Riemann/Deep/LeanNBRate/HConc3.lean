import Riemann.Deep.LeanNBRate.HConc2
import Riemann.Deep.LeanNBRate.Contour6

/-!
# Growth of the concrete H = (1 - E_low) Φ: HBound and the bound at the trivial zeros
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem norm_eS_le {ℓ : ℝ} (z : ℂ) {u : ℝ} (hu0 : 0 ≤ u) (hul : u ≤ ℓ) :
    ‖eS z u‖ ≤ Real.exp (ℓ * max 0 (1 / 2 - z.re)) := by
  unfold eS
  rw [Complex.norm_exp, Real.exp_le_exp]
  have : (-(z - 1 / 2) * (u : ℂ)).re = (1 / 2 - z.re) * u := by simp
  rw [this]
  have h1 : (1 / 2 - z.re) * u ≤ max 0 (1 / 2 - z.re) * u := mul_le_mul_of_nonneg_right (le_max_right _ _) hu0
  have h2 : max 0 (1 / 2 - z.re) * u ≤ max 0 (1 / 2 - z.re) * ℓ :=
    mul_le_mul_of_nonneg_left hul (le_max_left _ _)
  linarith

theorem norm_eS'_eq (z : ℂ) (u : ℝ) : ‖eS' z u‖ = ‖z - 1 / 2‖ * ‖eS z u‖ := by
  unfold eS' eS; rw [norm_mul, norm_neg]

theorem Nself_eS_le {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (z : ℂ) :
    Nself ℓ (eS z) (eS' z) ≤ (ℓ * (‖z - 1 / 2‖ ^ 2 + 1 / 4) + 1) * Real.exp (2 * ℓ * max 0 (1 / 2 - z.re)) := by
  set m := max 0 (1 / 2 - z.re)
  set E := Real.exp (2 * ℓ * m)
  have hE : ∀ u : ℝ, 0 ≤ u → u ≤ ℓ → ‖eS z u‖ ^ 2 ≤ E := by
    intro u h0 h1
    have := norm_eS_le z h0 h1 (ℓ := ℓ)
    calc ‖eS z u‖ ^ 2 ≤ Real.exp (ℓ * m) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) this 2
      _ = E := by rw [← Real.exp_nat_mul]; congr 1; push_cast; ring
  have hcont : Continuous fun u : ℝ => ‖eS' z u‖ ^ 2 + 1 / 4 * ‖eS z u‖ ^ 2 :=
    ((continuous_eS' z).norm.pow 2).add (continuous_const.mul ((continuous_eS z).norm.pow 2))
  have hint : ∫ u in (0 : ℝ)..ℓ, (‖eS' z u‖ ^ 2 + 1 / 4 * ‖eS z u‖ ^ 2) ≤
      ∫ u in (0 : ℝ)..ℓ, (‖z - 1 / 2‖ ^ 2 + 1 / 4) * E := by
    refine intervalIntegral.integral_mono_on hℓ (hcont.intervalIntegrable _ _) intervalIntegrable_const
      fun u hu => ?_
    rw [norm_eS'_eq, mul_pow]
    have := hE u hu.1 hu.2
    have h0 : 0 ≤ ‖z - 1 / 2‖ ^ 2 := by positivity
    nlinarith
  rw [intervalIntegral.integral_const, smul_eq_mul, sub_zero] at hint
  have h0 := hE 0 le_rfl hℓ
  have hl := hE ℓ hℓ le_rfl
  unfold Nself
  nlinarith

theorem sqrt_Nself_eS_le {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (z : ℂ) :
    Real.sqrt (Nself ℓ (eS z) (eS' z)) ≤
      Real.sqrt (ℓ + 1) * (‖z - 1 / 2‖ + 1) * Real.exp (ℓ * max 0 (1 / 2 - z.re)) := by
  have h := Nself_eS_le hℓ z
  set x := ‖z - 1 / 2‖
  have hx : 0 ≤ x := norm_nonneg _
  have hsq : (Real.sqrt (ℓ + 1) * (x + 1) * Real.exp (ℓ * max 0 (1 / 2 - z.re))) ^ 2 =
      (ℓ + 1) * (x + 1) ^ 2 * Real.exp (2 * ℓ * max 0 (1 / 2 - z.re)) := by
    rw [mul_pow, mul_pow, Real.sq_sqrt (by linarith), ← Real.exp_nat_mul]
    congr 2; push_cast; ring
  have hle : Nself ℓ (eS z) (eS' z) ≤ (Real.sqrt (ℓ + 1) * (x + 1) * Real.exp (ℓ * max 0 (1 / 2 - z.re))) ^ 2 := by
    rw [hsq]
    refine h.trans (mul_le_mul_of_nonneg_right ?_ (Real.exp_pos _).le)
    nlinarith
  calc Real.sqrt (Nself ℓ (eS z) (eS' z)) ≤
      Real.sqrt ((Real.sqrt (ℓ + 1) * (x + 1) * Real.exp (ℓ * max 0 (1 / 2 - z.re))) ^ 2) :=
        Real.sqrt_le_sqrt hle
    _ = _ := Real.sqrt_sq (by positivity)

theorem norm_one_sub_Elow_le {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (r : ExpSum) (z : ℂ) :
    ‖1 - Elow ℓ r.eval (r.dmod 0).eval z‖ ≤ 1 + Real.sqrt (ℓ + 1) * Real.sqrt (Nself ℓ r.eval (r.dmod 0).eval) *
      (‖z - 1 / 2‖ + 1) * Real.exp (ℓ * max 0 (1 / 2 - z.re)) := by
  have h1 := norm_Elow_le hℓ r.eval (r.dmod 0).eval r.continuous_eval (r.dmod 0).continuous_eval z
  have h2 := sqrt_Nself_eS_le hℓ z
  have h3 := norm_sub_le (1 : ℂ) (Elow ℓ r.eval (r.dmod 0).eval z)
  rw [norm_one] at h3
  have hN := Real.sqrt_nonneg (Nself ℓ r.eval (r.dmod 0).eval)
  have := mul_le_mul_of_nonneg_right h2 hN
  nlinarith

theorem le_two_div_add {c x y : ℝ} (hc : 0 < c) (hx : 0 ≤ x) (hy1 : y ≤ 1) (hyx : y * x ≤ c) :
    y ≤ 2 * c / (c + x) := by
  rw [le_div_iff₀ (by linarith)]
  nlinarith

theorem norm_Mw_le_two {a : ℝ} (ha : 0 < a) (x : ℂ) :
    ‖Mw a x‖ ≤ Real.exp (a * max 0 (-x.re)) * (2 * (2 / a) / (2 / a + ‖x‖)) := by
  have hE := Real.exp_pos (a * max 0 (-x.re))
  have h1 := norm_Mw_le ha.le x
  set y := ‖Mw a x‖ / Real.exp (a * max 0 (-x.re))
  have hy1 : y ≤ 1 := by rw [div_le_one hE]; exact h1
  have hyx : y * ‖x‖ ≤ 2 / a := by
    by_cases hx : x = 0
    · rw [hx, norm_zero, mul_zero]; positivity
    · have h2 := norm_Mw_le_inv ha hx
      have hxpos : 0 < ‖x‖ := norm_pos_iff.mpr hx
      have h3 : Real.exp (-(a * x.re)) ≤ Real.exp (a * max 0 (-x.re)) := by
        rw [Real.exp_le_exp]; nlinarith [le_max_right 0 (-x.re)]
      have h4 : 1 ≤ Real.exp (a * max 0 (-x.re)) := Real.one_le_exp (by positivity)
      have h5 : ‖Mw a x‖ * (a * ‖x‖) ≤ 2 * Real.exp (a * max 0 (-x.re)) := by
        have := (le_div_iff₀ (by positivity)).mp h2
        linarith
      rw [div_mul_eq_mul_div, div_le_div_iff₀ hE ha]
      nlinarith
  have := le_two_div_add (by positivity : (0 : ℝ) < 2 / a) (norm_nonneg x) hy1 hyx
  rwa [div_le_iff₀ hE, mul_comm] at this

theorem norm_Phi_le_two {a : ℝ} (ha : 0 < a) (n : ℕ) (z : ℂ) :
    ‖Phi a n z‖ ≤ Real.exp (n * a * max 0 (1 / 2 - z.re)) * 2 ^ n *
      ((2 / a) / (2 / a + ‖z - 1 / 2‖)) ^ n := by
  unfold Phi
  rw [norm_pow]
  have h := norm_Mw_le_two ha (z - 1 / 2)
  have hre : -(z - 1 / 2).re = 1 / 2 - z.re := by simp
  rw [hre] at h
  calc ‖Mw a (z - 1 / 2)‖ ^ n ≤ (Real.exp (a * max 0 (1 / 2 - z.re)) *
        (2 * (2 / a) / (2 / a + ‖z - 1 / 2‖))) ^ n := pow_le_pow_left₀ (norm_nonneg _) h n
    _ = _ := by
        rw [mul_pow, ← Real.exp_nat_mul, mul_div_assoc, mul_pow]
        ring_nf

/-- **HBound for the concrete H.** -/
theorem hBound_Hc {ℓ a : ℝ} {n : ℕ} (hℓ : 0 ≤ ℓ) (ha : 0 < a) (hn : 7 ≤ n) (r : ExpSum) :
    HBound (Hc ℓ a n r)
      (2 ^ n * (1 + Real.sqrt (ℓ + 1) * Real.sqrt (Nself ℓ r.eval (r.dmod 0).eval) * (2 / a + 1)))
      (ℓ + n * a) (2 / a) (n - 1) := by
  have hEd : Differentiable ℂ (Elow ℓ r.eval (r.dmod 0).eval) :=
    differentiable_Elow hℓ r.eval (r.dmod 0).eval ((r.dmod 0).dmod 0).eval
      (fun x => (r.dmod 0).hasDerivAt_eval x) r.continuous_eval (r.dmod 0).continuous_eval
      ((r.dmod 0).dmod 0).continuous_eval
  refine ⟨?_, by positivity, by positivity, by positivity, by omega, ?_⟩
  · unfold Hc
    exact ((differentiable_const _).sub hEd).mul (differentiable_Phi a n)
  · intro z
    set K := Real.sqrt (ℓ + 1) * Real.sqrt (Nself ℓ r.eval (r.dmod 0).eval)
    have hK : 0 ≤ K := by positivity
    set c := 2 / a
    have hc : 0 < c := by positivity
    set x := ‖z - 1 / 2‖
    have hx : 0 ≤ x := norm_nonneg _
    set m := max 0 (1 / 2 - z.re)
    have hm : 0 ≤ m := le_max_left _ _
    have h1 := norm_one_sub_Elow_le hℓ r z
    have h2 := norm_Phi_le_two ha n z
    have hq0 : 0 ≤ c / (c + x) := by positivity
    have hq1 : c / (c + x) ≤ 1 := (div_le_one (by linarith)).mpr (by linarith)
    have hpow : (c / (c + x)) ^ n = (c / (c + x)) * (c / (c + x)) ^ (n - 1) := by
      rw [← pow_succ']; congr 1; omega
    have hEl : 1 ≤ Real.exp (ℓ * m) := Real.one_le_exp (by positivity)
    have hkey : (1 + K * (x + 1) * Real.exp (ℓ * m)) * (c / (c + x)) ≤ (1 + K * (c + 1)) * Real.exp (ℓ * m) := by
      have h3 : (x + 1) * (c / (c + x)) ≤ c + 1 := by
        rw [mul_div_assoc', div_le_iff₀ (by linarith)]; nlinarith
      have h4 : K * Real.exp (ℓ * m) * ((x + 1) * (c / (c + x))) ≤ K * Real.exp (ℓ * m) * (c + 1) :=
        mul_le_mul_of_nonneg_left h3 (by positivity)
      nlinarith
    unfold Hc
    rw [norm_mul]
    have hexp : Real.exp ((ℓ + n * a) * m) = Real.exp (n * a * m) * Real.exp (ℓ * m) := by
      rw [← Real.exp_add]; congr 1; ring
    calc ‖1 - Elow ℓ r.eval (r.dmod 0).eval z‖ * ‖Phi a n z‖
        ≤ (1 + K * (x + 1) * Real.exp (ℓ * m)) * (Real.exp (n * a * m) * 2 ^ n * (c / (c + x)) ^ n) :=
          mul_le_mul h1 h2 (norm_nonneg _) (by positivity)
      _ = 2 ^ n * Real.exp (n * a * m) * (c / (c + x)) ^ (n - 1) *
          ((1 + K * (x + 1) * Real.exp (ℓ * m)) * (c / (c + x))) := by rw [hpow]; ring
      _ ≤ 2 ^ n * Real.exp (n * a * m) * (c / (c + x)) ^ (n - 1) * ((1 + K * (c + 1)) * Real.exp (ℓ * m)) :=
          mul_le_mul_of_nonneg_left hkey (by positivity)
      _ = 2 ^ n * (1 + K * (c + 1)) * Real.exp ((ℓ + n * a) * m) * (c / (c + x)) ^ (n - 1) := by
          rw [hexp]; ring

/-- **H at the trivial zeros.** -/
theorem norm_Hc_tz {ℓ a : ℝ} {n : ℕ} (hℓ : 0 ≤ ℓ) (ha : 0 < a) (r : ExpSum) (k : ℕ) :
    ‖Hc ℓ a n r (tz k)‖ ≤ (1 + Real.sqrt (ℓ + 1) * Real.sqrt (Nself ℓ r.eval (r.dmod 0).eval)) *
      Real.exp ((ℓ + 1 + n * a) * (2 * k + 5 / 2)) := by
  set K := Real.sqrt (ℓ + 1) * Real.sqrt (Nself ℓ r.eval (r.dmod 0).eval)
  have hK : 0 ≤ K := by positivity
  have hk0 : (0 : ℝ) ≤ k := k.cast_nonneg
  have hre : (tz k).re = -2 * k - 2 := tz_re k
  have hm : max 0 (1 / 2 - (tz k).re) = 2 * k + 5 / 2 := by rw [hre, max_eq_right (by linarith)]; ring
  have hx : ‖tz k - 1 / 2‖ = 2 * k + 5 / 2 := by
    have e : tz k - 1 / 2 = ((-(2 * k + 5 / 2) : ℝ) : ℂ) := by unfold tz; push_cast; ring
    rw [e, Complex.norm_real, Real.norm_eq_abs, abs_neg, abs_of_nonneg (by linarith)]
  have h1 := norm_one_sub_Elow_le hℓ r (tz k)
  have h2 := norm_Phi_le ha.le n (tz k)
  rw [hm, hx] at h1
  rw [hm] at h2
  have h3 : (2 * k + 5 / 2 : ℝ) + 1 ≤ Real.exp (2 * k + 5 / 2) := by
    have := Real.add_one_le_exp (2 * k + 5 / 2 : ℝ); linarith
  have hE1 : 1 ≤ Real.exp ((ℓ + 1) * (2 * k + 5 / 2)) := Real.one_le_exp (by positivity)
  have h4 : 1 + K * ((2 * k + 5 / 2) + 1) * Real.exp (ℓ * (2 * k + 5 / 2)) ≤
      (1 + K) * Real.exp ((ℓ + 1) * (2 * k + 5 / 2)) := by
    have e : Real.exp ((ℓ + 1) * (2 * k + 5 / 2)) = Real.exp (2 * k + 5 / 2) * Real.exp (ℓ * (2 * k + 5 / 2)) := by
      rw [← Real.exp_add]; congr 1; ring
    rw [e]
    have h5 : K * ((2 * k + 5 / 2) + 1) * Real.exp (ℓ * (2 * k + 5 / 2)) ≤
        K * Real.exp (2 * k + 5 / 2) * Real.exp (ℓ * (2 * k + 5 / 2)) := by
      have := Real.exp_pos (ℓ * (2 * k + 5 / 2))
      have h6 := mul_le_mul_of_nonneg_left h3 hK
      nlinarith
    rw [← e]
    nlinarith
  unfold Hc
  rw [norm_mul]
  have hexp : Real.exp ((ℓ + 1 + n * a) * (2 * k + 5 / 2)) =
      Real.exp ((ℓ + 1) * (2 * k + 5 / 2)) * Real.exp (n * a * (2 * k + 5 / 2)) := by
    rw [← Real.exp_add]; congr 1; ring
  rw [hexp]
  calc ‖1 - Elow ℓ r.eval (r.dmod 0).eval (tz k)‖ * ‖Phi a n (tz k)‖
      ≤ ((1 + K) * Real.exp ((ℓ + 1) * (2 * k + 5 / 2))) * Real.exp (n * a * (2 * k + 5 / 2)) :=
        mul_le_mul (h1.trans h4) h2 (norm_nonneg _) (by positivity)
    _ = _ := by ring

end Riemann.Deep.LeanNBRate


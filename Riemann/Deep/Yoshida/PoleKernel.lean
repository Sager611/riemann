/-
Riemann/Deep/Yoshida/PoleKernel.lean — the pole term as a Fourier-side weight.

For `f` supported in `[−a, a]` the autocorrelation `k` lives in `[−2a, 2a]`, so
`pole = ∫ k(x) 2cosh(x/2) dx = (1/2π) ∫ |F(r)|² J(r) dr` with
`J(r) = ∫_{−2a}^{2a} 2cosh(x/2) cos(rx) dx = 4((1/2) sinh a cos 2ar + r cosh a sin 2ar)/(1/4 + r²)`
(Fourier inversion of `k` + Fubini). Hence, below `log 2`,
`Re Q(f) = (1/2π) ∫ |F(r)|² (2π μ(r) + J(r)) dr`: positivity reduces to a one-variable weight.
-/
import Riemann.Deep.Yoshida.Small

noncomputable section

open Complex MeasureTheory Set Filter
open scoped ComplexConjugate

namespace Fayed.Yoshida

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF

/-- Closed form of `∫_{−2a}^{2a} 2 cosh(x/2) cos(rx) dx`. -/
def J (a r : ℝ) : ℝ :=
  4 * (1 / 2 * Real.sinh a * Real.cos (2 * a * r) + r * Real.cosh a * Real.sin (2 * a * r))
    / (1 / 4 + r ^ 2)

lemma hasDerivAt_antideriv (r x : ℝ) : HasDerivAt
    (fun y => 2 * (1 / 2 * Real.sinh (y / 2) * Real.cos (r * y)
      + r * Real.cosh (y / 2) * Real.sin (r * y)) / (1 / 4 + r ^ 2))
    (2 * Real.cosh (x / 2) * Real.cos (r * x)) x := by
  have h1 : HasDerivAt (fun y : ℝ => y / 2) (1 / 2) x := by
    simpa using (hasDerivAt_id x).div_const 2
  have h2 : HasDerivAt (fun y : ℝ => r * y) r x := by
    simpa using (hasDerivAt_id x).const_mul r
  have hsh := h1.sinh
  have hch := h1.cosh
  have hc := h2.cos
  have hs := h2.sin
  have := ((((hsh.const_mul (1 / 2)).mul hc).add ((hch.const_mul r).mul hs)).const_mul 2).div_const
    (1 / 4 + r ^ 2)
  convert this using 1
  have hpos : (1 / 4 + r ^ 2) ≠ 0 := by positivity
  field_simp
  ring

lemma integral_J (a r : ℝ) :
    ∫ x in (-(2 * a))..(2 * a), 2 * Real.cosh (x / 2) * Real.cos (r * x) = J a r := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hasDerivAt_antideriv r x)
    ((by fun_prop : Continuous fun x : ℝ => 2 * Real.cosh (x / 2) * Real.cos (r * x)).intervalIntegrable _ _)]
  rw [show 2 * a / 2 = a by ring, show -(2 * a) / 2 = -a by ring, show r * (2 * a) = 2 * a * r by ring,
    show r * -(2 * a) = -(2 * a * r) by ring, Real.sinh_neg, Real.cosh_neg, Real.cos_neg, Real.sin_neg]
  unfold J
  have hpos : (1 / 4 + r ^ 2) ≠ 0 := by positivity
  field_simp
  ring

lemma J_even (a r : ℝ) : J a (-r) = J a r := by
  unfold J
  rw [show 2 * a * -r = -(2 * a * r) by ring, Real.cos_neg, Real.sin_neg]
  ring

lemma continuous_J (a : ℝ) : Continuous (J a) := by
  unfold J
  refine Continuous.div (by fun_prop) (by fun_prop) (fun r => by positivity)

/-- Lipschitz bound `|J(r) − J(p)| ≤ 8a² cosh a · |r − p|` (from the integral form). -/
lemma abs_J_sub_le {a : ℝ} (ha : 0 ≤ a) (r p : ℝ) :
    |J a r - J a p| ≤ 8 * a ^ 2 * Real.cosh a * |r - p| := by
  have hcont : ∀ q : ℝ, Continuous fun x : ℝ => 2 * Real.cosh (x / 2) * Real.cos (q * x) :=
    fun q => by fun_prop
  rw [← integral_J, ← integral_J, ← intervalIntegral.integral_sub ((hcont r).intervalIntegrable _ _)
    ((hcont p).intervalIntegrable _ _)]
  have hb : ∀ x ∈ Set.Ioc (-(2 * a)) (2 * a),
      ‖2 * Real.cosh (x / 2) * Real.cos (r * x) - 2 * Real.cosh (x / 2) * Real.cos (p * x)‖
        ≤ 2 * Real.cosh a * |r - p| * |x| := by
    intro x hx
    have hxa : |x / 2| ≤ a := by rw [abs_le]; constructor <;> linarith [hx.1, hx.2]
    have hch : Real.cosh (x / 2) ≤ Real.cosh a := by
      rw [← Real.cosh_abs (x / 2)]
      exact Real.cosh_le_cosh.2 (by rw [abs_abs, abs_of_nonneg ha]; exact hxa)
    have hcos := Real.abs_cos_sub_cos_le (r * x) (p * x)
    rw [Real.norm_eq_abs, ← mul_sub, abs_mul, abs_of_pos (by positivity : 0 < 2 * Real.cosh (x / 2))]
    have e : |r * x - p * x| = |r - p| * |x| := by rw [← sub_mul, abs_mul]
    rw [e] at hcos
    have h0 : 0 ≤ |r - p| * |x| := by positivity
    calc 2 * Real.cosh (x / 2) * |Real.cos (r * x) - Real.cos (p * x)|
        ≤ 2 * Real.cosh a * (|r - p| * |x|) :=
          mul_le_mul (by linarith) hcos (abs_nonneg _) (by positivity)
      _ = 2 * Real.cosh a * |r - p| * |x| := by ring
  have hle := intervalIntegral.norm_integral_le_of_norm_le (show -(2 * a) ≤ 2 * a by linarith)
    (Eventually.of_forall hb) (μ := volume)
    ((by fun_prop : Continuous fun x : ℝ => 2 * Real.cosh a * |r - p| * |x|).intervalIntegrable _ _)
  rw [Real.norm_eq_abs] at hle
  refine hle.trans (le_of_eq ?_)
  rw [intervalIntegral.integral_const_mul]
  have habs : ∫ x in (-(2 * a))..(2 * a), |x| = 4 * a ^ 2 := by
    rw [← intervalIntegral.integral_add_adjacent_intervals (b := 0)
      (continuous_abs.intervalIntegrable _ _) (continuous_abs.intervalIntegrable _ _)]
    have e1 : ∫ x in (-(2 * a))..0, |x| = ∫ x in (-(2 * a))..0, -x :=
      intervalIntegral.integral_congr fun x hx => by
        rw [Set.uIcc_of_le (by linarith)] at hx
        exact abs_of_nonpos hx.2
    have e2 : ∫ x in (0 : ℝ)..(2 * a), |x| = ∫ x in (0 : ℝ)..(2 * a), x :=
      intervalIntegral.integral_congr fun x hx => by
        rw [Set.uIcc_of_le (by linarith)] at hx
        exact abs_of_nonneg hx.1
    rw [e1, e2, intervalIntegral.integral_neg, integral_id, integral_id]
    ring
  rw [habs]
  ring

/-- Crude tail bound `J(r) ≥ −4 (sinh a/(2R²) + cosh a/R)` for `|r| ≥ R ≥ 1`. -/
lemma J_ge_tail {a R r : ℝ} (ha : 0 ≤ a) (hR : 1 ≤ R) (hr : R ≤ |r|) :
    -4 * (Real.sinh a / (2 * R ^ 2) + Real.cosh a / R) ≤ J a r := by
  have hJ : J a r = J a |r| := by
    rcases le_total 0 r with h | h
    · rw [abs_of_nonneg h]
    · rw [abs_of_nonpos h, J_even]
  rw [hJ]
  set t := |r| with ht
  have ht1 : 1 ≤ t := le_trans hR hr
  have hsh : 0 ≤ Real.sinh a := Real.sinh_nonneg_iff.2 ha
  have hch : 0 < Real.cosh a := Real.cosh_pos a
  unfold J
  have hden : 0 < 1 / 4 + t ^ 2 := by positivity
  rw [le_div_iff₀ hden]
  have hc := Real.neg_one_le_cos (2 * a * t)
  have hs := Real.neg_one_le_sin (2 * a * t)
  have h1 : -(Real.sinh a) ≤ Real.sinh a * Real.cos (2 * a * t) := by nlinarith
  have h2 : -(t * Real.cosh a) ≤ t * Real.cosh a * Real.sin (2 * a * t) := by
    have : 0 ≤ t * Real.cosh a := by positivity
    nlinarith
  have hR0 : 0 < R := by linarith
  -- (sinh a / 2 + t cosh a) ≤ (sinh a/(2R²) + cosh a/R) (1/4 + t²)
  have k1 : Real.sinh a / 2 ≤ Real.sinh a / (2 * R ^ 2) * (1 / 4 + t ^ 2) := by
    rw [div_mul_eq_mul_div, le_div_iff₀ (by positivity)]
    have : R ^ 2 ≤ 1 / 4 + t ^ 2 := by nlinarith
    nlinarith
  have k2 : t * Real.cosh a ≤ Real.cosh a / R * (1 / 4 + t ^ 2) := by
    rw [div_mul_eq_mul_div, le_div_iff₀ hR0]
    have : t * R ≤ 1 / 4 + t ^ 2 := by nlinarith
    nlinarith
  nlinarith

/-! ### The pole term on the Fourier side -/

/-- `Re(pole) = (1/2π) ∫ |F(r)|² J(r) dr` (Fourier inversion of the autocorrelation, Fubini). -/
theorem pole_re_eq_J {a : ℝ} (ha : 0 ≤ a) {f : ℝ → ℂ} (hf : TestFn a f) :
    (poleTerm (autocorr f)).re = 1 / (2 * Real.pi) * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * J a r := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  set k := autocorr f with hkdef
  have hP : Integrable (fun r : ℝ => paperFT k r) := integrable_paperFT_ofReal hFk
  set S := Icc (-(2 * a)) (2 * a) with hS
  set c : ℝ → ℂ := fun x => ((2 * Real.cosh (x / 2) : ℝ) : ℂ) with hcdef
  have hcc : Continuous c := by fun_prop
  -- step 1: pole = ∫ k · 2cosh(x/2)
  have h1 : poleTerm k = ∫ x, k x * c x := by
    unfold poleTerm
    rw [ft_eq_paperFT, ft_eq_paperFT]
    unfold paperFT
    rw [← integral_add (Fayed.WeilSuff.integrable_mul_cexp hk.continuous hkc _)
      (Fayed.WeilSuff.integrable_mul_cexp hk.continuous hkc _)]
    congr 1
    funext x
    rw [← mul_add]
    congr 1
    have e1 : I * (I / 2) * (x : ℂ) = ((-(x / 2) : ℝ) : ℂ) := by
      push_cast
      linear_combination ((x : ℂ) / 2) * I_sq
    have e2 : I * (-(I / 2)) * (x : ℂ) = ((x / 2 : ℝ) : ℂ) := by
      push_cast
      linear_combination (-(x : ℂ) / 2) * I_sq
    rw [e1, e2, ← Complex.ofReal_exp, ← Complex.ofReal_exp, ← Complex.ofReal_add, hcdef]
    congr 1
    rw [Real.cosh_eq]
    ring
  -- step 2: restrict to [−2a, 2a]
  have h2 : ∫ x, k x * c x = ∫ x in S, k x * c x := by
    refine (setIntegral_eq_integral_of_forall_compl_eq_zero fun x hx => ?_).symm
    have hx' : 2 * a < |x| := by
      simp only [hS, mem_Icc, not_and_or, not_le] at hx
      rcases hx with h | h
      · exact lt_abs.2 (Or.inr (by linarith))
      · exact lt_abs.2 (Or.inl h)
    rw [hkdef, autocorr_eq_zero hf.2 hx', zero_mul]
  -- step 3: Fourier inversion inside
  set Fx : ℝ → ℝ → ℂ := fun x r => paperFT k r * cexp (-I * r * x) * c x with hFx
  have hinv := paper_inversion hk.continuous (hk.continuous.integrable_of_hasCompactSupport hkc) hFk
  have h3 : ∫ x in S, k x * c x = (1 / (2 * Real.pi) : ℂ) * ∫ x in S, ∫ r, Fx x r := by
    rw [← integral_const_mul]
    refine integral_congr_ae (Eventually.of_forall fun x => ?_)
    simp only [hFx]
    rw [hinv x, mul_assoc, ← integral_mul_const]
  -- step 4: Fubini
  have hcS : Integrable c (volume.restrict S) := hcc.integrableOn_Icc
  have hint : Integrable (Function.uncurry Fx) ((volume.restrict S).prod volume) := by
    have hprod := hcS.norm.mul_prod hP.norm
    refine hprod.mono' ?_ (Eventually.of_forall fun z => ?_)
    · have m1 : AEStronglyMeasurable (fun z : ℝ × ℝ => paperFT k z.2)
          ((volume.restrict S).prod volume) := hP.aestronglyMeasurable.comp_snd
      have m2 : AEStronglyMeasurable (fun z : ℝ × ℝ => cexp (-I * z.2 * z.1) * c z.1)
          ((volume.restrict S).prod volume) := (by fun_prop : Continuous fun z : ℝ × ℝ =>
            cexp (-I * z.2 * z.1) * c z.1).aestronglyMeasurable
      refine (m1.mul m2).congr (Eventually.of_forall fun z => ?_)
      simp only [Function.uncurry, hFx, Pi.mul_apply]
      ring
    · simp only [Function.uncurry, hFx]
      have hn : ‖cexp (-I * (z.2 : ℂ) * (z.1 : ℂ))‖ = 1 := by
        rw [show -I * (z.2 : ℂ) * (z.1 : ℂ) = ((-(z.2 * z.1) : ℝ) : ℂ) * I by push_cast; ring]
        exact Complex.norm_exp_ofReal_mul_I _
      rw [norm_mul, norm_mul, hn, mul_one, mul_comm]
  have h4 : ∫ x in S, ∫ r, Fx x r = ∫ r, ∫ x in S, Fx x r := integral_integral_swap hint
  have hint6 : Integrable (fun r => ∫ x in S, Fx x r) := hint.integral_prod_right
  -- step 5: real parts
  have hre : ∀ r : ℝ, (∫ x in S, Fx x r).re = ‖paperFT f r‖ ^ 2 * J a r := by
    intro r
    have hcx : Integrable (fun x => Fx x r) (volume.restrict S) :=
      (by fun_prop : Continuous fun x : ℝ => paperFT k r * cexp (-I * r * x) * c x).integrableOn_Icc
    have hr := integral_re hcx
    simp only [RCLike.re_to_complex] at hr
    rw [← hr]
    have hpt : ∀ x : ℝ, (Fx x r).re = ‖paperFT f r‖ ^ 2 * (2 * Real.cosh (x / 2) * Real.cos (r * x)) := by
      intro x
      simp only [hFx, hcdef]
      rw [paperFT_autocorr_real hf r,
        show -I * (r : ℂ) * (x : ℂ) = ((-(r * x)) : ℝ) * I by push_cast; ring,
        mul_comm, ← mul_assoc, ← Complex.ofReal_mul, Complex.re_ofReal_mul,
        Complex.exp_ofReal_mul_I_re, Real.cos_neg]
      ring
    simp only [hpt]
    rw [integral_const_mul, integral_Icc_eq_integral_Ioc, ← intervalIntegral.integral_of_le (by linarith),
      integral_J]
  calc (poleTerm (autocorr f)).re = ((1 / (2 * Real.pi) : ℂ) * ∫ r, ∫ x in S, Fx x r).re := by
        rw [← hkdef, h1, h2, h3, h4]
    _ = 1 / (2 * Real.pi) * (∫ r, ∫ x in S, Fx x r).re := by
        rw [show (1 / (2 * Real.pi) : ℂ) = ((1 / (2 * Real.pi) : ℝ) : ℂ) by push_cast; ring,
          Complex.re_ofReal_mul]
    _ = 1 / (2 * Real.pi) * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * J a r := by
        have hr := integral_re hint6
        simp only [RCLike.re_to_complex] at hr
        rw [← hr]
        simp only [hre]

end Fayed.Yoshida

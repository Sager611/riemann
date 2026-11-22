import Riemann.Deep.LeanNBRate.LowPart
import Mathlib.Analysis.InnerProductSpace.ProdL2
import Mathlib.Analysis.InnerProductSpace.Projection.FiniteDimensional
import Mathlib.MeasureTheory.Function.L2Space

/-!
# Dual problem: the Ornstein–Uhlenbeck form and its Riesz representer
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate
open scoped InnerProductSpace

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The OU form with endpoint terms: N(f,g) = ∫_0^ℓ (f′ḡ′ + ¼fḡ) + ½(f(0)ḡ(0) + f(ℓ)ḡ(ℓ)). -/
def NForm (ℓ : ℝ) (f f' g g' : ℝ → ℂ) : ℂ :=
  oubil ℓ f f' g g' + (1 / 2 : ℂ) * (f 0 * conj (g 0) + f ℓ * conj (g ℓ))

/-- The measure on (0, ℓ]. -/
abbrev muL (ℓ : ℝ) : Measure ℝ := volume.restrict (Ioc (0 : ℝ) ℓ)

/-- The Hilbert space carrying the OU form. -/
abbrev WSp (ℓ : ℝ) : Type := WithLp 2 (WithLp 2 (Lp ℂ 2 (muL ℓ) × Lp ℂ 2 (muL ℓ)) × WithLp 2 (ℂ × ℂ))

theorem memLp_muL_of_continuous {ℓ : ℝ} {F : ℝ → ℂ} (hF : Continuous F) : MemLp F 2 (muL ℓ) := by
  obtain ⟨C, hC⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := ℓ)).exists_bound_of_continuousOn hF.continuousOn
  refine MemLp.of_bound hF.aestronglyMeasurable C ?_
  rw [ae_restrict_iff' measurableSet_Ioc]
  exact Filter.Eventually.of_forall fun x hx => hC x ⟨hx.1.le, hx.2⟩

/-- The embedding f ↦ (f′, f/2, f(0)/√2, f(ℓ)/√2). -/
def Temb (ℓ : ℝ) (F F' : ℝ → ℂ) (hF : Continuous F) (hF' : Continuous F') : WSp ℓ :=
  WithLp.toLp 2 (WithLp.toLp 2 ((memLp_muL_of_continuous (ℓ := ℓ) hF').toLp F',
      (memLp_muL_of_continuous (ℓ := ℓ) (hF.mul continuous_const)).toLp (fun u => F u * (1 / 2 : ℂ))),
    WithLp.toLp 2 (F 0 * ((Real.sqrt 2)⁻¹ : ℂ), F ℓ * ((Real.sqrt 2)⁻¹ : ℂ)))

theorem inner_Lp_muL {ℓ : ℝ} (hℓ : 0 ≤ ℓ) {F G : ℝ → ℂ} (hF : MemLp F 2 (muL ℓ)) (hG : MemLp G 2 (muL ℓ)) :
    ⟪hG.toLp G, hF.toLp F⟫_ℂ = ∫ u in (0 : ℝ)..ℓ, F u * conj (G u) := by
  rw [MeasureTheory.L2.inner_def, intervalIntegral.integral_of_le hℓ]
  refine integral_congr_ae ?_
  filter_upwards [MemLp.coeFn_toLp hF, MemLp.coeFn_toLp hG] with u h1 h2
  rw [h1, h2, RCLike.inner_apply]

theorem inner_Temb {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (F F' G G' : ℝ → ℂ) (hF : Continuous F) (hF' : Continuous F')
    (hG : Continuous G) (hG' : Continuous G') :
    ⟪Temb ℓ G G' hG hG', Temb ℓ F F' hF hF'⟫_ℂ = NForm ℓ F F' G G' := by
  unfold Temb NForm oubil
  simp only [WithLp.prod_inner_apply]
  rw [inner_Lp_muL hℓ, inner_Lp_muL hℓ]
  simp only [RCLike.inner_apply]
  have hint1 : IntervalIntegrable (fun u => F' u * conj (G' u)) volume 0 ℓ :=
    (hF'.mul (Complex.continuous_conj.comp hG')).intervalIntegrable _ _
  have hint2 : IntervalIntegrable (fun u => F u * (1 / 2 : ℂ) * conj (G u * (1 / 2 : ℂ))) volume 0 ℓ :=
    ((hF.mul continuous_const).mul (Complex.continuous_conj.comp (hG.mul continuous_const))).intervalIntegrable _ _
  have hint3 : IntervalIntegrable (fun u => (1 / 4 : ℂ) * (F u * conj (G u))) volume 0 ℓ :=
    ((hF.mul (Complex.continuous_conj.comp hG)).intervalIntegrable _ _).const_mul _
  have e2 : (∫ u in (0 : ℝ)..ℓ, F u * (1 / 2 : ℂ) * conj (G u * (1 / 2 : ℂ))) =
      ∫ u in (0 : ℝ)..ℓ, (1 / 4 : ℂ) * (F u * conj (G u)) := by
    congr 1; ext u; simp only [map_mul, map_div₀, map_one, map_ofNat]; ring
  have hs2 : ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ = 1 / 2 := by
    rw [← mul_inv, ← Complex.ofReal_mul, Real.mul_self_sqrt (by norm_num)]
    norm_num
  have e3 : F 0 * ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * conj (G 0 * ((Real.sqrt 2 : ℝ) : ℂ)⁻¹) +
      F ℓ * ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * conj (G ℓ * ((Real.sqrt 2 : ℝ) : ℂ)⁻¹) =
      (1 / 2 : ℂ) * (F 0 * conj (G 0) + F ℓ * conj (G ℓ)) := by
    rw [map_mul, map_mul, map_inv₀, Complex.conj_ofReal, ← hs2]
    ring
  rw [e2, e3, intervalIntegral.integral_add hint1 hint3]

theorem norm_sq_Temb {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (F F' : ℝ → ℂ) (hF : Continuous F) (hF' : Continuous F') :
    ‖Temb ℓ F F' hF hF'‖ ^ 2 = (NForm ℓ F F' F F').re := by
  rw [← inner_Temb hℓ F F' F F' hF hF' hF hF', norm_sq_eq_re_inner (𝕜 := ℂ)]
  rfl

/-- The evaluation vector: ⟪e, T f⟫ = f(0). -/
def evec (ℓ : ℝ) : WSp ℓ := WithLp.toLp 2 (0, WithLp.toLp 2 (((Real.sqrt 2 : ℝ) : ℂ), 0))

theorem inner_evec {ℓ : ℝ} (F F' : ℝ → ℂ) (hF : Continuous F) (hF' : Continuous F') :
    ⟪evec ℓ, Temb ℓ F F' hF hF'⟫_ℂ = F 0 := by
  unfold evec Temb
  simp only [WithLp.prod_inner_apply, inner_zero_left, zero_add, RCLike.inner_apply, map_zero,
    Complex.conj_ofReal, mul_zero, add_zero]
  have h2 : (Real.sqrt 2 : ℂ) ≠ 0 := by exact_mod_cast (Real.sqrt_pos.mpr (by norm_num : (0 : ℝ) < 2)).ne'
  field_simp

end Riemann.Deep.LeanNBRate


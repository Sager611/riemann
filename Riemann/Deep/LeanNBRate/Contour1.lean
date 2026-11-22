import Riemann.Deep.LeanNBRate.Strip
import Riemann.Deep.LeanNBRate.ZetaBd
import Mathlib.Analysis.SpecialFunctions.JapaneseBracket

/-!
# The explicit construction (U1): definitions and generic bounds

H is an entire function with the bound HBound; F = H/zeta; g(u) = (2 pi)^{-1} int F(3+iy) e^{u(3+iy)} dy;
G(s) = int_0^L g(u) e^{-us} du; Psi_s(z) = F(z) e^{-L(s-z)}/(s-z).
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- H/zeta. -/
def Fq (H : ℂ → ℂ) (z : ℂ) : ℂ := H z / riemannZeta z

/-- The kernel e^{-L x}. -/
def kap (L : ℝ) (x : ℂ) : ℂ := cexp (-(L : ℂ) * x)

def Psi0 (H : ℂ → ℂ) (s z : ℂ) : ℂ := Fq H z / (s - z)

def Psi (H : ℂ → ℂ) (L : ℝ) (s z : ℂ) : ℂ := Fq H z * kap L (s - z) / (s - z)

/-- g(u) = (2 pi)^{-1} int F(3+iy) e^{u(3+iy)} dy. -/
def gF (H : ℂ → ℂ) (u : ℝ) : ℂ := (2 * π)⁻¹ * ∫ y : ℝ, Fq H (3 + y * I) * cexp (u * (3 + y * I))

/-- G(s) = int_0^L g(u) e^{-us} du. -/
def GF (H : ℂ → ℂ) (L : ℝ) (s : ℂ) : ℂ := ∫ u in (0 : ℝ)..L, gF H u * cexp (-(u * s))

/-- Residue coefficients of F at the trivial zeros -2(k+1). -/
def rk (H : ℂ → ℂ) (k : ℕ) : ℂ := H (-2 * ((k : ℂ) + 1)) / deriv riemannZeta (-2 * ((k : ℂ) + 1))

def Rtr (H : ℂ → ℂ) (L : ℝ) (s : ℂ) : ℂ :=
  ∑' k : ℕ, rk H k * kap L (s + 2 * ((k : ℂ) + 1)) / (s + 2 * ((k : ℂ) + 1))

def Jf (H : ℂ → ℂ) (L T₂ : ℝ) (s : ℂ) : ℂ :=
  (2 * π * I)⁻¹ * (UpperU (Psi H L s) (-1) 2 T₂ - LowerU (Psi H L s) (-1) 2 T₂)

/-- The growth/decay hypothesis on H. -/
structure HBound (H : ℂ → ℂ) (M β D : ℝ) (N : ℕ) : Prop where
  diff : Differentiable ℂ H
  M_nonneg : 0 ≤ M
  β_nonneg : 0 ≤ β
  D_pos : 0 < D
  N_ge : 4 ≤ N
  bound : ∀ z, ‖H z‖ ≤ M * (1 + ‖z‖) * Real.exp (β * max 0 (1 / 2 - z.re)) * (D / (D + ‖z - 1 / 2‖)) ^ N

theorem norm_kap (L : ℝ) (x : ℂ) : ‖kap L x‖ = Real.exp (-(L * x.re)) := by
  rw [kap, Complex.norm_exp]; congr 1; simp

theorem weight_mono {D : ℝ} (hD : 0 < D) (N : ℕ) {z : ℂ} :
    (D / (D + ‖z - 1 / 2‖)) ^ N ≤ (D / (D + |z.im|)) ^ N := by
  apply pow_le_pow_left₀ (by positivity)
  apply div_le_div_of_nonneg_left hD.le (by positivity)
  have : |z.im| = |(z - 1 / 2).im| := by simp
  rw [this]
  linarith [Complex.abs_im_le_norm (z - 1 / 2)]

theorem weight_le {D : ℝ} (hD : 0 < D) {N : ℕ} (hN : 4 ≤ N) (y : ℝ) :
    (1 + |y|) ^ 2 * (D / (D + |y|)) ^ N ≤ (D / min D 1) ^ N * ((1 + |y|) ^ 2)⁻¹ := by
  set c := min D 1 with hc
  have hc0 : 0 < c := lt_min hD one_pos
  have hcD : c ≤ D := min_le_left _ _
  have hc1 : c ≤ 1 := min_le_right _ _
  have hy : 0 ≤ |y| := abs_nonneg y
  have h1 : c * (1 + |y|) ≤ D + |y| := by nlinarith
  have hp : 0 < 1 + |y| := by linarith
  have h2 : (c * (1 + |y|)) ^ N ≤ (D + |y|) ^ N := pow_le_pow_left₀ (by positivity) h1 N
  have h3 : (1 + |y|) ^ 4 ≤ (1 + |y|) ^ N := pow_le_pow_right₀ (by linarith) hN
  have hcN : 0 < c ^ N := pow_pos hc0 N
  have hA : c ^ N * (1 + |y|) ^ 4 ≤ (D + |y|) ^ N := by
    rw [mul_pow] at h2
    exact le_trans (mul_le_mul_of_nonneg_left h3 hcN.le) h2
  rw [div_pow]
  have hB : 0 < c ^ N * (1 + |y|) ^ 4 := by positivity
  calc (1 + |y|) ^ 2 * (D ^ N / (D + |y|) ^ N) ≤ (1 + |y|) ^ 2 * (D ^ N / (c ^ N * (1 + |y|) ^ 4)) := by
        gcongr
    _ = (D / c) ^ N * ((1 + |y|) ^ 2)⁻¹ := by
        rw [div_pow]; first | (field_simp; done) | (field_simp; ring)

theorem integrable_inv_one_add_abs_sq : Integrable (fun y : ℝ => ((1 + |y|) ^ 2)⁻¹) := by
  have h := integrable_one_add_norm (E := ℝ) (μ := volume) (r := 2) (by norm_num)
  refine h.congr (Filter.Eventually.of_forall fun y => ?_)
  simp only [Real.norm_eq_abs]
  rw [Real.rpow_neg (by positivity), show (2 : ℝ) = ((2 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]

/-- (G1): a line function bounded by C/(1+|y|)^2 and measurable is integrable. -/
theorem integrable_of_le_inv_sq {g : ℝ → ℂ} (hg : AEStronglyMeasurable g) {C : ℝ}
    (h : ∀ y, ‖g y‖ ≤ C * ((1 + |y|) ^ 2)⁻¹) : Integrable g :=
  (integrable_inv_one_add_abs_sq.const_mul C).mono' hg (Filter.Eventually.of_forall h)

/-- (G2): uniform decay from a C/(1+|y|)^2 bound for |y| >= 1. -/
theorem decay_of_le_inv_sq {f : ℂ → ℂ} {a b C : ℝ}
    (h : ∀ z : ℂ, a ≤ z.re → z.re ≤ b → 1 ≤ |z.im| → ‖f z‖ ≤ C * ((1 + |z.im|) ^ 2)⁻¹) :
    ∀ ε > 0, ∀ᶠ T : ℝ in atTop, ∀ x ∈ Icc a b, ‖f (x + T * I)‖ ≤ ε ∧ ‖f (x + ((-T : ℝ) : ℂ) * I)‖ ≤ ε := by
  intro ε hε
  have hC : ∀ᶠ T : ℝ in atTop, C * ((1 + T) ^ 2)⁻¹ ≤ ε := by
    filter_upwards [Filter.eventually_ge_atTop (max 1 (C / ε))] with T hT
    have hT1 : 1 ≤ T := le_of_max_le_left hT
    have hT2 : C / ε ≤ T := le_of_max_le_right hT
    have hp : 0 < (1 + T) ^ 2 := by positivity
    rw [← div_eq_mul_inv, div_le_iff₀ hp]
    rcases le_or_gt C 0 with hC0 | hC0
    · nlinarith [mul_pos hε hp]
    · have : C ≤ ε * T := by rwa [div_le_iff₀ hε, mul_comm] at hT2
      nlinarith
  filter_upwards [hC, Filter.eventually_ge_atTop 1] with T hT hT1 x hx
  have hT0 : |T| = T := abs_of_pos (by linarith)
  constructor
  · have e : ((x : ℂ) + T * I).im = T := by simp
    have := h (x + T * I) (by simpa using hx.1) (by simpa using hx.2) (by rw [e, hT0]; exact hT1)
    rw [e, hT0] at this
    exact this.trans hT
  · have e : ((x : ℂ) + ((-T : ℝ) : ℂ) * I).im = -T := by simp
    have := h (x + ((-T : ℝ) : ℂ) * I) (by simpa using hx.1) (by simpa using hx.2)
      (by rw [e, abs_neg, hT0]; exact hT1)
    rw [e, abs_neg, hT0] at this
    exact this.trans hT

end Riemann.Deep.LeanNBRate


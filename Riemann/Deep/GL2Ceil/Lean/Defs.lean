import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.Calculus.ContDiff.Deriv

/-!
# GL2Ceil (gl2lean): Bernard's functional and the expansion (E)

Bernard arXiv 1404.6432, Thm th:mm2 (TeX l. 770-783): for P ∈ ℝ[x], P(0) = 0, P(1) = 1 (l. 161) and
Q with Q(0) = 1,
  c(P,Q,r,ξ) = 1 + (1/ξ) ∫₀¹∫₀¹ e^{2rs} [d/dx (e^{rξx} Q(s+ξx) P(x+u))|_{x=0}]² du ds,
the bracket being P'(u)Q(s) + ξP(u)Q'(s) + ξrP(u)Q(s) (bracket_eq; AT21 arXiv 2004.03581 l. 209-211).
Zero count (B l. 785-795, AT21 l. 214-218): κ_f ≥ 1 - inf (1/R) log c(P,Q,R,ν/2).
DERIVATION §1.2 (E): c = (1 + w(1)²)/2 + (C_P/ξ) W₀ + ξ B_P W₁ with w = e^{Ry} Q.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

/-- Bernard's c(P,Q,r,ξ) for real P, Q, bracket expanded (AT21 l. 210). -/
def c (P Q : ℝ → ℝ) (r ξ : ℝ) : ℝ :=
  1 + (1 / ξ) * ∫ v in (0:ℝ)..1, ∫ u in (0:ℝ)..1,
    Real.exp (2 * r * v) * (deriv P u * Q v + ξ * P u * deriv Q v + ξ * r * P u * Q v) ^ 2

/-- Bernard's literal c (B l. 779-781), bracket as a derivative at x = 0. -/
def cBernard (P Q : ℝ → ℝ) (r ξ : ℝ) : ℝ :=
  1 + (1 / ξ) * ∫ s in (0:ℝ)..1, ∫ u in (0:ℝ)..1,
    Real.exp (2 * r * s) *
      (deriv (fun x => Real.exp (r * ξ * x) * Q (s + ξ * x) * P (x + u)) 0) ^ 2

/-- κ(P,Q,R,ξ) = 1 - (1/R) log c(P,Q,R,ξ) (B l. 793, AT21 l. 216). -/
def kappa (P Q : ℝ → ℝ) (R ξ : ℝ) : ℝ := 1 - (1 / R) * Real.log (c P Q R ξ)

/-- Admissible P: C¹ (polynomials included), P(0) = 0, P(1) = 1 (B l. 161, AT21 l. 194). -/
structure PAdm (P : ℝ → ℝ) : Prop where
  smooth : ContDiff ℝ 1 P
  zero : P 0 = 0
  one : P 1 = 1

/-- Real Q class (DERIVATION §1.1): C¹, Q(0) = 1, Q(y) + Q(1-y) ≡ β on [0,1]. -/
structure QAdm (Q : ℝ → ℝ) : Prop where
  smooth : ContDiff ℝ 1 Q
  zero : Q 0 = 1
  symm : ∃ β : ℝ, ∀ y ∈ Icc (0:ℝ) 1, Q y + Q (1 - y) = β

def BP (P : ℝ → ℝ) : ℝ := ∫ u in (0:ℝ)..1, P u ^ 2
def CP (P : ℝ → ℝ) : ℝ := ∫ u in (0:ℝ)..1, deriv P u ^ 2
def w (R : ℝ) (Q : ℝ → ℝ) (y : ℝ) : ℝ := Real.exp (R * y) * Q y
def wd (R : ℝ) (Q : ℝ → ℝ) (y : ℝ) : ℝ := Real.exp (R * y) * (deriv Q y + R * Q y)
def W0 (R : ℝ) (Q : ℝ → ℝ) : ℝ := ∫ y in (0:ℝ)..1, w R Q y ^ 2
def W1 (R : ℝ) (Q : ℝ → ℝ) : ℝ := ∫ y in (0:ℝ)..1, wd R Q y ^ 2

/-- J = (1 + w(1)²)/2 + α ∫w² + γ ∫w'² (DERIVATION §1.4, §2). -/
def Jag (R α γ : ℝ) (Q : ℝ → ℝ) : ℝ := (1 + w R Q 1 ^ 2) / 2 + α * W0 R Q + γ * W1 R Q

lemma bracket_eq {P Q : ℝ → ℝ} (hP : Differentiable ℝ P) (hQ : Differentiable ℝ Q) (r ξ s u : ℝ) :
    deriv (fun x => Real.exp (r * ξ * x) * Q (s + ξ * x) * P (x + u)) 0
      = deriv P u * Q s + ξ * P u * deriv Q s + ξ * r * P u * Q s := by
  have hl : HasDerivAt (fun x : ℝ => r * ξ * x) (r * ξ * 1) 0 :=
    (hasDerivAt_id' (x := (0:ℝ))).const_mul (r * ξ)
  have h1 := hl.exp
  have hq : HasDerivAt (fun x : ℝ => s + ξ * x) (ξ * 1) 0 :=
    ((hasDerivAt_id' (x := (0:ℝ))).const_mul ξ).const_add s
  have h2 : HasDerivAt (fun x => Q (s + ξ * x)) (deriv Q (s + ξ * 0) * (ξ * 1)) 0 :=
    (hQ (s + ξ * 0)).hasDerivAt.comp (0:ℝ) hq
  have hp : HasDerivAt (fun x : ℝ => x + u) 1 0 := (hasDerivAt_id' (x := (0:ℝ))).add_const u
  have h3 : HasDerivAt (fun x => P (x + u)) (deriv P (0 + u) * 1) 0 :=
    (hP (0 + u)).hasDerivAt.comp (0:ℝ) hp
  rw [((h1.fun_mul h2).fun_mul h3).deriv]
  simp only [mul_zero, add_zero, zero_add, mul_one, Real.exp_zero, one_mul]
  ring

theorem cBernard_eq {P Q : ℝ → ℝ} (hP : Differentiable ℝ P) (hQ : Differentiable ℝ Q) (r ξ : ℝ) :
    cBernard P Q r ξ = c P Q r ξ := by
  unfold cBernard c
  simp_rw [bracket_eq hP hQ r ξ]

lemma hasDerivAt_w {R : ℝ} {Q : ℝ → ℝ} (hQ : Differentiable ℝ Q) (y : ℝ) :
    HasDerivAt (w R Q) (wd R Q y) y := by
  have h1 : HasDerivAt (fun y : ℝ => R * y) (R * 1) y := (hasDerivAt_id' (x := y)).const_mul R
  have h2 := h1.exp.fun_mul (hQ y).hasDerivAt
  show HasDerivAt (fun y => Real.exp (R * y) * Q y) _ y
  convert h2 using 1
  simp only [wd]
  ring

lemma w_cont {R : ℝ} {Q : ℝ → ℝ} (hQ : Continuous Q) : Continuous (w R Q) := by
  unfold w; fun_prop

lemma wd_cont {R : ℝ} {Q : ℝ → ℝ} (hQ : Continuous Q) (hQ' : Continuous (deriv Q)) :
    Continuous (wd R Q) := by
  unfold wd; fun_prop

lemma integral_lin3 (A B C : ℝ) {f g h : ℝ → ℝ} (hf : Continuous f) (hg : Continuous g)
    (hh : Continuous h) :
    ∫ u in (0:ℝ)..1, (A * f u + B * g u + C * h u) =
      A * (∫ u in (0:ℝ)..1, f u) + B * (∫ u in (0:ℝ)..1, g u) + C * ∫ u in (0:ℝ)..1, h u := by
  have i1 := (hf.intervalIntegrable (μ := volume) 0 1).const_mul A
  have i2 := (hg.intervalIntegrable (μ := volume) 0 1).const_mul B
  have i3 := (hh.intervalIntegrable (μ := volume) 0 1).const_mul C
  rw [intervalIntegral.integral_add (i1.add i2) i3, intervalIntegral.integral_add i1 i2,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul]

/-- ∫₀¹ P P' = 1/2. -/
lemma int_PdP {P : ℝ → ℝ} (hP : PAdm P) : ∫ u in (0:ℝ)..1, P u * deriv P u = 1 / 2 := by
  obtain ⟨hPd, hPc⟩ := contDiff_one_iff_deriv.mp hP.smooth
  have h : ∀ x ∈ uIcc (0:ℝ) 1, HasDerivAt (fun u => P u * P u / 2) (P x * deriv P x) x := by
    intro x _
    have := ((hPd x).hasDerivAt.fun_mul (hPd x).hasDerivAt).div_const 2
    convert this using 1
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt h
    ((hPd.continuous.fun_mul hPc).intervalIntegrable 0 1)]
  simp only [hP.zero, hP.one]
  norm_num

/-- The inner integral (over u) for fixed v, with q = Q v, q' = Q'(v). -/
lemma inner_expand {P : ℝ → ℝ} (hP : PAdm P) (R ξ q q' v : ℝ) :
    ∫ u in (0:ℝ)..1, Real.exp (2 * R * v) * (deriv P u * q + ξ * P u * q' + ξ * R * P u * q) ^ 2
      = CP P * (Real.exp (R * v) * q) ^ 2
        + ξ * ((Real.exp (R * v) * q) * (Real.exp (R * v) * (q' + R * q)))
        + ξ ^ 2 * BP P * (Real.exp (R * v) * (q' + R * q)) ^ 2 := by
  obtain ⟨hPd, hPc⟩ := contDiff_one_iff_deriv.mp hP.smooth
  have hP0 := hPd.continuous
  have hexp : Real.exp (2 * R * v) = Real.exp (R * v) ^ 2 := by
    rw [sq, ← Real.exp_add]; congr 1; ring
  calc (∫ u in (0:ℝ)..1, Real.exp (2 * R * v) * (deriv P u * q + ξ * P u * q' + ξ * R * P u * q) ^ 2)
      = ∫ u in (0:ℝ)..1, ((Real.exp (R * v) * q) ^ 2 * deriv P u ^ 2
          + (2 * ξ * (Real.exp (R * v) * q) * (Real.exp (R * v) * (q' + R * q))) * (P u * deriv P u)
          + (ξ ^ 2 * (Real.exp (R * v) * (q' + R * q)) ^ 2) * P u ^ 2) := by
        apply intervalIntegral.integral_congr
        intro u _
        simp only [hexp]
        ring
    _ = (Real.exp (R * v) * q) ^ 2 * (∫ u in (0:ℝ)..1, deriv P u ^ 2)
          + (2 * ξ * (Real.exp (R * v) * q) * (Real.exp (R * v) * (q' + R * q)))
            * (∫ u in (0:ℝ)..1, P u * deriv P u)
          + (ξ ^ 2 * (Real.exp (R * v) * (q' + R * q)) ^ 2) * ∫ u in (0:ℝ)..1, P u ^ 2 :=
        integral_lin3 _ _ _ (hPc.fun_pow 2) (hP0.fun_mul hPc) (hP0.fun_pow 2)
    _ = _ := by
        rw [int_PdP hP]
        unfold CP BP
        ring

/-- Expansion (E), DERIVATION §1.2. -/
theorem c_expand {P Q : ℝ → ℝ} (hP : PAdm P) (hQ : ContDiff ℝ 1 Q) (hQ0 : Q 0 = 1) (R : ℝ) {ξ : ℝ}
    (hξ : ξ ≠ 0) :
    c P Q R ξ = (1 + w R Q 1 ^ 2) / 2 + (CP P / ξ) * W0 R Q + ξ * BP P * W1 R Q := by
  obtain ⟨hQd, hQc⟩ := contDiff_one_iff_deriv.mp hQ
  have hw := w_cont (R := R) hQd.continuous
  have hwd := wd_cont (R := R) hQd.continuous hQc
  have hww : ∫ y in (0:ℝ)..1, w R Q y * wd R Q y = (w R Q 1 ^ 2 - 1) / 2 := by
    have h : ∀ x ∈ uIcc (0:ℝ) 1,
        HasDerivAt (fun y => w R Q y * w R Q y / 2) (w R Q x * wd R Q x) x := by
      intro x _
      have := ((hasDerivAt_w (R := R) hQd x).fun_mul (hasDerivAt_w (R := R) hQd x)).div_const 2
      convert this using 1
      ring
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt h ((hw.fun_mul hwd).intervalIntegrable 0 1)]
    simp only [w, hQ0, mul_zero, Real.exp_zero, mul_one]
    ring
  have hin : ∀ v, (∫ u in (0:ℝ)..1, Real.exp (2 * R * v) *
      (deriv P u * Q v + ξ * P u * deriv Q v + ξ * R * P u * Q v) ^ 2)
      = CP P * w R Q v ^ 2 + ξ * (w R Q v * wd R Q v) + ξ ^ 2 * BP P * wd R Q v ^ 2 := by
    intro v
    rw [inner_expand hP R ξ (Q v) (deriv Q v) v]
    rfl
  unfold c
  simp_rw [hin]
  rw [integral_lin3 (CP P) ξ (ξ ^ 2 * BP P) (hw.fun_pow 2) (hw.fun_mul hwd) (hwd.fun_pow 2), hww]
  unfold W0 W1
  field_simp
  ring

lemma W0_pos {R : ℝ} {Q : ℝ → ℝ} (hQ : Continuous Q) (hQ0 : Q 0 = 1) : 0 < W0 R Q := by
  unfold W0
  apply intervalIntegral.integral_pos (by norm_num : (0:ℝ) < 1) ((w_cont hQ).fun_pow 2).continuousOn
  · intro x _
    positivity
  · exact ⟨0, ⟨le_refl 0, zero_le_one⟩, by simp [w, hQ0]⟩

lemma W1_nonneg (R : ℝ) (Q : ℝ → ℝ) : 0 ≤ W1 R Q :=
  intervalIntegral.integral_nonneg zero_le_one (fun _ _ => sq_nonneg _)

lemma W0_nonneg (R : ℝ) (Q : ℝ → ℝ) : 0 ≤ W0 R Q :=
  intervalIntegral.integral_nonneg zero_le_one (fun _ _ => sq_nonneg _)

end GL2Ceil



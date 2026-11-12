import Riemann.Deep.GL2Ceil.Lean.Main
import Mathlib.Analysis.Complex.RealDeriv

/-!
# GL2Ceil (gl2lean): Bernard's complex Q (DERIVATION §1.1)

B l. 786-788: Q(x) = 1 + Σ_{n=1}^M i^{n+1} λ_n [(1-2x)^n - 1], λ ∈ ℝ^M; c with |·|² (AT21 l. 210).
For real P the integrand |A_r + iA_i|² = A_r² + A_i² splits, so c(P, Q) ≥ c(P, Re Q); Re Q is a
real polynomial with Re Q(0) = 1 and Re Q(y) + Re Q(1-y) constant (the even-n coefficients
Re i^{n+1} vanish, odd n are antisymmetric under y ↦ 1-y). Hence the ceiling holds for Bernard's
literal class.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

/-- c for complex Q, modulus squared (AT21 l. 210; B l. 779-781). -/
def cC (P : ℝ → ℝ) (Q : ℝ → ℂ) (r ξ : ℝ) : ℝ :=
  1 + (1 / ξ) * ∫ v in (0:ℝ)..1, ∫ u in (0:ℝ)..1, Real.exp (2 * r * v) *
    ‖((deriv P u : ℝ) : ℂ) * Q v + (ξ : ℂ) * (P u : ℂ) * deriv Q v
      + (ξ : ℂ) * (r : ℂ) * (P u : ℂ) * Q v‖ ^ 2

def kappaC (P : ℝ → ℝ) (Q : ℝ → ℂ) (R ξ : ℝ) : ℝ := 1 - (1 / R) * Real.log (cC P Q R ξ)

/-- Bernard's Q (B l. 786-788). -/
def QB (M : ℕ) (lam : ℕ → ℝ) (x : ℝ) : ℂ :=
  1 + ∑ n ∈ Finset.Icc 1 M, Complex.I ^ (n + 1) * (lam n : ℂ) * ((1 - 2 * (x : ℂ)) ^ n - 1)

lemma deriv_re {Q : ℝ → ℂ} (hQ : Differentiable ℝ Q) (v : ℝ) :
    deriv (fun x => (Q x).re) v = (deriv Q v).re := by
  have h := Complex.reCLM.hasFDerivAt.comp_hasDerivAt v (hQ v).hasDerivAt
  have e : (fun x => (Q x).re) = ⇑Complex.reCLM ∘ Q := by funext x; simp
  rw [e, h.deriv]
  simp

lemma deriv_im {Q : ℝ → ℂ} (hQ : Differentiable ℝ Q) (v : ℝ) :
    deriv (fun x => (Q x).im) v = (deriv Q v).im := by
  have h := Complex.imCLM.hasFDerivAt.comp_hasDerivAt v (hQ v).hasDerivAt
  have e : (fun x => (Q x).im) = ⇑Complex.imCLM ∘ Q := by funext x; simp
  rw [e, h.deriv]
  simp

lemma normsq_split (a b c d : ℝ) (q q' : ℂ) :
    ‖(a : ℂ) * q + (b : ℂ) * (c : ℂ) * q' + (b : ℂ) * (d : ℂ) * (c : ℂ) * q‖ ^ 2
      = (a * q.re + b * c * q'.re + b * d * c * q.re) ^ 2
        + (a * q.im + b * c * q'.im + b * d * c * q.im) ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp only [Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im, Complex.ofReal_re,
    Complex.ofReal_im]
  ring

/-- Inner integral with constants q = Q(v), q' = Q'(v). -/
def Gin (P : ℝ → ℝ) (r ξ q q' v : ℝ) : ℝ :=
  ∫ u in (0:ℝ)..1, Real.exp (2 * r * v) * (deriv P u * q + ξ * P u * q' + ξ * r * P u * q) ^ 2

lemma Gin_cont {P : ℝ → ℝ} (hP : PAdm P) (r ξ : ℝ) {f g : ℝ → ℝ} (hf : Continuous f)
    (hg : Continuous g) : Continuous (fun v => Gin P r ξ (f v) (g v) v) := by
  have e : (fun v => Gin P r ξ (f v) (g v) v) = fun v => CP P * (Real.exp (r * v) * f v) ^ 2
      + ξ * ((Real.exp (r * v) * f v) * (Real.exp (r * v) * (g v + r * f v)))
      + ξ ^ 2 * BP P * (Real.exp (r * v) * (g v + r * f v)) ^ 2 :=
    funext fun v => inner_expand hP r ξ (f v) (g v) v
  rw [e]
  fun_prop

lemma Gin_nonneg (P : ℝ → ℝ) (r ξ q q' v : ℝ) : 0 ≤ Gin P r ξ q q' v :=
  intervalIntegral.integral_nonneg zero_le_one (fun _ _ => by positivity)

/-- §1.1: complex Q never helps for real P: c(P, Re Q) ≤ c(P, Q). -/
theorem cC_ge {P : ℝ → ℝ} (hP : PAdm P) {Q : ℝ → ℂ} (hQ : Differentiable ℝ Q)
    (hr : ContDiff ℝ 1 (fun v => (Q v).re)) (hi : ContDiff ℝ 1 (fun v => (Q v).im))
    (r : ℝ) {ξ : ℝ} (hξ : 0 < ξ) : c P (fun v => (Q v).re) r ξ ≤ cC P Q r ξ := by
  obtain ⟨hrd, hrc⟩ := contDiff_one_iff_deriv.mp hr
  obtain ⟨hid, hic⟩ := contDiff_one_iff_deriv.mp hi
  obtain ⟨hPd, hPc⟩ := contDiff_one_iff_deriv.mp hP.smooth
  have hP0 := hPd.continuous
  have hdr : deriv (fun v => (Q v).re) = fun v => (deriv Q v).re := funext fun v => deriv_re hQ v
  have hdi : deriv (fun v => (Q v).im) = fun v => (deriv Q v).im := funext fun v => deriv_im hQ v
  have hcR : Continuous (fun v => (deriv Q v).re) := hdr ▸ hrc
  have hcI : Continuous (fun v => (deriv Q v).im) := hdi ▸ hic
  have hpt : ∀ v, (∫ u in (0:ℝ)..1, Real.exp (2 * r * v) *
      ‖((deriv P u : ℝ) : ℂ) * Q v + (ξ : ℂ) * (P u : ℂ) * deriv Q v
        + (ξ : ℂ) * (r : ℂ) * (P u : ℂ) * Q v‖ ^ 2)
      = Gin P r ξ (Q v).re (deriv Q v).re v + Gin P r ξ (Q v).im (deriv Q v).im v := by
    intro v
    unfold Gin
    have i1 : IntervalIntegrable (fun u => Real.exp (2 * r * v) *
        (deriv P u * (Q v).re + ξ * P u * (deriv Q v).re + ξ * r * P u * (Q v).re) ^ 2) volume 0 1 :=
      (by fun_prop : Continuous fun u => Real.exp (2 * r * v) *
        (deriv P u * (Q v).re + ξ * P u * (deriv Q v).re + ξ * r * P u * (Q v).re) ^ 2).intervalIntegrable 0 1
    have i2 : IntervalIntegrable (fun u => Real.exp (2 * r * v) *
        (deriv P u * (Q v).im + ξ * P u * (deriv Q v).im + ξ * r * P u * (Q v).im) ^ 2) volume 0 1 :=
      (by fun_prop : Continuous fun u => Real.exp (2 * r * v) *
        (deriv P u * (Q v).im + ξ * P u * (deriv Q v).im + ξ * r * P u * (Q v).im) ^ 2).intervalIntegrable 0 1
    rw [← intervalIntegral.integral_add i1 i2]
    apply intervalIntegral.integral_congr
    intro u _
    simp only
    rw [normsq_split]
    ring
  have hc : c P (fun v => (Q v).re) r ξ
      = 1 + (1 / ξ) * ∫ v in (0:ℝ)..1, Gin P r ξ (Q v).re (deriv Q v).re v := by
    unfold c Gin
    rw [hdr]
  have hint : (∫ v in (0:ℝ)..1, ∫ u in (0:ℝ)..1, Real.exp (2 * r * v) *
      ‖((deriv P u : ℝ) : ℂ) * Q v + (ξ : ℂ) * (P u : ℂ) * deriv Q v
        + (ξ : ℂ) * (r : ℂ) * (P u : ℂ) * Q v‖ ^ 2)
      = (∫ v in (0:ℝ)..1, Gin P r ξ (Q v).re (deriv Q v).re v)
        + ∫ v in (0:ℝ)..1, Gin P r ξ (Q v).im (deriv Q v).im v := by
    simp_rw [hpt]
    exact intervalIntegral.integral_add
      ((Gin_cont hP r ξ hrd.continuous hcR).intervalIntegrable 0 1)
      ((Gin_cont hP r ξ hid.continuous hcI).intervalIntegrable 0 1)
  have hnn : 0 ≤ ∫ v in (0:ℝ)..1, Gin P r ξ (Q v).im (deriv Q v).im v :=
    intervalIntegral.integral_nonneg zero_le_one (fun v _ => Gin_nonneg _ _ _ _ _ _)
  have h1 : 0 ≤ (1 / ξ) * ∫ v in (0:ℝ)..1, Gin P r ξ (Q v).im (deriv Q v).im v :=
    mul_nonneg (by positivity) hnn
  rw [hc]
  unfold cC
  rw [hint]
  linarith

lemma c_pos {P Q : ℝ → ℝ} (hP : PAdm P) (hQ : ContDiff ℝ 1 Q) (hQ0 : Q 0 = 1) (R : ℝ) {ξ : ℝ}
    (hξ : 0 < ξ) : 0 < c P Q R ξ := by
  rw [c_expand hP hQ hQ0 R hξ.ne']
  have hCP : 0 ≤ CP P := intervalIntegral.integral_nonneg zero_le_one (fun _ _ => sq_nonneg _)
  have hBP : 0 ≤ BP P := intervalIntegral.integral_nonneg zero_le_one (fun _ _ => sq_nonneg _)
  have h1 : 0 ≤ CP P / ξ * W0 R Q := mul_nonneg (div_nonneg hCP hξ.le) (W0_nonneg R Q)
  have h2 : 0 ≤ ξ * BP P * W1 R Q := mul_nonneg (mul_nonneg hξ.le hBP) (W1_nonneg R Q)
  have h3 : 0 < (1 + w R Q 1 ^ 2) / 2 := by positivity
  linarith

/-- Ceiling for complex Q whose real part is in the class (any C¹ imaginary part). -/
theorem ceiling_complex (cells : List Cell) (hcov : Cover cells) (hcert : CertCells cells)
    (hS : TailSmallR) (hL : TailLargeLam) (hR : TailLargeR)
    {P : ℝ → ℝ} (hP : PAdm P) {Q : ℝ → ℂ} (hQ : Differentiable ℝ Q)
    (hr : QAdm (fun v => (Q v).re)) (hi : ContDiff ℝ 1 (fun v => (Q v).im))
    {R ξ : ℝ} (hR0 : 0 < R) (hξ : 0 < ξ) (hξ0 : ξ ≤ 25 / 256) :
    kappaC P Q R ξ ≤ 0.07005633 := by
  have h1 := ceiling cells hcov hcert hS hL hR hP hr hR0 hξ hξ0
  have h2 := cC_ge hP hQ hr.smooth hi R hξ
  have hpos := c_pos hP hr.smooth hr.zero R hξ
  have hlog := Real.log_le_log hpos h2
  have h3 : 1 / R * Real.log (c P (fun v => (Q v).re) R ξ) ≤ 1 / R * Real.log (cC P Q R ξ) :=
    mul_le_mul_of_nonneg_left hlog (by positivity)
  unfold kappaC
  unfold kappa at h1
  linarith

/-! ## Bernard's literal class -/

lemma QB_term (n : ℕ) (a x : ℝ) :
    Complex.I ^ (n + 1) * (a : ℂ) * ((1 - 2 * (x : ℂ)) ^ n - 1)
      = Complex.I ^ (n + 1) * ((a * ((1 - 2 * x) ^ n - 1) : ℝ) : ℂ) := by
  push_cast
  ring

lemma QB_re (M : ℕ) (lam : ℕ → ℝ) (x : ℝ) :
    (QB M lam x).re = 1 + ∑ n ∈ Finset.Icc 1 M,
      ((Complex.I ^ (n + 1)).re * lam n) * ((1 - 2 * x) ^ n - 1) := by
  unfold QB
  rw [Complex.add_re, Complex.one_re, Complex.re_sum]
  congr 1
  apply Finset.sum_congr rfl
  intro n _
  rw [QB_term, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im]
  ring

lemma QB_im (M : ℕ) (lam : ℕ → ℝ) (x : ℝ) :
    (QB M lam x).im = ∑ n ∈ Finset.Icc 1 M,
      ((Complex.I ^ (n + 1)).im * lam n) * ((1 - 2 * x) ^ n - 1) := by
  unfold QB
  rw [Complex.add_im, Complex.one_im, zero_add, Complex.im_sum]
  apply Finset.sum_congr rfl
  intro n _
  rw [QB_term, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im]
  ring

/-- Re Q as a real polynomial. -/
def QBr (M : ℕ) (lam : ℕ → ℝ) : Polynomial ℝ :=
  Polynomial.C 1 + ∑ n ∈ Finset.Icc 1 M, Polynomial.C ((Complex.I ^ (n + 1)).re * lam n) *
    ((Polynomial.C 1 - Polynomial.C 2 * Polynomial.X) ^ n - Polynomial.C 1)

/-- Im Q as a real polynomial. -/
def QBi (M : ℕ) (lam : ℕ → ℝ) : Polynomial ℝ :=
  ∑ n ∈ Finset.Icc 1 M, Polynomial.C ((Complex.I ^ (n + 1)).im * lam n) *
    ((Polynomial.C 1 - Polynomial.C 2 * Polynomial.X) ^ n - Polynomial.C 1)

lemma QB_re_eval (M : ℕ) (lam : ℕ → ℝ) (x : ℝ) : (QB M lam x).re = (QBr M lam).eval x := by
  rw [QB_re]
  simp only [QBr, Polynomial.eval_add, Polynomial.eval_C, Polynomial.eval_finsetSum,
    Polynomial.eval_mul, Polynomial.eval_pow, Polynomial.eval_sub, Polynomial.eval_X]

lemma QB_im_eval (M : ℕ) (lam : ℕ → ℝ) (x : ℝ) : (QB M lam x).im = (QBi M lam).eval x := by
  rw [QB_im]
  simp only [QBi, Polynomial.eval_C, Polynomial.eval_finsetSum,
    Polynomial.eval_mul, Polynomial.eval_pow, Polynomial.eval_sub, Polynomial.eval_X]

lemma QB_eq (M : ℕ) (lam : ℕ → ℝ) :
    QB M lam = fun x => (((QBr M lam).eval x : ℝ) : ℂ) + (((QBi M lam).eval x : ℝ) : ℂ) * Complex.I := by
  funext x
  apply Complex.ext
  · simp [QB_re_eval]
  · simp [QB_im_eval]

lemma QB_diff (M : ℕ) (lam : ℕ → ℝ) : Differentiable ℝ (QB M lam) := by
  rw [QB_eq]
  intro x
  exact (((QBr M lam).hasDerivAt x).ofReal_comp.add
    (((QBi M lam).hasDerivAt x).ofReal_comp.mul_const Complex.I)).differentiableAt

lemma I_pow_re_even {n : ℕ} (hn : Even n) : (Complex.I ^ (n + 1)).re = 0 := by
  obtain ⟨k, rfl⟩ := hn
  rw [pow_succ, show k + k = 2 * k by ring, pow_mul, Complex.I_sq]
  have : ((-1 : ℂ) ^ k) = (((-1 : ℝ) ^ k : ℝ) : ℂ) := by push_cast; rfl
  rw [this, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
  ring

lemma term_symm (n : ℕ) (c y : ℝ) (hc : Even n → c = 0) :
    c * ((1 - 2 * y) ^ n - 1) + c * ((1 - 2 * (1 - y)) ^ n - 1)
      = c * ((1 - 2 * 0) ^ n - 1) + c * ((1 - 2 * 1) ^ n - 1) := by
  rcases Nat.even_or_odd n with h | h
  · rw [hc h]
    ring
  · have e1 : (1 - 2 * (1 - y)) ^ n = -(1 - 2 * y) ^ n := by
      rw [show (1 - 2 * (1 - y)) = -(1 - 2 * y) by ring, h.neg_pow]
    have e2 : ((1:ℝ) - 2 * 1) ^ n = -1 := by
      rw [show ((1:ℝ) - 2 * 1) = -1 by norm_num]
      exact h.neg_one_pow
    rw [e1, e2, show ((1:ℝ) - 2 * 0) = 1 by norm_num, one_pow]
    ring

lemma QB_re_symm (M : ℕ) (lam : ℕ → ℝ) (y : ℝ) :
    (QB M lam y).re + (QB M lam (1 - y)).re = (QB M lam 0).re + (QB M lam 1).re := by
  rw [QB_re, QB_re, QB_re, QB_re]
  have h : ∀ n ∈ Finset.Icc 1 M,
      ((Complex.I ^ (n + 1)).re * lam n) * ((1 - 2 * y) ^ n - 1)
        + ((Complex.I ^ (n + 1)).re * lam n) * ((1 - 2 * (1 - y)) ^ n - 1)
      = ((Complex.I ^ (n + 1)).re * lam n) * ((1 - 2 * 0) ^ n - 1)
        + ((Complex.I ^ (n + 1)).re * lam n) * ((1 - 2 * 1) ^ n - 1) :=
    fun n _ => term_symm n _ y (fun he => by rw [I_pow_re_even he, zero_mul])
  have := Finset.sum_congr rfl h
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib] at this
  linarith

lemma QB_re_zero (M : ℕ) (lam : ℕ → ℝ) : (QB M lam 0).re = 1 := by
  rw [QB_re]
  simp

lemma QBr_adm (M : ℕ) (lam : ℕ → ℝ) : QAdm (fun v => (QB M lam v).re) where
  smooth := by
    have e : (fun v => (QB M lam v).re) = fun v => (QBr M lam).eval v :=
      funext fun v => QB_re_eval M lam v
    rw [e]
    exact poly_C1 _
  zero := QB_re_zero M lam
  symm := ⟨(QB M lam 0).re + (QB M lam 1).re, fun y _ => QB_re_symm M lam y⟩

lemma QBi_C1 (M : ℕ) (lam : ℕ → ℝ) : ContDiff ℝ 1 (fun v => (QB M lam v).im) := by
  have e : (fun v => (QB M lam v).im) = fun v => (QBi M lam).eval v :=
    funext fun v => QB_im_eval M lam v
  rw [e]
  exact poly_C1 _

/-- Main theorem for Bernard's literal class (B l. 161, 770-795; AT21 l. 192-218, 227-231):
P ∈ ℝ[x] with P(0) = 0, P(1) = 1, Q = 1 + Σ_{n=1}^M i^{n+1} λ_n[(1-2x)^n - 1], λ real, any M,
R > 0, 0 < ξ = ν/2 ≤ 25/256: 1 - (1/R) log c(P,Q,R,ξ) ≤ 0.07005633. -/
theorem ceiling_bernard (cells : List Cell) (hcov : Cover cells) (hcert : CertCells cells)
    (hS : TailSmallR) (hL : TailLargeLam) (hR : TailLargeR)
    (P : Polynomial ℝ) (hP0 : P.eval 0 = 0) (hP1 : P.eval 1 = 1) (M : ℕ) (lam : ℕ → ℝ)
    {R ξ : ℝ} (hR0 : 0 < R) (hξ : 0 < ξ) (hξ0 : ξ ≤ 25 / 256) :
    kappaC (fun x => P.eval x) (QB M lam) R ξ ≤ 0.07005633 :=
  ceiling_complex cells hcov hcert hS hL hR ⟨poly_C1 P, hP0, hP1⟩ (QB_diff M lam)
    (QBr_adm M lam) (QBi_C1 M lam) hR0 hξ hξ0

end GL2Ceil



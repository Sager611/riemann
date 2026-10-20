/-
Riemann/Deep/Yoshida/Moment.lean — moment (Galerkin) certificate for Weil positivity.

Below `log 2`, `2π Re Q(f) = ∫ φ W` (`φ = |F|²`, `W = 2πμ + J`, `Hybrid.lean`). If `W ≥ s − t` with a
step function `t ≥ 0` (cells), then `∫ φ W ≥ 2πN s − ∫ φ t`. Instead of the cap `φ ≤ 2aN`, bound
`∫ φ t` through the Taylor polynomial of `F`:

* `F(r) = P(r) + E(r)`, `P(r) = Σ_{m<K} r^m B_m`, `B_m = ∫ f(u) (iu)^m/m!`, `|E(r)| ≤ 2(|r|a)^K/K! ∫|f|`;
* `∫ t |P|² = xᵀHx + yᵀHy` (`x = Re B`, `y = Im B`, `H_{mk}` = moments of `t`), exact;
* Bessel: `N ≥ 2(d·x + e·y) − dᵀĜd − eᵀĜe` for all real `d, e` (`0 ≤ ∫|f − ḡ|²`, `g` polynomial);
* a PSD certificate `H ≤ α(2Λ − ΛᵀĜΛ)` then gives `∫ t|P|² ≤ αN`, and `∫ φ t ≤ ((1+ε)α + (1+1/ε)β) N`.
-/
import Riemann.Deep.Yoshida.Hybrid

noncomputable section

open Complex MeasureTheory Set Filter Finset
open scoped ComplexConjugate

namespace Fayed.Yoshida

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF

/-! ### Polynomial integrals -/

lemma integral_poly_sq (γ : ℕ → ℝ) (K : ℕ) (lo hi : ℝ) :
    ∫ r in lo..hi, (∑ m ∈ range K, γ m * r ^ m) ^ 2
      = ∑ m ∈ range K, ∑ k ∈ range K,
          γ m * γ k * ((hi ^ (m + k + 1) - lo ^ (m + k + 1)) / ((m + k : ℕ) + 1 : ℝ)) := by
  have e : ∀ r : ℝ, (∑ m ∈ range K, γ m * r ^ m) ^ 2
      = ∑ m ∈ range K, ∑ k ∈ range K, γ m * γ k * r ^ (m + k) := by
    intro r
    rw [sq, Finset.sum_mul_sum]
    refine Finset.sum_congr rfl fun m _ => Finset.sum_congr rfl fun k _ => ?_
    ring
  simp_rw [e]
  rw [intervalIntegral.integral_finsetSum (fun m _ => ?_)]
  · refine Finset.sum_congr rfl fun m _ => ?_
    rw [intervalIntegral.integral_finsetSum (fun k _ => ?_)]
    · refine Finset.sum_congr rfl fun k _ => ?_
      rw [intervalIntegral.integral_const_mul, integral_pow]
    · exact (continuous_const.mul (continuous_pow _)).intervalIntegrable _ _
  · exact (continuous_finsetSum _ fun k _ => by fun_prop).intervalIntegrable _ _

/-- `∫ r^p` over `lo ≤ |r| ≤ hi` (both signs). -/
def cellMom (lo hi : ℝ) (p : ℕ) : ℝ :=
  ((hi ^ (p + 1) - lo ^ (p + 1)) + ((-lo) ^ (p + 1) - (-hi) ^ (p + 1))) / ((p : ℝ) + 1)

/-- Moments of the step function of a cell list. -/
def momSum : List (ℝ × ℝ × ℝ) → ℕ → ℝ
  | [], _ => 0
  | c :: cs, p => c.2.2 * cellMom c.1 c.2.1 p + momSum cs p

lemma cellMom_even_nonneg {lo hi : ℝ} (h0 : 0 ≤ lo) (h1 : lo ≤ hi) (K : ℕ) :
    0 ≤ cellMom lo hi (2 * K) := by
  unfold cellMom
  have e1 : (-lo) ^ (2 * K + 1) = -lo ^ (2 * K + 1) := by rw [neg_pow, pow_succ, pow_mul]; simp
  have e2 : (-hi) ^ (2 * K + 1) = -hi ^ (2 * K + 1) := by rw [neg_pow, pow_succ, pow_mul]; simp
  rw [e1, e2]
  have : lo ^ (2 * K + 1) ≤ hi ^ (2 * K + 1) := pow_le_pow_left₀ h0 h1 _
  apply div_nonneg _ (by positivity)
  linarith

/-! ### Taylor split of `F` -/

/-- Scaled moments `B_m = ∫ f(u) (iu)^m / m!`. -/
def momB (f : ℝ → ℂ) (m : ℕ) : ℂ := ∫ u : ℝ, f u * ((I * u) ^ m / (m.factorial : ℂ))

/-- Taylor polynomial `P(r) = Σ_{m<K} r^m B_m` of `F(r)`. -/
def taylorP (f : ℝ → ℂ) (K : ℕ) (r : ℝ) : ℂ := ∑ m ∈ range K, (r : ℂ) ^ m * momB f m

lemma testFn_integrable_mul {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) {g : ℝ → ℂ} (hg : Continuous g) :
    Integrable (fun u => f u * g u) :=
  (hf.continuous.mul hg).integrable_of_hasCompactSupport (hf.hasCompactSupport.mul_right)

lemma testFn_eq_zero {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) {u : ℝ}
    (hu : u ∉ Icc (-a) a) : f u = 0 := by
  by_contra h
  exact hu (hf.2 (Function.mem_support.mpr h))

lemma taylorP_eq_integral {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (K : ℕ) (r : ℝ) :
    taylorP f K r = ∫ u : ℝ, f u * ∑ m ∈ range K, (I * r * u) ^ m / (m.factorial : ℂ) := by
  unfold taylorP momB
  simp_rw [Finset.mul_sum]
  rw [integral_finsetSum _ (fun m _ => testFn_integrable_mul hf (by fun_prop))]
  refine Finset.sum_congr rfl fun m _ => ?_
  rw [← integral_const_mul]
  congr 1
  funext u
  rw [mul_pow, mul_pow]
  ring

lemma taylor_split {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (K : ℕ) {r : ℝ}
    (hr : |r| * a ≤ (K + 1) / 2) :
    ‖paperFT f r - taylorP f K r‖ ≤ 2 * (|r| * a) ^ K / K.factorial * ∫ u, ‖f u‖ := by
  rw [taylorP_eq_integral hf K r]
  unfold paperFT
  have hi1 : Integrable (fun u : ℝ => f u * Complex.exp (I * (r : ℂ) * u)) :=
    testFn_integrable_mul hf (by fun_prop)
  have hi2 : Integrable (fun u : ℝ => f u * ∑ m ∈ range K, (I * r * u) ^ m / (m.factorial : ℂ)) :=
    testFn_integrable_mul hf (by fun_prop)
  rw [← integral_sub hi1 hi2]
  have hfn : Integrable (fun u => ‖f u‖) := hf.continuous.norm.integrable_of_hasCompactSupport
    hf.hasCompactSupport.norm
  set C := 2 * (|r| * a) ^ K / K.factorial
  have hb : ∀ u, ‖f u * Complex.exp (I * (r : ℂ) * u)
      - f u * ∑ m ∈ range K, (I * r * u) ^ m / (m.factorial : ℂ)‖ ≤ C * ‖f u‖ := by
    intro u
    by_cases hu : u ∈ Icc (-a) a
    · rw [← mul_sub, norm_mul, mul_comm]
      gcongr
      have hx : ‖I * (r : ℂ) * u‖ = |r| * |u| := by
        rw [norm_mul, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Complex.norm_real,
          Real.norm_eq_abs, Real.norm_eq_abs]
      have hua : |u| ≤ a := abs_le.2 ⟨hu.1, hu.2⟩
      have hxa : ‖I * (r : ℂ) * u‖ ≤ |r| * a := by
        rw [hx]; exact mul_le_mul_of_nonneg_left hua (abs_nonneg r)
      have hcond : ‖I * (r : ℂ) * u‖ / (K.succ : ℝ) ≤ 1 / 2 := by
        rw [div_le_iff₀ (by positivity)]
        push_cast
        linarith
      have h := Complex.exp_bound' hcond
      refine h.trans ?_
      unfold C
      have : ‖I * (r : ℂ) * u‖ ^ K ≤ (|r| * a) ^ K := pow_le_pow_left₀ (norm_nonneg _) hxa K
      have hf0 : (0 : ℝ) < K.factorial := by positivity
      rw [div_mul_eq_mul_div, mul_comm 2 ((|r| * a) ^ K)]
      gcongr
    · rw [testFn_eq_zero hf hu]
      simp
  refine (norm_integral_le_of_norm_le (hfn.const_mul C) (Eventually.of_forall hb)).trans
    (le_of_eq ?_)
  rw [integral_const_mul]

/-! ### Cells -/

lemma norm_sq_taylorP (f : ℝ → ℂ) (K : ℕ) (r : ℝ) :
    ‖taylorP f K r‖ ^ 2 = (∑ m ∈ range K, (momB f m).re * r ^ m) ^ 2
      + (∑ m ∈ range K, (momB f m).im * r ^ m) ^ 2 := by
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply]
  unfold taylorP
  rw [Complex.re_sum, Complex.im_sum]
  have e1 : ∀ m, ((r : ℂ) ^ m * momB f m).re = (momB f m).re * r ^ m := by
    intro m
    rw [← Complex.ofReal_pow, Complex.re_ofReal_mul, mul_comm]
  have e2 : ∀ m, ((r : ℂ) ^ m * momB f m).im = (momB f m).im * r ^ m := by
    intro m
    rw [← Complex.ofReal_pow, Complex.im_ofReal_mul, mul_comm]
  simp_rw [e1, e2]
  ring

lemma sq_add_le_eps (z w : ℂ) {ε : ℝ} (hε : 0 < ε) :
    ‖z + w‖ ^ 2 ≤ (1 + ε) * ‖z‖ ^ 2 + (1 + 1 / ε) * ‖w‖ ^ 2 := by
  have h1 := norm_add_le z w
  have h0 : 0 ≤ ‖z + w‖ := norm_nonneg _
  have hz := norm_nonneg z
  have hw := norm_nonneg w
  have key : (‖z‖ + ‖w‖) ^ 2 ≤ (1 + ε) * ‖z‖ ^ 2 + (1 + 1 / ε) * ‖w‖ ^ 2 := by
    have e : (1 + ε) * ‖z‖ ^ 2 + (1 + 1 / ε) * ‖w‖ ^ 2 - (‖z‖ + ‖w‖) ^ 2
        = (ε * ‖z‖ - ‖w‖) ^ 2 / ε := by
      field_simp
      ring
    have : 0 ≤ (ε * ‖z‖ - ‖w‖) ^ 2 / ε := by positivity
    linarith
  nlinarith

/-- Pointwise bound on a Taylor-valid range. -/
lemma phi_le_taylor {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f) (K : ℕ) {ε : ℝ}
    (hε : 0 < ε) {r : ℝ} (hr : |r| * a ≤ (K + 1) / 2) :
    ‖paperFT f r‖ ^ 2 ≤ (1 + ε) * ((∑ m ∈ range K, (momB f m).re * r ^ m) ^ 2
        + (∑ m ∈ range K, (momB f m).im * r ^ m) ^ 2)
      + (1 + 1 / ε) * (4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2) * r ^ (2 * K) := by
  have hs := sq_add_le_eps (taylorP f K r) (paperFT f r - taylorP f K r) hε
  rw [add_sub_cancel] at hs
  rw [norm_sq_taylorP] at hs
  have hE := taylor_split hf K hr
  have hE2 : ‖paperFT f r - taylorP f K r‖ ^ 2
      ≤ 4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2 * r ^ (2 * K) := by
    have h0 := norm_nonneg (paperFT f r - taylorP f K r)
    refine (pow_le_pow_left₀ h0 hE 2).trans (le_of_eq ?_)
    have e : (|r| * a) ^ K * (|r| * a) ^ K = a ^ (2 * K) * r ^ (2 * K) := by
      rw [← pow_add, mul_pow, ← two_mul, pow_mul, pow_mul, sq_abs]
      ring
    field_simp
    rw [show ((|r| * a) ^ K) ^ 2 = (|r| * a) ^ K * (|r| * a) ^ K by ring, e]
    ring
  have hc : 0 ≤ 1 + 1 / ε := by positivity
  calc ‖paperFT f r‖ ^ 2 ≤ _ := hs
    _ ≤ _ := by
      have := mul_le_mul_of_nonneg_left hE2 hc
      linarith

/-- The majorant of `φ` on a Taylor-valid range. -/
def phiMaj (f : ℝ → ℂ) (K : ℕ) (ε Cβ : ℝ) (r : ℝ) : ℝ :=
  (1 + ε) * ((∑ m ∈ range K, (momB f m).re * r ^ m) ^ 2 + (∑ m ∈ range K, (momB f m).im * r ^ m) ^ 2)
    + (1 + 1 / ε) * Cβ * r ^ (2 * K)

lemma continuous_phiMaj (f : ℝ → ℂ) (K : ℕ) (ε Cβ : ℝ) : Continuous (phiMaj f K ε Cβ) := by
  unfold phiMaj
  fun_prop

lemma integral_phiMaj (f : ℝ → ℂ) (K : ℕ) (ε Cβ lo hi : ℝ) :
    ∫ r in lo..hi, phiMaj f K ε Cβ r
      = (1 + ε) * ∑ m ∈ range K, ∑ k ∈ range K, ((momB f m).re * (momB f k).re
          + (momB f m).im * (momB f k).im) * ((hi ^ (m + k + 1) - lo ^ (m + k + 1)) / ((m + k : ℕ) + 1 : ℝ))
        + (1 + 1 / ε) * Cβ * ((hi ^ (2 * K + 1) - lo ^ (2 * K + 1)) / ((2 * K : ℕ) + 1 : ℝ)) := by
  unfold phiMaj
  have c1 : Continuous fun r : ℝ => (∑ m ∈ range K, (momB f m).re * r ^ m) ^ 2 := by fun_prop
  have c2 : Continuous fun r : ℝ => (∑ m ∈ range K, (momB f m).im * r ^ m) ^ 2 := by fun_prop
  rw [intervalIntegral.integral_add
    (f := fun r : ℝ => (1 + ε) * ((∑ m ∈ range K, (momB f m).re * r ^ m) ^ 2
      + (∑ m ∈ range K, (momB f m).im * r ^ m) ^ 2))
    (g := fun r : ℝ => (1 + 1 / ε) * Cβ * r ^ (2 * K))
    ((continuous_const.mul (c1.add c2)).intervalIntegrable _ _)
    ((continuous_const.mul (continuous_pow _)).intervalIntegrable _ _)]
  rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_add (c1.intervalIntegrable _ _)
    (c2.intervalIntegrable _ _), integral_poly_sq, integral_poly_sq, intervalIntegral.integral_const_mul,
    integral_pow, ← Finset.sum_add_distrib]
  congr 2
  refine Finset.sum_congr rfl fun m _ => ?_
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  ring

lemma cellMom_split (lo hi : ℝ) (p : ℕ) :
    (hi ^ (p + 1) - lo ^ (p + 1)) / ((p : ℝ) + 1) + ((-lo) ^ (p + 1) - (-hi) ^ (p + 1)) / ((p : ℝ) + 1)
      = cellMom lo hi p := by
  unfold cellMom
  ring

/-- One cell: `∫ φ · cellFn ≤ τ · (bound)`. -/
lemma cell_bound {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f)
    (hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2)) (K : ℕ) {ε : ℝ} (hε : 0 < ε)
    {c : ℝ × ℝ × ℝ} (h0 : 0 ≤ c.1) (h1 : c.1 ≤ c.2.1) (hτ : 0 ≤ c.2.2)
    (hK : c.2.2 = 0 ∨ c.2.1 * a ≤ (K + 1) / 2) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * cellFn c r) ∧
    ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * cellFn c r ≤ c.2.2 * ((1 + ε) * ∑ m ∈ range K, ∑ k ∈ range K,
        ((momB f m).re * (momB f k).re + (momB f m).im * (momB f k).im) * cellMom c.1 c.2.1 (m + k)
      + (1 + 1 / ε) * (4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2)
        * cellMom c.1 c.2.1 (2 * K)) := by
  obtain ⟨lo, hi, τ⟩ := c
  simp only at h0 h1 hτ hK ⊢
  set φ := fun r : ℝ => ‖paperFT f r‖ ^ 2 with hφ
  have e : (fun r => φ r * cellFn (lo, hi, τ) r)
      = fun r => (Ico lo hi).indicator (fun r => τ * φ r) r
        + (Ioc (-hi) (-lo)).indicator (fun r => τ * φ r) r := by
    funext r
    unfold cellFn
    simp only [Set.indicator]
    split_ifs <;> ring
  have i1 : Integrable fun r => (Ico lo hi).indicator (fun r => τ * φ r) r :=
    (hφi.const_mul τ).indicator measurableSet_Ico
  have i2 : Integrable fun r => (Ioc (-hi) (-lo)).indicator (fun r => τ * φ r) r :=
    (hφi.const_mul τ).indicator measurableSet_Ioc
  refine ⟨by rw [e]; exact i1.add i2, ?_⟩
  rw [e, integral_add i1 i2, integral_indicator measurableSet_Ico,
    integral_indicator measurableSet_Ioc, integral_const_mul, integral_const_mul, ← mul_add]
  rcases hK with hK | hK
  · rw [hK]; simp
  refine mul_le_mul_of_nonneg_left ?_ hτ
  set Cβ := 4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2
  have hle : ∀ r : ℝ, |r| ≤ hi → φ r ≤ phiMaj f K ε Cβ r := by
    intro r hr
    have hr' : |r| * a ≤ (K + 1) / 2 :=
      le_trans (mul_le_mul_of_nonneg_right hr ha.le) hK
    exact phi_le_taylor ha hf K hε hr'
  have hP := continuous_phiMaj f K ε Cβ
  have j1 : ∫ r in Ico lo hi, φ r ≤ ∫ r in lo..hi, phiMaj f K ε Cβ r := by
    rw [integral_Ico_eq_integral_Ioc, ← intervalIntegral.integral_of_le h1]
    refine intervalIntegral.integral_mono_on h1 (hφi.intervalIntegrable) (hP.intervalIntegrable _ _)
      fun r hr => hle r ?_
    rw [abs_of_nonneg (by linarith [hr.1])]
    exact hr.2
  have j2 : ∫ r in Ioc (-hi) (-lo), φ r ≤ ∫ r in (-hi)..(-lo), phiMaj f K ε Cβ r := by
    rw [← intervalIntegral.integral_of_le (by linarith)]
    refine intervalIntegral.integral_mono_on (by linarith) (hφi.intervalIntegrable)
      (hP.intervalIntegrable _ _) fun r hr => hle r ?_
    rw [abs_of_nonpos (by linarith [hr.2])]
    linarith [hr.1]
  refine (add_le_add j1 j2).trans (le_of_eq ?_)
  rw [integral_phiMaj, integral_phiMaj]
  have e2 : ∀ p : ℕ, ((p : ℕ) : ℝ) + 1 = ((p : ℝ) + 1) := fun p => rfl
  rw [← cellMom_split lo hi (2 * K)]
  simp_rw [← cellMom_split lo hi]
  push_cast
  simp only [mul_add, Finset.sum_add_distrib]
  ring

/-- Bound for the whole step function. -/
lemma step_bound {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f)
    (hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2)) (K : ℕ) {ε : ℝ} (hε : 0 < ε) :
    ∀ cs : List (ℝ × ℝ × ℝ), (∀ c ∈ cs, 0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2) →
    (∀ c ∈ cs, c.2.2 = 0 ∨ c.2.1 * a ≤ (K + 1) / 2) →
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * stepFn cs r) ∧
    ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * stepFn cs r ≤ (1 + ε) * ∑ m ∈ range K, ∑ k ∈ range K,
        ((momB f m).re * (momB f k).re + (momB f m).im * (momB f k).im) * momSum cs (m + k)
      + (1 + 1 / ε) * (4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2)
        * momSum cs (2 * K)
  | [], _, _ => by
    have e : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * stepFn [] r) = fun _ => (0 : ℝ) := by
      funext r; simp [stepFn]
    rw [e]
    simp [momSum]
  | c :: cs, h, hK => by
    simp only [List.mem_cons, forall_eq_or_imp] at h hK
    obtain ⟨hc1, hc2⟩ := cell_bound ha hf hφi K hε h.1.1 h.1.2.1 h.1.2.2 hK.1
    obtain ⟨hs1, hs2⟩ := step_bound ha hf hφi K hε cs h.2 hK.2
    have e : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * stepFn (c :: cs) r)
        = fun r : ℝ => ‖paperFT f r‖ ^ 2 * cellFn c r + ‖paperFT f r‖ ^ 2 * stepFn cs r := by
      funext r
      simp only [stepFn]
      ring
    refine ⟨by rw [e]; exact hc1.add hs1, ?_⟩
    rw [e, integral_add hc1 hs1]
    refine (add_le_add hc2 hs2).trans (le_of_eq ?_)
    simp only [momSum, mul_add, Finset.sum_add_distrib, Finset.mul_sum]
    ring_nf

/-! ### Bessel inequality for the moments -/

/-- `σ_m = (−1)^{⌊m/2⌋}`. -/
def sgn2 (m : ℕ) : ℝ := (-1) ^ (m / 2)

/-- `Ĝ_{mk} = σ_m σ_k (∫_{−a}^a u^{m+k} du) / (m! k!)`. -/
def Ghat (a : ℝ) (m k : ℕ) : ℝ :=
  sgn2 m * sgn2 k * ((a ^ (m + k + 1) - (-a) ^ (m + k + 1)) / ((m + k : ℕ) + 1 : ℝ))
    / ((m.factorial : ℝ) * k.factorial)

lemma I_pow_eq (m : ℕ) : I ^ m = (sgn2 m : ℂ) * I ^ (m % 2) := by
  conv_lhs => rw [← Nat.div_add_mod m 2, pow_add, pow_mul, Complex.I_sq]
  unfold sgn2
  push_cast
  ring

lemma norm_sq_poly (κ : ℕ → ℂ) (K : ℕ) (u : ℝ) :
    ‖∑ m ∈ range K, (u : ℂ) ^ m * κ m‖ ^ 2
      = (∑ m ∈ range K, (κ m).re * u ^ m) ^ 2 + (∑ m ∈ range K, (κ m).im * u ^ m) ^ 2 := by
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply, Complex.re_sum, Complex.im_sum]
  have e1 : ∀ m, ((u : ℂ) ^ m * κ m).re = (κ m).re * u ^ m := by
    intro m
    rw [← Complex.ofReal_pow, Complex.re_ofReal_mul, mul_comm]
  have e2 : ∀ m, ((u : ℂ) ^ m * κ m).im = (κ m).im * u ^ m := by
    intro m
    rw [← Complex.ofReal_pow, Complex.im_ofReal_mul, mul_comm]
  simp_rw [e1, e2]
  ring

lemma ghat_term (a : ℝ) (d e : ℕ → ℝ) (m k : ℕ) :
    let κ := fun j : ℕ => ((d j : ℂ) - I * e j) * I ^ j / (j.factorial : ℂ)
    ((κ m).re * (κ k).re + (κ m).im * (κ k).im)
        * ((a ^ (m + k + 1) - (-a) ^ (m + k + 1)) / ((m + k : ℕ) + 1 : ℝ))
      = (d m * d k + e m * e k) * Ghat a m k := by
  intro κ
  have hfm : (m.factorial : ℝ) ≠ 0 := by positivity
  have hfk : (k.factorial : ℝ) ≠ 0 := by positivity
  have hκ : ∀ j, κ j = ((d j : ℂ) - I * e j) * (sgn2 j : ℂ) * I ^ (j % 2) / ((j.factorial : ℝ) : ℂ) := by
    intro j
    show ((d j : ℂ) - I * e j) * I ^ j / (j.factorial : ℂ) = _
    rw [I_pow_eq]
    push_cast
    ring
  have hodd : m % 2 ≠ k % 2 → (a ^ (m + k + 1) - (-a) ^ (m + k + 1)) = 0 := by
    intro h
    have hev : Even (m + k + 1) := Nat.even_iff.2 (by omega)
    rw [hev.neg_pow]
    ring
  unfold Ghat
  rcases Nat.mod_two_eq_zero_or_one m with hm | hm <;>
    rcases Nat.mod_two_eq_zero_or_one k with hk | hk
  · rw [hκ m, hκ k, hm, hk]
    simp only [pow_zero, mul_one, Complex.div_ofReal_re, Complex.div_ofReal_im, Complex.mul_re,
      Complex.mul_im, Complex.sub_re, Complex.sub_im, Complex.ofReal_re, Complex.ofReal_im,
      Complex.I_re, Complex.I_im]
    field_simp
    ring
  · rw [hodd (by omega)]
    simp
  · rw [hodd (by omega)]
    simp
  · rw [hκ m, hκ k, hm, hk]
    simp only [pow_one, Complex.div_ofReal_re, Complex.div_ofReal_im, Complex.mul_re,
      Complex.mul_im, Complex.sub_re, Complex.sub_im, Complex.ofReal_re, Complex.ofReal_im,
      Complex.I_re, Complex.I_im]
    field_simp
    ring

lemma norm_sub_conj_sq (z w : ℂ) : ‖z - conj w‖ ^ 2 = ‖z‖ ^ 2 - 2 * (z * w).re + ‖w‖ ^ 2 := by
  rw [← Complex.normSq_eq_norm_sq, ← Complex.normSq_eq_norm_sq, ← Complex.normSq_eq_norm_sq]
  simp only [Complex.normSq_apply, Complex.sub_re, Complex.sub_im, Complex.conj_re,
    Complex.conj_im, Complex.mul_re]
  ring

/-- **Bessel inequality** for the moments `B_m`: for all real `d, e`,
`2 Σ (d_m Re B_m + e_m Im B_m) − Σ (d_m d_k + e_m e_k) Ĝ_{mk} ≤ ∫ |f|²`. -/
theorem bessel {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f) (K : ℕ) (d e : ℕ → ℝ) :
    2 * ∑ m ∈ range K, (d m * (momB f m).re + e m * (momB f m).im)
      - ∑ m ∈ range K, ∑ k ∈ range K, (d m * d k + e m * e k) * Ghat a m k
      ≤ ∫ u, ‖f u‖ ^ 2 := by
  set κ := fun j : ℕ => ((d j : ℂ) - I * e j) * I ^ j / (j.factorial : ℂ) with hκ
  set g : ℝ → ℂ := fun u => ∑ m ∈ range K, (u : ℂ) ^ m * κ m with hg
  have hgc : Continuous g := by rw [hg]; fun_prop
  have hgw : ∀ u : ℝ, g u = ∑ m ∈ range K, ((d m : ℂ) - I * e m) * ((I * u) ^ m / (m.factorial : ℂ)) := by
    intro u
    refine Finset.sum_congr rfl fun m _ => ?_
    simp only [hκ]
    rw [mul_pow]
    ring
  -- (i) the linear term
  have hlin : (∫ u, f u * g u).re = ∑ m ∈ range K, (d m * (momB f m).re + e m * (momB f m).im) := by
    have e1 : (fun u => f u * g u)
        = fun u => ∑ m ∈ range K, ((d m : ℂ) - I * e m) * (f u * ((I * u) ^ m / (m.factorial : ℂ))) := by
      funext u
      rw [hgw, Finset.mul_sum]
      refine Finset.sum_congr rfl fun m _ => ?_
      ring
    rw [e1, integral_finsetSum _ (fun m _ => (testFn_integrable_mul hf (by fun_prop)).const_mul _),
      Complex.re_sum]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [integral_const_mul]
    unfold momB
    simp [Complex.mul_re]
  -- (ii) the quadratic term
  have hquad : ∫ u in (-a)..a, ‖g u‖ ^ 2
      = ∑ m ∈ range K, ∑ k ∈ range K, (d m * d k + e m * e k) * Ghat a m k := by
    have e1 : (fun u => ‖g u‖ ^ 2) = fun u => (∑ m ∈ range K, (κ m).re * u ^ m) ^ 2
        + (∑ m ∈ range K, (κ m).im * u ^ m) ^ 2 := by
      funext u
      exact norm_sq_poly κ K u
    rw [e1, intervalIntegral.integral_add (Continuous.intervalIntegrable (by fun_prop) _ _)
      (Continuous.intervalIntegrable (by fun_prop) _ _), integral_poly_sq,
      integral_poly_sq, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [← ghat_term a d e m k]
    ring
  -- (iii) positivity of ∫ |f − ḡ|²
  have hzero : ∀ u, u ∉ Icc (-a) a → f u = 0 := fun u hu => testFn_eq_zero hf hu
  have hfc := hf.continuous
  have hpt : ∀ u, ‖f u - conj (g u)‖ ^ 2 = ‖f u‖ ^ 2 - 2 * (f u * g u).re + ‖g u‖ ^ 2 :=
    fun u => norm_sub_conj_sq _ _
  have hI : IntegrableOn (fun u => ‖f u‖ ^ 2 - 2 * (f u * g u).re + ‖g u‖ ^ 2) (Icc (-a) a) :=
    (Continuous.integrableOn_Icc (by fun_prop))
  have hpos : 0 ≤ ∫ u in Icc (-a) a, (‖f u‖ ^ 2 - 2 * (f u * g u).re + ‖g u‖ ^ 2) := by
    refine setIntegral_nonneg measurableSet_Icc fun u _ => ?_
    rw [← hpt]
    positivity
  have i1 : IntegrableOn (fun u => ‖f u‖ ^ 2) (Icc (-a) a) := Continuous.integrableOn_Icc (by fun_prop)
  have i2 : IntegrableOn (fun u => 2 * (f u * g u).re) (Icc (-a) a) :=
    Continuous.integrableOn_Icc (by fun_prop)
  have i3 : IntegrableOn (fun u => ‖g u‖ ^ 2) (Icc (-a) a) := Continuous.integrableOn_Icc (by fun_prop)
  have i12 : IntegrableOn (fun u => ‖f u‖ ^ 2 - 2 * (f u * g u).re) (Icc (-a) a) := i1.sub i2
  rw [integral_add i12 i3, integral_sub i1 i2, integral_const_mul] at hpos
  have hN : ∫ u in Icc (-a) a, ‖f u‖ ^ 2 = ∫ u, ‖f u‖ ^ 2 :=
    setIntegral_eq_integral_of_forall_compl_eq_zero (fun u hu => by simp [hzero u hu])
  have hR : ∫ u in Icc (-a) a, (f u * g u).re = (∫ u, f u * g u).re := by
    have hfg : IntegrableOn (fun u => f u * g u) (Icc (-a) a) := Continuous.integrableOn_Icc (by fun_prop)
    have h := integral_re hfg
    simp only [RCLike.re_to_complex] at h
    rw [h]
    congr 1
    exact setIntegral_eq_integral_of_forall_compl_eq_zero (fun u hu => by simp [hzero u hu])
  have hG : ∫ u in Icc (-a) a, ‖g u‖ ^ 2 = ∫ u in (-a)..a, ‖g u‖ ^ 2 := by
    rw [intervalIntegral.integral_of_le (by linarith), integral_Icc_eq_integral_Ioc]
  rw [hN, hR, hG, hlin, hquad] at hpos
  linarith

lemma momSum_even_nonneg : ∀ cs : List (ℝ × ℝ × ℝ),
    (∀ c ∈ cs, 0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2) → ∀ K : ℕ, 0 ≤ momSum cs (2 * K)
  | [], _, _ => by simp [momSum]
  | c :: cs, h, K => by
    simp only [List.mem_cons, forall_eq_or_imp] at h
    simp only [momSum]
    have h1 := cellMom_even_nonneg h.1.1 h.1.2.1 K
    have h2 := momSum_even_nonneg cs h.2 K
    have h3 := mul_nonneg h.1.2.2 h1
    linarith

/-! ### The certificate theorem -/

/-- **Moment certificate ⇒ Weil positivity** on the window `2a`. -/
theorem weilPos_of_moment {a s R ε α β : ℝ} (ha : 0 < a) (hlog : 2 * a < Real.log 2)
    (cs : List (ℝ × ℝ × ℝ)) (hchain : Chain 0 cs R) (hok : ∀ c ∈ cs, CellOK a s c)
    (htail : ∀ r : ℝ, R ≤ |r| → s ≤ W a r) (K : ℕ)
    (hK : ∀ c ∈ cs, c.2.2 = 0 ∨ c.2.1 * a ≤ (K + 1) / 2) (Λ : ℕ → ℕ → ℝ)
    (hpsd : ∀ x : ℕ → ℝ, ∑ m ∈ range K, ∑ k ∈ range K, x m * x k * momSum cs (m + k)
      ≤ α * (2 * ∑ m ∈ range K, (∑ i ∈ range K, Λ m i * x i) * x m
        - ∑ m ∈ range K, ∑ k ∈ range K,
            ((∑ i ∈ range K, Λ m i * x i) * (∑ i ∈ range K, Λ k i * x i)) * Ghat a m k))
    (hε : 0 < ε) (hα : 0 ≤ α)
    (hβ : 8 * a ^ (2 * K + 1) * momSum cs (2 * K) / (K.factorial : ℝ) ^ 2 ≤ β)
    (hnum : (1 + ε) * α + (1 + 1 / ε) * β ≤ 2 * Real.pi * s) :
    WeilPos (2 * a) := by
  intro f hf
  have hf' : TestFn a f := by
    have e : 2 * a / 2 = a := by ring
    rwa [e] at hf
  obtain ⟨hφw, hre⟩ := weilQ_re_eq_W ha.le hlog hf'
  have hk := hf'.autocorr_contDiff
  have hkc := hf'.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  set N := ∫ u, ‖f u‖ ^ 2 with hN
  have hN0 : 0 ≤ N := integral_nonneg fun u => by positivity
  have hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2) := by
    refine (integrable_paperFT_ofReal hFk).re.congr (Eventually.of_forall fun r => ?_)
    simp only [paperFT_autocorr_real hf', RCLike.re_to_complex, Complex.ofReal_re]
  have hτ : ∀ c ∈ cs, 0 ≤ c.2.2 := fun c hc => (hok c hc).2.2.1
  have hw : ∀ r : ℝ, s - stepFn cs r ≤ W a r := by
    intro r
    rcases le_or_gt R |r| with hR | hR
    · have := stepFn_nonneg cs hτ r
      linarith [htail r hR]
    · obtain ⟨c, hc, h1, h2⟩ := chain_cover cs 0 R hchain |r| (abs_nonneg r) hR
      have := le_stepFn cs hτ c hc r h1 h2
      linarith [(hok c hc).2.2.2 r h1 h2]
  obtain ⟨hsi, hsb⟩ := step_bound ha hf' hφi K hε cs
    (fun c hc => ⟨(hok c hc).1, (hok c hc).2.1, (hok c hc).2.2.1⟩) hK
  -- lower bound by the weight
  have hlow : s * (2 * Real.pi * N) - ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * stepFn cs r
      ≤ ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W a r := by
    have hint : Integrable (fun r : ℝ => s * ‖paperFT f r‖ ^ 2 - ‖paperFT f r‖ ^ 2 * stepFn cs r) :=
      (hφi.const_mul s).sub hsi
    have h := integral_mono hint hφw (fun r => by
      have h1 := mul_le_mul_of_nonneg_left (hw r) (by positivity : (0 : ℝ) ≤ ‖paperFT f r‖ ^ 2)
      simp only
      nlinarith)
    rw [integral_sub (hφi.const_mul s) hsi, integral_const_mul, plancherel hf'] at h
    exact h
  -- the moment part
  set x : ℕ → ℝ := fun m => (momB f m).re
  set y : ℕ → ℝ := fun m => (momB f m).im
  have hX : ∑ m ∈ range K, ∑ k ∈ range K, ((momB f m).re * (momB f k).re
      + (momB f m).im * (momB f k).im) * momSum cs (m + k) ≤ α * N := by
    have hx := hpsd x
    have hy := hpsd y
    have hb := bessel ha hf' K (fun m => ∑ i ∈ range K, Λ m i * x i) (fun m => ∑ i ∈ range K, Λ m i * y i)
    have e1 : ∑ m ∈ range K, ∑ k ∈ range K, ((momB f m).re * (momB f k).re
        + (momB f m).im * (momB f k).im) * momSum cs (m + k)
        = ∑ m ∈ range K, ∑ k ∈ range K, x m * x k * momSum cs (m + k)
          + ∑ m ∈ range K, ∑ k ∈ range K, y m * y k * momSum cs (m + k) := by
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun m _ => ?_
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun k _ => ?_
      ring
    have e2 : 2 * ∑ m ∈ range K, ((∑ i ∈ range K, Λ m i * x i) * (momB f m).re
        + (∑ i ∈ range K, Λ m i * y i) * (momB f m).im)
        - ∑ m ∈ range K, ∑ k ∈ range K, ((∑ i ∈ range K, Λ m i * x i) * (∑ i ∈ range K, Λ k i * x i)
          + (∑ i ∈ range K, Λ m i * y i) * (∑ i ∈ range K, Λ k i * y i)) * Ghat a m k
        = (2 * ∑ m ∈ range K, (∑ i ∈ range K, Λ m i * x i) * x m
          - ∑ m ∈ range K, ∑ k ∈ range K,
            ((∑ i ∈ range K, Λ m i * x i) * (∑ i ∈ range K, Λ k i * x i)) * Ghat a m k)
        + (2 * ∑ m ∈ range K, (∑ i ∈ range K, Λ m i * y i) * y m
          - ∑ m ∈ range K, ∑ k ∈ range K,
            ((∑ i ∈ range K, Λ m i * y i) * (∑ i ∈ range K, Λ k i * y i)) * Ghat a m k) := by
      simp only [x, y, add_mul, Finset.sum_add_distrib, mul_add]
      ring
    rw [e1]
    rw [e2] at hb
    nlinarith
  have hA : (∫ u, ‖f u‖) ^ 2 ≤ 2 * a * N := sq_integral_norm_le ha hf'
  have hM0 : 0 ≤ momSum cs (2 * K) :=
    momSum_even_nonneg cs (fun c hc => ⟨(hok c hc).1, (hok c hc).2.1, (hok c hc).2.2.1⟩) K
  set A := ∫ u, ‖f u‖
  have hK0 : (0 : ℝ) < (K.factorial : ℝ) ^ 2 := by positivity
  have hβN : (1 + 1 / ε) * (4 * a ^ (2 * K) * A ^ 2 / (K.factorial : ℝ) ^ 2) * momSum cs (2 * K)
      ≤ (1 + 1 / ε) * β * N := by
    have hc : 0 ≤ 1 + 1 / ε := by positivity
    rw [mul_assoc, mul_assoc (1 + 1 / ε)]
    refine mul_le_mul_of_nonneg_left ?_ hc
    have h1 : 4 * a ^ (2 * K) * A ^ 2 / (K.factorial : ℝ) ^ 2 * momSum cs (2 * K)
        ≤ 4 * a ^ (2 * K) * (2 * a * N) / (K.factorial : ℝ) ^ 2 * momSum cs (2 * K) := by
      gcongr
    have h2 : 4 * a ^ (2 * K) * (2 * a * N) / (K.factorial : ℝ) ^ 2 * momSum cs (2 * K)
        = (8 * a ^ (2 * K + 1) * momSum cs (2 * K) / (K.factorial : ℝ) ^ 2) * N := by
      ring
    rw [h2] at h1
    exact h1.trans (mul_le_mul_of_nonneg_right hβ hN0)
  have hstep : ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * stepFn cs r ≤ (1 + ε) * α * N + (1 + 1 / ε) * β * N := by
    have h1 := mul_le_mul_of_nonneg_left hX (by linarith : (0 : ℝ) ≤ 1 + ε)
    linarith
  have hfin : 0 ≤ ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W a r := by
    have h3 : ((1 + ε) * α + (1 + 1 / ε) * β) * N ≤ (2 * Real.pi * s) * N :=
      mul_le_mul_of_nonneg_right hnum hN0
    nlinarith
  rw [hre]
  have hπ := Real.pi_pos
  positivity

end Fayed.Yoshida

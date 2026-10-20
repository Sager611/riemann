import Riemann.Deep.Inputs.PerronKernel
import Mathlib.NumberTheory.Chebyshev
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.NumberTheory.Harmonic.ZetaAsymp
import Mathlib.NumberTheory.ZetaValues

/-!
# Power-saving PNT, part 1: Perron's formula for ψ and the contour shift

* `perron_psi`: for `x ≥ 1` at distance `≥ 1/2` from every integer and `T > 0`,
  `‖ψ(x) - (1/2π) ∫_{-T}^{T} F(4+iy) dy‖ ≤ 4x⁴(x+1)/(πT) · π²/6`, where
  `F(s) = -ζ'(s)/ζ(s) · x^s/s` (truncated Perron on the line Re s = 4; the kernel bound is
  `perron_kernel` from Inputs/PerronKernel.lean).
* `shift_psi`: if ζ has no zero in Re s > σ₀ and σ₀ < σ₁ < 1, the rectangle
  [σ₁, 4] × [-T, T] contains only the pole s = 1 of F (residue x).
-/

open Complex Set MeasureTheory intervalIntegral Filter Topology Asymptotics
open scoped Real Interval LSeries.notation ArithmeticFunction.vonMangoldt

namespace Fayed.PNTError

open Fayed.InputsProof.PerronAux

noncomputable section

/-- ζ has no zero in the open half-plane Re s > σ₀. -/
def ZeroFree (σ₀ : ℝ) : Prop := ∀ s : ℂ, σ₀ < s.re → riemannZeta s ≠ 0

/-- The Perron integrand for ψ: -ζ'(s)/ζ(s) · x^s / s. -/
def Fψ (x : ℝ) (s : ℂ) : ℂ := (-deriv riemannZeta s / riemannZeta s) * (x : ℂ) ^ s / s

/-- The Perron frequency log x - log n. -/
def Lx (x : ℝ) (n : ℕ) : ℝ := Real.log x - Real.log n

lemma term_mul_eq {x : ℝ} (hx : 0 < x) (s : ℂ) (n : ℕ) :
    LSeries.term ↗Λ s n * (x : ℂ) ^ s / s = (Λ n : ℂ) * kern (Lx x n) s := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp
  · rw [LSeries.term_of_ne_zero hn, kern, Lx]
    have hn' : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn
    have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
    rw [cpow_def_of_ne_zero hn', cpow_def_of_ne_zero hx', ← Complex.natCast_log,
      ← Complex.ofReal_log hx.le]
    push_cast
    rw [sub_mul, Complex.exp_sub]
    ring

lemma Fψ_eq_tsum {x : ℝ} (hx : 0 < x) {s : ℂ} (hs : 1 < s.re) :
    Fψ x s = ∑' n : ℕ, (Λ n : ℂ) * kern (Lx x n) s := by
  rw [Fψ, ← ArithmeticFunction.LSeries_vonMangoldt_eq_deriv_riemannZeta_div hs, LSeries,
    ← tsum_mul_right, ← tsum_div_const]
  exact tsum_congr (fun n => term_mul_eq hx s n)

lemma kern_cont' (ℓ c : ℝ) (hc : c ≠ 0) : Continuous (fun y : ℝ => kern ℓ ((c : ℂ) + y * I)) := by
  unfold kern
  refine Continuous.div (by fun_prop) (by fun_prop) (fun y => ?_)
  intro h
  have := congrArg Complex.re h
  simp at this
  exact hc this

lemma norm_kern_four (ℓ y : ℝ) : ‖kern ℓ ((4 : ℝ) + y * I)‖ ≤ Real.exp (4 * ℓ) / 4 := by
  have h := @norm_kern_le_re ℓ 4 y (by norm_num)
  rwa [abs_of_pos (by norm_num : (0 : ℝ) < 4), mul_comm ℓ 4] at h

lemma exp_four_Lx {x : ℝ} (hx : 0 < x) {n : ℕ} (hn : n ≠ 0) :
    Real.exp (4 * Lx x n) = x ^ 4 / (n : ℝ) ^ 4 := by
  have hn' : (0 : ℝ) < n := Nat.cast_pos.mpr (Nat.pos_of_ne_zero hn)
  rw [Lx, mul_sub, Real.exp_sub]
  have e : ∀ z : ℝ, 0 < z → Real.exp (4 * Real.log z) = z ^ 4 := by
    intro z hz
    rw [show (4 : ℝ) * Real.log z = ((4 : ℕ) : ℝ) * Real.log z by norm_num, Real.exp_nat_mul,
      Real.exp_log hz]
  rw [e x hx, e n hn']

lemma vonMangoldt_le_self (n : ℕ) : Λ n ≤ n :=
  ArithmeticFunction.vonMangoldt_le_log.trans (Real.log_le_self (Nat.cast_nonneg n))

/-- Truncated Perron on Re s = 4 as a convergent sum of Perron kernels. -/
lemma hasSum_Pψ {x T : ℝ} (hx : 0 < x) (hT : 0 < T) :
    HasSum (fun n : ℕ => (Λ n : ℂ) *
        (((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, kern (Lx x n) ((4 : ℝ) + y * I)))
      (((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fψ x ((4 : ℝ) + y * I)) := by
  set c : ℂ := ((1 / (2 * π) : ℝ) : ℂ)
  set F : ℕ → ℝ → ℂ := fun n y => (Λ n : ℂ) * kern (Lx x n) ((4 : ℝ) + y * I) with hF
  have hint : ∀ n, Integrable (F n) (volume.restrict (Ioc (-T) T)) := fun n =>
    ((continuous_const.mul (kern_cont' _ 4 (by norm_num))).integrableOn_Icc).mono_set
      Ioc_subset_Icc_self
  have hbd : ∀ n, ∫ y in Ioc (-T) T, ‖F n y‖ ≤ 2 * T * x ^ 4 * (1 / (n : ℝ) ^ 3) := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · simp [F]
    · have hn0 : (0 : ℝ) < n := Nat.cast_pos.mpr (Nat.pos_of_ne_zero hn)
      have hn0' : (n : ℝ) ≠ 0 := hn0.ne'
      have h1 : ∀ y, ‖F n y‖ ≤ x ^ 4 * (1 / (n : ℝ) ^ 3) := by
        intro y
        simp only [F]
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs,
          abs_of_nonneg ArithmeticFunction.vonMangoldt_nonneg]
        have hk := norm_kern_four (Lx x n) y
        rw [exp_four_Lx hx hn] at hk
        calc Λ n * ‖kern (Lx x n) ((4 : ℝ) + y * I)‖ ≤ n * (x ^ 4 / (n : ℝ) ^ 4 / 4) :=
              mul_le_mul (vonMangoldt_le_self n) hk (norm_nonneg _) (Nat.cast_nonneg n)
          _ = x ^ 4 * (1 / (n : ℝ) ^ 3) / 4 := by field_simp
          _ ≤ x ^ 4 * (1 / (n : ℝ) ^ 3) := by
              have : 0 ≤ x ^ 4 * (1 / (n : ℝ) ^ 3) := by positivity
              linarith
      have hBint : IntegrableOn (fun _ : ℝ => x ^ 4 * (1 / (n : ℝ) ^ 3)) (Ioc (-T) T) :=
        (continuous_const.integrableOn_Icc).mono_set Ioc_subset_Icc_self
      calc ∫ y in Ioc (-T) T, ‖F n y‖ ≤ ∫ y in Ioc (-T) T, x ^ 4 * (1 / (n : ℝ) ^ 3) :=
            setIntegral_mono_on (hint n).norm hBint measurableSet_Ioc (fun y _ => h1 y)
        _ = (T - -T) * (x ^ 4 * (1 / (n : ℝ) ^ 3)) := by
            rw [setIntegral_const, Real.volume_real_Ioc_of_le (by linarith), smul_eq_mul]
        _ = 2 * T * x ^ 4 * (1 / (n : ℝ) ^ 3) := by ring
  have hsum : Summable (fun n => ∫ y in Ioc (-T) T, ‖F n y‖) :=
    Summable.of_nonneg_of_le (fun n => integral_nonneg (fun y => norm_nonneg _)) hbd
      ((Real.summable_one_div_nat_pow.mpr (by norm_num : 1 < 3)).mul_left (2 * T * x ^ 4))
  have hswap := integral_tsum_of_summable_integral_norm hint hsum
  have hsumF : Summable (fun n => ∫ y in Ioc (-T) T, F n y) :=
    Summable.of_norm_bounded hsum (fun n => norm_integral_le_integral_norm _)
  have hf : ∀ n, (Λ n : ℂ) * (c * ∫ y in (-T)..T, kern (Lx x n) ((4 : ℝ) + y * I)) =
      c * ∫ y in Ioc (-T) T, F n y := by
    intro n
    rw [intervalIntegral.integral_of_le (by linarith)]
    simp only [F]
    rw [MeasureTheory.integral_const_mul]
    ring
  have hFz : ∫ y in (-T)..T, Fψ x ((4 : ℝ) + y * I) = ∫ y in Ioc (-T) T, ∑' n, F n y := by
    rw [intervalIntegral.integral_of_le (by linarith)]
    congr 1
    funext y
    exact Fψ_eq_tsum hx (by norm_num)
  rw [hFz, ← hswap]
  simp_rw [hf]
  exact hsumF.hasSum.mul_left c

/-- |log x - log y| ≥ γ/(x+y) when |x - y| ≥ γ. -/
lemma log_sub_ge' {x y γ : ℝ} (hx : 0 < x) (hy : 0 < y) (hγ : 0 < γ) (h : γ ≤ |x - y|) :
    γ / (x + y) ≤ |Real.log x - Real.log y| := by
  rcases le_total y x with hxy | hxy
  · have hd : γ ≤ x - y := by rwa [abs_of_nonneg (by linarith)] at h
    have h1 := Real.log_le_sub_one_of_pos (div_pos hy hx)
    rw [Real.log_div hy.ne' hx.ne'] at h1
    have h2 : x * (Real.log y - Real.log x) ≤ y - x := by
      have e : x * (y / x - 1) = y - x := by field_simp
      calc x * (Real.log y - Real.log x) ≤ x * (y / x - 1) := mul_le_mul_of_nonneg_left h1 hx.le
        _ = y - x := e
    have h4 : 0 ≤ Real.log x - Real.log y := by linarith [Real.log_le_log hy hxy]
    rw [abs_of_nonneg h4, div_le_iff₀ (by linarith)]
    nlinarith
  · have hd : γ ≤ y - x := by rwa [abs_of_nonpos (by linarith), neg_sub] at h
    have h1 := Real.log_le_sub_one_of_pos (div_pos hx hy)
    rw [Real.log_div hx.ne' hy.ne'] at h1
    have h2 : y * (Real.log x - Real.log y) ≤ x - y := by
      have e : y * (x / y - 1) = x - y := by field_simp
      calc y * (Real.log x - Real.log y) ≤ y * (x / y - 1) := mul_le_mul_of_nonneg_left h1 hy.le
        _ = x - y := e
    have h4 : Real.log x - Real.log y ≤ 0 := by linarith [Real.log_le_log hx hxy]
    rw [abs_of_nonpos h4, div_le_iff₀ (by linarith)]
    nlinarith

/-- The per-term numerics of the Perron error. -/
lemma key_ineq {A ℓ x T : ℝ} {n : ℕ} (hAn : A ≤ n) (hn : (1 : ℝ) ≤ n)
    (hx : 1 ≤ x) (hT : 0 < T) (hL : 1 / 2 / (x + n) ≤ |ℓ|) :
    A * (2 * (x ^ 4 / (n : ℝ) ^ 4) / (π * |ℓ| * T)) ≤
      4 * x ^ 4 * (x + 1) / (π * T) * (1 / (n : ℝ) ^ 2) := by
  have hπ : 0 < π := Real.pi_pos
  have hn0 : (0 : ℝ) < n := by linarith
  have hn0' : (n : ℝ) ≠ 0 := hn0.ne'
  have hx0 : 0 < x := by linarith
  have hLpos : 0 < |ℓ| := lt_of_lt_of_le (by positivity) hL
  have hL0 : |ℓ| ≠ 0 := hLpos.ne'
  have hT0 : T ≠ 0 := hT.ne'
  have h1L : 1 / |ℓ| ≤ 2 * (x + n) := by
    rw [div_le_iff₀ hLpos]
    rw [div_div, div_le_iff₀ (by positivity)] at hL
    linarith
  have hE0 : 0 ≤ x ^ 4 / (n : ℝ) ^ 4 := by positivity
  have hq : 0 ≤ 2 / (π * T) := div_nonneg (by norm_num) (mul_pos hπ hT).le
  have hb0 : 0 ≤ (n : ℝ) * (x ^ 4 / (n : ℝ) ^ 4) * (2 / (π * T)) :=
    mul_nonneg (mul_nonneg hn0.le hE0) hq
  calc A * (2 * (x ^ 4 / (n : ℝ) ^ 4) / (π * |ℓ| * T))
        = A * (x ^ 4 / (n : ℝ) ^ 4) * (2 / (π * T)) * (1 / |ℓ|) := by
          field_simp
    _ ≤ (n : ℝ) * (x ^ 4 / (n : ℝ) ^ 4) * (2 / (π * T)) * (2 * (x + n)) := by
        apply mul_le_mul _ h1L (by positivity) hb0
        exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hAn hE0) hq
    _ = 4 * x ^ 4 * (x + n) / (π * T * (n : ℝ) ^ 3) := by
        field_simp
        ring
    _ ≤ 4 * x ^ 4 * ((x + 1) * n) / (π * T * (n : ℝ) ^ 3) := by
        apply div_le_div_of_nonneg_right _ (mul_pos (mul_pos hπ hT) (pow_pos hn0 3)).le
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        nlinarith
    _ = 4 * x ^ 4 * (x + 1) / (π * T) * (1 / (n : ℝ) ^ 2) := by
        field_simp

/-- **Truncated Perron formula for ψ** (line Re s = 4, height T). -/
theorem perron_psi {x T : ℝ} (hx : 1 ≤ x) (hT : 0 < T) (hfar : ∀ n : ℕ, 1 / 2 ≤ |x - n|) :
    ‖(Chebyshev.psi x : ℂ) - ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fψ x ((4 : ℝ) + y * I)‖ ≤
      4 * x ^ 4 * (x + 1) / (π * T) * (π ^ 2 / 6) := by
  have hx0 : 0 < x := by linarith
  have hS : HasSum (fun n : ℕ => (Λ n : ℂ) * (if 0 < Lx x n then 1 else 0))
      (Chebyshev.psi x : ℂ) := by
    have hfin : ∀ n ∉ Finset.Ioc 0 ⌊x⌋₊,
        (Λ n : ℂ) * (if 0 < Lx x n then 1 else 0) = 0 := by
      intro n hn
      rcases eq_or_ne n 0 with rfl | hn0
      · simp
      · have h1 : ⌊x⌋₊ < n := by
          simp only [Finset.mem_Ioc, not_and, not_le] at hn
          exact hn (Nat.pos_of_ne_zero hn0)
        have hxn : x < n := (Nat.floor_lt hx0.le).mp h1
        have h2 : Lx x n < 0 := by
          simp only [Lx]
          linarith [Real.log_lt_log hx0 hxn]
        simp [not_lt.mpr h2.le]
    have hh : HasSum (fun n : ℕ => (Λ n : ℂ) * (if 0 < Lx x n then 1 else 0))
        (∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, (Λ n : ℂ) * (if 0 < Lx x n then 1 else 0)) :=
      hasSum_sum_of_ne_finset_zero hfin
    have e : (∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, (Λ n : ℂ) * (if 0 < Lx x n then 1 else 0)) =
        (Chebyshev.psi x : ℂ) := by
      rw [Chebyshev.psi]
      push_cast
      refine Finset.sum_congr rfl (fun n hn => ?_)
      have hn' := Finset.mem_Ioc.mp hn
      have hnx : (n : ℝ) < x := by
        have h1 : (n : ℝ) ≤ x := (Nat.le_floor_iff hx0.le).mp hn'.2
        rcases h1.lt_or_eq with h | h
        · exact h
        · have := hfar n
          rw [h, sub_self, abs_zero] at this
          norm_num at this
      have h3 : 0 < Lx x n := by
        simp only [Lx]
        have hn0 : (0 : ℝ) < n := Nat.cast_pos.mpr hn'.1
        linarith [Real.log_lt_log hn0 hnx]
      simp [h3]
    rw [e] at hh
    exact hh
  have hd := hS.sub (hasSum_Pψ hx0 hT)
  have hg : HasSum (fun n : ℕ => 4 * x ^ 4 * (x + 1) / (π * T) * (1 / (n : ℝ) ^ 2))
      (4 * x ^ 4 * (x + 1) / (π * T) * (π ^ 2 / 6)) := hasSum_zeta_two.mul_left _
  refine hd.norm_le_of_bounded hg (fun n => ?_)
  rcases eq_or_ne n 0 with rfl | hn
  · simp
  · have hn1 : (1 : ℝ) ≤ n := Nat.one_le_cast.mpr (Nat.one_le_iff_ne_zero.mpr hn)
    have hn0 : (0 : ℝ) < n := by linarith
    have hLabs : 1 / 2 / (x + n) ≤ |Lx x n| := log_sub_ge' hx0 hn0 (by norm_num) (hfar n)
    have hLpos : 0 < |Lx x n| := lt_of_lt_of_le (by positivity) hLabs
    have hk := perron_kernel (abs_pos.mp hLpos) (by norm_num : (0 : ℝ) < 4) hT
    rw [← mul_sub, norm_mul, Complex.norm_real, Real.norm_eq_abs,
      abs_of_nonneg ArithmeticFunction.vonMangoldt_nonneg, norm_sub_rev]
    rw [exp_four_Lx hx0 hn] at hk
    exact (mul_le_mul_of_nonneg_left hk ArithmeticFunction.vonMangoldt_nonneg).trans
      (key_ineq (vonMangoldt_le_self n) hn1 hx hT hLabs)

/-! ### The contour shift -/

lemma analyticAt_zeta {w : ℂ} (hw : w ≠ 1) : AnalyticAt ℂ riemannZeta w :=
  DifferentiableOn.analyticAt (s := {1}ᶜ)
    (fun _ hv => (differentiableAt_riemannZeta hv).differentiableWithinAt)
    (isOpen_compl_singleton.mem_nhds hw)

lemma differentiableAt_deriv_zeta {w : ℂ} (hw : w ≠ 1) :
    DifferentiableAt ℂ (deriv riemannZeta) w :=
  (analyticAt_zeta hw).deriv.differentiableAt

lemma Fψ_differentiableAt {σ₀ : ℝ} (hZ : ZeroFree σ₀) {x : ℝ} (hx : 0 < x) {s : ℂ}
    (hs : σ₀ < s.re) (hs0 : s ≠ 0) (hs1 : s ≠ 1) : DifferentiableAt ℂ (Fψ x) s := by
  have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
  have h1 : DifferentiableAt ℂ (fun s => -deriv riemannZeta s / riemannZeta s) s :=
    ((differentiableAt_deriv_zeta hs1).neg).div (differentiableAt_riemannZeta hs1) (hZ s hs)
  have h2 : DifferentiableAt ℂ (fun s : ℂ => (x : ℂ) ^ s) s :=
    DifferentiableAt.const_cpow differentiableAt_id (Or.inl hx')
  exact (h1.mul h2).div differentiableAt_id hs0

/-- Near s = 1, F(s) - x/(s-1) is bounded. -/
lemma Fψ_near_one {x : ℝ} (hx : 0 < x) :
    (Fψ x - fun s => (x : ℂ) / (s - 1)) =O[𝓝[≠] 1] (1 : ℂ → ℂ) := by
  have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
  set g : ℂ → ℂ := fun s => (x : ℂ) ^ s / s with hg
  have hgd : DifferentiableAt ℂ g 1 :=
    (DifferentiableAt.const_cpow differentiableAt_id (Or.inl hx')).div differentiableAt_id
      one_ne_zero
  have hg1 : g 1 = x := by simp [hg]
  have hslope : Tendsto (slope g 1) (𝓝[≠] 1) (𝓝 (deriv g 1)) :=
    hasDerivAt_iff_tendsto_slope.mp hgd.hasDerivAt
  have h1 : (fun s => (g s - x) * (s - 1)⁻¹) =O[𝓝[≠] 1] (fun _ => (1 : ℂ)) := by
    refine (hslope.isBigO_one ℂ).congr' (Eventually.of_forall fun s => ?_) EventuallyEq.rfl
    rw [slope_def_field, hg1, div_eq_mul_inv]
  have hgb : g =O[𝓝[≠] 1] (fun _ => (1 : ℂ)) :=
    (hgd.continuousAt.tendsto.mono_left nhdsWithin_le_nhds).isBigO_one ℂ
  have h2 : (fun s => -(deriv riemannZeta s / riemannZeta s + (s - 1)⁻¹) * g s)
      =O[𝓝[≠] 1] (fun _ => (1 : ℂ)) := by
    have := log_deriv_riemannZeta_add_inv_sub_bounded.neg_left.mul hgb
    simpa only [mul_one] using this
  refine (h2.add h1).congr' (Eventually.of_forall fun s => ?_) (Eventually.of_forall fun s => rfl)
  simp only [Pi.sub_apply, Fψ, hg]
  ring

/-- **Contour shift**: the rectangle [σ₁, 4] × [-T, T] contains only the pole s = 1. -/
theorem shift_psi {σ₀ σ₁ x T : ℝ} (hZ : ZeroFree σ₀) (h01 : σ₀ < σ₁) (hσ₁0 : 0 < σ₁)
    (hσ₁ : σ₁ < 1) (hx : 0 < x) (hT : 0 < T) :
    ((1 / (2 * π) : ℝ) : ℂ) * (∫ y in (-T)..T, Fψ x ((4 : ℝ) + y * I)) -
      ((1 / (2 * π) : ℝ) : ℂ) * (∫ y in (-T)..T, Fψ x ((σ₁ : ℂ) + y * I)) - (x : ℂ) =
    ((1 / (2 * π) : ℝ) : ℂ) * (I * (∫ t in σ₁..4, Fψ x (t + (-T : ℝ) * I)) -
      I * (∫ t in σ₁..4, Fψ x (t + (T : ℝ) * I))) := by
  have hres : RectangleIntegral' (Fψ x) (⟨σ₁, -T⟩ : ℂ) ⟨4, T⟩ = x := by
    refine ResidueTheoremOnRectangleWithSimplePole' (p := 1) (A := x)
      (by simp only; linarith) (by simp only; linarith) ?_ ?_ (Fψ_near_one hx)
    · rw [rectangle_mem_nhds_iff]
      simp only [mem_reProdIm, one_re, one_im]
      rw [uIoo_of_le (by linarith), uIoo_of_le (by linarith)]
      exact ⟨⟨hσ₁, by norm_num⟩, ⟨by linarith, hT⟩⟩
    · intro z hz
      have hzre : z.re ∈ [[σ₁, 4]] := hz.1.1
      rw [uIcc_of_le (by linarith)] at hzre
      have hz0 : z ≠ 0 := by
        intro h
        have h0 : σ₁ ≤ z.re := hzre.1
        rw [h, zero_re] at h0
        linarith
      exact (Fψ_differentiableAt hZ hx (by linarith [hzre.1]) hz0 hz.2).differentiableWithinAt
  simp only [RectangleIntegral', RectangleIntegral, HIntegral, VIntegral, smul_eq_mul] at hres
  set Hb := ∫ t in σ₁..4, Fψ x (t + (-T : ℝ) * I)
  set Ht := ∫ t in σ₁..4, Fψ x (t + (T : ℝ) * I)
  set V4 := ∫ y in (-T)..T, Fψ x ((4 : ℝ) + y * I)
  set Vs := ∫ y in (-T)..T, Fψ x ((σ₁ : ℂ) + y * I)
  have hπ' : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_pos.ne'
  have h' : Hb - Ht + I * V4 - I * Vs = 2 * π * I * x := by
    field_simp at hres
    linear_combination hres
  push_cast
  field_simp
  linear_combination (-I) * h' + (V4 - Vs) * I_sq + (-2 * (π : ℂ) * x) * I_sq

/-! ### An elementary integral -/

lemma integral_inv_one_add_abs {T : ℝ} (hT : 0 ≤ T) :
    ∫ y in (-T)..T, (1 + |y|)⁻¹ = 2 * Real.log (1 + T) := by
  have hc : Continuous (fun y : ℝ => (1 + |y|)⁻¹) :=
    Continuous.inv₀ (by fun_prop) (fun y => by positivity)
  have h1 : ∫ y in (0 : ℝ)..T, (1 + |y|)⁻¹ = Real.log (1 + T) := by
    have e : ∫ y in (0 : ℝ)..T, (1 + |y|)⁻¹ = ∫ y in (0 : ℝ)..T, (fun u : ℝ => u⁻¹) (1 + y) := by
      refine intervalIntegral.integral_congr (fun y hy => ?_)
      rw [uIcc_of_le hT] at hy
      simp [abs_of_nonneg hy.1]
    rw [e, intervalIntegral.integral_comp_add_left (fun u : ℝ => u⁻¹) (1 : ℝ),
      integral_inv_of_pos (by norm_num) (by linarith)]
    simp
  have h2 : ∫ y in (-T)..0, (1 + |y|)⁻¹ = Real.log (1 + T) := by
    have e := intervalIntegral.integral_comp_neg (a := 0) (b := T) (fun y : ℝ => (1 + |y|)⁻¹)
    simp only [abs_neg, neg_zero] at e
    rw [← e, h1]
  rw [← intervalIntegral.integral_add_adjacent_intervals (b := 0) (hc.intervalIntegrable _ _)
    (hc.intervalIntegrable _ _), h1, h2]
  ring

end

end Fayed.PNTError

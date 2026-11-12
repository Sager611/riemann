import Mathlib.Analysis.Fourier.PoissonSummation
import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.MeasureTheory.Integral.DominatedConvergence

/-!
# Lamzouri's small-gap criterion has a floor: weak duality (explore34, N8)

Lamzouri, arXiv 2611.14486, Thm Main2 (SmallGapsPaper.tex l. 392-400): for 1/2 < θ < 4/7, 0 < λ ≤ 1 and
a nonnegative even compactly supported f ∈ C_c²(ℝ) with f^(0) = 1, f^(u) ≤ 0 for |u| > λ and
f(0) + 2∫_0^θ f(α) α dα > 1, one has D(λ) > 0. Fourier convention l. 307-309:
f^(ξ) = ∫ f(u) e^(-2πiξu) du, which is Mathlib's 𝓕.

* Admissible λ f: the f-hypotheses of Thm Main2 (Re f^ ≤ 0 is weaker than "f^ real and ≤ 0",
  so the barrier below covers every f Lamzouri allows).
* J θ f = f(0) + 2∫_0^θ f(α) α dα.
* weak_duality (Lemma 1): a lattice dual σ = Σ_{k≠0} h(1 + η_|k|) δ_{kh} with structure factor
  P ≥ |α|·1[|α| ≤ θ] on [0, 1/(2h)] forces J θ f ≤ 1 for every admissible f.
  Proof: Poisson summation (Mathlib), periodicity and evenness of P, σ ≥ 0 on |u| ≥ λ.
-/

open Real MeasureTheory FourierTransform Filter Asymptotics
open scoped Real

noncomputable section

namespace ExploreDPP

structure Admissible (lam : ℝ) (f : ℝ → ℝ) : Prop where
  contDiff : ContDiff ℝ 2 f
  hasCompactSupport : HasCompactSupport f
  even : ∀ x, f (-x) = f x
  nonneg : ∀ x, 0 ≤ f x
  fourier_zero : 𝓕 (fun x => (f x : ℂ)) 0 = 1
  fourier_nonpos : ∀ u, lam < |u| → (𝓕 (fun x => (f x : ℂ)) u).re ≤ 0

def J (θ : ℝ) (f : ℝ → ℝ) : ℝ := f 0 + 2 * ∫ α in (0:ℝ)..θ, f α * α

def FT (f : ℝ → ℝ) (u : ℝ) : ℝ := (𝓕 (fun x => (f x : ℂ)) u).re

variable {f : ℝ → ℝ}

lemma FT_eq_integral_cos (hc : Continuous f) (hcs : HasCompactSupport f) (u : ℝ) :
    FT f u = ∫ x, f x * Real.cos (2 * π * x * u) := by
  have hint : Integrable (fun v : ℝ => Complex.exp (↑(-2 * π * v * u) * Complex.I) • (f v : ℂ)) := by
    apply Continuous.integrable_of_hasCompactSupport
    · fun_prop
    · refine hcs.mono ?_
      intro x hx
      simp only [Function.mem_support] at hx ⊢
      intro h0; apply hx; simp [h0]
  rw [FT, Real.fourier_real_eq_integral_exp_smul, ← Complex.reCLM_apply,
    ← ContinuousLinearMap.integral_comp_comm _ hint]
  congr 1
  ext x
  simp only [Complex.reCLM_apply, smul_eq_mul, Complex.mul_re, Complex.exp_ofReal_mul_I_re,
    Complex.ofReal_re, Complex.ofReal_im, mul_zero, sub_zero]
  rw [mul_comm]
  congr 1
  rw [show -2 * π * x * u = -(2 * π * x * u) by ring, Real.cos_neg]

lemma continuous_FT (hi : Integrable f) : Continuous (FT f) := by
  have : Continuous (𝓕 (fun x => (f x : ℂ))) :=
    VectorFourier.fourierIntegral_continuous Real.continuous_fourierChar
      (innerSL ℝ).continuous₂ hi.ofReal
  exact Complex.continuous_re.comp this

lemma decay (hf : ContDiff ℝ 2 f) (hcs : HasCompactSupport f) :
    ∃ C, ∀ x : ℝ, (2 * π * |x|) ^ 2 * ‖𝓕 (fun x => (f x : ℂ)) x‖ ≤ C := by
  set g : ℝ → ℂ := fun x => (f x : ℂ) with hgdef
  have h0 : ContDiff ℝ 2 ((↑) : ℝ → ℂ) := Complex.ofRealCLM.contDiff
  have hg : ContDiff ℝ 2 g := h0.comp hf
  have hgcs : HasCompactSupport g := hcs.comp_left Complex.ofReal_zero
  have hint : ∀ n : ℕ, n ≤ (2 : ℕ∞) → Integrable (iteratedDeriv n g) := by
    intro n hn
    apply (hg.continuous_iteratedDeriv n (by exact_mod_cast hn)).integrable_of_hasCompactSupport
    refine (hgcs.iteratedFDeriv (𝕜 := ℝ) n).mono ?_
    intro x hx
    simp only [Function.mem_support] at hx ⊢
    intro h0; apply hx; simp [iteratedDeriv, h0]
  have hF := Real.fourier_iteratedDeriv (N := 2) (n := 2) (by exact_mod_cast hg) hint le_rfl
  refine ⟨∫ x, ‖iteratedDeriv 2 g x‖, fun x => ?_⟩
  have h1 := VectorFourier.norm_fourierIntegral_le_integral_norm 𝐞 volume (innerₗ ℝ)
    (iteratedDeriv 2 g) x
  change ‖𝓕 (iteratedDeriv 2 g) x‖ ≤ _ at h1
  rw [hF, norm_smul, norm_pow] at h1
  have hn : ‖(2 * π * Complex.I * (x:ℂ))‖ = 2 * π * |x| := by
    simp [abs_of_pos Real.pi_pos]
  rw [hn] at h1
  exact h1


lemma FT_neg (hc : Continuous f) (hcs : HasCompactSupport f) (u : ℝ) : FT f (-u) = FT f u := by
  rw [FT_eq_integral_cos hc hcs, FT_eq_integral_cos hc hcs]
  congr 1; funext x
  rw [show 2 * π * x * -u = -(2 * π * x * u) by ring, Real.cos_neg]

lemma abs_FT_le (u : ℝ) : |FT f u| ≤ ∫ x, |f x| := by
  have h1 := VectorFourier.norm_fourierIntegral_le_integral_norm 𝐞 volume (innerₗ ℝ)
    (fun x => (f x : ℂ)) u
  change ‖𝓕 (fun x => (f x : ℂ)) u‖ ≤ _ at h1
  simp only [Complex.norm_real, Real.norm_eq_abs] at h1
  exact (Complex.abs_re_le_norm _).trans h1

lemma FT_nonpos_of_ge {lam : ℝ} (hf : Admissible lam f) {u : ℝ} (hu : lam ≤ u) : FT f u ≤ 0 := by
  have hc : Continuous (FT f) := continuous_FT
    (hf.contDiff.continuous.integrable_of_hasCompactSupport hf.hasCompactSupport)
  have hcl : IsClosed {u | FT f u ≤ 0} := isClosed_le hc continuous_const
  have hsub : Set.Ioi lam ⊆ {u | FT f u ≤ 0} := fun u hu' =>
    hf.fourier_nonpos u (lt_of_lt_of_le hu' (le_abs_self u))
  have := closure_minimal hsub hcl
  rw [closure_Ioi] at this
  exact this hu

/-- Poisson summation on the lattice: sum f(m/h) = h * sum FT f (k h). -/
lemma poisson (hf : ContDiff ℝ 2 f) (hcs : HasCompactSupport f) {h : ℝ} (hh : 0 < h) :
    Summable (fun m : ℤ => f (m / h)) ∧ Summable (fun k : ℤ => FT f (k * h)) ∧
    ∑' m : ℤ, f (m / h) = h * ∑' k : ℤ, FT f (k * h) := by
  set g : ℝ → ℂ := fun x => (f x : ℂ) with hgdef
  set G : ℝ → ℂ := fun x => g (x / h) with hGdef
  have hgc : Continuous g := Complex.continuous_ofReal.comp hf.continuous
  have hGc : Continuous G := hgc.comp (continuous_id.div_const h)
  have hFG : ∀ ξ : ℝ, 𝓕 G ξ = (h : ℂ) * 𝓕 g (h * ξ) := by
    intro ξ
    rw [Real.fourier_real_eq, Real.fourier_real_eq]
    have := MeasureTheory.Measure.integral_comp_div
      (fun y : ℝ => 𝐞 (-((h * y) * ξ)) • g y) h
    rw [abs_of_pos hh] at this
    have e1 : (fun v : ℝ => 𝐞 (-(v * ξ)) • G v) =
        fun v => 𝐞 (-((h * (v / h)) * ξ)) • g (v / h) := by
      funext v; congr 3; field_simp
    rw [e1, this, Complex.real_smul]
    congr 1
    congr 1; funext y; congr 3; ring
  obtain ⟨R, hR⟩ := hcs.isCompact.isBounded.subset_closedBall 0
  have hGO : G =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
    have hev : G =ᶠ[cocompact ℝ] 0 := by
      have hmem := (isCompact_closedBall (0:ℝ) (|R| * h)).compl_mem_cocompact
      filter_upwards [hmem] with x hx
      simp only [Set.mem_compl_iff, Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs,
        not_le] at hx
      simp only [hGdef, hgdef, Pi.zero_apply, Complex.ofReal_eq_zero]
      apply image_eq_zero_of_notMem_tsupport
      intro hxs
      have h1 := hR hxs
      simp only [Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs] at h1
      rw [abs_div, abs_of_pos hh, div_le_iff₀ hh] at h1
      nlinarith [le_abs_self R]
    exact hev.isBigO.trans (isBigO_zero _ _)
  obtain ⟨C, hC⟩ := decay hf hcs
  have hFGO : (𝓕 G) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
    refine IsBigO.of_bound (C / (4 * π ^ 2 * h)) ?_
    have hmem : ({0}ᶜ : Set ℝ) ∈ cocompact ℝ := isCompact_singleton.compl_mem_cocompact
    filter_upwards [hmem] with x hx
    have hx0 : x ≠ 0 := hx
    rw [hFG x, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hh, Real.norm_eq_abs,
      abs_of_nonneg (by positivity : (0:ℝ) ≤ |x| ^ (-(2:ℝ))),
      Real.rpow_neg (abs_nonneg x), Real.rpow_two, sq_abs]
    have h2 := hC (h * x)
    have hhx : 0 < |h * x| := abs_pos.mpr (mul_ne_zero hh.ne' hx0)
    have hpos : 0 < (2 * π * |h * x|) ^ 2 := by positivity
    have h3 : ‖𝓕 g (h * x)‖ ≤ C / (2 * π * |h * x|) ^ 2 := by
      rw [le_div_iff₀ hpos]; linarith
    calc h * ‖𝓕 g (h * x)‖ ≤ h * (C / (2 * π * |h * x|) ^ 2) := by gcongr
      _ = C / (4 * π ^ 2 * h) * (x ^ 2)⁻¹ := by
        rw [abs_mul, abs_of_pos hh, show (2 * π * (h * |x|)) ^ 2 = 4 * π ^ 2 * h ^ 2 * x ^ 2 by
          rw [mul_pow, mul_pow, mul_pow, sq_abs]; ring]
        have := Real.pi_pos.ne'
        field_simp
  have hsG : Summable (fun n : ℤ => G n) :=
    summable_of_isBigO (Real.summable_abs_int_rpow one_lt_two)
      (hGO.comp_tendsto Int.tendsto_coe_cofinite)
  have hsFG : Summable (fun n : ℤ => 𝓕 G n) :=
    summable_of_isBigO (Real.summable_abs_int_rpow one_lt_two)
      (hFGO.comp_tendsto Int.tendsto_coe_cofinite)
  have hpois := Real.tsum_eq_tsum_fourier_of_rpow_decay hGc one_lt_two hGO hFGO 0
  simp only [zero_add, fourier_coe_apply, Complex.ofReal_zero, mul_zero, zero_div,
    Complex.exp_zero, mul_one] at hpois
  have hre : ∀ n : ℤ, (𝓕 G n).re = h * FT f (n * h) := by
    intro n
    rw [hFG, FT, Complex.re_ofReal_mul, mul_comm (h : ℝ) (n : ℝ)]
  refine ⟨?_, ?_, ?_⟩
  · exact Complex.summable_ofReal.mp hsG
  · have := (Complex.reCLM.summable hsFG)
    simp only [Complex.reCLM_apply, hre] at this
    exact (summable_mul_left_iff hh.ne').mp this
  · have e2 : ((∑' m : ℤ, f (m / h) : ℝ) : ℂ) = ∑' n : ℤ, G n := by
      rw [Complex.ofReal_tsum]
    have e3 := congrArg Complex.re (e2.trans hpois)
    rw [Complex.ofReal_re, Complex.re_tsum hsFG] at e3
    simp only [hre] at e3
    rw [e3, tsum_mul_left]

lemma tsum_int_fold {a : ℤ → ℝ} (hev : ∀ n, a (-n) = a n) (hs : Summable a) :
    ∑' n : ℤ, a n = a 0 + 2 * ∑' n : ℕ, a ((n : ℤ) + 1) := by
  rw [tsum_int_eq_zero_add_two_mul_tsum_pnat (f := a) (fun n => hev n) hs,
    tsum_pnat_eq_tsum_succ (f := fun n : ℕ => a n)]
  simp [nsmul_eq_mul]


/-- The structure factor P(α) = (1 - h) + 2h Σ_{k≥1} η_k cos(2π α k h), with k = n + 1. -/
def P (h : ℝ) (η : ℕ → ℝ) (α : ℝ) : ℝ :=
  (1 - h) + 2 * h * ∑' n : ℕ, η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))

section P
variable {h : ℝ} {η : ℕ → ℝ}

lemma summable_shift (hη : Summable η) : Summable (fun n : ℕ => |η (n + 1)|) :=
  ((summable_nat_add_iff 1).mpr hη).abs

lemma continuous_S (hη : Summable η) :
    Continuous (fun α : ℝ => ∑' n : ℕ, η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))) := by
  refine continuous_tsum (fun n => by fun_prop) (summable_shift hη) (fun n α => ?_)
  rw [Real.norm_eq_abs, abs_mul]
  exact mul_le_of_le_one_right (abs_nonneg _) (Real.abs_cos_le_one _)

lemma continuous_P (hη : Summable η) : Continuous (P h η) := by
  have := continuous_S (h := h) hη
  unfold P
  fun_prop

lemma P_neg (α : ℝ) : P h η (-α) = P h η α := by
  simp only [P]
  congr 2
  apply tsum_congr; intro n
  rw [show 2 * π * -α * ((n + 1) * h) = -(2 * π * α * ((n + 1) * h)) by ring, Real.cos_neg]

lemma P_add_period (hh : 0 < h) (α : ℝ) (m : ℤ) : P h η (α + m / h) = P h η α := by
  simp only [P]
  congr 2
  apply tsum_congr; intro n
  congr 1
  have : 2 * π * (α + m / h) * ((n + 1) * h) =
      2 * π * α * ((n + 1) * h) + ((m * (n + 1) : ℤ) : ℝ) * (2 * π) := by
    push_cast; field_simp
  rw [this, Real.cos_add_int_mul_two_pi]

lemma P_ge_of_base {θ : ℝ} (hh : 0 < h) (hθ : θ ≤ 1 / (2 * h))
    (hP : ∀ α, 0 ≤ α → α ≤ 1 / (2 * h) → (if α ≤ θ then α else 0) ≤ P h η α) :
    (∀ α, 0 ≤ P h η α) ∧ (∀ α, |α| ≤ θ → |α| ≤ P h η α) := by
  have habs : ∀ α, P h η α = P h η |α| := by
    intro α; rcases le_or_gt 0 α with h0 | h0
    · rw [abs_of_nonneg h0]
    · rw [abs_of_neg h0, P_neg]
  constructor
  · intro α
    set m : ℤ := round (α * h) with hm
    have hb : |α - m / h| ≤ 1 / (2 * h) := by
      have h1 := abs_sub_round (α * h)
      rw [← hm] at h1
      have e : α - m / h = (α * h - m) / h := by field_simp
      rw [e, abs_div, abs_of_pos hh, div_le_div_iff₀ hh (by positivity)]
      nlinarith [h1]
    have e2 : P h η α = P h η (α - m / h) := by
      conv_lhs => rw [← sub_add_cancel α ((m : ℝ) / h)]
      exact P_add_period hh _ m
    rw [e2, habs]
    have := hP |α - m / h| (abs_nonneg _) hb
    split_ifs at this <;> linarith [abs_nonneg (α - m / h)]
  · intro α hα
    rw [habs]
    have := hP |α| (abs_nonneg _) (hα.trans hθ)
    simpa [hα] using this

lemma integral_mul_P (hf : ContDiff ℝ 2 f) (hcs : HasCompactSupport f) (hη : Summable η) :
    ∫ α, f α * P h η α =
      (1 - h) * FT f 0 + 2 * h * ∑' n : ℕ, η (n + 1) * FT f ((n + 1) * h) := by
  have hc := hf.continuous
  have hi : Integrable f := hc.integrable_of_hasCompactSupport hcs
  set S : ℝ → ℝ := fun α => ∑' n : ℕ, η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h)) with hS
  have hSc : Continuous S := continuous_S hη
  have hfS : Integrable (fun α => f α * S α) :=
    (hc.mul hSc).integrable_of_hasCompactSupport hcs.mul_right
  have e1 : (fun α => f α * P h η α) = fun α => (1 - h) * f α + 2 * h * (f α * S α) := by
    funext α; simp only [P, hS]; ring
  rw [e1, integral_add (hi.const_mul _) (hfS.const_mul _), integral_const_mul, integral_const_mul]
  have e0 : ∫ α, f α = FT f 0 := by
    rw [FT_eq_integral_cos hc hcs]; simp
  rw [e0]
  congr 2
  have hFi : ∀ n : ℕ, Integrable
      (fun α => f α * (η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h)))) :=
    fun n => (hc.mul (by fun_prop)).integrable_of_hasCompactSupport hcs.mul_right
  have hFs : Summable (fun n : ℕ =>
      ∫ α, ‖f α * (η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h)))‖) := by
    refine Summable.of_nonneg_of_le (fun n => integral_nonneg (fun α => norm_nonneg _))
      (fun n => ?_) ((summable_shift hη).mul_right (∫ α, ‖f α‖))
    rw [← integral_const_mul]
    apply integral_mono (hFi n).norm (hi.norm.const_mul _)
    intro α
    simp only [norm_mul, Real.norm_eq_abs]
    have := Real.abs_cos_le_one (2 * π * α * ((n + 1) * h))
    nlinarith [mul_nonneg (abs_nonneg (f α)) (abs_nonneg (η (n + 1)))]
  have key := integral_tsum_of_summable_integral_norm hFi hFs
  have e3 : (fun α => f α * S α) = fun α =>
      ∑' n : ℕ, f α * (η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))) := by
    funext α; simp only [hS]; rw [tsum_mul_left]
  rw [e3, ← key]
  apply tsum_congr; intro n
  rw [FT_eq_integral_cos hc hcs, ← integral_const_mul]
  congr 1; funext α; ring

end P

/-- Lemma 1 (weak duality). A lattice dual η whose structure factor P dominates
|α|·1[|α| ≤ θ] on the base period [0, 1/(2h)] gives J_θ(f) ≤ 1 = f^(0) for every
admissible f. -/
theorem weak_duality {h lam θ : ℝ} {k0 : ℕ} {η : ℕ → ℝ}
    (hh : 0 < h) (hlam : (k0 : ℝ) * h = lam) (hη : Summable η)
    (hcore : ∀ k : ℕ, 1 ≤ k → k < k0 → η k = -1) (hge : ∀ k : ℕ, k0 ≤ k → -1 ≤ η k)
    (hθ0 : 0 ≤ θ) (hθ : θ ≤ 1 / (2 * h))
    (hP : ∀ α, 0 ≤ α → α ≤ 1 / (2 * h) → (if α ≤ θ then α else 0) ≤ P h η α)
    {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1 := by
  obtain ⟨hP0, hPθ⟩ := P_ge_of_base hh hθ hP
  have hcf : Continuous f := hf.contDiff.continuous
  have hPc : Continuous (P h η) := continuous_P hη
  have hint : Integrable (fun α => f α * P h η α) :=
    (hcf.mul hPc).integrable_of_hasCompactSupport hf.hasCompactSupport.mul_right
  have hnn : ∀ α, 0 ≤ f α * P h η α := fun α => mul_nonneg (hf.nonneg α) (hP0 α)
  -- Step A: J θ f ≤ f 0 + ∫ f P
  have hA : J θ f ≤ f 0 + ∫ α, f α * P h η α := by
    unfold J
    have h1 : ∫ α in (0:ℝ)..θ, f α * α ≤ ∫ α in (0:ℝ)..θ, f α * P h η α := by
      apply intervalIntegral.integral_mono_on hθ0
      · exact (hcf.mul continuous_id).intervalIntegrable _ _
      · exact (hcf.mul hPc).intervalIntegrable _ _
      · intro α hα
        apply mul_le_mul_of_nonneg_left _ (hf.nonneg α)
        have := hPθ α (by rw [abs_of_nonneg hα.1]; exact hα.2)
        rwa [abs_of_nonneg hα.1] at this
    have hsym : ∫ α in (-θ)..0, f α * P h η α = ∫ α in (0:ℝ)..θ, f α * P h η α := by
      have := intervalIntegral.integral_comp_neg (a := 0) (b := θ) (fun α => f α * P h η α)
      simp only [neg_zero, hf.even, P_neg] at this
      exact this.symm
    have hadd := intervalIntegral.integral_add_adjacent_intervals (μ := volume)
      ((hcf.mul hPc).intervalIntegrable (-θ) 0) ((hcf.mul hPc).intervalIntegrable 0 θ)
    simp only [Pi.mul_apply] at hadd
    have h3 : ∫ α in (-θ)..θ, f α * P h η α ≤ ∫ α, f α * P h η α := by
      rw [intervalIntegral.integral_of_le (by linarith)]
      exact setIntegral_le_integral hint (Eventually.of_forall hnn)
    linarith
  -- Step B: Poisson
  obtain ⟨hs1, hs2, hpois⟩ := poisson hf.contDiff hf.hasCompactSupport hh
  have hB : f 0 ≤ ∑' m : ℤ, f (m / h) := by
    have := hs1.le_tsum 0 (fun j _ => hf.nonneg _)
    simpa using this
  have hFTeven : ∀ k : ℤ, FT f ((-k : ℤ) * h) = FT f (k * h) := by
    intro k
    rw [Int.cast_neg, neg_mul, FT_neg hcf hf.hasCompactSupport]
  have hfold : ∑' k : ℤ, FT f (k * h) = FT f 0 + 2 * ∑' n : ℕ, FT f ((n + 1) * h) := by
    rw [tsum_int_fold (a := fun k : ℤ => FT f (k * h)) hFTeven hs2]
    push_cast
    simp
  -- Step C: the integral
  have hC := integral_mul_P (h := h) hf.contDiff hf.hasCompactSupport hη
  have hF0 : FT f 0 = 1 := by simp [FT, hf.fourier_zero]
  -- Step D: termwise sign
  have hsumF : Summable (fun n : ℕ => FT f ((n + 1) * h)) := by
    have hinj : Function.Injective (fun n : ℕ => (n : ℤ) + 1) := by
      intro a b hab; simpa using hab
    have := hs2.comp_injective hinj
    refine this.congr (fun n => ?_)
    simp
  have hsumηF : Summable (fun n : ℕ => η (n + 1) * FT f ((n + 1) * h)) := by
    refine Summable.of_norm_bounded ((summable_shift hη).mul_right (∫ x, |f x|)) (fun n => ?_)
    rw [Real.norm_eq_abs, abs_mul]
    exact mul_le_mul_of_nonneg_left (abs_FT_le _) (abs_nonneg _)
  have hD : ∑' n : ℕ, FT f ((n + 1) * h) + ∑' n : ℕ, η (n + 1) * FT f ((n + 1) * h) ≤ 0 := by
    rw [← hsumF.tsum_add hsumηF]
    apply tsum_nonpos
    intro n
    have e : FT f ((n + 1) * h) + η (n + 1) * FT f ((n + 1) * h) =
        (1 + η (n + 1)) * FT f ((n + 1) * h) := by ring
    rw [e]
    by_cases hn : n + 1 < k0
    · rw [hcore (n + 1) (by omega) hn]; simp
    · rw [not_lt] at hn
      apply mul_nonpos_of_nonneg_of_nonpos
      · linarith [hge (n + 1) hn]
      · apply FT_nonpos_of_ge hf
        rw [← hlam]
        have : (k0 : ℝ) ≤ (n : ℝ) + 1 := by exact_mod_cast hn
        nlinarith
  have hE : 2 * h * (∑' n : ℕ, FT f ((n + 1) * h) +
      ∑' n : ℕ, η (n + 1) * FT f ((n + 1) * h)) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (by positivity) hD
  calc J θ f ≤ f 0 + ∫ α, f α * P h η α := hA
    _ ≤ ∑' m : ℤ, f (m / h) + ∫ α, f α * P h η α := by linarith
    _ = FT f 0 + 2 * h * (∑' n : ℕ, FT f ((n + 1) * h) +
          ∑' n : ℕ, η (n + 1) * FT f ((n + 1) * h)) := by
      rw [hpois, hfold, hC]; ring
    _ ≤ FT f 0 := by linarith
    _ = 1 := hF0


/-! ## Monotonicity -/

lemma Admissible.mono {lam lam' : ℝ} (hf : Admissible lam' f) (h : lam' ≤ lam) :
    Admissible lam f :=
  { hf with fourier_nonpos := fun u hu => hf.fourier_nonpos u (lt_of_le_of_lt h hu) }

lemma J_mono (hnn : ∀ x, 0 ≤ f x) (hc : Continuous f) {θ' θ : ℝ} (h0 : 0 ≤ θ') (h : θ' ≤ θ) :
    J θ' f ≤ J θ f := by
  unfold J
  have hadd := intervalIntegral.integral_add_adjacent_intervals (μ := volume)
    ((hc.mul continuous_id).intervalIntegrable 0 θ') ((hc.mul continuous_id).intervalIntegrable θ' θ)
  simp only [Pi.mul_apply, id] at hadd
  have hnn2 : 0 ≤ ∫ α in θ'..θ, f α * α :=
    intervalIntegral.integral_nonneg h (fun α hα => mul_nonneg (hnn α) (by linarith [hα.1]))
  linarith

/-- f satisfies every f-hypothesis of Thm Main2 at (θ, λ), including J_θ(f) > 1. -/
def Main2Kernel (θ lam : ℝ) (f : ℝ → ℝ) : Prop := Admissible lam f ∧ 1 < J θ f

/-- All hypotheses of Lamzouri's Thm Main2 (l. 392-396) for the data (θ, λ, f). -/
def Main2Hyp (θ lam : ℝ) (f : ℝ → ℝ) : Prop :=
  1 / 2 < θ ∧ θ < 4 / 7 ∧ 0 < lam ∧ lam ≤ 1 ∧ Main2Kernel θ lam f

lemma Main2Kernel.mono {θ' θ lam' lam : ℝ} (hf : Main2Kernel θ' lam' f) (h0 : 0 ≤ θ')
    (hθ : θ' ≤ θ) (hlam : lam' ≤ lam) : Main2Kernel θ lam f :=
  ⟨hf.1.mono hlam, lt_of_lt_of_le hf.2 (J_mono hf.1.nonneg hf.1.contDiff.continuous h0 hθ)⟩

end ExploreDPP

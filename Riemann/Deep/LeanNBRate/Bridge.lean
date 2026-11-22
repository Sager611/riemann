import Riemann.Deep.LeanNBRate.Closure
import Riemann.Deep.LeanNBRate.U2aux

/-!
# The lower bound (Burnol), real versus complex spans, the constrained distance, the value
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- Complexification L²(0,∞; ℝ) → L²(0,∞; ℂ). -/
def cplx (f : L2pos) : L2posC := (Complex.ofRealCLM : ℝ →L[ℝ] ℂ).compLp f

theorem cplx_ae (f : L2pos) : (cplx f : ℝ → ℂ) =ᵐ[volPos] fun x => ((f x : ℝ) : ℂ) :=
  ContinuousLinearMap.coeFn_compLp' _ f

theorem cplx_sub (f g : L2pos) : cplx f - cplx g = cplx (f - g) := by
  unfold cplx
  exact ((Complex.ofRealCLM : ℝ →L[ℝ] ℂ).compLpₗ 2 volPos).map_sub f g |>.symm

theorem cplx_add (f g : L2pos) : cplx (f + g) = cplx f + cplx g := by
  unfold cplx
  exact ((Complex.ofRealCLM : ℝ →L[ℝ] ℂ).compLpₗ 2 volPos).map_add f g

theorem cplx_smul (r : ℝ) (f : L2pos) : cplx (r • f) = (r : ℂ) • cplx f := by
  unfold cplx
  rw [Complex.coe_smul]
  exact ((Complex.ofRealCLM : ℝ →L[ℝ] ℂ).compLpₗ 2 volPos).map_smul r f

theorem norm_cplx (f : L2pos) : ‖cplx f‖ = ‖f‖ := by
  have h : ∀ᵐ x ∂volPos, ‖(cplx f : ℝ → ℂ) x‖ = ‖(f : ℝ → ℝ) x‖ := by
    filter_upwards [cplx_ae f] with x hx
    rw [hx, Complex.norm_real]
  rw [Lp.norm_def (cplx f), Lp.norm_def f,
    eLpNorm_congr_norm_ae (Lp.aestronglyMeasurable _) (Lp.aestronglyMeasurable _) h]

theorem isometry_cplx : Isometry cplx := by
  refine Isometry.of_dist_eq fun f g => ?_
  rw [dist_eq_norm, dist_eq_norm, cplx_sub, norm_cplx]

theorem cplx_chiVec : cplx chiVec = chiVecC := by
  refine Lp.ext ?_
  have h1 := cplx_ae chiVec
  have h2 : (chiVec : ℝ → ℝ) =ᵐ[volPos] (Ioo (0 : ℝ) 1).indicator (fun _ => (1 : ℝ)) := indicatorConstLp_coeFn
  have h3 : (chiVecC : ℝ → ℂ) =ᵐ[volPos] (Ioc (0 : ℝ) 1).indicator (fun _ => (1 : ℂ)) := indicatorConstLp_coeFn
  have hne : ∀ᵐ x ∂volPos, x ≠ 1 := by
    refine ae_restrict_of_ae ?_
    rw [ae_iff]; simp [measure_singleton]
  filter_upwards [h1, h2, h3, hne] with x hx1 hx2 hx3 hx
  rw [hx1, hx2, hx3]
  by_cases h : x ∈ Ioo (0 : ℝ) 1
  · rw [indicator_of_mem h, indicator_of_mem (Ioo_subset_Ioc_self h)]; simp
  · rw [indicator_of_notMem h]
    have : x ∉ Ioc (0 : ℝ) 1 := fun h' => h ⟨h'.1, lt_of_le_of_ne h'.2 hx⟩
    rw [indicator_of_notMem this]; simp

theorem cplx_span (lam : ℝ) :
    cplx '' (Submodule.span ℝ (nbSet lam) : Set L2pos) ⊆ (Submodule.span ℂ (nbSetC lam) : Set L2posC) := by
  rintro _ ⟨f, hf, rfl⟩
  induction hf using Submodule.span_induction with
  | mem f hf =>
    obtain ⟨θ, hθ, hfθ⟩ := hf
    apply Submodule.subset_span
    refine ⟨θ, hθ, ?_⟩
    filter_upwards [cplx_ae f, hfθ] with x hx1 hx2
    rw [hx1, hx2]; rfl
  | zero =>
    have : cplx 0 = 0 := by
      unfold cplx; exact ((Complex.ofRealCLM : ℝ →L[ℝ] ℂ).compLpₗ 2 volPos).map_zero
    rw [this]; exact Submodule.zero_mem _
  | add f g _ _ hf hg => rw [cplx_add]; exact Submodule.add_mem _ hf hg
  | smul r f _ hf => rw [cplx_smul]; exact Submodule.smul_mem _ _ hf

/-- Burnol's complex distance is at most the real distance. -/
theorem nbDistC_le (lam : ℝ) : nbDistC lam ≤ nbDist lam := by
  unfold nbDistC nbDist
  rw [← cplx_chiVec, ← Metric.infDist_image isometry_cplx]
  exact Metric.infDist_le_infDist_of_subset (cplx_span lam)
    ⟨cplx 0, ⟨0, Submodule.zero_mem _, rfl⟩⟩

/-- The constrained set lies in the span. -/
theorem nymanSet_subset (lam : ℝ) (hlam : 0 < lam) : nymanSet lam ⊆ (Submodule.span ℝ (nbSet lam) : Set L2pos) := by
  rintro f ⟨n, c, θ, hθ, _, hf⟩
  have hθ' : ∀ k, θ k ∈ Ioc (0 : ℝ) 1 := fun k => ⟨by linarith [(hθ k).1], (hθ k).2⟩
  have hmem : ∀ k, MemLp (nbFun (θ k)) 2 volPos := fun k => memLp_nbFun (hθ' k)
  have heq : f = ∑ k, c k • (hmem k).toLp := by
    refine Lp.ext ?_
    have hs := Lp.coeFn_finsetSum (μ := volPos) (p := 2) Finset.univ (fun k => c k • (hmem k).toLp)
    have hk : ∀ k, ((c k • (hmem k).toLp : L2pos) : ℝ → ℝ) =ᵐ[volPos] fun x => c k * nbFun (θ k) x := by
      intro k
      filter_upwards [Lp.coeFn_smul (c k) (hmem k).toLp, MemLp.coeFn_toLp (hmem k)] with x hx1 hx2
      rw [hx1, Pi.smul_apply, hx2, smul_eq_mul]
    have hall : ∀ᵐ x ∂volPos, ∀ k, ((c k • (hmem k).toLp : L2pos) : ℝ → ℝ) x = c k * nbFun (θ k) x :=
      ae_all_iff.mpr hk
    filter_upwards [hf, hs, hall] with x hx1 hx2 hx3
    rw [hx1, hx2, Finset.sum_apply]
    exact (Finset.sum_congr rfl fun k _ => (hx3 k).symm)
  rw [heq]
  refine Submodule.sum_mem _ fun k _ => Submodule.smul_mem _ _ (Submodule.subset_span ?_)
  exact ⟨θ k, hθ k, MemLp.coeFn_toLp (hmem k)⟩

theorem nbDist_le_nymanDist (lam : ℝ) (hlam : 0 < lam) : nbDist lam ≤ nymanDist lam := by
  unfold nbDist nymanDist
  refine Metric.infDist_le_infDist_of_subset (nymanSet_subset lam hlam) ⟨0, ?_⟩
  refine ⟨0, Fin.elim0, Fin.elim0, fun k => Fin.elim0 k, by simp, ?_⟩
  filter_upwards [Lp.coeFn_zero ℝ 2 volPos] with x hx
  rw [hx]; simp

/-- **The lower bound** (Burnol 2002): under RH, L·D(λ)² is eventually ≥ S − ε. -/
theorem lower_bound (hB : Burnol2002_Thm) (hc : Titchmarsh_Thm9_2) (hRH : RiemannHypothesis) (ε : ℝ) (hε : 0 < ε) :
    ∀ᶠ lam in 𝓝[>] (0 : ℝ), zeroSum - ε ≤ Real.log (1 / lam) * nbDist lam ^ 2 := by
  have hS := summable_zeroSum hc hRH
  have hZ := zeroSum_nonneg
  have hev1 : ∀ᶠ lam in 𝓝[>] (0 : ℝ), lam < 1 := by
    filter_upwards [Ioo_mem_nhdsGT (show (0 : ℝ) < 1 by norm_num)] with lam hlam using hlam.2
  have hev0 : ∀ᶠ lam in 𝓝[>] (0 : ℝ), 0 < lam := self_mem_nhdsWithin
  by_cases hcase : zeroSum - ε ≤ 0
  · filter_upwards [hev1, hev0] with lam h1 h0
    have : 0 ≤ Real.log (1 / lam) := Real.log_nonneg (by rw [le_div_iff₀ h0]; linarith)
    have : 0 ≤ Real.log (1 / lam) * nbDist lam ^ 2 := by positivity
    linarith
  push Not at hcase
  set c := Real.sqrt (zeroSum - ε)
  have hc0 : 0 ≤ c := Real.sqrt_nonneg _
  have hlt : ENNReal.ofReal c <
      (∑' ρ : NTZ, ENNReal.ofReal ((mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2)) ^ (1 / 2 : ℝ) := by
    rw [← ENNReal.ofReal_tsum_of_nonneg (fun ρ => by positivity) hS,
      ENNReal.ofReal_rpow_of_nonneg (tsum_nonneg fun ρ => by positivity) (by norm_num)]
    have hZ' : 0 < ∑' ρ : NTZ, (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2 := by
      have : zeroSum = ∑' ρ : NTZ, (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2 := rfl
      rw [← this]; linarith
    rw [ENNReal.ofReal_lt_ofReal_iff (Real.rpow_pos_of_pos hZ' _)]
    rw [← Real.sqrt_eq_rpow]
    have hzs : zeroSum = ∑' ρ : NTZ, (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2 := rfl
    exact Real.sqrt_lt_sqrt (by linarith) (by rw [← hzs]; linarith)
  filter_upwards [hB hRH c hlt, hev1, hev0] with lam hlam h1 h0
  have hL : 0 ≤ Real.log (1 / lam) := Real.log_nonneg (by rw [le_div_iff₀ h0]; linarith)
  have h2 : c ≤ nbDist lam * Real.sqrt (Real.log (1 / lam)) :=
    hlam.trans (mul_le_mul_of_nonneg_right (nbDistC_le lam) (Real.sqrt_nonneg _))
  have h3 := pow_le_pow_left₀ hc0 h2 2
  rw [mul_pow, Real.sq_sqrt hL, Real.sq_sqrt (by linarith)] at h3
  linarith

/-- **The value**: with simple zeros, S = 2 + γ − log 4π (Davenport ch. 12, RH). -/
theorem zeroSum_value (hD : Davenport_Ch12) (hRH : RiemannHypothesis) (hsimple : ∀ ρ : NTZ, mult (ρ : ℂ) = 1) :
    zeroSum = 2 + Real.eulerMascheroniConstant - Real.log (4 * Real.pi) := by
  have e : ∀ ρ : NTZ, (mult (ρ : ℂ) : ℝ) * ((ρ : ℂ)⁻¹).re = (1 / 2) * ((mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2) := by
    intro ρ
    have hre := NTZ_re_half hRH ρ.2
    have hne : (ρ : ℂ) ≠ 0 := fun h => by rw [h] at hre; simp at hre
    rw [Complex.inv_re, hre, hsimple ρ, Complex.normSq_eq_norm_sq]
    push_cast
    have : ‖(ρ : ℂ)‖ ≠ 0 := norm_ne_zero_iff.mpr hne
    field_simp
  have h : HasSum (fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * ((ρ : ℂ)⁻¹).re)
      (1 + Real.eulerMascheroniConstant / 2 - Real.log (4 * Real.pi) / 2) := hD
  simp_rw [e] at h
  have h2 := h.mul_left 2
  have h3 : HasSum (fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2)
      (2 * (1 + Real.eulerMascheroniConstant / 2 - Real.log (4 * Real.pi) / 2)) := by
    convert h2 using 1
    funext ρ; ring
  unfold zeroSum
  rw [h3.tsum_eq]
  ring

end Riemann.Deep.LeanNBRate


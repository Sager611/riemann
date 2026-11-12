import Riemann.Deep.ExploreDPP.Lean.TheoremB
import Riemann.Deep.ExploreDPP.Lean.Triangle

/-!
# Decaying kernels (every Schwartz f): the barrier also covers SmallGap's Lemma T route

AdmissibleDecay replaces C_c² by: f continuous, integrable, f and f^ = O(|x|^(-2)). Every
Admissible f and every real Schwartz f with Lamzouri's sign conditions is AdmissibleDecay;
weak duality and Theorems A, B hold verbatim for this class (DERIVATION: "Schwartz f: same proof").
-/

open Real MeasureTheory FourierTransform Filter Asymptotics
open scoped Real SchwartzMap

noncomputable section

namespace ExploreDPP

structure AdmissibleDecay (lam : ℝ) (f : ℝ → ℝ) : Prop where
  continuous : Continuous f
  integrable : Integrable f
  f_decay : f =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ)))
  fourier_decay : 𝓕 (fun x => (f x : ℂ)) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ)))
  even : ∀ x, f (-x) = f x
  nonneg : ∀ x, 0 ≤ f x
  fourier_zero : 𝓕 (fun x => (f x : ℂ)) 0 = 1
  fourier_nonpos : ∀ u, lam < |u| → (𝓕 (fun x => (f x : ℂ)) u).re ≤ 0

variable {f : ℝ → ℝ}

lemma AdmissibleDecay.mono {lam lam' : ℝ} (hf : AdmissibleDecay lam' f) (h : lam' ≤ lam) :
    AdmissibleDecay lam f :=
  { hf with fourier_nonpos := fun u hu => hf.fourier_nonpos u (lt_of_le_of_lt h hu) }

lemma integrable_exp_smul (hc : Continuous f) (hi : Integrable f) (u : ℝ) :
    Integrable (fun v : ℝ => Complex.exp (↑(-2 * π * v * u) * Complex.I) • (f v : ℂ)) := by
  have hm : Continuous (fun v : ℝ => Complex.exp (↑(-2 * π * v * u) * Complex.I) • (f v : ℂ)) := by
    fun_prop
  refine hi.norm.mono' hm.aestronglyMeasurable (Eventually.of_forall fun v => ?_)
  rw [norm_smul, Complex.norm_exp_ofReal_mul_I, one_mul, Complex.norm_real]

lemma FT_eq_integral_cos' (hc : Continuous f) (hi : Integrable f) (u : ℝ) :
    FT f u = ∫ x, f x * Real.cos (2 * π * x * u) := by
  have hint := integrable_exp_smul hc hi u
  rw [FT, Real.fourier_real_eq_integral_exp_smul, ← Complex.reCLM_apply,
    ← ContinuousLinearMap.integral_comp_comm _ hint]
  congr 1
  ext x
  simp only [Complex.reCLM_apply, smul_eq_mul, Complex.mul_re, Complex.exp_ofReal_mul_I_re,
    Complex.ofReal_re, Complex.ofReal_im, mul_zero, sub_zero]
  rw [mul_comm]
  congr 1
  rw [show -2 * π * x * u = -(2 * π * x * u) by ring, Real.cos_neg]

lemma FT_neg' (hc : Continuous f) (hi : Integrable f) (u : ℝ) : FT f (-u) = FT f u := by
  rw [FT_eq_integral_cos' hc hi, FT_eq_integral_cos' hc hi]
  congr 1; funext x
  rw [show 2 * π * x * -u = -(2 * π * x * u) by ring, Real.cos_neg]

lemma FT_nonpos_of_ge' {lam : ℝ} (hf : AdmissibleDecay lam f) {u : ℝ} (hu : lam ≤ u) :
    FT f u ≤ 0 := by
  have hc : Continuous (FT f) := continuous_FT hf.integrable
  have hcl : IsClosed {u | FT f u ≤ 0} := isClosed_le hc continuous_const
  have hsub : Set.Ioi lam ⊆ {u | FT f u ≤ 0} := fun u hu' =>
    hf.fourier_nonpos u (lt_of_lt_of_le hu' (le_abs_self u))
  have := closure_minimal hsub hcl
  rw [closure_Ioi] at this
  exact this hu

lemma tendsto_div_cocompact {h : ℝ} (hh : 0 < h) :
    Tendsto (fun x : ℝ => x / h) (cocompact ℝ) (cocompact ℝ) := by
  rw [cocompact_eq_atBot_atTop]
  exact ((tendsto_id.atBot_div_const hh).mono_right le_sup_left).sup
    ((tendsto_id.atTop_div_const hh).mono_right le_sup_right)

lemma tendsto_mul_cocompact {h : ℝ} (hh : 0 < h) :
    Tendsto (fun x : ℝ => h * x) (cocompact ℝ) (cocompact ℝ) := by
  rw [cocompact_eq_atBot_atTop]
  exact ((tendsto_id.const_mul_atBot hh).mono_right le_sup_left).sup
    ((tendsto_id.const_mul_atTop hh).mono_right le_sup_right)

lemma rpow_scale {c : ℝ} (hc : 0 < c) :
    (fun x : ℝ => |x / c| ^ (-(2:ℝ))) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
  have e : (fun x : ℝ => |x / c| ^ (-(2:ℝ))) = fun x => c ^ 2 * |x| ^ (-(2:ℝ)) := by
    funext x
    rw [abs_div, abs_of_pos hc, Real.div_rpow (abs_nonneg x) hc.le, Real.rpow_neg hc.le,
      Real.rpow_two, div_inv_eq_mul, mul_comm]
  rw [e]
  exact isBigO_const_mul_self _ _ _

lemma rpow_scale' {c : ℝ} (hc : 0 < c) :
    (fun x : ℝ => |c * x| ^ (-(2:ℝ))) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
  have e : (fun x : ℝ => |c * x| ^ (-(2:ℝ))) = fun x => (c ^ 2)⁻¹ * |x| ^ (-(2:ℝ)) := by
    funext x
    rw [abs_mul, abs_of_pos hc, Real.mul_rpow hc.le (abs_nonneg x), Real.rpow_neg hc.le,
      Real.rpow_two]
  rw [e]
  exact isBigO_const_mul_self _ _ _

lemma poisson' (hc : Continuous f) (hfO : f =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))))
    (hFO : 𝓕 (fun x => (f x : ℂ)) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))))
    {h : ℝ} (hh : 0 < h) :
    Summable (fun m : ℤ => f (m / h)) ∧ Summable (fun k : ℤ => FT f (k * h)) ∧
    ∑' m : ℤ, f (m / h) = h * ∑' k : ℤ, FT f (k * h) := by
  set g : ℝ → ℂ := fun x => (f x : ℂ) with hgdef
  set G : ℝ → ℂ := fun x => g (x / h) with hGdef
  have hgc : Continuous g := Complex.continuous_ofReal.comp hc
  have hGc : Continuous G := hgc.comp (continuous_id.div_const h)
  have hFG : ∀ ξ : ℝ, 𝓕 G ξ = (h : ℂ) * 𝓕 g (h * ξ) := fun ξ => fourier_comp_div g hh ξ
  have hgO : g =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
    refine (IsBigO.of_bound 1 (Eventually.of_forall fun x => ?_)).trans hfO
    simp [hgdef, Complex.norm_real]
  have hGO : G =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) :=
    (hgO.comp_tendsto (tendsto_div_cocompact hh)).trans (rpow_scale hh)
  have hFGO : (𝓕 G) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
    have h1 : (𝓕 G) = fun ξ => (h : ℂ) * 𝓕 g (h * ξ) := funext hFG
    rw [h1]
    exact ((hFO.comp_tendsto (tendsto_mul_cocompact hh)).const_mul_left (h : ℂ)).trans
      (rpow_scale' hh)
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

lemma norm_S_le {h : ℝ} {η : ℕ → ℝ} (hη : Summable η) (α : ℝ) :
    ‖∑' n : ℕ, η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))‖ ≤ ∑' n : ℕ, |η (n + 1)| := by
  have hle : ∀ n : ℕ, ‖η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))‖ ≤ |η (n + 1)| := by
    intro n; rw [Real.norm_eq_abs, abs_mul]
    exact mul_le_of_le_one_right (abs_nonneg _) (Real.abs_cos_le_one _)
  have hs : Summable (fun n : ℕ => ‖η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))‖) :=
    Summable.of_nonneg_of_le (fun n => norm_nonneg _) hle (summable_shift hη)
  exact (norm_tsum_le_tsum_norm hs).trans (hs.tsum_le_tsum hle (summable_shift hη))

lemma integrable_mul_P' {h : ℝ} {η : ℕ → ℝ} (hc : Continuous f) (hi : Integrable f)
    (hη : Summable η) : Integrable (fun α => f α * P h η α) := by
  have hPc : Continuous (P h η) := continuous_P hη
  have hb : ∀ α, ‖P h η α‖ ≤ |1 - h| + |2 * h| * ∑' n : ℕ, |η (n + 1)| := by
    intro α
    unfold P
    refine (norm_add_le _ _).trans (add_le_add (le_of_eq (Real.norm_eq_abs _)) ?_)
    rw [norm_mul, Real.norm_eq_abs]
    exact mul_le_mul_of_nonneg_left (norm_S_le hη α) (abs_nonneg _)
  refine (hi.norm.mul_const (|1 - h| + |2 * h| * ∑' n : ℕ, |η (n + 1)|)).mono'
    (hc.mul hPc).aestronglyMeasurable (Eventually.of_forall fun α => ?_)
  rw [norm_mul]
  exact mul_le_mul_of_nonneg_left (hb α) (norm_nonneg _)

lemma integral_mul_P' {h : ℝ} {η : ℕ → ℝ} (hc : Continuous f) (hi : Integrable f)
    (hη : Summable η) :
    ∫ α, f α * P h η α =
      (1 - h) * FT f 0 + 2 * h * ∑' n : ℕ, η (n + 1) * FT f ((n + 1) * h) := by
  set S : ℝ → ℝ := fun α => ∑' n : ℕ, η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h)) with hS
  have hSc : Continuous S := continuous_S hη
  have hfS : Integrable (fun α => f α * S α) := by
    refine (hi.norm.mul_const (∑' n : ℕ, |η (n + 1)|)).mono' (hc.mul hSc).aestronglyMeasurable
      (Eventually.of_forall fun α => ?_)
    rw [norm_mul]; exact mul_le_mul_of_nonneg_left (norm_S_le hη α) (norm_nonneg _)
  have e1 : (fun α => f α * P h η α) = fun α => (1 - h) * f α + 2 * h * (f α * S α) := by
    funext α; simp only [P, hS]; ring
  rw [e1, integral_add (hi.const_mul _) (hfS.const_mul _), integral_const_mul, integral_const_mul]
  have e0 : ∫ α, f α = FT f 0 := by
    rw [FT_eq_integral_cos' hc hi]; simp
  rw [e0]
  congr 2
  have hFi : ∀ n : ℕ, Integrable
      (fun α => f α * (η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h)))) := by
    intro n
    have hm : Continuous (fun α => f α * (η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h)))) := by
      fun_prop
    refine (hi.norm.mul_const |η (n + 1)|).mono' hm.aestronglyMeasurable
      (Eventually.of_forall fun α => ?_)
    rw [norm_mul, norm_mul, Real.norm_eq_abs (η (n + 1)), Real.norm_eq_abs (Real.cos _)]
    exact mul_le_mul_of_nonneg_left (mul_le_of_le_one_right (abs_nonneg _)
      (Real.abs_cos_le_one _)) (norm_nonneg _)
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
  rw [FT_eq_integral_cos' hc hi, ← integral_const_mul]
  congr 1; funext α; ring

/-- Lemma 1 for decaying kernels (Schwartz f included). -/
theorem weak_duality_decay {h lam θ : ℝ} {k0 : ℕ} {η : ℕ → ℝ}
    (hh : 0 < h) (hlam : (k0 : ℝ) * h = lam) (hη : Summable η)
    (hcore : ∀ k : ℕ, 1 ≤ k → k < k0 → η k = -1) (hge : ∀ k : ℕ, k0 ≤ k → -1 ≤ η k)
    (hθ0 : 0 ≤ θ) (hθ : θ ≤ 1 / (2 * h))
    (hP : ∀ α, 0 ≤ α → α ≤ 1 / (2 * h) → (if α ≤ θ then α else 0) ≤ P h η α)
    {f : ℝ → ℝ} (hf : AdmissibleDecay lam f) : J θ f ≤ 1 := by
  obtain ⟨hP0, hPθ⟩ := P_ge_of_base hh hθ hP
  have hcf : Continuous f := hf.continuous
  have hPc : Continuous (P h η) := continuous_P hη
  have hint : Integrable (fun α => f α * P h η α) := integrable_mul_P' hcf hf.integrable hη
  have hnn : ∀ α, 0 ≤ f α * P h η α := fun α => mul_nonneg (hf.nonneg α) (hP0 α)
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
  obtain ⟨hs1, hs2, hpois⟩ := poisson' hcf hf.f_decay hf.fourier_decay hh
  have hB : f 0 ≤ ∑' m : ℤ, f (m / h) := by
    have := hs1.le_tsum 0 (fun j _ => hf.nonneg _)
    simpa using this
  have hFTeven : ∀ k : ℤ, FT f ((-k : ℤ) * h) = FT f (k * h) := by
    intro k
    rw [Int.cast_neg, neg_mul, FT_neg' hcf hf.integrable]
  have hfold : ∑' k : ℤ, FT f (k * h) = FT f 0 + 2 * ∑' n : ℕ, FT f ((n + 1) * h) := by
    rw [tsum_int_fold (a := fun k : ℤ => FT f (k * h)) hFTeven hs2]
    push_cast
    simp
  have hC := integral_mul_P' (h := h) hcf hf.integrable hη
  have hF0 : FT f 0 = 1 := by simp [FT, hf.fourier_zero]
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
      · apply FT_nonpos_of_ge' hf
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

/-- Every Admissible kernel is AdmissibleDecay. -/
lemma Admissible.toDecay {lam : ℝ} (hf : Admissible lam f) : AdmissibleDecay lam f where
  continuous := hf.contDiff.continuous
  integrable := hf.contDiff.continuous.integrable_of_hasCompactSupport hf.hasCompactSupport
  f_decay := by
    have hev : f =ᶠ[cocompact ℝ] 0 := by
      filter_upwards [hf.hasCompactSupport.isCompact.compl_mem_cocompact] with x hx
      exact image_eq_zero_of_notMem_tsupport hx
    exact hev.isBigO.trans (isBigO_zero _ _)
  fourier_decay := by
    obtain ⟨C, hC⟩ := decay hf.contDiff hf.hasCompactSupport
    refine IsBigO.of_bound (C / (4 * π ^ 2)) ?_
    have hmem : ({0}ᶜ : Set ℝ) ∈ cocompact ℝ := isCompact_singleton.compl_mem_cocompact
    filter_upwards [hmem] with x hx
    have hx0 : x ≠ 0 := hx
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity : (0:ℝ) ≤ |x| ^ (-(2:ℝ))),
      Real.rpow_neg (abs_nonneg x), Real.rpow_two, sq_abs]
    have hxp : 0 < |x| := abs_pos.mpr hx0
    have hpos : 0 < (2 * π * |x|) ^ 2 := by positivity
    have h3 : ‖𝓕 (fun x => (f x : ℂ)) x‖ ≤ C / (2 * π * |x|) ^ 2 := by
      rw [le_div_iff₀ hpos]; linarith [hC x]
    calc ‖𝓕 (fun x => (f x : ℂ)) x‖ ≤ C / (2 * π * |x|) ^ 2 := h3
      _ = C / (4 * π ^ 2) * (x ^ 2)⁻¹ := by
        rw [mul_pow, sq_abs]
        have := Real.pi_pos.ne'
        field_simp
        ring
  even := hf.even
  nonneg := hf.nonneg
  fourier_zero := hf.fourier_zero
  fourier_nonpos := hf.fourier_nonpos

/-- A real-valued Schwartz kernel (a complex Schwartz function g with real values f) with
Lamzouri's sign conditions is AdmissibleDecay. -/
lemma AdmissibleDecay.of_schwartz {lam : ℝ} (g : 𝓢(ℝ, ℂ)) (hfg : ∀ x, g x = (f x : ℂ))
    (heven : ∀ x, f (-x) = f x) (hnn : ∀ x, 0 ≤ f x) (h0 : 𝓕 (fun x => (f x : ℂ)) 0 = 1)
    (hneg : ∀ u, lam < |u| → (𝓕 (fun x => (f x : ℂ)) u).re ≤ 0) : AdmissibleDecay lam f := by
  have hgf : (fun x => (f x : ℂ)) = ⇑g := funext fun x => (hfg x).symm
  have hfre : f = fun x => (g x).re := funext fun x => by rw [hfg]; simp
  have hO : f =O[cocompact ℝ] (⇑g) :=
    IsBigO.of_bound 1 (Eventually.of_forall fun x => by rw [hfg, Complex.norm_real, one_mul])
  refine ⟨?_, ?_, hO.trans (g.isBigO_cocompact_rpow (-2)), ?_, heven, hnn, h0, hneg⟩
  · rw [hfre]; exact Complex.continuous_re.comp g.continuous
  · rw [hfre]; exact g.integrable.re
  · rw [hgf, ← SchwartzMap.fourier_coe]; exact (𝓕 g).isBigO_cocompact_rpow (-2)

/-- Theorem B for decaying kernels (SmallGap's Lemma T class), from CertB. -/
theorem theoremB_decay (hcert : CertB) {θ lam : ℝ} (hθ0 : 0 ≤ θ) (hθ : θ ≤ 4 / 7)
    (hlam : lam ≤ 31 / 40) {f : ℝ → ℝ} (hf : AdmissibleDecay lam f) : J θ f ≤ 1 := by
  have hf' : AdmissibleDecay (31 / 40) f := hf.mono hlam
  have h1 : J (4 / 7) f ≤ 1 := by
    refine weak_duality_decay (h := 31 / 280) (k0 := 7) (η := etaOf certB) (by norm_num)
      (by norm_num) (summable_etaOf _) (etaOf_core certB_core) (fun k _ => etaOf_ge certB_ge k)
      (by norm_num) (by norm_num) (fun α h0 h1 => ?_) hf'
    rw [P_etaOf]
    exact hcert α h0 (by linarith [show (1 : ℝ) / (2 * (31 / 280)) = 140 / 31 by norm_num])
  exact (J_mono hf'.nonneg hf'.continuous hθ0 hθ).trans h1

/-- Theorem A for decaying kernels, unconditional. -/
theorem theoremA_decay {θ lam : ℝ} (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hlam : lam ≤ 1 / (1 + θ))
    {f : ℝ → ℝ} (hf : AdmissibleDecay lam f) : J θ f ≤ 1 := by
  have ha0 : 0 < 1 / (1 + θ) := by positivity
  have hf' : AdmissibleDecay (1 / (1 + θ)) f := hf.mono hlam
  have hhalf : 1 / (2 * (1 / (1 + θ))) = (1 + θ) / 2 := by field_simp
  refine weak_duality_decay (h := 1 / (1 + θ)) (k0 := 1) (η := etaA (1 / (1 + θ))) ha0 (by simp)
    (summable_etaA ha0) (fun k hk1 hk => by omega) (fun k _ => etaA_ge _ k) hθ0.le
    (by rw [hhalf]; linarith) (fun α h0 h1 => ?_) hf'
  rw [hhalf] at h1
  rw [P_etaA (triangleSeries ha0), tri_sum_base hθ0 h0 h1]
  exact SA_ge hθ0 hθ1 h0 h1

/-- Theorem B (ii) for decaying kernels (θ = 1, RH case), from CertB1. -/
theorem theoremB1_decay (hcert : CertB1) {θ lam : ℝ} (hθ0 : 0 ≤ θ) (hθ : θ ≤ 1)
    (hlam : lam ≤ 119 / 200) {f : ℝ → ℝ} (hf : AdmissibleDecay lam f) : J θ f ≤ 1 := by
  have hf' : AdmissibleDecay (119 / 200) f := hf.mono hlam
  have h1 : J 1 f ≤ 1 := by
    refine weak_duality_decay (h := 119 / 600) (k0 := 3) (η := etaOf certB1) (by norm_num)
      (by norm_num) (summable_etaOf _) (etaOf_core certB1_core) (fun k _ => etaOf_ge certB1_ge k)
      (by norm_num) (by norm_num) (fun α h0 h1 => ?_) hf'
    rw [P_etaOf]
    exact hcert α h0 (by linarith [show (1 : ℝ) / (2 * (119 / 600)) = 300 / 119 by norm_num])
  exact (J_mono hf'.nonneg hf'.continuous hθ0 hθ).trans h1

end ExploreDPP


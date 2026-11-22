import Riemann.Deep.LeanNBRate.Mellin

/-!
# From the Mellin side to Nyman's constrained distance
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem integrableOn_Ioi_iff_exp_neg (g : ℝ → ℝ) :
    IntegrableOn g (Ioi 0) ↔ Integrable (fun v : ℝ => Real.exp (-v) * g (Real.exp (-v))) := by
  rw [← integrable_comp_exp g]
  constructor
  · intro h
    have := h.comp_neg
    simpa [smul_eq_mul] using this
  · intro h
    have := h.comp_neg
    simpa [smul_eq_mul] using this

theorem integrableOn_nbF_sq {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ) (hθ : ∀ k, θ k ∈ Ioc (0 : ℝ) 1) :
    IntegrableOn (fun x => ‖nbF c θ x‖ ^ 2) (Ioi 0) := by
  rw [integrableOn_Ioi_iff_exp_neg]
  have h := (memLp_two_iff_integrable_sq_norm (integrable_nbF_exp c θ hθ).aestronglyMeasurable).mp
    (memLp_nbF_exp c θ hθ)
  refine h.congr (Filter.Eventually.of_forall fun v => ?_)
  simp only [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (Real.exp_pos _), mul_pow]
  have hx : Real.exp (-v / 2) ^ 2 = Real.exp (-v) := by rw [sq, ← Real.exp_add]; ring_nf
  rw [hx]

theorem memLp_nbFun {θ : ℝ} (hθ : θ ∈ Ioc (0 : ℝ) 1) : MemLp (nbFun θ) 2 volPos := by
  have hm : AEStronglyMeasurable (nbFun θ) volPos :=
    (measurable_fract.comp (measurable_const.div measurable_id)).aestronglyMeasurable
  rw [memLp_two_iff_integrable_sq_norm hm]
  have hmeas : AEStronglyMeasurable (fun x => ‖nbFun θ x‖ ^ 2) volume :=
    ((measurable_fract.comp (measurable_const.div measurable_id)).norm.pow_const 2).aestronglyMeasurable
  have hsplit : Ioi (0 : ℝ) = Ioc 0 1 ∪ Ioi 1 := (Ioc_union_Ioi_eq_Ioi zero_le_one).symm
  show IntegrableOn (fun x => ‖nbFun θ x‖ ^ 2) (Ioi 0)
  rw [hsplit]
  refine IntegrableOn.union ?_ ?_
  · refine Measure.integrableOn_of_bounded (M := 1) (by simp) hmeas ?_
    refine Filter.Eventually.of_forall fun x => ?_
    rw [norm_pow, Real.norm_eq_abs, Real.norm_eq_abs, abs_abs, sq_abs]
    unfold nbFun
    have h0 := Int.fract_nonneg (θ / x)
    have h1 := Int.fract_lt_one (θ / x)
    nlinarith
  · have hc : IntegrableOn (fun x : ℝ => x ^ (-2 : ℝ)) (Ioi 1) := integrableOn_Ioi_rpow_of_lt (by norm_num) one_pos
    refine hc.mono' hmeas.restrict ?_
    rw [ae_restrict_iff' measurableSet_Ioi]
    refine Filter.Eventually.of_forall fun x hx => ?_
    have hx1 : (1 : ℝ) < x := hx
    have hx0 : 0 < x := by linarith
    unfold nbFun
    rw [fract_div_eq hθ.1 (by linarith [hθ.2]), norm_pow, Real.norm_eq_abs, Real.norm_eq_abs, abs_abs,
      sq_abs]
    rw [Real.rpow_neg hx0.le, Real.rpow_two]
    rw [div_pow, ← one_div]
    apply div_le_div_of_nonneg_right _ (by positivity)
    nlinarith [hθ.1, hθ.2]

theorem norm_sq_L2pos (u : L2pos) : ‖u‖ ^ 2 = ∫ x, (u x) ^ 2 ∂volPos := by
  rw [← real_inner_self_eq_norm_sq, MeasureTheory.L2.inner_def]
  congr 1; ext x
  rw [real_inner_self_eq_norm_sq, Real.norm_eq_abs, sq_abs]

/-- **Nyman's distance from a complex combination** (real parts; the constraint on the real part). -/
theorem nymanDist_sq_le {lam : ℝ} (hlam : 0 < lam) {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ)
    (hθ : ∀ k, θ k ∈ Icc lam 1) (hc : (∑ k, c k * (θ k : ℂ)).re = 0) :
    nymanDist lam ^ 2 ≤ ∫ x in Ioi (0 : ℝ), ‖nbF c θ x‖ ^ 2 := by
  have hθ' : ∀ k, θ k ∈ Ioc (0 : ℝ) 1 := fun k => ⟨by linarith [(hθ k).1], (hθ k).2⟩
  set fr : ℝ → ℝ := fun x => ∑ k, (c k).re * nbFun (θ k) x
  have hmem : MemLp fr 2 volPos :=
    memLp_finsetSum _ fun k _ => (memLp_nbFun (hθ' k)).const_mul _
  set f : L2pos := hmem.toLp fr
  have hf : f ∈ nymanSet lam := by
    refine ⟨n, fun k => (c k).re, θ, hθ, ?_, MemLp.coeFn_toLp hmem⟩
    have : (∑ k, c k * (θ k : ℂ)).re = ∑ k, (c k).re * θ k := by
      rw [Complex.re_sum]; refine Finset.sum_congr rfl fun k _ => ?_
      simp [Complex.mul_re]
    rw [← this]; exact hc
  have h1 : nymanDist lam ≤ ‖chiVec - f‖ := by
    rw [← dist_eq_norm]; exact Metric.infDist_le_dist_of_mem hf
  have h0 : 0 ≤ nymanDist lam := Metric.infDist_nonneg
  have h2 : ‖chiVec - f‖ ^ 2 = ∫ x in Ioi (0 : ℝ), ((Ioo (0 : ℝ) 1).indicator (fun _ => (1 : ℝ)) x - fr x) ^ 2 := by
    rw [norm_sq_L2pos]
    refine integral_congr_ae ?_
    have e1 := Lp.coeFn_sub chiVec f
    have e2 : (chiVec : ℝ → ℝ) =ᵐ[volPos] (Ioo (0 : ℝ) 1).indicator (fun _ => (1 : ℝ)) :=
      indicatorConstLp_coeFn
    have e3 : (f : ℝ → ℝ) =ᵐ[volPos] fr := MemLp.coeFn_toLp hmem
    filter_upwards [e1, e2, e3] with x hx1 hx2 hx3
    rw [hx1, Pi.sub_apply, hx2, hx3]
  have h3 : ∫ x in Ioi (0 : ℝ), ((Ioo (0 : ℝ) 1).indicator (fun _ => (1 : ℝ)) x - fr x) ^ 2 ≤
      ∫ x in Ioi (0 : ℝ), ‖nbF c θ x‖ ^ 2 := by
    refine integral_mono_ae ?_ (integrableOn_nbF_sq c θ hθ') ?_
    · -- integrability of the left side from the L² class
      have := (memLp_two_iff_integrable_sq_norm (Lp.aestronglyMeasurable (chiVec - f))).mp
        (Lp.memLp (chiVec - f))
      refine this.congr ?_
      have e1 := Lp.coeFn_sub chiVec f
      have e2 : (chiVec : ℝ → ℝ) =ᵐ[volPos] (Ioo (0 : ℝ) 1).indicator (fun _ => (1 : ℝ)) :=
        indicatorConstLp_coeFn
      have e3 : (f : ℝ → ℝ) =ᵐ[volPos] fr := MemLp.coeFn_toLp hmem
      filter_upwards [e1, e2, e3] with x hx1 hx2 hx3
      rw [hx1, Pi.sub_apply, hx2, hx3, Real.norm_eq_abs, sq_abs]
    · have hne : ∀ᵐ x ∂(volume.restrict (Ioi (0 : ℝ))), x ≠ 1 := by
        refine ae_restrict_of_ae ?_
        rw [ae_iff]; simp [measure_singleton]
      filter_upwards [hne] with x hx
      have hre : (Ioo (0 : ℝ) 1).indicator (fun _ => (1 : ℝ)) x - fr x = (nbF c θ x).re := by
        unfold nbF
        rw [Complex.sub_re, Complex.re_sum]
        have hind : ((Ioc (0 : ℝ) 1).indicator (fun _ => (1 : ℂ)) x).re =
            (Ioo (0 : ℝ) 1).indicator (fun _ => (1 : ℝ)) x := by
          by_cases h : x ∈ Ioo (0 : ℝ) 1
          · rw [indicator_of_mem h, indicator_of_mem (Ioo_subset_Ioc_self h)]; simp
          · rw [indicator_of_notMem h]
            have : x ∉ Ioc (0 : ℝ) 1 := fun h' => h ⟨h'.1, lt_of_le_of_ne h'.2 hx⟩
            rw [indicator_of_notMem this]; simp
        rw [hind]
        congr 1
        refine Finset.sum_congr rfl fun k _ => ?_
        simp [Complex.mul_re, nbFun]
      rw [hre]
      have := Complex.abs_re_le_norm (nbF c θ x)
      have h00 : 0 ≤ |(nbF c θ x).re| := abs_nonneg _
      calc (nbF c θ x).re ^ 2 = |(nbF c θ x).re| ^ 2 := (sq_abs _).symm
        _ ≤ ‖nbF c θ x‖ ^ 2 := pow_le_pow_left₀ h00 this 2
  calc nymanDist lam ^ 2 ≤ ‖chiVec - f‖ ^ 2 := pow_le_pow_left₀ h0 h1 2
    _ = _ := h2
    _ ≤ _ := h3

/-- **Nyman's distance and the critical line.** -/
theorem nymanDist_sq_le_line (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) {lam : ℝ} (hlam : 0 < lam)
    {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ) (hθ : ∀ k, θ k ∈ Icc lam 1) (hc : (∑ k, c k * (θ k : ℂ)).re = 0) :
    nymanDist lam ^ 2 ≤ (2 * Real.pi)⁻¹ * ∫ t : ℝ,
      ‖(1 + riemannZeta (1 / 2 + t * I) * ∑ k, c k * (θ k : ℂ) ^ (1 / 2 + t * I)) / (1 / 2 + t * I)‖ ^ 2 := by
  have hθ' : ∀ k, θ k ∈ Ioc (0 : ℝ) 1 := fun k => ⟨by linarith [(hθ k).1], (hθ k).2⟩
  rw [← integral_nbF_sq hP h215 c θ hθ']
  exact nymanDist_sq_le hlam c θ hθ hc

end Riemann.Deep.LeanNBRate


import Riemann.Deep.LeanNBRate.Contour2

/-!
# U1 contour argument, part 3: the residue at s (Cauchy on the right half-plane) and the strip [2,3]
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

/-- (G2') uniform decay from a C/(1+|y|)^2 bound for |y| >= Y0. -/
theorem decay_of_le_inv_sq' {f : ℂ → ℂ} {a b C Y₀ : ℝ}
    (h : ∀ z : ℂ, a ≤ z.re → z.re ≤ b → Y₀ ≤ |z.im| → ‖f z‖ ≤ C * ((1 + |z.im|) ^ 2)⁻¹) :
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
  filter_upwards [hC, Filter.eventually_ge_atTop (max 1 Y₀)] with T hT hTY x hx
  have hT1 : 1 ≤ T := le_of_max_le_left hTY
  have hTY' : Y₀ ≤ T := le_of_max_le_right hTY
  have hT0 : |T| = T := abs_of_pos (by linarith)
  constructor
  · have e : ((x : ℂ) + T * I).im = T := by simp
    have := h (x + T * I) (by simpa using hx.1) (by simpa using hx.2) (by rw [e, hT0]; exact hTY')
    rw [e, hT0] at this
    exact this.trans hT
  · have e : ((x : ℂ) + ((-T : ℝ) : ℂ) * I).im = -T := by simp
    have := h (x + ((-T : ℝ) : ℂ) * I) (by simpa using hx.1) (by simpa using hx.2)
      (by rw [e, abs_neg, hT0]; exact hTY')
    rw [e, abs_neg, hT0] at this
    exact this.trans hT

theorem continuousAt_Psi0 (hH : HBound H M β D N) {s z : ℂ} (hz : 2 ≤ z.re) (hzs : z ≠ s) :
    ContinuousAt (Psi0 H s) z :=
  (continuousAt_Fq hH (ne_one_of_two_le hz) (riemannZeta_ne_zero_of_one_lt_re (by linarith))).div
    (continuousAt_const.sub continuousAt_id) (sub_ne_zero.mpr (Ne.symm hzs))

theorem differentiableAt_Psi0 (hH : HBound H M β D N) {s z : ℂ} (hz : 2 ≤ z.re) (hzs : z ≠ s) :
    DifferentiableAt ℂ (Psi0 H s) z :=
  (differentiableAt_Fq hH (ne_one_of_two_le hz) (riemannZeta_ne_zero_of_one_lt_re (by linarith))).div
    ((differentiableAt_const _).sub differentiableAt_id) (sub_ne_zero.mpr (Ne.symm hzs))

theorem norm_Psi0_le (hH : HBound H M β D N) {s z : ℂ} (hz : 2 ≤ z.re) {δ : ℝ} (hδ : 0 < δ)
    (hd : δ ≤ ‖s - z‖) : ‖Psi0 H s z‖ ≤ 3 * M * (D / min D 1) ^ 2 / δ * ((1 + |z.im|) ^ 2)⁻¹ := by
  rw [Psi0, norm_div]
  have h := norm_Fq_right_le hH hz
  have hM := hH.M_nonneg
  calc ‖Fq H z‖ / ‖s - z‖ ≤ (3 * M * (D / min D 1) ^ 2 * ((1 + |z.im|) ^ 2)⁻¹) / δ :=
        div_le_div₀ (by positivity) h hδ hd
    _ = _ := by ring

theorem integrable_Psi0_line (hH : HBound H M β D N) {s : ℂ} {σ : ℝ} (hσ : 2 ≤ σ) (hσs : σ ≠ s.re) :
    Integrable fun y : ℝ => Psi0 H s (σ + y * I) := by
  have hδ : 0 < |s.re - σ| := abs_pos.mpr (sub_ne_zero.mpr (Ne.symm hσs))
  have hz : ∀ y : ℝ, 2 ≤ ((σ : ℂ) + y * I).re := fun y => by simp; exact hσ
  have hzs : ∀ y : ℝ, ((σ : ℂ) + y * I) ≠ s := fun y h => hσs (by rw [← h]; simp)
  refine integrable_of_le_inv_sq (C := 3 * M * (D / min D 1) ^ 2 / |s.re - σ|) ?_ (fun y => ?_)
  · refine (continuous_iff_continuousAt.mpr fun y => ?_).aestronglyMeasurable
    exact ContinuousAt.comp (f := fun y : ℝ => (σ : ℂ) + y * I) (x := y)
      (continuousAt_Psi0 hH (hz y) (hzs y)) (by fun_prop : Continuous fun y : ℝ => (σ : ℂ) + y * I).continuousAt
  · have hd : |s.re - σ| ≤ ‖s - ((σ : ℂ) + y * I)‖ := by
      have := Complex.abs_re_le_norm (s - ((σ : ℂ) + y * I)); simpa using this
    have := norm_Psi0_le hH (hz y) hδ hd
    have him : ((σ : ℂ) + y * I).im = y := by simp
    rwa [him] at this

/-- **Cauchy integral formula on the right half-plane**: int Psi0(3+iy) dy = 2 pi F(s) for Re s > 3. -/
theorem integral_Psi0_line (hH : HBound H M β D N) {s : ℂ} (hs : 3 < s.re) :
    ∫ y : ℝ, Psi0 H s (3 + y * I) = 2 * π * Fq H s := by
  have hs2 : 2 ≤ s.re := by linarith
  have hFd : DifferentiableAt ℂ (Fq H) s :=
    differentiableAt_Fq hH (ne_one_of_two_le hs2) (riemannZeta_ne_zero_of_one_lt_re (by linarith))
  -- strip [3, c'] for every c' > Re s
  have hstrip : ∀ c' : ℝ, s.re < c' → I * VI (Psi0 H s) c' - I * VI (Psi0 H s) 3 =
      2 * π * I * (-Fq H s) := by
    intro c' hc'
    have h := strip_residues (f := Psi0 H s) (a := 3) (b := c') (by linarith) {s} (fun _ => -Fq H s)
      (by intro p hp; rw [Finset.mem_singleton] at hp; subst hp; exact ⟨hs, hc'⟩)
      (fun x hx1 _ hxP => continuousAt_Psi0 hH (by linarith) (by simpa using hxP))
      (fun x hx1 _ hxP => differentiableAt_Psi0 hH (by linarith) (by simpa using hxP))
      ?_ (integrable_Psi0_line hH (by norm_num) (by linarith)) (integrable_Psi0_line hH (by linarith)
        (by linarith)) ?_
    · rw [Finset.sum_singleton] at h; exact h
    · intro p hp
      rw [Finset.mem_singleton] at hp; subst hp
      refine ⟨fun z => -dslope (Fq H) p z, (continuousAt_dslope_same.mpr hFd).neg, ?_⟩
      filter_upwards [self_mem_nhdsWithin] with z hz
      have hz' : z - p ≠ 0 := sub_ne_zero.mpr hz
      rw [dslope_of_ne _ hz, slope_def_field, Psi0, show p - z = -(z - p) by ring]
      field_simp
      ring
    · refine decay_of_le_inv_sq' (Y₀ := |s.im| + 1) (C := 3 * M * (D / min D 1) ^ 2 / 1) fun z hz1 hz2 hzY => ?_
      refine norm_Psi0_le hH (by linarith) one_pos ?_
      have := Complex.abs_im_le_norm (s - z)
      simp only [sub_im] at this
      have h2 : |z.im| ≤ |s.im| + |s.im - z.im| := by
        calc |z.im| = |s.im - (s.im - z.im)| := by ring_nf
          _ ≤ |s.im| + |s.im - z.im| := abs_sub _ _
      linarith
  -- c' -> infinity
  have hlim : Tendsto (fun c' : ℝ => VI (Psi0 H s) c') atTop (𝓝 0) := by
    have hK : 0 ≤ 3 * M * (D / min D 1) ^ 2 := by
      have := hH.M_nonneg; have := hH.D_pos; positivity
    have h := tendsto_integral_filter_of_dominated_convergence (l := (atTop : Filter ℝ))
      (F := fun c' : ℝ => fun y : ℝ => Psi0 H s (c' + y * I)) (f := fun _ => (0 : ℂ))
      (fun y => 3 * M * (D / min D 1) ^ 2 / 1 * ((1 + |y|) ^ 2)⁻¹) ?_ ?_
      (integrable_inv_one_add_abs_sq.const_mul _) ?_
    · simpa [VI] using h
    · filter_upwards [Filter.eventually_ge_atTop (s.re + 1)] with c' hc'
      exact (integrable_Psi0_line hH (by linarith) (by linarith)).aestronglyMeasurable
    · filter_upwards [Filter.eventually_ge_atTop (s.re + 1)] with c' hc'
      refine Filter.Eventually.of_forall fun y => ?_
      have hd : 1 ≤ ‖s - ((c' : ℂ) + y * I)‖ := by
        have := Complex.abs_re_le_norm (s - ((c' : ℂ) + y * I))
        simp at this
        have : |s.re - c'| = c' - s.re := by rw [abs_sub_comm]; exact abs_of_nonneg (by linarith)
        linarith
      have := norm_Psi0_le hH (z := (c' : ℂ) + y * I) (by simp; linarith) one_pos hd
      simpa using this
    · refine Filter.Eventually.of_forall fun y => ?_
      rw [tendsto_zero_iff_norm_tendsto_zero]
      have hb : ∀ c' : ℝ, s.re + 1 ≤ c' → ‖Psi0 H s (c' + y * I)‖ ≤ 3 * M * (D / (D + (c' - 1 / 2))) ^ N := by
        intro c' hc'
        have hz : 2 ≤ ((c' : ℂ) + y * I).re := by simp; linarith
        have hd : 1 ≤ ‖s - ((c' : ℂ) + y * I)‖ := by
          have := Complex.abs_re_le_norm (s - ((c' : ℂ) + y * I))
          simp at this
          have : |s.re - c'| = c' - s.re := by rw [abs_sub_comm]; exact abs_of_nonneg (by linarith)
          linarith
        rw [Psi0, norm_div]
        refine (div_le_self (norm_nonneg _) hd).trans ((norm_Fq_right hH hz).trans ?_)
        have hM := hH.M_nonneg; have hD := hH.D_pos
        have hc0 : 0 < D + (c' - 1 / 2) := by linarith
        gcongr
        have := Complex.abs_re_le_norm ((c' : ℂ) + y * I - 1 / 2)
        simp at this
        have h' : c' - 1 / 2 ≤ |c' - 1 / 2| := le_abs_self _
        linarith
      have hD := hH.D_pos
      have ht : Tendsto (fun c' : ℝ => 3 * M * (D / (D + (c' - 1 / 2))) ^ N) atTop (𝓝 0) := by
        have h1 : Tendsto (fun c' : ℝ => D + (c' - 1 / 2)) atTop atTop := by
          have : Tendsto (fun c' : ℝ => c' + (D - 1 / 2)) atTop atTop :=
            tendsto_atTop_add_const_right _ _ tendsto_id
          exact this.congr fun c' => by ring
        have h2 : Tendsto (fun c' : ℝ => D * (D + (c' - 1 / 2))⁻¹) atTop (𝓝 (D * 0)) :=
          h1.inv_tendsto_atTop.const_mul D
        have h3 := (h2.pow N).const_mul (3 * M)
        have hN : N ≠ 0 := by have := hH.N_ge; omega
        simp only [mul_zero, zero_pow hN] at h3
        refine h3.congr fun c' => ?_
        rw [← div_eq_mul_inv]
      refine squeeze_zero' (Filter.Eventually.of_forall fun _ => norm_nonneg _) ?_ ht
      filter_upwards [Filter.eventually_ge_atTop (s.re + 1)] with c' hc'
      exact hb c' hc'
  have hev : ∀ᶠ c' : ℝ in atTop, VI (Psi0 H s) c' + 2 * π * Fq H s = VI (Psi0 H s) 3 := by
    filter_upwards [Filter.eventually_gt_atTop s.re] with c' hc'
    have h := hstrip c' hc'
    have hI : (I : ℂ) ≠ 0 := Complex.I_ne_zero
    have : I * (VI (Psi0 H s) c' + 2 * π * Fq H s - VI (Psi0 H s) 3) = 0 := by
      rw [mul_sub, mul_add]; rw [show I * (2 * π * Fq H s) = -(2 * π * I * (-Fq H s)) by ring, ← h]; ring
    have := (mul_eq_zero.mp this).resolve_left hI
    exact sub_eq_zero.mp this
  have hl2 : Tendsto (fun c' : ℝ => VI (Psi0 H s) c' + 2 * π * Fq H s) atTop (𝓝 (2 * π * Fq H s)) := by
    simpa using hlim.add_const (2 * π * Fq H s)
  have hfin := tendsto_nhds_unique (tendsto_const_nhds.congr' (hev.mono fun c' h => h.symm)) hl2
  have e : VI (Psi0 H s) 3 = ∫ y : ℝ, Psi0 H s (3 + y * I) := by simp [VI]
  rw [← e]
  exact hfin

end Riemann.Deep.LeanNBRate

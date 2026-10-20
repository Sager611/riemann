import Mathlib

/-!
# Fayed, arXiv:2209.01890v55 — Appendix, Theorem 1 — REFUTED

Theorem 1 (L'Hôpital rule for a disconnected domain). F, G differentiable on ℝ⁺;
D = ⋃ₘ (a_M, b_M), b_M < a_{M+1}, b_M → ∞; along D:
(i) F → 0, G → 0; (ii) F'/G' → L (finite or infinite);
(iii) |F(u)/G(u) - F(b_M⁻)/G(b_M⁻)| ≤ δ_M → 0 on each D_M.
Claimed conclusion: F/G → L along D.

FALSE. Counterexample on the paper's own buffered domain
D_γ = ⋃_M (M + 1/2 + γ, M + 1 - γ), 0 < γ < 1/4:
F u = cos (2πu) / u  (first Fourier mode of the paper's own F ≈ (1-2σ₀) B̄₂(u)/(2u)),
G u = u^(-1/2)       (the paper's G = u^(-p), p = 1/2).
Then F/G = cos(2πu)/√u → 0, but F'/G' = 2(2πu sin 2πu + cos 2πu)/√u → -∞ on D_γ
(there sin 2πu ≤ -4γ). The paper's proof applies L'Hôpital at the finite endpoint b_M,
where F, G are not 0/0; across the gaps F undoes what it gains on D_γ.

Consequence: step (55) of the paper fails; (39) F/G → 0 and (50) F'/G' → ±∞ coexist.
-/

open Filter Topology Set Real

namespace Fayed

/-- u → ∞ inside D. -/
def along (D : Set ℝ) : Filter ℝ := atTop ⊓ 𝓟 D

/-- D = ⋃ₘ (a M, b M). -/
def unionIoo (a b : ℕ → ℝ) : Set ℝ := ⋃ M, Ioo (a M) (b M)

/-- Admissible limits in Theorem 1: finite or ±∞. -/
def IsLimitTarget (l : Filter ℝ) : Prop := l = atTop ∨ l = atBot ∨ ∃ L : ℝ, l = 𝓝 L

/-- Theorem 1, real-valued special case. -/
def Theorem1 : Prop :=
  ∀ (F G : ℝ → ℝ) (a b δ : ℕ → ℝ) (l : Filter ℝ),
    IsLimitTarget l →
    (∀ u, 0 < u → DifferentiableAt ℝ F u) →
    (∀ u, 0 < u → DifferentiableAt ℝ G u) →
    0 < a 0 → (∀ M, a M < b M) → (∀ M, b M < a (M + 1)) → Tendsto b atTop atTop →
    Tendsto F (along (unionIoo a b)) (𝓝 0) →
    Tendsto G (along (unionIoo a b)) (𝓝 0) →
    Tendsto (fun u => deriv F u / deriv G u) (along (unionIoo a b)) l →
    (∀ M, ∀ u ∈ Ioo (a M) (b M), |F u / G u - F (b M) / G (b M)| ≤ δ M) →
    Tendsto δ atTop (𝓝 0) →
    Tendsto (fun u => F u / G u) (along (unionIoo a b)) l

/-- Theorem 1 as printed: complex-valued F, G; L finite or ∞ (cobounded ℂ). -/
def Theorem1C : Prop :=
  ∀ (F G : ℝ → ℂ) (a b δ : ℕ → ℝ) (l : Filter ℂ),
    (l = Bornology.cobounded ℂ ∨ ∃ L : ℂ, l = 𝓝 L) →
    (∀ u, 0 < u → DifferentiableAt ℝ F u) →
    (∀ u, 0 < u → DifferentiableAt ℝ G u) →
    0 < a 0 → (∀ M, a M < b M) → (∀ M, b M < a (M + 1)) → Tendsto b atTop atTop →
    Tendsto F (along (unionIoo a b)) (𝓝 0) →
    Tendsto G (along (unionIoo a b)) (𝓝 0) →
    Tendsto (fun u => deriv F u / deriv G u) (along (unionIoo a b)) l →
    (∀ M, ∀ u ∈ Ioo (a M) (b M), ‖F u / G u - F (b M) / G (b M)‖ ≤ δ M) →
    Tendsto δ atTop (𝓝 0) →
    Tendsto (fun u => F u / G u) (along (unionIoo a b)) l

/-! ## The buffered domain D_γ (paper eq. (10)) -/

/-- Left ends a_M = M + 1/2 + γ (indices shifted by M₀). -/
noncomputable def aγ (γ : ℝ) (M₀ : ℕ) (M : ℕ) : ℝ := ((M₀ + M : ℕ) : ℝ) + 1 / 2 + γ
/-- Right ends b_M = M + 1 - γ. -/
noncomputable def bγ (γ : ℝ) (M₀ : ℕ) (M : ℕ) : ℝ := ((M₀ + M : ℕ) : ℝ) + 1 - γ
/-- D_γ = ⋃_{M ≥ M₀} (M + 1/2 + γ, M + 1 - γ). -/
abbrev Dγ (γ : ℝ) (M₀ : ℕ) : Set ℝ := unionIoo (aγ γ M₀) (bγ γ M₀)

lemma along_le_atTop (D : Set ℝ) : along D ≤ atTop := inf_le_left

lemma eventually_mem_along (D : Set ℝ) : ∀ᶠ u in along D, u ∈ D :=
  mem_inf_of_right (mem_principal_self D)

lemma along_Dγ_neBot {γ : ℝ} (hγ' : γ < 1 / 4) (M₀ : ℕ) : (along (Dγ γ M₀)).NeBot := by
  rw [along, inf_principal_neBot_iff]
  intro U hU
  obtain ⟨x, hx⟩ := mem_atTop_sets.1 hU
  obtain ⟨n, hn⟩ := exists_nat_ge x
  refine ⟨((M₀ + n : ℕ) : ℝ) + 3 / 4, hx _ ?_, mem_iUnion.2 ⟨n, ?_, ?_⟩⟩
  · push_cast
    have : (0 : ℝ) ≤ M₀ := Nat.cast_nonneg M₀
    linarith
  · simp only [aγ]; linarith
  · simp only [bγ]; linarith

lemma tendsto_aγ (γ : ℝ) (M₀ : ℕ) : Tendsto (aγ γ M₀) atTop atTop := by
  have h : Tendsto (fun M : ℕ => (M : ℝ) + ((M₀ : ℝ) + 1 / 2 + γ)) atTop atTop :=
    tendsto_atTop_add_const_right _ _ tendsto_natCast_atTop_atTop
  refine h.congr (fun M => ?_)
  simp only [aγ]; push_cast; ring

lemma tendsto_bγ (γ : ℝ) (M₀ : ℕ) : Tendsto (bγ γ M₀) atTop atTop := by
  have h : Tendsto (fun M : ℕ => (M : ℝ) + ((M₀ : ℝ) + 1 - γ)) atTop atTop :=
    tendsto_atTop_add_const_right _ _ tendsto_natCast_atTop_atTop
  refine h.congr (fun M => ?_)
  simp only [bγ]; push_cast; ring

/-- On D_γ: sin (2πu) ≤ -4γ. -/
lemma sin_le_on_buffer {γ u : ℝ} {n : ℕ} (hγ : 0 < γ) (_hγ' : γ ≤ 1 / 4)
    (hu : u ∈ Ioo ((n : ℝ) + 1 / 2 + γ) ((n : ℝ) + 1 - γ)) : sin (2 * π * u) ≤ -(4 * γ) := by
  have hπ := Real.pi_pos
  set θ := 2 * π * (u - n - 1 / 2) with hθ
  have h1 : 2 * π * γ < θ := by
    have : γ < u - n - 1 / 2 := by linarith [hu.1]
    have := mul_lt_mul_of_pos_left this (by positivity : (0 : ℝ) < 2 * π)
    linarith
  have h2 : θ < π - 2 * π * γ := by
    have : u - n - 1 / 2 < 1 / 2 - γ := by linarith [hu.2]
    have := mul_lt_mul_of_pos_left this (by positivity : (0 : ℝ) < 2 * π)
    linarith
  have hsin : sin (2 * π * u) = -sin θ := by
    have : 2 * π * u = θ + π + (n : ℕ) * (2 * π) := by rw [hθ]; ring
    rw [this, Real.sin_add_nat_mul_two_pi, Real.sin_add_pi]
  have hγ0 : 4 * γ = 2 / π * (2 * π * γ) := by field_simp; ring
  have hθ0 : 0 < 2 * π * γ := by positivity
  have key : 4 * γ ≤ sin θ := by
    rcases le_total θ (π / 2) with h | h
    · have hJ := Real.mul_le_sin (by nlinarith) h
      have := mul_le_mul_of_nonneg_left h1.le (by positivity : (0 : ℝ) ≤ 2 / π)
      linarith
    · have hJ := Real.mul_le_sin (x := π - θ) (by linarith) (by linarith)
      rw [Real.sin_pi_sub] at hJ
      have h3 : 2 * π * γ ≤ π - θ := by linarith
      have := mul_le_mul_of_nonneg_left h3 (by positivity : (0 : ℝ) ≤ 2 / π)
      linarith
  linarith

/-! ## The counterexample -/

/-- F u = cos (2πu) / u. -/
noncomputable def Fc (u : ℝ) : ℝ := cos (2 * π * u) / u
/-- G u = u^(-1/2). -/
noncomputable def Gc (u : ℝ) : ℝ := (√u)⁻¹

lemma hasDerivAt_Fc {u : ℝ} (hu : 0 < u) :
    HasDerivAt Fc ((-sin (2 * π * u) * (2 * π * 1) * u - cos (2 * π * u) * 1) / u ^ 2) u := by
  have h1 : HasDerivAt (fun x : ℝ => 2 * π * x) (2 * π * 1) u :=
    (hasDerivAt_id' u).const_mul (2 * π)
  exact h1.cos.div (hasDerivAt_id' u) hu.ne'

lemma hasDerivAt_Gc {u : ℝ} (hu : 0 < u) :
    HasDerivAt Gc (-(1 / (2 * √u)) / (√u) ^ 2) u :=
  (Real.hasDerivAt_sqrt hu.ne').inv (Real.sqrt_ne_zero'.2 hu)

lemma ratio_eq {u : ℝ} (hu : 0 < u) :
    deriv Fc u / deriv Gc u = 2 * (2 * π * u * sin (2 * π * u) + cos (2 * π * u)) / √u := by
  rw [(hasDerivAt_Fc hu).deriv, (hasDerivAt_Gc hu).deriv]
  obtain ⟨r, hr, rfl⟩ : ∃ r, 0 < r ∧ u = r ^ 2 :=
    ⟨√u, Real.sqrt_pos.2 hu, (Real.sq_sqrt hu.le).symm⟩
  rw [Real.sqrt_sq hr.le]
  field_simp
  ring

lemma abs_FG_le {u : ℝ} (hu : 0 < u) : |Fc u / Gc u| ≤ (√u)⁻¹ := by
  obtain ⟨r, hr, rfl⟩ : ∃ r, 0 < r ∧ u = r ^ 2 :=
    ⟨√u, Real.sqrt_pos.2 hu, (Real.sq_sqrt hu.le).symm⟩
  rw [Fc, Gc, Real.sqrt_sq hr.le]
  have hc := abs_cos_le_one (2 * π * r ^ 2)
  have : cos (2 * π * r ^ 2) / r ^ 2 / r⁻¹ = cos (2 * π * r ^ 2) / r := by
    field_simp
  rw [this, abs_div, abs_of_pos hr]
  calc |cos (2 * π * r ^ 2)| / r ≤ 1 / r := by gcongr
    _ = r⁻¹ := one_div r

lemma tendsto_Fc : Tendsto Fc atTop (𝓝 0) := by
  refine squeeze_zero_norm' ?_ tendsto_inv_atTop_zero
  filter_upwards [eventually_gt_atTop 0] with u hu
  have hc := abs_cos_le_one (2 * π * u)
  rw [Fc, norm_div, Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos hu]
  calc |cos (2 * π * u)| / u ≤ 1 / u := by gcongr
    _ = u⁻¹ := one_div u

lemma tendsto_Gc : Tendsto Gc atTop (𝓝 0) := Real.tendsto_sqrt_atTop.inv_tendsto_atTop

lemma tendsto_FG : Tendsto (fun u => Fc u / Gc u) atTop (𝓝 0) := by
  refine squeeze_zero_norm' ?_ tendsto_Gc
  filter_upwards [eventually_gt_atTop 0] with u hu
  rw [Real.norm_eq_abs]
  exact abs_FG_le hu

/-- F'/G' ≤ 2 - 16πγ√u on D_γ ∩ [1, ∞). -/
lemma ratio_le {γ : ℝ} (hγ : 0 < γ) (hγ' : γ ≤ 1 / 4) {M₀ : ℕ} {u : ℝ} (hu : u ∈ Dγ γ M₀) (h1 : 1 ≤ u) :
    deriv Fc u / deriv Gc u ≤ 2 - 16 * π * γ * √u := by
  have hu0 : 0 < u := by linarith
  obtain ⟨M, hM⟩ := mem_iUnion.1 hu
  have hS : sin (2 * π * u) ≤ -(4 * γ) :=
    sin_le_on_buffer (n := M₀ + M) hγ hγ' (by simpa [aγ, bγ] using hM)
  have hC := cos_le_one (2 * π * u)
  rw [ratio_eq hu0]
  have hr0 : 0 < √u := Real.sqrt_pos.2 hu0
  have hr1 : 1 ≤ √u := Real.one_le_sqrt.2 h1
  have hsq : √u * √u = u := Real.mul_self_sqrt hu0.le
  rw [div_le_iff₀ hr0]
  have hπ := Real.pi_pos
  have hmul := mul_le_mul_of_nonneg_left hS (by positivity : (0 : ℝ) ≤ 2 * π * u)
  have h16 : 16 * π * γ * (√u * √u) = 16 * π * γ * u := by rw [hsq]
  nlinarith

lemma tendsto_ratio {γ : ℝ} (hγ : 0 < γ) (hγ' : γ ≤ 1 / 4) (M₀ : ℕ) :
    Tendsto (fun u => deriv Fc u / deriv Gc u) (along (Dγ γ M₀)) atBot := by
  have hev : ∀ᶠ u in along (Dγ γ M₀), deriv Fc u / deriv Gc u ≤ 2 - 16 * π * γ * √u := by
    filter_upwards [(eventually_ge_atTop 1).filter_mono (along_le_atTop _),
      eventually_mem_along _] with u hu1 huD
    exact ratio_le hγ hγ' huD hu1
  have hlim : Tendsto (fun u => 2 - 16 * π * γ * √u) atTop atBot := by
    have h := Real.tendsto_sqrt_atTop.const_mul_atTop (by positivity : (0 : ℝ) < 16 * π * γ)
    have h' := tendsto_atBot_add_const_left atTop 2 (tendsto_neg_atTop_atBot.comp h)
    refine h'.congr (fun u => ?_)
    simp only [Function.comp]
    ring
  exact tendsto_atBot_mono' _ hev (hlim.mono_left (along_le_atTop _))

lemma osc_bound {γ : ℝ} {M₀ : ℕ} (hγ : 0 < γ) (M : ℕ) {u : ℝ}
    (hu : u ∈ Ioo (aγ γ M₀ M) (bγ γ M₀ M)) :
    |Fc u / Gc u - Fc (bγ γ M₀ M) / Gc (bγ γ M₀ M)| ≤ 2 * (√(aγ γ M₀ M))⁻¹ := by
  have ha : 0 < aγ γ M₀ M := by
    simp only [aγ]; have : (0 : ℝ) ≤ ((M₀ + M : ℕ) : ℝ) := Nat.cast_nonneg _; linarith
  have hu0 : 0 < u := ha.trans hu.1
  have hb0 : 0 < bγ γ M₀ M := hu0.trans hu.2
  have hpos : 0 < √(aγ γ M₀ M) := Real.sqrt_pos.2 ha
  have e1 : (√u)⁻¹ ≤ (√(aγ γ M₀ M))⁻¹ := by
    have := Real.sqrt_le_sqrt hu.1.le
    gcongr
  have e2 : (√(bγ γ M₀ M))⁻¹ ≤ (√(aγ γ M₀ M))⁻¹ := by
    have := Real.sqrt_le_sqrt (hu.1.trans hu.2).le
    gcongr
  calc |Fc u / Gc u - Fc (bγ γ M₀ M) / Gc (bγ γ M₀ M)|
      ≤ |Fc u / Gc u| + |Fc (bγ γ M₀ M) / Gc (bγ γ M₀ M)| := abs_sub _ _
    _ ≤ (√u)⁻¹ + (√(bγ γ M₀ M))⁻¹ := add_le_add (abs_FG_le hu0) (abs_FG_le hb0)
    _ ≤ 2 * (√(aγ γ M₀ M))⁻¹ := by linarith

lemma tendsto_δ (γ : ℝ) (M₀ : ℕ) :
    Tendsto (fun M => 2 * (√(aγ γ M₀ M))⁻¹) atTop (𝓝 0) := by
  have h := ((Real.tendsto_sqrt_atTop.comp (tendsto_aγ γ M₀)).inv_tendsto_atTop).const_mul 2
  simpa [Function.comp_def] using h

/-- The counterexample satisfies every hypothesis of Theorem 1, yet F/G → 0 ≠ -∞. -/
theorem counterexample {γ : ℝ} (hγ : 0 < γ) (hγ' : γ < 1 / 4) (M₀ : ℕ) :
    (∀ u, 0 < u → DifferentiableAt ℝ Fc u) ∧ (∀ u, 0 < u → DifferentiableAt ℝ Gc u) ∧
    0 < aγ γ M₀ 0 ∧ (∀ M, aγ γ M₀ M < bγ γ M₀ M) ∧ (∀ M, bγ γ M₀ M < aγ γ M₀ (M + 1)) ∧
    Tendsto (bγ γ M₀) atTop atTop ∧
    Tendsto Fc (along (Dγ γ M₀)) (𝓝 0) ∧ Tendsto Gc (along (Dγ γ M₀)) (𝓝 0) ∧
    Tendsto (fun u => deriv Fc u / deriv Gc u) (along (Dγ γ M₀)) atBot ∧
    (∀ M, ∀ u ∈ Ioo (aγ γ M₀ M) (bγ γ M₀ M),
      |Fc u / Gc u - Fc (bγ γ M₀ M) / Gc (bγ γ M₀ M)| ≤ 2 * (√(aγ γ M₀ M))⁻¹) ∧
    Tendsto (fun M => 2 * (√(aγ γ M₀ M))⁻¹) atTop (𝓝 0) ∧
    Tendsto (fun u => Fc u / Gc u) (along (Dγ γ M₀)) (𝓝 0) ∧
    ¬ Tendsto (fun u => Fc u / Gc u) (along (Dγ γ M₀)) atBot := by
  have := along_Dγ_neBot hγ' M₀
  have hFG := tendsto_FG.mono_left (along_le_atTop (Dγ γ M₀))
  refine ⟨fun u hu => (hasDerivAt_Fc hu).differentiableAt,
    fun u hu => (hasDerivAt_Gc hu).differentiableAt, ?_, ?_, ?_, tendsto_bγ γ M₀,
    tendsto_Fc.mono_left (along_le_atTop _), tendsto_Gc.mono_left (along_le_atTop _),
    tendsto_ratio hγ hγ'.le M₀, fun M u hu => osc_bound hγ M hu, tendsto_δ γ M₀, hFG,
    hFG.not_tendsto (disjoint_nhds_atBot 0)⟩
  · simp only [aγ]; push_cast; linarith
  · intro M; simp only [aγ, bγ]; linarith
  · intro M; simp only [aγ, bγ]; push_cast; linarith

/-- Theorem 1 of the paper is false (real-valued case). -/
theorem theorem1_false : ¬ Theorem1 := by
  intro h
  obtain ⟨hF, hG, ha0, hab, hba, hb, hF0, hG0, hR, hosc, hδ, -, hnot⟩ :=
    counterexample (γ := 1 / 8) (by norm_num) (by norm_num) 0
  exact hnot (h Fc Gc _ _ _ atBot (Or.inr (Or.inl rfl)) hF hG ha0 hab hba hb hF0 hG0 hR hosc hδ)

/-- Theorem 1 as printed (complex-valued, L = ∞) is false. -/
theorem theorem1C_false : ¬ Theorem1C := by
  intro h
  obtain ⟨hF, hG, ha0, hab, hba, hb, hF0, hG0, hR, hosc, hδ, hFG, -⟩ :=
    counterexample (γ := 1 / 8) (by norm_num) (by norm_num) 0
  have hne : (along (unionIoo (aγ (1 / 8) 0) (bγ (1 / 8) 0))).NeBot :=
    along_Dγ_neBot (γ := 1 / 8) (by norm_num) 0
  have dF : ∀ u, 0 < u → HasDerivAt (fun y => (Fc y : ℂ)) ((deriv Fc u : ℝ) : ℂ) u :=
    fun u hu => (hF u hu).hasDerivAt.ofReal_comp
  have dG : ∀ u, 0 < u → HasDerivAt (fun y => (Gc y : ℂ)) ((deriv Gc u : ℝ) : ℂ) u :=
    fun u hu => (hG u hu).hasDerivAt.ofReal_comp
  have hpos : ∀ᶠ u in along (Dγ (1 / 8) 0), 0 < u :=
    (eventually_gt_atTop 0).filter_mono (along_le_atTop _)
  have hR' : Tendsto (fun u => deriv (fun y => (Fc y : ℂ)) u / deriv (fun y => (Gc y : ℂ)) u)
      (along (Dγ (1 / 8) 0)) (Bornology.cobounded ℂ) := by
    rw [← tendsto_norm_atTop_iff_cobounded]
    refine (tendsto_abs_atBot_atTop.comp hR).congr' ?_
    filter_upwards [hpos] with u hu
    simp [(dF u hu).deriv, (dG u hu).deriv, Function.comp, abs_div]
  have hF0' : Tendsto (fun y => (Fc y : ℂ)) (along (Dγ (1 / 8) 0)) (𝓝 0) := by
    have := (Complex.continuous_ofReal.tendsto 0).comp hF0
    rwa [Complex.ofReal_zero] at this
  have hG0' : Tendsto (fun y => (Gc y : ℂ)) (along (Dγ (1 / 8) 0)) (𝓝 0) := by
    have := (Complex.continuous_ofReal.tendsto 0).comp hG0
    rwa [Complex.ofReal_zero] at this
  have hconc := h (fun y => (Fc y : ℂ)) (fun y => (Gc y : ℂ)) _ _ _ (Bornology.cobounded ℂ)
    (Or.inl rfl) (fun u hu => (dF u hu).differentiableAt) (fun u hu => (dG u hu).differentiableAt)
    ha0 hab hba hb hF0' hG0' hR'
    (fun M u hu => by
      have := hosc M u hu
      rw [← Complex.ofReal_div, ← Complex.ofReal_div, ← Complex.ofReal_sub, Complex.norm_real,
        Real.norm_eq_abs]
      exact this) hδ
  rw [← tendsto_norm_atTop_iff_cobounded] at hconc
  have h0 : Tendsto (fun u => ‖(Fc u : ℂ) / (Gc u : ℂ)‖)
      (along (unionIoo (aγ (1 / 8) 0) (bγ (1 / 8) 0))) (𝓝 0) := by
    have := hFG.norm
    simpa [Complex.norm_real] using this
  exact not_tendsto_nhds_of_tendsto_atTop hconc 0 h0

end Fayed

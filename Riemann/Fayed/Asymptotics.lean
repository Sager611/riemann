import Riemann.Fayed.Defs

/-!
# Fayed (39),(42),(50),(54): the asymptotic lemmas, proved from Inputs

Route: perron33 + em17 give u^s J(u,s) = -u^s EMcorr(u,s) + O(1/u) at a zero s.
F(u) = u^{s₀}J(u,s₀) - u^{1-s̄₀}J(u,1-s̄₀) (combine the integrands), so the B̄₁ terms cancel
and F = O(1/u). F'(u) = u⁻¹ (s₀u^{s₀}J(u,s₀) - (1-s̄₀)u^{1-s̄₀}J(u,1-s̄₀)) + boundary, and
the B̄₁ terms now add up to (1-2σ₀)B̄₁(u); the boundary terms are O(u^{-2}) by convexity.
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real ComplexConjugate

namespace Fayed

noncomputable section

/-! ### Bernoulli functions -/

lemma Bbar_bounded (k : ℕ) : ∃ c : ℝ, 0 ≤ c ∧ ∀ x, |Bbar k x| ≤ c := by
  have hc : Continuous fun x : ℝ => Polynomial.aeval x (Polynomial.bernoulli k) :=
    Polynomial.continuous_aeval _
  obtain ⟨c, hc'⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    hc.continuousOn
  refine ⟨max c 0, le_max_right _ _, fun x => ?_⟩
  have hx : Int.fract x ∈ Icc (0 : ℝ) 1 := ⟨Int.fract_nonneg x, (Int.fract_lt_one x).le⟩
  have h := hc' _ hx
  rw [Real.norm_eq_abs] at h
  exact h.trans (le_max_left _ _)

lemma Bbar_one (x : ℝ) : Bbar 1 x = Int.fract x - 1 / 2 := by
  simp [Bbar, Polynomial.bernoulli_one]

/-- On D_γ, B̄₁(u) = {u} - 1/2 > γ. -/
lemma Bbar_one_gt {γ : ℝ} (hγ : 0 < γ) {M₀ : ℕ} {u : ℝ} (hu : u ∈ Dγ γ M₀) : γ < Bbar 1 u := by
  obtain ⟨M, hM⟩ := mem_iUnion.1 hu
  simp only [aγ, bγ, mem_Ioo] at hM
  have hf : Int.fract u = u - ((M₀ + M : ℕ) : ℝ) := by
    rw [Int.fract_eq_iff]
    exact ⟨by linarith, by linarith, ((M₀ + M : ℕ) : ℤ), by push_cast; ring⟩
  rw [Bbar_one, hf]
  linarith

lemma aγ_pos {γ : ℝ} (hγ : 0 < γ) (M₀ M : ℕ) : 0 < aγ γ M₀ M := by
  simp only [aγ]
  have := Nat.cast_nonneg (α := ℝ) (M₀ + M)
  linarith

/-! ### Small algebra -/

lemma inv_div_rpow {u : ℝ} (hu : 0 < u) (a : ℝ) : u⁻¹ / u ^ a = u ^ (-1 - a) := by
  rw [Real.rpow_sub hu, Real.rpow_neg_one]

lemma norm_FG (F : ℂ) {p u : ℝ} (hu : 0 < u) : ‖F / Gfun p u‖ = ‖F‖ / u ^ (-p) := by
  rw [Gfun, norm_div, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos (Real.rpow_pos_of_pos hu _)]

lemma cpow_mul_EMcorr {u : ℝ} (hu : 0 < u) (s : ℂ) :
    (u : ℂ) ^ s * EMcorr u s = (Bbar 1 u : ℂ) + s * (Bbar 2 u : ℂ) / (2 * (u : ℂ)) +
      s * (s + 1) * (Bbar 3 u : ℂ) / (6 * (u : ℂ) ^ 2) := by
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  have h0 : (u : ℂ) ^ s ≠ 0 := by
    rw [Ne, Complex.cpow_eq_zero_iff]; exact fun h => hu' h.1
  have e1 : (u : ℂ) ^ (-s) = ((u : ℂ) ^ s)⁻¹ := Complex.cpow_neg _ _
  have e2 : (u : ℂ) ^ (s + 1) = (u : ℂ) ^ s * u := by
    rw [Complex.cpow_add _ _ hu', Complex.cpow_one]
  have e2' : (u : ℂ) ^ (2 : ℂ) = (u : ℂ) ^ 2 := by exact_mod_cast Complex.cpow_natCast (u : ℂ) 2
  have e3 : (u : ℂ) ^ (s + 2) = (u : ℂ) ^ s * (u : ℂ) ^ 2 := by
    rw [Complex.cpow_add _ _ hu', e2']
  rw [EMcorr, e1, e2, e3]
  field_simp

/-- ‖α b₂/(2u) + β b₃/(6u²)‖ ≤ (‖α‖+‖β‖)c/u. -/
lemma norm_tail_le {α β : ℂ} {b2 b3 c u : ℝ} (hu : 1 ≤ u) (h2 : |b2| ≤ c) (h3 : |b3| ≤ c) :
    ‖α * (b2 : ℂ) / (2 * (u : ℂ)) + β * (b3 : ℂ) / (6 * (u : ℂ) ^ 2)‖ ≤ (‖α‖ + ‖β‖) * c / u := by
  have hu0 : 0 < u := by linarith
  have hc : 0 ≤ c := (abs_nonneg _).trans h2
  have hα := norm_nonneg α
  have hβ := norm_nonneg β
  have n1 : ‖α * (b2 : ℂ) / (2 * (u : ℂ))‖ = ‖α‖ * |b2| / (2 * u) := by
    rw [norm_div, norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos hu0]
    norm_num
  have n2 : ‖β * (b3 : ℂ) / (6 * (u : ℂ) ^ 2)‖ = ‖β‖ * |b3| / (6 * u ^ 2) := by
    rw [norm_div, norm_mul, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs,
      Complex.norm_real, Real.norm_eq_abs, abs_of_pos hu0]
    norm_num
  refine (norm_add_le _ _).trans ?_
  rw [n1, n2]
  have t1 : ‖α‖ * |b2| / (2 * u) ≤ ‖α‖ * c / u := by
    rw [div_le_div_iff₀ (by linarith) hu0]
    have := mul_le_mul_of_nonneg_left h2 hα
    nlinarith [mul_nonneg hα hc, abs_nonneg b2, mul_nonneg hα (abs_nonneg b2)]
  have t2 : ‖β‖ * |b3| / (6 * u ^ 2) ≤ ‖β‖ * c / u := by
    rw [div_le_div_iff₀ (by positivity) hu0]
    have := mul_le_mul_of_nonneg_left h3 hβ
    have h6 : u ≤ 6 * u ^ 2 := by nlinarith
    nlinarith [mul_nonneg hβ hc, mul_nonneg hβ (abs_nonneg b3),
      mul_le_mul_of_nonneg_left h6 (mul_nonneg hβ hc)]
  calc _ ≤ ‖α‖ * c / u + ‖β‖ * c / u := add_le_add t1 t2
    _ = _ := by ring

/-! ### The integrals -/

/-- ζ-part of all integrands: u^{1-η+i(t₀+y)} ζ(1-η+i(t₀+y)). -/
def Kf (η t₀ u y : ℝ) : ℂ :=
  (u : ℂ) ^ (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I)

lemma cpow_mul_Φ {u : ℝ} (hu : 0 < u) (η : ℝ) (s : ℂ) (y : ℝ) :
    (u : ℂ) ^ s * Φ η u s y = Kf η s.im u y / (((1 - s.re - η : ℝ) : ℂ) + y * I) := by
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  have hexp : s + ((((1 - s.re - η : ℝ) : ℂ)) + y * I) = ((1 - η : ℝ) : ℂ) + (s.im + y) * I := by
    apply Complex.ext <;> simp <;> ring
  rw [Φ, Kf, ← hexp, Complex.cpow_add s _ hu']
  ring

lemma cpow_mul_J {u : ℝ} (hu : 0 < u) (η : ℝ) (s : ℂ) :
    (u : ℂ) ^ s * J η u s = ((1 / (2 * π) : ℝ) : ℂ) *
      ∫ y in (-(u ^ 8))..(u ^ 8), Kf η s.im u y / (((1 - s.re - η : ℝ) : ℂ) + y * I) := by
  rw [J, mul_left_comm, ← intervalIntegral.integral_const_mul]
  congr 1
  exact intervalIntegral.integral_congr (fun y _ => cpow_mul_Φ hu η s y)

lemma s0_re (σ₀ t₀ : ℝ) : ((σ₀ : ℂ) + t₀ * I).re = σ₀ := by simp
lemma s0_im (σ₀ t₀ : ℝ) : ((σ₀ : ℂ) + t₀ * I).im = t₀ := by simp
lemma s1_re (σ₀ t₀ : ℝ) : (1 - conj ((σ₀ : ℂ) + t₀ * I)).re = 1 - σ₀ := by simp
lemma s1_im (σ₀ t₀ : ℝ) : (1 - conj ((σ₀ : ℂ) + t₀ * I)).im = t₀ := by simp

lemma cpow_mul_J0 {u : ℝ} (hu : 0 < u) (η σ₀ t₀ : ℝ) :
    (u : ℂ) ^ ((σ₀ : ℂ) + t₀ * I) * J η u ((σ₀ : ℂ) + t₀ * I) = ((1 / (2 * π) : ℝ) : ℂ) *
      ∫ y in (-(u ^ 8))..(u ^ 8), Kf η t₀ u y / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) := by
  rw [cpow_mul_J hu, s0_re, s0_im]

lemma cpow_mul_J1 {u : ℝ} (hu : 0 < u) (η σ₀ t₀ : ℝ) :
    (u : ℂ) ^ (1 - conj ((σ₀ : ℂ) + t₀ * I)) * J η u (1 - conj ((σ₀ : ℂ) + t₀ * I)) = ((1 / (2 * π) : ℝ) : ℂ) *
      ∫ y in (-(u ^ 8))..(u ^ 8), Kf η t₀ u y / (((σ₀ - η : ℝ) : ℂ) + y * I) := by
  rw [cpow_mul_J hu, s1_re, s1_im, show (1 - (1 - σ₀) - η : ℝ) = σ₀ - η by ring]

lemma continuous_Kf_div (η t₀ u d : ℝ) (hu : 0 < u) (hη : 0 < η) (hd : 0 < d) :
    Continuous fun y : ℝ => Kf η t₀ u y / (((d : ℝ) : ℂ) + y * I) := by
  unfold Kf
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  have hg : Continuous fun y : ℝ => ((1 - η : ℝ) : ℂ) + (t₀ + y) * I := by fun_prop
  have hζ : Continuous fun y : ℝ => riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) := by
    refine continuous_iff_continuousAt.2 fun y => ?_
    refine ((differentiableAt_riemannZeta ?_).continuousAt).comp hg.continuousAt
    intro h
    have := congrArg Complex.re h
    simp at this
    linarith
  refine Continuous.div ((hg.const_cpow (Or.inl hu')).mul hζ) (by fun_prop) (fun y => ?_)
  intro h
  have := congrArg Complex.re h
  simp at this
  linarith

/-- Identity (a): F(u) = u^{s₀}J(u,s₀) - u^{1-s̄₀}J(u,1-s̄₀). -/
lemma Ffun_eq {η σ₀ t₀ u : ℝ} (hu : 0 < u) (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) :
    Ffun η σ₀ t₀ u = (u : ℂ) ^ ((σ₀ : ℂ) + t₀ * I) * J η u ((σ₀ : ℂ) + t₀ * I) - (u : ℂ) ^ (1 - conj ((σ₀ : ℂ) + t₀ * I)) * J η u (1 - conj ((σ₀ : ℂ) + t₀ * I)) := by
  rw [cpow_mul_J0 hu, cpow_mul_J1 hu, ← mul_sub, ← intervalIntegral.integral_sub
    ((continuous_Kf_div η t₀ u (1 - σ₀ - η) hu hη (by linarith)).intervalIntegrable _ _)
    ((continuous_Kf_div η t₀ u (σ₀ - η) hu hη (by linarith)).intervalIntegrable _ _), Ffun]
  congr 1
  refine intervalIntegral.integral_congr (fun y _ => ?_)
  simp only [ψ, Kf]
  ring

/-- (44) under the integral: F' = u⁻¹(s₀u^{s₀}J(s₀) - s₁u^{s₁}J(s₁)) + boundary. -/
lemma Fder_eq {η σ₀ t₀ u : ℝ} (hu : 0 < u) (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) :
    Fder η σ₀ t₀ u = (u : ℂ)⁻¹ * (((σ₀ : ℂ) + t₀ * I) * ((u : ℂ) ^ ((σ₀ : ℂ) + t₀ * I) * J η u ((σ₀ : ℂ) + t₀ * I)) -
        (1 - conj ((σ₀ : ℂ) + t₀ * I)) * ((u : ℂ) ^ (1 - conj ((σ₀ : ℂ) + t₀ * I)) * J η u (1 - conj ((σ₀ : ℂ) + t₀ * I)))) +
      ((1 / (2 * π) : ℝ) : ℂ) * (8 * (u : ℂ) ^ 7 * (ψ η σ₀ t₀ u (u ^ 8) + ψ η σ₀ t₀ u (-(u ^ 8)))) := by
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  have hint : ∀ y : ℝ, (1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)) *
      ((((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * (u : ℂ) ^ (((-η : ℝ) : ℂ) + (t₀ + y) * I)) *
      riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) =
      (u : ℂ)⁻¹ * (((σ₀ : ℂ) + t₀ * I) * (Kf η t₀ u y / (((1 - σ₀ - η : ℝ) : ℂ) + y * I)) -
        (1 - conj ((σ₀ : ℂ) + t₀ * I)) * (Kf η t₀ u y / (((σ₀ - η : ℝ) : ℂ) + y * I))) := by
    intro y
    have hc : (u : ℂ) ^ (((-η : ℝ) : ℂ) + (t₀ + y) * I) =
        (u : ℂ)⁻¹ * (u : ℂ) ^ (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) := by
      have he : ((1 - η : ℝ) : ℂ) + (t₀ + y) * I = 1 + (((-η : ℝ) : ℂ) + (t₀ + y) * I) := by
        push_cast; ring
      rw [he, Complex.cpow_add (1 : ℂ) _ hu', Complex.cpow_one, ← mul_assoc, inv_mul_cancel₀ hu',
        one_mul]
    rw [hc]
    have h44 := eq44_algebra (t₀ := t₀) h1 h2 y
    simp only [Kf]
    linear_combination ((u : ℂ)⁻¹ * (u : ℂ) ^ (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) *
      riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I)) * h44
  have hI : (∫ y in (-(u ^ 8))..(u ^ 8), (1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)) *
      ((((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * (u : ℂ) ^ (((-η : ℝ) : ℂ) + (t₀ + y) * I)) *
      riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I)) =
      ∫ y in (-(u ^ 8))..(u ^ 8), (u : ℂ)⁻¹ * (((σ₀ : ℂ) + t₀ * I) * (Kf η t₀ u y / (((1 - σ₀ - η : ℝ) : ℂ) + y * I)) -
        (1 - conj ((σ₀ : ℂ) + t₀ * I)) * (Kf η t₀ u y / (((σ₀ - η : ℝ) : ℂ) + y * I))) :=
    intervalIntegral.integral_congr (fun y _ => hint y)
  have hiA : IntervalIntegrable (fun y : ℝ => ((σ₀ : ℂ) + t₀ * I) * (Kf η t₀ u y / (((1 - σ₀ - η : ℝ) : ℂ) + y * I)))
      volume (-(u ^ 8)) (u ^ 8) :=
    (continuous_const.mul (continuous_Kf_div η t₀ u (1 - σ₀ - η) hu hη (by linarith))).intervalIntegrable _ _
  have hiB : IntervalIntegrable (fun y : ℝ => (1 - conj ((σ₀ : ℂ) + t₀ * I)) * (Kf η t₀ u y / (((σ₀ - η : ℝ) : ℂ) + y * I)))
      volume (-(u ^ 8)) (u ^ 8) :=
    (continuous_const.mul (continuous_Kf_div η t₀ u (σ₀ - η) hu hη (by linarith))).intervalIntegrable _ _
  rw [Fder, hI, intervalIntegral.integral_const_mul, intervalIntegral.integral_sub hiA hiB,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul, cpow_mul_J0 hu,
    cpow_mul_J1 hu]
  ring

/-! ### (34) at a zero: u^s J(u,s) + u^s EMcorr(u,s) = O(1/u) -/

lemma J_at_zero (hin : Inputs) {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) (hz : riemannZeta s = 0)
    {η : ℝ} (hη : 0 < η) (h1 : η < s.re) (h2 : η < 1 - s.re) {γ : ℝ} (hγ : 0 < γ)
    (hγ' : γ < 1 / 4) :
    ∃ C U : ℝ, 1 ≤ U ∧ 0 ≤ C ∧ ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖(u : ℂ) ^ s * J η u s + (u : ℂ) ^ s * EMcorr u s‖ ≤ C / u := by
  obtain ⟨C1, hC1⟩ := hin.em17 s hs0 hs1
  obtain ⟨C2, U, hC2⟩ := hin.perron33 s hs0 hs1 η hη h1 h2 γ hγ hγ' (1 / 2) (by norm_num)
  refine ⟨|C1| + |C2|, max U 1, le_max_right _ _, by positivity, fun M₀ u hu hUu => ?_⟩
  have hu1 : 1 ≤ u := le_trans (le_max_right _ _) hUu
  have hu0 : 0 < u := by linarith
  have e1 := hC1 u hu1
  have e2 := hC2 M₀ u hu (le_trans (le_max_left _ _) hUu)
  have hA := (norm_sub_le _ _).trans (add_le_add e1 e2)
  have hsplit : (u : ℂ) ^ s * J η u s + (u : ℂ) ^ s * EMcorr u s =
      (u : ℂ) ^ s * ((S u s - (riemannZeta s + (u : ℂ) ^ (1 - s) / (1 - s) - EMcorr u s)) -
        (S u s - (J η u s + (u : ℂ) ^ (1 - s) / (1 - s)))) := by
    rw [hz]; ring
  rw [hsplit, norm_mul, Complex.norm_cpow_eq_rpow_re_of_pos hu0, div_eq_mul_inv]
  have hpos : 0 ≤ u ^ s.re := Real.rpow_nonneg hu0.le _
  have ha : 0 ≤ u ^ (-(s.re + 3)) := Real.rpow_nonneg hu0.le _
  have hb : 0 ≤ u ^ (-4 + 1 / 2 : ℝ) := Real.rpow_nonneg hu0.le _
  have hp1 : u ^ s.re * u ^ (-(s.re + 3)) ≤ u⁻¹ := by
    rw [← Real.rpow_add hu0, ← Real.rpow_neg_one]
    exact Real.rpow_le_rpow_of_exponent_le hu1 (by linarith)
  have hp2 : u ^ s.re * u ^ (-4 + 1 / 2 : ℝ) ≤ u⁻¹ := by
    rw [← Real.rpow_add hu0, ← Real.rpow_neg_one]
    exact Real.rpow_le_rpow_of_exponent_le hu1 (by linarith)
  have k1 := mul_le_mul_of_nonneg_left hp1 (abs_nonneg C1)
  have k2 := mul_le_mul_of_nonneg_left hp2 (abs_nonneg C2)
  have k3 : C1 * u ^ (-(s.re + 3)) ≤ |C1| * u ^ (-(s.re + 3)) :=
    mul_le_mul_of_nonneg_right (le_abs_self _) ha
  have k4 : C2 * u ^ (-4 + 1 / 2 : ℝ) ≤ |C2| * u ^ (-4 + 1 / 2 : ℝ) :=
    mul_le_mul_of_nonneg_right (le_abs_self _) hb
  calc _ ≤ u ^ s.re * (C1 * u ^ (-(s.re + 3)) + C2 * u ^ (-4 + 1 / 2 : ℝ)) :=
        mul_le_mul_of_nonneg_left hA hpos
    _ ≤ u ^ s.re * (|C1| * u ^ (-(s.re + 3)) + |C2| * u ^ (-4 + 1 / 2 : ℝ)) :=
        mul_le_mul_of_nonneg_left (add_le_add k3 k4) hpos
    _ = |C1| * (u ^ s.re * u ^ (-(s.re + 3))) + |C2| * (u ^ s.re * u ^ (-4 + 1 / 2 : ℝ)) := by
        ring
    _ ≤ (|C1| + |C2|) * u⁻¹ := by linarith

/-- (34),(35) at s₀ and at 1-s̄₀ (a zero too, by reflection). -/
lemma J_pair (hin : Inputs) (X : Setting) {γ : ℝ} (hγ : 0 < γ) (hγ' : γ < 1 / 4) :
    ∃ C U : ℝ, 1 ≤ U ∧ 0 ≤ C ∧ ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖(u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * J X.η u ((X.σ₀ : ℂ) + X.t₀ * I) + (u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)‖ ≤ C / u ∧ ‖(u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * J X.η u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + (u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))‖ ≤ C / u := by
  obtain ⟨hσa, hσb⟩ := X.hσ₀
  obtain ⟨hη, hη1, hη2⟩ := X.hη
  have z0 := X.hzero
  have z1 := hin.reflect _ z0 (by rw [s0_re]; exact hσa) (by rw [s0_re]; exact hσb)
  obtain ⟨C0, U0, hU0, hC0, h0⟩ := J_at_zero hin (by rw [s0_re]; exact hσa)
    (by rw [s0_re]; exact hσb) z0 hη (by rw [s0_re]; exact hη1) (by rw [s0_re]; exact hη2) hγ hγ'
  obtain ⟨C1, U1, hU1, hC1, h1⟩ := J_at_zero hin (by rw [s1_re]; linarith)
    (by rw [s1_re]; linarith) z1 hη (by rw [s1_re]; exact hη2) (by rw [s1_re]; linarith) hγ hγ'
  refine ⟨C0 + C1, max U0 U1, le_trans hU0 (le_max_left _ _), by linarith,
    fun M₀ u hu hUu => ⟨?_, ?_⟩⟩
  · have hu0 : 0 < u := by linarith [le_max_left U0 U1]
    exact (h0 M₀ u hu (le_trans (le_max_left _ _) hUu)).trans
      (div_le_div_of_nonneg_right (by linarith) hu0.le)
  · have hu0 : 0 < u := by linarith [le_max_left U0 U1]
    exact (h1 M₀ u hu (le_trans (le_max_right _ _) hUu)).trans
      (div_le_div_of_nonneg_right (by linarith) hu0.le)

/-- F = O(1/u) on D_γ (the B̄₁ terms cancel). -/
lemma F_bound (hin : Inputs) (X : Setting) {γ : ℝ} (hγ : 0 < γ) (hγ' : γ < 1 / 4) :
    ∃ K U : ℝ, 1 ≤ U ∧ 0 ≤ K ∧ ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖Ffun X.η X.σ₀ X.t₀ u‖ ≤ K / u := by
  obtain ⟨C, U, hU, hC, hJ⟩ := J_pair hin X hγ hγ'
  obtain ⟨c2, hc2, hb2⟩ := Bbar_bounded 2
  obtain ⟨c3, hc3, hb3⟩ := Bbar_bounded 3
  have hc : 0 ≤ max c2 c3 := le_trans hc2 (le_max_left _ _)
  refine ⟨2 * C + (‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I)‖ + ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1)‖) * max c2 c3,
    U, hU, add_nonneg (by linarith) (mul_nonneg (add_nonneg (norm_nonneg _) (norm_nonneg _)) hc),
    fun M₀ u hu hUu => ?_⟩
  have hu1 : 1 ≤ u := hU.trans hUu
  have hu0 : 0 < u := by linarith
  obtain ⟨j0, j1⟩ := hJ M₀ u hu hUu
  have ht : ‖(u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - (u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)‖ ≤
      (‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I)‖ + ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1)‖) * max c2 c3 / u := by
    have hE : (u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - (u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I) = ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I)) * (Bbar 2 u : ℂ) / (2 * (u : ℂ)) +
        ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1)) * (Bbar 3 u : ℂ) / (6 * (u : ℂ) ^ 2) := by
      rw [cpow_mul_EMcorr hu0, cpow_mul_EMcorr hu0]; ring
    rw [hE]
    exact norm_tail_le hu1 ((hb2 u).trans (le_max_left _ _)) ((hb3 u).trans (le_max_right _ _))
  have key : ∀ a b e : ℂ, ‖a - b + e‖ ≤ ‖a‖ + ‖b‖ + ‖e‖ := fun a b e =>
    (norm_add_le _ _).trans (add_le_add (norm_sub_le _ _) le_rfl)
  have hsplit : Ffun X.η X.σ₀ X.t₀ u = ((u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * J X.η u ((X.σ₀ : ℂ) + X.t₀ * I) + (u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)) - ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * J X.η u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + (u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))) + ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - (u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)) := by
    rw [Ffun_eq hu0 X.hη.1 X.hη.2.1 X.hη.2.2]; ring
  rw [hsplit]
  refine (key _ _ _).trans ?_
  calc _ ≤ C / u + C / u +
        (‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I)‖ + ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1)‖) * max c2 c3 / u :=
        add_le_add (add_le_add j0 j1) ht
    _ = _ := by ring

/-- s₀u^{s₀}J(s₀) - s₁u^{s₁}J(s₁) = (1-2σ₀)B̄₁(u) + O(1/u) on D_γ. -/
lemma X_bound (hin : Inputs) (X : Setting) {γ : ℝ} (hγ : 0 < γ) (hγ' : γ < 1 / 4) :
    ∃ K U : ℝ, 1 ≤ U ∧ 0 ≤ K ∧ ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖(((X.σ₀ : ℂ) + X.t₀ * I) * ((u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * J X.η u ((X.σ₀ : ℂ) + X.t₀ * I)) - (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * J X.η u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)))) - ((1 - 2 * X.σ₀ : ℝ) : ℂ) * (Bbar 1 u : ℂ)‖ ≤ K / u := by
  obtain ⟨C, U, hU, hC, hJ⟩ := J_pair hin X hγ hγ'
  obtain ⟨c2, hc2, hb2⟩ := Bbar_bounded 2
  obtain ⟨c3, hc3, hb3⟩ := Bbar_bounded 3
  have hc : 0 ≤ max c2 c3 := le_trans hc2 (le_max_left _ _)
  refine ⟨‖((X.σ₀ : ℂ) + X.t₀ * I)‖ * C + ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))‖ * C + (‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I) * ((X.σ₀ : ℂ) + X.t₀ * I)‖ +
      ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1)) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1))‖) * max c2 c3, U, hU,
    add_nonneg (add_nonneg (mul_nonneg (norm_nonneg _) hC) (mul_nonneg (norm_nonneg _) hC))
      (mul_nonneg (add_nonneg (norm_nonneg _) (norm_nonneg _)) hc), fun M₀ u hu hUu => ?_⟩
  have hu1 : 1 ≤ u := hU.trans hUu
  have hu0 : 0 < u := by linarith
  obtain ⟨j0, j1⟩ := hJ M₀ u hu hUu
  have hs1 : (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) = ((1 - X.σ₀ : ℝ) : ℂ) + X.t₀ * I := by apply Complex.ext <;> simp
  have ht : ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))) - ((X.σ₀ : ℂ) + X.t₀ * I) * ((u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)) - ((1 - 2 * X.σ₀ : ℝ) : ℂ) * (Bbar 1 u : ℂ)‖ ≤
      (‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I) * ((X.σ₀ : ℂ) + X.t₀ * I)‖ +
        ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1)) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1))‖) * max c2 c3 / u := by
    have hE : (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))) - ((X.σ₀ : ℂ) + X.t₀ * I) * ((u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)) - ((1 - 2 * X.σ₀ : ℝ) : ℂ) * (Bbar 1 u : ℂ) =
        ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I) * ((X.σ₀ : ℂ) + X.t₀ * I)) * (Bbar 2 u : ℂ) / (2 * (u : ℂ)) +
        ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1)) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1))) * (Bbar 3 u : ℂ) /
          (6 * (u : ℂ) ^ 2) := by
      rw [cpow_mul_EMcorr hu0, cpow_mul_EMcorr hu0, hs1]; push_cast; ring
    rw [hE]
    exact norm_tail_le hu1 ((hb2 u).trans (le_max_left _ _)) ((hb3 u).trans (le_max_right _ _))
  have key : ∀ x y a b e : ℂ, ‖x * a - y * b + e‖ ≤ ‖x‖ * ‖a‖ + ‖y‖ * ‖b‖ + ‖e‖ := by
    intro x y a b e
    calc _ ≤ ‖x * a - y * b‖ + ‖e‖ := norm_add_le _ _
      _ ≤ ‖x * a‖ + ‖y * b‖ + ‖e‖ := add_le_add (norm_sub_le _ _) le_rfl
      _ = _ := by rw [norm_mul, norm_mul]
  have hsplit : (((X.σ₀ : ℂ) + X.t₀ * I) * ((u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * J X.η u ((X.σ₀ : ℂ) + X.t₀ * I)) - (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * J X.η u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)))) - ((1 - 2 * X.σ₀ : ℝ) : ℂ) * (Bbar 1 u : ℂ) =
      ((X.σ₀ : ℂ) + X.t₀ * I) * ((u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * J X.η u ((X.σ₀ : ℂ) + X.t₀ * I) + (u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)) - (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * J X.η u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + (u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))) +
        ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((u : ℂ) ^ (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * EMcorr u (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))) - ((X.σ₀ : ℂ) + X.t₀ * I) * ((u : ℂ) ^ ((X.σ₀ : ℂ) + X.t₀ * I) * EMcorr u ((X.σ₀ : ℂ) + X.t₀ * I)) - ((1 - 2 * X.σ₀ : ℝ) : ℂ) * (Bbar 1 u : ℂ)) := by
    ring
  rw [hsplit]
  refine (key _ _ _ _ _).trans ?_
  calc _ ≤ ‖((X.σ₀ : ℂ) + X.t₀ * I)‖ * (C / u) + ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I))‖ * (C / u) +
        (‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * (1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) - ((X.σ₀ : ℂ) + X.t₀ * I) * ((X.σ₀ : ℂ) + X.t₀ * I)‖ +
          ‖(1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) * ((1 - conj ((X.σ₀ : ℂ) + X.t₀ * I)) + 1)) - ((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) * (((X.σ₀ : ℂ) + X.t₀ * I) + 1))‖) * max c2 c3 / u :=
        add_le_add (add_le_add (mul_le_mul_of_nonneg_left j0 (norm_nonneg _))
          (mul_le_mul_of_nonneg_left j1 (norm_nonneg _))) ht
    _ = _ := by ring

/-! ### Boundary terms of (44): O(u^{-2}) by convexity -/

lemma norm_ψ_le (hin : Inputs) {η σ₀ : ℝ} (t₀ : ℝ) (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) :
    ∃ C U : ℝ, 1 ≤ U ∧ 0 ≤ C ∧ ∀ u : ℝ, U ≤ u → ∀ y : ℝ, |y| = u ^ 8 →
      ‖ψ η σ₀ t₀ u y‖ ≤ C / u ^ 11 := by
  obtain ⟨C, hC⟩ := hin.convexity η hη (by linarith) (1 / 4) (by norm_num)
  refine ⟨2 * |C|, |t₀| + 2, by linarith [abs_nonneg t₀], by positivity, fun u hu y hy => ?_⟩
  have hu1 : 1 ≤ u := by linarith [abs_nonneg t₀]
  have hu0 : 0 < u := by linarith
  have hu0' : u ≠ 0 := hu0.ne'
  have hu8 : u ≤ u ^ 8 := le_self_pow₀ hu1 (by norm_num)
  have hT1 : 1 ≤ |t₀ + y| := by
    have := abs_sub (t₀ + y) t₀
    rw [add_sub_cancel_left] at this
    linarith
  have hT2 : |t₀ + y| ≤ 4 * u ^ 8 := by
    have := abs_add_le t₀ y
    linarith
  have hz := hC (t₀ + y) hT1
  rw [Complex.ofReal_add] at hz
  have hr : |t₀ + y| ^ (η / 2 + 1 / 4) ≤ 2 * u ^ 4 := by
    calc |t₀ + y| ^ (η / 2 + 1 / 4) ≤ |t₀ + y| ^ (1 / 2 : ℝ) :=
          Real.rpow_le_rpow_of_exponent_le hT1 (by linarith)
      _ = √|t₀ + y| := (Real.sqrt_eq_rpow _).symm
      _ ≤ 2 * u ^ 4 := Real.sqrt_le_iff.2 ⟨by positivity, by nlinarith⟩
  have hZ : ‖riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I)‖ ≤ |C| * (2 * u ^ 4) := by
    calc _ ≤ C * |t₀ + y| ^ (η / 2 + 1 / 4) := hz
      _ ≤ |C| * |t₀ + y| ^ (η / 2 + 1 / 4) :=
          mul_le_mul_of_nonneg_right (le_abs_self C) (Real.rpow_nonneg (abs_nonneg _) _)
      _ ≤ |C| * (2 * u ^ 4) := mul_le_mul_of_nonneg_left hr (abs_nonneg C)
  have hP : ‖(u : ℂ) ^ (((1 - η : ℝ) : ℂ) + (t₀ + y) * I)‖ ≤ u := by
    have hre : (((1 - η : ℝ) : ℂ) + (t₀ + y) * I).re = 1 - η := by simp
    rw [Complex.norm_cpow_eq_rpow_re_of_pos hu0, hre]
    calc u ^ (1 - η) ≤ u ^ (1 : ℝ) := Real.rpow_le_rpow_of_exponent_le hu1 (by linarith)
      _ = u := Real.rpow_one u
  have hA : ((1 - σ₀ - η : ℝ) : ℂ) + y * I ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  have hB : ((σ₀ - η : ℝ) : ℂ) + y * I ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  have hAn : u ^ 8 ≤ ‖((1 - σ₀ - η : ℝ) : ℂ) + y * I‖ := by
    have := Complex.abs_im_le_norm (((1 - σ₀ - η : ℝ) : ℂ) + y * I)
    have him : (((1 - σ₀ - η : ℝ) : ℂ) + y * I).im = y := by simp
    rw [him, hy] at this
    exact this
  have hBn : u ^ 8 ≤ ‖((σ₀ - η : ℝ) : ℂ) + y * I‖ := by
    have := Complex.abs_im_le_norm (((σ₀ - η : ℝ) : ℂ) + y * I)
    have him : (((σ₀ - η : ℝ) : ℂ) + y * I).im = y := by simp
    rw [him, hy] at this
    exact this
  have hdiff : 1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I) =
      ((2 * σ₀ - 1 : ℝ) : ℂ) / ((((1 - σ₀ - η : ℝ) : ℂ) + y * I) * (((σ₀ - η : ℝ) : ℂ) + y * I)) := by
    rw [div_sub_div _ _ hA hB]
    congr 1
    push_cast
    ring
  have hD : ‖1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)‖ ≤
      1 / u ^ 16 := by
    rw [hdiff, norm_div, norm_mul, Complex.norm_real, Real.norm_eq_abs]
    have hab : |2 * σ₀ - 1| ≤ 1 := abs_le.2 ⟨by linarith, by linarith⟩
    have hAB := mul_le_mul hAn hBn (by positivity) (norm_nonneg _)
    rw [div_le_div_iff₀ (mul_pos (norm_pos_iff.2 hA) (norm_pos_iff.2 hB)) (pow_pos hu0 16)]
    calc |2 * σ₀ - 1| * u ^ 16 ≤ 1 * u ^ 16 :=
          mul_le_mul_of_nonneg_right hab (by positivity)
      _ ≤ _ := by linarith
  rw [ψ, norm_mul, norm_mul]
  calc _ ≤ (1 / u ^ 16) * u * (|C| * (2 * u ^ 4)) :=
        mul_le_mul (mul_le_mul hD hP (norm_nonneg _) (by positivity)) hZ (norm_nonneg _)
          (mul_nonneg (by positivity) hu0.le)
    _ = 2 * |C| / u ^ 11 := by field_simp

lemma bdry_bound (hin : Inputs) {η σ₀ : ℝ} (t₀ : ℝ) (hη : 0 < η) (h1 : η < σ₀)
    (h2 : η < 1 - σ₀) :
    ∃ K U : ℝ, 1 ≤ U ∧ 0 ≤ K ∧ ∀ u : ℝ, U ≤ u →
      ‖8 * (u : ℂ) ^ 7 * (ψ η σ₀ t₀ u (u ^ 8) + ψ η σ₀ t₀ u (-(u ^ 8)))‖ ≤ K / u ^ 2 := by
  obtain ⟨C, U, hU, hC, hψ⟩ := norm_ψ_le hin t₀ hη h1 h2
  refine ⟨16 * C, U, hU, by linarith, fun u hu => ?_⟩
  have hu1 : 1 ≤ u := hU.trans hu
  have hu0 : 0 < u := by linarith
  have hu0' : u ≠ 0 := hu0.ne'
  have h8 : (0 : ℝ) ≤ u ^ 8 := (pow_pos hu0 8).le
  have e1 := hψ u hu (u ^ 8) (abs_of_nonneg h8)
  have e2 := hψ u hu (-(u ^ 8)) (by rw [abs_neg, abs_of_nonneg h8])
  have hs := (norm_add_le _ _).trans (add_le_add e1 e2)
  have h8n : ‖(8 : ℂ)‖ = 8 := by norm_num
  rw [norm_mul, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hu0, h8n]
  calc _ ≤ 8 * u ^ 7 * (C / u ^ 11 + C / u ^ 11) :=
        mul_le_mul_of_nonneg_left hs (mul_nonneg (by norm_num) (pow_pos hu0 7).le)
    _ = 16 * C / u ^ 4 := by field_simp; ring
    _ ≤ 16 * C / u ^ 2 := by
      rw [div_le_div_iff₀ (pow_pos hu0 4) (pow_pos hu0 2)]
      have : u ^ 2 ≤ u ^ 4 := pow_le_pow_right₀ hu1 (by norm_num)
      exact mul_le_mul_of_nonneg_left this (by linarith)

/-- ‖u⁻¹X + c·bd‖ ≥ u⁻¹ m/2 when ‖X - T‖ ≤ K₁/u, ‖T‖ ≥ m, ‖bd‖ ≤ K₂/u², u large. -/
lemma lower_aux {Xe T bd : ℂ} {u m K1 K2 c : ℝ} (hu0 : 0 < u) (hc0 : 0 ≤ c) (hc1 : c ≤ 1)
    (hT : m ≤ ‖T‖) (hX : ‖Xe - T‖ ≤ K1 / u) (hB : ‖bd‖ ≤ K2 / u ^ 2)
    (hU : 2 * (K1 + K2) / m ≤ u) (hm : 0 < m) :
    u⁻¹ * (m / 2) ≤ ‖(u : ℂ)⁻¹ * Xe + (c : ℂ) * bd‖ := by
  have h1 : ‖(u : ℂ)⁻¹ * Xe‖ = u⁻¹ * ‖Xe‖ := by
    rw [norm_mul, norm_inv, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hu0]
  have h2 : ‖(c : ℂ) * bd‖ ≤ K2 / u ^ 2 := by
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hc0]
    calc c * ‖bd‖ ≤ 1 * ‖bd‖ := mul_le_mul_of_nonneg_right hc1 (norm_nonneg _)
      _ ≤ K2 / u ^ 2 := by rw [one_mul]; exact hB
  have h3 : ‖T‖ ≤ ‖Xe‖ + ‖Xe - T‖ := by
    have := norm_sub_le Xe (Xe - T)
    rwa [sub_sub_cancel] at this
  have h4 : ‖(u : ℂ)⁻¹ * Xe‖ ≤ ‖(u : ℂ)⁻¹ * Xe + (c : ℂ) * bd‖ + ‖(c : ℂ) * bd‖ := by
    have := norm_sub_le ((u : ℂ)⁻¹ * Xe + (c : ℂ) * bd) ((c : ℂ) * bd)
    rwa [add_sub_cancel_right] at this
  have h5 : (K1 + K2) / u ≤ m / 2 := by
    rw [div_le_iff₀ hu0]
    rw [div_le_iff₀ hm] at hU
    linarith
  have hv : 0 < u⁻¹ := inv_pos.2 hu0
  have h6 : m - K1 / u ≤ ‖Xe‖ := by linarith
  have h7 : u⁻¹ * (m - K1 / u) - K2 / u ^ 2 = u⁻¹ * (m - (K1 + K2) / u) := by ring
  have h8 : u⁻¹ * (m / 2) ≤ u⁻¹ * (m - (K1 + K2) / u) :=
    mul_le_mul_of_nonneg_left (by linarith) hv.le
  have h9 : u⁻¹ * (m - K1 / u) ≤ u⁻¹ * ‖Xe‖ := mul_le_mul_of_nonneg_left h6 hv.le
  calc u⁻¹ * (m / 2) ≤ u⁻¹ * (m - (K1 + K2) / u) := h8
    _ = u⁻¹ * (m - K1 / u) - K2 / u ^ 2 := h7.symm
    _ ≤ u⁻¹ * ‖Xe‖ - ‖(c : ℂ) * bd‖ := sub_le_sub h9 h2
    _ = ‖(u : ℂ)⁻¹ * Xe‖ - ‖(c : ℂ) * bd‖ := by rw [h1]
    _ ≤ _ := by linarith

/-- F/G = O(u^{p-1}) on D_γ. -/
lemma FG_bound (hin : Inputs) (X : Setting) {γ : ℝ} (hγ : 0 < γ) (hγ' : γ < 1 / 4) :
    ∃ K U : ℝ, 1 ≤ U ∧ 0 ≤ K ∧ ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖Ffun X.η X.σ₀ X.t₀ u / Gfun X.p u‖ ≤ K * u ^ (-(1 - X.p)) := by
  obtain ⟨K, U, hU, hK, hF⟩ := F_bound hin X hγ hγ'
  refine ⟨K, U, hU, hK, fun M₀ u hu hUu => ?_⟩
  have hu0 : 0 < u := by linarith
  rw [norm_FG _ hu0]
  calc ‖Ffun X.η X.σ₀ X.t₀ u‖ / u ^ (-X.p) ≤ K / u / u ^ (-X.p) :=
        div_le_div_of_nonneg_right (hF M₀ u hu hUu) (Real.rpow_nonneg hu0.le _)
    _ = K * (u⁻¹ / u ^ (-X.p)) := by ring
    _ = K * u ^ (-(1 - X.p)) := by rw [inv_div_rpow hu0, show (-1 - -X.p) = -(1 - X.p) by ring]

end

end Fayed

import Riemann.Deep.LeanNBRate.HighR
import Riemann.Deep.LeanNBRate.Legendre
import Mathlib.Analysis.Calculus.Deriv.Star

/-!
# Low zeros: Christoffel transfer, Markov bound, diagonal and cross terms ([R] Claim F)

For an exponential sum r with all nodes at 0 (a polynomial in u of degree < M):
|r(0)|² ≤ (M²/ℓ) ∫_0^ℓ |r|² (Christoffel), sup|r′| ≤ (C_M/ℓ) sup|r| (Bernstein with w = 0).
-/

open MeasureTheory Filter Topology Set Complex Polynomial

noncomputable section

namespace Riemann.Deep.LeanNBRate

namespace ExpSum

/-- All nodes at 0. -/
def IsPoly (r : ExpSum) : Prop := ∀ j, r.node j = 0

theorem IsPoly.dmod {r : ExpSum} (h : r.IsPoly) (c : ℝ) : (r.dmod c).IsPoly := fun j => h j

/-- The polynomial P with r(ℓ x) = P(x). -/
def polyOf (r : ExpSum) (ℓ : ℝ) : ℂ[X] :=
  ∑ j, ∑ k ∈ Finset.range (r.mult j), C (r.coef j k * (I * ℓ) ^ k) * X ^ k

theorem polyOf_eval {r : ExpSum} (h : r.IsPoly) (ℓ x : ℝ) :
    (r.polyOf ℓ).eval (x : ℂ) = r.eval (ℓ * x) := by
  unfold polyOf eval
  rw [eval_finsetSum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [eval_finsetSum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [h j]
  simp only [eval_mul, eval_C, eval_pow, eval_X, ofReal_zero, mul_zero, zero_mul, Complex.exp_zero,
    mul_one]
  push_cast; ring

theorem mult_le_M (r : ExpSum) (j : Fin r.J) : r.mult j ≤ r.M :=
  Finset.single_le_sum (f := r.mult) (fun _ _ => Nat.zero_le _) (Finset.mem_univ j)

theorem polyOf_degree (r : ExpSum) (ℓ : ℝ) : (r.polyOf ℓ).degree < r.M := by
  rw [degree_lt_iff_coeff_zero]
  intro n hn
  unfold polyOf
  rw [finset_sum_coeff]
  refine Finset.sum_eq_zero fun j _ => ?_
  rw [finset_sum_coeff]
  refine Finset.sum_eq_zero fun k hk => ?_
  rw [Finset.mem_range] at hk
  rw [coeff_C_mul, coeff_X_pow, if_neg (by have := r.mult_le_M j; omega), mul_zero]

end ExpSum

/-- **Christoffel transfer**: |r(0)|² ≤ (M²/ℓ) ∫_0^ℓ |r|² for r with all nodes at 0. -/
theorem christoffel_expSum {r : ExpSum} (h : r.IsPoly) {ℓ : ℝ} (hℓ : 0 < ℓ) :
    ‖r.eval 0‖ ^ 2 ≤ (r.M : ℝ) ^ 2 / ℓ * ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2 := by
  have hC := christoffel r.M (r.polyOf ℓ) (r.polyOf_degree ℓ)
  have h0 : (r.polyOf ℓ).eval 0 = r.eval 0 := by
    have := ExpSum.polyOf_eval h ℓ 0
    simpa using this
  rw [h0] at hC
  have hint : ∫ x in (0 : ℝ)..1, ‖(r.polyOf ℓ).eval (x : ℂ)‖ ^ 2 =
      ℓ⁻¹ * ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2 := by
    simp_rw [ExpSum.polyOf_eval h ℓ]
    rw [intervalIntegral.integral_comp_mul_left (fun u => ‖r.eval u‖ ^ 2) hℓ.ne', smul_eq_mul,
      mul_zero, mul_one]
  rw [hint] at hC
  calc ‖r.eval 0‖ ^ 2 ≤ (r.M : ℝ) ^ 2 * (ℓ⁻¹ * ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2) := hC
    _ = _ := by rw [div_eq_mul_inv]; ring

/-- Markov-type bound from Bernstein (w = 0): sup|r′| ≤ (markovC c M/ℓ) sup|r|. -/
def markovC (c : ℝ) (M : ℕ) : ℝ := (2 * c) ^ (M - 1) * CB * ((M : ℝ) + 1) ^ 3 * (4 * A₀ * ((M : ℝ) + 1) ^ 2)

theorem markov (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ (r : ExpSum) (ℓ K : ℝ), r.IsPoly → 0 < ℓ → 1 ≤ r.M → (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ K) →
      ∀ x ∈ Icc 0 ℓ, ‖(r.dmod 0).eval x‖ ≤ markovC c r.M / ℓ * K := by
  obtain ⟨c, hc, hB⟩ := bernstein hT
  refine ⟨c, hc, fun r ℓ K hpoly hℓ hM hK x hx => ?_⟩
  have hA0 := A₀_pos
  set η := 4 * A₀ * ((r.M : ℝ) + 1) ^ 2 / ℓ
  have hη : 0 < η := by positivity
  have hnodes : ∀ j, r.node j ∈ Icc (-(0 : ℝ)) 0 := fun j => by rw [hpoly j]; simp
  have hℓη : 4 * A₀ / η * ((r.M : ℝ) + 1) ^ 2 ≤ ℓ := by
    simp only [η]; field_simp; rfl
  have h := hB r 0 η ℓ K le_rfl hη hℓ hM hnodes hℓη hK x hx
  refine h.trans (le_of_eq ?_)
  simp only [markovC, η, zero_add]
  field_simp

/-- Oscillatory integration by parts: |∫_0^ℓ e^{iΩu} F| ≤ (|F(0)| + |F(ℓ)| + ∫_0^ℓ |F′|)/|Ω|. -/
theorem ibp_osc_bound {ℓ Ω : ℝ} (hℓ : 0 ≤ ℓ) (hΩ : Ω ≠ 0) {F F' : ℝ → ℂ}
    (hF : ∀ u, HasDerivAt F (F' u) u) (hF' : Continuous F') :
    ‖∫ u in (0 : ℝ)..ℓ, cexp (I * Ω * u) * F u‖ ≤
      (‖F 0‖ + ‖F ℓ‖ + ∫ u in (0 : ℝ)..ℓ, ‖F' u‖) / |Ω| := by
  have hIΩ : (I * Ω : ℂ) ≠ 0 := mul_ne_zero I_ne_zero (by exact_mod_cast hΩ)
  set v : ℝ → ℂ := fun u => cexp (I * Ω * u) / (I * Ω)
  have hv : ∀ u : ℝ, HasDerivAt v (cexp (I * Ω * u)) u := by
    intro u
    have h1 : HasDerivAt (fun u : ℝ => (I * Ω * (u : ℂ))) (I * Ω) u := by
      simpa using ((hasDerivAt_id (u : ℂ)).const_mul (I * (Ω : ℂ))).comp_ofReal
    have := (h1.cexp).div_const (I * Ω)
    convert this using 1
    rw [mul_div_assoc, div_self hIΩ, mul_one]
  have hFc : Continuous F := continuous_iff_continuousAt.mpr fun u => (hF u).continuousAt
  have hibp := intervalIntegral.integral_mul_deriv_eq_deriv_mul (a := (0 : ℝ)) (b := ℓ)
    (fun u _ => hF u) (fun u _ => hv u) (hF'.intervalIntegrable _ _)
    ((by fun_prop : Continuous fun u : ℝ => cexp (I * Ω * u)).intervalIntegrable _ _)
  have e : ∫ u in (0 : ℝ)..ℓ, cexp (I * Ω * u) * F u = ∫ u in (0 : ℝ)..ℓ, F u * cexp (I * Ω * u) := by
    congr 1; ext u; ring
  rw [e, hibp]
  have hvn : ∀ u : ℝ, ‖v u‖ = 1 / |Ω| := by
    intro u
    simp only [v, norm_div, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs]
    rw [show (I * (Ω : ℂ) * u) = ((Ω * u : ℝ) : ℂ) * I by push_cast; ring, Complex.norm_exp_ofReal_mul_I]
  have hint : ‖∫ u in (0 : ℝ)..ℓ, F' u * v u‖ ≤ (∫ u in (0 : ℝ)..ℓ, ‖F' u‖) / |Ω| := by
    refine (intervalIntegral.norm_integral_le_integral_norm hℓ).trans (le_of_eq ?_)
    simp_rw [norm_mul, hvn]
    rw [intervalIntegral.integral_mul_const]; ring
  have hΩp : 0 < |Ω| := abs_pos.mpr hΩ
  calc ‖F ℓ * v ℓ - F 0 * v 0 - ∫ u in (0 : ℝ)..ℓ, F' u * v u‖
      ≤ ‖F ℓ * v ℓ‖ + ‖F 0 * v 0‖ + ‖∫ u in (0 : ℝ)..ℓ, F' u * v u‖ := by
        refine (norm_sub_le _ _).trans ?_
        gcongr
        exact norm_sub_le _ _
    _ ≤ ‖F ℓ‖ / |Ω| + ‖F 0‖ / |Ω| + (∫ u in (0 : ℝ)..ℓ, ‖F' u‖) / |Ω| := by
        rw [norm_mul, norm_mul, hvn, hvn]
        gcongr <;> exact le_of_eq (by ring)
    _ = _ := by ring

open ComplexConjugate in
/-- The OU bilinear form on [0, ℓ] without endpoint terms: ∫_0^ℓ (f′ ḡ′ + ¼ f ḡ). -/
def oubil (ℓ : ℝ) (f f' g g' : ℝ → ℂ) : ℂ :=
  ∫ u in (0 : ℝ)..ℓ, (f' u * conj (g' u) + (1 / 4 : ℂ) * (f u * conj (g u)))

open ComplexConjugate in
theorem conj_cexp_I_mul (γ u : ℝ) : conj (cexp (I * γ * u)) = cexp (-(I * γ * u)) := by
  rw [← Complex.exp_conj]; congr 1; simp [Complex.conj_ofReal]

open ComplexConjugate in
/-- Generic cross-term bound for f = e^{iγ₁u} r₁, g = e^{iγ₂u} r₂ (γ₁ ≠ γ₂). -/
theorem cross_bound (r₁ r₂ : ExpSum) {γ₁ γ₂ ℓ : ℝ} (hℓ : 0 ≤ ℓ) (hγ : γ₁ ≠ γ₂)
    {a₁ a₂ a₁' a₂' s₁ s₂ s₁' s₂' : ℝ}
    (hA₁ : ∀ u ∈ Icc 0 ℓ, ‖(r₁.dmod γ₁).eval u‖ ≤ a₁) (hA₂ : ∀ u ∈ Icc 0 ℓ, ‖(r₂.dmod γ₂).eval u‖ ≤ a₂)
    (hA₁' : ∀ u ∈ Icc 0 ℓ, ‖((r₁.dmod γ₁).dmod 0).eval u‖ ≤ a₁')
    (hA₂' : ∀ u ∈ Icc 0 ℓ, ‖((r₂.dmod γ₂).dmod 0).eval u‖ ≤ a₂')
    (hs₁ : ∀ u ∈ Icc 0 ℓ, ‖r₁.eval u‖ ≤ s₁) (hs₂ : ∀ u ∈ Icc 0 ℓ, ‖r₂.eval u‖ ≤ s₂)
    (hs₁' : ∀ u ∈ Icc 0 ℓ, ‖(r₁.dmod 0).eval u‖ ≤ s₁') (hs₂' : ∀ u ∈ Icc 0 ℓ, ‖(r₂.dmod 0).eval u‖ ≤ s₂') :
    ‖oubil ℓ (fun u => cexp (I * γ₁ * u) * r₁.eval u) (fun u => cexp (I * γ₁ * u) * (r₁.dmod γ₁).eval u)
        (fun u => cexp (I * γ₂ * u) * r₂.eval u) (fun u => cexp (I * γ₂ * u) * (r₂.dmod γ₂).eval u)‖ ≤
      (2 * (a₁ * a₂ + s₁ * s₂ / 4) + ℓ * (a₁' * a₂ + a₁ * a₂' + (s₁' * s₂ + s₁ * s₂') / 4)) /
        |γ₁ - γ₂| := by
  set A₁ := (r₁.dmod γ₁).eval
  set A₂ := (r₂.dmod γ₂).eval
  set F : ℝ → ℂ := fun u => A₁ u * conj (A₂ u) + (1 / 4 : ℂ) * (r₁.eval u * conj (r₂.eval u))
  set F' : ℝ → ℂ := fun u => ((r₁.dmod γ₁).dmod 0).eval u * conj (A₂ u) +
    A₁ u * conj (((r₂.dmod γ₂).dmod 0).eval u) +
    (1 / 4 : ℂ) * ((r₁.dmod 0).eval u * conj (r₂.eval u) + r₁.eval u * conj ((r₂.dmod 0).eval u))
  have hF : ∀ u, HasDerivAt F (F' u) u := by
    intro u
    have h1 := (r₁.dmod γ₁).hasDerivAt_eval u
    have h2 := ((r₂.dmod γ₂).hasDerivAt_eval u).star
    have h3 := r₁.hasDerivAt_eval u
    have h4 := (r₂.hasDerivAt_eval u).star
    have := (h1.mul h2).add ((h3.mul h4).const_mul (1 / 4 : ℂ))
    convert this using 1
    · ext x; simp only [F, A₁, A₂, starRingEnd_apply, Pi.add_apply, Pi.mul_apply]
    · simp only [F', A₁, A₂, starRingEnd_apply]
  have hF'c : Continuous F' := by
    simp only [F', A₁, A₂]
    have c1 := ((r₁.dmod γ₁).dmod 0).continuous_eval
    have c2 := (r₂.dmod γ₂).continuous_eval
    have c3 := (r₁.dmod γ₁).continuous_eval
    have c4 := ((r₂.dmod γ₂).dmod 0).continuous_eval
    have c5 := (r₁.dmod 0).continuous_eval
    have c6 := r₂.continuous_eval
    have c7 := r₁.continuous_eval
    have c8 := (r₂.dmod 0).continuous_eval
    fun_prop
  have hint : oubil ℓ (fun u => cexp (I * γ₁ * u) * r₁.eval u)
      (fun u => cexp (I * γ₁ * u) * (r₁.dmod γ₁).eval u)
      (fun u => cexp (I * γ₂ * u) * r₂.eval u) (fun u => cexp (I * γ₂ * u) * (r₂.dmod γ₂).eval u) =
      ∫ u in (0 : ℝ)..ℓ, cexp (I * ((γ₁ - γ₂ : ℝ) : ℂ) * u) * F u := by
    unfold oubil
    congr 1; ext u
    simp only [F, A₁, A₂, map_mul, conj_cexp_I_mul]
    have e : cexp (I * ((γ₁ - γ₂ : ℝ) : ℂ) * u) = cexp (I * γ₁ * u) * cexp (-(I * γ₂ * u)) := by
      rw [← Complex.exp_add]; push_cast; ring_nf
    rw [e]; ring
  rw [hint]
  refine (ibp_osc_bound hℓ (sub_ne_zero.mpr hγ) hF hF'c).trans ?_
  have h0 : (0 : ℝ) ∈ Icc 0 ℓ := ⟨le_rfl, hℓ⟩
  have hl : ℓ ∈ Icc 0 ℓ := ⟨hℓ, le_rfl⟩
  have hFb : ∀ u ∈ Icc 0 ℓ, ‖F u‖ ≤ a₁ * a₂ + s₁ * s₂ / 4 := by
    intro u hu
    simp only [F]
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, norm_mul, norm_mul, Complex.norm_conj, Complex.norm_conj]
    have : ‖(1 / 4 : ℂ)‖ = 1 / 4 := by norm_num
    rw [this]
    have e1 := mul_le_mul (hA₁ u hu) (hA₂ u hu) (norm_nonneg _) ((norm_nonneg _).trans (hA₁ u hu))
    have e2 := mul_le_mul (hs₁ u hu) (hs₂ u hu) (norm_nonneg _) ((norm_nonneg _).trans (hs₁ u hu))
    linarith
  have hF'b : ∀ u ∈ Icc 0 ℓ, ‖F' u‖ ≤ a₁' * a₂ + a₁ * a₂' + (s₁' * s₂ + s₁ * s₂') / 4 := by
    intro u hu
    have e1 := mul_le_mul (hA₁' u hu) (hA₂ u hu) (norm_nonneg _) ((norm_nonneg _).trans (hA₁' u hu))
    have e2 := mul_le_mul (hA₁ u hu) (hA₂' u hu) (norm_nonneg _) ((norm_nonneg _).trans (hA₁ u hu))
    have e3 := mul_le_mul (hs₁' u hu) (hs₂ u hu) (norm_nonneg _) ((norm_nonneg _).trans (hs₁' u hu))
    have e4 := mul_le_mul (hs₁ u hu) (hs₂' u hu) (norm_nonneg _) ((norm_nonneg _).trans (hs₁ u hu))
    have n1 : ‖((r₁.dmod γ₁).dmod 0).eval u * conj (A₂ u) + A₁ u * conj (((r₂.dmod γ₂).dmod 0).eval u)‖ ≤
        ‖((r₁.dmod γ₁).dmod 0).eval u‖ * ‖A₂ u‖ + ‖A₁ u‖ * ‖((r₂.dmod γ₂).dmod 0).eval u‖ := by
      refine (norm_add_le _ _).trans (le_of_eq ?_)
      rw [norm_mul, norm_mul, Complex.norm_conj, Complex.norm_conj]
    have n2 : ‖(1 / 4 : ℂ) * ((r₁.dmod 0).eval u * conj (r₂.eval u) + r₁.eval u * conj ((r₂.dmod 0).eval u))‖ ≤
        (1 / 4) * (‖(r₁.dmod 0).eval u‖ * ‖r₂.eval u‖ + ‖r₁.eval u‖ * ‖(r₂.dmod 0).eval u‖) := by
      rw [norm_mul, show ‖(1 / 4 : ℂ)‖ = 1 / 4 by norm_num]
      gcongr
      refine (norm_add_le _ _).trans (le_of_eq ?_)
      rw [norm_mul, norm_mul, Complex.norm_conj, Complex.norm_conj]
    have n3 : ‖F' u‖ ≤ ‖((r₁.dmod γ₁).dmod 0).eval u * conj (A₂ u) + A₁ u * conj (((r₂.dmod γ₂).dmod 0).eval u)‖ +
        ‖(1 / 4 : ℂ) * ((r₁.dmod 0).eval u * conj (r₂.eval u) + r₁.eval u * conj ((r₂.dmod 0).eval u))‖ :=
      norm_add_le _ _
    linarith
  have hI : ∫ u in (0 : ℝ)..ℓ, ‖F' u‖ ≤ ℓ * (a₁' * a₂ + a₁ * a₂' + (s₁' * s₂ + s₁ * s₂') / 4) := by
    have := intervalIntegral.integral_mono_on (μ := volume) hℓ (hF'c.norm.intervalIntegrable 0 ℓ)
      (continuous_const.intervalIntegrable 0 ℓ) (fun u hu => hF'b u hu)
    rw [intervalIntegral.integral_const, smul_eq_mul, sub_zero] at this
    exact this
  have hΩ : 0 < |γ₁ - γ₂| := abs_pos.mpr (sub_ne_zero.mpr hγ)
  gcongr
  have := hFb 0 h0
  have := hFb ℓ hl
  linarith

end Riemann.Deep.LeanNBRate




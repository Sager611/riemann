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

/-- The sup bound √(2(2c)^{2(M-1)} ℓ⁻¹ ∫_0^ℓ |r|²) (Nikolskii on [0, ℓ]). -/
def supS (c ℓ : ℝ) (r : ExpSum) : ℝ :=
  Real.sqrt (2 * (2 * c) ^ (2 * (r.M - 1)) / ℓ * ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2)

theorem supS_nonneg (c ℓ : ℝ) (r : ExpSum) : 0 ≤ supS c ℓ r := Real.sqrt_nonneg _

theorem supS_sq {c ℓ : ℝ} (hc : 0 < c) (hℓ : 0 < ℓ) (r : ExpSum) :
    supS c ℓ r ^ 2 = 2 * (2 * c) ^ (2 * (r.M - 1)) / ℓ * ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2 := by
  unfold supS
  rw [Real.sq_sqrt]
  have : 0 ≤ ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2 :=
    intervalIntegral.integral_nonneg hℓ.le fun u _ => by positivity
  positivity

theorem markovC_mono {c c' : ℝ} (hc : 0 < c) (hcc : c ≤ c') (M : ℕ) : markovC c M ≤ markovC c' M := by
  unfold markovC
  have := CB_nonneg
  have := A₀_pos
  have h1 : (2 * c) ^ (M - 1) ≤ (2 * c') ^ (M - 1) := pow_le_pow_left₀ (by positivity) (by linarith) _
  gcongr

theorem markovC_nonneg {c : ℝ} (hc : 0 < c) (M : ℕ) : 0 ≤ markovC c M := by
  unfold markovC; have := CB_nonneg; have := A₀_pos; positivity

/-- Sup bounds for a polynomial exponential sum and its modulated derivatives. -/
theorem poly_bounds (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ (r : ExpSum) (ℓ γ : ℝ), r.IsPoly → 0 < ℓ → 1 ≤ r.M →
      (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ supS c ℓ r) ∧
      (∀ u ∈ Icc 0 ℓ, ‖(r.dmod 0).eval u‖ ≤ markovC c r.M / ℓ * supS c ℓ r) ∧
      (∀ u ∈ Icc 0 ℓ, ‖(r.dmod γ).eval u‖ ≤ (|γ| + markovC c r.M / ℓ) * supS c ℓ r) ∧
      (∀ u ∈ Icc 0 ℓ, ‖((r.dmod γ).dmod 0).eval u‖ ≤
        markovC c r.M / ℓ * ((|γ| + markovC c r.M / ℓ) * supS c ℓ r)) := by
  obtain ⟨c₁, hc₁, hN⟩ := nikolskii hT
  obtain ⟨c₂, hc₂, hMk⟩ := markov hT
  set c := max c₁ c₂
  have hc : 0 < c := lt_max_of_lt_left hc₁
  refine ⟨c, hc, fun r ℓ γ hpoly hℓ hM => ?_⟩
  have hS : ∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ supS c ℓ r := by
    intro u hu
    apply Real.le_sqrt_of_sq_le
    have h := hN r 0 ℓ hM hℓ u hu
    rw [sub_zero, integral_Icc_eq_integral_Ioc, ← intervalIntegral.integral_of_le hℓ.le] at h
    have hI : 0 ≤ ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2 :=
      intervalIntegral.integral_nonneg hℓ.le fun u _ => by positivity
    have hp : (2 * c₁) ^ (2 * (r.M - 1)) ≤ (2 * c) ^ (2 * (r.M - 1)) :=
      pow_le_pow_left₀ (by positivity) (by linarith [le_max_left c₁ c₂]) _
    calc ‖r.eval u‖ ^ 2 ≤ 2 * (2 * c₁) ^ (2 * (r.M - 1)) * (ℓ⁻¹ * ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2) := h
      _ ≤ 2 * (2 * c) ^ (2 * (r.M - 1)) * (ℓ⁻¹ * ∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2) := by gcongr
      _ = _ := by rw [div_eq_mul_inv]; ring
  have hmk : markovC c₂ r.M ≤ markovC c r.M := markovC_mono hc₂ (le_max_right c₁ c₂) _
  have hS0 := supS_nonneg c ℓ r
  have hS' : ∀ u ∈ Icc 0 ℓ, ‖(r.dmod 0).eval u‖ ≤ markovC c r.M / ℓ * supS c ℓ r := by
    intro u hu
    refine (hMk r ℓ (supS c ℓ r) hpoly hℓ hM hS u hu).trans ?_
    gcongr
  have hA : ∀ u ∈ Icc 0 ℓ, ‖(r.dmod γ).eval u‖ ≤ (|γ| + markovC c r.M / ℓ) * supS c ℓ r := by
    intro u hu
    rw [r.dmod_eval_add γ u]
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs]
    have := hS u hu
    have := hS' u hu
    nlinarith [abs_nonneg γ]
  refine ⟨hS, hS', hA, fun u hu => ?_⟩
  have hpolyA : (r.dmod γ).IsPoly := hpoly.dmod γ
  have hMA : 1 ≤ (r.dmod γ).M := by rw [r.dmod_M]; exact hM
  refine (hMk (r.dmod γ) ℓ ((|γ| + markovC c r.M / ℓ) * supS c ℓ r) hpolyA hℓ hMA hA u hu).trans ?_
  rw [r.dmod_M]
  have : 0 ≤ (|γ| + markovC c r.M / ℓ) * supS c ℓ r := by
    have := markovC_nonneg hc r.M
    positivity
  gcongr

/-- Diagonal lower bound: ∫_0^ℓ (|(r.dmod γ)|² + ¼|r|²) ≥ (γ² + ¼) ∫|r|² - 2|γ| ℓ S S'. -/
theorem diag_lower (r : ExpSum) {γ ℓ S S' : ℝ} (hℓ : 0 ≤ ℓ)
    (hS : ∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ S) (hS' : ∀ u ∈ Icc 0 ℓ, ‖(r.dmod 0).eval u‖ ≤ S') :
    (γ ^ 2 + 1 / 4) * (∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2) - 2 * |γ| * ℓ * (S * S') ≤
      ∫ u in (0 : ℝ)..ℓ, (‖(r.dmod γ).eval u‖ ^ 2 + 1 / 4 * ‖r.eval u‖ ^ 2) := by
  have hc1 := r.continuous_eval
  have hc2 := (r.dmod γ).continuous_eval
  have hpt : ∀ u ∈ Icc 0 ℓ, (γ ^ 2 + 1 / 4) * ‖r.eval u‖ ^ 2 - 2 * |γ| * (S * S') ≤
      ‖(r.dmod γ).eval u‖ ^ 2 + 1 / 4 * ‖r.eval u‖ ^ 2 := by
    intro u hu
    rw [r.dmod_eval_add γ u]
    set a := I * (γ : ℂ) * r.eval u
    set b := (r.dmod 0).eval u
    have ha : ‖a‖ = |γ| * ‖r.eval u‖ := by
      simp only [a, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs]
    have hab : ‖a‖ ^ 2 - 2 * ‖a‖ * ‖b‖ ≤ ‖a + b‖ ^ 2 := by
      have h1 : ‖a‖ - ‖b‖ ≤ ‖a + b‖ := by
        have := norm_sub_norm_le a (-b); rw [norm_neg, sub_neg_eq_add] at this; exact this
      rcases le_total ‖b‖ ‖a‖ with h | h
      · have := pow_le_pow_left₀ (by linarith) h1 2
        nlinarith [norm_nonneg b]
      · have : 0 ≤ ‖a + b‖ ^ 2 := by positivity
        nlinarith [norm_nonneg a, norm_nonneg b]
    have hr := hS u hu
    have hb := hS' u hu
    have hγ0 := abs_nonneg γ
    have hrb : ‖r.eval u‖ * ‖b‖ ≤ S * S' :=
      mul_le_mul hr hb (norm_nonneg _) ((norm_nonneg _).trans hr)
    rw [ha] at hab
    have hsq : |γ| ^ 2 = γ ^ 2 := sq_abs γ
    nlinarith
  have hi1 : IntervalIntegrable (fun u => (γ ^ 2 + 1 / 4) * ‖r.eval u‖ ^ 2 - 2 * |γ| * (S * S'))
      volume 0 ℓ :=
    ((continuous_const.mul (hc1.norm.pow 2)).sub continuous_const).intervalIntegrable 0 ℓ
  have hi2 : IntervalIntegrable (fun u => ‖(r.dmod γ).eval u‖ ^ 2 + 1 / 4 * ‖r.eval u‖ ^ 2) volume 0 ℓ :=
    ((hc2.norm.pow 2).add (continuous_const.mul (hc1.norm.pow 2))).intervalIntegrable 0 ℓ
  have hint : ∫ u in (0 : ℝ)..ℓ, ((γ ^ 2 + 1 / 4) * ‖r.eval u‖ ^ 2 - 2 * |γ| * (S * S')) ≤
      ∫ u in (0 : ℝ)..ℓ, (‖(r.dmod γ).eval u‖ ^ 2 + 1 / 4 * ‖r.eval u‖ ^ 2) :=
    intervalIntegral.integral_mono_on hℓ hi1 hi2 hpt
  have e : ∫ u in (0 : ℝ)..ℓ, ((γ ^ 2 + 1 / 4) * ‖r.eval u‖ ^ 2 - 2 * |γ| * (S * S')) =
      (γ ^ 2 + 1 / 4) * (∫ u in (0 : ℝ)..ℓ, ‖r.eval u‖ ^ 2) - 2 * |γ| * ℓ * (S * S') := by
    have j1 : IntervalIntegrable (fun u => (γ ^ 2 + 1 / 4) * ‖r.eval u‖ ^ 2) volume 0 ℓ :=
      (continuous_const.mul (hc1.norm.pow 2)).intervalIntegrable 0 ℓ
    have j2 : IntervalIntegrable (fun _ : ℝ => 2 * |γ| * (S * S')) volume 0 ℓ :=
      continuous_const.intervalIntegrable 0 ℓ
    rw [intervalIntegral.integral_sub j1 j2, intervalIntegral.integral_const_mul,
      intervalIntegral.integral_const, smul_eq_mul, sub_zero]
    ring
  linarith

theorem part_singleton_isPoly (P : ExpSum) (γ : ℝ) : (P.part {γ} γ).IsPoly := by
  classical
  intro j
  show (if P.node j ∈ ({γ} : Set ℝ) then P.node j - γ else 0) = 0
  split_ifs with h
  · rw [Set.mem_singleton_iff.mp h, sub_self]
  · rfl

theorem part_eval_singletons (P : ExpSum) (T : Set ℝ) (u : ℝ) :
    ∑ γ ∈ nodeSet P T, cexp (I * γ * u) * (P.part {γ} γ).eval u = (P.part T 0).eval u := by
  classical
  have h := P.part_eval_sum T ((nodeSet P T).image fun γ => ({γ} : Finset ℝ))
    (fun j hj => ⟨{P.node j}, Finset.mem_image.mpr ⟨P.node j, mem_nodeSet.mpr ⟨⟨j, rfl⟩, hj⟩, rfl⟩,
      Finset.mem_singleton_self _⟩)
    (fun g hg x hx => by
      obtain ⟨γ, hγ, rfl⟩ := Finset.mem_image.mp hg
      rw [Finset.mem_singleton.mp hx]; exact (mem_nodeSet.mp hγ).2)
    (fun g hg g' hg' x hx hx' => by
      obtain ⟨γ, _, rfl⟩ := Finset.mem_image.mp hg
      obtain ⟨γ', _, rfl⟩ := Finset.mem_image.mp hg'
      rw [Finset.mem_singleton.mp hx] at hx'
      rw [Finset.mem_singleton.mp hx']) u
  rw [Finset.sum_image (fun a _ b _ hab => Finset.singleton_injective hab)] at h
  rw [← h]
  refine Finset.sum_congr rfl fun γ _ => ?_
  rw [P.part_eval_shift, Finset.coe_singleton]

theorem part_dmod_eval_singletons (P : ExpSum) (T : Set ℝ) (u : ℝ) :
    ∑ γ ∈ nodeSet P T, cexp (I * γ * u) * ((P.part {γ} γ).dmod γ).eval u =
      ((P.part T 0).dmod 0).eval u := by
  have hF : (fun v : ℝ => ∑ γ ∈ nodeSet P T, cexp (I * γ * v) * (P.part {γ} γ).eval v) =
      (P.part T 0).eval := by
    funext v; exact part_eval_singletons P T v
  have h1 : HasDerivAt (fun v : ℝ => ∑ γ ∈ nodeSet P T, cexp (I * γ * v) * (P.part {γ} γ).eval v)
      (∑ γ ∈ nodeSet P T, cexp (I * γ * u) * ((P.part {γ} γ).dmod γ).eval u) u :=
    HasDerivAt.fun_sum fun γ _ => (P.part {γ} γ).hasDerivAt_mod_eval γ u
  rw [hF] at h1
  exact h1.unique ((P.part T 0).hasDerivAt_eval u)

open ComplexConjugate in
theorem oubil_sum {ι : Type*} (s : Finset ι) (f f' g g' : ι → ℝ → ℂ)
    (hf : ∀ i, Continuous (f i)) (hf' : ∀ i, Continuous (f' i))
    (hg : ∀ i, Continuous (g i)) (hg' : ∀ i, Continuous (g' i)) (ℓ : ℝ) :
    oubil ℓ (fun u => ∑ i ∈ s, f i u) (fun u => ∑ i ∈ s, f' i u) (fun u => ∑ j ∈ s, g j u)
        (fun u => ∑ j ∈ s, g' j u) =
      ∑ i ∈ s, ∑ j ∈ s, oubil ℓ (f i) (f' i) (g j) (g' j) := by
  unfold oubil
  have e : ∀ u : ℝ, (∑ i ∈ s, f' i u) * conj (∑ j ∈ s, g' j u) +
      (1 / 4 : ℂ) * ((∑ i ∈ s, f i u) * conj (∑ j ∈ s, g j u)) =
      ∑ i ∈ s, ∑ j ∈ s, (f' i u * conj (g' j u) + (1 / 4 : ℂ) * (f i u * conj (g j u))) := by
    intro u
    rw [map_sum, map_sum, Finset.sum_mul_sum, Finset.sum_mul_sum, Finset.mul_sum,
      ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  simp_rw [e]
  rw [intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun i _ => ?_
    rw [intervalIntegral.integral_finsetSum]
    intro j _
    exact (((hf' i).mul ((hg' j).star)).add (continuous_const.mul ((hf i).mul (hg j).star))).intervalIntegrable _ _
  · intro i _
    exact (continuous_finsetSum _ fun j _ => (((hf' i).mul ((hg' j).star)).add
      (continuous_const.mul ((hf i).mul (hg j).star)))).intervalIntegrable _ _

open ComplexConjugate in
theorem oubil_diag (r : ExpSum) (γ ℓ : ℝ) :
    oubil ℓ (fun u => cexp (I * γ * u) * r.eval u) (fun u => cexp (I * γ * u) * (r.dmod γ).eval u)
      (fun u => cexp (I * γ * u) * r.eval u) (fun u => cexp (I * γ * u) * (r.dmod γ).eval u) =
      ((∫ u in (0 : ℝ)..ℓ, (‖(r.dmod γ).eval u‖ ^ 2 + 1 / 4 * ‖r.eval u‖ ^ 2) : ℝ) : ℂ) := by
  unfold oubil
  rw [← intervalIntegral.integral_ofReal]
  congr 1; ext u
  have h1 : ∀ z : ℂ, cexp (I * γ * u) * z * conj (cexp (I * γ * u) * z) = ((‖z‖ ^ 2 : ℝ) : ℂ) := by
    intro z
    rw [map_mul, conj_cexp_I_mul]
    have : cexp (I * γ * u) * z * (cexp (-(I * γ * u)) * conj z) =
        (cexp (I * γ * u) * cexp (-(I * γ * u))) * (z * conj z) := by ring
    rw [this, ← Complex.exp_add, add_neg_cancel, Complex.exp_zero, one_mul, Complex.mul_conj,
      Complex.normSq_eq_norm_sq]
  rw [h1, h1]; push_cast; ring

/-- The constant of the Claim F core. -/
def claimFC (c : ℝ) (M₀ : ℕ) (H₀ g₀ : ℝ) (J : ℕ) : ℝ :=
  8 * H₀ * markovC c M₀ * (2 * (2 * c) ^ (2 * (M₀ - 1))) +
    4 * (2 * ((H₀ + markovC c M₀) ^ 2 + 1 / 4) + 2 * markovC c M₀ * (H₀ + markovC c M₀) ^ 2 +
      markovC c M₀ / 2) * (2 * (2 * c) ^ (2 * (M₀ - 1))) * J / g₀

theorem markovC_mono2 {c₀ c : ℝ} (hc₀ : 0 < c₀) (hc : c₀ ≤ c) (hc1 : 1 ≤ c) {M M₀ : ℕ} (hM : M ≤ M₀) :
    markovC c₀ M ≤ markovC c M₀ := by
  unfold markovC
  have := CB_nonneg
  have := A₀_pos
  have h1 := pow_mono_two hc₀ hc hc1 (Nat.sub_le_sub_right hM 1)
  have h2 : ((M : ℝ) + 1) ≤ ((M₀ : ℝ) + 1) := by
    have : (M : ℝ) ≤ M₀ := by exact_mod_cast hM
    linarith
  gcongr

open ComplexConjugate in
/-- **[R] Claim F core** (finite low configuration, ℓ ≥ 1): the OU form dominates the diagonal mass
up to 1 - C/ℓ, and |p(0)|² ≤ (Σ M_γ²/(γ²+¼))·(diagonal mass)/ℓ. -/
theorem claimF_core (hT : Nazarov1993_Turan) : ∃ c : ℝ, 1 ≤ c ∧
    ∀ (P : ExpSum) (T : Set ℝ) (ℓ g₀ H₀ : ℝ) (M₀ : ℕ), 1 ≤ ℓ → 0 < g₀ → 0 ≤ H₀ →
      (∀ γ ∈ nodeSet P T, ∀ γ' ∈ nodeSet P T, γ ≠ γ' → g₀ ≤ |γ - γ'|) →
      (∀ γ ∈ nodeSet P T, |γ| ≤ H₀) →
      (∀ γ ∈ nodeSet P T, 1 ≤ (P.part {γ} γ).M ∧ (P.part {γ} γ).M ≤ M₀) →
      (1 - claimFC c M₀ H₀ g₀ (nodeSet P T).card / ℓ) *
          (∑ γ ∈ nodeSet P T, (γ ^ 2 + 1 / 4) * ∫ u in (0 : ℝ)..ℓ, ‖(P.part {γ} γ).eval u‖ ^ 2) ≤
        (oubil ℓ (P.part T 0).eval ((P.part T 0).dmod 0).eval (P.part T 0).eval
          ((P.part T 0).dmod 0).eval).re ∧
      ‖(P.part T 0).eval 0‖ ^ 2 ≤
        (∑ γ ∈ nodeSet P T, ((P.part {γ} γ).M : ℝ) ^ 2 / (γ ^ 2 + 1 / 4)) *
          ((∑ γ ∈ nodeSet P T, (γ ^ 2 + 1 / 4) * ∫ u in (0 : ℝ)..ℓ, ‖(P.part {γ} γ).eval u‖ ^ 2) / ℓ) := by
  obtain ⟨c₀, hc₀, hPB⟩ := poly_bounds hT
  set c := max c₀ 1
  have hc1 : 1 ≤ c := le_max_right _ _
  have hcc : c₀ ≤ c := le_max_left _ _
  refine ⟨c, hc1, fun P T ℓ g₀ H₀ M₀ hℓ hg₀ hH₀ hsep hH hM => ?_⟩
  have hℓ0 : 0 < ℓ := by linarith
  set S := nodeSet P T
  set r : ℝ → ExpSum := fun γ => P.part {γ} γ
  set Ig : ℝ → ℝ := fun γ => ∫ u in (0 : ℝ)..ℓ, ‖(r γ).eval u‖ ^ 2
  set Sg : ℝ → ℝ := fun γ => supS c₀ ℓ (r γ)
  set σ := 2 * (2 * c) ^ (2 * (M₀ - 1))
  set mk := markovC c M₀
  have hmk0 : 0 ≤ mk := markovC_nonneg (by linarith) M₀
  have hσ0 : 0 ≤ σ := by positivity
  have hIg0 : ∀ γ, 0 ≤ Ig γ := fun γ => intervalIntegral.integral_nonneg hℓ0.le fun u _ => by positivity
  have hSg0 : ∀ γ, 0 ≤ Sg γ := fun γ => supS_nonneg _ _ _
  have hSg2 : ∀ γ ∈ S, Sg γ ^ 2 ≤ σ / ℓ * Ig γ := by
    intro γ hγ
    simp only [Sg]
    rw [supS_sq hc₀ hℓ0]
    have hp := pow_mono_two hc₀ hcc hc1 (Nat.mul_le_mul_left 2 (Nat.sub_le_sub_right (hM γ hγ).2 1))
    have := hIg0 γ
    simp only [σ, Ig]
    gcongr
  have hmkg : ∀ γ ∈ S, markovC c₀ (r γ).M ≤ mk := fun γ hγ =>
    markovC_mono2 hc₀ hcc hc1 (hM γ hγ).2
  -- sup bounds per node
  have hB : ∀ γ ∈ S, (∀ u ∈ Icc 0 ℓ, ‖(r γ).eval u‖ ≤ Sg γ) ∧
      (∀ u ∈ Icc 0 ℓ, ‖((r γ).dmod 0).eval u‖ ≤ mk / ℓ * Sg γ) ∧
      (∀ u ∈ Icc 0 ℓ, ‖((r γ).dmod γ).eval u‖ ≤ (H₀ + mk) * Sg γ) ∧
      (∀ u ∈ Icc 0 ℓ, ‖(((r γ).dmod γ).dmod 0).eval u‖ ≤ mk / ℓ * ((H₀ + mk) * Sg γ)) ∧
      (∀ u ∈ Icc 0 ℓ, ‖((r γ).dmod 0).eval u‖ ≤ markovC c₀ (r γ).M / ℓ * Sg γ) := by
    intro γ hγ
    obtain ⟨b1, b2, b3, b4⟩ := hPB (r γ) ℓ γ (part_singleton_isPoly P γ) hℓ0 (hM γ hγ).1
    have hmg := hmkg γ hγ
    have hq : markovC c₀ (r γ).M / ℓ ≤ mk / ℓ := by gcongr
    have hq2 : mk / ℓ ≤ mk := div_le_self hmk0 hℓ
    have hγH := hH γ hγ
    have hS0 := hSg0 γ
    have hq0 : 0 ≤ markovC c₀ (r γ).M / ℓ := div_nonneg (markovC_nonneg hc₀ _) hℓ0.le
    refine ⟨b1, fun u hu => (b2 u hu).trans (by gcongr), fun u hu => (b3 u hu).trans ?_,
      fun u hu => (b4 u hu).trans ?_, b2⟩
    · gcongr; linarith
    · have : 0 ≤ (|γ| + markovC c₀ (r γ).M / ℓ) * Sg γ := by positivity
      calc markovC c₀ (r γ).M / ℓ * ((|γ| + markovC c₀ (r γ).M / ℓ) * Sg γ)
          ≤ mk / ℓ * ((|γ| + markovC c₀ (r γ).M / ℓ) * Sg γ) := by gcongr
        _ ≤ mk / ℓ * ((H₀ + mk) * Sg γ) := by gcongr; linarith
  -- the functions
  set f : ℝ → ℝ → ℂ := fun γ u => cexp (I * γ * u) * (r γ).eval u
  set f' : ℝ → ℝ → ℂ := fun γ u => cexp (I * γ * u) * ((r γ).dmod γ).eval u
  have hfc : ∀ γ, Continuous (f γ) := fun γ => by
    have := (r γ).continuous_eval
    simp only [f]; fun_prop
  have hf'c : ∀ γ, Continuous (f' γ) := fun γ => by
    have := ((r γ).dmod γ).continuous_eval
    simp only [f']; fun_prop
  have hp : (P.part T 0).eval = fun u => ∑ γ ∈ S, f γ u := by
    funext u; exact (part_eval_singletons P T u).symm
  have hp' : ((P.part T 0).dmod 0).eval = fun u => ∑ γ ∈ S, f' γ u := by
    funext u; exact (part_dmod_eval_singletons P T u).symm
  set D := ∑ γ ∈ S, (γ ^ 2 + 1 / 4) * Ig γ
  have hsumI : ∑ γ ∈ S, Ig γ ≤ 4 * D := by
    simp only [D, Finset.mul_sum]
    refine Finset.sum_le_sum fun γ _ => ?_
    have := hIg0 γ
    nlinarith [sq_nonneg γ]
  refine ⟨?_, ?_⟩
  · -- the lower bound for the form
    rw [hp, hp', oubil_sum S f f' f f' hfc hf'c hfc hf'c ℓ]
    set Ξ := 2 * ((H₀ + mk) ^ 2 + 1 / 4) + 2 * mk * (H₀ + mk) ^ 2 + mk / 2
    have hΞ0 : 0 ≤ Ξ := by positivity
    have hrow : ∀ γ ∈ S, (γ ^ 2 + 1 / 4) * Ig γ - 2 * H₀ * mk * (σ / ℓ * Ig γ) -
        ∑ γ' ∈ S, Ξ * (Sg γ * Sg γ') / g₀ ≤ (∑ γ' ∈ S, oubil ℓ (f γ) (f' γ) (f γ') (f' γ')).re := by
      intro γ hγ
      obtain ⟨b1, b2, b3, b4, b5⟩ := hB γ hγ
      rw [← Finset.add_sum_erase S (fun γ' => oubil ℓ (f γ) (f' γ) (f γ') (f' γ')) hγ, Complex.add_re]
      have hdiag := diag_lower (r γ) (γ := γ) hℓ0.le b1 b5
      have hd : (oubil ℓ (f γ) (f' γ) (f γ) (f' γ)).re =
          ∫ u in (0 : ℝ)..ℓ, (‖((r γ).dmod γ).eval u‖ ^ 2 + 1 / 4 * ‖(r γ).eval u‖ ^ 2) := by
        simp only [f, f']; rw [oubil_diag]; simp
      have hdef : 2 * |γ| * ℓ * (Sg γ * (markovC c₀ (r γ).M / ℓ * Sg γ)) ≤
          2 * H₀ * mk * (σ / ℓ * Ig γ) := by
        have e : 2 * |γ| * ℓ * (Sg γ * (markovC c₀ (r γ).M / ℓ * Sg γ)) =
            2 * |γ| * markovC c₀ (r γ).M * Sg γ ^ 2 := by field_simp
        rw [e]
        have h1 := hSg2 γ hγ
        have h2 := hmkg γ hγ
        have h3 := hH γ hγ
        have h4 : 0 ≤ markovC c₀ (r γ).M := markovC_nonneg hc₀ _
        have h5 : 0 ≤ Sg γ ^ 2 := sq_nonneg _
        calc 2 * |γ| * markovC c₀ (r γ).M * Sg γ ^ 2 ≤ 2 * H₀ * mk * Sg γ ^ 2 := by gcongr
          _ ≤ 2 * H₀ * mk * (σ / ℓ * Ig γ) := by gcongr
      have hoff : -(∑ γ' ∈ S.erase γ, Ξ * (Sg γ * Sg γ') / g₀) ≤
          (∑ γ' ∈ S.erase γ, oubil ℓ (f γ) (f' γ) (f γ') (f' γ')).re := by
        rw [Complex.re_sum, ← Finset.sum_neg_distrib]
        refine Finset.sum_le_sum fun γ' hγ' => ?_
        obtain ⟨hne, hγ'S⟩ := Finset.mem_erase.mp hγ'
        obtain ⟨c1, c2, c3, c4, _⟩ := hB γ' hγ'S
        have hcb := cross_bound (r γ) (r γ') hℓ0.le (Ne.symm hne) b3 c3 b4 c4 b1 c1 b2 c2
        have hsep' := hsep γ hγ γ' hγ'S (Ne.symm hne)
        have hss : 0 ≤ Sg γ * Sg γ' := mul_nonneg (hSg0 γ) (hSg0 γ')
        have hnum : 2 * ((H₀ + mk) * Sg γ * ((H₀ + mk) * Sg γ') + Sg γ * Sg γ' / 4) +
            ℓ * (mk / ℓ * ((H₀ + mk) * Sg γ) * ((H₀ + mk) * Sg γ') +
              (H₀ + mk) * Sg γ * (mk / ℓ * ((H₀ + mk) * Sg γ')) +
              (mk / ℓ * Sg γ * Sg γ' + Sg γ * (mk / ℓ * Sg γ')) / 4) = Ξ * (Sg γ * Sg γ') := by
          simp only [Ξ]; field_simp; ring
        rw [hnum] at hcb
        have hpos : 0 < |γ - γ'| := abs_pos.mpr (sub_ne_zero.mpr (Ne.symm hne))
        have hb2 : ‖oubil ℓ (f γ) (f' γ) (f γ') (f' γ')‖ ≤ Ξ * (Sg γ * Sg γ') / g₀ :=
          hcb.trans (div_le_div_of_nonneg_left (by positivity) hg₀ hsep')
        have := (abs_le.mp (Complex.abs_re_le_norm (oubil ℓ (f γ) (f' γ) (f γ') (f' γ')))).1
        linarith
      have hfull : ∑ γ' ∈ S.erase γ, Ξ * (Sg γ * Sg γ') / g₀ ≤ ∑ γ' ∈ S, Ξ * (Sg γ * Sg γ') / g₀ :=
        Finset.sum_le_sum_of_subset_of_nonneg (Finset.erase_subset γ S) fun γ' _ _ => by
          have := mul_nonneg (hSg0 γ) (hSg0 γ'); positivity
      rw [hd]
      have hγ2 : (γ ^ 2 + 1 / 4) * Ig γ - 2 * |γ| * ℓ * (Sg γ * (markovC c₀ (r γ).M / ℓ * Sg γ)) ≤
          ∫ u in (0 : ℝ)..ℓ, (‖((r γ).dmod γ).eval u‖ ^ 2 + 1 / 4 * ‖(r γ).eval u‖ ^ 2) := hdiag
      linarith
    rw [Complex.re_sum]
    have hsum := Finset.sum_le_sum hrow
    have hSS : ∑ γ ∈ S, ∑ γ' ∈ S, Ξ * (Sg γ * Sg γ') / g₀ ≤ Ξ * ((S.card : ℝ) * (σ / ℓ * ∑ γ ∈ S, Ig γ)) / g₀ := by
      have e : ∑ γ ∈ S, ∑ γ' ∈ S, Ξ * (Sg γ * Sg γ') / g₀ = Ξ * (∑ γ ∈ S, Sg γ) ^ 2 / g₀ := by
        rw [sq, Finset.sum_mul_sum, Finset.mul_sum, Finset.sum_div]
        refine Finset.sum_congr rfl fun γ _ => ?_
        rw [Finset.mul_sum, Finset.sum_div]
      rw [e]
      have h1 : (∑ γ ∈ S, Sg γ) ^ 2 ≤ (S.card : ℝ) * ∑ γ ∈ S, Sg γ ^ 2 := sq_sum_le_card_mul_sum_sq
      have h2 : ∑ γ ∈ S, Sg γ ^ 2 ≤ σ / ℓ * ∑ γ ∈ S, Ig γ := by
        rw [Finset.mul_sum]; exact Finset.sum_le_sum hSg2
      have h3 : (∑ γ ∈ S, Sg γ) ^ 2 ≤ (S.card : ℝ) * (σ / ℓ * ∑ γ ∈ S, Ig γ) :=
        h1.trans (by gcongr)
      gcongr
    have hlin : ∑ γ ∈ S, ((γ ^ 2 + 1 / 4) * Ig γ - 2 * H₀ * mk * (σ / ℓ * Ig γ) -
        ∑ γ' ∈ S, Ξ * (Sg γ * Sg γ') / g₀) =
        D - 2 * H₀ * mk * (σ / ℓ) * ∑ γ ∈ S, Ig γ - ∑ γ ∈ S, ∑ γ' ∈ S, Ξ * (Sg γ * Sg γ') / g₀ := by
      rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib, Finset.mul_sum]
      congr 1; congr 1
      refine Finset.sum_congr rfl fun γ _ => ?_; ring
    rw [hlin] at hsum
    have hfin : (1 - claimFC c M₀ H₀ g₀ S.card / ℓ) * D ≤
        D - 2 * H₀ * mk * (σ / ℓ) * ∑ γ ∈ S, Ig γ - Ξ * ((S.card : ℝ) * (σ / ℓ * ∑ γ ∈ S, Ig γ)) / g₀ := by
      have hI4 := hsumI
      have hIs : 0 ≤ ∑ γ ∈ S, Ig γ := Finset.sum_nonneg fun γ _ => hIg0 γ
      have e : (1 - claimFC c M₀ H₀ g₀ S.card / ℓ) * D =
          D - (8 * H₀ * mk * σ + 4 * Ξ * σ * S.card / g₀) / ℓ * D := by
        simp only [claimFC, mk, σ, Ξ]; ring
      rw [e]
      have hk1 : 0 ≤ 2 * H₀ * mk * (σ / ℓ) := by positivity
      have hk2 : 0 ≤ Ξ * (S.card : ℝ) * (σ / ℓ) / g₀ := by positivity
      have e2 : Ξ * ((S.card : ℝ) * (σ / ℓ * ∑ γ ∈ S, Ig γ)) / g₀ =
          Ξ * (S.card : ℝ) * (σ / ℓ) / g₀ * ∑ γ ∈ S, Ig γ := by ring
      have e3 : (8 * H₀ * mk * σ + 4 * Ξ * σ * S.card / g₀) / ℓ * D =
          2 * H₀ * mk * (σ / ℓ) * (4 * D) + Ξ * (S.card : ℝ) * (σ / ℓ) / g₀ * (4 * D) := by
        field_simp; ring
      rw [e2, e3]
      have := mul_le_mul_of_nonneg_left hI4 hk1
      have := mul_le_mul_of_nonneg_left hI4 hk2
      linarith
    linarith
  · -- the endpoint value
    have h0 : (P.part T 0).eval 0 = ∑ γ ∈ S, (r γ).eval 0 := by
      rw [← part_eval_singletons P T 0]
      refine Finset.sum_congr rfl fun γ _ => ?_
      simp only [ofReal_zero, mul_zero, Complex.exp_zero, one_mul]
      rfl
    rw [h0]
    have hval : ∀ γ ∈ S, ‖(r γ).eval 0‖ ^ 2 ≤
        (((r γ).M : ℝ) ^ 2 / (γ ^ 2 + 1 / 4)) * ((γ ^ 2 + 1 / 4) * Ig γ / ℓ) := by
      intro γ hγ
      have := christoffel_expSum (part_singleton_isPoly P γ) hℓ0
      have hq : 0 < γ ^ 2 + 1 / 4 := by positivity
      calc ‖(r γ).eval 0‖ ^ 2 ≤ ((r γ).M : ℝ) ^ 2 / ℓ * Ig γ := this
        _ = _ := by field_simp
    have hCS := norm_sum_sq_le_of_le S (fun γ => (r γ).eval 0)
      (fun γ => ((r γ).M : ℝ) ^ 2 / (γ ^ 2 + 1 / 4)) (fun γ => (γ ^ 2 + 1 / 4) * Ig γ / ℓ)
      (fun γ _ => by positivity) (fun γ _ => by have := hIg0 γ; positivity) hval
    refine hCS.trans (le_of_eq ?_)
    congr 1
    rw [Finset.sum_div]

end Riemann.Deep.LeanNBRate







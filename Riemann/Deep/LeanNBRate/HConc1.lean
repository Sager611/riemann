import Riemann.Deep.LeanNBRate.Dual5
import Mathlib.Analysis.Analytic.Order

/-!
# Interpolation: 1 - E_low vanishes to order m_ρ at the configuration's zeros

The k-th derivative of the Green form E_low(s) = N(e_s, p) is N((-u)^k e_s, p); at a configuration
zero ρ the function (-u)^k e_ρ is a configuration element (k < m), so the representer property gives δ_{k0}.
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem hasDerivAt_laplace {ℓ : ℝ} (hℓ : 0 ≤ ℓ) {q : ℝ → ℂ} (hq : Continuous q) (s₀ : ℂ) :
    HasDerivAt (fun s : ℂ => ∫ u in (0 : ℝ)..ℓ, cexp (-(s - 1 / 2) * u) * q u)
      (∫ u in (0 : ℝ)..ℓ, (-(u : ℂ)) * cexp (-(s₀ - 1 / 2) * u) * q u) s₀ := by
  obtain ⟨M, hM⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := ℓ)).exists_bound_of_continuousOn hq.continuousOn
  set B : ℝ := ℓ * Real.exp ((‖s₀‖ + 2) * ℓ) * max M 0
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume) (a := 0) (b := ℓ)
    (F := fun s u => cexp (-(s - 1 / 2) * u) * q u)
    (F' := fun s u => (-(u : ℂ)) * cexp (-(s - 1 / 2) * u) * q u) (x₀ := s₀)
    (s := Metric.ball s₀ 1) (bound := fun _ => B) (Metric.ball_mem_nhds s₀ one_pos)
    (Filter.Eventually.of_forall fun s => (by fun_prop : Continuous fun u : ℝ =>
      cexp (-(s - 1 / 2) * u) * q u).aestronglyMeasurable)
    ((by fun_prop : Continuous fun u : ℝ => cexp (-(s₀ - 1 / 2) * u) * q u).intervalIntegrable _ _)
    ((by fun_prop : Continuous fun u : ℝ => (-(u : ℂ)) * cexp (-(s₀ - 1 / 2) * u) * q u).aestronglyMeasurable)
    ?_ intervalIntegrable_const ?_
  · exact key.2
  · refine Filter.Eventually.of_forall fun u hu s hs => ?_
    rw [Set.uIoc_of_le hℓ] at hu
    have hu0 : 0 ≤ u := hu.1.le
    have hul : u ≤ ℓ := hu.2
    have hs' : ‖s‖ ≤ ‖s₀‖ + 1 := by
      have := Metric.mem_ball.mp hs
      have h2 := norm_le_norm_add_norm_sub' s s₀
      rw [dist_eq_norm] at this
      linarith
    rw [norm_mul, norm_mul, norm_neg, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hu0, Complex.norm_exp]
    have hre : (-(s - 1 / 2) * (u : ℂ)).re ≤ (‖s₀‖ + 2) * ℓ := by
      have : (-(s - 1 / 2) * (u : ℂ)).re = -(s.re - 1 / 2) * u := by simp
      rw [this]
      have h1 : |s.re| ≤ ‖s‖ := Complex.abs_re_le_norm s
      have h2 : -(s.re - 1 / 2) ≤ ‖s₀‖ + 2 := by have := neg_abs_le s.re; linarith
      calc -(s.re - 1 / 2) * u ≤ (‖s₀‖ + 2) * u := mul_le_mul_of_nonneg_right h2 hu0
        _ ≤ (‖s₀‖ + 2) * ℓ := mul_le_mul_of_nonneg_left hul (by positivity)
    have hqu : ‖q u‖ ≤ max M 0 := (hM u ⟨hu0, hul⟩).trans (le_max_left _ _)
    have hfin : u * Real.exp (-(s - 1 / 2) * (u : ℂ)).re * ‖q u‖ ≤
        ℓ * Real.exp ((‖s₀‖ + 2) * ℓ) * max M 0 := by
      gcongr
    exact hfin
  · refine Filter.Eventually.of_forall fun u _ s _ => ?_
    have h1 : HasDerivAt (fun y : ℂ => -(y - 1 / 2) * (u : ℂ)) (-(u : ℂ)) s := by
      have := ((hasDerivAt_id s).sub_const (1 / 2 : ℂ)).const_mul (-(u : ℂ))
      convert this using 1
      · ext y; simp only [id]; ring
      · ring
    have h := (h1.cexp).mul_const (q u)
    convert h using 1
    ring

/-- The k-th derivative of the Green form. -/
def Ek (ℓ : ℝ) (q : ℝ → ℂ) (A B : ℂ) (k : ℕ) (s : ℂ) : ℂ :=
  (∫ u in (0 : ℝ)..ℓ, cexp (-(s - 1 / 2) * u) * ((-(u : ℂ)) ^ k * q u)) + (if k = 0 then A else 0) +
    (-(ℓ : ℂ)) ^ k * cexp (-(s - 1 / 2) * ℓ) * B

theorem hasDerivAt_Ek {ℓ : ℝ} (hℓ : 0 ≤ ℓ) {q : ℝ → ℂ} (hq : Continuous q) (A B : ℂ) (k : ℕ) (s : ℂ) :
    HasDerivAt (Ek ℓ q A B k) (Ek ℓ q A B (k + 1) s) s := by
  have h1 := hasDerivAt_laplace hℓ (q := fun u => (-(u : ℂ)) ^ k * q u) (by fun_prop) s
  have h2 : HasDerivAt (fun s : ℂ => (-(ℓ : ℂ)) ^ k * cexp (-(s - 1 / 2) * ℓ) * B)
      ((-(ℓ : ℂ)) ^ (k + 1) * cexp (-(s - 1 / 2) * ℓ) * B) s := by
    have h3 : HasDerivAt (fun y : ℂ => -(y - 1 / 2) * (ℓ : ℂ)) (-(ℓ : ℂ)) s := by
      have := ((hasDerivAt_id s).sub_const (1 / 2 : ℂ)).const_mul (-(ℓ : ℂ))
      convert this using 1
      · ext y; simp only [id]; ring
      · ring
    have h4 := ((h3.cexp).const_mul ((-(ℓ : ℂ)) ^ k)).mul_const B
    convert h4 using 1
    ring
  have h := (h1.add_const (if k = 0 then A else 0)).add h2
  unfold Ek
  convert h using 1
  rw [if_neg (Nat.succ_ne_zero k), add_zero]
  congr 1
  refine intervalIntegral.integral_congr fun u _ => ?_
  ring

theorem iteratedDeriv_one_sub_Ek {ℓ : ℝ} (hℓ : 0 ≤ ℓ) {q : ℝ → ℂ} (hq : Continuous q) (A B : ℂ) (k : ℕ) :
    iteratedDeriv k (fun s => 1 - Ek ℓ q A B 0 s) = fun s => (if k = 0 then 1 else 0) - Ek ℓ q A B k s := by
  induction k with
  | zero => funext s; simp
  | succ k ih =>
    rw [iteratedDeriv_succ, ih]
    funext s
    rw [((hasDerivAt_Ek hℓ hq A B k s).const_sub (if k = 0 then (1 : ℂ) else 0)).deriv,
      if_neg (Nat.succ_ne_zero k)]
    ring

theorem Elow_eq_Ek {ℓ : ℝ} (p p' p'' : ℝ → ℂ) (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcp : Continuous p)
    (hcp' : Continuous p') (hcp'' : Continuous p'') :
    Elow ℓ p p' = Ek ℓ (fun u => conj (gDens p p' p'' u)) (conj (gA p p')) (conj (gB ℓ p p')) 0 := by
  funext s
  rw [Elow_green p p' p'' hp' hcp hcp' hcp'' s]
  unfold Ek
  simp

theorem Ek_eq_NForm {ℓ : ℝ} (p p' p'' : ℝ → ℂ) (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcp : Continuous p)
    (hcp' : Continuous p') (hcp'' : Continuous p'') {f f' : ℝ → ℂ} (hf : ∀ x, HasDerivAt f (f' x) x)
    (hcf : Continuous f) (hcf' : Continuous f') {ρ : ℂ} {k : ℕ}
    (hfe : ∀ u : ℝ, f u = (-(u : ℂ)) ^ k * cexp (-(ρ - 1 / 2) * u)) :
    Ek ℓ (fun u => conj (gDens p p' p'' u)) (conj (gA p p')) (conj (gB ℓ p p')) k ρ = NForm ℓ f f' p p' := by
  rw [NForm_green f f' p p' p'' hf hp' hcf hcf' hcp hcp' hcp'']
  unfold Ek
  have h0 : f 0 = if k = 0 then 1 else 0 := by
    rw [hfe]
    simp only [Complex.ofReal_zero, neg_zero, mul_zero, Complex.exp_zero, mul_one]
    rcases Nat.eq_zero_or_pos k with h | h
    · subst h; simp
    · rw [zero_pow (by omega), if_neg (by omega)]
  rw [h0, hfe ℓ]
  have hint : ∫ u in (0 : ℝ)..ℓ, cexp (-(ρ - 1 / 2) * u) * ((-(u : ℂ)) ^ k * conj (gDens p p' p'' u)) =
      ∫ u in (0 : ℝ)..ℓ, f u * conj (gDens p p' p'' u) := by
    refine intervalIntegral.integral_congr fun u _ => ?_
    rw [hfe]; ring
  rw [hint]
  have hif : (if k = 0 then conj (gA p p') else 0) = (if k = 0 then (1 : ℂ) else 0) * conj (gA p p') := by
    split_ifs <;> simp
  rw [hif]

/-- The single-coefficient vector: I^k at (j₀, k). -/
def sc (P : ExpSum) (j₀ : Fin P.J) (k : ℕ) : Fin P.J → ℕ → ℂ :=
  fun j k' => if j = j₀ ∧ k' = k then I ^ k else 0

theorem eval_sc (P : ExpSum) (j₀ : Fin P.J) {k : ℕ} (hk : k < P.mult j₀) (u : ℝ) :
    (P.withCoef (sc P j₀ k)).eval u = (-(u : ℂ)) ^ k * cexp (I * P.node j₀ * u) := by
  unfold ExpSum.eval ExpSum.withCoef sc
  simp only
  rw [Finset.sum_eq_single j₀]
  · rw [Finset.sum_eq_single k]
    · simp only [and_self, if_true]
      have e : I ^ k * (I * (u : ℂ)) ^ k = (-(u : ℂ)) ^ k := by
        rw [← mul_pow]; congr 1; rw [← mul_assoc, I_mul_I]; ring
      rw [e]
    · intro k' _ hk'
      simp [hk']
    · intro h
      exact absurd (Finset.mem_range.mpr hk) h
  · intro j _ hj
    apply Finset.sum_eq_zero
    intro k' _
    simp [hj]
  · intro h
    exact absurd (Finset.mem_univ j₀) h

/-- **Interpolation**: all derivatives of order < m of 1 - E_low vanish at a configuration zero. -/
theorem iteratedDeriv_one_sub_Elow {P : ExpSum} {ℓ : ℝ} (hℓ : 0 ≤ ℓ) {c : Fin P.J → ℕ → ℂ}
    (hc : ∀ d : Fin P.J → ℕ → ℂ, NForm ℓ (P.withCoef d).eval ((P.withCoef d).dmod 0).eval
      (P.withCoef c).eval ((P.withCoef c).dmod 0).eval = (P.withCoef d).eval 0)
    {j₀ : Fin P.J} {k : ℕ} (hk : k < P.mult j₀) {ρ : ℂ} (hρ : ρ - 1 / 2 = -((P.node j₀ : ℝ) : ℂ) * I) :
    iteratedDeriv k (fun s => 1 - Elow ℓ (P.withCoef c).eval ((P.withCoef c).dmod 0).eval s) ρ = 0 := by
  obtain ⟨r, hr⟩ : ∃ r : ExpSum, r = P.withCoef c := ⟨_, rfl⟩
  rw [← hr]
  rw [← hr] at hc
  have hp' : ∀ x, HasDerivAt (r.dmod 0).eval (((r.dmod 0).dmod 0).eval x) x := fun x => (r.dmod 0).hasDerivAt_eval x
  have hcp : Continuous r.eval := r.continuous_eval
  have hcp' : Continuous (r.dmod 0).eval := (r.dmod 0).continuous_eval
  have hcp'' : Continuous ((r.dmod 0).dmod 0).eval := ((r.dmod 0).dmod 0).continuous_eval
  have hq : Continuous fun u => conj (gDens r.eval (r.dmod 0).eval ((r.dmod 0).dmod 0).eval u) := by
    unfold gDens; fun_prop
  rw [Elow_eq_Ek r.eval (r.dmod 0).eval ((r.dmod 0).dmod 0).eval hp' hcp hcp' hcp'',
    iteratedDeriv_one_sub_Ek hℓ hq]
  simp only
  obtain ⟨d, hd⟩ : ∃ d, d = sc P j₀ k := ⟨_, rfl⟩
  have hfe : ∀ u : ℝ, (P.withCoef d).eval u = (-(u : ℂ)) ^ k * cexp (-(ρ - 1 / 2) * u) := by
    intro u; rw [hd, eval_sc P j₀ hk u, hρ]; congr 2; ring
  rw [Ek_eq_NForm r.eval (r.dmod 0).eval ((r.dmod 0).dmod 0).eval hp' hcp hcp' hcp''
    (P.withCoef d).hasDerivAt_eval (P.withCoef d).continuous_eval ((P.withCoef d).dmod 0).continuous_eval hfe,
    hc d, hfe 0]
  simp only [Complex.ofReal_zero, neg_zero, mul_zero, Complex.exp_zero, mul_one]
  rcases Nat.eq_zero_or_pos k with h | h
  · subst h; simp
  · rw [zero_pow (by omega), if_neg (by omega), sub_zero]

end Riemann.Deep.LeanNBRate


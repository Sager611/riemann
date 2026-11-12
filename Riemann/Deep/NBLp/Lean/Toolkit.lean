import Riemann.Deep.NBLp.Lean.PropN
import Mathlib.Analysis.SpecialFunctions.Gamma.Deligne

/-!
# NBLp in Lean, part 3: the §4–§5 toolkit as named hypotheses; Lemma M proved

Named hypotheses (PROOF.md line numbers in LEAN.md):
* MellinToolkit HM — §4 Mellin–Plancherel (HM σ F G: F ∈ L²_σ and 𝓜_σ F = G a.e.).
* LemmaT HM, LemmaPW HM — §4 Lemma T, Lemma PW.
* LemmaGb, LemmaH, LemmaA, LemmaP — §5 Lemma Γ(b), Lemma H, Lemma A, Lemma P.
Proved: lemmaM (§4 Lemma M, Hölder).
-/

open Complex MeasureTheory Set Filter Topology
open scoped ENNReal

namespace NBLp

/-- §0: squared norm in L²_σ = L²((0,∞), x^{2σ−1} dx). -/
noncomputable def L2w (σ : ℝ) (F : ℝ → ℂ) : ℝ := ∫ x in Ioi (0:ℝ), ‖F x‖ ^ 2 * x ^ (2 * σ - 1)

/-- §0: F ∈ L²_σ. -/
def MemL2w (σ : ℝ) (F : ℝ → ℂ) : Prop :=
  AEStronglyMeasurable F (volume.restrict (Ioi 0)) ∧
    IntegrableOn (fun x => ‖F x‖ ^ 2 * x ^ (2 * σ - 1)) (Ioi 0)

/-- §0: χ₊ = 1_(0,1] on (0,∞). -/
noncomputable def chiP (x : ℝ) : ℂ := if x ≤ 1 then 1 else 0

/-- §4 toolkit (named hypothesis; [C7] Plancherel). HM σ F G: F ∈ L²_σ and 𝓜_σ F = G a.e.,
where 𝓜_σ F = (u ↦ F(e^{−u})e^{−σu})^∧ (Fourier–Plancherel). Satisfied by
HM σ F G := MemL2w σ F ∧ G =ᵐ 𝓜_σ F. -/
structure MellinToolkit (HM : ℝ → (ℝ → ℂ) → (ℝ → ℂ) → Prop) : Prop where
  mem : ∀ {σ F G}, HM σ F G → MemL2w σ F
  plancherel : ∀ {σ F G}, HM σ F G →
    Integrable (fun t => ‖G t‖ ^ 2) ∧ L2w σ F = (2 * Real.pi)⁻¹ * ∫ t, ‖G t‖ ^ 2
  add : ∀ {σ F G F' G'}, HM σ F G → HM σ F' G' → HM σ (F + F') (G + G')
  smul : ∀ {σ F G} (c : ℂ), HM σ F G → HM σ (c • F) (c • G)
  congr : ∀ {σ F G F' G'}, HM σ F G → F =ᵐ[volume.restrict (Ioi 0)] F' → G =ᵐ[volume] G' →
    HM σ F' G'
  l1 : ∀ {σ F}, MemL2w σ F →
    IntegrableOn (fun x : ℝ => (x : ℂ) ^ ((σ : ℂ) - 1) * F x) (Ioi 0) →
    HM σ F (fun t => ∫ x in Ioi (0:ℝ), (x : ℂ) ^ ((σ : ℂ) + t * I - 1) * F x)

/-- §4 Lemma T (named hypothesis). -/
def LemmaT (HM : ℝ → (ℝ → ℂ) → (ℝ → ℂ) → Prop) : Prop :=
  (∀ σ : ℝ, 0 < σ → σ < 1 → ∀ k : ℕ, 1 ≤ k →
    IntegrableOn (fun x : ℝ => x ^ (σ - 1) * rho k x) (Ioi 0) ∧
    MemL2w σ (fun x => (rho k x : ℂ)) ∧
    L2w σ (fun x => (rho k x : ℂ)) = (k : ℝ) ^ (-2 * σ) * L2w σ (fun x => (rho 1 x : ℂ)) ∧
    L2w σ (fun x => (rho 1 x : ℂ)) ≤ 1 / (2 * σ) + 1 / (2 - 2 * σ) ∧
    HM σ (fun x => (rho k x : ℂ))
      (fun t => -((k : ℂ) ^ (-((σ : ℂ) + t * I))) * riemannZeta (σ + t * I) / (σ + t * I))) ∧
  (∀ σ : ℝ, 0 < σ → HM σ chiP (fun t => 1 / ((σ : ℂ) + t * I)))

/-- §4 Lemma PW (named hypothesis). -/
def LemmaPW (HM : ℝ → (ℝ → ℂ) → (ℝ → ℂ) → Prop) : Prop :=
  ∀ (a b : ℝ) (Ψ : ℂ → ℂ) (U : Set ℂ), a < b → IsOpen U →
    {s : ℂ | a ≤ s.re ∧ s.re ≤ b} ⊆ U → DifferentiableOn ℂ Ψ U →
    (∃ M : ℝ, ∀ σ ∈ Icc a b, Integrable (fun t : ℝ => ‖Ψ (σ + t * I)‖ ^ 2) ∧
      ∫ t : ℝ, ‖Ψ (σ + t * I)‖ ^ 2 ≤ M) →
    ∀ F : ℝ → ℂ, HM b F (fun t => Ψ (b + t * I)) → HM a F (fun t => Ψ (a + t * I))

/-- §0: γ(s) = ½ s(s−1) π^{−s/2} Γ(s/2). -/
noncomputable def gam (s : ℂ) : ℂ := 1 / 2 * s * (s - 1) * Gammaℝ s

/-- §0: Γr_ε(s) = γ(s+ε)/γ(s). -/
noncomputable def Gr (ε : ℝ) (s : ℂ) : ℂ := gam (s + ε) / gam s

/-- §0: ξ(s) = γ(s)ζ(s), entire: ½s(s−1)Λ₀(s) + ½ (Λ = Λ₀ − 1/s − 1/(1−s)). -/
noncomputable def xi (s : ℂ) : ℂ := 1 / 2 * s * (s - 1) * completedRiemannZeta₀ s + 1 / 2

/-- §5 Lemma Γ(b) (named hypothesis). -/
def LemmaGb : Prop :=
  ∀ (ε : ℝ) (s : ℂ), 0 < ε → ε ≤ 1 → 1 / 2 ≤ s.re →
    (∀ u ∈ Icc (0:ℝ) ε, s.re + u < 1 ∨ 1 ≤ |s.im|) →
    ‖Gr ε s‖ ≤ Real.exp (5 * ε) * (2 + ‖s‖) ^ (ε / 2) ∧
    (s.re + ε < 1 → ‖Gr ε s‖ ^ 2 ≤ Real.exp (8 * ε) * (2 + |s.im|) ^ ε)

/-- §5 Lemma H (named hypothesis). -/
def LemmaH : Prop :=
  ∀ (θs ε : ℝ) (s : ℂ), (∀ ρ : ℂ, riemannZeta ρ = 0 → 0 < ρ.re → ρ.re ≤ θs) →
    0 < ε → θs - ε / 2 ≤ s.re → s ≠ 1 → 0 < s.re →
    ‖xi s‖ ≤ ‖xi (s + ε)‖ ∧
      (s + ε ≠ 1 → ‖riemannZeta s‖ ≤ ‖riemannZeta (s + ε)‖ * ‖Gr ε s‖)

/-- §5 Lemma A (named hypothesis). -/
def LemmaA : Prop :=
  ∀ θ : ℝ, 1 / 2 ≤ θ → θ < 1 → ZeroFree θ → ∀ δ η' : ℝ, 0 < δ → δ ≤ (1 - θ) / 2 → 0 < η' →
    ∃ C : ℝ, ∀ s : ℂ, θ + 2 * δ ≤ s.re →
      ‖riemannZeta s‖⁻¹ ≤ C * (2 + |s.im|) ^ η' ∧
        (1 ≤ |s.im| → ‖riemannZeta s‖ ≤ C * (2 + |s.im|) ^ η')

/-- §0: P_N(s) = Σ_{k≤N} μ(k)(1 − k/N) k^{−s}. -/
noncomputable def PN (N : ℕ) (s : ℂ) : ℂ :=
  ∑ k ∈ Finset.Icc 1 N, (ArithmeticFunction.moebius k : ℂ) * (1 - (k : ℂ) / N) * (k : ℂ) ^ (-s)

/-- 1/ζ with 1/ζ(1) := 0 (Lemma P). -/
noncomputable def invZeta (s : ℂ) : ℂ := if s = 1 then 0 else (riemannZeta s)⁻¹

/-- §5 Lemma P (named hypothesis). -/
def LemmaP : Prop :=
  (∀ θ : ℝ, 1 / 2 ≤ θ → θ < 1 → ZeroFree θ → ∀ η δ η' : ℝ, 0 < η → η < 1 → 0 < δ →
    δ ≤ (1 - θ) / 2 → 0 < η' → η' ≤ 1 / 2 →
    ∃ C : ℝ, ∀ s : ℂ, θ + 2 * δ ≤ s.re - η → ∀ N : ℕ, 1 ≤ N →
      ‖PN N s - invZeta s‖ ≤ C * (N : ℝ) ^ (-η) * (2 + |s.im|) ^ η') ∧
  (∀ η : ℝ, 0 < η → η < 1 → ∀ s : ℂ, 1 < s.re - η → ∀ N : ℕ, 1 ≤ N →
      ‖PN N s - invZeta s‖ ≤ (1 / (1 - η) + 1 / η) * (N : ℝ) ^ (-η))

/-- §4 Lemma M (Hölder), proved: ‖F‖_{L^p(0,1)} ≤ (1 − r(a − ½))^{−1/r} ‖x^{a−½}F‖_{L²(0,1)}. -/
theorem lemmaM {p r a : ℝ} (hp1 : 1 ≤ p) (hr0 : 0 < r) (hr : 1 / p = 1 / 2 + 1 / r)
    (ha : a < 1 / p) {F : ℝ → ℂ} (hF : AEStronglyMeasurable F vol01) :
    eLpNorm F (ENNReal.ofReal p) vol01 ≤
      ENNReal.ofReal ((1 - r * (a - 1 / 2)) ^ (-1 / r)) *
        eLpNorm (fun x => ((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x) 2 vol01 := by
  have hp0 : 0 < p := by linarith
  set c : ℝ := 1 / 2 - a with hc
  have hrc : -1 < r * c := by
    have e1 : 1 / 2 - 1 / p = -(1 / r) := by linarith
    have e2 : r * (1 / r) = 1 := by field_simp
    have h3 : r * (1 / 2 - 1 / p) < r * c := mul_lt_mul_of_pos_left (by rw [hc]; linarith) hr0
    rw [e1, mul_neg, e2] at h3
    exact h3
  set φ : ℝ → ℂ := fun x => ((x ^ c : ℝ) : ℂ) with hφ
  set f : ℝ → ℂ := fun x => ((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x with hf
  have hφm : AEStronglyMeasurable φ vol01 :=
    (Complex.measurable_ofReal.comp (measurable_id.pow_const c)).aestronglyMeasurable
  have hfm : AEStronglyMeasurable f vol01 :=
    (Complex.measurable_ofReal.comp (measurable_id.pow_const (a - 1 / 2))).aestronglyMeasurable.mul
      hF
  have hFeq : F =ᵐ[vol01] φ • f := by
    refine (ae_restrict_iff' measurableSet_Ioc).mpr (Eventually.of_forall fun x hx => ?_)
    show F x = φ x * f x
    simp only [hφ, hf]
    have e0 : c + (a - 1 / 2) = 0 := by rw [hc]; ring
    rw [← mul_assoc, ← Complex.ofReal_mul, ← Real.rpow_add hx.1, e0, Real.rpow_zero,
      Complex.ofReal_one, one_mul]
  have : ENNReal.HolderTriple (ENNReal.ofReal r) 2 (ENNReal.ofReal p) := ⟨by
    rw [← ENNReal.ofReal_inv_of_pos hr0, ← ENNReal.ofReal_inv_of_pos hp0,
      show (2:ℝ≥0∞)⁻¹ = ENNReal.ofReal 2⁻¹ by rw [ENNReal.ofReal_inv_of_pos two_pos]; simp,
      ← ENNReal.ofReal_add (by positivity) (by positivity)]
    congr 1
    simp only [← one_div]
    linarith [hr]⟩
  have hint : IntegrableOn (fun x : ℝ => x ^ (r * c)) (Ioc 0 1) := by
    have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := 1) (r := r * c) hrc
    rwa [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
  have hval : ∫ x in Ioc (0:ℝ) 1, x ^ (r * c) = 1 / (r * c + 1) := by
    rw [← intervalIntegral.integral_of_le zero_le_one, integral_rpow (Or.inl hrc)]
    simp [Real.zero_rpow (by linarith : r * c + 1 ≠ 0)]
  have hφnorm : eLpNorm φ (ENNReal.ofReal r) vol01 =
      ENNReal.ofReal ((1 - r * (a - 1 / 2)) ^ (-1 / r)) := by
    rw [eLpNorm_eq_lintegral_rpow_enorm_toReal (by simp [hr0]) ENNReal.ofReal_ne_top hφm,
      ENNReal.toReal_ofReal hr0.le]
    have heq : ∀ x ∈ Ioc (0:ℝ) 1, ‖φ x‖ₑ ^ r = ENNReal.ofReal (x ^ (r * c)) := by
      intro x hx
      rw [← ofReal_norm, hφ]
      simp only [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hx.1.le c)]
      rw [ENNReal.ofReal_rpow_of_nonneg (Real.rpow_nonneg hx.1.le c) hr0.le,
        ← Real.rpow_mul hx.1.le, mul_comm]
    rw [setLIntegral_congr_fun measurableSet_Ioc heq,
      ← ofReal_integral_eq_lintegral_ofReal hint
        (ae_restrict_of_forall_mem measurableSet_Ioc fun x hx => Real.rpow_nonneg hx.1.le _),
      hval, ENNReal.ofReal_rpow_of_nonneg (div_pos one_pos (by linarith)).le
        (div_pos one_pos hr0).le]
    congr 1
    have e : r * c + 1 = 1 - r * (a - 1 / 2) := by rw [hc]; ring
    rw [e, one_div, Real.inv_rpow (by linarith), ← Real.rpow_neg (by linarith), neg_div]
  calc eLpNorm F (ENNReal.ofReal p) vol01 = eLpNorm (φ • f) (ENNReal.ofReal p) vol01 :=
        eLpNorm_congr_ae hFeq
    _ ≤ eLpNorm φ (ENNReal.ofReal r) vol01 * eLpNorm f 2 vol01 :=
        eLpNorm_smul_le_mul_eLpNorm hφm hfm
    _ = _ := by rw [hφnorm]

end NBLp

#print axioms NBLp.lemmaM



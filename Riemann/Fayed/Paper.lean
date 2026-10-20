import Riemann.Fayed.Theorem1

/-!
# Fayed, arXiv:2209.01890v55 — the whole argument in Lean

s = σ + it, 0 < σ < 1; p ∈ (0,1); η ∈ (0, min(σ,1-σ)); γ ∈ (0,1/4); T = u⁸.

* Literature facts the paper cites (all true, standard) → fields of Inputs.
* Paper-internal steps → lemmas eqNN, numbered as in the paper.
  Algebraic identities: proved. Asymptotic steps (39),(42),(50),(54): true, stated,
  proofs pending (sorry, marked HARDEN).
* fayed_RH : Inputs → Theorem1C → RiemannHypothesis — the paper's logic, checked.
* theorem1C_false : ¬ Theorem1C — the argument breaks exactly at (55).
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real ComplexConjugate

namespace Fayed

noncomputable section

/-- B̄ₖ(x) = Bₖ({x}). -/
def Bbar (k : ℕ) (x : ℝ) : ℝ := Polynomial.aeval (Int.fract x) (Polynomial.bernoulli k)

/-- ∑_{1 ≤ n ≤ u} n^{-s}. -/
def S (u : ℝ) (s : ℂ) : ℂ := ∑ n ∈ Finset.Icc 1 ⌊u⌋₊, (n : ℂ) ^ (-s)

/-- Line integrand of (33): u^{1-σ-η+iy} ζ(1-η+i(t+y)) / (1-σ-η+iy). -/
def Φ (η u : ℝ) (s : ℂ) (y : ℝ) : ℂ :=
  (u : ℂ) ^ (((1 - s.re - η : ℝ) : ℂ) + y * I) *
    riemannZeta (((1 - η : ℝ) : ℂ) + (s.im + y) * I) / (((1 - s.re - η : ℝ) : ℂ) + y * I)

/-- (33): J(u,s) = (1/2π) ∫_{-u⁸}^{u⁸} Φ dy. -/
def J (η u : ℝ) (s : ℂ) : ℂ := ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-(u ^ 8))..(u ^ 8), Φ η u s y

/-- Euler–Maclaurin correction of (17)/(34): B̄₁u^{-s} + sB̄₂/(2u^{s+1}) + s(s+1)B̄₃/(6u^{s+2}). -/
def EMcorr (u : ℝ) (s : ℂ) : ℂ :=
  (Bbar 1 u : ℂ) * (u : ℂ) ^ (-s) + s * (Bbar 2 u : ℂ) / (2 * (u : ℂ) ^ (s + 1)) +
    s * (s + 1) * (Bbar 3 u : ℂ) / (6 * (u : ℂ) ^ (s + 2))

/-- (36): I(u,s) = u^{s+p} ζ(s) - u^{1-s̄+p} ζ(1-s̄). -/
def Ifun (p u : ℝ) (s : ℂ) : ℂ :=
  (u : ℂ) ^ (s + p) * riemannZeta s - (u : ℂ) ^ (1 - conj s + p) * riemannZeta (1 - conj s)

/-- Integrand of F, (41). -/
def ψ (η σ₀ t₀ u y : ℝ) : ℂ :=
  (1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)) *
    (u : ℂ) ^ (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I)

/-- (41): F(u) = (1/2π) ∫_{-u⁸}^{u⁸} ψ dy. -/
def Ffun (η σ₀ t₀ u : ℝ) : ℂ := ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-(u ^ 8))..(u ^ 8), ψ η σ₀ t₀ u y

/-- (41): G(u) = u^{-p}. -/
def Gfun (p u : ℝ) : ℂ := ((u ^ (-p) : ℝ) : ℂ)

/-- (44): Leibniz form of F'(u) (integrand derivative + both moving endpoints). -/
def Fder (η σ₀ t₀ u : ℝ) : ℂ :=
  ((1 / (2 * π) : ℝ) : ℂ) * ((∫ y in (-(u ^ 8))..(u ^ 8),
      (1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)) *
      ((((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * (u : ℂ) ^ (((-η : ℝ) : ℂ) + (t₀ + y) * I)) *
      riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I)) +
    8 * (u : ℂ) ^ 7 * (ψ η σ₀ t₀ u (u ^ 8) + ψ η σ₀ t₀ u (-(u ^ 8))))

/-- Literature facts used by the paper. All are true theorems. -/
structure Inputs : Prop where
  /-- (17): Euler–Maclaurin, k = 4 [Montgomery–Vaughan; Ivić]. -/
  em17 : ∀ s : ℂ, 0 < s.re → s.re < 1 → ∃ C : ℝ, ∀ u : ℝ, 1 ≤ u →
    ‖S u s - (riemannZeta s + (u : ℂ) ^ (1 - s) / (1 - s) - EMcorr u s)‖ ≤ C * u ^ (-(s.re + 3))
  /-- (25)–(33): truncated Perron (c = 2, T = u⁸) + contour shift to Re z = 1-σ-η
  [Apostol Lemma 11.12; convexity bound for the horizontal sides]. -/
  perron33 : ∀ s : ℂ, 0 < s.re → s.re < 1 → ∀ η : ℝ, 0 < η → η < s.re → η < 1 - s.re →
    ∀ γ : ℝ, 0 < γ → γ < 1 / 4 → ∀ ε : ℝ, 0 < ε → ∃ C U : ℝ, ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖S u s - (J η u s + (u : ℂ) ^ (1 - s) / (1 - s))‖ ≤ C * u ^ (-4 + ε)
  /-- (44): Leibniz rule (integrand smooth in (u,y) for u > 0). -/
  leibniz44 : ∀ η σ₀ t₀ : ℝ, 0 < η → η < σ₀ → η < 1 - σ₀ → ∀ u : ℝ, 0 < u →
    HasDerivAt (Ffun η σ₀ t₀) (Fder η σ₀ t₀ u) u
  /-- (30),(45): convexity bound on Re s = 1-η [Titchmarsh–Heath-Brown §5.1]. -/
  convexity : ∀ η : ℝ, 0 < η → η < 1 / 2 → ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ T : ℝ, 1 ≤ |T| →
    ‖riemannZeta (((1 - η : ℝ) : ℂ) + T * I)‖ ≤ C * |T| ^ (η / 2 + ε)
  /-- §2: nontrivial zeros lie in 0 < Re s < 1 [Hadamard; de la Vallée Poussin]. -/
  strip : ∀ s : ℂ, riemannZeta s = 0 → (¬∃ n : ℕ, s = -2 * (n + 1)) → s ≠ 1 →
    0 < s.re ∧ s.re < 1
  /-- §2: ζ(s) = 0 in the strip ⇒ ζ(1 - s̄) = 0 (functional equation (8) + conjugation). -/
  reflect : ∀ s : ℂ, riemannZeta s = 0 → 0 < s.re → s.re < 1 → riemannZeta (1 - conj s) = 0

/-- §3 at a hypothetical zero s₀ = σ₀ + it₀ with parameters p, η, γ. -/
structure Setting where
  σ₀ : ℝ
  t₀ : ℝ
  p : ℝ
  η : ℝ
  γ : ℝ
  hσ₀ : 0 < σ₀ ∧ σ₀ < 1
  hzero : riemannZeta ((σ₀ : ℂ) + t₀ * I) = 0
  hp : 0 < p ∧ p < 1
  hη : 0 < η ∧ η < σ₀ ∧ η < 1 - σ₀
  hγ : 0 < γ ∧ γ < 1 / 4

/-! ## Algebraic identities of the paper (proved) -/

/-- (37): B̄₂ coefficient, s - (1 - s̄) = 2σ - 1. -/
lemma eq37_B2 (s : ℂ) : s - (1 - conj s) = 2 * (s.re : ℂ) - 1 := by
  apply Complex.ext <;> simp <;> ring

/-- (38): B̄₃ coefficient, s(s+1) - (1-s̄)(2-s̄) = 2(2σ-1)(1+it). -/
lemma eq38_B3 (s : ℂ) :
    s * (s + 1) - (1 - conj s) * (2 - conj s) = 2 * (2 * (s.re : ℂ) - 1) * (1 + s.im * I) := by
  apply Complex.ext <;> simp [Complex.mul_re, Complex.mul_im] <;> ring

/-- (48): B̄₁ coefficient, (1 - s̄) - s = 1 - 2σ. -/
lemma eq48_B1 (s : ℂ) : (1 - conj s) - s = 1 - 2 * (s.re : ℂ) := by
  apply Complex.ext <;> simp <;> ring

/-- (37): the B̄₁ terms of u^{s+p}ζ(s) and u^{1-s̄+p}ζ(1-s̄) cancel. -/
lemma eq37_B1_cancel {u : ℝ} (hu : 0 < u) (p : ℝ) (s x : ℂ) :
    (u : ℂ) ^ (s + p) * (x * (u : ℂ) ^ (-s)) -
      (u : ℂ) ^ (1 - conj s + p) * (x * (u : ℂ) ^ (-(1 - conj s))) = 0 := by
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  have e1 : (u : ℂ) ^ (s + p) * (u : ℂ) ^ (-s) = (u : ℂ) ^ (p : ℂ) := by
    rw [← Complex.cpow_add _ _ hu']; ring_nf
  have e2 : (u : ℂ) ^ (1 - conj s + p) * (u : ℂ) ^ (-(1 - conj s)) = (u : ℂ) ^ (p : ℂ) := by
    rw [← Complex.cpow_add _ _ hu']; ring_nf
  linear_combination x * e1 - x * e2

/-- w(1/A - 1/B) = a/A - b/B when w = a + A = b + B. -/
lemma key_alg {K : Type*} [Field K] {w a b A B : K} (hA : A ≠ 0) (hB : B ≠ 0)
    (h1 : w = a + A) (h2 : w = b + B) : w * (1 / A - 1 / B) = a / A - b / B := by
  have e1 : w / A = a / A + 1 := by rw [h1, add_div, div_self hA]
  have e2 : w / B = b / B + 1 := by rw [h2, add_div, div_self hB]
  calc w * (1 / A - 1 / B) = w / A - w / B := by ring
    _ = a / A - b / B := by rw [e1, e2]; ring

/-- (44), second line: (1-η+i(t₀+y)) [1/(1-σ₀-η+iy) - 1/(σ₀-η+iy)]
 = s₀/(1-σ₀-η+iy) - (1-s̄₀)/(σ₀-η+iy). -/
lemma eq44_algebra {σ₀ t₀ η : ℝ} (h1 : η < σ₀) (h2 : η < 1 - σ₀) (y : ℝ) :
    (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) *
        (1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)) =
      ((σ₀ : ℂ) + t₀ * I) / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) -
        (1 - conj ((σ₀ : ℂ) + t₀ * I)) / (((σ₀ - η : ℝ) : ℂ) + y * I) := by
  have hA : ((1 - σ₀ - η : ℝ) : ℂ) + y * I ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  have hB : ((σ₀ - η : ℝ) : ℂ) + y * I ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  refine key_alg hA hB ?_ ?_
  · push_cast; ring
  · simp only [map_add, map_mul, Complex.conj_ofReal, Complex.conj_I]; push_cast; ring

/-! ## G = u^{-p} -/

lemma hasDerivAt_Gfun (p : ℝ) {u : ℝ} (hu : 0 < u) :
    HasDerivAt (Gfun p) (((-p) * u ^ (-p - 1) : ℝ) : ℂ) u :=
  (Real.hasDerivAt_rpow_const (Or.inl hu.ne')).ofReal_comp

/-- (43): G → 0. -/
lemma eq43 {p : ℝ} (hp : 0 < p) : Tendsto (Gfun p) atTop (𝓝 0) := by
  have := (Complex.continuous_ofReal.tendsto 0).comp (tendsto_rpow_neg_atTop hp)
  rwa [Complex.ofReal_zero] at this

/-! ## Asymptotic steps (true; proofs pending) -/

/-- (39): F/G = u^p F(u) → 0 along D_γ. From (34),(35) and the identity
u^{s+p}J(u,s) - u^{1-s̄+p}J(u,1-s̄) = u^p F(u): u^p F = (1-2σ₀)B̄₂ u^{p-1}/2 + O(u^{p-2}). -/
lemma eq39 (hin : Inputs) (X : Setting) (M₀ : ℕ) :
    Tendsto (fun u => Ffun X.η X.σ₀ X.t₀ u / Gfun X.p u) (along (Dγ X.γ M₀)) (𝓝 0) := by
  sorry -- HARDEN

/-- (42): F → 0 along D_γ. -/
lemma eq42 (hin : Inputs) (X : Setting) (M₀ : ℕ) :
    Tendsto (Ffun X.η X.σ₀ X.t₀) (along (Dγ X.γ M₀)) (𝓝 0) := by
  sorry -- HARDEN

/-- (50): σ₀ ≠ 1/2 ⇒ |F'/G'| → ∞ along D_γ; indeed F'/G' = (1-2σ₀)u^p B̄₁(u)/(-p) + o(1)
and B̄₁ ≥ γ on D_γ. From (44)-(49). -/
lemma eq50 (hin : Inputs) (X : Setting) (M₀ : ℕ) (hσ : X.σ₀ ≠ 1 / 2) :
    Tendsto (fun u => deriv (Ffun X.η X.σ₀ X.t₀) u / deriv (Gfun X.p) u) (along (Dγ X.γ M₀))
      (Bornology.cobounded ℂ) := by
  sorry -- HARDEN

/-- (51)–(54): oscillation of F/G on each D_M is ≤ δ_M → 0. -/
lemma eq54 (hin : Inputs) (X : Setting) (M₀ : ℕ) :
    ∃ δ : ℕ → ℝ, Tendsto δ atTop (𝓝 0) ∧ ∀ M, ∀ u ∈ Ioo (aγ X.γ M₀ M) (bγ X.γ M₀ M),
      ‖Ffun X.η X.σ₀ X.t₀ u / Gfun X.p u -
        Ffun X.η X.σ₀ X.t₀ (bγ X.γ M₀ M) / Gfun X.p (bγ X.γ M₀ M)‖ ≤ δ M := by
  sorry -- HARDEN

/-! ## The paper's conclusion, (55)–(56) -/

/-- The paper's argument, checked: given the literature inputs AND its Theorem 1, RH follows. -/
theorem fayed_RH (hin : Inputs) (thm1 : Theorem1C) : RiemannHypothesis := by
  intro s hs hnt hs1
  obtain ⟨hσ0, hσ1⟩ := hin.strip s hs hnt hs1
  by_contra hne
  have hm : 0 < min s.re (1 - s.re) := lt_min hσ0 (by linarith)
  have hmle1 : min s.re (1 - s.re) ≤ s.re := min_le_left _ _
  have hmle2 : min s.re (1 - s.re) ≤ 1 - s.re := min_le_right _ _
  let X : Setting :=
    { σ₀ := s.re, t₀ := s.im, p := 1 / 2, η := min s.re (1 - s.re) / 2, γ := 1 / 8,
      hσ₀ := ⟨hσ0, hσ1⟩, hzero := by rw [Complex.re_add_im]; exact hs,
      hp := by norm_num, hη := ⟨by linarith, by linarith, by linarith⟩, hγ := by norm_num }
  obtain ⟨-, -, ha0, hab, hba, hb, -⟩ := counterexample (γ := 1 / 8) (by norm_num) (by norm_num) 0
  obtain ⟨δ, hδ, hosc⟩ := eq54 hin X 0
  have hneBot : (along (Dγ X.γ 0)).NeBot := along_Dγ_neBot (by norm_num) 0
  -- (55): Theorem 1 turns (42),(43),(50),(54) into F/G → ∞
  have h55 := thm1 (Ffun X.η X.σ₀ X.t₀) (Gfun X.p) (aγ X.γ 0) (bγ X.γ 0) δ
    (Bornology.cobounded ℂ) (Or.inl rfl)
    (fun u hu => (hin.leibniz44 _ _ _ X.hη.1 X.hη.2.1 X.hη.2.2 u hu).differentiableAt)
    (fun u hu => (hasDerivAt_Gfun X.p hu).differentiableAt)
    ha0 hab hba hb (eq42 hin X 0) ((eq43 X.hp.1).mono_left (along_le_atTop _))
    (eq50 hin X 0 hne) hosc hδ
  -- (56): contradiction with (39)
  rw [← tendsto_norm_atTop_iff_cobounded] at h55
  have h39 := (eq39 hin X 0).norm
  rw [norm_zero] at h39
  exact not_tendsto_nhds_of_tendsto_atTop h55 0 h39

/-- **Where the paper breaks.** Its argument is valid given its Theorem 1 (fayed_RH),
and Theorem 1 is false (theorem1C_false). The failing step is (55). -/
theorem fayed_gap : (Inputs → Theorem1C → RiemannHypothesis) ∧ ¬ Theorem1C :=
  ⟨fayed_RH, theorem1C_false⟩

end

end Fayed

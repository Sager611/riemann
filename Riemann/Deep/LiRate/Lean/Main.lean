import Riemann.Deep.LiRate.Lean.Trend
import Riemann.Deep.LiRate.Lean.Cor3
import Riemann.Deep.LiRate.Lean.Lemma0
import Riemann.Deep.LiRate.Lean.Fejer

/-!
# LiRate: the Li-remainder theorems from named inputs (summary module)

`R(n) = λ_n - t(n)`. `P` is any proposition standing for RH (instantiate with Mathlib's
`RiemannHypothesis`); the analytic inputs are hypotheses, see `LEAN.md`.
-/

open Finset Real Filter

namespace LiRate

/-- The inputs of Theorem 1: (H_Li) Li's formula with real phases, `t(0) = 0`, (H_N) RvM,
(H_w)+(H_τ) weights with the trend/tail estimate. -/
structure Thm1Inputs {κ : Type*} (θ : κ → ℝ) (lam t : ℕ → ℝ) (Z : ℝ → Finset κ) : Prop where
  li : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)
  t0 : t 0 = 0
  rvm : RvM Z
  weights : ∀ η > 0, ∃ w : ℕ → ℕ → ℝ, GoodWeights w η ∧ TrendTail w θ lam t Z η

/-- The finer inputs of Theorem 1 (the lemmas of proofs.md): (H_Li), `t(0) = 0`, (H_N), (H7), and
for every `η > 0` weights with (H_w), (H_tiny), (H4). -/
structure Thm1Lemmas {κ : Type*} [DecidableEq κ] (θ : κ → ℝ) (lam t : ℕ → ℝ) (Z : ℝ → Finset κ) :
    Prop where
  li : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)
  t0 : t 0 = 0
  rvm : RvM Z
  tauT : TauT θ Z
  weights : ∀ η > 0, ∃ w : ℕ → ℕ → ℝ, ∃ Xs : ℕ → Finset κ,
    GoodWeights w η ∧ TinyTerms w θ lam Z Xs η ∧ TrendPt w θ t Z

variable {κ : Type*} {θ : κ → ℝ} {lam t : ℕ → ℝ} {Z : ℝ → Finset κ}

/-- The finer inputs give the inputs of `thm1`. -/
theorem Thm1Lemmas.toInputs [DecidableEq κ] (h : Thm1Lemmas θ lam t Z) : Thm1Inputs θ lam t Z :=
  ⟨h.li, h.t0, h.rvm, fun η hη => by
    obtain ⟨w, Xs, hg, ⟨hsub, htiny⟩, h4⟩ := h.weights η hη
    exact ⟨w, hg, trendTail_of w θ lam t Z η hg.2.2.1 Xs hsub htiny
      (htrend_of w θ t Z h.rvm h4 h.tauT η hη)⟩⟩

/-- Theorem 1. -/
theorem li_thm1 (h : Thm1Inputs θ lam t Z) : Thm1Bound (fun n => lam n - t n) :=
  thm1 θ lam t h.li h.t0 Z h.rvm h.weights

/-- Corollary 1: `limsup |R(n)|/(n log n)^{1/4} ≥ 0.289`. -/
theorem li_cor1 (h : Thm1Inputs θ lam t Z) :
    ∃ᶠ n : ℕ in atTop, (0.289 : ℝ) * √(√(n * Real.log n)) ≤ |lam n - t n| :=
  cor1_num _ (li_thm1 h)

/-- Corollary 2: `Σ m^{1/2} y_m² = ∞`. -/
theorem li_cor2 (h : Thm1Inputs θ lam t Z) :
    ¬ Summable (fun m : ℕ => √(m : ℝ) * ((lam m - t m) / m) ^ 2) :=
  cor2 _ (li_thm1 h)

/-- Corollary 3: with Theorem 2's upper bound, `Σ m^δ y_m² < ∞ ⟺ δ < 1/2`. -/
theorem li_cor3 (h : Thm1Inputs θ lam t Z) (h2 : Thm2Upper (fun n => lam n - t n)) (δ : ℝ) :
    Summable (fun m : ℕ => (m : ℝ) ^ δ * ((lam m - t m) / m) ^ 2) ↔ δ < 1 / 2 :=
  cor3 _ (li_thm1 h) h2 δ

/-- Corollary 1, unconditional: `R(n) = Ω((n log n)^{1/4})`. -/
theorem li_omega (P : Prop) (h1 : P → Thm1Inputs θ lam t Z)
    (h0 : ¬P → Lemma0 (fun n => lam n - t n)) :
    ∃ c > 0, ∃ᶠ n : ℕ in atTop, c * √(√(n * Real.log n)) ≤ |lam n - t n| :=
  omega_uncond _ P (fun hP => li_thm1 (h1 hP)) h0

/-- Corollary 4: Merlini–Sala–Sala's `|λ̃_n| ≤ a log n` fails for every `a`, and `λ̃_n = O(n^ε)`
fails for every `ε < 1/4`, where `λ̃_n - R(n)` is bounded. -/
theorem li_cor4 (P : Prop) (h1 : P → Thm1Inputs θ lam t Z)
    (h0 : ¬P → Lemma0 (fun n => lam n - t n)) (lt : ℕ → ℝ)
    (hlt : ∃ C, ∀ n, |lt n - (lam n - t n)| ≤ C) :
    (∀ a : ℝ, ¬ ∀ᶠ n : ℕ in atTop, |lt n| ≤ a * Real.log n) ∧
      ∀ ε < (1 / 4 : ℝ), ∀ C' : ℝ, ¬ ∀ᶠ n : ℕ in atTop, |lt n| ≤ C' * (n : ℝ) ^ ε :=
  ⟨fun a => cor4 _ lt P (fun hP => li_thm1 (h1 hP)) h0 hlt a,
    fun _ hε C' => cor4_rpow _ lt P (fun hP => li_thm1 (h1 hP)) h0 hlt hε C'⟩

/-- Corollary 3, criterion: `RH ⟺ Σ m^δ y_m² < ∞` for some / every `δ ∈ [0, 1/2)`. -/
theorem li_rh_iff (P : Prop) (h2 : P → Thm2Upper (fun n => lam n - t n))
    (h0 : ¬P → Lemma0 (fun n => lam n - t n)) :
    (P ↔ ∃ δ ∈ Set.Ico (0 : ℝ) (1 / 2),
        Summable (fun m : ℕ => (m : ℝ) ^ δ * ((lam m - t m) / m) ^ 2)) ∧
      (P ↔ ∀ δ ∈ Set.Ico (0 : ℝ) (1 / 2),
        Summable (fun m : ℕ => (m : ℝ) ^ δ * ((lam m - t m) / m) ^ 2)) :=
  ⟨rh_iff_exists _ P h2 h0, rh_iff_forall _ P h2 h0⟩

/-- Corollary 1, unconditional, for the actual trend `t = liTrend`, with Lemma 0 from the
radius input (¬RH ⇒ `λ_n r^n` unbounded for some `r < 1`). -/
theorem li_omega_radius (P : Prop) (h1 : P → Thm1Inputs θ lam liTrend Z)
    (hrad : ¬P → ∃ r : ℝ, 0 < r ∧ r < 1 ∧ ¬ ∃ B : ℝ, ∀ n, |lam n| * r ^ n ≤ B) :
    ∃ c > 0, ∃ᶠ n : ℕ in atTop, c * √(√(n * Real.log n)) ≤ |lam n - liTrend n| :=
  li_omega P h1 fun hP => by
    obtain ⟨r, hr0, hr1, hunb⟩ := hrad hP
    exact lemma0_of_unbounded lam liTrend hr0 hr1 hunb 4 liTrend_bound

/-- Corollary 3, criterion, for the actual trend, with Lemma 0 from the radius input. -/
theorem li_rh_iff_radius (P : Prop) (h2 : P → Thm2Upper (fun n => lam n - liTrend n))
    (hrad : ¬P → ∃ r : ℝ, 0 < r ∧ r < 1 ∧ ¬ ∃ B : ℝ, ∀ n, |lam n| * r ^ n ≤ B) :
    (P ↔ ∃ δ ∈ Set.Ico (0 : ℝ) (1 / 2),
        Summable (fun m : ℕ => (m : ℝ) ^ δ * ((lam m - liTrend m) / m) ^ 2)) ∧
      (P ↔ ∀ δ ∈ Set.Ico (0 : ℝ) (1 / 2),
        Summable (fun m : ℕ => (m : ℝ) ^ δ * ((lam m - liTrend m) / m) ^ 2)) :=
  li_rh_iff P h2 fun hP => by
    obtain ⟨r, hr0, hr1, hunb⟩ := hrad hP
    exact lemma0_of_unbounded lam liTrend hr0 hr1 hunb 4 liTrend_bound

end LiRate

#print axioms LiRate.core_bound
#print axioms LiRate.tail_eq_tsum_compl
#print axioms LiRate.trendTail_of
#print axioms LiRate.thm1
#print axioms LiRate.htrend_of
#print axioms LiRate.thm1_of_lemmas
#print axioms LiRate.Thm1Lemmas.toInputs
#print axioms LiRate.li_thm1
#print axioms LiRate.li_cor1
#print axioms LiRate.li_cor2
#print axioms LiRate.li_cor3
#print axioms LiRate.li_omega
#print axioms LiRate.li_cor4
#print axioms LiRate.li_rh_iff
#print axioms LiRate.lemma0_of_unbounded
#print axioms LiRate.li_omega_radius
#print axioms LiRate.li_rh_iff_radius
#print axioms LiRate.fejer_goodWeights

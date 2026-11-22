module

public import Riemann.Deep.LeanLi.Lit

/-!
# LeanLi: the Claim of `discoveries/li-remainder.md`, as Lean propositions

`R n = λ_n - (n/2)(log n - 1 - log 2π + γ)` with `λ_n = lim_{T→∞} Σ_{|ρ|≤T} m_ρ [1 - (1-1/ρ)ⁿ]`
(`Defs.lean`). "Under RH" is Mathlib's `RiemannHypothesis`; `≪`, `O`, `Ω`, `o(1)` are explicit.

Scope. Not formalized: Claim 5's last sentence "Numerically a = 2 first fails at n = 5083" (a numerical
computation of `λ_n`, not a theorem; out of scope). The parenthetical "(Li: RH ⟺ λ_n > 0 for all n)" is
cited background, not part of the Claim. Zeros are counted with multiplicity throughout (no simplicity).
Sums `Σ_{n≤N}`, `Σ_{N<n≤2N}` run over integers `N` (`∀ᶠ N : ℕ`); Claim 4's scale `N` is real.
-/

@[expose] public section

open Complex Filter Topology
open scoped ContDiff

namespace LeanLi

/-- The definition of `λ_n` is sound: the partial sums converge, and the limit is real. -/
def Part0 : Prop :=
  ∀ n : ℕ, Tendsto (liPartial n) atTop (𝓝 (liCoeff n)) ∧ (liCoeff n).im = 0

/-- Claim 1. Under RH: `Σ_{n≤N} R(n)² ≥ (c₁ - o(1)) N^{3/2} (log N)^{1/2}` with
`c₁ = (3/5)^{1/2}(4/5)²/(2π√2) = 0.0558`; `Σ_{N<n≤2N} R(n)² ≪ N^{3/2}(log N)^{3/2}`; and (with
Montgomery's pair-correlation theorem, a literature input) `≪ N^{3/2}(log N)^{1/2} log log N`. -/
def Part1 : Prop :=
  |c₁ - 0.0558| < 0.00005 ∧
  (RiemannHypothesis →
    (∀ ε > (0 : ℝ), ∀ᶠ N : ℕ in atTop,
      (c₁ - ε) * ((N : ℝ) ^ (3 / 2 : ℝ) * Real.log N ^ (1 / 2 : ℝ)) ≤
        ∑ n ∈ Finset.Icc 1 N, R n ^ 2) ∧
    (∃ C : ℝ, ∀ᶠ N : ℕ in atTop,
      ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤
        C * ((N : ℝ) ^ (3 / 2 : ℝ) * Real.log N ^ (3 / 2 : ℝ))) ∧
    (∃ C : ℝ, ∀ᶠ N : ℕ in atTop,
      ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤
        C * ((N : ℝ) ^ (3 / 2 : ℝ) * Real.log N ^ (1 / 2 : ℝ) * Real.log (Real.log N))))

/-- Claim 2. Unconditionally `R(n) = Ω((n log n)^{1/4})`; under RH
`limsup |R(n)|/(n log n)^{1/4} ≥ 0.289` (even: `≥ 0.289` infinitely often); if RH fails,
`|R(n)| ≥ e^{cn}` infinitely often. -/
def Part2 : Prop :=
  (∃ c > (0 : ℝ), ∃ᶠ n : ℕ in atTop, c * ((n : ℝ) * Real.log n) ^ (1 / 4 : ℝ) ≤ |R n|) ∧
  (RiemannHypothesis →
    (∃ᶠ n : ℕ in atTop, (0.289 : ℝ) * ((n : ℝ) * Real.log n) ^ (1 / 4 : ℝ) ≤ |R n|) ∧
    ((0.289 : ℝ) : EReal) ≤
      limsup (fun n : ℕ => ((|R n| / ((n : ℝ) * Real.log n) ^ (1 / 4 : ℝ) : ℝ) : EReal)) atTop) ∧
  (¬ RiemannHypothesis → ∃ c > (0 : ℝ), ∃ᶠ n : ℕ in atTop, Real.exp (c * n) ≤ |R n|)

/-- Claim 3. `RH ⟺ Σ_m m^δ y_m² < ∞` for some (equivalently every) `δ ∈ [0, 1/2)`, Keiper's `y_m`;
under RH the weight `m^{1/2}` is never summable. -/
def Part3 : Prop :=
  (RiemannHypothesis ↔
    ∃ δ ∈ Set.Ico (0 : ℝ) (1 / 2), Summable (fun m : ℕ => (m : ℝ) ^ δ * keiperY m ^ 2)) ∧
  (RiemannHypothesis ↔
    ∀ δ ∈ Set.Ico (0 : ℝ) (1 / 2), Summable (fun m : ℕ => (m : ℝ) ^ δ * keiperY m ^ 2)) ∧
  (RiemannHypothesis → ¬ Summable (fun m : ℕ => (m : ℝ) ^ (1 / 2 : ℝ) * keiperY m ^ 2))

/-- Claim 4. Under RH, for `φ ∈ C_c^∞(0, ∞)`:
`Σ_n φ(n/N) R(n) = (7/4) Σ_n φ(n/N) + O(log N)`, for every real scale `N ≥ 2` (the sums are finite;
`O(log N)` is uniform in `N`, its constant depends on `φ`). -/
def Part4 : Prop :=
  RiemannHypothesis → ∀ φ : ℝ → ℝ, ContDiff ℝ ∞ φ → HasCompactSupport φ →
    tsupport φ ⊆ Set.Ioi 0 →
      ∃ C : ℝ, ∀ N : ℝ, 2 ≤ N →
        |∑' n : ℕ, φ (n / N) * R n - (7 / 4) * ∑' n : ℕ, φ (n / N)| ≤ C * Real.log N

/-- Claim 5. Merlini–Sala–Sala's conjecture `|λ_tiny(n)| ≤ a log n` is false for every `a`, and
`λ_tiny = O(n^ε)` fails for every `ε < 1/4`: for every trend within `O(1)` of
`(n/2)(log n - 1 - log 2π + γ)` (their eq. (2): `(n/2) log n - 1.13…·n`), in particular for
Maślanka's trend `S̄_n` (Voros's `vorosTrend`). -/
def Part5 : Prop :=
  (∀ trend : ℕ → ℝ, (∃ C : ℝ, ∀ n, |trend n - liTrend n| ≤ C) →
    (∀ a : ℝ, ¬ ∀ᶠ n : ℕ in atTop, |lam n - trend n| ≤ a * Real.log n) ∧
    (∀ ε : ℝ, ε < 1 / 4 →
      ¬ (fun n : ℕ => lam n - trend n) =O[atTop] (fun n : ℕ => (n : ℝ) ^ ε))) ∧
  (∀ a : ℝ, ¬ ∀ᶠ n : ℕ in atTop, |lam n - vorosTrend n| ≤ a * Real.log n) ∧
  (∀ ε : ℝ, ε < 1 / 4 →
    ¬ (fun n : ℕ => lam n - vorosTrend n) =O[atTop] (fun n : ℕ => (n : ℝ) ^ ε))

/-- The Claim of `discoveries/li-remainder.md` (parts 1–5), with the soundness of `λ_n`. -/
def LiRemainderLaws : Prop := Part0 ∧ Part1 ∧ Part2 ∧ Part3 ∧ Part4 ∧ Part5

end LeanLi

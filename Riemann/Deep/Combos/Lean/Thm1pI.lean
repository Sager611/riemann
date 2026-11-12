module

public import Riemann.Deep.Combos.Lean.Count
public import Riemann.Deep.Combos.Lean.Expo1pB
public import Riemann.Deep.Combos.Lean.Expo1pC
public import Riemann.Deep.Combos.Lean.Expo1pD

@[expose] public section

/-!
# Theorem 1′(a) (RCApps2 NOTE Prop 1 + Thm 1′(a), l.28–68): interface and named inputs

Ring class family of one D < 0: X = Ĝ (any finite type), zero data (Count.lean), |D|, the analytic
conductors A(χ*) ≤ |D|, and the two detectors of Lemma 5 as real numbers:
detI V Y χ ρ = |Σ_{N𝔞>V} a_𝔞 χ(𝔞) N𝔞^{−ρ} e^{−N𝔞/Y}|, detII V χ γ = ∫|L(½+i(γ+u),χ)M_V(½+i(γ+u),χ)|e^{−|u|}du.
-/

open Real Finset

namespace RingClass

/-- analytic data of the ring class family of one D < 0 -/
structure RCZeros (X : Type*) extends ZeroData X where
  /-- |D| -/
  absD : ℝ
  absD_ge : 3 ≤ absD
  /-- A(χ*) = |d_K| N𝔣_χ (NOTE l.35) -/
  A : X → ℝ
  A_ge : ∀ χ, 1 ≤ A χ
  A_le : ∀ χ, A χ ≤ absD
  /-- Lemma 5 (I) quantity, mollifier length V, smoothing Y (NOTE l.139) -/
  detI : ℝ → ℝ → X → ℂ → ℝ
  /-- Lemma 5 (II) quantity, mollifier length V, at height γ (NOTE l.140) -/
  detII : ℝ → X → ℝ → ℝ

/-- the named inputs of Theorem 1′(a), with constants c₀, c₁, C_I(ε₀), C_II(ε₀), 𝔔₀(ε₀) that do not
depend on D. 𝔔 = 9|D| (= |D|(M₀ + 3)², M₀ = 0), V = 𝔔^{1/2}, T = |t| + 3. -/
structure RCZeros.Hyp {X : Type*} [Fintype X] (d : RCZeros X) (c₀ c₁ : ℝ)
    (CI CII Q0 : ℝ → ℝ) : Prop where
  /-- (F3) NOTE l.52: |Ĝ| = h(D) ≤ 3|D|^{1/2} log²(3|D|) -/
  F3 : (Fintype.card X : ℝ) ≤ 3 * d.absD ^ (1 / 2 : ℝ) * Real.log (3 * d.absD) ^ 2
  /-- (F4) NOTE l.54–55 (Iwaniec–Kowalski Prop 5.7) -/
  F4 : ∀ χ t, (∑ ρ ∈ d.win χ t, (d.mult χ ρ : ℝ)) ≤ c₀ * Real.log (d.A χ * (|t| + 3) ^ 2)
  c₀_nonneg : 0 ≤ c₀
  /-- Lemma 5 (NOTE l.134–149): the class I / class II split, any mollifier length 2 ≤ V ≤ Y -/
  L5 : ∀ σ, 0.6 ≤ σ → σ ≤ 1 → ∀ V Y, 2 ≤ V → V ≤ Y → ∀ χ, χ ≠ d.one → ∀ t, ∀ ρ ∈ d.win χ t,
    σ ≤ ρ.re → 1 / 4 ≤ d.detI V Y χ ρ ∨ c₁ * Y ^ (σ - 1 / 2) ≤ d.detII V χ ρ.im
  /-- class II, first moment (RCApps2 l.40–48: Cauchy–Schwarz, [D2] (2.1) = Lemma 4 at mollifier
  length 1, Lemma 2, Lemma 3), at Y = (𝔔T)^{1/(3−2σ)}: |𝒳_II| ≤ C(ε₀) Y^{1/2−σ} 𝔔^{1/2+ε₀} T^{1/2+ε₀} -/
  classII : ∀ ε₀, 0 < ε₀ → ε₀ ≤ 1 / 10 → ∀ σ, 0.6 ≤ σ → σ ≤ 1 → ∀ t,
    ∀ (S : Finset X) (ρ : X → ℂ), (∀ χ ∈ S, χ ≠ d.one ∧ ρ χ ∈ d.win χ t ∧ σ ≤ (ρ χ).re ∧
      c₁ * ((9 * d.absD * (|t| + 3)) ^ (1 / (3 - 2 * σ))) ^ (σ - 1 / 2) ≤
        d.detII ((9 * d.absD) ^ (1 / 2 : ℝ)) χ (ρ χ).im) →
    (S.card : ℝ) ≤ CII ε₀ * ((9 * d.absD * (|t| + 3)) ^ (1 / (3 - 2 * σ))) ^ (1 / 2 - σ) *
      (9 * d.absD) ^ (1 / 2 + ε₀) * (|t| + 3) ^ (1 / 2 + ε₀)
  /-- class I, no k-th powers (RCApps2 l.50–55: tail beyond Z, dyadic blocks, removal of β, (LS) =
  Lemma 2, Sobolev, Lemma 3), for 𝔔 ≥ 𝔔₀(ε₀), at Y = (𝔔T)^{1/(3−2σ)}:
  |𝒳_I| ≤ C(ε₀) (𝔔T)^{ε₀} Z^{2−2σ+2ε₀}, Z = 2Y log²Y -/
  classI : ∀ ε₀, 0 < ε₀ → ε₀ ≤ 1 / 10 → Q0 ε₀ ≤ 9 * d.absD → ∀ σ, 0.6 ≤ σ → σ ≤ 1 → ∀ t,
    ∀ (S : Finset X) (ρ : X → ℂ), (∀ χ ∈ S, χ ≠ d.one ∧ ρ χ ∈ d.win χ t ∧ σ ≤ (ρ χ).re ∧
      1 / 4 ≤ d.detI ((9 * d.absD) ^ (1 / 2 : ℝ)) ((9 * d.absD * (|t| + 3)) ^ (1 / (3 - 2 * σ)))
        χ (ρ χ)) →
    (S.card : ℝ) ≤ CI ε₀ * (9 * d.absD * (|t| + 3)) ^ ε₀ *
      (2 * (9 * d.absD * (|t| + 3)) ^ (1 / (3 - 2 * σ)) *
        Real.log ((9 * d.absD * (|t| + 3)) ^ (1 / (3 - 2 * σ))) ^ 2) ^ (2 - 2 * σ + 2 * ε₀)
  CI_nonneg : ∀ ε₀, 0 ≤ CI ε₀
  CII_nonneg : ∀ ε₀, 0 ≤ CII ε₀

end RingClass

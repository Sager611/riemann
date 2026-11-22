module

public import Riemann.Deep.LeanLi.Defs

/-!
# LeanLi: literature inputs (theorems of others, as named `Prop`s)

Each input is a published theorem (or a classical textbook fact missing from Mathlib), stated no
stronger than its source: asymptotic statements are written in their eventual form
(`∃ C T₀, ∀ T ≥ T₀`). They are bundled in `LiteratureInputs`, an explicit hypothesis of the main
theorem, so `#print axioms` stays standard.
-/

@[expose] public section

open Complex Filter Topology

namespace LeanLi

/-- **ζ has no zeros on the real segment `(0, 1)`** (classical: `ζ(σ) < 0` for `0 < σ < 1`;
T. Nakamura, arXiv 1405.1504v3, TeX l.308–312, for the Hurwitz zeta function with `a = 1`:
"we obtain ζ(σ,a) < 0 for all 0 < σ < 1 ... Therefore ζ(σ,a) does not vanish in the interval
(0,1) when a ≥ 1/2", via Apostol, Introduction to Analytic Number Theory, Thm 12.21). -/
def ZetaNoZerosUnitInterval : Prop :=
  ∀ x : ℝ, 0 < x → x < 1 → riemannZeta x ≠ 0

/-- **Stirling asymptotics of the Riemann–Siegel theta function** (arXiv 2004.00926, TeX l.70–74:
"The well-known asymptotic expansion for ϑ(t) is given by
ϑ(t) ∼ ½ t (log(t/2π) − 1) − π/8 + Σ_{k=1}^∞ (1 − 2^{1−2k})|B_{2k}| / (4k(2k−1) t^{2k−1})
+ ½ arctan(e^{−πt}), valid as t → +∞"; Edwards, Riemann's Zeta Function, §6.5).
Its first terms: `θ(t) = (t/2)(log(t/2π) - 1) - π/8 + O(1/t)`. -/
def RSThetaAsymp : Prop :=
  ∃ C t₀ : ℝ, ∀ t : ℝ, t₀ ≤ t →
    |rsTheta t - (t / 2 * (Real.log (t / (2 * Real.pi)) - 1) - Real.pi / 8)| ≤ C / t

/-- **Asymptotics of the digamma function** (DLMF 5.11.2: "As z → ∞ in the sector
|ph z| ≤ π − δ, ψ(z) ∼ ln z − 1/(2z) − Σ_{k=1}^∞ B_{2k}/(2k z^{2k})"). Here the sector
`|ph z| ≤ π/2` and the remainder after two terms: `ψ(z) = ln z - 1/(2z) + O(|z|^{-2})`. -/
def DigammaAsymp : Prop :=
  ∃ C R₀ : ℝ, ∀ z : ℂ, R₀ ≤ ‖z‖ → 0 ≤ z.re →
    ‖digamma z - (Complex.log z - 1 / (2 * z))‖ ≤ C / ‖z‖ ^ 2

/-- **Riemann–von Mangoldt: `S(T) = O(log T)`** (Titchmarsh, Thm 9.4, with the exact formula
`N(T) = θ(T)/π + 1 + S(T)` of §9.3; quoted in Brent–Platt–Trudgian, arXiv 2009.13791,
TeX l.263–292 and l.313–325: "From [Titch, Thm. 9.4], S(T) ≪ log T"; "The Riemann–von Mangoldt
formula states, in its most precise form, N(t) = θ(t)/π + 1 + S(t)"). -/
def RvMSBound : Prop :=
  ∃ C T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → |zetaS T| ≤ C * Real.log T

/-- **`S₁(T) = ∫_0^T S(t) dt = O(log T)`** (Titchmarsh, Thm 9.9(A); quoted in Brent–Platt–Trudgian,
arXiv 2009.13791, TeX l.373–379: "Define S₁(T) := ∫_0^T S(t) dt. We know that S₁(T) ≪ log T"). -/
def S1Bound : Prop :=
  ∃ C T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → |∫ t in (0 : ℝ)..T, zetaS t| ≤ C * Real.log T

/-- **Selberg's mean value of `S(t)²`** (A. Selberg, Contributions to the theory of the Riemann
zeta-function, Arch. Math. Naturvid. 48 (1946), no. 5, 89–155, Theorems 6 and 7, unconditional;
quoted in arXiv 2006.08503, TeX l.227–230: "∫_0^T |S(t)|² dt = (T/(2π²)) log log T +
O(T √(log log T))"). Used here only as an upper bound. -/
def SelbergS2 : Prop :=
  ∃ C T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
    ∫ t in (0 : ℝ)..T, zetaS t ^ 2 ≤ C * T * Real.log (Real.log T)

/-- **Montgomery's pair-correlation theorem** (H. L. Montgomery, The pair correlation of zeros of the
zeta function, Analytic Number Theory, Proc. Sympos. Pure Math. 24 (1973), 181–193; restated in
D. A. Goldston, arXiv math/0412313, TeX l.257–264: "Assume the Riemann Hypothesis. Then F(α) is
real, even, and non-negative. Further, uniformly for 0 ≤ α ≤ 1 − ε, we have
F(α) = α + o(1) + (1 + o(1)) T^{−2α} log T"). Used for `0 ≤ α ≤ 1/2` (`ε = 1/2`). -/
def MontgomeryPC : Prop :=
  RiemannHypothesis → ∀ ε > (0 : ℝ), ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → ∀ α : ℝ, 0 ≤ α → α ≤ 1 / 2 →
    ‖montgomeryF α T - ((α : ℂ) + ((T ^ (-2 * α) * Real.log T : ℝ) : ℂ))‖ ≤
      ε * (1 + T ^ (-2 * α) * Real.log T)

/-- **Voros's asymptotic expansion (TS) of the trend `S̄_n`** (A. Voros, Sharpenings of Li's
criterion for the Riemann Hypothesis, Math. Phys. Anal. Geom. 9 (2006) 53–63, arXiv math/0506326,
TeX l.494–499, eq. (TS): "S̄_n ∼ ½ n (log n − 1 + γ − log 2π) + 3/4 − Σ_{k=1}^∞ B_{2k}/(4k)
n^{1−2k} (n → ∞) unconditionally"). Its first terms: `S̄_n = t(n) + 3/4 + O(1/n)`. -/
def VorosTS : Prop :=
  ∃ C : ℝ, ∃ n₀ : ℕ, ∀ n : ℕ, n₀ ≤ n → |vorosTrend n - liTrend n - 3 / 4| ≤ C / n

/-- The literature inputs except Montgomery's pair-correlation theorem (all unconditional theorems).
Everything except Claim 1's pair-correlation bound follows from these. -/
structure LitCore : Prop where
  zeta_unit_interval : ZetaNoZerosUnitInterval
  rs_theta : RSThetaAsymp
  digamma : DigammaAsymp
  rvm : RvMSBound
  s1 : S1Bound
  selberg : SelbergS2
  voros : VorosTS

/-- The literature inputs, bundled: `LitCore` and Montgomery's pair-correlation theorem. -/
structure LiteratureInputs : Prop extends LitCore where
  montgomery : MontgomeryPC

instance : Coe LiteratureInputs LitCore := ⟨LiteratureInputs.toLitCore⟩

end LeanLi

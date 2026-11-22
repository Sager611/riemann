import Riemann.Deep.LeanNBRate.Defs
import Mathlib.MeasureTheory.Integral.Bochner.Set

/-!
# Literature inputs (theorems of others, stated no stronger than their sources)

Each input is a `def … : Prop`; the main theorem takes them as explicit hypotheses bundled in
`LitInputs`, so the main theorem uses only propext, Classical.choice and Quot.sound. Every docstring gives the citation and the source
statement verbatim (TeX copied from the arXiv e-print where one exists).
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-! ## Objects used only to state the inputs -/

/-- Complex `L²(0,∞)` (Burnol's `K = L²(]0,∞[, dt)` over the complex numbers). -/
abbrev L2posC : Type := Lp ℂ 2 volPos

theorem volPos_Ioc_ne_top : volPos (Ioc (0 : ℝ) 1) ≠ ⊤ := by
  refine ne_of_lt (lt_of_le_of_lt (Measure.restrict_apply_le (Ioi (0 : ℝ)) (Ioc (0 : ℝ) 1)) ?_)
  rw [Real.volume_Ioc]
  exact ENNReal.ofReal_lt_top

/-- Burnol's `χ`, the indicator function of `]0,1]`, in complex `L²(0,∞)`. -/
def chiVecC : L2posC := indicatorConstLp 2 measurableSet_Ioc volPos_Ioc_ne_top (1 : ℂ)

/-- `t ↦ ρ(θ/t)` (`ρ` = fractional part), complex-valued. -/
def nbFunC (θ : ℝ) (x : ℝ) : ℂ := ((Int.fract (θ / x) : ℝ) : ℂ)

/-- The classes of `t ↦ ρ(θ/t)`, `λ ≤ θ ≤ 1`, in complex `L²(0,∞)`. -/
def nbSetC (lam : ℝ) : Set L2posC := {f | ∃ θ ∈ Icc lam 1, (f : ℝ → ℂ) =ᵐ[volPos] nbFunC θ}

/-- Burnol's `D(λ) = inf_{f ∈ B_λ} ‖χ − f‖`, `B_λ` the finite complex linear combinations. -/
def nbDistC (lam : ℝ) : ℝ := Metric.infDist chiVecC (Submodule.span ℂ (nbSetC lam) : Set L2posC)

/-- Nontrivial zeros with `a < Im ρ ≤ b`. -/
def zerosIn (a b : ℝ) : Set ℂ := {ρ | ρ ∈ NTZ ∧ a < ρ.im ∧ ρ.im ≤ b}

/-- `N(b) − N(a)`: nontrivial zeros with ordinate in `(a, b]`, counted with multiplicity. -/
def Ncount (a b : ℝ) : ℕ := ∑ᶠ ρ ∈ zerosIn a b, mult ρ

/-! ## The inputs -/

/-- **Burnol 2002.** J.-F. Burnol, *A lower bound in an approximation problem involving the zeros of
the Riemann zeta function*, Adv. Math. 170 (2002) 56–70, arXiv:math/0103058 (TeX l.169–206).
Setting (l.169–175): "Let $K = L^2(]0,\infty[,dt)$ (over the complex numbers), let $\chi$ be the
indicator function of the interval $]0,1]$, and let $\rho$ be the function ``fractional part'' […]
Let $0<\lambda<1$ and let $\cB_\lambda$ be the sub-vector space of $K$ consisting of the finite
linear combinations of the functions $t\mapsto \rho({\theta\over t})$, for $\lambda\leq\theta\leq1$."
(l.189) "$D(\la)$ for the Hilbert-space distance $\inf_{f\in\cB_\lambda} \| \chi -f \|$."
(l.197–206) "If the Riemann Hypothesis fails this result is true but trivial […] So we will assume
that the Riemann Hypothesis holds. […] counted \emph{only once} independently of their
multiplicities $m_\rho$. We prove the following: \begin{thm} We have:
$$\liminf_{\la\to0} D(\la)\sqrt{\log({1\over\la})} \geq \sqrt{\sum_\rho {m_\rho^2\over |\rho|^2}}$$"
Lean form: under RH, every real `c` below the right side is eventually `≤ D(λ)·√log(1/λ)`
(the zero sum taken in `[0, ∞]`, so no summability is presupposed). -/
def Burnol2002_Thm : Prop :=
  RiemannHypothesis → ∀ c : ℝ,
    ENNReal.ofReal c <
        (∑' ρ : NTZ, ENNReal.ofReal ((mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2)) ^ (1 / 2 : ℝ) →
      ∀ᶠ lam in 𝓝[>] (0 : ℝ), c ≤ nbDistC lam * Real.sqrt (Real.log (1 / lam))

/-- **Turán–Nazarov lemma** (F. Nazarov, *Local estimates for exponential polynomials and their
applications to inequalities of the uncertainty principle type*, Algebra i Analiz 5 (1993);
St. Petersburg Math. J. 5 (1994) 663–717), as stated in O. Friedland, Y. Yomdin, *An observation
on Turán–Nazarov inequality*, Studia Math. 218 (2013), arXiv:1107.0039, Theorem 1.1 (TeX l.72–78):
"Let $p(t) = \sum_{k=0}^m c_k e^{\lambda_k t}$ be an exponential polynomial, where
$c_k ,\lambda_k \in \C$. Let $B \subset \R$ be an interval, and let $\Omega \subset B$ be a
measurable set. Then $$\sup_{B} |p| \le e^{\mu_1(B) \cdot \max|\Re \lambda_k|} \cdot \left(
\frac{c \mu_1(B)}{\mu_1(\Omega) }\right)^{m} \cdot \sup_{\Omega} |p|$$ where $\mu_1$ is the
Lebesgue measure on $\R$ and $c>0$ is an absolute constant."
Lean form: `B = [a, b]` (a bounded interval), `Ω` of positive measure, and `K ≥ sup_Ω |p|`. -/
def Nazarov1993_Turan : Prop :=
  ∃ c : ℝ, 0 < c ∧ ∀ (m : ℕ) (coef lam : Fin (m + 1) → ℂ) (a b : ℝ) (Ω : Set ℝ) (K : ℝ),
    a < b → MeasurableSet Ω → Ω ⊆ Icc a b → 0 < volume Ω →
    (∀ t ∈ Ω, ‖∑ k, coef k * Complex.exp (lam k * t)‖ ≤ K) →
    ∀ t ∈ Icc a b, ‖∑ k, coef k * Complex.exp (lam k * t)‖ ≤
      Real.exp ((b - a) * ⨆ k, |(lam k).re|) * (c * (b - a) / (volume Ω).toReal) ^ m * K

/-- **Goldston–Gonek 2007, Theorem 1** (D. A. Goldston, S. M. Gonek, *A note on S(t) and the zeros
of the Riemann zeta-function*, Bull. London Math. Soc. 39 (2007) 482–486, arXiv:math/0511092,
TeX l.85–87 and l.124–133): "We assume the Riemann Hypothesis (RH) throughout this paper. Let $N(t)$
denote the number of zeros $\rho = \frac12 +i\gamma$ of the Riemann zeta-function with ordinates in
the interval $(0, t]$. […] \begin{theorem} Assume the Riemann Hypothesis. Let $t$ be large and
$0 < h \leq \sqrt{t}$. Then we have $$\Big| N(t + h) - N( t) - \frac{h}{2\pi} \log
\frac{t}{2\pi}\Big| \le \big(\frac{1}{2}+o(1)\big)\frac{\log t}{\log\log t}.$$"
(`N(t+h) − N(t)` counts the zeros with ordinates in `(t, t+h]`, with multiplicity.) -/
def GoldstonGonek2007_Thm1 : Prop :=
  RiemannHypothesis → ∀ ε : ℝ, 0 < ε → ∃ t₀ : ℝ, ∀ t : ℝ, t₀ ≤ t → ∀ h : ℝ, 0 < h → h ≤ Real.sqrt t →
    |(Ncount t (t + h) : ℝ) - h / (2 * Real.pi) * Real.log (t / (2 * Real.pi))| ≤
      (1 / 2 + ε) * (Real.log t / Real.log (Real.log t))

/-- **Titchmarsh, Theorem 9.2** (E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*,
2nd ed., rev. D. R. Heath-Brown, Oxford 1986, §9.2): "$N(T+1) - N(T) = O(\log T)$", `N(T)` the number
of zeros with `0 < γ ≤ T`, counted with multiplicity. -/
def Titchmarsh_Thm9_2 : Prop :=
  ∃ A T₀ : ℝ, 0 < T₀ ∧ ∀ T : ℝ, T₀ ≤ T → (Ncount T (T + 1) : ℝ) ≤ A * Real.log T

/-- **Titchmarsh, Theorem 9.7** (same book, §9.7): "There is a constant $A$ such that each interval
$[T, T+1]$ contains a value of $t$ for which $|\zeta(\sigma+it)|^{-1} = O(t^A)$ $(-1 \le \sigma \le 2)$."
(The implied constant is absolute; `ζ(σ+it) ≠ 0` is part of the bound.) -/
def Titchmarsh_Thm9_7 : Prop :=
  ∃ A C T₀ : ℝ, 0 < T₀ ∧ ∀ T : ℝ, T₀ ≤ T → ∃ t ∈ Icc T (T + 1), ∀ σ ∈ Icc (-1 : ℝ) 2,
    riemannZeta (σ + t * I) ≠ 0 ∧ ‖riemannZeta (σ + t * I)‖⁻¹ ≤ C * t ^ A

/-- **Paley–Wiener** (W. Rudin, *Real and Complex Analysis*, 3rd ed., McGraw-Hill 1987, Theorem 19.3):
"Suppose $A$ and $C$ are positive constants and $f$ is an entire function such that
$|f(z)| \le C e^{A|z|}$ for all $z$, and $\int_{-\infty}^{\infty} |f(x)|^2\,dx < \infty$. Then there
exists $F \in L^2(-A, A)$ such that $f(z) = \int_{-A}^{A} F(t) e^{itz}\,dt$ for all $z$." -/
def Rudin1987_Thm19_3 : Prop :=
  ∀ (A C : ℝ) (f : ℂ → ℂ), 0 < A → 0 < C → Differentiable ℂ f →
    (∀ z, ‖f z‖ ≤ C * Real.exp (A * ‖z‖)) → MemLp (fun x : ℝ => f x) 2 volume →
    ∃ F : ℝ → ℂ, MemLp F 2 (volume.restrict (Ioo (-A) A)) ∧
      ∀ z : ℂ, f z = ∫ t in Ioo (-A) A, F t * Complex.exp (I * t * z)

/-- **Davenport, ch. 12** (H. Davenport, *Multiplicative Number Theory*, 3rd ed., GTM 74, Springer
2000, ch. 12, eqs. (10)–(11)): for the Hadamard product
$\xi(s) = e^{A+Bs}\prod_\rho (1-s/\rho)e^{s/\rho}$ (zeros with multiplicity),
"$B = -\sum_\rho \Re(\rho^{-1})$" and "$B = -\tfrac12\gamma - 1 + \tfrac12\log 4\pi$" ($\gamma$
Euler's constant). Lean form: the absolutely convergent sum of `m_ρ·Re(1/ρ)` over the distinct
nontrivial zeros equals `1 + γ/2 − ½ log 4π`. -/
def Davenport_Ch12 : Prop :=
  HasSum (fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * ((ρ : ℂ)⁻¹).re)
    (1 + Real.eulerMascheroniConstant / 2 - Real.log (4 * Real.pi) / 2)

/-- **Titchmarsh, eq. (2.1.5)** (same book, §2.1): "$\zeta(s) = s\int_0^\infty \frac{[x]-x}{x^{s+1}}\,dx$
$(0 < \sigma < 1)$". -/
def Titchmarsh_2_1_5 : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 →
    riemannZeta s = s * ∫ x in Ioi (0 : ℝ), ((⌊x⌋ : ℂ) - x) * (x : ℂ) ^ (-s - 1)

/-- All literature inputs, bundled. -/
structure LitInputs : Prop where
  burnol : Burnol2002_Thm
  turan : Nazarov1993_Turan
  gg : GoldstonGonek2007_Thm1
  count : Titchmarsh_Thm9_2
  heights : Titchmarsh_Thm9_7
  paleyWiener : Rudin1987_Thm19_3
  davenport : Davenport_Ch12
  mellinFract : Titchmarsh_2_1_5

end Riemann.Deep.LeanNBRate

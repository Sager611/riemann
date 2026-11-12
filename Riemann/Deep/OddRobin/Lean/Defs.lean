import Mathlib.Basic.Real.Basic
import Mathlib.NumberTheory.ArithmeticFunction.Misc

/-!
# The sharp odd Robin criterion: definitions (light file)

Washington–Yang, arXiv 2008.04787 (Int. J. Number Theory 17 (2021) 843–870).
The logarithm L and the constant E = e^gamma are parameters; Instance.lean sets L = Real.log,
E = Real.exp Real.eulerMascheroniConstant, which turns ORf into the brief's OR.
-/

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

/-- N = 3^4 * 5^3 * 7^2 * 11 * 13 * ... * 61 (WY, proof of Prop. 3.5, TeX l. 401). -/
def NN : ℕ := 277093800961598968947975375

/-- c0 = 3^4 * 5^3 * 7^2 * 11 * 13 * ... * 67 = 67 N (WY Theorem 1.1, TeX l. 92). -/
def c0 : ℕ := 18565284664427130919514350125

/-- (OR) sigma(n)/n < (E/2) L(L n); with L = Real.log, E = e^gamma this is the brief's OR n. -/
def ORf (L : ℝ → ℝ) (E : ℝ) (n : ℕ) : Prop := (σ 1 n : ℝ) / n < E / 2 * L (L n)

/-- WY's form sigma(n)/(n L(L n)) < E/2 (TeX eq. (RobinOdd), l. 89–91). -/
def WYIneqF (L : ℝ → ℝ) (E : ℝ) (n : ℕ) : Prop := (σ 1 n : ℝ) / (n * L (L n)) < E / 2

/-- Washington–Yang, Theorem 1.1 (label MainThm1, TeX l. 87–93; restated and proved l. 440–460):
RH is equivalent to sigma(n)/(n log log n) < e^gamma/2 for all odd n ≥ c0. RH, L, E are parameters;
the intended instance is RH = Mathlib's RiemannHypothesis, L = Real.log, E = Real.exp gamma. -/
def WYTheoremFor (RH : Prop) (L : ℝ → ℝ) (E : ℝ) : Prop :=
  RH ↔ ∀ n : ℕ, Odd n → c0 ≤ n → WYIneqF L E n

theorem WYIneqF_iff {L : ℝ → ℝ} {E : ℝ} {n : ℕ} (hn : 0 < n) (hll : 0 < L (L n)) :
    WYIneqF L E n ↔ ORf L E n := by
  unfold WYIneqF ORf
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : (0 : ℝ) < n * L (L n) := mul_pos hn' hll
  rw [div_lt_iff₀ h1, div_lt_iff₀ hn']
  constructor <;> intro h <;> linarith

end OddRobin


import Riemann.Deep.OddRobin.Lean.Main
import Riemann.Deep.OddRobin.Lean.LogSeries
import Mathlib.NumberTheory.Harmonic.EulerMascheroni
import Mathlib.NumberTheory.LSeries.RiemannZeta

/-!
# Instantiation with Mathlib's Real.log, Real.exp, γ and RiemannHypothesis

NOT COMPILED: its imports exceed the 290 MB memory cap (Mathlib.Analysis.SpecialFunctions.Log.Basic alone
peaks at 326 MB, Mathlib.NumberTheory.LSeries.RiemannZeta alone at 342 MB). Every field below is a
Mathlib theorem or a LogSeries.lean theorem; the file only plugs them into Main.lean.
-/

namespace OddRobin

open Real ArithmeticFunction
open scoped ArithmeticFunction.sigma

/-- The brief's (OR): sigma(n)/n < (e^gamma/2) log log n. -/
noncomputable def OR (n : ℕ) : Prop :=
  (σ 1 n : ℝ) / n < Real.exp Real.eulerMascheroniConstant / 2 * Real.log (Real.log n)

theorem OR_eq (n : ℕ) : OR n = ORf Real.log (Real.exp eulerMascheroniConstant) n := rfl

theorem logLike_log : LogLike Real.log where
  map_mul x y hx hy := Real.log_mul hx.ne' hy.ne'
  mono x y hx hxy := Real.log_le_log hx hxy
  one_sub_inv_le x hx := Real.one_sub_inv_le_log_of_pos hx
  series_lo a ha m := lsum_le_log ha m
  series_hi a ha m := log_le_lsum ha m

theorem mathlibFacts :
    MathlibFacts Real.log (Real.exp eulerMascheroniConstant) eulerMascheroniConstant where
  logLike := logLike_log
  exp_pos := Real.exp_pos _
  log_exp := Real.log_exp _
  seq_lt n := by
    have := eulerMascheroniSeq_lt_eulerMascheroniConstant n
    simpa [eulerMascheroniSeq] using this
  lt_seq' n hn := by
    have := eulerMascheroniConstant_lt_eulerMascheroniSeq' n
    simpa [eulerMascheroniSeq', hn.ne'] using this

theorem gap_real : ∀ n : ℕ, Odd n → NN < n → n < c0 → OR n := gap_theorem mathlibFacts

theorem sharp_real : ¬ OR NN := sharp mathlibFacts

theorem main_real
    (hWY : WYTheoremFor RiemannHypothesis Real.log (Real.exp eulerMascheroniConstant)) :
    RiemannHypothesis ↔ ∀ n : ℕ, Odd n → NN < n → OR n :=
  main mathlibFacts RiemannHypothesis hWY

end OddRobin


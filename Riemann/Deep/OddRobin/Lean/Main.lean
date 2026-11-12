import Riemann.Deep.OddRobin.Lean.Gap
import Riemann.Deep.OddRobin.Lean.GammaL

/-!
# Main theorems: the gap, sharpness, and RH ⟺ (OR) for every odd n > N

All statements are for a logarithm L, a constant E (= e^gamma) and g (= gamma) satisfying MathlibFacts:
every field is a Mathlib theorem for L = Real.log, E = Real.exp γ, g = γ (see Instance.lean).
-/

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

/-- The facts about Real.log, Real.exp and γ that the proof uses; for L = Real.log,
E = Real.exp γ, g = γ each field is a Mathlib theorem (series_lo/hi: LogSeries.lean, from
Real.hasSum_log_one_add; seq_lt: Real.eulerMascheroniSeq_lt_eulerMascheroniConstant;
lt_seq': Real.eulerMascheroniConstant_lt_eulerMascheroniSeq'). -/
structure MathlibFacts (L : ℝ → ℝ) (E g : ℝ) : Prop where
  logLike : LogLike L
  exp_pos : 0 < E
  log_exp : L E = g
  seq_lt : ∀ n : ℕ, (harmonic n : ℝ) - L ((n : ℝ) + 1) < g
  lt_seq' : ∀ n : ℕ, 0 < n → g < (harmonic n : ℝ) - L n

/-- THE GAP THEOREM (PROOF.md §1): every odd n with N < n < c0 satisfies (OR). -/
theorem gap_theorem {E g : ℝ} (h : MathlibFacts L E g) :
    ∀ n : ℕ, Odd n → NN < n → n < c0 → ORf L E n := fun _ hodd h1 h2 =>
  gap_thm h.logLike h.exp_pos h.log_exp (gamma_bounds h.logLike h.seq_lt h.lt_seq').1 hodd h1 h2

/-- Sharpness: N itself violates (OR). -/
theorem sharp {E g : ℝ} (h : MathlibFacts L E g) : ¬ ORf L E NN :=
  not_OR_NN h.logLike h.exp_pos h.log_exp (gamma_bounds h.logLike h.seq_lt h.lt_seq').2

lemma LL_pos {E g : ℝ} (h : MathlibFacts L E g) {n : ℕ} (hn : NN ≤ n) : 0 < L (L (n : ℝ)) := by
  have hL := h.logLike
  have h1 := LNN_gt_one hL
  have h2 : L (NN : ℝ) ≤ L n := hL.mono _ _ (by norm_num [NN]) (by exact_mod_cast hn)
  exact hL.pos_of_one_lt (by linarith)

/-- PROOF.md §5(b): from Washington–Yang's Theorem 1.1, RH ⟺ (OR) for every odd n > N. -/
theorem main {E g : ℝ} (h : MathlibFacts L E g) (RH : Prop) (hWY : WYTheoremFor RH L E) :
    RH ↔ ∀ n : ℕ, Odd n → NN < n → ORf L E n := by
  have hNc : NN < c0 := by norm_num [NN, c0]
  constructor
  · intro hRH n hodd hn
    rcases lt_or_ge n c0 with hc | hc
    · exact gap_theorem h n hodd hn hc
    · have hw := (hWY.1 hRH) n hodd hc
      exact (WYIneqF_iff (by omega) (LL_pos h (by omega))).1 hw
  · intro hall
    refine hWY.2 (fun n hodd hn => ?_)
    have hNn : NN < n := lt_of_lt_of_le hNc hn
    exact (WYIneqF_iff (by omega) (LL_pos h hNn.le)).2 (hall n hodd hNn)

end OddRobin


/-
Riemann/Deep/Yoshida/MuBounds.lean — certified lower bounds for `2π μ(r) = Re ψ(1/4 + ir/2) − log π`.

* `two_pi_mu_ge`: partial sums of the digamma series (all terms `≥ 0` on the quarter line) give
  `2π μ(r) ≥ −γ̄ − ℓ̄ − (1/4)/((1/4)² + (r/2)²) + Σ_{n<M} (1/(n+1) − (n+5/4)/((n+5/4)² + (r/2)²))`.
* `euler_lt`: `γ < 0.5792` (`H₂₅₆ − log 256`); `log_pi_lt`: `log π < 1.1448`.
* `exp_mem`: two-sided Taylor bounds for `exp` on `|x| ≤ 1`.
-/
import Riemann.Deep.Yoshida.Small

noncomputable section

open Finset

namespace Fayed.Yoshida

open Zeta23

/-- Series lower bound for `2π μ(r)`. -/
theorem two_pi_mu_ge (r : ℝ) (M : ℕ) {g l : ℝ}
    (hg : Real.eulerMascheroniConstant ≤ g) (hl : Real.log Real.pi ≤ l) :
    -g - l - (1 / 4) / ((1 / 4) ^ 2 + (r / 2) ^ 2)
      + ∑ n ∈ range M, (1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2
        + (r / 2) ^ 2)) ≤ 2 * Real.pi * mu r := by
  have hπ := Real.pi_pos
  rw [MuFields.mu_eq, MuFields.re_digamma_vertical (a := 1 / 4) (by norm_num) (by norm_num)]
  have hnn : ∀ n : ℕ, 0 ≤ 1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2
      + (r / 2) ^ 2) := by
    intro n
    have hα : (0 : ℝ) < (n : ℝ) + 1 + 1 / 4 := by positivity
    have h1 : ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (r / 2) ^ 2)
        ≤ 1 / ((n : ℝ) + 1 + 1 / 4) := by
      rw [div_le_div_iff₀ (by positivity) hα]
      nlinarith [sq_nonneg (r / 2)]
    have h2 : 1 / ((n : ℝ) + 1 + 1 / 4) ≤ 1 / ((n : ℝ) + 1) :=
      one_div_le_one_div_of_le (by positivity) (by linarith)
    linarith
  have hsum := (MuFields.summable_re_terms (a := 1 / 4) (by norm_num) (by norm_num) (r / 2)).sum_le_tsum
    (range M) (fun n _ => hnn n)
  have e : 2 * Real.pi * (1 / (2 * Real.pi) * (-Real.eulerMascheroniConstant
      - 1 / 4 / ((1 / 4) ^ 2 + (r / 2) ^ 2) + ∑' n : ℕ, (1 / ((n : ℝ) + 1)
        - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (r / 2) ^ 2)))
      - Real.log Real.pi / (2 * Real.pi))
      = -Real.eulerMascheroniConstant - 1 / 4 / ((1 / 4) ^ 2 + (r / 2) ^ 2)
        + ∑' n : ℕ, (1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2
          + (r / 2) ^ 2)) - Real.log Real.pi := by
    field_simp
  rw [e]
  linarith

theorem euler_lt : Real.eulerMascheroniConstant < 0.5792 := by
  have h := Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' 256
  have e : Real.eulerMascheroniSeq' 256 = (harmonic 256 : ℝ) - 8 * Real.log 2 := by
    simp only [Real.eulerMascheroniSeq', OfNat.ofNat_ne_zero, ite_false]
    rw [show ((256 : ℕ) : ℝ) = 2 ^ 8 by norm_num, Real.log_pow]
    push_cast
    ring
  rw [e] at h
  have hl := Real.log_two_gt_d9
  have hH : (harmonic 256 : ℝ) < 6.124345 := by
    have h2 : harmonic 256 < (6124345 / 1000000 : ℚ) := by
      norm_num [harmonic, Finset.sum_range_succ]
    have h3 := (Rat.cast_lt (K := ℝ)).2 h2
    push_cast at h3
    linarith
  linarith

theorem log_pi_lt : Real.log Real.pi < 1.1448 := by
  rw [Real.log_lt_iff_lt_exp Real.pi_pos]
  have h := Real.sum_le_exp_of_nonneg (x := 1.1448) (by norm_num) 9
  have hs : (3.1417 : ℝ) < ∑ i ∈ range 9, (1.1448 : ℝ) ^ i / (i.factorial : ℝ) := by
    norm_num [Finset.sum_range_succ, Nat.factorial]
  linarith [Real.pi_lt_d4]

/-- Two-sided Taylor bounds for `exp` on `|x| ≤ 1`. -/
theorem exp_mem {x : ℝ} (hx : |x| ≤ 1) (n : ℕ) (hn : 0 < n) :
    ∑ m ∈ range n, x ^ m / (m.factorial : ℝ) - |x| ^ n * ((n.succ : ℝ) / ((n.factorial : ℝ) * n))
      ≤ Real.exp x ∧
    Real.exp x ≤ ∑ m ∈ range n, x ^ m / (m.factorial : ℝ)
      + |x| ^ n * ((n.succ : ℝ) / ((n.factorial : ℝ) * n)) := by
  have h := abs_le.1 (Real.exp_bound hx hn)
  constructor <;> linarith [h.1, h.2]

end Fayed.Yoshida

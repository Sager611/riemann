import Mathlib

/-!
# Round 21 — learning from failed RH proofs

| id | claim | flaw | class |
|----|-------|------|-------|
| R21.0 | Fayed 2022–26 (root) | L'Hôpital on gaps false; identities σ-blind | false lemma + σ-blindness |
| R21.1 | X.-J. Li 2008 (withdrawn) | positive-definiteness assumed | false positivity lemma |
| R21.2 | de Branges 2004– | Conrey–Li negative value | false positivity lemma |
| R21.3 | Atiyah 2018 (Todd function) | real-zero claim ≡ RH, unproved | false lemma |
| R21.4 | Eswaran 2016– | deterministic signs treated as independent | randomness as proof |
| R21.5 | Bender–Brody–Müller 2017 | self-adjointness/domain not shown | false self-adjointness |
| R21.6 | Wong 2023 (Nyman–Beurling inequality) | finite-N to limit leap | false limiting inference |
| R21.7 | numerical Báez-Duarte / Li evidence | finite checks | numerics as proof |
| R21.8 | "FE symmetry forces Re = 1/2" | pairs zeros, does not fix them (symmetric_offline_example) | symmetry ≠ fixedness |

Salvage: none of these yields an RH route; each failure class is now a filter for new branches.
-/

namespace Fayed.R21

/-- R21.8: a polynomial with ζ's symmetry s ↦ 1 - s (and real coefficients) whose zeros
1/3 ± i are off the line. -/
theorem symmetric_offline_example :
    let P : ℂ → ℂ := fun s => ((s - 1 / 3) ^ 2 + 1) * ((s - 2 / 3) ^ 2 + 1)
    (∀ s, P (1 - s) = P s) ∧ P (1 / 3 + Complex.I) = 0 ∧ (1 / 3 + Complex.I).re ≠ 1 / 2 := by
  refine ⟨fun s => by ring, ?_, ?_⟩
  · simp only
    have : (1 / 3 + Complex.I - 1 / 3) ^ 2 + 1 = 0 := by
      rw [show (1 / 3 + Complex.I - 1 / 3 : ℂ) = Complex.I by ring, Complex.I_sq]; ring
    rw [this, zero_mul]
  · simp

end Fayed.R21

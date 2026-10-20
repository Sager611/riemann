import Mathlib

/-!
# Round 47 — horizontal monotonicity / modulus criteria: local-looking, globally RH

| criterion | proven region (unconditional) | enlarging it | verdict |
|-----------|-------------------------------|--------------|---------|
| Sondow–Dumitrescu: |ξ| increasing in σ > 1/2 | near σ ≥ 1 | ⇒ zero-free region | CIRCULAR |
| Matiyasevich–Saidak–Zvengrowski: |ζ| monotone for σ < 1/2 | σ < 0, large |t| | stronger than zero-freeness | DEAD |
| Spira / Dixon–Schoenfeld: |ζ(1-s)| > |ζ(s)| | strict where ζ(s) ≠ 0 (|χ| < 1) | strictness everywhere ⇔ RH | CIRCULAR |
| Speiser / Levinson–Montgomery: ζ' ≠ 0 in 0 < σ < 1/2 | bounded heights only | ⇔ RH | CIRCULAR |
| Lagarias Re ξ'/ξ > 0 | Re s > 1 | positivity ⇒ zero-free | CIRCULAR |
| log|ξ| convexity (Hinkkanen, Garunkštis) | σ ≥ 1 | ⇒ zero-free | CIRCULAR |
| |ξ| minimised on the line (Saidak–Zvengrowski) | σ ≥ 1 | ⇒ zero-free | DEAD |
| arg monotonicity / Gram's law | fails infinitely often | — | DEAD |

Each pointwise criterion fails at a zero, so proving it in a region proves that region zero-free:
"local" criteria are zero-free-region statements in disguise. The unconditional half (Lean):
the functional equation with |χ(s)| < 1 gives |ζ(s)| < |ζ(1-s)| wherever ζ(s) ≠ 0.
-/

namespace Fayed.R47

/-- If z = c·w with ‖c‖ < 1 and z ≠ 0, then ‖z‖ < ‖w‖. -/
theorem modulus_lt {z w c : ℂ} (h : z = c * w) (hc : ‖c‖ < 1) (hz : z ≠ 0) : ‖z‖ < ‖w‖ := by
  have hw : w ≠ 0 := by
    rintro rfl
    exact hz (by simp [h])
  rw [h, norm_mul]
  exact mul_lt_of_lt_one_left (norm_pos_iff.mpr hw) hc

end Fayed.R47

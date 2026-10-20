/-!
# Round 28 — 16 micro-variants of the root argument: all σ₀-blind

| id | variant | outcome |
|----|---------|---------|
| R28.1a | complement domain (B̄₁ < 0) | sign of F'/G' flips; F/G → 0 unchanged |
| R28.1b | oscillating G ~ u^{-1} cos 2πu | F/G → (1-2σ₀)/2: consistent, no second computation |
| R28.2a | p > 1 | F/G oscillates, unbounded |
| R28.2b | p = 1 | F/G bounded oscillating |
| R28.3a | T = u^k | error control only |
| R28.3b | η → 0 | contour hits the pole; undefined |
| R28.4a | period averages | O((1-2σ₀)M^{-2}) since ∫B̄₂ = 0 |
| R28.4b | integers + Stolz | ΔF/ΔG ~ (1-2σ₀)M^{p-1}/(12p) → 0 |
| R28.5a | all Euler–Maclaurin orders | same leading terms |
| R28.5b | ∂/∂s at s₀ | ζ'(s₀) absorbs the identity |
| R28.6a | double zero | extra s-derivative identity only |
| R28.6b | conjugate zero s̄₀ | redundant |
| R28.7a | Perron line past z = 0 | residue ζ(s₀) = 0 |
| R28.7b | unequal weights p ≠ q | B̄₁(u)(u^p - u^q), σ₀-free |
| R28.8a | Poisson summation | rewrites B̄ₖ only |
| R28.8b | eta function instead of ζ | nonzero factor, same zero |

The root's neighbourhood is exhausted: every local modification is consistent with an off-line zero.
-/

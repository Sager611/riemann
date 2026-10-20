/-!
# Round 90 — triage of older arXiv claims

| target | key step | flaw |
|--------|----------|------|
| Wünsche, arXiv 1607.04551 (2016) | Ξ(z) = ∫₀^∞ Ω(u) cosh(uz) du with Ω ≥ 0 non-increasing, rapidly decreasing ⇒ all zeros on z = iy (second mean value theorem) | the general theorem is FALSE: K(u) = e^{-u²}(1 + u⁴/4) is positive, non-increasing on [0, ∞) (K' = u e^{-u²}(−2 + u² − u⁴/2) < 0) and rapidly decreasing, yet its transform has four zeros off the axis (Fayed.R38.quartic_factor_pos; Gaussian-moment identity cited) |
| Stenger, arXiv 1708.01209 | G(z) = Γ(z)(1 − 2^{1−z})ζ(z) symmetry + zero correspondence | symmetry ≠ fixedness; Davenport–Heilbronn |
| "Bernoulli Sums, Euler–Maclaurin, and RH" (2022) | — | no such arXiv record (misattributed); nearest: Mantzakouras–López 2112.11228 (Bernoulli series for ξ, corrected) |

Next: prove the Gaussian-moment identity in Lean to make the counterexample (and the R38 log-concavity
counterexample) fully machine-checked.
-/

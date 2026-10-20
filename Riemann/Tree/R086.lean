/-!
# Round 86 — certify zeros of other functions with our machinery

| target | cost | value |
|--------|------|-------|
| DH off-line zero near 0.8085 + 85.6993i (rectangle count = 1) — would make dh_offline unconditional | high | high → queued (Orca) |
| certified DH evaluation via L(s, χ) mod 5 at height 85.7 | high | high (prerequisite) |
| first zero of Dirichlet beta L(s, χ₋₄) ≈ 1/2 + 6.0209i | medium | high |
| first zero of L(s, χ₋₃) ≈ 1/2 + 8.0397i | medium | high |
| off-line zero of the Epstein zeta for x² + 5y² | high | high |
| first DH zeros on the line | medium | medium |
| off-line zero of a partial sum Σ_{n≤N} n^{-s} | medium | high (no FE, but Dirichlet series) |
| interval widths at height 85.7 vs 33 | medium | high (feasibility check first) |

Elementary remark: 1 + 2^{-s} vanishes at s = (2k+1)πi / log 2 (Re s = 0), a Dirichlet polynomial with exact
off-line zeros; it has no functional equation, so it only kills arguments using neither FE nor Euler product.
-/

/-!
# Round 87 — machine-checkable RH analogs

| analog | RH-type statement | Mathlib has | cost |
|--------|-------------------|-------------|------|
| Hasse bound | |#E(F_q) − q − 1| ≤ 2√q | WeierstrassCurve, points | high |
| Gauss sums | |g(χ)| = √p | gaussSum, gaussSum_mul_gaussSum_eq_card, star_gaussSum_eq | low (identity) |
| Ihara zeta, Ramanujan graphs | RH ⇔ Ramanujan | adjacency matrices, charpoly | high (no Ihara API) |
| elliptic curve over F₂ | roots of 1 − aT + 2T² have modulus √2 | localPolynomial | high |
| Kloosterman sums | |K(a,b;p)| ≤ 2√p | Gauss/Jacobi sums | high |
| Fermat curves via Jacobi sums | Jacobi sums have modulus √q | jacobiSum | high |
| Hashimoto edge matrix | spectral modulus | charpoly, spectrum | high |
| function-field zeta of F_q[T] | trivial RH | RatFunc, Dedekind domains | high |

The cheapest proven RH analog in Lean is the Gauss-sum modulus |g(χ)|² = p (weight 1/2 purity).
-/

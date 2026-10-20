/-!
# Round 102 — RH for curves over finite fields as a formal target

| item | formalized? | Lean cost |
|------|-------------|-----------|
| Hasse's theorem |#E(F_q) − q − 1| ≤ 2√q | none known (Lean, Isabelle, Coq, HOL) | high |
| Manin 1956 elementary proof (Chahal–Osserman, Monthly 2008) | none known | high |
| Stepanov–Bombieri–Schmidt (Weil bound for curves) | none known | high |
| Kloosterman bound |K(a,b;p)| ≤ 2√p | none known | high |
| Hasse–Weil zeta of E/F_q (rationality, FE) | none known | high |
| Mathlib gaps | Frobenius endomorphism, degree map, isogenies, dual isogeny | — |
| brute force p ≤ 50 by kernel decide | feasible but uninformative | med |
| Weil via Riemann–Roch (Bombieri 1973) | Mathlib lacks curves / divisors / Riemann–Roch | high |

Mathlib already has the Gauss-sum analog (|g|² = p). Next big formal target (deferred until a slot frees):
Hasse's theorem via Manin / Chahal–Osserman, or Stepanov's method plus the tensor-power trick.
-/

/-!
# Round 97 — semi-local and small-support positivity theorems

| result | statement | Lean feasibility |
|--------|-----------|------------------|
| Yoshida 1992 | Weil positivity when the autocorrelation is supported in (−log 2, log 2) (our L < log 2) | medium (target; we have L ≤ 0.55) |
| Connes–Consani 2021 (archimedean place) | archimedean positivity on the same range, via the semi-local trace formula | medium |
| Connes–Consani, Spectral triples and ζ-cycles | no proved window; numerical truncations | low |
| Burnol 2001 | semi-local positivity for S = {∞, p}, support in [−log p, log p] | medium |
| Bombieri 2000 | positivity for support ≤ log 2 (prime-free), via L¹ minimisation | medium |
| windows log 2 < L < log 3 | only the n = 2 prime term appears (−2 log 2 · g(log 2)); no sign in general | decomposition easy, positivity open |
| Connes–Moscovici 2022 | archimedean form positive for support in [−2π, 2π] (prolate functions) | medium |
| numerical windows | Zhu 2026, Mori (log 15), ours 1.75 | medium (finite PSD + interval bounds) |

Note: Connes–Moscovici's archimedean-only positivity up to 2π does not include the prime terms, which
start at log 2; full Weil positivity beyond log 2 needs control of −2 log 2 · g(log 2) and later terms.
-/

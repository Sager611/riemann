/-!
# Round 80 — the Hadamard / 1/ζ bottleneck, broken into lemmas

| lemma | Mathlib has | effort |
|-------|-------------|--------|
| Jensen ⇒ n(r) ≪ r log r for ξ | Jensen formula (MeromorphicOn.circleAverage_log_norm), isolated zeros | 300–1000 |
| Σ 1/|ρ|^{1+ε} < ∞ | summability comparison lemmas | 30–180 |
| genus-1 canonical product entire, order ≤ 1 | multipliable_one_add_of_summable, locally uniform limits | 150–300 |
| Hadamard factorization (order 1) | none | 1000–2000 |
| ξ'/ξ = B + Σ (1/(s−ρ) + 1/ρ) | HasProdLocallyUniformlyOn, deriv of limits | 300–600 |
| Littlewood: log ζ ≪ (log t)^{2−2σ+ε} under RH | three-lines (Hadamard.lean), Jensen | 1500–3000 |
| |1/ζ(σ+it)| ≪ t^ε for σ ≥ 1/2 + δ under RH | Phragmén–Lindelöf | 1000–3000 |
| Perron for μ ⇒ M(x) = O(x^{1/2+ε}) | moebius, LSeries identities | 800–1500 |

Total ≈ 5k–12k lines: a multi-agent formalization project; queued after the Yoshida extension.
-/

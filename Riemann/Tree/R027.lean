/-!
# Round 27 — Mathlib infrastructure gap map for RH work

| item | Mathlib status | missing | effort (lines) |
|------|----------------|---------|----------------|
| Hadamard factorization | partial: analyticOrderAt, tprod/log tools | finite order, canonical products, theorem | 1500–3000 |
| argument principle | partial: meromorphicOrderAt, logDeriv, Cauchy integrals, Jensen | ∮ f'/f = Σ orders | 300–600 (circles) |
| Jensen / Borel–Carathéodory / Landau lemma | Jensen ✓, Complex.borelCaratheodory ✓ | Landau log-derivative lemma | 200–400 |
| Phragmén–Lindelöf in strips | PhragmenLindelof.vertical_strip ✓ | ζ convexity bound itself | 200–900 (with Γ bounds) |
| Perron's formula | LSeries/Mellin identities | truncated Perron with errors | 500–900 |
| explicit formula for ψ | ψ, Λ, zero set discreteness, logDeriv FE ✓ | the formula | 1000–2500 |
| AFE / Riemann–Siegel / Hardy Z | completed ζ, FE ✓ | AFE, θ, Z | 1500–3000 |
| complex Stirling in strips | real Stirling only | complex log Γ asymptotics | 500–1000 |

Next formal target chosen: ζ convexity bound (would discharge Fayed's Inputs.convexity).
-/

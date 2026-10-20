/-!
# Round 105 — formalizing the consequences of a zero-free half-plane Re s > θ (θ = 7/8 now known)

| consequence | reference | Lean cost |
|-------------|-----------|-----------|
| log ζ, ζ'/ζ ≪ (log t)^A for σ ≥ θ + ε | Titchmarsh §14.2 (Borel–Carathéodory, three circles) | high |
| ψ(x) = x + O(x^{θ+ε}) | Titchmarsh §17.3, Ingham III.4 (truncated Perron + contour shift) | high |
| M(x) = O(x^{θ+ε}) | Littlewood, Titchmarsh §14.25 (Perron for 1/ζ) | high |
| π(x) = li(x) + O(x^{θ+ε}); a prime in [x, x + x^{θ+ε}] | partial summation | high |
| ζ(σ+it) ≪ t^ε for σ > (1+θ)/2 = 15/16 | Titchmarsh §5.1 / 14.2 | high |

With θ = 7/8 these are the first unconditional power savings (x^{7/8+ε}) in the prime number theorem and for M(x),
assuming the OpenAI theorem. Formal target after the port: ψ(x) = x + O(x^{7/8+ε}) (reuse perron33, convexity_holds).
Also: arXiv 2609.29898 (Kazin–Kadyrov) disproves "Yang's conjecture" on a tempered ξ, not RH — no audit needed.
-/

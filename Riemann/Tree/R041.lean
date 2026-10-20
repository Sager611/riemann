/-!
# Round 41 — next formal targets on the critical line

| target | route | missing pieces | effort (Lean lines) |
|--------|-------|----------------|---------------------|
| Hardy–Littlewood N₀(T) ≥ cT | Z moments | effective moment asymptotics | 1.5k–3k |
| Selberg: positive proportion | mollifier | approximate FE, mean values | 8k–15k |
| Levinson: ≥ 1/3 | mollifier + Levinson's device | shifted convolutions | 4k–8k |
| Riemann–von Mangoldt N(T) | argument principle | zero-counting API, Stirling | 1.5k–3k |
| Riemann–Siegel / Hardy Z | approximate FE | Stirling, Euler–Maclaurin | 2k–4k |
| Turing's method up to T₀ | certified numerics | interval arithmetic for ζ | 15k–30k |
| first zero γ₁ ≈ 14.1347 | certified numerics + counting | same | 3k–8k |
| de la Vallée Poussin zero-free region | 3-4-1 + growth bounds | log-derivative estimates | 2k–4k |

Prior art found this round: Alpöge–Furman (arXiv 2608.13637, 2026) prove that more than two thirds of
the zeros are simple and on the critical line, formalized in Lean 4 (github.com/anthropics/formal-math,
zeta23). It implies Hardy's theorem and the proportion results above, so those targets are not new.
Open formal targets that remain: certified numerics (Turing's method, the first zero) and zero-free regions.
-/

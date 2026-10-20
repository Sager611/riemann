/-!
# Round 24 — positive results on the critical line, as Lean targets

(Scout provider was degraded: 8 of 16 scouts returned empty; 5 answers below suffice.)

| id | result | key input | Lean size (est.) | status |
|----|--------|-----------|------------------|--------|
| R24.1 | Hardy 1914: infinitely many zeros on the line | Ξ(t) = ∫Φ(u)cos(ut)du, positivity + moments of Φ | 800–1500 lines | PROMISING formal target → Orca deep dive |
| R24.2 | Hardy–Littlewood 1921: N₀(T) ≫ T | Z(t), AFE, moments, sign changes | 1500–3000 | HARD |
| R24.3 | Selberg 1942: positive proportion | mollifiers, mean values | 20k–50k | HARD |
| R24.4 | Levinson/Conrey ≥ 1/3, 2/5 | mollified moments of ζ' | 20k–50k | HARD |
| R24.5 | de la Vallée Poussin region 1 - c/log t | Landau lemma, Borel–Carathéodory | 2k–5k | HARD |
| R24.6–8 | N(T) formula, Turing method, Pólya route | — | — | no answer (provider) |

Formal filter object built this round: Davenport–Heilbronn function with its Riemann-type
functional equation, machine-checked (Riemann/Deep/DH/FE.lean).
-/

# Proof tree

Root: RH via Fayed (arXiv 2209.01890v55).

- R0 paper as written — DEAD at (55): Thm 1 false. [Theorem1.lean](Riemann/Fayed/Theorem1.lean), [Paper.lean](Riemann/Fayed/Paper.lean) `fayed_gap`
- R1 repair (55) — EXHAUSTED: 5 dead, 5 circular. [R001.lean](Riemann/Tree/R001.lean)
- R2 global σ-sensitive input passing the Davenport–Heilbronn filter — EXHAUSTED: 7 dead, 3 circular. [R002.lean](Riemann/Tree/R002.lean)
- R3 push partial theorems / finite faces — 6 dead, 3 circular, 1 promising. [R003.lean](Riemann/Tree/R003.lean)
  - R3.1 Weil positivity beyond window log 2 (prime 2 enters) — FRONTIER (R4 + Orca deep dive)

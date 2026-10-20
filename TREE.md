# Proof tree

Root: RH via Fayed (arXiv 2209.01890v55).

- R0 paper as written — DEAD at (55): Thm 1 false. [Theorem1.lean](Riemann/Fayed/Theorem1.lean), [Paper.lean](Riemann/Fayed/Paper.lean) `fayed_gap`
- R1 repair (55) — EXHAUSTED: 5 dead, 5 circular. [R001.lean](Riemann/Tree/R001.lean)
- R2 global σ-sensitive input passing the Davenport–Heilbronn filter — EXHAUSTED: 7 dead, 3 circular. [R002.lean](Riemann/Tree/R002.lean)
- R3 push partial theorems / finite faces — 6 dead, 3 circular, 1 promising. [R003.lean](Riemann/Tree/R003.lean)
  - R3.1 Weil positivity beyond window log 2 — OPEN: R4 (7 dead, 1 promising) [R004.lean](Riemann/Tree/R004.lean); Orca deep dive [Deep/Weil](Riemann/Deep/Weil/)
- R5 breadth: frontiers, disproof, meta — 5 dead, 2 circular, 1 open-insufficient. [R005.lean](Riemann/Tree/R005.lean)
- R6 "Fayed + X": root identity + one global input — EXHAUSTED: 8 dead, 2 circular. [R006.lean](Riemann/Tree/R006.lean)
- R7 transfer from proven RH analogs — 9 dead, 1 open on a weaker goal (FUP strip). [R007.lean](Riemann/Tree/R007.lean): RH ⟺ QuasiRH ½ proved
- R8 invention: 10 new mechanisms — 5 dead, 5 circular. [R008.lean](Riemann/Tree/R008.lean)
- R9 gap map of RH-equivalents — no gap closable by known methods. [R009.lean](Riemann/Tree/R009.lean); formal backbone: Mertens ⇒ zero-free (Orca deep dive, Riemann/Deep/Mertens)
- R10 radical reformulations — 6 dead, 4 circular. [R010.lean](Riemann/Tree/R010.lean); Redheffer det formalization (ultraworker)
- R11 attack the obstructions — all heuristic/theorem; none broken. [R011.lean](Riemann/Tree/R011.lean)
- Formal backbone (machine-checked): Redheffer det R_n = M(n) [Det](Riemann/Deep/Redheffer/Det.lean); Mertens ⇒ zero-free [ZeroFree](Riemann/Deep/Mertens/ZeroFree.lean); RH ⟸ Mertens / Redheffer bounds [Backbone](Riemann/Deep/Backbone.lean); Weil windows (in progress)
- R12 unusual angles — 7 dead, 3 circular. [R012.lean](Riemann/Tree/R012.lean)
- R13 exotic equivalents — 3 dead, 7 circular. [R013.lean](Riemann/Tree/R013.lean); Nyman–Beurling sufficiency formalization (ultraworker)
- R14 FRONTIER (next)

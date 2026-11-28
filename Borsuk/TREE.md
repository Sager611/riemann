# Borsuk in ℝ⁴: approach tree

Root B: every bounded X ⊂ ℝ⁴ with diam X > 0 splits into 5 parts of diameter < diam X. Known: 5 ≤ b(4) ≤ 6 (arXiv:2611.17375, 2611.16775; R001).
Status: open | active | dead (one line in NOGO.md). Rank = chance × payoff (1–5 each).

## Disproof: finite X whose diameter graph G(X) is not 5-colourable (≥ 11 points: R001 L4)
| id | node | status | rank | rounds |
|---|---|---|---|---|
| D1 | subsets of uniform 4-polytope vertex sets and root systems (A4, B4, F4, H4 families; F4/B4/C4 roots) | active | 1×5 | R002, R005 |
| D2 | lattice subsets (ℤ⁴, D4, A4, A4*) in bounded regions, every squared-distance threshold | active | 1×5 | R005 |
| D3 | two-distance sets (≤ 10 points) | dead | – | R001, R002 |
| D4a | spherical sets, circumradius > D/√2 | dead | – | R001 |
| D4b | spherical sets, circumradius ≤ D/√2 | dead | – | R006 |
| D5 | algebraic embeddings (strongly regular graphs, Frankl–Wilson counting) at d = 4 | dead | – | R006 |
| D6 | OpenAI projector construction pushed to d = 4 | dead | – | R001 |
| D7 | small counterexamples (11–16 points; includes three-distance sets) | dead | – | R001, R006 |
| D8 | complete joins A ∨ B (Lenz-type) | dead | – | R001 |
| D9 | optimisation of β₅/diam over finite sets (Gale-type multistart search) | dead | – | R004, R006 |

## Proof
| id | node | status | rank | rounds |
|---|---|---|---|---|
| P5 | push the b(4) ≤ 6 proofs to 5 parts (5 simplex cones for r < 3/10; Tolmachev's outer bodies H(u) above) | dead | – | R002, R003, R006 |
| P1 | universal cover of diameter-1 sets in ℝ⁴ cut into 5 parts of diameter < 1 | dead | – | R006 |
| P2 | constant-width bodies: extend the smooth / centrally symmetric / revolution / belt cases | dead | – | R004, R006 |
| P3 | reduction to finitely many configurations + exact computation (= P5's method) | dead | – | R006 |
| P4 | topological (Lyusternik–Schnirelmann, genus, equivariant maps) | dead | – | R001, R004, R006 |

## Rounds
| round | node | verdict |
|---|---|---|
| R001 | root: literature | b(4) ∈ {5, 6}; P5 opened; D4a, D6, D8 dead; D3 dead pending g(4) = 10 |
| R002 | P5: audit of the three proofs | no stated obstruction to 5; same contact-triangle reduction; g(4) = 10 ⇒ D3 dead |
| R003 | P5: obstructions in H(u) | no clique, topological or volume obstruction; 40k-point far graph of the hardest-looking H(u) 5-colourable |
| R004 | root: new routes + regular-simplex case | no new route; sets containing a regular unit simplex split into 5 (0.874) |
| R005 | D1, D2: exhaustive colouring checks | see file |
| R006 | P5 barrier; P1–P4, D4b, D5, D7, D9 | β₅(H(1,0,0,0)) ≥ 0.9699 certified; all dead |

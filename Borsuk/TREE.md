# Borsuk in ℝ⁴: approach tree

Root B: every bounded X ⊂ ℝ⁴ with diam X > 0 splits into 5 parts of diameter < diam X. Known: 5 ≤ b(4) ≤ 6 (arXiv:2611.17375, 2611.16775; R001).
Status: open | active | dead (one line in NOGO.md). Rank = chance × payoff (1–5 each).

## Disproof: finite X whose diameter graph G(X) is not 5-colourable (≥ 11 points: R001 L4)
| id | node | status | rank | rounds |
|---|---|---|---|---|
| D1 | subsets of uniform 4-polytope vertex sets and root systems (A4, B4, F4, H4 families; F4/B4/C4 roots) | active | 1×5 | R002 |
| D2 | lattice subsets (ℤ⁴, D4, A4, A4*, icosians) in bounded regions, every squared-distance threshold | open | 1×5 | |
| D3 | two-distance sets (≤ 10 points) | dead | – | R001, R002 |
| D4a | spherical sets, circumradius > D/√2 | dead | – | R001 |
| D4b | spherical sets, circumradius ≤ D/√2 (vertex links in 3-flats, edge links in 2-flats, triangle links ≤ 2 points) | open | 1×5 | |
| D5 | algebraic embeddings (strongly regular graphs, Frankl–Wilson counting: Kahn–Kalai, Bondarenko, Jenrich–Brouwer) at d = 4 | open | 1×5 | |
| D6 | OpenAI projector construction pushed to d = 4 | dead | – | R001 |
| D7 | small counterexamples (11–16 points; includes three-distance sets): 6-critical graphs obeying L1–L3, then realizability | open | 1×5 | |
| D8 | complete joins A ∨ B (Lenz-type) | dead | – | R001 |

## Proof
| id | node | status | rank | rounds |
|---|---|---|---|---|
| P5 | push the b(4) ≤ 6 proofs to 5 parts: 5 simplex cones for r < 3/10 (done); 5-partitions of Tolmachev's outer bodies H(u) for ρ ∈ [3/10, 2/5] | active | 2×5 | R002, R003 |
| P1 | universal cover of diameter-1 sets in ℝ⁴ cut into 5 parts of diameter < 1 (Langi ≤ 14-facet covers) | open | 1×5 | |
| P2 | constant-width bodies: extend smooth / centrally symmetric / revolution / belt cases | open | 1×5 | |
| P3 | reduction to finitely many configurations + exact computation | open | 1×5 | |
| P4 | topological (Lyusternik–Schnirelmann, equivariant maps; Musin arXiv:2609.09688) | open | 1×5 | |


## Rounds
| round | node | verdict |
|---|---|---|
| R001 | root: literature | b(4) ∈ {5, 6}; P5 opened; D4a, D6, D8 dead; D3 dead pending g(4) = 10 |
| R002 | P5: audit of the three proofs | no stated obstruction to 5; same contact-triangle reduction; g(4) = 10 ⇒ D3 dead |
| R003 | P5: obstructions in H(u) | no clique, topological or volume obstruction; 40k-point far graph of the hardest-looking H(u) 5-colourable |

# Borsuk in ℝ⁴: approach tree

Root B: every bounded X ⊂ ℝ⁴ with diam X > 0 splits into 5 parts of diameter < diam X. Known: 5 ≤ b(4) ≤ 6 (arXiv:2611.17375, 2611.16775; R001).
Status: open | active | dead (one line in NOGO.md). Rank = chance × payoff (1–5 each).

## Disproof: finite X whose diameter graph G(X) is not 5-colourable (≥ 11 points: R001 L4)
| id | node | status | rank | rounds |
|---|---|---|---|---|
| D1 | subsets of uniform 4-polytope vertex sets and root systems (A4, B4, F4, H4 families; F4/B4/C4 roots) | dead | – | R005 |
| D2 | lattice subsets (ℤ⁴, D4, A4, A4*) in bounded regions, every squared-distance threshold | dead | – | R005 |
| D3 | two-distance sets (≤ 10 points) | dead | – | R001, R002 |
| D4a | spherical sets, circumradius > D/√2 | dead | – | R001 |
| D4b | spherical sets, circumradius ≤ D/√2 | dead | – | R006 |
| D5 | algebraic embeddings (strongly regular graphs, Frankl–Wilson counting) at d = 4 | dead | – | R006 |
| D6 | OpenAI projector construction pushed to d = 4 | dead | – | R001 |
| D7 | small counterexamples (11–16 points; includes three-distance sets) | dead | – | R001, R006 |
| D8 | complete joins A ∨ B (Lenz-type) | dead | – | R001 |
| D9 | optimisation of β₅/diam over finite sets (Gale-type multistart search) | dead | – | R004, R006 |
| D10 | D1/D2 remainders with a stronger 2-level solver (incremental CaDiCaL proposer, kissat checks, stabilizer symmetry breaking) | open | 1×5 | |
| D11 | colour-forcing rigidity gadgets (K5/K4 sharing patterns; non-topological, non-counting obstructions) | open | 1×5 | |
| D12 | 2-level search inside the tightest outer bodies: diameter-≤1 subsets with large β₅ (fuzzy graphs, CEGAR over compatible sets) | open | 1×5 | R022 |

## Proof
| id | node | status | rank | rounds |
|---|---|---|---|---|
| P5 | push the b(4) ≤ 6 proofs to 5 parts (5 simplex cones for r < 3/10; Tolmachev's outer bodies H(u) above) | dead | – | R002, R003, R006 |
| P1 | universal cover of diameter-1 sets in ℝ⁴ cut into 5 parts of diameter < 1 | dead | – | R006 |
| P2 | constant-width bodies: extend the smooth / centrally symmetric / revolution / belt cases | dead | – | R004, R006 |
| P3 | reduction to finitely many configurations + exact computation (= P5's method) | dead | – | R006 |
| P4 | topological (Lyusternik–Schnirelmann, genus, equivariant maps) | dead | – | R001, R004, R006 |
| P6 | case analysis by contact type (3, 4, 5 contacts) with explicit 5-partitions of the all-contact outer bodies Ω = B(0,√ρ) ∩ ⋂ B(pᵢ,1) | open | 1×5 | R007, R020–R022 |
| P7 | analytic large-r argument: contact-aligned cones/power cells + moment bound M ⪯ (1/2 − r)I (Chalmers/ABPR Lemma 6) | open | 2×5 | |
| P8 | fast C certifier + adaptive branch-and-bound that adds contacts only where the 3-contact form is near-critical | active | 2×5 | R015–R022 |
| P9 | rotation-avoidance near ρ = 3/10: rotate the 5 simplex cones away from the near-contact diametral pairs (thin-shell lemma) | open | 1×4 | R008, R009 |
| P10 | Hadwiger partition with canonical normals + repair of bad (singular–singular) diametral pairs | open | 2×5 | R010 |
| P11 | X-dependent partitions: case splits on whether X meets critical cells of the outer body (adds balls B(c, 1+δ) or removes cells), leaf partitions certified (cert2) | active | 2×5 | R023 |
| P12 | fan conjecture: an X-dependent simplex fan (rotation Q, apex c) splits every diameter-1 set; prove via the diametral-chord family of the completion (blocking measure, adversarial game) | active | 2×5 | R023 |

## Rounds
| round | node | verdict |
|---|---|---|
| R001 | root: literature | b(4) ∈ {5, 6}; P5 opened; D4a, D6, D8 dead; D3 dead pending g(4) = 10 |
| R002 | P5: audit of the three proofs | no stated obstruction to 5; same contact-triangle reduction; g(4) = 10 ⇒ D3 dead |
| R003 | P5: obstructions in H(u) | no clique, topological or volume obstruction; 40k-point far graph of the hardest-looking H(u) 5-colourable |
| R004 | root: new routes + regular-simplex case | no new route; sets containing a regular unit simplex split into 5 (0.874) |
| R005 | D1, D2: exhaustive colouring checks | no counterexample where complete (regular polytopes, root systems, A4/B4/F4 families, 8 H4 orbits, lattices to s ≤ 12–32); rest stalls |
| R006 | P5 barrier; P1–P4, D4b, D5, D7, D9 | β₅(H(1,0,0,0)) ≥ 0.9849 certified; all dead |
| R007 | P6: all-contact margins; no-unit-pair relaxation; Petersen critical pairs | room 3–6% near ρ = 0.3 for all contact types |
| R008 | P6: density bounds on S³; thin-shell lemma | no density obstruction; cone failures need near-contact diametral pairs |
| R009 | P9: normal cones, rotation avoidance | at most a small analytic extension 3/10 + ε₀; middle range needs adapted partitions |
| R010 | P10: Hadwiger colouring with canonical normals | monochromatic edges only among bad edges (θ_x + θ_y ≥ 48.2°); necessary condition for counterexamples |
| R011 | P6/P8: relaxed partitions on samples | Jung corner passes the band test; u = (0.6,1,0,0) fails it (tabu) |
| R012 | P8: exact SAT for the band graph | SAT but cluster-inconsistent (components 1.02–1.06); needs a consistent encoding |
| R013 | P8: cluster-consistent relaxed SAT | Jung-corner certificate stays UNSAT with 37 far region pairs |
| R014 | P8: a fourth contact | room stays 3–6% (8k samples); certified computation needs adaptive boxes + ms-cost C certifier |
| R015 | P8: certifier prototype | exact pivoting enclosure works; too loose (1.0199 vs 0.9667) and slow (279 s) |
| R016 | P8: certifier v2 (support functions + duality) | rigorous in principle; Python 28 ms/box too slow; port to C |
| R017 | P8: certifier v2 in C | ball 0.29 ∩ simplex cone certified < 0.985 in 48k boxes (0.4 ms/box) |
| R018 | P8: Jung-corner piece | Ω₅ aligned cone certified < 0.88 |
| R019 | P8: end-to-end cones + certifier | 2 of 5 cells certified; sampled diameters underestimate by ≈ 1% |
| R020 | P8: extreme contact configurations | type 3 at ρ = 1/3 is the crux (β₅ sample (0.97, 0.99]); regular tetrahedron ≤ 0.95; cones fail for type 3 |
| R021 | P8: SAT-certified brackets | tabu NOT ≠ UNSAT (corrected); extreme type 3 ≥ 0.97, type 4 at 0.36 ∈ [0.95, 0.97] |
| R022 | P8: relaxation, explicit partitions, Jung stability | relaxation useless at type-3 extreme; LP-fitted cell partitions 1.03–1.10; type-5 near-regular only for ρ ≥ 0.395 |
| R023 | P12: fans on diameter-1 sets | optimised simplex fans give 0.88 (4D) on sets whose outer bodies need ≥ 0.97; fan conjecture; blocking measure ≤ 0.31% per chord at ρ ≤ 1/3 |

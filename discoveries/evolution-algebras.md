# Evolution algebras of graphs: when is A(G) ≅ A_RW(G)?

[Back to DISCOVERIES](../DISCOVERIES.md)

**Definitions.** For a finite connected graph G, the evolution algebra A(G) has a basis (e_v) with e_v² = Σ_{w∼v} e_w; its random-walk version A_RW(G) has e_v² = Σ_{w∼v} e_w / deg v.

**Claim.** Over every field 𝕜 of characteristic 0, A(G) ≅ A_RW(G) if and only if G is regular, or G is bipartite biregular with part degrees d₁, d₂ such that d₁/d₂ is a cube in 𝕜. Characteristic p is also classified.

**Consequences.**
- Solves an open problem of arXiv 2610.22248 (Oct 2026). As posed (over any field), its statement holds over R and C but fails over Q (P₃, K_{1,3}, K_{2,3}).
- The Cadavid–Rodiño Montoya–Rodriguez conjecture (2021, over R; known before only for non-singular graphs) holds; its "moreover" clause fails over C for every singular G that is neither regular nor biregular.

**Key steps.** The natural bases of a non-degenerate evolution algebra are exactly unions of orthogonal bases of the twin blocks (via Boudi–Cabrera–Siles Molina 2022); an all-ones identity shows twin rotations never change the scalars.

**Checks.** All 995 connected graphs with at most 7 vertices, over C, R and Q.

**Status.** Paper proof; not refereed. Lean (evolean, 12 Nov 2026): FORMALIZED for characteristic 0 — EvoAlg.theoremA, both directions (G finite connected simple, n ≥ 2: A(G) ≅ A_RW(G) ⟺ G regular or bipartite biregular with d₁/d₂ ∈ (𝕜^×)³), with Theorem I, Lemma 2, Corollary 4.1, Theorem N (BCS22 Thm 2.4(i) proved, not assumed) and Corollary A (over ℚ the cube condition is needed: K_{1,2}, K_{1,3}, K_{2,3}); no named hypotheses, sorry-free, standard axioms, every compile under the 290 MB watchdog; the ℝ/ℂ instances (imports above the cap) and the characteristic-p classification are not formalized; reviewer FAITHFUL ([LEAN.md](../Riemann/Deep/EvoAlg/LEAN.md)).

**Review.** Proofs correct (the quotation of Boudi–Cabrera–Siles Molina 2022 verified against the source).

**Novelty.** Nothing states any criterion for singular graphs ([novelty](../Riemann/Deep/EvoAlg/NOVELTY.md)).

**Files.** [report](../Riemann/Deep/EvoAlg/REPORT.md), [note](../Riemann/Deep/EvoAlg/NOTE.md), [verify.py](../Riemann/Deep/EvoAlg/verify.py)

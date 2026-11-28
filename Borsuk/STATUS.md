WORKING
UTC: 2026-11-28 11:21
Node: P8 middle range: find convex (or union) 5-cell partitions of Ω(contacts) for ρ ∈ [0.3, 0.39]; certifier ready (cert2.c); D10 in background
Last round: R018 (Jung-corner piece certified < 0.88)
Next: todo #19: partition search in the middle range with the certifier in the loop (optimise power-diagram sites against cert2's worst direction, warm-started), starting at the all-contact type-3 bodies near ρ = 0.31 (room 3–6%); then types 4, 5

## Decisions
- 2026-11-28 00:00 Tree seeded with 8 disproof nodes (D1–D8) and 4 proof nodes (P1–P4); ranks set after R001.
- 2026-11-28 00:00 Colouring checks: own exact-integer code + CaDiCaL via python-sat (uv); a 5-colouring found is checked edge by edge; UNSAT only → DRAT + drat-trim (build from source when needed).
- 2026-11-28 00:00 Scouts: concurrency 4, sequential web searches each (brief: ≤ 4 searches at once).
- 2026-11-28 00:40 R001 found b(4) ≤ 6 (Chalmers, Tolmachev, Nov 2026): new top node P5; the question is now 5 vs 6.
- 2026-11-28 00:40 Killed D4a, D6, D8 by Kupavskii's S³ theorem and own lemmas L1–L5 (R001); D3 pending g(4) = 10.
- 2026-11-28 02:45 P5 split: r < 3/10 done by 5 simplex cones; the open part is ρ ∈ [3/10, 2/5] with Tolmachev's normal form H(u) (his verification pipeline is the template).
- 2026-11-28 02:45 Falsify first: test whether some H(u) needs 6 parts (CEGAR on boundary samples) before writing any verification code.
- 2026-11-28 02:45 Own C tabu colourer (Borsuk/code/tabucol.c); numpy without BLAS for 4-dim dot products and fixed array shapes (Accelerate and macOS malloc kept growing buffers past 300 MB).
- 2026-11-28 03:00 P5 declared dead: certified ≥ 0.9699 at the Jung corner of the 4-parameter normal form, no convex 5-cell partition below 1 found; 6–10-parameter refinements judged beyond reach (R006).
- 2026-11-28 03:00 macOS: run every numpy job with MallocLargeCache=0 (the large-allocation cache kept growing-size buffers past 300 MB).
- 2026-11-28 03:00 Route check before stopping: 8 scouts (R004) found no new route; the only concrete candidate (4D Meissner-type bodies) is 5-partitionable (0.874).
- 2026-11-28 04:03 D1, D2 closed with explicit unresolved remainders (2-level search stalls on large local sets); declared dead because the remaining checks cannot finish with our solver and carry no structure beyond D7.
- 2026-11-28 04:03 STOPPED.
- 2026-11-28 08:42 Reopened (Adrian, 28 Nov): STOPPED is not allowed; dead nodes spawn new nodes. New nodes P6 (contact-type case analysis with explicit 5-partitions), P7 (analytic large-r argument), P8 (C certifier + adaptive branch-and-bound adding contacts only where needed), D10 (stronger 2QBF solver for D1/D2 remainders), D11 (rigidity gadgets). Todo goal set; one task in progress at a time.
- 2026-11-28 09:00 Relaxation adopted for all future partition searches: parts may contain pairs at distance > 1 (unions of mutually far clusters); only unit pairs are forbidden.
- 2026-11-28 10:15 D10 (background): CaDiCaL proposer + 3600 s limit now resolves H4:01 thresholds that timed out before (s = 181.4, 185.4: UNSAT); pipeline continues (d10_cegar.out).

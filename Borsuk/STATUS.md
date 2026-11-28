WORKING
UTC: 2026-11-28 04:10
Node: D1, D2 (last computations); all other nodes dead
Last round: R006 (P5 barrier certified: β₅(H(1,0,0,0)) ≥ 0.9699)
Next: finish d1v2 + d1_cegar (H4 orbits) and d2 (A4 to s = 22); write R005; then STOPPED unless a check fails

## Decisions
- 2026-11-28 00:00 Tree seeded with 8 disproof nodes (D1–D8) and 4 proof nodes (P1–P4); ranks set after R001.
- 2026-11-28 00:00 Colouring checks: own exact-integer code + CaDiCaL via python-sat (uv); a 5-colouring found is checked edge by edge; UNSAT only → DRAT + drat-trim (build from source when needed).
- 2026-11-28 00:00 Scouts: concurrency 4, sequential web searches each (brief: ≤ 4 searches at once).
- 2026-11-28 00:40 R001 found b(4) ≤ 6 (Chalmers, Tolmachev, Nov 2026): new top node P5; the question is now 5 vs 6.
- 2026-11-28 00:40 Killed D4a, D6, D8 by Kupavskii's S³ theorem and own lemmas L1–L5 (R001); D3 pending g(4) = 10.
- 2026-11-28 02:45 P5 split: r < 3/10 done by 5 simplex cones; the open part is ρ ∈ [3/10, 2/5] with Tolmachev's normal form H(u) (his verification pipeline is the template).
- 2026-11-28 02:45 Falsify first: test whether some H(u) needs 6 parts (CEGAR on boundary samples) before writing any verification code.
- 2026-11-28 02:45 Own C tabu colourer (Borsuk/code/tabucol.c); numpy without BLAS for 4-dim dot products and fixed array shapes (Accelerate and macOS malloc kept growing buffers past 300 MB).
- 2026-11-28 04:10 P5 declared dead: certified ≥ 0.9699 at the Jung corner of the 4-parameter normal form, no convex 5-cell partition below 1 found; 6–10-parameter refinements judged beyond reach (R006).
- 2026-11-28 04:10 macOS: run every numpy job with MallocLargeCache=0 (the large-allocation cache kept growing-size buffers past 300 MB).
- 2026-11-28 04:10 Route check before stopping: 8 scouts (R004) found no new route; the only concrete candidate (4D Meissner-type bodies) is 5-partitionable (0.874).

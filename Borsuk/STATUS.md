WORKING
UTC: 2026-11-28 17:50
Node: P5 barrier certificate for H(0.8,0,0,0) at t = 0.9899 (2,488 exact points kissat-UNSAT; drat-trim needs more than 300 MB → shrinking the set: Borsuk/code/p8_shrink.py → Borsuk/certs/h08_0989.txt, log Borsuk/.scratch/hbar3.out); P12 fan conjecture (spherical version next)
Last round: R028 (merging lemma; 6-piece merging fails)
Next: finish the certificate (DRAT verified) and write R029; then spherical fan conjecture (all-contact sets X = R·D, D ⊂ S³ with angular diameter θ_c(ρ), 0 ∈ conv D) as the 3-dimensional core of P12

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
- 2026-11-28 12:01 Type 3 (0 inside a contact triangle on a great circle, 3 parameters) is the crux of the contact-type split; work on its extreme body first, both as proof bottleneck and as disproof lead.
- 2026-11-28 12:46 Correction recorded (R021): CEGAR "NOT" verdicts were tabu failures; all future brackets use the kissat fallback.
- 2026-11-28 14:09 Pivot to P12 (fan conjecture): diameter-1 sets are ≥ 9% easier than their outer bodies, so the proof must use X's own diametral structure; fans are the simplest X-dependent partitions.
- 2026-11-28 17:50 Barrier certificates use symmetry breaking by a clique of the far graph (sound: any colouring can be permuted to colour the clique 1..|Q|); shrink point sets before DRAT so that drat-trim fits in 300 MB.

# prefixlean STATE
Start 2026-11-12 17:02 UTC. Stop rule: 19:32 write-up; hard end 20:02.
Plan: (0) defs (1) B.1 n-cycle (2) parity (3) kernel counterexamples (4) blocks/orbits necessity (5) conditional sufficiency.
Status: reading inputs; measuring import memory next.
17:25 Defs, CycleA/B/C compile (peaks 249/269/264/254 MB). B.1 n-cycle (transitive form), Lemma 1, Lemma 2 proved. Next: Parity.lean.
17:30 Parity compiles (278 MB): sign r_i = (-1)^(i/2), G ≤ A_n iff ..., eq_top_of_alt, parity_line. Next: IsCycleP (Closure import) + Small (kernel words).
17:35 IsCycleP (267 MB) + Small (279 MB) compile: P.IsCycle, support univ, G_eq_top_of_swap; G 7 4 2 = G 11 6 4 = G 13 10 2 = ⊤ by decide +kernel; not_conjecture4Literal. Next: fetch TeX statements; Blocks (necessity).
17:45 TeX fetched to /tmp/prl (delete at end): Conj4 = 2511.16959v3 l.884-886; LZZW Thm 3.1 = 2610.30892 l.221-223; Jones Thm 1.2 = arXiv 1209.5169 TeX l.54-72. Writing BlocksA (gcd invariant set, interval blocks).
17:43 BlocksA (261 MB), BlocksB (272 MB) compile: gcd≠1 intransitive, interval blocks (α|m), n even parity blocks, mod-4 ⇒ ≤ A_n, line_necessary (no hypotheses). Next: Main.lean (Theorem A from named hyps OffLineLZZW, Prop5Primitive [ours, open], JonesLemma7).
17:47 Main (269 MB) compiles: theoremA (from OffLineLZZW, Prop5Primitive, JonesLemma7), theoremA_necessary (from OffLineLZZW only), line_B4. Next: Prim.lean = Lemma 3 (blocks are residue classes of psi) reducing Prop5Primitive to arithmetic Prop5Core (★_d ⇒ α|n).
17:58 PrimA/PrimB1/PrimB2 compile (≤270 MB; one PrimB draft killed at 303 MB, split): Lemma 3 (blocks = residue classes of psi), star_of_block, Prop5Core (open, arithmetic Steps 1-4) ⇒ Prop5Primitive; theoremA_core. Next: split Jones hypothesis (JonesK0 with Coll predicate + Lemma7) + case (a) affine proof + n=11,23 words.
18:04 Small23 (n=23 words) compiled rc=0 but /usr/bin/time footprint 466 MB (memcap sampled 250, missed spike) = CAP BREACH; file deleted, not part of deliverable. Instead: Jones.lean with JonesK0(Coll) [Jones Thm 1.2 k=0, (a),(c) with fixed-point facts] + Lemma7(Coll).
18:06 Jones.lean (274 MB) compiles: JonesK0 Coll, Lemma7 Coll, line_alt_le_jones, theoremA_jones, line_B4_jones. Next: CheckA/CheckB (#print axioms), then LEAN.md, reviewer.
18:21 Reviewer round 1 ISSUES (evidence gap: section hypotheses not shown; docstring overstated) → fixed docstring, supplied #check signatures; round 2 FAITHFUL. LEAN.md: header, files, statements, hypotheses written. Clean rebuild (verify.sh) running in background → verify.log.
19:08 Clean rebuild done (verify.log VERIFY_DONE): 17/17 rc=0, max footprint 281 MB, 49 #print axioms standard, scan clean. LEAN.md assembled (PARTIAL; order per brief). Oleans and test files deleted; /tmp/prl deleted. DONE.

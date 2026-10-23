# Quaternary covering codes: K_4(10,3) <= 768 and K_4(11,3) <= 1920

**Statement.** There is a quaternary code of length 10 with 768 words, and one of length 11 with 1920 words, such that every word
is within Hamming distance 3 of a codeword.

**Before.** Keri's tables (rev. 2009, update log to 2011): 313 <= K_4(10,3) <= 832 and 849 <= K_4(11,3) <= 2048. Both upper bounds
are computer results of Bertolo, Di Pasquale and Santisi (2004-06). arXiv 2610.03760 (Oct 2026) improved K_4(8,2) to 256 with the
octacode and K_4(6,2) to 48, and reports no other improvement among the open quaternary cells with n <= 11.

**Construction.** O = octacode, the Z4-linear code [I_4 | A] of arXiv 2610.03760 (256 words, length 8, covering radius 2).
K = {c in O : c_8 = 0} is Z4-linear with 64 words. Pad K with m zero coordinates (m = 2 or 3) and take a union of its cosets:
- length 10: 12 cosets, 12 * 64 = 768 words; length 11: 30 cosets, 30 * 64 = 1920 words.
- The 768 code also has the form C = {(c_1..c_7, a, y) : c in O, c_8 = a, (a, y) in B} with |B| = 12 (header of the code file).
A union of cosets x_i + K' has covering radius 3 iff the representatives S satisfy S + Sigma_3 = Z_4^n / K', where Sigma_3 is the
image of the radius-3 ball; the quotient has 16384 (n=10) or 65536 (n=11) elements.

**Method.** Tabu search over the coset representatives with the kernel fixed, in C (Ostergard-style moves: pick an uncovered
syndrome, move one representative by a unit vector so that it covers it, exact delta). 12 cosets take about 1 s, 30 cosets about
30 s. Random Z4-linear kernels of the same size, and random linear extensions of the octacode, do much worse (150-240 uncovered
syndromes at 12 cosets), so the octacode structure is essential.

**Limits inside the construction (not lower bounds for K_4).**
- Octacode split on its last coordinate plus 2 new coordinates (this includes the amalgamated direct sum with any code of
  length 3): an exact hitting-set search shows that 12 blocks are needed, so 768 is optimal there; the ADS with an 8-word code
  (512 words) fails.
- All 120 distinct 64-word subcodes {c in O : v.c = 0 mod 4} as kernels: 11 cosets (704 words) were never reached.
- Length 11 with 29 cosets (1856 words): best 19 uncovered syndromes.

**Review.** verify.py (a pure-Python bitset closure over all 4^n words) gives covering radius exactly 3 for both codes in under 1 s.
An independent breadth-first search agrees; words at distance 0/1/2/3: 768/21504/253184/773120 (n=10) and
1920/62656/843648/3286080 (n=11). Both codes are invariant under translation by K x 0^m (checked).

**Novelty.** Newer tables still list the old values: Florath's covering-codes-lean data (Sep 2026) has 832 and 2048, and its
Lean-proved bounds are only 4096; the thiagopatzdorf/Matematica ledger (Oct 2026) has 832 and 2048. arXiv API searches for
covering-code papers of 2025-2026 find only 2610.03760 for q = 4. NOVELTY.md records a deeper check of the 768 code (Keri's live
site, OpenAlex, Zenodo, GitHub, evidence API); the 1920 code got the lighter checks above.

**Significance.** Modest table improvements (-7.7% and -6.3%), a fairly direct consequence of the octacode of 2610.03760 plus
computer search.

**Files.** codes/K4_10_3_768.txt, codes/K4_11_3_1920.txt, codes/ref_K4_6_2_48_maharaj.txt, codes/ref_K4_8_2_256_octacode.txt,
verify.py, REPORT.md, NOVELTY.md.

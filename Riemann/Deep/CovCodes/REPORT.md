RECORD: K_4(10,3) <= 768 (was 832); K_4(11,3) <= 1920 (was 2048)
Certified: verify.py (pure Python, checks all 4^n words, < 1 s per code) gives covering radius exactly 3 for both codes;
an independent breadth-first search agrees. Words at distance 0/1/2/3: 768/21504/253184/773120 (n=10) and
1920/62656/843648/3286080 (n=11).
Previous bounds: Keri's tables (4-5_tables.pdf, 2009; log to 2011): 313 <= K_4(10,3) <= 832 and 849 <= K_4(11,3) <= 2048. Both
upper bounds are by Bertolo-Di Pasquale-Santisi (2004-06). Newer tables still list them: Florath's covering-codes-lean data
(Sep 2026) has 832 and 2048, its Lean-proved bounds only 4096; the thiagopatzdorf/Matematica ledger (Oct 2026) has 832 and 2048.
arXiv 2610.03760 (Oct 2026) changes only K_4(6,2) and K_4(8,2); arXiv API searches (2025-2026) show nothing newer for q=4.
NOVELTY.md holds a separate, deeper novelty check of the 768 code; the 1920 code got the table and arXiv checks above.
Construction: K = {c in octacode : c_8 = 0} (64 words, Z4-linear; octacode = [I|A] of 2610.03760), padded with zeros.
768 = 12 cosets of K x 0^2 in Z_4^10; 1920 = 30 cosets of K x 0^3 in Z_4^11. Coset representatives: tabu search in C.
Limits of the construction (not lower bounds for K_4):
- octacode split on its last coordinate + 2 new ones (includes every amalgamated direct sum with a length-3 code): an exact
  search shows 12 blocks are needed, so 768 is the best this form gives.
- all 120 distinct 64-word subcodes of the octacode as kernels: 704 words never reached (good kernels stall at 33 uncovered).
- 1856 words for (11,3): best 19 uncovered syndromes; 1919 words by plain tabu started from the 1920 code: best 147.
NO RECORD for the brief's main targets: K_4(6,2) stays 48 and K_4(8,2) stays 256.
- K_4(6,2) = 47: orbit, coset and plain tabu searches never closed the gap (plain tabu: >= 139 uncovered words).
- K_4(8,2) = 255: octacode minus one word plus repair fails; the needed density is below the near-perfect octacode's
  (lower bound 251). 47 and 255 are odd, so neither can be a union of cosets of a nontrivial additive code.
Other cells tried, no improvement: (10,2) 4096; (10,4) published 208, 192 claimed in the thiagopatzdorf/Matematica ledger,
our 192-word coset code misses 576 words; (11,2) 15872; (11,4) 512; (11,5) 128; (9,2) 1024; (9,3) 256; (7,2) 128;
(7,3) 32; (6,3) 14.
Significance: modest table improvements (-7.7% and -6.3%) that follow from the 2610.03760 octacode plus computer search.
Files: codes/K4_10_3_768.txt, codes/K4_11_3_1920.txt, codes/ref_* (reference codes of 2610.03760), verify.py, NOTE.md, NOVELTY.md.

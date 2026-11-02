NO RECORD: no tabulated cell of K_4(n,R) (n <= 11) or K_3(n,R) (n <= 13), R <= 4, was beaten. Side result: K_4(12,4) <= 1536 (length 12 is not tabulated for q=4; the best direct sum of published codes gives 1792).
Certified: verify.py (bitset closure over all 4^12 words, 1.7 s) gives covering radius exactly 4 for codes/K4_12_4_1536.txt; an independent C check (ball enumeration) finds 0 uncovered words and no redundant codeword.

Where the tables stand (upper bounds; * = changed in 2026; sources: Keri's tables, thiagopatzdorf ledger, CovCodes):
- q=4, n=6..11: R=1 256, 992, 3456, 12288, 49152, 131072 | R=2 48*, 128, 256*, 1024, 4096, 15872 | R=3 14, 32, 96, 256, 768*, 1920* | R=4 (n=7..11) 10, 27*, 64, 192*, 512.
- q=4, n=12 (not tabulated; product bounds): R=1 524288, R=2 63488, R=3 6144, R=4 1536 (ours; was 1792).
- q=3 (Keri 2011; the 2026 ledger shows no ternary change): R=1 n=6..12: 73, 186, 486, 1269, 3645, 9477, 27702 | R=2 n=6..13: 17, 34, 81, 219, 555, 729 (exact), 2187, 6561 | R=3 n=7..13: 12, 27, 54, 105, 243, 657, 1215 | R=4 n=9..13: 18, 36, 81, 175, 335. All other cells in range are exact.

What was tried:
1. Product rules over all pairs of table entries and the 2026 codes: a plain direct sum never beats the table. The amalgamated direct sum (ADS) formula would beat it in K_4(10,3), (11,3), (11,4), (8,4) and in several ternary cells, but the octacode and the 48-word code are not normal, so each case needs an exact search.
2. Exact split searches (the code is the base code split into classes, each class joined with a small code B; B chosen by an exact set cover over distance profiles; SCIP, or local search when too big):
   - octacode split on 1 coordinate: (10,3) needs 12 blocks (768, optimal); (9,3) 256 optimal; (11,4) needs >= 8 blocks = 512 (hand proof); (11,3) needs 26-30 blocks (best 30 = 1920). All 135 labellings of the octacode fall into 5 classes; none helps (10,3), (11,3) or (11,4). Split on 2 coordinates for (11,4): best 624.
   - Maharaj 48-word code: any such split needs at least 12*K_4(1+m, R-2) words, which is at least the table value in every cell ((8,3) 96, (9,3) 288, (10,3) 768, (9,4) 84, (10,4) 192, (11,4) 576).
   - 192-, 27-, 10-, 768- and 1920-word codes plus coordinates at the same radius: SCIP proves every block is needed (for example 192 -> 768 for (11,4), 27 -> 108 for (9,4)).
3. Ternary: subcode unions of the Golay [11,6] code and of the Hamming [4,2] code, with coordinates kept, dropped or added: 62 exact instances, 32 more ruled out by a structural lower bound. Nothing reached the table (for example K_3(12,3) via Golay: optimum 729 > 657).
4. Short local searches: random-kernel unions of cosets (50 cell/kernel pairs; F_3-linear and additive F_4 kernels) reach only table values or worse (K_4(7,3) 32, K_4(7,2) 128, K_3(8,3) 27). A length-11 linear extension of the octacode plus greedy patching gives 802 words for (11,4), far above 512.
Not done (out of time): extended Golay, Hamming [13,10], hexacode and [5,3]_4 bases; the ADS search for K_4(8,4) built on K_4(4,1) (24 < 27 would be a record); exact 2-coordinate splits (local search too weak, and SCIP needs more than 300 MB).
Significance: nothing for the tables. The 1536 bound is in a cell nobody tabulates (no novelty search done beyond the tables and the ledger).
Files: NOTE.md, codes/K4_12_4_1536.txt, verify.py (q = 3 or 4).


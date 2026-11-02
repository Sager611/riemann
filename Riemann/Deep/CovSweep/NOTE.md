# CovSweep: the q=3 and q=4 covering-code tables against this year's new codes

**Statement.** No tabulated cell of K_4(n,R) (n <= 11) or K_3(n,R) (n <= 13) with R <= 4 was improved. There is a quaternary code of
length 12 with 1536 words in which every word lies within Hamming distance 4 of a codeword: K_4(12,4) <= 1536.

**Before.** Keri's q=4 table stops at n = 11, and so does the thiagopatzdorf/Matematica ledger (which also imports Florath's data).
The best bound for K_4(12,4) from direct sums of published codes is 1792 (octacode (+) K_4(4,2) = 256*7). Other product rules give
1920 (K_4(11,3)) and 2048 (4*K_4(11,4)). An amalgamated direct sum would need normal factors, and the octacode is not normal.
Lower bound: sphere covering, 4^12 / V(12,4) = 16777216 / 46666, so at least 360.

**Construction.** O = octacode (generator [I|A] of arXiv 2610.03760; ../CovCodes/codes/ref_K4_8_2_256_octacode.txt).
label(x) = x_7 + x_8 mod 4 splits O into 4 classes of 64 words. C = {(x, y) : x in O, y in Z_4^4, (label(x), y) in B}, where
|B| = 24 (6 per label). This gives 24 * 64 = 1536 words. B is listed in the code file header.

**Method.** For a base code A with classes A_l and m new coordinates, the code {(x, y) : x in A_l, (l, y) in B} has covering radius
<= R iff, for every distance profile P(w) = (d(w, A_l))_l and every z, some (l, y) in B has P_l(w) + d(z, y) <= R. Only maximal
profiles matter. This is an exact set cover on B (SCIP, or local search with row weights when large). Splitting on coordinates
that stay in the codeword gives the amalgamated form (ADS and its generalisations). The same tool, with structural lower bounds,
gave the negative results in REPORT.md:
- octacode ADS cells: (10,3) 768 optimal; (11,4) at least 512; (11,3) at least 1664.
- Maharaj 48: cannot beat the ADS bound in any cell.
- the R=4 codes: every block is needed.
- ternary Golay/Hamming subcode unions: nothing below the table.

**Review.** verify.py (pure Python bitset closure over all 4^12 = 16,777,216 words) gives covering radius exactly 4 in 1.7 s.
An independent C ball enumeration finds 0 uncovered words, and removing any single codeword uncovers a word. The file was rebuilt
from O and B, and it equals {(x, y) : x in O, (label(x), y) in B} exactly.

**Novelty.** The cell is not in Keri's tables or in the 2026 ledger; no literature search beyond those was done. This is not a table
record: it only beats a product bound in an untabulated cell.

**Files.** codes/K4_12_4_1536.txt, verify.py (q = 3 or 4; header n=, R=, q=), REPORT.md.


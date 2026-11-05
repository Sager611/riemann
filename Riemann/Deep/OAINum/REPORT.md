NO FLAW FOUND: OpenAI's five number-theory preprints (deterministic factorization, Gaussian moat, Ostmann, prime-argument polynomial differences, totient fibers): every finitely checkable lemma we tested holds, and the two Lean proofs scan clean.

# Audit: OpenAI's number-theory proofs (agent "oainum", Nov 2026)
Source: github.com/openai/math, preprints/<name>/build (TeX, single-file downloads) and lean/ (single files, GitHub code search). Method as in ../OAIProofs/REPORT.md: each lemma with a finite instance retested with our own code, exact arithmetic, one single-threaded job at a time. Scripts are in this folder; their output lines are in OUTPUT.txt.

| preprint | lemmas tested | ranges | verdict |
|---|---|---|---|
| Deterministic factorization (4 Oct) | Sec. 3 table: the Gauss-period field U_q has degree q; omega_2, omega_q have the exact primary orders | 54 entries, p <= 103, q = 2, 3, 5, 7 | holds |
| | Lemma 2.2 root extraction (incl. the q dividing d eigenvector branch); distinct-degree splitting | 720 random q-th powers, q = 3, 5, 7, p <= 10111, moduli of degree <= 8; 1000 moduli | holds |
| | Prop. 2.3 even-degree splitter (orientation masks, row counts) | 4425 split F: all of p <= 17, N <= 4; random N <= 8, p <= 65537 | always a proper factor; exactly one orientation per pair |
| | Sec. 4 geometry on Y^3 = F(X) (genus 1, m = 3): Lemma 4.3 rationality, Lemma 4.2 label formula, Prop. 4.6 forced separation | every cubic over F_7 (35) and F_13 (286); 40 over F_19, 12 over F_31; 10 with K = F_25 | all 27 lambda^3-divisions of [3e] rational over k; separation for all 10,341 admissible d_3 (first split at t = 1, 2, 3 in ratio 18:6:3); formula = brute-force labels 14,937 times |
| | Thm 6.1 norm solver, implemented as written (local maximal orders, section ansatz, fiber at infinity, trace kernel, rank-one idempotent, b^q, Ve = he) | 28 random promised norms with q = 3 (p = 331, 1009, 10009; bundle shapes 3, 2+1, 1+2, 1+1+1), 3 with q = 5 | N(h) = G exactly every time; non-norms rejected; section space unchanged for M + 6 |
| | Prop. 7.1 promise: the normalized pushforward is a norm iff the class is lambda-divisible | genus 1 over F_p, p = 37, 43, 61, 67: 12 curves, 705 classes | 235/235 divisible classes solved, 470/470 others rejected |
| Gaussian moat (26 Sep) | Lean proof lean/OAI/NumberTheory/GaussianMoat (41 files, Mathlib only) | comparator statement; full scan | faithful; no sorry, axiom, native_decide, extern, unsafe; not built |
| | Lemma 4.1 (A <= 4 pi D^2 abs(E-E), n+1 <= (9/pi)(R+2/3)(2W+2/3)), cube edge-isoperimetry, Lemmas 3.2, 3.3 | 3000 walks, D^2 <= 9, n <= 320; all sets r <= 4, random r <= 10; 20k laws each | hold (smallest abs(E-E)/A = 0.72) |
| Ostmann (24 Sep) | Lean proof lean/OAI/NumberTheory/Ostmann (7716 files) | comparator statements; code search over all 7720 indexed files | faithful; final theorems carry no hypotheses; cited inputs (zero density, quadratic large sieve, Siegel bound, ...) are proved props; no sorry, native_decide, implemented_by, unsafe, opaque; "axiom" only in comments; not built |
| Prime-argument polynomial differences (5 Oct) | Lemma 2.1 auxiliary family (integrality, T_l(Dx-i) = T_lD, h_l(Dx-i) = lambda(D) h_lD, D divides lambda(D) <= D^k); Lemma 2.3 uniform lifting classes | 7 prime-intersective h, incl. (x^2-2)(x^2-17)(x^2-34); l <= 48, D <= 12; p <= 7, l <= 30 | hold |
| Totient fibers (24 Sep) | Type II exact determinant parametrization (bijection); block-sieve polynomials (V <= 1 <= U, level <= z^(4h+2)) | e <= 12, m, h <= 60; z = 13, 30, h = 2, 4, n <= 60000 | hold; without coprimality 4258 converse candidates fail, as the paper says |

## Hand check
- Factorization: all eleven sections read. The odd-degree mechanism is sound: lifted ramification classes sum to a lambda^m-division of [Ne], and constant labels at every test would force the lattice divisor Ne one lambda-power deeper than it is. The norm solver's local orders, vector-bundle argument and fiber idempotent are correct, and the max-of-partial-sums step gives (1-sigma)V = E. Not run: the generic ideal-based divisor arithmetic (Sec. 5; standard Hess/Khuri-Makdisi techniques) and the complexity bookkeeping (Sec. 10). The unconditional theorem rests on the companion's uniform Hecke zero-free strip (delta = 1e-6, "Primitive roots for every admissible integer base", Thm 1.2), which we audited earlier (NOGO.md: no false step). The reduction alone gives deterministic polynomial-time factoring under GRH.
- Prime-argument paper: its analytic input is OpenAI's 7/8 half-plane (Lean-verified; ported in Riemann/Deep/QRH). Tuple-law section: the transposition-distance and no-opposite-edge claims were checked by hand.
- Totient paper: the analytic core (transference, ideal kernels, Type II) was read but not re-derived. openai/math's Lean covers only the trivial case eps >= 1.
- Ostmann and Gaussian moat: introductions and statements read; correctness rests on the Lean proofs (not built, per brief).

## Prior reports
openai/math has issues and discussions disabled; no public refutation of these papers found.

## Caveat
One q = 5 norm-solver run peaked at 597 MB (cap 500 MB). Row reduction was then made in place; later runs stayed at or below 310 MB, the others below 200 MB.

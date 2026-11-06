# Imaginary quadratic fields with class number at most 5000 (node N11)

**Claim.** Every imaginary quadratic field K = Q(√−d) (d = |disc K|) with class number h(K) ≤ 5000 is known, without RH or GRH: there are 102,746,940 of them, and the largest d is 6,912,348,787 (h = 4967). Every h ≤ 5000 occurs. Counts, largest d and SHA-256 of the sorted lists for N = 1000, 2000, 2500, …, 5000 are in [REPORT](../Riemann/Deep/ImqClass/REPORT.md); the per-h table (h, number of fields, largest d) is perh_h5000.txt.

**Before.** h ≤ 2000 unconditionally (arXiv 2610.21511, Oct 2026, with a Lean-checked zero gap); h ≤ 3000 conditional on Kimura's real-zero bound, Lean-checked by its author but unrefereed (M. Jaswal, Zenodo 23547020, 21 Oct 2026; found by harvest2); h ≤ 2000 under GRH (Garcia–Grenié–Lee–Molteni, arXiv 2610.14257); counts for odd h < 10⁶ under GRH (Holmin–Jones–Kurlberg–McLeman–Petersen, Exp. Math. 2019); h ≤ 100 (Watkins 2004).

**Method.** Three ranges of d.
- d ≥ e^25.97 ≈ 1.9·10^11: h(−d) > 5000 from our explicit Littlewood bound ([littlewood-half-plane](littlewood-half-plane.md)), which rests on the Lean-checked zero-free half-plane Re s > 0.87338196661320752 for every Dirichlet L-function; certified cells of log d in Arb (threshold.py, loading Littlewood78/cert.py unchanged).
- 56·2^28 = 15,032,385,536 ≤ d ≤ 1.9·10^11: no field has h ≤ 5000. For fundamental −d, distinct primitive ideals of norm a ≤ √(d/4) give distinct reduced forms, so h(−d) ≥ M(d) = Σ_{a ≤ √(d/4)} r(a) with r multiplicative (2 at split, 0 at inert, 1 at ramified primes, 0 at ramified squares). A bit-parallel sieve (isieve.c) certified M(d) > 5000 for all 65,612,855,424 candidates (every d ≡ 3, 7, 11, 15, 4, 8 mod 16), none left over; about 60 min on one core.
- d < 56·2^28: the Mosunov–Jacobson tables (Math. Comp. 85 (2016), unconditional, verified with the Eichler–Selberg trace formula; LMFDB files k = 0..55), extracted with checks (each file's d are exactly the fundamental discriminants of its class and slab; invariant factors multiply to h; gzip CRCs).

**Checks.** The lists for h ≤ 1000 and h ≤ 2000 have the same SHA-256 as 2610.21511's lists; counts for all 2500 odd h ≤ 4999 equal Holmin et al.'s GRH counts; positive-case tests of the sieve (2610.21511's own slab; a window at production size). The orchestrator reran (6 Nov 2026, on copies): two sieve slices ([15,032,385,536, 15,532,385,535] and [1.896·10^11, 1.9·10^11]: all lanes certified), threshold.py (25.97 / 26.11, identical), PARI for every row's largest d (fundamental; qfbclassno = quadclassunit = the stated h), 400 random fundamental d in the sieve range (minimum h 20,020), and the per-h counts and largest d for h ≤ 100 against OEIS A046125 and A038552 (all equal).

**Review.** A reviewer read isieve.c against the bound (correct; an unchecked command-line overflow at N = INT_MAX is irrelevant). A second reviewer (verdict only, argument inline): CORRECT as a deduction from the stated inputs (the ranges overlap; the reduced-form count is a valid lower bound for every fundamental d; non-fundamental candidates are harmless extras; d = 3, 4 lie in the tables); SIGNIFICANT (h ≤ 2000 → h ≤ 5000, including even class numbers).

**Novelty.** arXiv math.NT listings Oct–Nov 2026 (only 2610.21511 classifies), Zenodo ("imaginary quadratic": no classification found then; harvest2 later found Jaswal's conditional h ≤ 3000, Zenodo 23547020), HAL (nothing beyond small h), github.com/iwaokimura/cln-imq-2k (tags to v1.4, last commit 23 Oct 2026, no run beyond 2000), our files (R370 candidate only). New as far as we can see.

**Caveat.** The analytic step is our reviewed paper proof (Littlewood78), not Lean-checked; the tables are trusted as published; the full list (1.4 GB) is not stored (rerun.md recomputes it in about 2 h plus a 26.5 GB streamed download).

**RH relevance:** N11 (the class number problem is the effective shadow of Landau–Siegel zeros; here made effective by the Lean-checked zero-free half-plane).

**Files.** [Riemann/Deep/ImqClass/](../Riemann/Deep/ImqClass/) (REPORT.md, rerun.md, isieve.c, mjx.c, threshold.py, perh_h5000.txt, lists_sha256.txt, logs).

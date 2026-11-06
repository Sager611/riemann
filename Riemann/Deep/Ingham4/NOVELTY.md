# NOVELTY — Ingham4 (6 Nov 2026, 09:40–11:20 UTC)

**Verdict: new (Theorem A).** A prime in (x, x + x^θ] for every real x ≥ 1 was known (explicitly, unconditionally) only with θ = 7/8 (our ShortPrimes, R323; Visser 2610.20058 states the same interval). θ = 0.856 is the first all-x exponent below 7/8 and below 6/7. Theorem B (conditional on an unrefereed density with a corrected constant) is reported as conditional only.

## Previous results found
| statement | source | status |
|---|---|---|
| (x, x + x^{7/8}] for all real x ≥ 1 | our ShortPrimes (7/8 half-plane + Fiori–Jaskari); Visser, arXiv 2610.20058 ("for all real x > 1 there are primes in (x, x + x^{7/8})", found again by the GXL search) | ours / preprint |
| (n^7, (n+1)^7) for all n ≥ 1, i.e. intervals ≈ 7x^{6/7} at x = n^7 only | our Ingham3 (6 Nov 2026) | ours |
| (x, x(1 + 1/(111 log² x))], x ≥ 2898239; (x, x(1 + 1/(5000 log² x))], x ≥ 468991632 | Trudgian 2016; Dusart 2018 | refereed; no power saving |
| (x, x + 3x^{2/3}) for x ≥ exp(exp(33.3)) etc. (cubes) | Dudek 2016, Cully-Hugill 2023 | refereed; huge thresholds |
| [n − n^{0.518}, n] "for sufficiently large n" (abstract: "always") | Lu–Yuan, OSF preprint 10.31219/osf.io/7fyrn (2025), Thm 1 (PDF p. 1: "for sufficiently large n") | unrefereed; asymptotic sieve result, no explicit range |
| x^{0.525} for x ≥ x₀ (x₀ not computed) | Baker–Harman–Pintz 2001 | refereed; not explicit |
| g_n/p_n ≤ 2a(log p_n)^b e^{−c√log p_n} | Visser, arXiv 2211.06469 | no power saving |

## Searches
- **arXiv monthly listings** 2026-10 and 2026-11, math.NT (python3 .scratch/lst.py 'short interval|primes in short|prime in|between consecutive|Ingham|prime gap|gaps between primes|zero.density|explicit' math.NT --zenodo '"short intervals"'): 26 title hits; related: 2610.15989 (short-strip density, Guth–Maynard type, non-explicit), 2610.27830 (RH-conditional short intervals), 2610.34525 (frangible gaps), 2610.10582; none gives an all-x power-saving interval or an explicit density near σ = 0.85.
- **arXiv API** (.scratch/axq.py): 'abs:"zero density" AND abs:explicit' (newest: 2610.15989, 2609.23222, 2405.12545, 2311.05136, 2204.01980, 2101.12263, 2412.02068 = Chourasiya's Carlson estimate), 'abs:"short intervals" AND abs:explicit AND abs:prime' (2206.00433 Cully-Hugill–Dudek mean value, 1902.05065 conditional, 1611.07251 Dudek's thesis, 0810.2113), 'abs:"for all x" AND …', 'abs:"consecutive" AND abs:"sixth powers"' (none), 'ti:"primes between" OR ti:"prime gaps" AND abs:explicit' (2610.20058, 2602.14340 Lee, 2601.15564, 2512.24230, 2501.18048, 2211.06469, 2107.14468, 1401.4233); authors Ramaré–Zuniga Alterman (2603.25961 = the b₄ input), Chourasiya/Simonič (2507.15184 v2, 2412.02068).
- **Explicit densities refereed after KLN** (TeX checked): Johnston–Yang 2204.01980 Table (σ ≥ 0.98 only) and Cully-Hugill–Johnston 2111.10001 Lemma 7.1 (σ ≥ 0.98 only): no gain at σ ≈ 0.85; 2405.12545 (Bellotti) has unsupported constants (REFUTATIONS.md).
- **web_search (arXiv mode)**: "Chourasiya Simonič explicit form of Ingham's zero density estimate journal": only arXiv v2, no journal version.
- **GXL evidence API** (2 claims): all-x exponent below 7/8 → only Visser 2610.20058 (7/8), Dudek's cubes, BHP, Lu–Yuan (asymptotic); Ingham-type explicit density → no explicit Ingham-type result besides 2507.15184.
- **Zenodo** (lst.py, '"short intervals"', 222 hits, newest 15): square-interval and semiprime notes, nothing on all-x prime intervals.
- **HAL** (API): '"short intervals" AND prime AND explicit' (1, unrelated), '"primes in short intervals"' (5, old/asymptotic), '"consecutive sixth powers" OR "consecutive powers"' (0), '"zero density" AND explicit' (unrelated).
- **Our files**: rg over TREE.md, NOGO.md, REFUTATIONS.md, Riemann/Tree, Riemann/Deep (*/NOVELTY.md, ShortPrimes, Ingham3, Harvest*): all-x record 7/8 (ShortPrimes), k = 7 (Ingham3), 2507.15184 only as R405's candidate.

## Theorem B (conditional) against the literature
- Consecutive powers for every n: k ≥ 90 (Cully-Hugill–Johnston 2402.04272), k ≥ 86 (Lee 2602.14340), k ≥ 8 (7/8 half-plane), k = 7 (our Ingham3). Nothing found on k = 6 for every n, conditional or not (arXiv API 'abs:"consecutive" AND abs:"sixth powers"': 0; HAL: 0).
- Thresholds for k = 3: before Ingham3 doubly exponential (Dudek, Cully-Hugill); Ingham3 e^{627.01}; conditional here e^{468.43}.

## What is not new
- The method is I3's (Ingham's smoothed explicit formula + zero-free half-plane + explicit density), with δ = log(1 + x^{θ−1}) instead of k log(1 + 1/n); the small-x part is the standard gap-table argument (as in ShortPrimes).
- The two refinements are routine: a trapezoid weight (convolution of two indicators, as in many smoothing arguments) and KLN's lemma at the newer RH height (done before by Cully-Hugill–Johnston and Johnston–Yang for σ ≥ 0.98).
- Theorem B rests on 2507.15184 (unrefereed) and on Ramaré–Zuniga Alterman's preprint; the correction of b₄ is ours (a 1.6% increase of B₁ near σ = 0.86).

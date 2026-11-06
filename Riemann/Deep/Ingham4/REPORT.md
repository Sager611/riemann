RECORD: a prime in (x, x + x^{0.856}] for every real x ≥ 1 (before: x^{7/8}); conditionally on the unrefereed Chourasiya–Simonič zero density (its constant corrected): x^{0.842}, a prime between n^6 and (n+1)^6 for every n ≥ 1, and n₀(5) = e^{45.584}, n₀(4) = e^{91.18}, n₀(3) = e^{468.43}.

# Ingham4 — primes in (x, x + x^θ] for every x ≥ 1 from the Lean-checked half-plane and explicit zero density (node N10)

## Claim
**Theorem A.** For every real x ≥ 1 the interval (x, x + x^{0.856}] contains a prime. 0.856 is the least multiple of 0.001 certified: at 0.855 the bound reaches 1.152 near log x = 190.6. Ingham3's pipeline unchanged gives 0.857. Theorem A implies Ingham3's k = 7 result, since (n+1)^7 − n^7 > 7n^6 > (n^7)^{0.856}.
Inputs: ζ(s) ≠ 0 for Re s > 0.8733819666132076 (Bebel, Lean); RH up to 3·10^{12} (Platt–Trudgian); Kadiri–Lumley–Ng's explicit zero density (JMAA 2018, Lemma 4.14), evaluated by us at every σ with their RH height replaced by 3·10^{12}, as Cully-Hugill–Johnston (IJNT 2023) did; N(t) ≤ t log t (Hasanalizade–Shen–Wong); prime gaps ≤ 1476 below 4·10^{18} (Oliveira e Silva–Herzog–Pardi) and a sieve to 10^7. No RH above 3·10^{12} and no unrefereed zero count.

**Theorem B (conditional).** If the zero-density bound of Chourasiya–Simonič (arXiv 2507.15184 v2, Cor. 1.1; unrefereed) holds with its B₁ raised for the constant that Ramaré–Zuniga Alterman's preprint (2603.25961) actually proves, then: (i) (x, x + x^{0.842}] contains a prime for every x ≥ 1; (ii) (n^6, (n+1)^6) contains a prime for every n ≥ 1; (iii) (n^k, (n+1)^k) does for n ≥ e^{45.584} (k = 5), e^{91.18} (k = 4), e^{468.43} (k = 3) (Ingham3: e^{57.746}, e^{117.21}, e^{627.01}).

## Before
- Every x: x^{7/8} (our ShortPrimes; Visser 2610.20058 states the same interval). Ingham3's "every n, k = 7" gives intervals ≈ 7x^{6/7}, but only at x = n^7. Log-power intervals: x/(111 log² x) (Trudgian), x/(5000 log² x) (Dusart). Not explicit: x^{0.525} for large x (Baker–Harman–Pintz); Lu–Yuan's OSF preprint claims n^{0.518} for large n.
- R405 estimated 2507.15184 at 0.77× KLN's count (σ = 0.85, T = 10^{13}); that compares with the improved KLN in 2507.15184's own appendix. Against KLN as used here it is 0.085× (9.0·10^9 against 1.06·10^{11}), mostly because its power of log T is about one lower (the Ramaré–Zuniga Alterman input removes a log X from KLN's mean-value lemma).
- k = 6 for every n: open (Ingham3 needed about 6× fewer zeros than KLN allows near Re s = B; this density allows about 12× fewer).

## Method (DERIVATION.md)
1. Ingham3's pipeline with δ = log(1 + x^{θ−1}): the weight sits exactly on [x, x + x^θ]; absolutely convergent explicit formula, ρ ↔ 1 − ρ̄ pairing, KLN in 144 σ-cells, cells in log x (δ decreasing in x), a new tail lemma for log x ≥ 300 (proved non-increasing).
2. Two changes take 0.857 to 0.856: a trapezoid weight (the convolution of indicators of [0, 1/10] and [0, 9/10]) instead of the triangle, which lowers the weight of zeros just above 3·10^{12} (−8%); and KLN's lemma evaluated with RH(3·10^{12}) instead of RH(3.06·10^{10}) (−5%). KLN state RH(H₀) as a hypothesis and set H₀ = 3.06·10^{10} only for their numbers; one of their 20 parameter rows becomes inadmissible and is dropped; our value at σ = 0.98 matches Cully-Hugill–Johnston's within 0.06%.
3. Small x: for x < e^{42} ≈ 1.7·10^{18} the next prime is within the gap after the largest prime ≤ x; the 22 record gaps below 10^7 satisfy g ≤ p^{0.856} (worst: 4 after 7, ratio log 4/log 7 = 0.712), and beyond them every gap below 4·10^{18} is ≤ 1476 ≤ x^{0.856}. Every real x is covered.
4. Theorem B: the Chourasiya–Simonič bound becomes a second bound in every σ-cell (the smaller is used); its exponents 3(1−σ)/(2−σ), (7−5σ)/(2−σ) are concave, so tangent lines give exact σ-integrals.
5. Audit of 2507.15184 (TeX, AUDIT-2507.15184.md): not refereed; its Lemma 4.2 takes b₄ = 0.54 from a private communication; the published result (Ramaré–Zuniga Alterman, Thm 1.2) gives only 0.5501 for 10^{12} ≤ X < 10^{16}. We raise B₁ by (b₄″/0.54)^{(2σ₂−1)/(2−σ₂)} with b₄″ = 0.5501(1 + (1−γ)ε) (+1.6% near σ = 0.86) and reproduce their table from their formulas within 0.1% (cs_repro.py).

## Certificate (one command)
uv run --no-project --with python-flint --with mpmath --with sympy python3 Riemann/Deep/Ingham4/cert.py
(Arb, 128 bits; 632 s, one thread; four parts run one after another in separate processes, peak 104, 92, 104, 164 MB; nothing left running; writes logs/cert.log and logs/cert_result.json.) θ = 0.856: 4932 cells on log x ∈ [42, 300] (max bound 0.9711 at 191.82), tail from 300 (0.0142), small x: CERTIFIED; θ = 0.855 fails. Ingham3's pipeline: θ = 0.857 (0.9319), 0.856 fails. Conditional runs: θ = 0.842 (max 0.8666; 0.841 fails), k = 6 every n (5448 cells on [6 log 100, 300], max 0.7032; n ≤ 100 direct), k = 5, 4, 3 thresholds.

## Checks
- Independent float quadrature (check.py; exact KLN and CS shapes, true kernel envelope, no moment tables): certificate / direct = 1.012–1.207 at 19 points (Ingham3's pipeline: 1.015–1.126 at 8 points); random test of both kernel bounds (6·10⁴ points each, max ratio 0.9999999995); logs/check.log, logs/check_A0.log.
- Fine cells (width 0.004) at θ = 0.855: max 1.152 at log x = 190.61, so 0.856 is not a cell-width artefact.
- δ = 1 in every KLN row (KLN Thm 1.1 as printed): θ = 0.861 certified (0.860 fails).
- Joints: small x below e^{42} (gap records to 10^7, then 1476 up to 4·10^{18}); cells from log x = 42; tail from 300; k = 6 direct check to n = 100 overlaps the first cell.
- Monotonicity: δ(x) decreases and h = x^θ increases in x; cell bounds use the worst end of each cell (asserted log x (1 − B) > 1); the tail majorant is non-increasing (conditions checked in Arb).
- No RH above 3·10^{12} (T₀ = 3·10^{12} in both densities). Theorem A uses no unrefereed zero count.
- cs_repro.py: Chourasiya–Simonič's B₁, B₂, B₃ recomputed from their lemmas: table/ours = 1.0000–1.0009 in all 13 rows used.
- Reviewer subagents on the tail lemma, the small-x step, the B₁ correction, the generalization to real x and the CS-in-cells step: two wording slips and three missing assertions found and fixed (Bertrand for the gap below 4·10^{18}, log x (1 − B) > 1 per cell, CS rows contiguous).

## Novelty
NOVELTY.md: arXiv listings Oct–Nov 2026 (math.NT), arXiv API, web search (arXiv mode), GXL evidence API, Zenodo, HAL, OSF, our files. The only all-x power-saving results found are x^{7/8} (ours; Visser). Lu–Yuan (OSF 2025) is asymptotic ("for sufficiently large n", Thm 1). Refereed post-KLN explicit densities (Johnston–Yang, Cully-Hugill–Johnston) tabulate only σ ≥ 0.98.

## Caveat
- Theorem A rests on Bebel's half-plane (Lean-kernel-checked, unrefereed) and on KLN (refereed; its a₁ = 0.63 is restored by Hiary–Patel–Yang's 0.618), evaluated by us with RH(3·10^{12}) (the same substitution as Cully-Hugill–Johnston's refereed Lemma 7.1; with KLN's own height the exponent is 0.857). KLN's Theorem 1.1 is printed with δ ≥ 1, its Lemma 4.14 and Table 1 with δ > 0; under δ = 1 the exponent is 0.861.
- θ = 0.856 is the limit of this pipeline, not of the method: about 13% fewer zeros at β ∈ [0.83, B] and heights 3·10^{12}–10^{15} would give 0.855.
- Theorem B is conditional: 2507.15184 is unrefereed, its b₄ = 0.54 is unsupported (we use 0.5501 from 2603.25961, itself an unrefereed preprint), and its explicit fourth moment (Mathematica, randomized parameter search) was not rechecked; cs_repro.py checks only the arithmetic from their lemmas to their table.
- The explicit formula is the classical one; no Lean formalization here.
- Side finding (for REFUTATIONS.md): 2507.15184's Table 1 uses b₄ = 0.54 from a private communication, but the published source (2603.25961, Thm 1.2) proves only 0.5501 for 10^{12} ≤ X < 10^{16}; its B₁ must rise by up to 1.6% (σ ≤ 0.875). The theorem is not disproved. Details: AUDIT-2507.15184.md.

RH relevance: N10

## Files
- REPORT.md (this), DERIVATION.md (all steps and constants), NOVELTY.md, AUDIT-2507.15184.md (the TeX audit), CODE-AUDIT.md (rounding of every table and bound, with an independent test)
- cert.py (the certificate), check.py (independent float cross-check; PIPE=A0 for Ingham3's pipeline), cs_repro.py (reproduction of 2507.15184's Table 1), code_audit.py (256-bit table and table-free cell test)
- logs/cert.log, logs/cert_result.json, logs/cert_time.txt, logs/check.log, logs/check_A0.log, logs/cs_repro.log, logs/code_audit.log

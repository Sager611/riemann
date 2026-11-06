RECORD: J^{d,n} (Jensen polynomials of Riemann's Xi) is hyperbolic whenever n^3 log^2(n+2) >= 7000 d^5 (every d >= 1, using the all-n range d <= 9.0e24); with that range this is the first explicit hyperbolic region unbounded in d.

# Explicit constant in Holland's hyperbolicity wedge (node N6)

## Theorem

Let `D_all >= 9.0·10^24` be a degree up to which J^{d,n} is hyperbolic for all n (Chasse's sector lemma with RH to
height T = 3,000,175,332,800 gives T^2 ≈ 9.0·10^24; ours in TREE N6 is d0 = 16,140,857,488,599,731,839,205,676).

1. For all integers d > D_all and n >= 0 with **n^3 log^2(n+2) >= 7000 d^5**, J^{d,n} has d distinct negative zeros.
2. Hence J^{d,n} is hyperbolic for every (d, n) with **d <= D_all or n^3 log^2(n+2) >= 7000 d^5**; in particular
   "n^3 log^2(n+2) >= 7000 d^5 => hyperbolic" holds for every d >= 1.
3. Self-contained form of Holland's Theorem 1.1 (no all-n range used): K = 2.43·10^114 works for every d >= 1
   (certificate with L* = 80, certify_self.log; there N* = 40(pi e^80 - 9/4) ≈ 7.0e36 and K = N*^3 log^2(N*+2)).

Before: Holland (arXiv 2608.08682, 9 Aug 2026, TeX l. 81-89) proves the wedge with an absolute K that is never
computed (l. 1814-1823, 2042-2043). Every other large-n result is ineffective (below). The only explicit region was the
all-n range, which is bounded in d.

## Method

Holland's chain is kept: Laguerre -> Jacobi -> finite-free (Jacobi x Jacobi) comparison model matching R_0..R_4,
Gershgorin root localization, fifth-order residual multiplier, perturbed-Jacobi derivative ratios, stability
proposition. Its analytic input, the complex form of the GORTTW saddle asymptotic (Holland l. 508-520), is a
one-sentence claim without proof or constant in GORTTW (1910.01227, TeX l. 514-518). We replace it by our own explicit
complex saddle bound on the disk |z - n| <= n/10 (relative error <= 5.33e-20 for n >= N* = 45(pi e^90 - 9/4) ≈ 1.725e41).
Every "sufficiently large n" becomes n >= N*, every O(.) an Arb ball enclosure over 1000 cells of eps = 1/L,
L = 2u_0(n) ≈ log n - log log n. Sharpenings: Brouwer instead of the 1/30 contraction; Gershgorin constant 2 instead
of 8; a signed fifth-derivative bound (12.769 instead of about 84 by the triangle inequality); Cauchy radius 1.4r instead of
2r in the stability proposition; a + b + theta < 1 instead of each <= 1/4. Full derivation: DERIVATION.md.

## Constants (Holland TeX line -> explicit value)

| Holland | TeX line | here |
|---|---|---|
| K (main theorem) | 82 | 7000 (d > D_all >= 9e24); 2.43e114 for all d |
| sector theta, x_0 (Prop. 4.1) | 473 | disk |z-n| <= n/10, n >= N* ≈ 1.725e41 |
| O(|z|^{-1+eps0}) (complex GORTTW (3.2)) | 519 | |delta| <= 5.33e-20; log-error G: sup|G - G(n)| <= 0.053287 |
| eta, C (h'''', h^(5) bounds) | 578 | eta = 1/10; n^4 L |Psi^(5)| <= 14.968 (|Q5|) |
| "sufficiently large n" | 624, 653, 683, 730, 926, 1015, 1646 | n >= N* |
| q0 > 1, D_n > 0, Q_n > 0 | 624, 653, 926 | n^2 L D_n >= 1.9129, n^3 L Q_n >= 1.9116 |
| A_J ~ nL | 742 | A_J/(nL) <= 1.0226 |
| K_r >= 32 | 1014 | 32 (B, D >= 32d checked) |
| eps_n <= 1/30 | 1175 | Brouwer box [1.10,1.30]x[0.25,0.42]x[0.015,0.07]; c >= 0.04012 |
| K_4 and the finite range of n | 1032, 1215 | absorbed: d > D_all and the wedge force n > N* |
| C_J = 8 | 1455 | 2.00000002 |
| C_loc | 1379 | 2.0874 |
| K_6 | 1531, 1645 | 2.6 |
| C_* (fifth derivative) | 1588-1621 | n^4 L sup|E_F^(5)| <= 12.769 |
| C_* (V_m) | 1723-1734 | C_V <= 1.00000001 |
| a_m, b_m, theta <= 1/4 | 1793 | a + b + theta <= 0.97084 |
| C_def, K_5 | 1524-1527 | Cfin <= 70.600; Kres = (Cfin/thetaE)^2 (1+lamL)^2 <= 6416.59 (thetaE >= 0.929840, 1+lamL <= 1.05502) |
| eps < 16 (stability) | 233 | e <= e^{thetaE} - 1 < 1.4^4 * 0.4 = 1.53664 |
| n_0 | 1815 | N* >= 1.7253e41; KN = N*^3 log^2(N*+2)/D_all^5 <= 784.13 |

## Certificate (one command, 4 seconds)

```
cd Riemann/Deep/HollandK && uv run --no-project --with python-flint --with sympy python3 certify.py 90 2.6 1000 1.4
```
It prints "K = max(Kres, KN(9e24)) = 6416.58.. <= 7000: PASS" and a block of documented constants, every one rounded
outward (upper bounds up, lower bounds down; certify.log):
* A: beta <= 0.0533474, c3 <= 7.63354, phi2 <= 1.01113, cos(2beta) >= 0.994313, rel_core <= 3.36e-20, rel_am <= 1.97e-20,
  the other four terms < 1e-300, deltaA <= 5.33e-20, b_rk <= -1.2224, gnet >= 0.61746, etaG <= 0.053287; |Q5| <= 14.968.
* cells: tQ >= 1.9116, n^2 L D_n >= 1.9129, aJ <= 1.0226, xc >= 0.04012, CJ1 <= 2.00000002, Cloc <= 2.0874,
  C_V <= 1.00000001, a + b + theta <= 0.97084, ce5 <= 12.769, Cfin <= 70.600.
* final: thetaE >= 0.929840, 1 + lamL <= 1.05502, Kres <= 6416.59, KN(d0) <= 42.27, KN(9e24) <= 784.13, K <= 6416.59.
The same command with 80 instead of 90 also passes (certify_self.log): beta <= 0.0534309, c3 <= 7.63382, deltaA <= 7.91e-18,
etaG <= 0.053366, Cloc <= 2.0930, a + b + theta <= 0.97568, Cfin <= 72.118, Kres <= 6764.23 <= 7000, KSELF <= 2.430e114.

## Checks (sanity, not part of the proof; sanity.log)

* check_identities.py: the three exact identities of Holland's Lemma 9.1 and the 3F2 equation hold (symbolic, exact).
* check_jacobi.py: the Jacobi-matrix entries and the diagonal-difference identity reproduce the roots (error < 1e-36).
* sanity.py: kernel normalization (moment ratios 4, 16, 64: a pure 4^n rescaling, irrelevant for h'', h''', ...); the
  saddle main term matches quadrature to 2e-3 at x0 = 2000, also for complex z with |Im z| = x0/10.
* sanity2.py b: the pipeline at n ≈ 1.7e41 (L = 90) from the main term at 600 digits: n^2 L D_n = 1.91305,
  tQ = 1.91282, aJ = 1.02247, c = 0.04190, inside the certified enclosures (>= 1.91298, >= 1.91168, <= 1.02256, >= 0.0401).
* sanity3.py: at n = 1e43 and d = 6.75e25 on the wedge boundary, the true max |E_F| on the boundary of Omega is 0.784;
  the certified bound there is <= 0.8887 (Cfin <= 70.600); the requirement is <= thetaE (>= 0.929840). The bound is
  within 14% of the truth.
* sanity2.py a: J^{d,n} for d = 2..30 and n in {0,1,2,3,5,10,20,50,100,200,300} (319 polynomials, gamma by
  50-digit quadrature of the Polya integral): all hyperbolic, real roots certified by Arb root isolation (python-flint
  complex_roots returns certified real roots with exactly zero imaginary part).

## Review

Reviewer subagent, three rounds on the argument pasted inline. Round 1: all replacements valid (lambda = 1.4 variant,
Brouwer, contour/concavity, Cauchy for G, wedge logic); gap = holomorphy/branches/uniformity not stated. Round 2: those
closed except the enclosure of the variable s in Q5 and three wording points. Round 3: (b)-(e) closed; the last point
(t-ball centred at twice the eps-cell) is "merely wording" — the certificate does centre it at 2 eps. No mathematical
error was found.

## Repair (orchestrator review)

The orchestrator's reviewer found two constants in DERIVATION.md rounded the wrong way. certify.py never used the
printed values: it computes c3, beta and every other constant in Arb and takes ball endpoints (up/lo). So the earlier
K, deltaA and Kres were valid; only the text was wrong. certify.py now prints every quoted constant itself, rounded
outward by exact rational arithmetic, and both certificates were rerun (K still <= 7000).

| constant | old text | new (outward) | where used |
|---|---|---|---|
| c3 (cubic Taylor constant, A2) | 7.63 | <= 7.63354 (L* = 90), 7.63382 (L* = 80); exceeds (4pi/3)e^{0.6} = 7.63247.. | C3 s1 <= kappa0/4, rel_core, deltaA |
| beta (saddle shift, A1) | <= 0.0533 | <= 0.0533474 (L* = 90), 0.0534309 (L* = 80) | cos(2beta), alo, phi2, c3, all of Part A, etaG |
| deltaA (saddle error) | 5.4e-20 / 5.3e-20 | <= 5.33e-20 (L* = 90), 7.91e-18 (L* = 80) | etaG, G-terms |
| etaG (log error) | <= 0.0533 | <= 0.053287 (L* = 90), 0.053366 (L* = 80) | Cauchy bounds g_k, G^{(5)} term |
| abs(Q5) | 14.56 (DERIVATION) | <= 14.968 (L* = 90), 15.526 (L* = 80) | h^(5) on Omega |
| n^2 L D_n, n^3 L Q_n (lower) | >= 1.913, >= 1.912 / 1.9130, 1.9117 | >= 1.9129, >= 1.9116 | defect signs, T > q0, matching target |
| a + b + theta | <= 0.971 / 0.9708 | <= 0.97084 (L* = 90), 0.97568 (L* = 80) | derivative-ratio contradiction |
| Cloc | <= 2.09 / 2.088 | <= 2.0874 (L* = 90), 2.0930 (L* = 80) | a_m |
| CJ1, C_V | 2 + 1e-8, 1 + 1e-8 | <= 2.00000002, <= 1.00000001 | Cloc, theta |
| ce5, Cfin | 13.3, 73.2 (DERIVATION) | <= 12.769, <= 70.600 (L* = 90); 12.981, 72.118 (L* = 80) | Kres |
| thetaE, 1 + lamL | 0.9298, 1.0550 | >= 0.929840, <= 1.05502 | Kres |
| lambda^4 (lambda - 1) | 1.5366 | 1.53664 (exact) | stability threshold |
| Kres, KN(9e24) | 6416.6, 784.1 / 784.2 | <= 6416.59, <= 784.13 | K = max <= 7000 |
| Kres (L* = 80) | 6764.2 | <= 6764.23 | self-contained run |
| KSELF | 2.43e114 | <= 2.430e114 | Theorem item 3 |
| sanity3 bound / requirement | 0.921 / <= 0.930 | <= 0.8887 / <= 0.92984 | check only |

Recertified: K <= 6416.59 <= 7000 (L* = 90), deltaA <= 5.33e-20, Kres <= 6416.59; KSELF <= 2.430e114 (L* = 80,
Kres <= 6764.23). The wedge implication was confirmed by the reviewer.

## Novelty

* arXiv listings 2026-10/11 (math.NT, CA, CV; .scratch/lst.py 'Jensen|hyperbolic|Turán|Laguerre'): nothing on explicit
  Jensen regions (2610.38870 is a GRH criterion).
* Zenodo '"Jensen polynomials"': Savva 22293642 (3 Sep 2026, unrefereed, AI-assisted) claims a stronger two-thirds wedge
  n^2 log(n+2) >= K d^3 with K "effective in the usual proof-theoretic sense" but explicitly leaves two analytic constants
  symbolic and prints no K (PDF text l. 186-193, 2000-2013). Chen 22951677 and Chung 23987495 are RH claims; Schieber
  22936342 is conditional on hyperbolicity; Dumas 22897068 transfers published bounds.
* HAL: Planat hal-05675741 (1 Jul 2026) states "there exists C0 > 0 ... n >= C0 d^4" with C0 ≈ 0.020 obtained
  "numerically" from fitted constants (K_tau^eff fitted on n <= 23; a bound "verified numerically for n = 10..50"); not
  an explicit rigorous constant.
* Earlier: GORTTW Thm 1.1 (c e^{d/2}, c from a compactness argument, TeX l. 306, 584-587), GORZ and Kim-Lee
  (ineffective, R313), Holland (non-explicit). arXiv API and arXiv web search for citing papers: HTTP 429 (throttled).
* Verdict: the explicit K (7000, and 2.43e114 without an all-n range) is new as far as found.

RH relevance: N6 (RH iff all J^{d,n} hyperbolic; first explicit region unbounded in d)

## Files

certify.py (the certificate), certify.log, certify_self.log, DERIVATION.md, check_identities.py, check_jacobi.py, sanity.py,
sanity2.py, sanity3.py, sanity.log.

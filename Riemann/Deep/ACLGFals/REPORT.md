HOLDS: nothing fails on the computable range: Galerkin spectra lambda^2 = 2...20 (eps(j)/delta_4j 1.67-11.96, odd 1.67-8.39, interlacing, counts lambda^2-1 / lambda^2-2), certified eps(1) >= 3.63 delta_4 and positivity for lambda^2 = 3...31, 45, 47, Gram constant >= 0.151, (H2ii) for S = {inf, 2}: Tr(M_chiR N^(n0)) ~ 0.006 R^0.8/n0 (for S = {inf, 2, 3} the decay starts only past n0 ~ 50); smallest certified margin eps(2) >= 4.9e-4 delta_8 (lambda^2 = 20).

# ACLGFals: falsification of the hypotheses of Andrews-Connes-La Gatta (arXiv 2611.13129) at small lambda (aclgfals, 2026-11-11)

Statements with definitions and line numbers: [STATEMENTS.md](STATEMENTS.md). Models: [DERIVATION.md](DERIVATION.md). Source: e-print
TeX `prolate-weil-fine15.tex` (1527 lines, v1, 9 Nov 2026). Notation: c = lambda^2; delta_n(lambda) = 1 - |chi_n(lambda)|, the exact rung.

## 1. Statements
| # | statement | TeX lines | class | test | result |
|---|---|---|---|---|---|
| 1 | Hypothesis LA: c delta_4j <= eps(j) <= C delta_4j log lambda; odd sector with delta_4j+2 | 588-599 | (iii) as stated; (ii) operative | ratios, log-gaps, interlacing, counts, lambda^2 = 2...20; certified j = 1 bounds up to 47 | holds, no trend (2.1, 2.2) |
| 2 | Positivity QW_lambda >= 0 | 898, 904-915, 1035, 1357 | (i) per lambda; RH for all lambda | sign of every Galerkin eigenvalue; certificates | holds: certified for lambda^2 <= 47 (Chekmachev, audited; c = 45, 47 on his unpublished code); all Galerkin eigenvalues > 0 |
| 3 | Gap (Lemma gap): eigenvalues other than eps >= g = c delta_6 (full) or c delta_8 (even) | 904-920 | (ii) | certified lower bounds against delta_6, delta_8 | holds: certified c >= 0.0707 (full), >= 4.9e-4 (even), lambda^2 <= 20 |
| 4 | (H2'): sum_{n>n0} kappa_n ||theta(f) a_n||^2 <= eps_t(n0) Lambda(f), eps_t(n0) -> 0 | 1221-1240 | (iii) | none possible as stated; its mechanism | kappa_n sqrt(n) -> about 0.9 (order 1/2, as Lemma H2 says) |
| 5 | Localisation (H2ii): Tr(M_chiR N_t^(n0)) -> 0 polynomially in R | 1242-1262 | (ii) | Nystrom model of the semilocal pair, S = {inf, 2} and {inf, 2, 3}, t = 1 | S = {inf, 2}: holds on 20 <= n0 <= 100, R = 5, 20, about 0.006 R^0.8/n0; S = {inf, 2, 3}: low-frequency mass extends to n ~ 50, rate not reached (2.4) |
| 6 | (R2), Gram part: inf ||E f||^2/||f||^2 over H_K(lambda) bounded below up to K ~ lambda^2 | 1364-1373, 1395 | (ii) | g_K^Herm, K <= 40; g_K(lambda), K <= lambda^2 + 2 at lambda^2 = 5, 12, 20, 30, 50 | holds: >= 0.1510 (2.3) |

## 2. Computations
### 2.1 Galerkin spectra, lambda^2 = 2...20 (code/spec_small.py, run_spec.sh, analyze_spec.py; logs/la_tables.md)
Matrix: weilsimple's closed-form Arb matrix (same space E_N as ACLG App. D), N = 120 and 240, 192-512 bits. Galerkin eigenvalues are
min-max UPPER bounds: they test the upper halves and the shape; the lower halves are tested by the certified bounds of 2.2.
Validation: lambda^2 = 5, N = 240: -log eps(1) = 39.219, as ACLG Table conv; our exact rungs (code/prolate.py) give ACLG's D_1 = 41.185 and
-log delta_6, delta_8, delta_10 = 33.833, 27.235, 21.261 at lambda^2 = 5 to every printed digit.
| c | eps(1)/delta_4 | eps(1)/(delta_4 log lambda) | eps(j)/delta_4j, all minuscule j (#) | eps_odd(j)/delta_4j+2 (#) | log-gap shortfall even / odd | interlacing |
|---|---|---|---|---|---|---|
| 2 | 4.51 | 13.00 | 4.51 (1) | - (0) | 1.128 / 0.675 | yes |
| 3 | 6.36 | 11.58 | 2.29-6.36 (2) | 4.09 (1) | 1.022 / 0.443 | yes |
| 5 | 7.15 | 8.88 | 1.77-7.15 (4) | 3.05-5.05 (3) | 0.776 / 0.348 | yes |
| 8 | 8.44 | 8.12 | 2.07-8.44 (7) | 2.09-5.41 (6) | 0.672 / 0.446 | yes |
| 12 | 8.85 | 7.13 | 1.82-8.85 (11) | 1.88-6.71 (10) | 0.596 / 0.277 | yes |
| 16 | 10.41 | 7.51 | 1.67-10.41 (15) | 2.16-7.52 (14) | 0.573 / 0.325 | yes |
| 20 | 11.96 | 7.98 | 1.74-11.96 (19) | 1.67-8.39 (18) | 0.565 / 0.354 | yes |
(c = 4, 6, 10 in logs/la_tables.md, same pattern.) Below 1e-2 there are exactly lambda^2 - 1 even and lambda^2 - 2 odd eigenvalues for
every c = 2...20, one fewer than U_lambda in each sector. At the threshold 1e-10 the difference is 2 at c = 8 (even: 4 vs 6;
eps(5) = 1.44e-10, delta_20 = 5.5e-11) and c = 20 (odd: 15 vs 17): threshold effects, as ACLG's Fact count item 3 says; c = 8, 20 are
outside their table. Fact logbound's constant 8.8 is exceeded at c = 2, 3 (13.0, 11.6), outside its stated range 5 <= lambda^2 <= 300.
N-convergence: log(eps1(N = 120)/eps1(N = 240)) = 0.0003 (c = 2) to 0.62 (c = 20); with increments shrinking by 0.3-0.4 per doubling
(ACLG), the ratios at c = 16, 20 are overestimates by about 15-40%.

### 2.2 Certified finite-lambda instances (code/rungs_cert.py, cert_la.py; logs/cert_la.json, rungs_cert.json)
Inputs: Chekmachev's eps_c >= 0.97 Theta_c, 3 <= c <= 31 and c = 45, 47 (Zenodo 23848597, a lower bound of QW_lambda on all of L^2,
audited sound in ../ChekAudit); weilsimple's certified U >= eps(1), L_odd <= eps_odd(1), L_2 <= eps(2) (../WeilSimple); our exact rungs.
| c | 3 | 5 | 10 | 16 | 20 | 25 | 31 | 45 | 47 |
|---|---|---|---|---|---|---|---|---|---|
| eps(1)/delta_4 >= (0.97 Theta_c) | 4.15 | 3.63 | 3.98 | 4.60 | 4.99 | 5.46 | 5.97 | 7.03 | 7.17 |
| eps(1)/delta_4 <= (U) | 6.86 | 9.41 | 11.14 | 13.69 | 14.76 | - | - | - | - |
| eps_odd(1)/delta_6 >= (L_odd) | 0.272 | 0.196 | 0.106 | 0.085 | 0.071 | - | - | - | - |
| eps(2)/delta_8 >= (L_2) | 0.0148 | 0.0062 | 0.0015 | 0.0007 | 0.0005 | - | - | - | - |
The minimum of the certified lower ratio over every integer c in [3, 31] and c = 45, 47 is 3.63 (c = 5). Chekmachev's odd-sector
constants give eps_odd(1) >= 2.81 delta_6 (c = 5) and 3.89 delta_6 (c = 16). L_odd and L_2 fall with c because our certificate method
loses accuracy, not the form (Galerkin eps_odd(1)/delta_6 = 2.3-8.4, eps(2)/delta_8 = 1.5-6.8). So at every certified window the lower
half of LA (j = 1) and the gap hypothesis hold with explicit constants. LA itself (lambda -> infinity) is not touched.
Side remark (logs/theta_vs_aclg.md): Chekmachev's empirical scale is Theta_c = (1.03 + o(1)) lambda delta_4(lambda) (leading terms of
Phi(log lambda) and Fuchs). So eps_c >= 0.97 Theta_c cannot hold for all large c: it contradicts eps_lambda <= C delta_4 (log lambda)^3
(Kastoris, Connes-Kastoris). ACLG's min-max upper bound at lambda^2 = 250 (16.08 delta_4, N = 5000, still falling) is only 1.2% above
0.97 Theta_c = 15.90 delta_4. His theorem claims the bound only for c <= 47, so nothing he proved is affected.

### 2.3 Gram constant, input of (R2) (code/gram_herm.py, gram_lambda.py; logs/gram_*.jsonl)
Hermite limit, from ||E h||^2 = (1/2pi) int |zeta(1/2+it)|^2 |Mh(1/2+it)|^2 dt (Mellin transforms by FFT; quadrature check
||G_2 - I|| = 2e-15): g_1 = 0.2192, g_5 = 0.1878, g_10 = 0.1610, g_20 = 0.1518, g_30 = 0.15105, g_40 = 0.15096 (the minimiser's
Mellin mass sits at t ~ 2.5). Finite lambda (prolate h_4k,lambda, E restricted to the window, every K through the plunge):
min_K g_K(lambda) = 0.2059 (c = 5, K <= 7), 0.1637 (c = 12), 0.1538 (c = 20), 0.1514 (c = 30), 0.1510 (c = 50, K <= 52); g_1(lambda) =
0.2227, 0.2204, 0.2199, 0.2197, 0.2195 tends to g_1^Herm = 0.2192, which cross-checks the two codes. So on this range the Gram constant
stays >= 0.151 through the whole ladder, and (R2) loses nothing from it. g_K^Herm is non-increasing in K (nested spaces); it can tend
to 0 only by concentrating the Mellin mass |Gamma_R P|^2 at zeros of zeta, and the weight |Gamma_R(1/2+it)|^2 ~ e^{-pi |t|/2} is the
Carleman-borderline case for polynomial approximation, where such concentration is at best logarithmic in K (heuristic, not proved).
Observed: -9e-5 from K = 30 to 40.

### 2.4 Semilocal pair and (H2ii) (code/hankel.py, semi_pair.py; logs/hk_*.json, sp_*.json)
Two Galerkin models (DERIVATION.md sec. 1, 3): (a) Laguerre: Gamma_t is the Hankel matrix (phi_{j+l+1}) of the symbol on the circle;
(b) Nystrom: Gamma_1 is unitarily equivalent to (A G)(w) = 2 int_0^1 C_S(uw) G(u) du on L^2(0,1) with
C_S(z) = sum_n gamma(n) cos(2 pi n z); for S = {inf, 2}, C_S = (1/2) sum_{k>=0} cos(2 pi 2^k z) - (1/2) cos(pi z), cut at n <= nmax.
- Validation, S = {inf}: kappa_n = 1 - delta_n(1) (0.999971376, 0.979484735, 0.524085896, 0.058976589) to 2e-8 in (a) and 1e-9 in (b);
  <e_n, f_n> = kappa_n to 6e-6 (n <= 4); the low-frequency masses of (a) and (b) agree to 2e-4.
- S = {inf, 2}: (a) gives kappa_0 = 1 - 2.6e-6, kappa_1 = 0.9927, kappa_2 = 0.6195 (M = 1500) but resolves only about 30 pairs. In (b),
  nmax = 64...1024: kappa_n sqrt(n) at fixed n/nmax is exactly self-similar (0.442 at (n, nmax) = (50, 64), (100, 128), (200, 256),
  (400, 512)) and at n = 50 it rises 0.44, 0.60, 0.71, 0.78, 0.82: kappa_n ~ C n^{-1/2} with C about 0.9, the order 1/2 of Lemma H2;
  Sigma kappa_n^2 grows by about 0.5 per doubling (4.01, 4.54, 5.05, 5.54, 5.92; the last capped at 800 vectors), as
  Sigma kappa^2 = infinity requires.
- The localisation quantity, as sums of kappa_n m_n(R) over windows of n (even boundaries: couples kappa_{2k+1} ~ kappa_{2k+2}):
  | window | (10,20] | (20,50] | (50,100] | (100,200] | (200,400] |
  |---|---|---|---|---|---|
  | R = 5, nmax = 256 / 512 / 1024 | .00238 / .00244 / .00248 | .00064 / .00060 / .00058 | .00037 / .00025 / .00022 | .0373 / .00028 / .00013 | .081 / .030 / .015 |
  | R = 20, nmax = 256 / 512 / 1024 | .0232 / .0235 / .0237 | .0024 / .0022 / .0022 | .0012 / .0008 / .0007 | .249 / .0007 / .0004 | .302 / .200 / .029 |
  A window is large while the cut is near it (vectors at n ~ 0.5-2 nmax with kappa 0.01-0.05 and low-frequency mass up to 0.95, an
  artifact of the cut), then drops by orders of magnitude and settles once nmax >= 5 b. Converged windows fall by a factor 1.7-3 per
  doubling of n. With the last window as tail estimate, T(n0, R) = Tr(M_chiR N^(n0)) is (0.021-0.026)/n0 (R = 5) and 0.07/n0 (R = 20)
  for n0 = 20, 50, 100: a polynomial rate in n0 with R-dependence about R^0.8 (0.7-0.9), which is what (H2ii) asks for; the values
  still fall with nmax, so these are overestimates. The first ten pairs carry 1.46 (R = 5) and 3.87 (R = 20), bounded as Lemma
  localization proves. Possible tool (novelty against CCqi/CCsel unchecked): the semilocal pair at t = 1 is the symmetric kernel
  2 C_S(uw) on L^2(0,1), C_S the lacunary cosine series of the local factors.
- Two primes, S = {inf, 2, 3} (window 3 <= lambda^2 < 4; 37-52 Laurent terms, nmax = 128, 256, 512, peak 198 MB): kappa_n is larger
  (kappa_n sqrt(n) = 1.24 at n = 50, nmax = 512, still rising) and the low-frequency content reaches much higher index. Window sums at
  nmax = 128 / 256 / 512: R = 20: (10,20] 1.18 / 1.21 / 1.25, (20,50] 0.51 / 0.96 / 0.98 (one prime: 0.0022), (50,100] 0.54 / 0.18 /
  0.058 (edge passing); R = 5: (10,20] 0.05 / 0.14 / 0.16, (20,50] 0.009 / 0.018 / 0.019. So the index n0 beyond which T(n0, R) is small
  grows with |S| (here from about 10 to about 50 at R = 20), as ACLG's caveat on the counting function of the angles (l. 1259) says;
  whether the growth stays polynomial in |S| ~ lambda^2/log lambda is exactly what (H2ii) with Theorem condpos needs, and it is not
  measurable at these sizes.
- Scope: S = {inf, 2} and {inf, 2, 3} (windows lambda^2 < 4), scale t = 1 only, n <= 200, a cut Laurent series; Theorem condpos needs
  S with all p <= lambda^2 and scales sigma > 1. Not a proof of anything; no sign of failure.

## 3. Novelty and prior tests
- arXiv listings 2026-10/11 (math.NT, CA, CV, SP, PR, FA, CO, OA, math-ph, NA; titles matching Connes|Consani|Moscovici|La Gatta|
  Weil quadratic|prolate|ground state|Weil positivity|uncertainty operator|Sonin|semilocal): only 2611.13129 itself (v1). arXiv search
  (three queries): nothing else.
- Zenodo: Spark 24126058 (10 Nov; at least m eigenvalues of the semilocal Weil operator in [-delta_m, delta_m]: an upper-bound statement),
  Wiktorowicz 23983773 (under RH and |zeta'(rho)|^{-1} <= |gamma|^{theta+eps}: eps_lambda >= e^{-(4 pi + eps) lambda^2}, the exponent of
  LA's lower half, not its constant), Chekmachev 23848597 (certificates, compared with Theta_c, not with delta_4), Botha 23852617,
  Kramarenko-Byrd 23007100, Mori 24047155, Connes-Kastoris 22864295 (upper bounds). HAL: nothing relevant. Our files: R419, R432 and
  Harvest6 only list the hypotheses.
- No numerical test of (H2'), (H2ii), (R2) or the Gram constant, and no certified comparison with delta_4, was found. Gate: the
  certified table of 2.2 is a routine corollary of audited certificates (COROLLARIES); the Gram values and the (H2ii) profile are
  records of a computation (RECORDS, N13); the remark on Theta_c is a routine deduction (COROLLARIES). No DISCOVERIES entry.

## 4. Verdict in plain words
Nothing in Andrews-Connes-La Gatta's hypotheses breaks where it can be checked. Up to lambda^2 = 20 the tiny eigenvalues of the Weil
form follow the prolate ladder within factors 1.7 to 12, interlace, and come in the predicted numbers; at the certified windows up
to lambda^2 = 47 the lowest eigenvalue is provably at least 3.6 times the first rung and the form is provably positive. The Gram
constant that their open step (R2) needs stays above 0.15 through the whole ladder up to lambda^2 = 50. For their main open step,
the localisation (H2ii), a model of the semilocal operator with the prime 2 shows the predicted behaviour: the near-self-dual
directions of high index have little low-frequency content, and the quantity falls like 1/n0. With the primes 2 and 3 the same
quantity stays large up to about five times higher index, so the onset of the decay moves with the number of primes; how fast it
moves is the open point. These are checks on small cases; the asymptotic statements (LA for lambda -> infinity, (H2') uniformly
in f) cannot be falsified by finite computation.
Trend: eps(1)/delta_4 rises slowly (4.5 at lambda^2 = 2, 12 at 20; 19 at 300 in ACLG) while eps(1)/(delta_4 log lambda) stays at 7-9
for lambda^2 >= 4;
the log-gap shortfalls settle near 0.57 (even) and 0.35 (odd). The weakest spot is the smallest ratio over all minuscule j, near the
plunge: 1.67-2.3 here (lambda^2 = 5...20), 1.40, 1.36, 1.23 at lambda^2 = 50, 100, 150 (ACLG), and these are upper bounds. LA with fixed J does not need
it, but a lower bound uniform in the index (what (R2) and Prop. reduction use near the transition block) would fail first if the
drift continues.

Process note: five runs broke the 300 MB cap: semi_pair.py first version (879 MB, 4 s), nmax = 256 by SVD (459 MB, 9 s),
nmax = 512 (313 MB, 34 s), a first nmax = 1024 attempt (killed at 517 MB) and the kept nmax = 1024 run (427 MB peak per
/usr/bin/time, missed by the 2 s RSS polling; 59 s). The code was then fixed (row-block kernel, Fortran-order in-place eigensolver,
real matrix products); its nmax = 256 rerun peaks at 160 MB. All other runs stayed below 270 MB.

RH relevance: N13 (hypotheses of a conditional RH theorem)

## Files
- STATEMENTS.md (the six statements), DERIVATION.md (models), this REPORT.md.
- code/: prolate.py (exact rungs), spec_small.py, run_spec.sh, analyze_spec.py (Galerkin spectra), rungs_cert.py, cert_la.py
  (certified ratios), theta_vs_aclg.py, gram_herm.py, gram_lambda.py (Gram constant), hankel.py, coeff_probe.py (Laguerre/Hankel
  model), semi_pair.py (Nystrom model; usage: semi_pair.py primes nmax R1,R2 [npw] [taper]; env SP_EIGH, SP_F32, SP_CONSISTENT).
- logs/: la_tables.md, la.json, spec.jsonl, cert_la.json, rungs_cert.json, theta_vs_aclg.md, gram_herm.jsonl, gram_lambda_c*.jsonl,
  hk_*.json, sp_*.json, *.time (wall time, peak memory).

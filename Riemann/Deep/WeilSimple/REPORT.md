RESULT: For every lambda in [sqrt 2, sqrt 20] (window c = lambda^2 in [2, 20]) the smallest eigenvalue of the Weil form QW_lambda of Connes-Consani-Moscovici (arXiv 2511.22755) is simple with an even eigenvector (computer-assisted, Arb). At 52 explicit windows (c = 2, 3, ..., 20 and the 38 points of a chain) the hypothesis of their Theorem 5.10 holds for every truncation N >= N_0 (30 to 160), so its conclusion does too; on the integer windows the upper bound U of the ground eigenvalue is 7.7 (c = 2) to 150 (c = 20) times below the lower bounds L of the next eigenvalues, along the chain at least 2 times. Known before: the single window lambda = e^0.8 (Zhu, arXiv 2608.24827).

# WeilSimple: the Weil form QW_lambda has a simple, even ground state for all lambda in [sqrt 2, sqrt 20] (weilsimple, 2026-11-05)

RH relevance: N13 (spectral realization: CCM's self-adjoint approximants) and N1 (Weil positivity).

## 1. Statement
Notation: c = lambda^2, L = ln c = 2 log lambda; QW_lambda on L^2([lambda^-1, lambda], d*u) is, via u = e^x, the Weil form Q of
../WeilZeta on L^2(-L/2, L/2) (DERIVATION.md sec. 1, term by term); "even" = invariant under u -> 1/u (x -> -x).

Theorem (computer-assisted: python-flint 0.9.0 Arb balls, plus the proved-rounding node scan of ../WeilZeta/RIGOR.md).
(i) For every lambda in [sqrt 2, sqrt 20] (c in [2, 20]), the self-adjoint operator A_lambda with QW_lambda(f, f) = <A_lambda f, f> (CCM sec. 3.3) has
    a simple smallest eigenvalue mu_lambda, and its eigenvector xi_lambda is even. More precisely spec(A_lambda) meets
    (-infinity, L(c)) in the single simple eigenvalue mu_lambda <= U(c), with U, L as in the tables (for c between two chain points,
    U of the left point and L of the right point). This is the first "missing step" of CCM sec. 8 (TeX l. 1390-1393) for these lambda.
    (It does not use CCM Thm 3.6 on discreteness: the three inequalities alone force a one-dimensional even spectral subspace below L.)
(ii) At each window of the tables (the 19 integers c = 2, ..., 20 and the 38 points of the chain, rational c; 52 distinct), for every N >= N_0(c)
    the truncation QW_lambda^N (restriction to E_N = span{V_n : |n| <= N}) is even-simple (CCM Def. 5.3). Hence CCM Thm 5.10
    holds for (lambda, N): with eps_N in [mu_lambda, U(c)] its smallest eigenvalue and xi its even eigenvector normalized by
    delta_N(xi) = 1, (a) D_{lambda,N} = D_log^(lambda) - |D_log^(lambda) xi><delta_N| is self-adjoint on E'_N + E_N^perp, the inner
    product on E'_N = E_N/C xi being QW_lambda^N - eps_N <|>; (b) det_reg(D_{lambda,N} - z) = -i lambda^{-iz} xi^(z); (c) xi^ is entire
    and all its zeros are real and equal the spectrum of D_{lambda,N}.
    CCM's own numerical cases lambda = 3, sqrt 12, sqrt 13, sqrt 14 with N = 120 are covered (N_0 = 80, 80, 100, 100): the reality
    of the zeros they computed now follows from their Theorem 5.10.
(iii) Not covered: lambda > sqrt 20, and the second missing step (k_lambda close to xi_lambda; lambda -> infinity). Between chain points
    the truncation statement holds only for N >= some N_0(c) that is not explicit (E is a form core, CCM Prop. 3.3).

## 2. The paper (feasibility, sec. 1 of the brief)
- QW_lambda (CCM eq. (3.10)): archimedean term int |f^|^2 2 theta'(t) dt/2 pi, 2 theta' = -log pi + Re psi(1/4 + it/2); pole term
  2 Re(f^(i/2) conj f^(-i/2)); primes 1 < n <= lambda^2 with <f|T(n)f> = n^{-1/2}((f^* * f)(n) + (f^* * f)(1/n)). Same form as ours
  with c = lambda^2, so our window c is lambda = sqrt c as the brief guessed.
- Thm 5.10 = Theorem "finmain" (sec. 5): hypothesis "smallest eigenvalue eps_N of the truncated matrix QW_lambda^N simple, eigenvector
  even"; conclusion (a)-(c) above. Sec. 8 (l. 1390-1396): to use it for QW_lambda itself one must prove that the smallest eigenvalue of
  QW_lambda (exists by Thm 3.6) is simple with even eigenvector; known for the prolate operator PW_lambda; listed as missing step one.
- Their numerics: lambda = 3, sqrt 12, sqrt 13, sqrt 14, N = 120, 200 digits; figure of log eps_lambda against mu = lambda^2 (about
  -12 mu); simplicity is assumed, not checked rigorously.

## 3. Method (DERIVATION.md)
Three certified inequalities per window c:
- (a) U = Q(g)/||g||^2 >= lambda_1^even for g = the (inverse-iteration) ground vector of the even block of QW_lambda^NG, rounded to exact
  dyadics; the block is assembled in Arb from closed forms (weil_matrix.py: pole and prime parts elementary, archimedean part by
  digamma and Hurwitz-zeta values at 1/4 + i pi j/L plus geometric c^{-b_k} sums with an enclosed tail).
- (b) Q(f) >= L_odd ||f||^2 for all odd f: ../WeilZeta's cert_zeta.py (Chekmachev's periodized lattice minorant, Woodbury bound,
  Cholesky positivity proof), copied unchanged into cert_simple.py.
- (c) Q(f) - L_2 ||f||^2 + <f, g>^2 >= 0 for all even f (new "defl" sector: g added as one more positive rank-one column), so even f
  orthogonal to g have Q(f) >= L_2 ||f||^2 and lambda_2^even >= L_2 (Courant-Fischer).
U < min(L_odd, L_2) gives (i) at c; g in E_NG gives (ii) for all N >= NG. For a continuum: every lambda_k of each parity sector is
non-increasing in c (L^2(I_c) sits in L^2(I_c'), same form), so U(c_i) < L(c_{i+1}) covers all of [c_i, c_{i+1}]; chain.py builds 30
such steps from 2 to 16 and 7 more to 20 (37 in all), with L(c_{i+1}) = 2.00001 U(c_i).

## 4. Results
Integer windows (logs/int.jsonl; "seconds": the two certificates, in brackets including the floating-point parameter search):
| c = lambda^2 | lambda | U | L_odd | L_2 | NG | bits | seconds: certificates (with parameter search) | peak MB | margins (defl, odd) |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 1.4142 | 1.335e-3 | 1.026e-2 | 2.971e-2 | 30 | 192 | 1 (1) | 48 | 3.4e-1, 8.2e-2 |
| 3 | 1.7321 | 5.967e-8 | 1.023e-6 | 9.431e-6 | 30 | 192 | 1 (2) | 48 | 9.9e-4, 1.6e-5 |
| 4 | 2.0000 | 9.637e-13 | 2.290e-11 | 4.215e-10 | 30 | 192 | 1 (2) | 48 | 9.8e-8, 5.5e-10 |
| 5 | 2.2361 | 1.222e-17 | 3.979e-16 | 9.178e-15 | 30 | 192 | 1 (2) | 48 | 2.9e-12, 7.8e-15 |
| 6 | 2.4495 | 9.575e-23 | 4.174e-21 | 1.232e-19 | 40 | 192 | 1 (2) | 47 | 7.8e-17, 1.9e-19 |
| 7 | 2.6458 | 1.006e-27 | 4.779e-26 | 1.596e-24 | 40 | 256 | 2 (4) | 58 | 8.8e-22, 1.0e-24 |
| 8 | 2.8284 | 5.211e-33 | 2.887e-31 | 1.233e-29 | 60 | 256 | 2 (5) | 48 | 1.3e-26, 1.0e-29 |
| 9 | 3.0000 | 3.229e-38 | 1.987e-36 | 8.893e-35 | 80 | 320 | 3 (12) | 60 | 1.1e-31, 6.2e-35 |
| 10 | 3.1623 | 1.944e-43 | 1.350e-41 | 7.348e-40 | 80 | 320 | 3 (14) | 54 | 9.6e-37, 7.2e-40 |
| 11 | 3.3166 | 1.371e-48 | 1.047e-46 | 5.706e-45 | 80 | 320 | 10 (24) | 134 | 7.8e-42, 3.5e-45 |
| 12 | 3.4641 | 6.271e-54 | 5.571e-52 | 3.583e-50 | 80 | 384 | 8 (23) | 102 | 7.7e-47, 5.6e-50 |
| 13 | 3.6056 | 3.721e-59 | 3.526e-57 | 2.328e-55 | 100 | 384 | 12 (44) | 97 | 5.3e-52, 1.5e-55 |
| 14 | 3.7417 | 1.642e-64 | 1.749e-62 | 1.309e-60 | 100 | 448 | 18 (51) | 155 | 4.0e-57, 3.3e-60 |
| 15 | 3.8730 | 8.446e-70 | 9.177e-68 | 7.240e-66 | 100 | 448 | 16 (53) | 146 | 1.8e-62, 8.0e-66 |
| 16 | 4.0000 | 3.743e-75 | 4.581e-73 | 4.086e-71 | 120 | 512 | 24 (96) | 148 | 1.1e-67, 2.9e-71 |
| 17 | 4.1231 | 2.095e-80 | 2.674e-78 | 2.496e-76 | 120 | 512 | 44 (73) | 161 | 7.2e-73, 1.4e-76 |
| 18 | 4.2426 | 6.775e-86 | 9.348e-84 | 9.476e-82 | 160 | 576 | 96 (159) | 312 | 4.8e-78, 1.6e-81 |
| 19 | 4.3589 | 4.223e-91 | 5.809e-89 | 5.768e-87 | 160 | 576 | 126 (207) | 255 | 2.5e-83, 5.5e-87 |
| 20 | 4.4721 | 1.657e-96 | 2.497e-94 | 2.841e-92 | 160 | 576 | 246 (315) | 676 | 1.6e-88, 6.6e-92 |

Chain over [2, 20] (logs/chain.jsonl; at c_{i+1} both certificates at L = L_odd = L_2 = 2.00001 U(c_i); U(c_{i+1}) from the defl run;
margins = smallest eigenvalue of Id - W_up; seconds for both certificates). Chain: 37 records, 37 ok; covered interval: [2, 20] = lambda in [1.4142, 4.4721], 37 steps; chain compute: 1548 s; max RSS 671 MB:
| c_i -> c_{i+1} | U(c_i) | L(c_{i+1}) = eps | U(c_{i+1}) | margins (defl, odd) | NG | bits | s |
|---|---|---|---|---|---|---|---|
| 2 -> 9/4 | 1.34e-3 | 2.67e-3 | 1.45e-4 | 2.2e-1, 2.9e-2 | 30 | 192 | 2.5 |
| 9/4 -> 5/2 | 1.45e-4 | 2.90e-4 | 1.05e-5 | 6.1e-2, 2.5e-3 | 30 | 192 | 2.6 |
| 5/2 -> 67/24 | 1.05e-5 | 2.10e-5 | 4.95e-7 | 6.7e-3, 8.9e-5 | 30 | 192 | 2.5 |
| 67/24 -> 25/8 | 4.95e-7 | 9.90e-7 | 1.52e-8 | 5.4e-4, 1.6e-5 | 30 | 192 | 2.6 |
| 25/8 -> 83/24 | 1.52e-8 | 3.03e-8 | 4.22e-10 | 2.3e-5, 2.8e-7 | 30 | 192 | 2.5 |
| 83/24 -> 23/6 | 4.22e-10 | 8.45e-10 | 6.12e-12 | 6.1e-7, 3.6e-9 | 30 | 192 | 2.5 |
| 23/6 -> 101/24 | 6.12e-12 | 1.22e-11 | 1.08e-13 | 1.2e-8, 7.7e-11 | 30 | 192 | 2.5 |
| 101/24 -> 55/12 | 1.08e-13 | 2.17e-13 | 1.28e-15 | 3.0e-10, 1.1e-12 | 30 | 192 | 2.5 |
| 55/12 -> 5 | 1.28e-15 | 2.57e-15 | 1.22e-17 | 3.1e-12, 7.4e-15 | 30 | 192 | 2.6 |
| 5 -> 65/12 | 1.22e-17 | 2.44e-17 | 9.26e-20 | 5.0e-14, 2.8e-16 | 40 | 192 | 2.7 |
| 65/12 -> 47/8 | 9.26e-20 | 1.85e-19 | 4.56e-22 | 3.4e-16, 8.6e-19 | 40 | 192 | 2.8 |
| 47/8 -> 19/3 | 4.56e-22 | 9.12e-22 | 2.44e-24 | 1.5e-18, 2.9e-21 | 40 | 256 | 2.8 |
| 19/3 -> 163/24 | 2.44e-24 | 4.88e-24 | 1.03e-26 | 9.2e-21, 9.0e-24 | 40 | 256 | 2.8 |
| 163/24 -> 29/4 | 1.03e-26 | 2.05e-26 | 3.68e-29 | 8.4e-23, 1.9e-25 | 60 | 256 | 3.9 |
| 29/4 -> 31/4 | 3.68e-29 | 7.35e-29 | 1.08e-31 | 2.4e-25, 3.5e-28 | 60 | 256 | 3.7 |
| 31/4 -> 33/4 | 1.08e-31 | 2.16e-31 | 2.96e-34 | 8.3e-28, 1.4e-30 | 80 | 320 | 5.1 |
| 33/4 -> 35/4 | 2.96e-34 | 5.91e-34 | 5.46e-37 | 2.4e-30, 1.4e-33 | 80 | 320 | 4.5 |
| 35/4 -> 223/24 | 5.46e-37 | 1.09e-36 | 9.13e-40 | 5.0e-33, 4.5e-36 | 80 | 320 | 6.0 |
| 223/24 -> 59/6 | 9.13e-40 | 1.83e-39 | 1.43e-42 | 9.2e-36, 4.8e-39 | 80 | 320 | 5.9 |
| 59/6 -> 83/8 | 1.43e-42 | 2.85e-42 | 2.49e-45 | 1.3e-38, 7.3e-42 | 80 | 320 | 6.6 |
| 83/8 -> 131/12 | 2.49e-45 | 4.98e-45 | 3.33e-48 | 2.1e-41, 7.1e-45 | 80 | 320 | 5.6 |
| 131/12 -> 275/24 | 3.33e-48 | 6.65e-48 | 4.50e-51 | 4.6e-44, 3.5e-47 | 80 | 384 | 10.5 |
| 275/24 -> 12 | 4.50e-51 | 8.99e-51 | 6.27e-54 | 7.7e-47, 3.6e-50 | 80 | 384 | 7.3 |
| 12 -> 151/12 | 6.27e-54 | 1.25e-53 | 4.66e-57 | 7.2e-50, 1.6e-53 | 100 | 384 | 11.9 |
| 151/12 -> 79/6 | 4.66e-57 | 9.32e-57 | 4.05e-60 | 9.3e-53, 4.4e-56 | 100 | 448 | 25.5 |
| 79/6 -> 55/4 | 4.05e-60 | 8.11e-60 | 3.63e-63 | 7.3e-56, 3.2e-59 | 100 | 448 | 18.1 |
| 55/4 -> 43/3 | 3.63e-63 | 7.26e-63 | 3.12e-66 | 6.9e-59, 2.6e-62 | 100 | 448 | 29.1 |
| 43/3 -> 179/12 | 3.12e-66 | 6.23e-66 | 2.31e-69 | 4.7e-62, 1.2e-65 | 100 | 448 | 27.9 |
| 179/12 -> 31/2 | 2.31e-69 | 4.63e-69 | 1.74e-72 | 4.5e-65, 8.2e-69 | 120 | 512 | 31.1 |
| 31/2 -> 16 | 1.74e-72 | 3.48e-72 | 3.74e-75 | 1.1e-67, 2.5e-71 | 120 | 512 | 22.3 |
| 16 -> 133/8 | 3.74e-75 | 7.49e-75 | 2.00e-78 | 6.3e-71, 7.0e-75 | 120 | 512 | 41.2 |
| 133/8 -> 69/4 | 2.00e-78 | 4.00e-78 | 6.96e-82 | 4.1e-74, 6.2e-78 | 160 | 576 | 119.6 |
| 69/4 -> 143/8 | 6.96e-82 | 1.39e-81 | 3.33e-85 | 2.2e-77, 5.0e-81 | 160 | 576 | 86.1 |
| 143/8 -> 37/2 | 3.33e-85 | 6.66e-85 | 1.62e-88 | 1.0e-80, 1.8e-84 | 160 | 576 | 131.7 |
| 37/2 -> 153/8 | 1.62e-88 | 3.25e-88 | 7.84e-92 | 5.3e-84, 7.6e-88 | 160 | 576 | 286.1 |
| 153/8 -> 79/4 | 7.84e-92 | 1.57e-91 | 3.81e-95 | 3.4e-87, 1.0e-90 | 160 | 576 | 317.6 |
| 79/4 -> 20 | 3.81e-95 | 7.62e-95 | 1.66e-96 | 1.6e-88, 6.7e-92 | 160 | 576 | 307.2 |

Two-sided bounds for the ground eigenvalue mu_lambda (= CCM's eps_lambda) at the integer windows (bracket.py, logs/bracket.jsonl):
L_1 from the plain even certificate at eps = U/2 (U/5 where U/2 failed: c = 11, 12, 14-17, 19, 20).
| c | L_1 <= mu_lambda | mu_lambda <= U | U/L_1 |
|---|---|---|---|
| 2 | 6.6759e-04 | 1.3352e-03 | 2 |
| 3 | 2.9833e-08 | 5.9666e-08 | 2 |
| 4 | 4.8185e-13 | 9.6369e-13 | 2 |
| 5 | 6.1089e-18 | 1.2218e-17 | 2 |
| 6 | 4.7874e-23 | 9.5748e-23 | 2 |
| 7 | 5.0282e-28 | 1.0056e-27 | 2 |
| 8 | 2.6056e-33 | 5.2112e-33 | 2 |
| 9 | 1.6143e-38 | 3.2286e-38 | 2 |
| 10 | 9.7225e-44 | 1.9445e-43 | 2 |
| 11 | 2.7411e-49 | 1.3705e-48 | 5 |
| 12 | 1.2541e-54 | 6.2705e-54 | 5 |
| 13 | 1.8604e-59 | 3.7209e-59 | 2 |
| 14 | 3.2842e-65 | 1.6421e-64 | 5 |
| 15 | 1.6891e-70 | 8.4456e-70 | 5 |
| 16 | 7.4860e-76 | 3.7430e-75 | 5 |
| 17 | 4.1909e-81 | 2.0954e-80 | 5 |
| 18 | 3.3875e-86 | 6.7751e-86 | 2 |
| 19 | 8.4454e-92 | 4.2227e-91 | 5 |
| 20 | 3.3149e-97 | 1.6574e-96 | 5 |

Node sets: K2 (all nodes k <= K2 kept in the minorant) was cert_zeta's default 201 in the integer rows c = 2-6, 8, 10 (run before
the failure at c = 7, 9 was seen), max(401, ceil(5c^2/3)) in all other rows and in the whole chain; at c = 11 the odd certificate
needed K2 = 803 (second try). The exact arguments of every run are in the logs.

Approximate spectrum (eigs_quick.py, floating point, orientation only): the low spectrum alternates
even < odd < even < odd with ratios 60-1e4 (c = 2: 1.33e-3, 0.077, 0.66; c = 16: 3.2e-75, 4.7e-71, 3.7e-67, 2.0e-63 at N = 160).

## 5. Checks
- Entries: every entry with i, j <= 2 of both blocks at c = 5 (logs/test_entries_c5.log) and six entries with indices up to 5 at
  c = 25/8 (logs/test_entries_c25o8.log) agree with independent mpmath quadrature in t (archimedean) and x (pole, primes) to 1e-6..1e-7,
  the quadrature's truncation error; logs/test_hp.log: ten entries (both blocks, indices up to 4, c = 5 and
  25/8) equal a 50-digit mpmath quadrature of the x-space form of the archimedean term (DERIVATION.md sec. 2; autocorrelations, pole and
  prime parts also by quadrature) to 3.5e-49, the quadrature's precision: the digamma, Hurwitz-zeta and c^{-b_k} closed forms are exact.
- Spectrum against independent code: U(c) is at or above weilzeta's certified lower bounds of lambda_1^even (lattice minorant, no
  shared formulas): c = 5: 1.2218e-17 >= 7.447e-18; c = 10: 1.9445e-43 >= 9.134e-44; c = 16: 3.7430e-75 >= 1.5128e-75.
- cert_simple.py reproduces weilzeta's c = 5 margins to all 12 digits (even 4.00280318265e-19, odd 8.28212749750e-16).
- Negative controls at c = 5 (logs/controls_c5.log): defl at eps = 1e-11 and 6e-12 (above lambda_2^even ~ 6e-12) fails; odd at
  2e-14 (above lambda_1^odd ~ 1.1e-14) fails; even without deflation at 1e-15 fails, with deflation at 1e-15 passes (margin 2.9e-12).
- A failure found and fixed: with cert_zeta's default K2 = 201 the odd minorant at c = 7 and c = 9 has a negative eigenvalue that does
  not depend on eps (margin -2.6e-24, -2.27e-28; logs/int_k201_fail.jsonl); K2 = max(401, 5c^2/3) removes it (the minorant drops
  positive terms beyond K2). Not a soundness issue (a failing certificate proves nothing), but ../WeilZeta's default K2 is not safe at
  small c.

## 6. Rerun (smallest window, under 10 s)
From Riemann/Deep/WeilSimple, one single-threaded process:
    sh rerun.sh
(cc -O2 -ffp-contract=off rscan.c, then run_windows.py on c = 2: prints U = 1.335e-3 < L_odd = 1.026e-2, L_2 = 2.971e-2, even_simple
true; measured 5.8 s with a fresh scan cache; margins differ in the 4th digit from logs/int.jsonl because the logged row used K2 = 201.) A chain segment: uv run --no-project --with mpmath --with numpy --with python-flint python3 chain.py logs/int.jsonl /tmp/ch.jsonl 2 5
(about 30 s); python3 table.py re-checks every chain inequality from the logs.

## 7. Novelty pass (2026-11-05)
- python3 .scratch/lst.py 'prolate|Weil|spectral triple|Connes' math.NT math.SP math.OA --zenodo '"Weil quadratic form"': 9 arXiv
  title hits (2026-10, 2026-11), none on simplicity (2610.32704 assumes RH; 2610.32082 prolate matrices). Zenodo: see below.
- arXiv API: all:"Weil quadratic form" 9 papers (newest 2610.41115, 2610.32704), abs:"zeta spectral triples" 0, all:"2511.22755" 0
  (the API has no citation search). Citations of 2511.22755: Semantic Scholar 9, OpenAlex 57 works (titles screened); Zenodo searches
  "even-simplicity", "zeta spectral triples", "compact windows". Relevant:
  - Zhu, arXiv 2608.24827 v2 (2 Sep 2026), Theorem "Simple even ground state, support 1.6": the same statement for ONE window,
    supp f in [-0.8, 0.8] (lambda = e^0.8 = 2.2255, c = 4.953): 8.9e-18 <= lambda_1^even <= 2.523e-16, lambda_2^even >= 2.085e-12,
    lambda_1^odd >= 8.206e-15, by the same three inequalities (Courant-Fischer on a Householder complement of the approximate ground
    vector of a 200-mode head, then a block-perturbation lift). Our result contains this window and extends it to the continuum
    lambda in [sqrt 2, sqrt 20] (half-width 0.347 to 1.498). Our second-eigenvalue step (rank-one deflation inside the lattice minorant) is different.
  - Suzuki, arXiv 2606.09096 Thm 4: simple and even for sufficiently small half-width a (perturbative, not explicit).
  - de Andrade Silva, Zenodo 20887027 (25 Jun 2026): even-simplicity of three FINITE matrices (c, N) = (13, 120), (29, 160), (53, 200);
    states it does not treat the continuum operator. Our c = 13 row covers every N >= 100, so includes (13, 120).
  - Andrews, Zenodo 23905446 v2.6 (Nov 2026): certified finite controls at C = 13, N = 10 and 120; cites Zhu as the continuum certificate.
  - Mondal, Zenodo 23187644 (6 Oct 2026): infinite-tail estimates for 17 <= c <= 19, intermediate modes not bounded, no simplicity theorem.
  - Luekens, Zenodo 23752720 (1 Nov 2026): Galerkin measurements c = 9..150 (gap to the second even state ~ 1900 c^4 eps_1, which our
    eigenvalues match: c = 16 ratio 1.2e8); no continuum simplicity.
  - Mehennaoui, Zenodo 23441946/23442026: "simplicity implies evenness" given the normalization delta_N(xi) = 1 (which presupposes it);
    Franco Olguin, Zenodo 23228212 ("analytic proof of even-simplicity" in the title): a rank-one pole identity plus numerics;
    "even dominance" papers (Zenodo 20011910, 20291994): finite values for 100 <= lambda <= 1.3e6, "not theorem-level because the
    odd-sector lower-bound construction remains open"; McDermott 20651102: N = 1 criterion under RH.
- Our files: R359 records the problem as open; TREE N1 lists Zhu's positivity theorem but not its parity theorem; ../WeilZeta, ../ChekAudit:
  lower bounds only.
- Verdict: new as far as found (first certificate on a continuum of windows; 52 explicit windows with all-N truncation statements);
  the method is Zhu's three-inequality scheme with a different second-eigenvalue certificate; not a step toward lambda -> infinity.

## 8. Limits and trust
- Trust: Arb/python-flint 0.9.0 enclosures; IEEE binary64 +, -, * in rscan.c (RIGOR.md of ../WeilZeta, flag margin 1e-6, flagged nodes
  re-checked in Arb); this Python, unchecked by a proof assistant. No row rests on an unproved float budget: g and the Cholesky factors
  are floating point but enter only as exact inputs to enclosed residuals.
- Closed forms of sec. 2 of DERIVATION.md are our derivation; checked numerically (sec. 5), not formally.
- Scope: cost grows like the certificates of ../WeilZeta (c = 16: 24 s, 150 MB; c = 20: 4 min, 0.68 GB); about c = 25 is the limit
  for 1.5 GB in this layout (../WeilZeta: c = 25 at 576 bits used 1.08 GB; c >= 30 does not fit). Nothing here bears on lambda -> infinity (mu_lambda >= 0 for every lambda is equivalent to RH).

## 9. Files
weil_matrix.py (closed-form blocks of QW_lambda^N); cert_simple.py (cert_zeta.py + defl sector; marked WS); rscan.c (copied unchanged
from ../WeilZeta, md5 bd77d55f0a0bff7e7955e6c734aab9a9; cert_zeta.py md5 8e6e503f680d9c56a40d0e87b0b63803 at copy time);
patch_cert.py (the exact edit cert_zeta.py -> cert_simple.py); run_windows.py, chain.py (drivers); table.py (tables + chain check);
eigs_quick.py, test_matrix.py, test_hp.py (checks); rerun.sh; DERIVATION.md; logs/.

## 10. A401172 extension (2026-11-06, follow-up)
OEIS A401172: a(n) = floor(-log10 kappa(n)), kappa(n) = inf Q(f)/||f||^2 over f supported in (-log(n)/2, log(n)/2). Its %C formula
(pole 2 Re(F(i/2) conj F(-i/2)); archimedean (1/2 pi) int (Re psi(1/4 + it/2) - log pi) |F|^2 dt; primes sum 2 log p p^{-k/2}
Re int f(u + k log p) conj f(u) du) is term by term our Q at c = n: the pole identity is DERIVATION.md sec. 1, and u = y - k log p
turns int f(y) conj f(y - log n) dy (n = p^k, Lambda(n) n^{-1/2} = log p p^{-k/2}) into its prime term. So kappa(n) = mu_lambda at
c = n, and since lambda_1^odd(n) >= L_odd(n) > U(n) (sec. 4), kappa(n) = lambda_1^even(n).
Lower bounds: the plain even certificate (cert_simple.py even sector = cert_zeta.py) at eps proves kappa(n) >= eps (bracket.py,
abracket.py). Upper bounds: Rayleigh quotients of the ground vector of the even block of QW_lambda^N in Arb (int rows, uimprove.py up
to N = 400). a(n) is certified when both ends of the interval give the same floor. Table from the logs: python3 a401172.py.
| n | certified interval for kappa(n) | -log10 kappa(n) in | a(n) | OEIS data | sources (lower; upper) |
|---|---|---|---|---|---|
| 2 | [0.00100001, 0.00132983] | [2.87620496, 2.99999566] | 2 | 2 | even cert, K2 = 401; g in E_120 |
| 3 | [2.9833e-8, 5.96656e-8] | [7.22427622, 7.52530307] | 7 | 7 | even cert, K2 = max(401, 5n^2/3); g in E_30 |
| 4 | [4.8185e-13, 9.63694e-13] | [12.0160607, 12.3170881] | 12 | 12 | even cert, K2 = max(401, 5n^2/3); g in E_30 |
| 5 | [6.1089e-18, 9.28104e-18] | [17.0324033, 17.214037] | 17 | 17 | even cert, K2 = max(401, 5n^2/3); g in E_240 |
| 6 | [4.7874e-23, 9.57482e-23] | [22.0188695, 22.3199003] | 22 | 22 | even cert, K2 = max(401, 5n^2/3); g in E_40 |
| 7 | [5.0282e-28, 7.4023e-28] | [27.1306332, 27.2985875] | 27 | 27 | even cert, K2 = max(401, 5n^2/3); g in E_160 |
| 8 | [2.6056e-33, 5.21118e-33] | [32.2830639, 32.5840923] | 32 | 32 | even cert, K2 = max(401, 5n^2/3); g in E_60 |
| 9 | [1.6143e-38, 3.22857e-38] | [37.4909894, 37.7920158] | 37 | 37 | even cert, K2 = max(401, 5n^2/3); g in E_80 |
| 10 | [1.00001e-43, 1.57234e-43] | [42.8034542, 42.9999957] | 42 | - | even cert, K2 = 401; g in E_240 |
| 11 | [2.7411e-49, 9.45177e-49] | [48.024487, 48.5620751] | 48 | - | even cert, K2 = max(401, 5n^2/3); g in E_320 |
| 12 | [1.2541e-54, 6.27053e-54] | [53.2026956, 53.9016678] | 53 | - | even cert, K2 = max(401, 5n^2/3); g in E_80 |
| 13 | [1.8604e-59, 3.7209e-59] | [58.429352, 58.7303937] | 58 | - | even cert, K2 = max(401, 5n^2/3); g in E_100 |
| 14 | [3.2842e-65, 1.18064e-64] | [63.9278843, 64.4835704] | 63 or 64 (open) | - | even cert, K2 = max(401, 5n^2/3); g in E_320 |
| 15 | [1.6891e-70, 8.44558e-70] | [69.0733705, 69.7723446] | 69 | - | even cert, K2 = max(401, 5n^2/3); g in E_100 |
| 16 | [1.00001e-75, 3.74299e-75] | [74.4267811, 74.9999957] | 74 | - | even cert, K2 = 427; g in E_120 |
| 17 | [4.1909e-81, 1.40047e-80] | [79.8537248, 80.3776927] | 79 or 80 (open) | - | even cert, K2 = max(401, 5n^2/3); g in E_240 |
| 18 | [3.3875e-86, 6.77506e-86] | [85.1690867, 85.4701207] | 85 | - | even cert, K2 = max(401, 5n^2/3); g in E_160 |
| 19 | [1.00001e-91, 4.22271e-91] | [90.374409, 90.9999957] | 90 | - | even cert, K2 = 603; g in E_160 |
| 20 | [3.3149e-97, 1.19228e-96] | [95.9236215, 96.4795296] | 95 or 96 (open) | - | even cert, K2 = max(401, 5n^2/3); g in E_400 |

Certified: a(n) for n = 2..13, 15, 16, 18, 19 (16 values). n = 2..9 equal the 8 data terms of A401172; new: a(10) = 42, a(11) = 48,
a(12) = 53, a(13) = 58, a(15) = 69, a(16) = 74, a(18) = 85, a(19) = 90.
Open: n = 14, 17, 20. In each case kappa(n) lies within a few percent of a power of 10 (estimates below are fits, not proofs):
- n = 14: U(N) = 1.6420e-64 (N = 100), 1.2653e-64 (200), 1.1806e-64 (320); fits a + b log(N)/N and a + b/N through the last two give
  kappa(14) = 1.00e-64 and 1.04e-64; the even certificate fails at 8.21e-65 (K2 = 401). Deciding needs both bounds to about 1%.
- n = 17: the even certificate at eps = 1.00001e-80 fails with margin -6.07e-82 (K2 = 483) and -2.61e-82 (K2 = 967); a 1/K2 fit of the
  margin turns positive only near K2 = 4000-8000 (about 40 min or more in this layout), i.e. kappa(17) is just above 1e-80 (then
  a(17) = 79); U(240) = 1.4005e-80.
- n = 20: U(N) = 1.6574e-96 (160), 1.3425e-96 (240), 1.1923e-96 (400); fits give kappa(20) = 0.91e-96 to 0.97e-96 (then a(20) = 96);
  a certificate needs U < 1e-96, about N = 1500-2700 (hours here). The run at eps = 1.00001e-96 (logs/abracket.jsonl, rigorous null)
  was stopped by hand: by the margins of logs/bracket.jsonl the K2 = 667 minorant bottoms out near 6e-97, below 1e-96.
Why the bounds are loose: U(N) converges slowly in N (the cosine basis resolves the ground state's boundary behaviour badly; the
relative error falls roughly like log(N)/N), and the minorant loses the nodes beyond K2 at a similar rate.
Rerun of one decision (n = 10, about 15 s; WZ_RSCAN and WZ_SCAN as in rerun.sh): uv run --no-project --with mpmath --with numpy
--with python-flint python3 abracket.py /tmp/a.jsonl 10:1.00001e-43:320:401, then the same with uimprove.py /tmp/u.jsonl 10:240:320.
Files: a401172.py, uimprove.py, abracket.py; logs/uimprove.jsonl, logs/abracket.jsonl, logs/a401172_table.md. All runs: one process,
nice 19, single thread, at most 300 MB.

DONE2

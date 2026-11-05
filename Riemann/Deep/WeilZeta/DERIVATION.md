# DERIVATION: the Weil form of zeta, its lattice identity and the certificate (weilzeta, 2026-11-05)

Port of ../WeilChi/DERIVATION.md (Dirichlet L-functions) back to zeta. Differences from L(s, chi): two pole terms (s = 0, 1),
constant -log pi instead of log(q/pi), gamma factor Gamma(s/2) (b_n = 2n + 1/2), weights Lambda(n)/sqrt n (no character).
The same chain was audited in Chekmachev's paper (Zenodo 23848597; ../ChekAudit/REPORT.md, verdict SOUND).

## 1. Explicit formula (the one the code uses)
Xi(s) = pi^{-s/2} Gamma(s/2) zeta(s) has simple poles at s = 0 and s = 1, zeros exactly at the nontrivial zeros rho = 1/2 + i gamma,
and Xi(s) = Xi(1 - s). Test function g in C_c^2(R), supp g in [-A, A], h(t) = int g(x) e^{itx} dx, H(s) = h(-i(s - 1/2)).
|H(sigma + it)| <= C (1 + |t|)^{-2} e^{A |sigma - 1/2|}; with |zeta'/zeta(sigma + i T_m)| = O(log^2 T_m) on suitable T_m -> infinity
and N(T) = O(T log T) the horizontal integrals vanish (Davenport ch. 15-17; Iwaniec-Kowalski Thm 5.12). Integrating Xi'/Xi(s) H(s)
over the boundary of the strip 1 - c <= Re s <= c (c > 1) and using Xi'/Xi(s) = -Xi'/Xi(1 - s) on the left line:
  sum_rho H(rho) = H(0) + H(1) + (1/2 pi i) int_(c) Xi'/Xi(s) [H(s) + H(1 - s)] ds,
  Xi'/Xi(s) = -(1/2) log pi + (1/2) psi(s/2) - sum_n Lambda(n) n^{-s}   (Re s > 1).
- Poles: residue -1 of Xi'/Xi at s = 0 and at s = 1 (as poles of Xi) give + H(0) + H(1) = h(i/2) + h(-i/2).
- Primes: (1/2 pi i) int_(c) n^{-s} H(s) ds = n^{-1/2} g(log n); H(1 - s) is the H-transform of g(-x). Total -sum Lambda(n) n^{-1/2} [g(log n) + g(-log n)].
- Archimedean: shift to Re s = 1/2 (psi(s/2) has poles only at s = 0, -2, ...; the pole at s = 0 is not crossed for 1/2 < Re s):
  (1/2 pi) int h(t) [-log pi + Re psi(1/4 + it/2)] dt.
  (EF)  sum_rho h(gamma_rho) = h(i/2) + h(-i/2) + (1/2 pi) int h(t) phi(t) dt - sum_n Lambda(n) n^{-1/2} [g(log n) + g(-log n)],
        phi(t) = -log pi + Re psi(1/4 + it/2).
This is the Guinand-Weil formula (Iwaniec-Kowalski (5.44) with the gamma factor of zeta).

## 2. The Weil form and Weil's criterion
For f in L^2(R) with compact support: g(x) = int f(y) conj f(y - x) dy, F(t) = int f(x) e^{itx} dx, so h = |F|^2 on R,
h(z) = F(z) conj F(conj z) for complex z, and (1/2 pi) int |F|^2 = g(0) = ||f||^2. Then h(i/2) + h(-i/2) = 2 Re[F(i/2) conj F(-i/2)]
= 2 |<f, cosh(x/2)>|^2 - 2 |<f, sinh(x/2)>|^2 (write F(+-i/2) = a -+ b with a = int f cosh(x/2), b = int f sinh(x/2)). Define
  Q(f) = 2 |<f, cosh(x/2)>|^2 - 2 |<f, sinh(x/2)>|^2 + (1/2 pi) int phi |F|^2 - 2 sum_n Lambda(n) n^{-1/2} Re g(log n).
By (EF), Q(f) = sum_rho F(gamma) conj F(conj gamma) for f in C_c^infty; under RH this is sum |F(gamma)|^2 >= 0, and an off-line zero
pair makes Q(f) < 0 for a suitable f (Weil). So RH <=> Q(f) >= 0 for all f in C_c^infty(R).
If supp f in [-a, a] only n < e^{2a} = c enter. The form is real symmetric and commutes with x -> -x: Q(u + iv) = Q(u) + Q(v) for real
u, v, and Q(f) = Q(f_even) + Q(f_odd). The pole part is +2 <f, cosh(x/2)>^2 on even f and -2 <f, sinh(x/2)>^2 on odd f.

## 3. Lattice identity (Chekmachev Thm 3.1)
c > 1, L = ln c, I = [-L/2, L/2], r_k = pi k/L. For every f in L^2(I), as an identity in (-infinity, +infinity]:
  (LAT) Q(f) = pole(f) + (1/2L) sum_{k in Z} mu(k) |F(r_k)|^2,
  mu(k) = -log pi + Re psi(1/4 + i r_k/2) + 2 (-1)^k sum_{n>=0} c^{-b_n} b_n/(b_n^2 + r_k^2) - 2 sum_{m<c} Lambda(m) m^{-1/2} cos(r_k ln m),
  b_n = 2n + 1/2.
Proof: ../WeilChi/DERIVATION.md sec. 3 and 3b with q = 1, a = 0, chi = 1 (Fourier series of the 2L-periodized g with coefficients
|F(r_k)|^2 >= 0; Gauss's series psi(z) = -gamma_E + sum [1/(n+1) - 1/(n+z)] and Re 1/(n + 1/4 + ir/2) = 2 b_n/(b_n^2 + r^2), the transform
of e^{-b_n |x|}; Tonelli for the nonnegative parts, absolute convergence for the corrections). The pole part is a bounded form on L^2(I)
(cosh, sinh in L^2(I)) and stays outside the lattice sum. The atom m = c, if c is a prime power, sits at g(+-L) = 0 and is dropped.
Check (ef_check_zeta.py, A = 1.2, so c = e^2.4): (LAT) with 3000 nodes equals the direct form to 3e-31 (even f), 1e-32 (odd f).
Monotonicity in c: f supported in I_{c'} (c' < c) is supported in I_c with the same Q, so a bound at c holds at every c' <= c.

## 4. Sectors and node systems (Chekmachev Thm 3.2)
Even f: phi_j = N_j cos(2 pi j x/L) (N_0 = L^{-1/2}, N_j = (2/L)^{1/2}), orthonormal; e_k = cos(pi k x/L), k odd, ||e_k||^2 = L/2;
pole vector v = cosh(x/2) with weight +2. Odd f: phi_j = (2/L)^{1/2} sin(2 pi j x/L) (j >= 1), e_k = sin(pi k x/L); v = sinh(x/2), weight -2.
  Q(f) = (1/2) sum_j mu(2j) <f, phi_j>^2 + (1/L) sum_{k odd >= 1} mu(k) <f, e_k>^2 +- 2 <f, v>^2,   (1/L) sum_{k odd} <f, e_k>^2 = ||f||^2/2.
Closed forms (re-derived, as in cert2.py of ChekAudit): with A_j = N_j (L/pi)(-1)^j (even) or N_j (-1)^j 4j (L/pi) (odd),
<phi_j, e_k> = A_j beta_k/(k^2 - 4j^2), beta_k = 2k s_k (even) or s_k (odd), s_k = (-1)^{(k-1)/2};
<phi_j, v> = A_j beta_v/(y_v - 4j^2), y_v = -L^2/(4 pi^2), beta_v = -(L/pi) sinh(L/4) (even) or sinh(L/4) (odd);
<e_k, v> = 2 r_k s_k cosh(L/4)/(r_k^2 + 1/4) (even) or s_k cosh(L/4)/(r_k^2 + 1/4) (odd);
||v||^2 = sinh(L/2) + L/2 (even) or sinh(L/2) - L/2 (odd).

## 5. Minorant, dual check, Woodbury bound, positivity
Choose finite K (odd) and E (j) with mu(k) >= t1 = 1 for odd k not in K and mu(2j) >= 0 for j not in E. For eps < 1/2:
  Q(f) - eps ||f||^2 >= sum_j dl_j <f, phi_j>^2 + sum_{k in K} w_k <f, e_k>^2 +- 2 <f, v>^2,  dl_j = (mu(2j) + 1)/2 - eps,  w_k = (mu(k) - 1)/L,
with dl_j >= deff = 1/2 - eps for j not in E. Write the right side as P - N:
  P = P0 + sum_{K+} w_k e_k e_k^T (+ 2 v v^T, even f), P0 = deff Id + sum_{j in E, dl_j > 0} (dl_j - deff) phi_j phi_j^T,
  N = B B^T, columns sqrt(deff - dl_j) phi_j (dl_j <= 0), sqrt(-w_k) e_k (w_k <= 0) (and sqrt 2 v, odd f).
P - N >= 0 iff Id - W >= 0, W = B^T P^{-1} B = X - Y^T S^{-1} Y (Woodbury) with X = B^T P0^{-1} B, Y = U^T P0^{-1} B, S = D^{-1} + U^T P0^{-1} U
(U the K+ columns and, for even f, v; D their weights). All entries of X, Y, S are closed forms: P0^{-1} is diagonal in phi_j, and the sums
over j in E are Cauchy sums F(y) = sum (1/dl_j - 1/deff) A_j^2/(y - 4j^2), F2(y) = sum (...)/(y - 4j^2)^2 (partial fractions).
- Woodbury upper bound (Chekmachev Lemma 5.3): for ANY real Z, (Z - S^{-1}Y)^T S (Z - S^{-1}Y) >= 0 gives
  Y^T S^{-1} Y >= Z^T Y + Y^T Z - Z^T S Z, so W <= W_up := X - Y^T Z - Z^T Y + Z^T S Z = X - (Z^T Y)^T + Z^T (S Z - Y).
  cert_zeta.py takes Z = an approximate floating solve of S Z = Y (exact entries, no error bound needed) and encloses W_up in Arb.
  Id - W >= Id - W_up, so it suffices to prove Id - W_up > 0. The loss W_up - W is quadratic in the error of Z. (WeilChi's
  cert_chi.py used a rigorous arb_mat solve here; that solve was its memory wall.)
- Positivity (Chekmachev Lemma 5.5): R = mid(Id - W_up); lambda from Rayleigh-Ritz on 6 block inverse-iteration vectors; s0 = lambda/2
  (then lambda/8, lambda/100), nu = s0/2. An approximate Cholesky factor L of R - s0 Id is computed in floating point at the working
  precision (block version: panel solves and trailing updates by arb_mat products, no error bounds needed). The residual
  E = (Id - W_up) - L L^T - s0 Id is enclosed in Arb; if max row sum of |E| <= s0 - nu then lambda_min(Id - W_up) >= nu > 0.
- Scan: for k > kscan, mu(k) >= h(r_k) - 2 sigma - 2 tail(r_k) - 2 trem with h increasing and tail decreasing (DERIVATION of WeilChi sec. 5),
  so kscan = first k with that bound >= t1 (odd) or >= 0 (even), found in Arb. Below kscan every node is checked by rscan.c (or rscan2.c,
  with early exit) with the proved error bound of RIGOR.md (E < margin/10; flag margin 1e-6, and 1e-5 for the c = 45, 47 scans where
  E = 1.8e-7; flagged nodes re-evaluated in Arb at 96 bits).
  Neighbour rule (Chekmachev Def. 6.1): E gets j in [(k-1)/2 - 32, (k-1)/2 + 32] around every selected odd k; every j with mu(2j) < 0
  adds odd k = 2j +- 1 and even j +- 3. K and E contain all odd k <= K2 and j <= K2, K2 = max(201, ceil(5c^2/3)) by default.

## 6. Conventions checked
- Pole sign and weight: + 2 (int f cosh(x/2))^2 (even), - 2 (int f sinh(x/2))^2 (odd); confirmed against zeros (ef_check_zeta.py:
  pole terms 1.0534 and -0.00105, the form 7.93e-9 and 4.94e-8, zero sums equal to 4e-22).
- Factor 2 between f and g: (1/2 pi) int |F|^2 = ||f||^2 = g(0); sum_k |F(r_k)|^2 = 2L ||f||^2.
- Fourier sign: Chekmachev's f^(r) = int f e^{-irx} = F(-r); every term is even in r for real f.
- Validation: cert_zeta.py reproduces ChekAudit's independent margins at c = 5, 10, 16 to all printed digits (logs/val.log).

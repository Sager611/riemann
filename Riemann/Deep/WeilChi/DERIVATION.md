# DERIVATION: the Weil form of a Dirichlet L-function and its lattice minorant (weilchi, 2026-11-05)

## 1. Explicit formula for L(s, chi)
chi primitive mod q > 1 (so not principal), a = 0 if chi(-1) = 1, a = 1 if chi(-1) = -1.
Lambda(s, chi) = (q/pi)^{(s+a)/2} Gamma((s+a)/2) L(s, chi) is entire, its zeros are the nontrivial zeros rho = 1/2 + i gamma of L(s, chi),
and Lambda(s, chi) = eps(chi) Lambda(1-s, conj chi), |eps| = 1.
Test function: g in C_c^infty(R) (C_c^2 suffices), h(t) = int g(x) e^{itx} dx, H(s) = h(-i(s - 1/2)) = int g(x) e^{(s-1/2)x} dx.
If supp g in [-A, A], H is entire and, by two integrations by parts, |H(sigma + it)| <= C (1 + |t|)^{-2} e^{A |sigma - 1/2|}. With
|Lambda'/Lambda(sigma + iT_m, chi)| = O(log^2 T_m) on horizontal segments T_m -> infinity chosen away from zeros (Davenport, Multiplicative
Number Theory, ch. 15-17) and N(T, chi) = O(T log qT), the horizontal integrals vanish and sum_rho H(rho) converges absolutely. This is
the standard proof (Iwaniec-Kowalski Thm 5.12); continuity of g alone would not suffice.
Contour argument (lines Re s = c > 1 and 1 - c; on the left line use Lambda'/Lambda(s, chi) = -Lambda'/Lambda(1-s, conj chi), then s -> 1-s):
  sum_rho H(rho) = (1/2 pi i) int_(c) [ (Lambda'/Lambda)(s, chi) H(s) + (Lambda'/Lambda)(s, conj chi) H(1-s) ] ds.
On Re s = c: Lambda'/Lambda(s, chi) = (1/2) log(q/pi) + (1/2) psi((s+a)/2) - sum_n chi(n) Lambda(n) n^{-s}.
- Prime terms: (1/2 pi i) int_(c) n^{-s} H(s) ds = n^{-1/2} g(log n) (shift to Re s = 1/2, Fourier inversion); H(1-s) is the H-transform of g(-x).
  Total: - sum_n Lambda(n) n^{-1/2} [ chi(n) g(log n) + conj chi(n) g(-log n) ].
- Archimedean terms: shift to Re s = 1/2 (the poles of psi((s+a)/2) are at s = -a - 2k < 1/2): with s = 1/2 + it,
  (1/2 pi) int h(t) [ log(q/pi) + Re psi(1/4 + a/2 + it/2) ] dt.
- No pole terms: L(s, chi) is entire (chi not principal). (For zeta the same computation gives the extra h(i/2) + h(-i/2).)
Result (gamma = (rho - 1/2)/i, complex if rho is off the line):
  (EF)  sum_rho h(gamma_rho) = (1/2 pi) int h(t) phi_chi(t) dt - sum_n Lambda(n) n^{-1/2} [chi(n) g(log n) + conj chi(n) g(-log n)],
  phi_chi(t) = log(q/pi) + Re psi(1/4 + a/2 + it/2).
Check (ef_check.py, logs/ef_check.log): chi_-4 and chi_-3, even and odd test functions f = (1-(x/A)^2)^8 and (x/A)(1-(x/A)^2)^8
on [-1.2, 1.2] (prime powers < e^2.4 = 11.02), zeros found by sign changes of the real function e^{i theta(t)} L(1/2+it, chi) up to
t = 90 (44 zeros for chi_-4, first 6.020948905; 40 for chi_-3, first 8.039737156): zero sum = right side of (EF) to <= 6.7e-22
(values 5.3e-2 ... 1.8e-3), i.e. relative <= 4e-19.

## 2. The Weil form and GRH (smooth f) and its extension to L^2 (sec. 3)
For f in L^2(R) with compact support put g(x) = int f(y) conj f(y - x) dy (so g(-x) = conj g(x), supp g in supp f - supp f) and
F(t) = int f(x) e^{itx} dx, so h = |F|^2 and (1/2 pi) int |F|^2 = ||f||^2. Define
  Q_chi(f) := (1/2 pi) int phi_chi(t) |F(t)|^2 dt - 2 sum_n Lambda(n) n^{-1/2} Re[ chi(n) g(log n) ].
By (EF), Q_chi(f) = sum_rho F(gamma) conj F(conj gamma). Zeros of L(s, chi) are symmetric under rho -> 1 - conj rho, i.e. gamma -> conj gamma,
so under GRH(chi) Q_chi(f) = sum |F(gamma)|^2 >= 0; conversely an off-line pair (gamma, conj gamma) makes Q_chi negative for a suitable f
(Weil's criterion for one L-function). So GRH for L(s, chi) <=> Q_chi(f) >= 0 for all f in C_c^infty(R).
If supp f in [-A, A], only n < e^{2A} enter (g(+-log n) = 0 for log n >= 2A).
Operator form: Q_chi(f) = (1/2 pi) int phi_chi |F|^2 - <f, C_chi f>, C_chi = sum_{n} Lambda(n) n^{-1/2} (conj chi(n) S_{log n} + chi(n) S_{-log n}),
S_s f(x) = f(x - s), <u, v> = int u conj v; C_chi is Hermitian. For real chi: C_chi = sum chi(n) Lambda(n) n^{-1/2} (S_{log n} + S_{-log n}),
real symmetric and commuting with x -> -x, so Q_chi splits into even and odd f with no cross term, and real f suffice.

## 3. Lattice identity (port of Chekmachev, Zenodo 23848597, Thm 3.1)
c > 1, L = ln c, I = [-L/2, L/2], f in L^2(I), r_k = pi k / L. g is continuous, supported in [-L, L], g(+-L) = 0, so its 2L-periodic
extension has Fourier coefficients |F(r_k)|^2 >= 0 (sum = 2L ||f||^2, absolutely convergent), g(x) = (1/2L) sum_k |F(r_k)|^2 e^{-i r_k x}.
- Archimedean part: psi(z) = -gamma_E + sum_{n>=0} [1/(n+1) - 1/(n+z)], z = 1/4 + a/2 + ir/2, and Re 1/(n+z) = 2 b_n/(b_n^2 + r^2) with
  b_n = 2n + 1/2 + a. 2b/(b^2 + r^2) = int e^{-b|x|} e^{irx} dx, so the x-space kernel of phi_chi is (delta-type constants at 0) - sum_n e^{-b_n |x|}.
  int_{-L}^{L} e^{-b|x|} e^{-i r_k x} dx = 2b(1 - (-1)^k c^{-b})/(b^2 + r_k^2). Constants and the singular part at 0 periodize exactly.
- Prime part: 2 Re[chi(n) g(log n)] = (1/2L) sum_k |F(r_k)|^2 2 Re[chi(n) e^{-i r_k log n}]; real chi: 2 chi(n) cos(r_k log n).
  All certified characters (chi_-4, chi_-3, chi_5, chi_8) are real. For complex chi the prime multiplier would be 2 Re[chi(m) e^{-i r_k log m}], mu would
  not be even in k, and the even and odd sectors would couple; that case is not treated here.
Hence, for real chi (all f in L^2(I), the form possibly +infinity):
  (LAT) Q_chi(f) = (1/2L) sum_{k in Z} mu(k) |F(r_k)|^2,
  mu(k) = log(q/pi) + Re psi(1/4 + a/2 + i r_k/2) + 2 (-1)^k sum_{n>=0} c^{-b_n} b_n/(b_n^2 + r_k^2) - 2 sum_{m<c} chi(m) Lambda(m) m^{-1/2} cos(r_k log m).
  Differences from zeta: no pole (rank-one) term; b_n = 2n + 3/2 (odd chi) instead of 2n + 1/2; constant log(q/pi) instead of -log pi;
  prime weights chi(m) Lambda(m)/sqrt m. The term m = c (if c is a prime power) sits at g(+-L) = 0 and may be dropped.
Check: (LAT) with K = 3000 nodes = direct (EF) right side to <= 1.1e-31 (ef_check.py, same functions as in sec. 1).

### 3b. Q_chi on all of L^2(I), and (LAT) there (proof)
The certificate is a statement about the functional Q_chi(f) := (1/2 pi) int phi_chi |F|^2 - 2 sum_{n<c} Lambda(n) n^{-1/2} Re[chi(n) g(log n)]
for f in L^2(I). It is well defined with values in (-infinity, +infinity]: phi_chi is bounded below, and |g(s)| <= ||f||^2.
(EF) is used only for smooth f (sec. 1-2); the certificate needs only (LAT) on L^2(I), proved directly:
- For f in L^2(I), g = f * f~ is continuous with supp g in [-L, L] and g(+-L) = 0. Its 2L-periodic Fourier coefficients are |F(r_k)|^2 >= 0
  with sum 2L ||f||^2 < infinity. So the Fourier series converges absolutely and uniformly to g (Fejer means converge to the continuous g,
  and a series with nonnegative summable coefficients converges absolutely). This gives the prime terms.
- Write phi_chi(t) = C0 + sum_{n>=0} p_n(t), with C0 = log(q/pi) + psi(1/4 + a/2) and p_n(t) = 2/b_n - 2 b_n/(b_n^2 + t^2) = 2t^2/(b_n (b_n^2 + t^2)) >= 0.
  By Tonelli, (1/2 pi) int phi_chi |F|^2 = C0 ||f||^2 + sum_n (1/2 pi) int p_n |F|^2, each term in [0, infinity].
- For each n, |F|^2 is in L^1, the transform of e^{-b_n |x|} is 2 b_n/(b_n^2 + t^2) in L^1, and g is in L^1 and continuous. So
  (1/2 pi) int p_n |F|^2 = (2/b_n) g(0) - int e^{-b_n |x|} g(x) dx.
  Insert the absolutely convergent series for g; the integral over [-L, L] of a bounded function may be swapped with the sum. This gives
  (1/2L) sum_k |F(r_k)|^2 [p_n(r_k) + 2 (-1)^k c^{-b_n} b_n/(b_n^2 + r_k^2)].
- Sum over n. The p_n(r_k) >= 0 parts may be exchanged by Tonelli. The corrections are absolutely summable:
  sum_n sum_k |F(r_k)|^2 c^{-b_n} b_n/(b_n^2 + r_k^2) <= 2L ||f||^2 sum_n c^{-b_n}/b_n < infinity.
  So (LAT) holds for every f in L^2(I) as an identity in (-infinity, +infinity], with the same mu(k).
- Consequence. The node-system and minorant steps (sec. 4-5) are inequalities between such series. If the right side of (LAT) is +infinity, Q(f) >= eps ||f||^2 is trivial.
  So the certificate gives Q_chi(f) >= eps ||f||^2 for every f in L^2(I), in particular for all f in C_c^infty(I). On those, Q_chi is Weil's form of
  sec. 2 and equals the zero sum.
Monotonicity in c: a function supported in I_{c'} (c' < c) is supported in I_c and Q_chi is the same functional, so a bound at c
holds for every c' <= c.

## 4. Sectors and node systems (Chekmachev Thm 3.2)
Real even f: F(r) = int f cos(rx). Even nodes: phi_j = N_j cos(2 pi j x/L), j >= 0 (N_0 = L^{-1/2}, N_j = (2/L)^{1/2}), orthonormal basis
of even L^2(I). Odd nodes: e_k = cos(pi k x/L), k odd >= 1, orthogonal basis, ||e_k||^2 = L/2. (LAT) becomes
  Q(f) = (1/2) sum_{j>=0} mu(2j) <f, phi_j>^2 + (1/L) sum_{k odd} mu(k) <f, e_k>^2,   (1/L) sum_{k odd} <f, e_k>^2 = ||f||^2/2.
Real odd f: phi_j = (2/L)^{1/2} sin(2 pi j x/L) (j >= 1), e_k = sin(pi k x/L); the same formula.

## 5. Minorant and dual check (Chekmachev Thm 4.1, Lemma 4.2, Thm 5.1, Prop. 5.2, Lemma 5.3, Thm 5.10)
Choose finite K (odd k) and E (j) with mu(k) >= 1 for odd k not in K and mu(2j) >= 0 for j not in E. Then, for eps < 1/2,
  Q(f) - eps ||f||^2 >= sum_j dl_j <f, phi_j>^2 + sum_{k in K} w_k <f, e_k>^2,  dl_j = (mu(2j) + 1)/2 - eps,  w_k = (mu(k) - 1)/L,
and dl_j >= deff := 1/2 - eps for j not in E. Write the right side as P - N with
  P = deff I + sum_{j in E, dl_j > 0} (dl_j - deff) phi_j phi_j^T + sum_{k in K, w_k > 0} w_k e_k e_k^T  (> 0),
  N = sum_{j in E, dl_j <= 0} (deff - dl_j) phi_j phi_j^T + sum_{k in K, w_k <= 0} (-w_k) e_k e_k^T = B B^T.
P - N >= 0 iff lambda_max(B^T P^{-1} B) <= 1 (B^T P^{-1} B by Woodbury; all inner products in closed form, the sums over j by
partial fractions, cert_chi.py). margin := 1 - lambda_max > 0 (rigorous) proves Q_chi >= eps on L^2(I).
Scan: for k > kscan, mu(k) >= monotone(k) := h(r_k) - 2 sum |chi(m)| Lambda(m)/sqrt m - 2 sum c^{-b_n} b_n/(b_n^2 + r_k^2) - (remainder),
since h is increasing in r > 0 (d/dy Re psi(x + iy) = -Im psi'(x+iy) > 0) and b/(b^2 + r^2) decreases; kscan is the first k with
monotone(k+1) >= 1 (odd) or >= 0 (even), found in Arb. Nodes below kscan:
- rigorous scan (rscan.c, proof in RIGOR.md): every node gets mu_c(k) with mu(k) >= mu_c(k) - E, E <= 1e-7, computed with integer phases,
  an Arb-built cosine table and IEEE +, -, * only. Every k with mu_c(k) < threshold + 1e-6 is re-evaluated in Arb (96 bits), and selected
  if its Arb lower endpoint is below the threshold. So K and E provably contain every node below the thresholds.
- float scan (fscan/cscan, the first runs): the same flags from float64 with libm cos and an estimated error budget (not a proof).
Selected nodes and Chekmachev's neighbour rule (Def. 6.1) define K and E.
Rigorous positivity (Chekmachev's Lemma 5.5): take an approximate smallest eigenvalue lambda of the midpoint of Id - W, a shift s0 < lambda
and nu < s0. Compute a float Cholesky factor L of mid(Id - W) - s0 Id, enclose the residual E = (Id - W) - L L^T - s0 Id in Arb, and check
max row sum |E| <= s0 - nu. Then lambda_min(Id - W) >= s0 - ||E||_2 >= nu > 0 (E is symmetric, so ||E||_2 <= max row sum).
Early runs: lambda from 4 inverse iterations, s0 = 0.995 lambda, nu = 0.99 lambda. Final code: lambda = smallest Rayleigh-Ritz value of
6 block inverse iterations, s0 = lambda/2 (then lambda/8, lambda/100 if the Cholesky fails), nu = s0/2. Any nu > 0 proves positivity.

## 6. Sign and convention checks (self-review)
- chi(n) vs conj chi(n): the n^{-s} terms of -L'/L(s, chi) carry chi(n) and pair with g(log n); the reflected line carries conj chi(n)
  with g(-log n). For real chi both are chi(n) and the term is 2 chi(n) g(log n) (g even for real f). Verified numerically (sec. 1).
- Parity: a = 1 for odd chi (chi_-4, chi_-3): Gamma((s+1)/2); a = 0 for even chi (chi_5, chi_8). Verified via the zero sums.
- Factor 2 between f and g = f * f~: (1/2 pi) int |F|^2 = ||f||^2 = g(0); sum_k |F(r_k)|^2 = 2L ||f||^2.
- Fourier sign: Chekmachev uses f^(r) = int f e^{-irx}; ours F(r) = f^(-r). mu is even in r for real chi, so nothing changes.

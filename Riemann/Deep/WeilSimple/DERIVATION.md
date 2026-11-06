# DERIVATION: simple, even ground state of the Weil form QW_lambda (weilsimple, 2026-11-05)

Notation: c = lambda^2, L = ln c = 2 log lambda, I = [-L/2, L/2], w_j = 2 pi j/L, b_k = 2k + 1/2, w_n = Lambda(n)/sqrt n.

## 1. QW_lambda is the Weil form Q of ../WeilZeta on I
CCM (arXiv 2511.22755, eq. (3.10) "quadratsemi"): on L^2([lambda^-1, lambda], d*u),
  QW_lambda(f, f) = int |f^(t)|^2 2 theta'(t) dt/(2 pi) + 2 Re(f^(i/2) conj f^(-i/2)) - sum_{1<n<=lambda^2} Lambda(n) <f | T(n) f>,
  f^(t) = int f(u) u^{-it} d*u, 2 theta'(t) = -log pi + Re psi(1/4 + it/2), <f|T(n)f> = n^{-1/2}((f^* * f)(n) + (f^* * f)(1/n)).
With u = e^x (d*u = dx) this is, term by term, Q(f) = 2|<f,cosh(x/2)>|^2 - 2|<f,sinh(x/2)>|^2 + (1/2 pi) int phi |F|^2
- sum_n w_n (C_f(log n) + C_f(-log n)) of ../WeilZeta/DERIVATION.md sec. 2 on L^2(I) (F(t) = f^(-t); the pole term
2 Re(F(i/2) conj F(-i/2)) = 2|<f,cosh>|^2 - 2|<f,sinh>|^2; CCM's W_R = -W_infty, eq. (3.7)). The atom n = lambda^2 sits at the
endpoint of supp(f * f^*) and contributes 0. So lambda = sqrt c, and "even" (gamma xi = xi, u -> 1/u) is x -> -x.
CCM's basis: U_n = L^{-1/2} e^{2 pi i n x/L} on [0, L], V_n(u) = U_n(log(lambda u)) = (-1)^n L^{-1/2} e^{i w_n x}. Hence the even
part of E_N = span{V_n : |n| <= N} is span{phi_j = N_j cos(w_j x) : 0 <= j <= N} (N_0 = L^{-1/2}, N_j = (2/L)^{1/2}) and the odd
part is span{psi_j = (2/L)^{1/2} sin(w_j x) : 1 <= j <= N}; QW_lambda^N is the direct sum of the two blocks.

## 2. Closed-form matrix entries (weil_matrix.py)
Bilinear form for real f1, f2 of equal parity, C12(x) = int f1(y) f2(y - x) dy (= C12(-x) for equal parity):
  Q(f1,f2) = +-2<f1,v><f2,v> + A(f1,f2) - 2 sum_{n<c} w_n C12(log n),   v = cosh(x/2) (even, +) or sinh(x/2) (odd, -).
Archimedean part: psi(z) = -gamma + int_0^inf (e^{-u} - e^{-zu})/(1 - e^{-u}) du with u = 2x and (1/2pi) int F1 F2 cos(tx) dt = C12(x):
  A(f1,f2) = -(log pi + gamma) <f1,f2> + int_0^inf [2 e^{-x} <f1,f2> - 2 e^{x/2} C12(x)]/(e^x - e^{-x}) dx.
(Check: minus CCM's W_R of eq. (3.5); the difference of the two expressions is (2 log 2 - log 4) C(0) = 0.)
Autocorrelations on 0 <= x <= L (elementary; s = (-1)^{i+j}):
  even, i != j: N_i N_j s (w_j sin w_j x - w_i sin w_i x)/(w_i^2 - w_j^2); i = j >= 1: (1 - x/L) cos w x - sin(w x)/(w L); i = j = 0: 1 - x/L;
  odd,  i != j: (2/L) s (w_i sin w_j x - w_j sin w_i x)/(w_i^2 - w_j^2);   i = j: (1 - x/L) cos w x + sin(w x)/(w L).
Pole vectors: <phi_j, cosh(x/2)> = N_j (-1)^j sinh(L/4)/(w_j^2 + 1/4); <psi_j, sinh(x/2)> = -2 N_j (-1)^j w_j sinh(L/4)/(w_j^2 + 1/4).
Archimedean integrals: expand 1/(e^x - e^{-x}) = sum_k e^{-(2k+1)x}; e^{i w_j L} = 1 makes every boundary term a multiple of
c^{-b_k}; the k-sums are digamma values (psi(z) = -gamma + sum [1/(k+1) - 1/(k+z)], sum 1/(k+z)^2 = psi'(z)), z_j = 1/4 + i w_j/2:
  S_j := int_0^L e^{x/2} sin(w_j x)/(e^x - e^{-x}) dx = (1/2) Im psi(z_j) - sum_k c^{-b_k} w_j/(b_k^2 + w_j^2);
  A_ij (i != j) = -2 N_i N_j s (w_j S_j - w_i S_i)/(w_i^2 - w_j^2) (even),  -2 (2/L) s (w_i S_j - w_j S_i)/(w_i^2 - w_j^2) (odd);
  A_00 = -log pi + psi(1/4) + psi'(1/4)/(2L) - (2/L) sum_k c^{-b_k}/b_k^2;
  A_jj (even, j >= 1) = -log pi + Re psi(z_j) + (Im psi(z_j)/w_j + Re psi'(z_j)/2)/L - (4/L) sum_k c^{-b_k} b_k^2/(b_k^2 + w_j^2)^2;
  A_jj (odd)          = -log pi + Re psi(z_j) - (Im psi(z_j)/w_j - Re psi'(z_j)/2)/L + (4 w_j^2/L) sum_k c^{-b_k}/(b_k^2 + w_j^2)^2.
(Per k: int_0^L e^{-bx}(1 - x/L) e^{iwx} dx = 1/z - (1 - c^{-b})/(L z^2), z = b - iw; int_0^L e^{-bx} sin(wx) dx = (1 - c^{-b}) w/(b^2 + w^2).)
As L -> infinity, A_jj -> phi(w_j) (the symbol). Rigor: Arb digamma/Hurwitz zeta; the c^{-b_k} sums are cut when c^{-b_K} < 2^-(prec+40)
and every term is <= 16 c^{-b_k} (b_k >= 1/2, w/(b^2+w^2) <= 1/(2b), 4b^2/(b^2+w^2)^2 <= 4/b^2, 4w^2/(b^2+w^2)^2 <= 1/b^2), so the tail
is enclosed by a ball of radius 16 c^{-b_K}/(1 - c^{-2}).
Checks (test_matrix.py, logs/test_entries_c5.log): every entry with i, j <= 2 of both blocks at c = 5 equals an independent mpmath
computation (x-quadrature for pole and primes, t-quadrature of (1/2pi) int phi F_i F_j for the archimedean part, F_j checked against
quadrature) to 1e-7 (the quadrature's own truncation error); six entries with indices up to 5 at c = 25/8 to 5e-6 (test_entries_c25o8.log);
ten entries at c = 5 and 25/8 equal a 50-digit quadrature of the x-space formula above to 3.5e-49 (test_hp.log). Spectrum check: our U(c) >= weilzeta's certified lower bound of the even
ground state (lattice minorant, a code that shares nothing with these closed forms): c = 5: 1.2218e-17 vs 7.447e-18; c = 10:
1.9445e-43 vs 9.134e-44; c = 16: 3.7430e-75 vs 1.5128e-75.

## 3. The criterion (Courant-Fischer)
A_lambda (CCM Thm 3.6 "thmsmallest") is self-adjoint, bounded below, with discrete spectrum, and commutes with parity. Let
lambda_1^e <= lambda_2^e <= ... and lambda_1^o <= ... be its eigenvalues on even and odd functions (min-max values of Q on L^2(I)).
  (a) U := Q(g)/||g||^2 for an explicit even g (finite cosine series)            =>  lambda_1^e <= U;
  (b) Q(f) >= L_odd ||f||^2 for all odd f in L^2(I)                                =>  lambda_1^o >= L_odd;
  (c) Q(f) - L_2 ||f||^2 + M <f,g>^2 >= 0 for all even f (some M >= 0)            =>  Q(f) >= L_2 ||f||^2 for even f orthogonal to g
                                                                                        =>  lambda_2^e >= L_2 (Courant-Fischer).
Complex f = u + iv (u, v real): Q(f) = Q(u) + Q(v) and |<f,g>|^2 = <u,g>^2 + <v,g>^2 (g real), so (b), (c) hold for complex f.
If U < L_odd and U < L_2, the smallest eigenvalue mu_lambda of A_lambda equals lambda_1^e, has multiplicity one, and its eigenvector is
even: the spectral subspace of (-inf, min(L_odd, L_2)) has no odd vector and its even part has dimension <= 1 (a 2-dim even subspace
contains a vector orthogonal to g) and >= 1 (U < L).
Truncations: if g lies in the even part of E_NG then for every N >= NG the matrix QW_lambda^N is even-simple (CCM Def. 5.3):
its odd block is >= L_odd (compression of (b)); its even block has lowest eigenvalue <= Q(g)/||g||^2 = U (g in E_N) and second
eigenvalue >= L_2 (an even f in E_N orthogonal to g is orthogonal to g in L^2, so (c) applies). Hence CCM Thm 5.10 applies to
QW_lambda^N for all N >= NG, with eps_N in [mu_lambda, U].

## 4. The deflated certificate (cert_simple.py, sector defl)
cert_zeta.py's even certificate proves Q(f) - eps ||f||^2 >= (P - N)(f) >= 0 via the minorant P - N and Woodbury
(../WeilZeta/DERIVATION.md sec. 5). Adding M g g^T to the positive part: Q(f) - eps ||f||^2 + M <f,g>^2 >= (P + M g g^T - N)(f), and
P' - N >= 0 iff Id - B^T P'^{-1} B >= 0, the same Woodbury reduction with one more column u = g, weight M (here M = 1), in U and D.
Inner products of g = sum_l g_l phi_l (l <= NG <= K2, so every l is in the node set E): <phi_j, g> = g_j, <e_k, g> = sum_l g_l <phi_l, e_k>,
<v, g> = sum_l g_l <phi_l, v> (closed forms of sec. 4 there), and since P0 phi_l = p_l phi_l (p_l = dl_l on G, deff otherwise):
<g, P0^{-1} b> = sum_l (g_l/p_l) <phi_l, b>, <g, P0^{-1} g> = sum_l g_l^2/p_l. Everything else (scan, minorant, Lemma 5.3 bound,
Cholesky positivity proof) is cert_zeta.py unchanged; it reproduces weilzeta's c = 5 margins exactly (even 4.00280318265e-19, odd
8.28212749750e-16). g = inverse-iteration ground vector of mid(even block of QW_lambda^NG), rounded to exact dyadics; U is enclosed
in Arb from the Arb matrix (sec. 2). K2 (all nodes k <= K2 kept) = max(401, ceil(5c^2/3)): with cert_zeta's default 201 the odd
minorant at c = 7 and c = 9 has a negative eigenvalue independent of eps (logs/int_k201_fail.jsonl); 401 removes it. (The integer
rows c = 2-6, 8, 10 ran with 201 and succeeded; any K2 gives a valid minorant, a larger one a tighter one.)

## 5. All windows in an interval
For c < c' the space L^2(I_c) sits isometrically in L^2(I_c') (extension by 0) and Q_{c'} = Q_c there (only n < c can meet
supp(f * f^*)), parity is preserved, so by min-max every lambda_k^e(c), lambda_k^o(c) is non-increasing in c. Therefore if
U(c_i) < L_odd(c_{i+1}) and U(c_i) < L_2(c_{i+1}), then for every c in [c_i, c_{i+1}]:
lambda_1^e(c) <= lambda_1^e(c_i) <= U(c_i) < L(c_{i+1}) <= lambda_1^o(c), lambda_2^e(c). chain.py builds such a chain on [2, 20]
(37 steps: 30 from 2 to 16, 7 more to 20; logs/chain.jsonl, all ok; table.py rechecks every step; eps at c_{i+1} = 2.00001 U(c_i)). Between grid points the truncation statement holds for N >= some N_0(c) (E is a form core, CCM
Prop. 3.3, so the lowest even truncated eigenvalue tends to lambda_1^e(c) < L), not with an explicit N_0.

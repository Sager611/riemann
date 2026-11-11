# DERIVATION (aclgfals): Galerkin models used for the tests

## 1. Semilocal pair as a Hankel matrix (statements 4, 5)
ACLG l. 1097-1113: on H = L^2(R_+^*, d*u), Mg(s) = int g u^s d*u, ||g||^2 = (1/2pi) int |Mg(i tau)|^2 dtau; P_t = 1_(0,t];
F_S = I o U (I: u -> 1/u; U = Mellin multiplier by u(s) = prod_{v in S} u_v(s)); hat P_t = U^*(1 - P_{1/t})U (eq. hankel l. 1173).
- Local factors. u_infty(s) = Gamma_R(1/2 - s)/Gamma_R(1/2 + s) (l. 1111), = exp(-2 i vartheta(tau)) at s = i tau (vartheta = Riemann-Siegel
  theta; ACLG l. 1111). Tate: int hat f |y|^{1-w} d*y = rho_v(w) int f |x|^w d*x with rho_v(w) = L_v(1-w)/L_v(w); the archimedean case of
  this, at w = 1/2 + s, is ACLG's u_infty, hence u_p(s) = rho_p(1/2 + s) = (1 - p^{-1/2-s})/(1 - p^{-1/2+s}): unimodular on the line,
  2 pi i/log p-periodic, Laurent coefficients c_{-1} = -p^{-1/2}, c_m = (1 - 1/p) p^{-m/2} (m >= 0), = O(p^{-|k|/2}) as ACLG l. 1240 says.
- Cayley. z = (tau - i)/(tau + i) = e^{i theta}, tau = -cot(theta/2), dtheta = 2 dtau/(1 + tau^2). F_k(tau) = pi^{-1/2} z^k/(tau + i) is an
  ONB of L^2(R, dtau); k >= 0: analytic in Im tau > 0; k <= -1: analytic in Im tau < 0. For g supported in (0,1], Mg(i tau) =
  int_{x<=0} G(x) e^{i tau x} dx is analytic in Im tau < 0; so Range P = span{F_k : k <= -1}, Range(1 - P) = span{F_k : k >= 0}.
  Multiplication by phi(tau) has matrix phi_{k-l}, phi_m = (1/2pi) int phi(e^{i theta}) e^{-i m theta} dtheta; phi(-tau) = conj phi(tau)
  (true for every u_v) makes all phi_m real.
- Gamma_1 = (1 - P) M_u P has matrix <phi F_{-l-1}, F_j> = phi_{j+l+1} (j, l >= 0): a Hankel matrix. P_t = M_{t^{i tau}} P M_{t^{-i tau}}
  gives Gamma_t = M_{t^{-i tau}} [(1 - P) M_{phi_t} P] M_{t^{-i tau}}, phi_t = t^{2 i tau} u (ACLG's phi_t of l. 1182). Hence
  kappa_n(t) = singular values of (phi_{t, j+l+1}). Check: S = {infinity} gives |chi_n(t)| (prolate, gamma = 2 pi t^2), n even.
- Pair vectors. Gamma e = kappa g, f = U^* g, a = (e + f)/sqrt(2 + 2 kappa) (<e, f> = <Gamma e, g> = kappa). In the conjugated frame
  |M a(i tau)|^2 = |e~ + conj(phi_t) g~|^2/(2 + 2 kappa) (the multipliers t^{+-i tau} are unimodular), so for a real cut-off chi(tau),
  ||chi^{1/2} a||^2 = (v^T T_chi v + u^T T_chi u + 2 v^T X u)/(2 + 2 kappa), T_chi = (chi_{a-b}) Toeplitz, X_{lj} = (chi phi_t)_{l+j+1},
  using |phi_t| = 1 and <chi conj(phi) g, e> = sum v_l u_j (chi conj phi)_{-(l+1+j)} = sum v_l u_j (chi phi)_{l+j+1} (real coefficients).
- Truncation. H_M = (phi_{j+l+1})_{j,l < M}: the compression of Gamma to the first M Laguerre functions on each half-line (F_k <-> Laguerre
  functions in x = log u). Finite sections of a compact Hankel converge in norm; singular values of H_M increase with M (interlacing of
  compressions), so kappa_n^(M) <= kappa_n. Tested at two M.

## 2. Gram constant (statement 6)
Prop. xipoly (l. 506-510) and eq. mellinE (l. 422): for h in H_K^Herm (int h = h(0) = 0), M(E h)(s) = zeta(s + 1/2) Mh(s + 1/2), so
||E h||^2 = (1/2pi) int |zeta(1/2 + it)|^2 |Mh(1/2 + it)|^2 dt, while ||h||^2_{L^2(R)} = (2/2pi) int |Mh(1/2 + it)|^2 dt. In the
orthonormal Hermite basis the Gram matrix of ||.||_2 is the identity; g_K^Herm = lambda_min of the |zeta|^2-weighted matrix on the
constraint hyperplane. Hence g_K^Herm = (1/2) min over the K-dimensional family of probability densities |Mh|^2 / int |Mh|^2 of the
mean of |zeta(1/2 + it)|^2: small iff that family can concentrate near zeros of zeta.
Finite lambda: E(f) restricted to [1/lambda, lambda], f in H_K(lambda), prolate h_{4k,lambda} orthonormal on [-lambda, lambda].

## 3. Nystrom model of the semilocal pair at t = 1 (statements 4, 5; code/semi_pair.py)
- S = {infinity}: F_infty on even f is cos-transform; on g(u) = u^{1/2} f(u) in L^2(d*u): (F g)(v) = int 2 (uv)^{1/2} cos(2 pi u v) g(u) d*u;
  U = I o F: (U g)(v) = int k(u/v) g(u) d*u, k(w) = 2 w^{1/2} cos(2 pi w) (a convolution; Mellin multiplier u_infty).
- U_p = sum_m c_m theta(p^m) (theta(v) g(u) = g(u/v), Mellin multiplier v^s; c_m from sec. 1), so U = U_infty prod_p U_p has kernel
  sum_m c_m k(p^m u/v) (one prime; products of the Laurent series for several).
- Gamma = (1 - P) U P with P = 1_(0,1]: for g = u^{1/2} G on (0,1] and w = 1/v in (0,1),
  (Gamma g)(1/w) = w^{1/2} (A G)(w),  (A G)(w) = 2 int_0^1 C_S(u w) G(u) du,  C_S(z) = sum_n gamma(n) cos(2 pi n z),
  n = prod_{p in S} p^{k_p} (k_p >= -1), gamma(n) = prod_p g_p(k_p), g_p(k) = 1 - 1/p (k >= 0), -1/p (k = -1)
  (c_m p^{m/2} = 1 - 1/p for m >= 0 and -1/p for m = -1). ||Gamma g|| = ||A G||_{L^2(0,1)}, ||g|| = ||G||: Gamma is unitarily
  equivalent to the symmetric operator A on L^2(0,1); kappa_n = |eigenvalues of A|. S = {infinity}: A = finite cosine transform with
  gamma = 2 pi, kappa_n = |chi_n(1)| (checked to 1e-9).
- For S = {infinity, 2} the series is lacunary (frequencies 2^k, weights 1/2, plus -cos(pi z)/2): each term 2^k adds about 2^{k+1}
  singular values of size 2^{-k/2} (the "p^k values of size p^{-k/2}" of ACLG l. 1240). Cut: n <= nmax (then U is replaced by a
  non-unitary U_K; f_n = U_K^* g_n is used consistently). Observed exact self-similarity kappa_{2n}(2 nmax) sqrt(2n) = kappa_n(nmax) sqrt(n)
  (0.442 at (50, 64), (100, 128), (200, 256), (400, 512)), so kappa_n sqrt(n) -> C (about 0.85-0.95) as nmax -> infinity: order 1/2.
- Mellin data: M e(i tau) = int_0^1 G(u) u^{-1/2 + i tau} du, M g_out(i tau) = int_0^1 H(w) w^{-1/2 - i tau} dw (A G = kappa H),
  M f = conj(u_K(i tau)) M g_out, a = (e + f)/sqrt(2 + 2 kappa); m_n(R) = (1/2pi) int chi_R |M a_n|^2 dtau (tau-symmetric, Gauss panels).
  Product integration: G interpolated (degree 19 per panel in s = sqrt u), u^{-1/2 + i tau} du = 2 s^{2 i tau} ds integrated on
  16-point Gauss panels in log s. Discretisation: composite Gauss-Legendre in s, panels <= 3 wavelengths of cos(2 pi nmax s^2);
  npw = 5 (float32) and npw = 7 (float64) agree to 4 digits in every window sum at nmax = 256.
- Near-degenerate couples (kappa_{2k+1} ~ kappa_{2k+2}) exchange their vectors between runs: window boundaries at even n.

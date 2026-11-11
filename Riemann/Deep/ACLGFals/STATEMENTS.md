# ACLG 2611.13129: the six statements, verbatim content, classified (aclgfals, 2026-11-11)

Source: arXiv e-print of 2611.13129 (Andrews, Connes, La Gatta), file `prolate-weil-fine15.tex` (1527 lines; tarball dated
11 Nov 2026 01:14). Line numbers below are of this file. (R419's l. 8488-9225 refer to another line breaking; R432's l. 1255-1369
match this file.)

## Common definitions
- QW(g,h) = Psi(F), F(x) = int g(u)h(ux) d*u, Psi = W_{0,2} - W_R - sum_p W_p (l. 403-411, eqs. Psi, Wp, WR; CCM normalisation).
  QW_lambda = restriction to L^2([1/lambda, lambda], d*u); commutes with iota(u) = 1/u; eigenvalues eps_lambda(1) <= eps_lambda(2) <= ...
  (even sector), eps^odd_lambda(1) <= ... (odd); eps_lambda = eps_lambda(1) (l. 412-413). Primes enter for p^m <= lambda^2.
- Prolate side (l. 427-446): gamma = 2 pi lambda^2; h_{n,lambda} (n even) orthonormal prolate functions on [-lambda, lambda];
  hat h_{n,lambda} = chi_n(lambda) h_{n,lambda} on [-lambda, lambda]; delta_n(lambda) := 1 - |chi_n(lambda)| (EXACT rung, eq. chi l. 436);
  Fuchs: delta_n ~ (2 sqrt pi 8^n (2pi)^{n+1/2}/n!) lambda^{2n+1} e^{-4 pi lambda^2} (eq. fuchsintro, l. 377; l. 468-481).
- E(f)(u) = u^{1/2} sum_{n>=1} f(nu) (l. 356); H_K(lambda) = {sum_{k<=K} c_k h_{4k,lambda} : int h = 0}, dim K; Hermite limit H_K^Herm
  (l. 504); Gram constant g_K^Herm := inf{||E(h)||^2 : h in H_K^Herm, ||h||_2 = 1} > 0 with
  ||E(h)||^2 = (1/2pi) int Xi(t)^2 |q_h(-t^2)|^2 dt (Prop. xipoly, l. 506-510, eq. gramherm l. 508).
- Semilocal pair (l. 1097-1113, 1170-1222): S finite set of places, {infinity} and all p <= lambda^2 in S; H = L^2(X_S)^{K_S} = L^2(R_+^*, d*u);
  theta(f) g = f * g (multiplicative convolution; Mellin multiplier by Mf); P_t = multiplication by 1_(0,t]; F_S = I o U, I(u) = 1/u,
  U = Mellin multiplier by u(s) = prod_{v in S} u_v(s), u_infty(s) = Gamma_R(1/2 - s)/Gamma_R(1/2 + s) (l. 1111), u_v = ratio of local
  factors; hat P_t = F_S P_t F_S^{-1} = U^*(1 - P_{1/t})U (eqs. IUI l. 1107, hankel l. 1173). kappa_n(t) = singular values of
  Gamma = (1 - P_{1/t}) U P_t (Lemma H2, l. 1191-1199), Gamma e_n = kappa_n g_n, f_n = U^* g_n, a_n ~ e_n + f_n (eigenvector of
  1 - P_t - hat P_t for -kappa_n; eqs. Kn l. 1214, decompK l. 1218, l. 1239).
  N_t^(n0) := sum_{n > n0} kappa_n(t) |a_n><a_n| (l. 1221); Lambda(f) := (1/2pi) int |Mf(is)|^2 log(2 + |s|) ds (eq. Lambda l. 1222);
  rho^(n0)_N(s) = sum_{n>n0} kappa_n |M a_n(is)|^2 (l. 1252, 1259).

## Statements

| # | Statement (exact content) | TeX lines | Class | Test here |
|---|---|---|---|---|
| 1 | **Hypothesis LA.** There are 0 < c <= C (and lambda_0, J) with, for all lambda >= lambda_0 and 1 <= j <= J: c delta_{4j}(lambda) <= eps_lambda(j) <= C delta_{4j}(lambda) log lambda and c delta_{4j+2}(lambda) <= eps^odd_lambda(j) <= C delta_{4j+2}(lambda) log lambda. Upper halves = Thm upperK with (log lambda)^3; the lower halves are the hypothesis ("a statement about constants, not only exponents"; sharp test: log(eps(2)/eps(1)) vs log(delta_8/delta_4), log(eps^odd(1)/eps(1)) vs log(delta_6/delta_4), eq. gaptest). | 588-599 (eq. LA 590, gaptest 596) | (iii) as stated (existence of constants for lambda >= lambda_0, J free: no finite computation can violate it); (ii) in its operative form: the ratios eps(j)/delta_{4j} and the log-gap offsets must not drift with lambda | ratios, offsets, interlacing, counts at lambda^2 = 2...20 (below the paper's 5 <= lambda^2 <= 300); certified one-sided bounds |
| 2 | **Positivity.** QW_lambda >= 0, i.e. eps_lambda >= 0 (Weil positivity on the window), assumed in all of Sec. 6; for all lambda it is RH (Prop. monotone; Remark "on the logic": cannot be used in an unconditional proof). | 898; 904-915; 1035-1037; 1357-1363 | (i) at each lambda (a sign); for all lambda = RH | falsifier = any Galerkin eigenvalue < 0 (min-max makes it a certified counterexample); certified odd-sector positivity |
| 3 | **Gap (Lemma gap).** Assume QW_lambda >= 0 and every eigenvalue other than eps_lambda is >= g_lambda > 0. Then ||omega_lambda||^2 <= C delta_4 (log lambda)^3/g_lambda, ||k_lambda - eta/||eta|| || <= ||omega||(1+||omega||^2)^{1/2}; if ||omega_lambda|| < 1/2, simple and iota-even. Versions: full gap g_lambda = c delta_6(lambda) (lower half of LA for eps^odd(1), eps(2)); even gap g_lambda = c delta_8(lambda) (second even eigenvalue). | 904-920 (eq. omega 906) | (ii): fixed c, uniform in lambda; at each lambda a numerical inequality | g/delta_6 and eps(2)/delta_8 against lambda; certified lower bounds of weilsimple (L_odd, L_2) |
| 4 | **(H2').** For f real, smooth, supported in [1/lambda, lambda]: Tr(theta(f) N_t^(n0) theta(f)^*) = sum_{n>n0} kappa_n(t) ||theta(f) a_n||^2 <= eps_t(n0) Lambda(f), eps_t(n0) -> 0 (n0 -> infinity), uniformly in f. Open for the semilocal pair (l. 1240); equivalent (Remark H2reduced) to the density bound rho^(n0)_N(s) <= eps_t(n0) log(2+|s|) together with (5). Predicted mechanism: eps_t(n0) ~ |S| kappa_{n0}(t), strong form with ||f||^2 (l. 1240). | 1221-1240 | (iii): a limit with no rate over all f; the predicted size |S| kappa_{n0} and the density bound are (ii) | density rho^(n0)_N(s) at bounded s on Hankel truncations |
| 5 | **Localisation (H2ii).** Tr(M_{chi_R} N_t^(n0)) = sum_{n>n0} kappa_n(t) ||chi_R^{1/2} a_n||^2 <= eps_2(n0, R) -> 0 (n0 -> infinity), polynomially in R; chi_R in C_c^infty, 0 <= chi_R <= 1, = 1 on [-R, R], supp in [-2R, 2R], M_chi = Mellin multiplier. Proved halves: <= C R^A for every n0 (Lemma localization, eq. localmass l. 1244); -> 0 for fixed R without rate (Lemma H2q). "The precise statement to establish from u = u_infty prod_p u_p"; gives (H2') on the smoothed Galerkin spaces theta_w E_N. | 1242-1262 (eq. H2ii 1257); 1392-1394 item 1 | (ii): a rate in n0 with polynomial R-dependence, computable for given S, t | T(n0, R) on Galerkin (Laguerre = Hankel-matrix) truncations of the pair, S = {inf}, {inf,2}, {inf,2,3}; fitted rates |
| 6 | **(R2).** For k in V_low (conditioning vectors attached to modes with delta_m(lambda sigma) <= tau = lambda^{-A}): k = K - r, K in Rad(QW), ||r|| <= C tau^{1/2} (log lambda)^C ||k||; hence |QW_lambda(k, v)| <= C tau^{1/2} (log lambda)^C ||k|| (QW_lambda(v,v) + C_0 lambda ||v||^2)^{1/2} and |QW_lambda(k, k')| <= C tau (log lambda)^C ||k|| ||k'||. "Would follow from an extension of [CK] uniform in the index up to tau, together with the Gram lower bound of Prop. xipoly beyond the Hermite regime; we leave this open." Needed in Prop. reduction with alpha, beta = O(tau^{1/2} lambda^{C'}), the power of lambda coming partly "from the Gram constant". | 1364-1373; 1376-1378; 1395 item 2 | (ii): uniform in the index; the Gram input g_K must stay >= (log lambda)^{-C} (or lambda^{-C'}) up to K ~ lambda^2 | g_K^Herm for K <= 40 from eq. gramherm; g_K(lambda) for K <= lambda^2 - 1 at small lambda; ratio QW_lambda(E f, E f)/(delta_{4K} ||E f||^2) |

## How they combine (for the report)
- Sec. 6: (2) + (3) (+ upper bound of Kastoris) => ground state within O(lambda^{-4}(log lambda)^{3/2}) of k_lambda, distance to the target
  tau = theta/theta(1) O(lambda^{-2}) in ||.||_1; Prop. hurwitz (l. 1026-1033): ||.||_1-convergence to tau + real zeros (Thm realzeros)
  => RH. Circular by (2) (their Remark l. 1035).
- Sec. 8: (4) => Thm condpos (QW_lambda >= (log sigma)/2 on an explicit finite-codimension space); (4) + R1 (proved) + (6) + the sign of a
  Schur complement on the transition-and-pole block => eps_lambda >= -eta(lambda) => (Prop. monotone) RH (Prop. reduction, l. 1376).
- Only (2) at a fixed lambda and the finite-lambda instances of (1), (3) are single-number checks; (4) is not; (5) and (6) are rate
  statements: they can fail visibly on a range, which is what this audit tests.

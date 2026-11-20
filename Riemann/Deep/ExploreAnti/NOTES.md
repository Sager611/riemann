# explore55 working notes (anti-concentration -> N8)

## Choice (02:05 UTC)
- Taken routes (do not repeat): triples/RS n-level (explore50/52, npafals); T^inf concentration (ExploreConcentration); Nevanlinna (explore54).
- Tool: second-moment small-ball bound (Paley–Zygmund / Christoffel–Markov–Stieltjes) in its REAL 2D form:
  for a random vector W in R^d, P(W != 0) >= m^T E[WW^T]^{-1} m, m = E W (sharp given first and second moments).
  CGG (Conrey–Ghosh–Gonek 1998, Bui–Heath-Brown 1302.5018: 19/27 under RH) use the complex 1D form
  P(Y != 0) >= |EY|^2 / E|Y|^2 with Y = zeta'(rho) M(rho); this ignores the pseudo-covariance E[Y^2].
- Under RH: zeta'(rho) = -i e^{-i theta(gamma)} Z'(gamma), zeta'(rho)^2 = -chi(rho)|zeta'(rho)|^2,
  conj(Y) = zeta'(rho) * (-chi(1-rho) M(1-rho)); so the real 2D bound = CGG with the two-piece mollifier
  psi(s) = M1(s) + chi(1-s) M2(1-s) (Bui–Conrey–Young-type chi-twisted piece) and needs the holomorphic
  moment sum_gamma zeta'(rho)^2 M_i(rho) M_j(rho) (computable type: Landau–Gonek with chi, arXiv 2601.18025).
- Control: Milinovich–Ng 1106.1160 (J_{-1} >= (3/(2pi^3))T, half of Gonek's 3/pi^3): their factor 1/2 = theta/(1+theta) at theta->1
  (mollifier length), and the conjugate (negative-frequency) features have only O(T log T) cross terms: predicted no 2D gain there.

## Numerics, zeta zeros (Z' by Riemann-Siegel, validated vs mpmath: rel err <= 2e-3 at gamma_1, < 1e-5 beyond zero 1000; at 10^12 |Z(gamma)| <= 7.5e-7)
Y = zeta'(rho) M(rho), M = sum_{n<=y} mu(n) x^j n^{-s}, x = log(y/n)/log y, y = T^theta. CGG_1D = max complex-linear |sum c.Y|^2/(N sum|c.Y|^2);
2D = max real-linear (sum a.W)^2/(N sum (a.W)^2), W = (Re Y, Im Y). par = Var(parallel to mean)/Var(total), P = x.
| zeros | L | theta | CGG J=1 | 2D J=1 | CGG J=3 | 2D J=3 | mu/L | A/L^2 | ReB/L^2 | par |
|---|---|---|---|---|---|---|---|---|---|---|
| 50001-1e5 | 9.39 | 1/4 | 0.6105 | 0.8177 | 0.6493 | 0.8402 | 0.666 | 0.727 | 0.359 | 0.349 |
| 50001-1e5 | 9.39 | 1/2 | 0.7856 | 0.9275 | 0.8106 | 0.9348 | 0.805 | 0.826 | 0.573 | 0.286 |
| 1e12+1..+1e4 | 24.48 | 1/4 | 0.5299 | 0.7451 | 0.5478 | 0.7588 | 0.640 | 0.7725 | 0.326 | 0.386 |
CGG asymptotics (BHB 1302.5018 l.240): mu/L = 1/2 + theta*intP, A/L^2 = 1/3 + theta*intP + theta^2 (intP)^2 + intP'^2/(12 theta):
P = x, theta = 1/4: 0.625, 0.8073 (CGG 0.484). Data converge like c/L (mu: 0.666 -> 0.640, predicted 0.640 from c = 0.37).
Two-point 1/L extrapolation, theta = 1/4: ReB/L^2 -> 0.306, par -> 0.40, 2D(J=1) -> ~0.70 vs CGG 0.484 (heuristic, needs B's main term).

## Recipe (CFKRS ratios) main term of the holomorphic moment — derived 02:25 UTC
Sum_gamma G(rho) = (1/2pi i) contour of (zeta'/zeta) G; left edge: zeta'/zeta(s) = chi'/chi(s) - zeta'/zeta(1-s); M(s) = (1/2pi i) int y^z K(z)/zeta(s+z) dz,
K(z) = sum_j p_j j!/((log y)^j z^{j+1}). Only the swap (alpha <-> delta1) survives d_alpha d_beta; at leading order (zeta(1+x) ~ 1/x, arithmetic factor 1):
D(alpha,beta) = [beta^2 g(alpha) - alpha^2 g(beta)]/(alpha beta (beta-alpha)), g(x) = e^{-x l}(1 - x lam)^2, lam = theta L intP, l = log(t/2pi);
d_a d_b D|_0 = -g_3 = l^3/6 + l^2 lam + l lam^2. Hence (per zero, units L^2):
  B := E[zeta'(rho)^2 M(rho)^2] = 1/6 + theta I + theta^2 I^2 = u^2 - 1/12,   u := 1/2 + theta I = mu/L (BHB),  I = intP.
Checks: same recipe gives Sum zeta'(rho) ~ (T/4pi)L^2 and BHB's S_1 = (T L^2/2pi)u exactly. Numerics: B/L^2 extrapolated (1/L, L = 9.39, 24.48)
0.4758 vs 0.4792 (theta 1/2), 0.306 vs 0.3073 (theta 1/4); par fraction -> 1/3, 0.4 predicted, data 0.316, 0.386 (rising).
Consequence: Var(Y) = A - u^2 = 1/12 + J'/(12 theta) (J' = int P'^2), pseudo-variance B - u^2 = -1/12:
  Var(Re Y) = J'/(24 theta), Var(Im Y) = 1/12 + J'/(24 theta); the real projection removes the 1/12 and half of J'/(12 theta).
  kappa_2D = u^2/(u^2 + J'/(24 theta)); optimum over P (P(0)=0): Riesz representer R = (1/2+theta)x - theta x^2/2,
  sup u^2/J' = int R'^2 = (3 + 6 theta + 4 theta^2)/12 =: Q/12, so kappa_2D(theta) = 2 theta Q/(1 + 2 theta Q); theta -> 1/2: 7/8 (CGG/BHB: 19/27).

## Third height + CUE model (02:40 UTC)
zeta, P = x: B - u^2 = -0.0757, -0.0799, -0.0838 (theta 1/2; L = 9.39, 12.10, 24.48), -0.0854, -0.0864, -0.0830 (theta 1/4): -> -1/12 = -0.0833.
Var A - u^2 = 0.177, 0.203, 0.228 (theta 1/2, pred 0.25), 0.283, 0.323, 0.363 (theta 1/4, pred 0.417): slow, from below (so finite-T bounds overshoot).
CUE / function-field model (L(u) = prod(1 - a_j u), Haar U(N), mollifier from h_k(a), K = theta N; 60,000 zeros per row), theta = 1/2:
u = 0.7505-0.7529, A = 0.804-0.816, B = 0.462-0.485, B - u^2 = -0.082 +- 0.002 (N = 10..200), CGG(P=x) 0.68-0.70 (pred 0.6923), 2D 0.862-0.872 (pred 0.871), par 0.334-0.340 (pred 1/3);
theta = 1/4: u 0.598-0.628 (0.625), A 0.805-0.857 (0.807), B - u^2 = -0.079..-0.084, 2D 0.63-0.71 (0.702), par 0.400-0.416 (0.4). Recipe confirmed in the model.
Memory: fit.py on zeros6 hit 436 MB max RSS once (rs.zprime chunk 4000; memcap sampled 18 MB); fixed chunk = 500 (35 MB); see logs/MEMORY_NOTE.txt.
- 02:50 hi_cross (10^12, theta .35/.45): all 15 pairwise B_ij - u_i u_j in [-0.0865,-0.0828]: bilinear H supported.

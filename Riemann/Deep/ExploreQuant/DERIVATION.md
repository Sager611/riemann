# ExploreQuant — derivations (explore29)

Setting. I = [A, B], p >= 0 continuous on I, positive mass on every cell. MSE quantizer y_1 < ... < y_N in I, thresholds
m_n = (y_n + y_{n+1})/2, cells V_1 = [A, m_1], V_n = [m_{n-1}, m_n], V_N = [m_{N-1}, B]. Lloyd–Max necessary conditions (stationarity of
D(y) = sum_n int_{V_n} (x - y_n)^2 p): y_n = int_{V_n} x p / int_{V_n} p. Gaps g_n = y_{n+1} - y_n.

## L1 centroid–gap lemma
Interior cell V_n (2 <= n <= N-1) has midpoint (m_{n-1} + m_n)/2 = y_n + (g_n - g_{n-1})/4. If p is nondecreasing on V_n, Chebyshev's
integral inequality (x and p comonotone) gives centroid >= midpoint, strict unless p is a.e. constant on V_n. Stationarity
(centroid = y_n) then gives g_n <= g_{n-1} (equality iff p a.e. constant on V_n). Nonincreasing p: g_n >= g_{n-1}.

## L2 peaks need valleys
Strict peak at j (2 <= j <= N-2): g_{j-1} < g_j > g_{j+1}. By L1, p is not nondecreasing on V_j: x1 < x2 in V_j with p(x1) > p(x2); and
not nonincreasing on V_{j+1}: x3 < x4 in V_{j+1} with p(x3) < p(x4). Since x2 <= m_j <= x3, the minimum of p on [x1, x4] is attained at
some a in (x1, x4) with p(a) < min(p(x1), p(x4)): a dip inside V_j u V_{j+1} (not necessarily a strict local minimum: flat
bottoms are allowed). Peaks are >= 2 indices apart, so dips of distinct peaks lie in disjoint open intervals. Hence #strict peaks of
(g_n) <= #disjoint dips of p. A unimodal p (every log-concave p) has no dip: its stationary quantizers have quasi-convex gap sequences.

## L3 oscillation bound
Rescale V_n to [-1, 1]. The centroid must sit o_n = |g_n - g_{n-1}|/(g_n + g_{n-1}) half-lengths from the midpoint. Among densities with
values in [c, c e^eps], the largest centroid offset is attained by the bang-bang density c(1 + K 1_{x > 1-u}), K = e^eps - 1
(linear-fractional program), offset K u(2-u) / (2(2 + K u)), maximised over u in (0, 2). So
  osc_{V_n} log p >= eps(o_n), eps(o) = log(1 + K(o)), K(o) solving max_u K u(2-u)/(2(2+Ku)) = o.
Closed form (reviewer): the maximiser is u* = 2/(r+1), r = sqrt(1+K), with value (r-1)/(r+1); so r = (1+o)/(1-o) and
  eps(o) = 2 log((1+o)/(1-o)) = 4 artanh(o) = 2 |log(g_n/g_{n-1})|
(eps(0.1) = 0.401, eps(0.3) = 1.238, eps(0.5) = 2.197, eps(0.8) = 4.39, as the numerical table). The extremiser is discontinuous, so for
continuous p the bound is a supremum bound (still valid); p = 0 somewhere on V_n gives osc = infinity. Cells are disjoint, so
  TV_I(log p) >= 2 sum_n |log(g_n/g_{n-1})|.

## Z  Zador / Fejes Toth efficiency
Unfold x_n = theta(gamma_n)/pi + 1 (uniform density). Voronoi MSE = sum_n g_n^3/12 + O(boundary); the N-point minimum on the same
interval is N (L/N)^3/12 (Fejes Toth; Zador's constant 1/12 in d = 1). Ratio kappa_2 = <delta^3>/<delta>^3. In the t-variable the
same ratio arises for p ∝ log^3(t/2 pi): Zador's point-density law lambda ∝ p^{1/3} is then the zero density, and
D_zeros / D_opt -> <delta^3>. General r: kappa_r = <delta^{r+1}>/<delta>^{r+1}; r -> infinity is the covering radius max delta / 2
(the large-gap constant Lambda of N8). Sine kernel: <s^2> = 2 int E(0;s) ds, <s^3> = 6 int s E(0;s) ds, <s^4> = 12 int s^2 E(0;s) ds,
E(0;s) = det(I - K_sine on [0,s]) by Nystrom (48 Gauss nodes): 1.179994, 1.577939 (1.981 dB), 2.333092.

## F  Fleischer bridge (the only RH-sensitive piece)
RH => Xi in the Laguerre–Polya class => Xi'^2 - Xi Xi'' >= 0 on R => |Xi| log-concave on every (gamma_n, gamma_{n+1}) => for every
alpha > 0, k >= 1 and every sub-interval J, the MSE k-point quantizer of |Xi|^alpha on J is unique (Fleischer 1964, Trushkin 1982).
Converse direction fails twice: (i) uniqueness holds for many non-log-concave densities (it is strictly weaker than log-concavity);
(ii) log-concavity between real zeros is the first Laguerre inequality, an RH consequence already verified wherever RH is
(|x| <= T - 1/2 for all orders, arXiv 2611.02944). Model Xi_delta = Xi ((t-g0)^2 + delta^2)/((t-ga)(t-gb)) (real zeros ga, gb replaced by t = g0 +- i delta):
(log|Xi_delta|)''(t) = 2(delta^2 - (t-g0)^2)/((t-g0)^2 + delta^2)^2 - sum_{other real zeros} (t-gamma)^-2 + (log F)'', so at t = g0 it is
2/delta^2 - O(gbar^-2): Laguerre fails for delta small against the mean gap gbar (measured: fails for every delta <= gbar tested);
for delta >> gbar the other zeros' -sum (t-gamma)^-2 wins and Laguerre is blind as well. Along (t-g0)^2 + delta^2 with delta^2 real,
the value at g0 is -2/|delta^2| for delta^2 < 0 (two real zeros) and +2/delta^2 for delta^2 > 0: the jump that quantizer functionals lack.

## RD  rate–distortion screen (N5, N1)
Gaussian source with spectral density S on a window W: reverse water-filling D(theta) = |W|^-1 int min(theta, S),
R(theta) = |W|^-1 int (1/2) log+(S/theta); both depend only on the distribution function of S (rearrangement invariant), and are defined
for every S >= 0. N5: the Nyman–Beurling space is the log-scale process with spectral density |zeta(1/2+it)|^2/|s|^2; RH <=> zeta(s)/s
outer in H^2(Re s > 1/2) (Beurling); the Kolmogorov–Szego/Jensen defect is the Balazard–Saias–Yor sum
(1/2 pi) int log|zeta(1/2+it)| dt/(1/4+t^2) = sum_{Re rho > 1/2} log|rho/(1-rho)| (BSY 1999): a phase statement, invisible to R(D).
N1: a Gaussian process with covariance = Weil distribution exists iff Weil positivity (RH); rate–distortion presupposes it.

# explore31 derivations (learning curves -> N5)

Notation: rho_theta(x) = {theta/x}, chi = chi_(0,1], d_Theta^2 = dist(chi, span{rho_theta : theta in Theta})^2 in L^2(0,inf);
Baez-Duarte d_N = d_Theta for Theta = {1/k : k <= N}; L = log N; C = 2 + gamma - log 4 pi; G = Vasyunin Gram of rho_{1/k};
t_k = <chi, rho_{1/k}> = (log k + 1 - gamma)/k; mu = dt/(2 pi |s|^2) on s = 1/2 + it.

## D1 The natural design has uniform leverage in u = log(1/theta)
Mellin: rho_theta -> -zeta(s) theta^s / s; on the line |theta^s|^2 = theta = e^{-u}. The points u_k = log k have density e^u.
Feature covariance of the natural design: Sum_k phi_k phi_k^* has kernel zeta(s) conj(zeta(s')) Sum_k e^{-u_k} e^{-i u_k (t-t')}
~ zeta(s) conj(zeta(s')) Int_0^L e^{-iu(t-t')} du (density x weight = e^u e^{-u} = 1). Its nonzero spectrum is that of the
truncated Wiener-Hopf operator W_L on L^2(0,L) with symbol w(t) = |zeta(1/2+it)|^2/(2 pi (1/4+t^2)), i.e. the Gram operator of
the continuous family {rho_theta : theta in [1/N, 1]} with log-uniform weight. Any formula that sees a design only through its
density and its feature covariance therefore predicts d_N^2 ~ delta(L)^2 := dist(chi, span{rho_theta : theta in [e^{-L}, 1]})^2.

## D2 The learning-curve (omniscient-risk) formula
Deterministic equivalent for n features sampled from a pool with Gram Gam and target vector beta = (<chi, phi_j>):
E[P_S] ~ T (T + kappa)^{-1}, T = pool covariance / M, tr T(T+kappa)^{-1} = n (ridgeless sketching / random-feature equivalent,
Bach 2023, Atanasov-Zavatone-Veth-Pehlevan 2024; same self-consistent equation as Bordelon-Canatar-Pehlevan and Simon et al.).
Push-through (Woodbury): E_LC(n) = |chi|^2 - beta^T (Gam + eps)^{-1} beta with tr Gam (Gam + eps)^{-1} = n, eps = M kappa.
Pinned version (rho_1 always in the design): Schur complement Gam' = Gam_11 - g10 g10^T/g00, beta' = beta_1 - g10 beta_0/g00,
|y'|^2 = 1 - beta_0^2/g00, then the formula with n = N - 1.

## D3 q-designs (exact Gram for a dense rational family)
D_{q,N} = {q/k : q <= k <= qN}. Scaling x -> x/q: <rho_{q/j}, rho_{q/k}> = q G_{jk}, <chi, rho_{q/k}> = q t_k with
t_k = (log(k/q) + 1 - gamma)/k (theta = q/k <= 1). So d_{q,N}^2 = 1 - q t^T G[q..qN]^{-1} t. D_{1,N} = {1/k}, D_{q,N} subset
D_{qq',N}, the union over q is dense in [1/N, 1] and theta -> rho_theta is L^2-continuous, hence d_{q,N} decreases to delta(log N)
along divisibility chains. Numerically d_{q,N}^2 - delta^2 ~ a/q; extrapolation by a + b/q (+ c/q^2) fits.
Mixed designs (natural plus head points) use the same index picture: dilation Q/k <-> index k, Gram Q G_{jk} with Vasyunin's
formula evaluated on the fly for arbitrary (j, k) (explore13's general form <rho_{a/b}, rho_{c/d}> = a c G_{bc, ad} agrees).
DH: A(x) = Sum_{n<=x} a_n vanishes on (0,1), so <chi, A(theta/t)> = theta beta, beta = Int_1^inf A(y) y^-2 dy; Gram q G^f_{jk}.

## D4 nb2nd's head functional for real dilations
l(x) = <zeta'/zeta, zeta x^{-s}>_mu (x >= 1 real) by the same residues at s = 0 and s = -2k as nb2nd eq. (2):
l(x) = -[log(2pi)(1 + log x - gamma) + a1]/x - Sum_k (zeta(2k+1) - 1) x^{-2k-1}/(2k(2k+1)) - [-(y/2) log(1-y^2) - atanh y + y],
y = 1/x, a1 = pi^2/12 - gamma^2 - 2 gamma_1 (the last bracket is Sum_k y^{2k+1}/(2k(2k+1)); its x = 1 value is 1 - log 2).
Check: l(1) = -1.7545539385, l(1)^2/G_11 = 2.441940 = nb2nd's Sigma_1. Fixed-head prediction (nb2nd step 3 with B ranging over the
span of a fixed head H): L^2 (d_nat^2 - d_{nat U H}^2) -> Sigma_{nat U H} - Sigma_nat, Sigma_H = l_H^T Gram_H^{-1} l_H.

## D5 Why density-based formulas fail: alignment
Sum_k c_k rho(1/(kx)) jumps at x = 1/m by an amount proportional to Sum_{k | m} c_k (each {1/(kx)} jumps by 1 where kx m' = 1);
Truncated Moebius c_k = mu(k) (k <= N) cancels the jumps at x = 1/m for 2 <= m <= N (Sum_{k|m} mu(k) = 0) and leaves the jump at x = 1
(chi's jump) plus jumps at 1/m, m > N, of size Sum_{k|m, k<=N} mu(k) (1 for primes m > N); the optimum is Moebius-like, not exactly
Moebius. Every natural feature has its jumps on the common lattice {1/m}, so they can cancel each other's jumps. A dilation theta not of the form 1/k jumps at {theta/n}, a set no other feature
shares; its jumps can only be smoothed, at L^2 cost of order the gap to the nearest aligned point. At the Baez-Duarte density
(theta prop. 1/theta^2) the interval [1/(k+1), 1/k] holds about one design point for every N, so a random design keeps an O(1)
head misalignment for every N: its error does not follow the density's continuum limit. (This is an explanation, consistent with
the head-localisation data of s3; the data establish the ensemble gap, not the mechanism.) For DH the jumps at 1/n have size
a_n in {1, 0.284, -0.284, -1, 0}, so the first misaligned jumps are smaller and the penalty milder (1.2-1.5x instead of 4-10x).

## D6 Side: liminf d_N^2 log N >= C exactly, unconditionally
d_N is non-increasing, so d_inf = lim d_N exists; RH <=> d_inf = 0 (Baez-Duarte 2003). If RH fails, d_N^2 log N >= d_inf^2 log N -> inf.
If RH holds, every zero is on the line and Burnol's Thm 1.3 gives liminf >= Sum_rho m_rho^2/|rho|^2 >= Sum_rho m_rho/|rho|^2
= Sum_rho (with multiplicity) 1/(rho(1-rho)) = C. Either way liminf >= C; equality forces RH and simple zeros. (No RH
verification height or tail bound needed; COROLLARIES/TREE N5 state C - 6e-12.)
Known: Burnol states the dichotomy himself (math/0103058 TeX l. 196-197: "If the Riemann Hypothesis fails this result is true but
trivial"). His theorem is for D(lambda) = dist(chi, span{rho(theta/t) : lambda <= theta <= 1}) in L^2(0,inf) (TeX l. 168-191), the
continuum delta(L) of D1-D3 with lambda = e^{-L}; Baez-Duarte's span is contained in B_{1/N}, so d_N >= D(1/N) and the bound transfers.

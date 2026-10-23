# Ramanujan's inequality for all x >= 38,358,837,683: a repaired, re-certified proof (RamRepair)

Ramanujan's inequality is pi(x)^2 < (e x/log x) pi(x/e). Padhi (Zenodo 10.5281/zenodo.23440770 v2, 2026-10-17;
code github.com/jaideepsaipadhi/ramanujans-inequality) claimed it for every x >= 38,358,837,683, unconditionally. Our audit
(Riemann/Deep/Ramanujan/) found a factor-2 undercount in his new zero-density bound; with it fixed his certificate fails near
log x = 875. This note repairs the argument, re-derives every remaining input and reruns the whole pipeline in ball arithmetic.
All computations: python-flint 0.9 (Arb balls), single-threaded; every number below is an enclosure or a certified bound.

## 0. Result
Theorem (unconditional). For every real x >= 38,358,837,683, pi(x)^2 < (e x/log x) pi(x/e); so N_R = 38,358,837,683.
The range e^58 <= x <= e^3159 (the gap between Johnston 2022 and Axler 2022) is proved here: all 317 intervals of the cover of
[57, 3159) have certified q < 1, max q = 0.92768 on [847, 857]; the transfer to pi gives max ratio 0.926851 at log y = 857 (< 1, 7.3% margin).
Inputs of the zeros-above-H0 part are all re-derived and certified here: zero-free region (own table, replaces Padhi Thm 4.6),
density cells (corrected Thm 6.1), Barban-Vehov mean square (own constants and remainders, certified k1, k0), (B1) box bounds;
Padhi Thm 5.2 (repulsion) and Bellotti's tables are not used. Buthe/smoothing constants and the Gamma integrals were rerun.
Published theorems used as black boxes are listed in 3.7. Section 5 states the theorem with credit.

## 1. What has to be shown, and how q encodes it
- Ramanujan at x follows from |pi(x) - li(x)| <= eps x/L and |pi(x/e) - li(x/e)| <= eps (x/e)/(L-1) with eps = eps_R(L),
  L = log x, the positive root of eps^2 + b(L) eps + g(L) = 0 (Padhi Lemma 2.1; algebra re-checked). eps_R(L) ~ 1/(3 L^4):
  3 L^4 eps_R = 1.2706 (L = 58), 1.0145 (907), 1.0041 (3159). Lemma 2.2 (evaluation without cancellation, A(L) = P(1/L) + R(L),
  -6e^{-L/2} <= R <= 16! L^-16/(1 - 34/L) + 6e^{-L/2}) re-derived: 16-fold integration by parts from sqrt x; the boundary term is
  at most 3 sqrt(x)/(L/2) since 0 < A, P <= 1.1 above 28.5; the integral is bounded by e^L L^-17/(1 - 34/L) because
  (log of e^s s^-17)' >= 1 - 34/L on [L/2, L].
- So with y = x/e one needs |pi(y) - li(y)| < eps_R(log y + 1) y/log y (57 <= log y <= 3158) and < eps_R(log y) y/log y
  (58 <= log y <= 3159). These come from psi via pi - li = (theta - y)/log y + int_2^y (theta - t)/(t log^2 t) dt + 2/log 2 - li(2),
  0 <= psi - theta <= 1.4262 sqrt t (Rosser-Schoenfeld), |psi - t| <= 0.94 sqrt t for t <= 1e19 (Buthe 2018 Thm 2).
- q: on each interval [L1, L2] of a cover of [43, 3159] the certified bound is
  |psi(t) - t| <= q * eps_R^-(L2 + 1)/(1 + 2/L1) * t + a(t), a(t) = 2 + e^{2e-4}(log t + 1e-4)(2.01e-4 sqrt t + (1/2) loglog(2t^2)),
  eps_R^- = certified lower end. q < 1 means the psi error uses the fraction q of the budget; the factor 1/(1 + 2/L1) is room for
  the integral term (it adds about 1/L + 6/L^2). q < 1 alone is not enough: the transfer computation (Padhi Prop. 2.4,
  transfer_arb.py, 0.05-steps in log t, no monotonicity of eps_R used) must give ratio < 1 for every log y in [57, 3159).
- q is the sum of four certified parts (Padhi Prop. 3.5 + Lemma 3.6, evaluated by his q_arb): smoothing (= f, the kernel
  fraction), truncation at T1 = c/eps (Buthe tail), zeros up to H0 = 3e12 (on the line), and the zeros above H0 (zero-free region +
  zero density). Only the last one involves the parts that were wrong or unchecked.

## 2. The false step, with the full original derivation
Bellotti, arXiv 2405.12545 (source final.tex), far off-diagonal term, lines 858-874 ((4.13)-(4.14) in the audit's numbering; TeX labels singdeltaT, boundc4), verbatim:
    858: \item  If $1<|\gamma_r-\gamma_s|\le T$, similarly to what we said for the previous case, by the construction \eqref{rectangles}, for every fixed $\gamma_r$ and for each $n=\delta^{-1},\dots,T\delta^{-1}-1$, there is a unique representative zero with imaginary part $\gamma_s$ such that $n\delta\le |\gamma_r-\gamma_s|\le (n+1)\delta$. Hence, by \eqref{asymboundgamma}, \eqref{sumrs} restricted to the case $1< |\gamma_r-\gamma_s|\le T$ is an application of Lemma \ref{stirling} gives
    859:  \begin{equation}\label{singdeltaT}
    860:         \begin{aligned}
    861:         &\frac{2d_5}{w\log T}\sum_{\substack{r\ne s\le J\\1<|\gamma_r-\gamma_s|\le T}} \left|\Gamma\left(2-\rho_r-\overline{\rho_k}\right)X^{2-\rho_r-\overline{\rho_k}}\right|\\
    862:             &\le \frac{2d_5}{w\log T}\cdot X^{2-2\alpha}\sum_{r=1}^J\sum_{1<|\gamma_r-\gamma_s|\le T}\left|\Gamma\left(2-\rho_r-\overline{\rho_k}\right)\right|\\&\le  \frac{2d_5}{w\log T}\cdot J\cdot X^{2-2\alpha}\cdot (2\pi)^{1/2}\cdot 2.4\cdot \int_{1/\delta}^{T/\delta}e^{1/6(n\delta)}e(n\delta)^{1.5-2\cdot\alpha}\cdot e^{-\pi n\delta/2}\text{d}n\\&\le  J\cdot \frac{X^{2-2\alpha}}{\delta \log T}\frac{2\cdot 2.4\cdot d_5\cdot (2\pi)^{1/2}}{w} \int_{1}^{T}e^{1/6y}y^{1.5-2\cdot\alpha_0}\cdot e^{-\pi y/2}\text{d}y.
    863:         \end{aligned}
    864:     \end{equation}
    865:     As before, multiplying \eqref{singdeltaT} by the quantity $d_4$ found in Lemma \ref{an2bn} we have
    866:     \begin{equation}\label{boundc4}
    867:          J\cdot c_4\cdot \frac{X^{2-2\alpha}}{\delta\log T},
    868:     \end{equation}
    869: where
    870:     \[
    871:     c_4=\frac{2\cdot (2\pi)^{1/2}\cdot 2.4\cdot d_5\cdot d_4}{w}\int_{1}^{T}e^{1/6y}y^{1.5-2\cdot\alpha_0}\cdot e^{-\pi y/2}\text{d}y.
    872:     \]
    873: \end{itemize}
    874: Summing \eqref{boundc2},\eqref{boundc3},\eqref{boundc4} together we can conclude that \eqref{sumrs} becomes
The sentence "for every fixed gamma_r and for each n ... there is a unique representative zero with imaginary part gamma_s such
that n delta <= |gamma_r - gamma_s| <= (n+1) delta" is false: the band n delta <= |gamma_r - gamma_s| <= (n+1) delta consists of two
intervals of heights (above and below gamma_r), each of which can hold a representative. The sum over n therefore counts one side
only; the correct bound has 2 J ... in place of J (her near term, line 843, "delta < |gamma_r - gamma_s| < 2", has the same wording; Padhi replaced that one by a two-sided LP).
Padhi keeps the one-sided count. His Lemma 6.4 (main.tex lines 736-753) bounds the number of representatives in a window of height
omega by n(omega) and replaces Bellotti's integral by
    744: \[
    745: I' = \delta\,n(\omega)\Big(\sum_{i\geq0}g_\alpha(1 + i\omega)\Big), \qquad g_\alpha(y) = e^{y/6}y^{3/2 - 2\alpha}e^{-\pi y/2}, \qquad \omega = \frac14,
    746: \]
and Lemma 6.2(v) (main.tex line 648) uses it once:
    648: \item[(v)] $\mathcal C = 1.0067\big(1 + 2\delta_k\mathcal L(\kappa_0\delta_k)\big) + 2\sqrt{2\pi}\cdot2.4\cdot I'$, where $\kappa_0 \geq 1$ is arbitrary, $\Lambda$ is defined in Lemma~\ref{lem:pair} and $I'$ in Lemma~\ref{lem:crowd}.
His notes/LEMMA_C.md says "tiling [1, oo) by windows [y_i, y_i + w] on each side" but sums one window per distance band, and
certify/arb/cells_arb.py line 262 adds 2 sqrt(2 pi) 2.4 Iprime_arb(k, lk) once. Correct: windows [1 + i omega, 1 + (i+1) omega]
above and below gamma_r, at most 2 n(omega) representatives per band. Effect (audit): far share of C 9-23%, log N+ up 0.09-0.20
per cell, q on [872, 877] 0.99137 -> 1.0726 (1.0327 at best kernel): Padhi's Theorem 1.1 is not proved as published.

## 3. The repaired pipeline: every input of the zeros-above-H0 part, re-derived
### 3.1 Zero-free region (replaces Padhi Thm 4.6; zfr.py, zfr_run.py, zfr_mp.py, zfr_lemma5.py, check_w*.py)
Method of Bellotti-Trudgian-Yang (BTY, arXiv 2603.21490): theta = 1.1338, w = sec^4(theta) (g*g) with g(v) = cos(kv - theta) - cos(theta)
on [0, D], k = tan(theta), D = 2 theta cot(theta) (identity proved symbolically, check_w3.py: w >= 0, even extension positive definite,
w decreasing on [0, D] since the autocorrelation of a symmetric log-concave function is log-concave); kappa_m = [1, -851, 780, -525,
171, 28, -29]/859 (kappa = 433/859); P(x) = |sum c_k e^{ikx}|^2/14912370 >= 0 (a0 = 1, a1 = 865534/497079, a = 2919857/828465, exact);
K = 16, T0 = 1e10. For a zero beta0 + it, t >= H0, eta = 1 - beta0, lambda0 = eta log t, l = log t, and the zero-free input
"no zero with beta >= 1 - A/log gamma, H0 <= gamma <= Kt + T0", sigma = 1 - A/log(Kt + T0), BTY's identity (their Lemma 6:
explicit formula for sum Lambda(n) n^-sigma f(log n) P(t log n), f(u) = eta w(eta u) sum kappa_m e^{-(2 sigma - 1) m u}) gives
  M1 := a1 (F(sigma-1+eta) + F(sigma-eta)) - [eta w0((a kappa/2) l - G - c10) + a0 F(sigma-1) + E2 (a lambda0 eta^2 + S2 eta^3) + E3 eta^3 + 1e-7] <= 0,
F(z) = sum kappa_m W((z + (2 sigma - 1) m)/eta), W the Laplace transform of w (closed form, checked against rigorous quadrature).
So M1 > 0 on a box of (l, lambda0) excludes zeros there. Our constants (all in Arb, per box at eta0 = eta rounded up, sigma0 rounded down):
- C_w(nu, r) = |w''(0)| + (|w'''(0)| + e^{-nu D}|w'''(D)| + int |w''''| e^{-nu u})/r, from four integrations by parts (w(D) = w'(D) =
  w''(D) = w'(0) = 0): |W(z) - w0/z| <= C_w |z|^-3 for Re z >= nu, |z| >= r (Ford's lemma, re-derived).
- Lemma 5 (pairs rho, 1 - conj(rho) contribute <= 0 when 0 <= Re z <= 2 sigma - 1): reduced to B(y) = sum_{m<=7} b_m (m - sgn(b_m) eps0/m)/
  (m^2 + y^2) > 0, eps0 = C_w eta0^2/(w0 (2 sigma0 - 1)^2); certified for eps0 = 1/500 by exact real-root isolation (sympy), and
  eps0 <= 1/500 asserted in every box.
- c10 = kappa sum_{15 log p <= u0} m_p log p. We certify f(u) >= kappa eta w0 = f(0) on [0, u0] (phi(u) = w(eta0 u) sum kappa_m e^{-d_m u} -
  kappa w0 >= 0 by Taylor bound near 0, phi(0) = 0, phi'(0) > 0, then ball bisection; valid for eta <= eta0, sigma >= sigma0 since w
  decreases and d_m is chosen per sign of kappa_m). Padhi verifies the same function but states only f >= kappa f(0) and uses kappa^2:
  a valid but weaker constant; this is why our table is about 0.5% above his. m_p = certified lower bounds of min_x G_p(1, x),
  G_p(1, x) = sum_{m<=15} p^-m P(mx) (bisection with derivative enclosure; m_5 = 0.0073564772, i.e. BTY's 0.00735709 is too large,
  as Padhi also noted).
- E2 = Sigma_kappa C1/(pi x0^2) multiplying a lambda0 eta^2 + S2 eta^3, S2 = sum_{k>=1} a_k log(k/2) (exact): the k >= 1 error terms
  (BTY Lemma 8) with |Re psi(1/4 + iY/2)| <= log(Y/2) for Y >= 10 (Stirling with remainder: Re psi <= log(Y/2) - 1/(24Y^2) + 1/(30|z|^4))
  and the symmetric-weight argument int log(kt + y)(x0^2 + y^2)^{-3/2} dy <= log(kt) 2/x0^2; eta log(kt/2) = lambda0 + eta log(k/2) exactly.
- E3 = a0 Sigma_kappa (C1 I0(x0)/(2 pi) + C2/sigma0^3): k = 0 terms; I0(x0) = int |Re psi(1/4 + iy/2)|(x0^2 + y^2)^{-3/2} dy by a
  certified upper sum (Re psi(1/4 + iy/2) is monotone in |y|); |W0(s/eta)| <= C2 eta^3/|s|^3 with |s| >= sigma0 (BTY and Padhi
  divide by x0^3: valid, weaker). For k >= 1 that term is C2 eta^3/(kt)^3 < 1e-30.
- G(sigma0): archimedean terms (BTY Lemma 12) with digamma values in Arb, + 1e-10 for the O(1/t) errors.
- 1e-7 covers: zeros above Kt + T0 (BTY bounds, re-evaluated for eta <= 0.0145, log t <= 46: <= 1.2e-9 per k, weight a0 + a = 4.52),
  the pole terms k >= 1 (BTY Lemma 13: <= 1e-10 eta), the digamma tails; total < 1e-8.
Iteration (zfr_run.py): grid L in {28.75, 29.00, ..., 40} step 1/4, then to 46 step 1/2; lambda(l) = value at the smallest grid point
>= l; start lambda = 1/5.558691 everywhere (Mossinghoff-Trudgian-Yang 2024, valid for t >= 2) and above 46. Step j assumes step j-1 at
heights <= Kt + T0 (A = lambda_{j-1}(log(K e^L + T0)) at the top of the cell; lambda/l non-increasing is enforced) and certifies
M1 > 0 on (l, lambda0) in (L_prev, L] x [lambda_{j-1}(L), lambda_j(L)] by adaptive ball bisection. 5 steps, last change 4.4e-6.
Result (out/zfr_table.json; 137,828 certified boxes; an independent serial rerun, zfr_run1.py, reproduces the table exactly):
| l | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 40 | 43 | 46 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ours lambda(l) | 0.21494 | 0.21447 | 0.21402 | 0.21360 | 0.21320 | 0.21283 | 0.21247 | 0.21213 | 0.21180 | 0.21146 | 0.21089 | 0.20979 | 0.20614 |
| Padhi Thm 4.6 | 0.21384 | 0.21345 | 0.21308 | 0.21272 | 0.21239 | 0.21206 | 0.21175 | 0.21145 | 0.21117 | 0.21084 | 0.21033 | 1/4.896 above | |
lambda(28.75) = 0.21506, i.e. sigma > 1 - 1/(4.650 log t) at log t = 28.75; 1/5.558691 above log t = 46. lambda(l)/l is non-increasing.
Constants from the self-test (out/zfr_selftest.log): w0 = 5.6727876, |w''(0)| = 50.06970, |w'''(0)| = |w'''(D)| = 117.26478,
C_w(69.4, 69.4) = 51.804, I0(0.493) <= 23.664, G(0.993) = 1.568399, G(0.9935164) = 1.568882, S2 = -0.882692, Sigma_kappa = 3.775320;
W closed form = rigorous quadrature at five test points; Lemma 5 numerator has no real root (eps0 = 1/500 and 1/400).

### 3.2 Repulsion (Padhi Thm 5.2): not used
All cells are evaluated with d* = 0 (the near-pair LP of Padhi Lemma 6.3 without excluded rectangles). Thm 5.2 is not needed.
### 3.3 Density cells (replaces Padhi Thm 6.1 / Lemma 6.2; cells_mine.py, cells_final.py)
N(1 - delta_k, e^l) <= nt * (d4 d5/w) * kap2 * C/(c1^2 - d4 c5') * e^{2 x delta_k l}/(delta_k l), for the same 3320 cells (k <= 300,
l in {29, 29.5, 30, 30.5, 31, 32, ..., 40}) and Padhi's parameters (u, tau, w, x) (any admissible choice is valid), nu = lambda(l)/l from 3.1.
- C = 1 + 2 delta Lp + 4 delta n(omega) sum_{i>=0} Gr(1 + i omega)/(1 + i omega): diagonal 1 (Gamma(1+s) <= 1, X^s/s increasing on
  [2 nu, 2 delta] as 2 nu x l >= 1); near pairs: LP over rectangles j >= 1 on both sides, b_j = Gr(m_j)/m_j, m_j = max(2 delta, (j-1)h),
  costs cost((j+1)h), budget M(eps, l) - cost(0), weak-duality upper bound; far pairs: both sides, at most n(omega) per window
  (Padhi Lemma 6.4's count), Gr(y) = |Gamma(1 + 2 delta + iy)|/Gamma(1 + 2 delta) (exact Gamma instead of Stirling; Gr(y)/y decreasing),
  terms y >= 60 by e^{-1.4 y}; best of omega in {0.25, 0.35, 0.5}, kappa_0 = h/delta in {1, 1.5, 2, 2.5, 3, 4}.
  M(eps, l) = min{(1/2) log(T/2 pi) + Z(eps), min_a [(1/2) log(T/2 pi) + (a^2 + 1/2)/(2a) Z(eps) + log(T/pi)/(8a)]}, Z = -zeta'/zeta(1+eps)
  (Hadamard + Stirling; the second from (a + cos phi)^2 >= 0 applied to -zeta'/zeta at 1 + eps + it and 1 + eps + 2it), + 1e-11 for O(1/t).
- kap2 = (1 + (Y/X)^{2 nu})/2, Y = U/l^2: in the residue Gamma(w)(X^w - Y^w), Re w = 2 - beta_r - beta_s >= 2 nu and Y < X give
  |X^w - Y^w| <= (1 + (Y/X)^{2 nu}) X^{Re w}; Bellotti and Padhi use 2 X^{Re w}. kap2 = 0.61-0.70 in the cells.
- c5' = c5/d4 with |G(1/2 + i tau)| <= sum |theta_j||theta_k| [j,k]^{-1/2} = sum_r h(r)(sum_{r|j} |theta_j| j^{-1/2})^2, h(p) = sqrt p - 1 >= 0,
  inner sum <= 4 sqrt(W)/(r log W) (sum_{m<=y} m^{-1/2} log(y/m) <= 4 sqrt y), so |G| <= 16 prod_p(1 + (sqrt p - 1)/p^2) W/log^2 W
  (product exact to 1e6, tail e^{0.002}); Bellotti's (4.15) W G(1) drops the signs of theta.
- c1 = e^{-l^2/U} - B1 - 2 e^{-l(l - 3u)}/l^2 (B1 as Padhi: 10.45 sqrt(VW), Hiary-Patel-Yang 0.618 |t|^{1/6} log|t|, I_Gamma <= 5.922;
  tail: Psi(n) = 0 for 1 < n <= U, |Psi(n) Theta(n)| <= d(n)^2 <= 4n and |n^-rho| <= n^-1/2, and n^{1/2} e^{-n l^2/U} decreases for n > U, so the
  tail is <= 4 int_U^oo n^{1/2} e^{-n l^2/U} dn <= 4.01 U^{3/2} e^{-l^2}/l^2 <= 2 U^3 e^{-l^2}/l^2 = 2 e^{-l(l-3u)}/l^2 (< e^-700); the dropped residue
  M(1) Gamma(1 - rho) X^{1-rho} (|Gamma| ~ e^{-pi gamma/2}, gamma > 3e12) is covered by an extra 1e-100).
- d5 = 1.0061(1 + (1.333 + 4 e^{-wl/2})/(wl)): Bellotti Lemma 4.3 = Graham's diagonalisation (j,k) = sum_{r|(j,k)} phi(r) with Ramare's
  |sum_{m<=x,(m,r)=1} mu(m) log(x/m)/m| <= 1.00303 r/phi(r) and sum_{r<=W} mu^2(r)/phi(r) <= log W + 1.333 + 3.95/sqrt W (checked).
- nt = max{1, floor(M/cost(h/2))} (zeros in one rectangle). Cells with c1^2 - d4 c5' < 0.2 or delta_k <= nu are dropped.
Example cells (out/cell_detail.log; u, tau, w, x = Padhi's parameters, nt, kp = kappa_0, then the d4 split, then the cell factors):
    80,31 u t w x [2.0744, 1.671, 0.6019, 5.0482] nt 1 kp 1.0
       d4 <= [1.7818294 +/- 9.44e-9]  D4 main [1.5783856 +/- 4.07e-8]  REM [0.1763 +/- 4.09e-5]  S# [6.674e-5 +/- 1.30e-10]
       R_i {'R1': '1.56e-03', 'R2': '1.75e-01', 'R3': '7.98e-14', 'R4': '8.51e-19', 'R5': '1.30e-05', 'R6': '3.61e-05', 'R7': '9.47e-05', 'R8': '2.11e-19'}
       C [2.76639 +/- 1.84e-6]  kap2 [0.627341 +/- 2.36e-7]  d5 [1.07800 +/- 8.50e-7]  c1 [0.97070109 +/- 1.42e-9]  c5/d4 [0.0183085 +/- 1.33e-8]  logN <= 5.704693658844008
    100,31 u t w x [2.0073, 1.6202, 0.5569, 4.7947] nt 1 kp 1.0
       d4 <= [1.8907263 +/- 4.74e-8]  D4 main [1.6712613 +/- 4.73e-9]  REM [0.1893 +/- 1.29e-5]  S# [7.917e-5 +/- 4.27e-9]
       R_i {'R1': '4.00e-03', 'R2': '1.85e-01', 'R3': '9.22e-14', 'R4': '9.48e-19', 'R5': '9.01e-05', 'R6': '2.27e-04', 'R7': '1.70e-04', 'R8': '1.80e-19'}
       C [3.05046 +/- 2.49e-6]  kap2 [0.637916 +/- 4.66e-7]  d5 [1.08383 +/- 3.76e-6]  c1 [0.96611700 +/- 2.53e-9]  c5/d4 [0.0181327 +/- 6.61e-9]  logN <= 6.218823613445571
    150,31 u t w x [1.9289, 1.5788, 0.4907, 4.5446] nt 1 kp 1.0
       d4 <= [2.0061108 +/- 3.71e-8]  D4 main [1.7637447 +/- 1.49e-8]  REM [0.2090 +/- 2.85e-5]  S# [9.294e-5 +/- 4.27e-9]
       R_i {'R1': '1.17e-02', 'R2': '1.95e-01', 'R3': '1.07e-13', 'R4': '1.06e-18', 'R5': '4.95e-04', 'R6': '1.15e-03', 'R7': '3.17e-04', 'R8': '1.54e-19'}
       C [3.73173 +/- 6.35e-7]  kap2 [0.648433 +/- 2.63e-7]  d5 [1.09440 +/- 4.59e-6]  c1 [0.95557618 +/- 2.44e-9]  c5/d4 [0.0163336 +/- 3.98e-8]  logN <= 7.501729501179532
    100,30 u t w x [2.0701, 1.621, 0.5808, 4.9413] nt 1 kp 1.0
       d4 <= [1.8838857 +/- 2.76e-8]  D4 main [1.6655863 +/- 3.67e-8]  REM [0.1882 +/- 2.47e-5]  S# [7.882e-5 +/- 1.15e-9]
       R_i {'R1': '4.18e-03', 'R2': '1.84e-01', 'R3': '9.18e-14', 'R4': '9.44e-19', 'R5': '9.05e-05', 'R6': '2.28e-04', 'R7': '1.84e-04', 'R8': '1.79e-19'}
       C [2.99163 +/- 3.18e-6]  kap2 [0.632395 +/- 1.20e-7]  d5 [1.08311 +/- 2.01e-6]  c1 [0.96605117 +/- 4.27e-9]  c5/d4 [0.0175696 +/- 2.05e-8]  logN <= 6.168115868810073
    200,32 u t w x [1.8286, 1.5615, 0.4236, 4.2953] nt 2 kp 1.0
       d4 <= [2.0785122 +/- 4.00e-8]  D4 main [1.8170946 +/- 4.42e-8]  REM [0.2263 +/- 2.42e-5]  S# [0.0001006 +/- 1.46e-8]
       R_i {'R1': '2.07e-02', 'R2': '2.02e-01', 'R3': '1.15e-13', 'R4': '1.13e-18', 'R5': '1.06e-03', 'R6': '2.39e-03', 'R7': '3.96e-04', 'R8': '1.44e-19'}
       C [4.49528 +/- 1.30e-6]  kap2 [0.658903 +/- 2.93e-7]  d5 [1.10537 +/- 3.40e-6]  c1 [0.94543431 +/- 3.27e-9]  c5/d4 [0.0150930 +/- 2.72e-8]  logN <= 9.561998942473506
    120,29 u t w x [2.0976, 1.5984, 0.5745, 4.9615] nt 1 kp 1.0
       d4 <= [1.9363205 +/- 3.43e-9]  D4 main [1.7082225 +/- 1.25e-8]  REM [0.1964 +/- 4.60e-5]  S# [8.576e-5 +/- 1.19e-9]
       R_i {'R1': '7.30e-03', 'R2': '1.88e-01', 'R3': '9.87e-14', 'R4': '9.97e-19', 'R5': '2.27e-04', 'R6': '5.49e-04', 'R7': '2.72e-04', 'R8': '1.64e-19'}
       C [3.19131 +/- 1.55e-6]  kap2 [0.632117 +/- 1.37e-8]  d5 [1.08666 +/- 1.13e-6]  c1 [0.96168758 +/- 2.41e-10]  c5/d4 [0.0162180 +/- 2.52e-8]  logN <= 6.621379731959285

### 3.4 Barban-Vehov mean square d4 (replaces Padhi Thm 6.5 / App. A; d4_consts.py, kconst_cert2.py, d4_mine.py, b1_cert.py)
d4 >= sum_n Psi(n)^2 n^-1 (e^{-n/X} - e^{-n l^2/U}). Following Padhi App. A: Gaussian smoothing chi = chi_flat + chi_sharp (width 0.02),
Minkowski, d4 <= S_flat + 2 sqrt(S_flat S_sharp) + S_sharp - e^{-l^2/U}, S_flat <= e^{-1/X} + D4 + REM.
- Main term (residues at w = 0 of Gamma(w) e^{Yw}[...], re-derived): D4 = Delta^-2 [(Y2^2 - Y1^2)/2 - gamma Delta + k1 (Y1 + Y2 - 2 gamma) + 2 k0],
  using Gamma(w) e^{Yw} = 1/w + (Y - gamma) + (Y^2/2 - gamma Y + gamma^2/2 + pi^2/12) w + ..., K_c = zeta(1+w) H_c = k1/w + k0 + ...,
  E(w)/w^2 = w^-2 + 0.0002 + ...; the pi^2/12 and varsigma^2/2 terms cancel between the U and V parts.
- k1 = H_c(0), k0 = H_c'(0) + gamma H_c(0), H_c(w) = (1/2 pi i) int_P A(w; -w-s, s) E(w+s) E(s)/(zeta(1-s) zeta(1+w+s) (w+s)^2 s^2) ds.
  Certified (kconst_cert2.py): on Re s = 0.38, |t| <= 80, rigorous acb.integral with the Euler product to p <= 400; truncation
  (e^tau - 1) int|g| with tau = sum_{p>400} sup|A_p - 1| = 0.0062642 (exact to 1e6 + Rosser-Schoenfeld tail) and int_0^80 |g| <= 3.6915,
  int |gD| <= 9.0078 (adaptive ball subdivision); d_w log A only through x2 = p^{-1-w-s}: |d_w log A_p| <= log p r2 (r1 + r3)/((1 - r2)(1 + r3 - r1 - r2));
  80 < |t| <= 400: |g| <= A* M_I iota(1.38) e^{-0.0004 t^2}/t^4 with A* = A*(0.62, 1.38, 1) = 10.9956; horizontal pieces e^{-64}.
  Result: k1 in -0.6055 +- 0.0075, k1 <= -0.59807; k0 <= 0.37066 (Padhi: k1 <= -0.5853, k0 <= 0.4214).
- A* bounds (d4_consts.py): sup over the 3-torus of phases of prod_{p<=100} |A_p| by branch and bound with derivative enclosures,
  100 < p <= 1e6 by log|A_p| <= |A_p - 1|, p > 1e6 by pi(t) < 1.25506 t/log t: A1 = 148.013, Ax = 7.7687, AV = 10.744, AU = 10.352,
  A_K = 17.664 (for e2 in [1.08, 1.40] it suffices to take e2 = 1.08: A_p is analytic in x2 on |x2| < 1, so its sup over |x2| <= r
  increases with r); boxes (crude: (1 + r1 + r2 + r3)(1 + r3)/((1 - r1)(1 - r2))) only multiply e^{-eta^2 T2^2/2}-small terms.
  Padhi: A1 <= 1202.27, Ax <= 14.56.
- (B1) certified (b1_cert.py, max modulus on the boundary, 225000 balls per edge): min |zeta| on the boundary of
  {0.6 <= sigma <= 1.3, |t| <= 450} >= 0.06998, so |1/zeta| <= 14.29 <= M_I = 14.57; |zeta| <= 10.593 <= M_Z = 11.89 off D(1, 0.1).
- Remainders R1-R8: Padhi's Lemma A-rem (label lem:A-rem, main.tex 1049-1110), integrand bounded factor by factor on each contour piece; all constants
  recomputed (Gamma integrals G_0 = 9.1053, G_1 = 11.5201, J-integrals, zeta and iota values, A* above). Where Padhi uses
  pi/(b^2(b +- c)) we use the exact int dt/((t^2 + b^2)((t + v)^2 + (b +- c)^2)) <= pi/(b (b +- c)(2b +- c)) (r_U, r_V).
- S_sharp: Padhi's Lemmas A-bsum and A-St (main.tex 1114-1160; Omega bounds for the Gaussian tail, Ramare-Akhilesh c0 < 1.3326 for sum mu^2/phi), eps = a/l_X optimised.
In the cells: d4 = D4 + REM + 2 sqrt(S_flat S_sharp) + S_sharp - e^{-l^2/U} with D4 = 1.58-1.82, REM = 0.18-0.23 (R2 dominates, see
section 4), S_sharp ~ 1e-4 (cell examples in 3.3). k1, k0 enclosures: [-0.6 +/- 0.0129], [0.3 +/- 0.0707] (upper ends used: -0.598066, 0.370655).

### 3.5 Explicit formula (Buthe) and smoothing constants
- Padhi Prop. 3.2 = Buthe 2016 (arXiv 1410.7015, zsum.tex lines 4-8) / Buthe 2018 (arXiv 1511.02032, method.tex lines 104-108):
  sum_{|gamma| > c/eps} |l(...) x^rho/(l(i/2) rho)| <= 0.16 (x+1)/sinh(c) e^{0.71 sqrt(c eps)} log(3c) log(c/eps) for x > 1, eps <= 1e-3,
  c >= 3. The scan uses 0.32/sinh(c) ... log(T1) (0.16 (x+1)/x <= 0.32) with c in [30, 110], eps <= 1e-4.
- |Theta| <= 2 (Buthe 2018 Prop. 2 proof: 0.025 + 1.87 + 0.006 + 0.008 = 1.909); |psi - t| <= 0.94 sqrt t, 11 < t <= 1e19 (Buthe 2018 Thm 2).
- Smoothing (Padhi Lemma 3.1): K(c) = 2 sqrt(c) e^{eps/2}(nu_c/2 + I0(c)/(200 sinh c)) <= 0.80319 < 0.7946*1.02 = 0.81049 for all c in
  {30, ..., 110} (rerun: max K = 0.8031848276113143  bound used 0.7946*1.02 = 0.810492  OK); proof re-read (Brun-Titchmarsh 2Y/log Y of Montgomery-Vaughan, prime powers, mass near 0).
- Logan kernel |l((rho - 1/2)/i)| <= e^{eps(c+1)/6} l(gamma), l(i/2) >= 1 (Padhi Lemma 3.3), Lemmas 3.4, 3.6, Prop. 3.5: re-derived in
  the audit; the evaluator q_arb is Padhi's code, unchanged apart from the three replacements of section 4.
### 3.6 RH height
H0 = 3e12: Platt and Trudgian, "The Riemann hypothesis is true up to 3*10^12", Bull. London Math. Soc. 53 (2021) 792-797. Used for:
zeros below H0 on the line (q's H-term e^{-L1/2} log^2 H0), the start of the zero-free table (log H0 = 28.7287 < 28.75), the cells
(zeros above H0 only), (B1) (|t| <= 450). KLN's own H0 = 3.0610046e10 is lower (their Table 1 needs only T >= that).
### 3.7 Inputs used as published theorems (not re-derived)
Mossinghoff-Trudgian-Yang 2024 (1/5.558691, start of 3.1); BTY Lemma 6 (explicit-formula identity) and Lemma 13 (pole terms);
Kadiri-Lumley-Ng 2018 Table 1 (N2; checked against the arXiv source table, all 20 rows); Hasanalizade-Shen-Wong 2022 (N1); Hiary-Patel-Yang
and Patel-Yang (zeta on the half line), Hiary (1.461); Ramare (two Moebius/phi sums), Ramare-Akhilesh (c0); Rosser-Schoenfeld; Dusart;
Montgomery-Vaughan; Buthe 2016/2018; Platt-Trudgian 2021; Johnston 2022 (x <= e^103); Axler 2022 (x >= e^3158.442).

## 4. Certified rerun (cell by cell, interval by interval)
- Cells (out/cells_final.json, 3320 cells, 487 s): 3316 admissible, 4 dropped because delta_k <= nu (those cells are empty by the
  zero-free region). Against Padhi's published cells, log N+ changes by -0.658 ... +0.041 (mean -0.416): the doubled far term is
  outweighed by kap2, the exact Gamma bound, our d4 and our larger nu; the few increases (<= 0.041) are all at k <= 65, next to the zero-free edge, where Padhi's repulsion
  bound d* > 0 (not used here) helps most.
  d4 in 1.653 ... 2.254, c1^2 - d4 c5' in 0.806 ... 0.936 (>= 0.2 required).
- Scan (out/scan_mine_0.txt; Padhi's q_arb with our cells, our zero-free table, N3 removed, N1 and N2 kept; Padhi's cover of
  [43, 3159] (332 intervals) and his kernel choices (c, f); 1701 s): every interval with L1 >= 57 has q < 1.
  max q on [57, 3159) = 0.927682 on [847, 857] (margin 7.2%); published max 0.99137 on [872, 877] (now 0.90817). Selected intervals:
  | [L1, L2] | q (ours) | parts: smooth / zeros > H0 | kernel | q published (flawed) |
  |---|---|---|---|---|
  | [57, 67] | 0.41880 | 0.4000 / 0.0000 | c=60 f=0.40 | 0.41880 |
  | [497, 507] | 0.40000 | 0.4000 / 0.0000 | c=60 f=0.40 | 0.40000 |
  | [797, 807] | 0.83137 | 0.6500 / 0.1814 | c=75 f=0.65 | 0.87070 |
  | [837, 847] | 0.91960 | 0.7000 / 0.2196 | c=110 f=0.70 | 0.96582 |
  | [847, 857] | 0.92768 | 0.6500 / 0.2777 | c=105 f=0.65 | 0.98559 |
  | [857, 862] | 0.91030 | 0.6000 / 0.3103 | c=105 f=0.60 | 0.97435 |
  | [872, 877] | 0.90817 | 0.5000 / 0.4082 | c=110 f=0.50 | 0.99137 |
  | [877, 882] | 0.89132 | 0.4000 / 0.4913 | c=95 f=0.40 | 0.99007 |
  | [1000, 1010] | 0.81149 | 0.4000 / 0.4115 | c=70 f=0.40 | 0.88862 |
  | [2000, 2010] | 0.40538 | 0.4000 / 0.0054 | c=55 f=0.40 | 0.40636 |
  | [3150, 3159] | 0.40001 | 0.4000 / 0.0000 | c=70 f=0.40 | 0.40001 |
  The truncation and H-terms are below 1e-30 in the binding range; q = f + (zeros above H0) there.

- Coverage (coverage_check.py): [57, 3159) covered by certified intervals with q < 1; [43, 57) bounded (needed only for the
  integral term). Transfer (transfer_arb.py, Padhi's code, ball arithmetic, 0.05-steps): max over log y in [57, 3159) of
  bound/(eps_R y/log y) = 0.926851 at log y = 857: TRANSFER OK. So |pi(y) - li(y)| < eps_R(log y + 1) y/log y for
  57 <= log y <= 3158 and < eps_R(log y) y/log y for 58 <= log y <= 3159, with 7.3% to spare, and Lemma 2.1 gives Ramanujan's
  inequality for every x in [e^58, e^3159].
- Not needed but available: the crude box bound A1box = 1.7e8 in R2 (REM ~ 0.18-0.23, mostly R2) could be replaced by A1 = 148
  (A_p is holomorphic on the closed polydisc, so its sup over a box of real parts is attained at the smallest real parts), which
  would lower d4 by about 10%; kernels were not re-optimised. Both would only enlarge the margin.

## 5. The theorem
Theorem (unconditional). For every real x >= 38,358,837,683, pi(x)^2 < (e x/log x) pi(x/e). Hence N_R = 38,358,837,683
(x = 38,358,837,682 violates the inequality: Dudek-Platt 2015).
Proof: 38,358,837,683 <= x <= e^103: Johnston 2022 (arXiv 2109.02249) Thm 5.1. e^58 <= x <= e^3159: sections 1-4 (Padhi's method,
repaired and re-certified). x >= e^3158.442: Axler 2022 (arXiv 2207.02486) Thm 1.1 (via Fiori-Kadiri-Swidinsky).
Credit: the framework (reduction, smoothed explicit formula, height-dependent zero-free region, log-free density cells with a
Barban-Vehov mean square, scan and transfer code) is J. S. Padhi's (Zenodo 10.5281/zenodo.23440770 v2). Our corrections: far term
counted on both sides; factor (1 + (Y/X)^{2 nu})/2 in the residue; rigorous bound for |G(1/2 + i tau)|; exact Gamma factors in the
pair sums; own zero-free region (c10 with kappa, sharper E2, E3, start from Mossinghoff-Trudgian-Yang); repulsion theorem not used;
Bellotti's tables (N3) not used; own Barban-Vehov constants and certified k1, k0; corrected c1 tail.
Novelty: before this, the best proved unconditional bound was N_R <= e^3158.442 (Axler 2022); Padhi's claim of N_R = 38,358,837,683
was not proved as published (audit). This is the first complete proof (subject to the cited theorems of 3.7).

## 6. Reproduce
code/ (python 3.14, python-flint 0.9, mpmath, sympy, numpy; run from code/ with the upstream repo unpacked as
code/repo/ramanujans-inequality-main from https://codeload.github.com/jaideepsaipadhi/ramanujans-inequality/zip/refs/heads/main):
zfr.py (constants, W, margins; self-test), zfr_mp.py (m_p), zfr_lemma5.py, check_w.py, check_w3.py, zfr_run.py / zfr_run1.py (table),
b1_cert.py, d4_consts.py (A*, Gamma/J integrals, zeta values), kconst_cert2.py (k1, k0), d4_mine.py (d4 bound), cells_mine.py +
cells_run.py (cell parts, uses Padhi's parameter tables), cells_final.py (cells), cells_compare.py, cell_detail.py, scan_mine.py
(scan), post.sh (coverage, transfer and the short certificates; set PYTHON to the interpreter), to.sh (hard time limit). Outputs in code/out/ (logs and JSON).
Total single-core time about 1.5 h.


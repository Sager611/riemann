# Explicit constant in Savva's two-thirds wedge — derivation

Source architecture: J. Savva, Zenodo 22293642 v1.1 (TeX l. = PAPER_SOURCE/JENSEN_TWO_THIRDS_MAIN.tex of GitHub
aimathpapers/Jensen-2-Thirds v1.1.0). Analytic input: our explicit complex saddle bound, Riemann/Deep/HollandK/DERIVATION.md
Sec. 2-3. Certificate: certify6.py (Arb ball arithmetic, python-flint; every inequality below is a chk() there).

## 0. Statement

J^{d,n}(X) = sum_{j<=d} binom(d,j) gamma(n+j) X^j, xi(1/2+w) = sum_n gamma(n) w^{2n}/n!.

**Theorem.** (i) Let D_all >= 9.0e24 be a degree up to which every J^{d,n} is hyperbolic (Chasse's sector lemma with RH to
height 3.0001753328e12; our TREE N6 range d0 = 1.614e25). For all d >= 1, n >= 0: **n^2 log(n+2) >= 24 d^3** implies that
J^{d,n} is hyperbolic; when d > D_all it has d distinct negative zeros. (The all-n range gives hyperbolicity only, not
simplicity, so for d <= D_all we claim hyperbolicity only.)
(ii) Without any all-n range (separate run, L* = 24, N* = Z(12) = 9.986e11): **n^2 log(n+2) >= 2.76e25 d^3** implies that J^{d,n}
has d distinct negative zeros, for every d >= 1.

Regime. Z(u) := u(pi e^{2u} - 9/4), u_0 = Z^{-1}, L := 2u_0(n), eps := 1/L, nu := 1/n, N* := Z(L*/2) exactly (an Arb ball in the
certificate; 6.9625e36 <= N* for L* = 80), so n >= N* <=> eps <= 1/L*. log(n+2) <= (1 + lamL) L, lamL = (log(pi L*/2) + 2/N*)/L*
(n < (pi L/2) e^L; the bound is decreasing in L): 1 + lamL <= 1.06043 at L* = 80 and <= 1.15124 at L* = 24 (outward-rounded).
(i): L* = 80, N* = 6.96e36. If d <= D_all, hyperbolic by hypothesis. If d > D_all: n^2 log(n+2) >= 24 d^3 > K_N D_all^3 with
K_N = N*^2 log(N*+2)/D_all^3 <= 5.65, so n > N*; and (d/n)^3 <= log(n+2)/(K n) <= log(N*+2)/N* (log(x+2)/x decreasing,
K >= 1): d/n <= omega^2 := (log(N*+2)/N*)^{1/3} = 2.3e-12.
(ii): L* = 24, N* = 9.99e11. K >= N*^2 log(N*+2) forces n >= N* for every d >= 1 and gives d/n <= 1/N* =: omega^2 (flag -self).

## 1. What is kept from Savva and what is replaced

Kept: six matches (Lemma 9.1, l. 615), finite-free factorization (l. 784), MMP Props 2.7(iii), 2.17 and the MSS max-root
product theorem in both orientations (l. 803-871), the terminating 3F2 recurrence in the decomposed form (13.6) (l. 905-926),
the global-maximum radius argument (l. 943-969), six-node Hermite-Genocchi (l. 1014), multiplier stability (l. 1064).
Replaced: the analytic input (Thm 7.1 + Cauchy transport + C_B6 + N_analytic) by HollandK Part A; the fixed-point box of
Sec. 11 by a per-cell Krawczyk certificate; the constants C_J = 8, C_loc = 12 + 8 sqrt 6, K_r = 4096, Omega radius 2 rho,
|z - j| <= 3 rho by certified sharp ones; stability with Cauchy radius lambda r, lambda = 1.33, using all six matches.

## 2. Analytic input (HollandK Part A, unchanged)

On D(n) = {|z - n| <= n/10}: h := log Mt = const + Psi(z) + G(z), |G - G(n)| <= etaG (0.05337 at L* = 80; 0.05523 at L* = 24),
Mt != 0; gamma(z) = C 4^{cz} Gamma(z+1) Mt(z)/Gamma(2z+1). Psi^{(k)}(z) = Q_k(t,s)/(z^{k-1} 2u_0(z)), t = 1/u_0(z),
s = 9/(4 pi e^{2u_0(z)}). New here: Q_6. With d/du t = -t^2, d/du s = -2s, Z'(u) = (9/4)(t + 2 - ts)/(ts),
Psi'' = 2t/Z', Psi^{(k+1)} = (1/Z')(-t^2 d_t - 2s d_s) Psi^{(k)}, Q_k = Psi^{(k)} Z^{k-1} 2/t (sympy, exact):
Q_2 = 2 - t + .., Q_3 = -2 + t^2 + .., Q_4 = 4 + t + .., Q_5 = -12 - 5t + .., **Q_6 = 48 + 26t - 28t^2 + ..** (t = 2 eps).
Cauchy: |G^{(k)}(w)| <= k! etaG/(n/10 - |w - n|)^k.

## 3. The six-match system in power-sum form

Model coefficient ratios m_j = (A)_j (C)_j D^j/((B)_j (D)_j (AC)^j); multiplier c(j) = gamma(n+j)/(gamma(n) S^j m_j),
S = B gamma(n+1)/gamma(n), so P(y) := J^{d,n}(-y/S)/gamma(n) = sum_j p_j c(j) y^j, p_j the coefficients of
p_F(y) = 3F2(-d, A, C; B, D; Dy/(AC)). c(0) = c(1) = 1, and c(2..5) = 1 iff (Lemma 9.1)
  Delta^2 log gamma(n+k) = Delta^2 log m_k, k = 0..3,
i.e. (duplication formula; linear terms drop; f(U) := log(1 + 1/U))
  Delta^2 h(n+k) = f(A+k) + f(C+k) - f(B+k) + f(n + 1/2 + k) - f(D+k).
Parametrize A = n/(a eps), C = t n, B = (t + w eps) n, D = n + 1/2 + delta eps n. Take Delta^j at k = 0 (j = 0..3),
use Delta^j f(U) = (-1)^j j! Phi_{j+1}(U), Phi_q(U) = int_{[0,1]^q} (U + sum u_i)^{-q} du, n^q Phi_q(U) = (n/U)^q (1 - theta q^2/U),
and divide by (-1)^j j! (j+1) eps / n^{j+1}. With k := j + 1:
  a^k eps^{k-1}/k + w R_k(t, w eps) + delta R_k(1 + nu/2, delta eps) = Y_{k-1}  (+- 6k L/n),
  R_k(T,U) := [T^{-k} - (T+U)^{-k}]/(kU) = int_0^1 (T + sU)^{-k-1} ds,
  Y_j := (-1)^j X_j/(j+1)!,  X_j := n^{j+1} L Delta^{j+2} h(n).
Since Delta^{j+2} h(n) is an average of h^{(j+2)} over [n, n+j+2],
  X_j in hull{ Q_{j+2}(t_w, s_w) (n/w)^{j+1} u_0(n)/u_0(w) : w in [n, n+j+2] } +- (j+2)! etaG (10/(1 - 50 nu))^{j+2} L/n,
with t_w in [2 eps (1 - 6 nu), 2 eps], s_w in [0, (9/4pi) e^{-1/eps}]; Q is evaluated in the mean-value form in t.
As eps -> 0: Y -> (2, 1, 2/3, 1/2) and the system becomes Savva's F = 0 (l. 689-702), root (a,t,w,delta) = (1/3, 2, 16/3, 1/3):
the Gauss-Radau rule for 2x dx on [0,1] with fixed node 1 (nodes 1, 1/2; weights 1/3, 2/3).
R_k: series in U/T with 24 terms and a geometric tail bound (ratio < 1/2 checked). Certificate per eps-cell (existence by Brouwer): let G_approx be the left side above without the O(L/n) term and
G_exact = G_approx + rho the exact matching map. For fixed n, G_exact(y) = sum_U s_U n^k Phi_k(U(y))/(k eps) - Y_{k-1} is continuous
(Phi_k is smooth in U > 0 and U(y) = n/(a eps), t n, (t + w eps) n, n + 1/2 + delta eps n are continuous and positive on Y).
On the box Y the certificate checks delta > 0, t > 1, a > 0 and a eps < 1 (for the cell [0, eh] this is 0 <= a eps < 1; every
actual n has eps > 0, hence 0 < a eps < 1 and A finite), so every n/U <= 1 and
|rho_k(y)| <= (k/eps) sum_U (n/U)^{k+1}/n <= 5k L/n <= 6k log(N*)/N* uniformly in y (L/n <= log n/n <= log(N*)/N*; no nu/eps
appears). Newton at the cell midpoint gives y_m and P ~ DG(y_m)^{-1}. For y in Y, the componentwise mean value theorem for
G_approx gives T(y) := y - P G_exact(y) in K(Y) := y_m - P [G_approx(y_m) +- 6k log(N*)/N*] + (I - P DG_approx(Y))(Y - y_m), where
G_approx(y_m) and DG_approx(Y) are enclosed with every datum (eps-cell, nu in [0, 1/N*], the Y-enclosures) as an interval. The
certificate checks K(Y) subset int(Y); then T(Y) subset Y, Brouwer gives a fixed point, and P is nonsingular (Rump), so
G_exact has a zero in K(Y) cap Y for every n in the cell. Uniqueness is not needed. Cells failing the test are bisected. Certified (L* = 80): a >= 0.32346, t in [1.95571, 2.00055],
w in [5.1802, 5.3370], delta in [0.33320, 0.34194], t - D/n >= 0.95144.

## 4. Positivity and simple roots of p_F

A - B - d > 0 (1 - (a eps)(B/n + d/n) > 0) and C - D - d > 0 (t - D/n - omega^2 > 0) are certified. Then both factors of
p_F = Q^(1)_{A,B} [x]_d Q^(D)_{C,D} (Savva l. 784; checked coefficientwise) are rescaled Jacobi polynomials
P_d^{(alpha,beta)}, alpha = V - 1, beta = U - V - d >= 0 ((U,V) = (A,B), (C,D)): d simple roots in (0, U) (resp. (0, C/D)).
MMP 2.7(iii) and 2.17 give d simple positive roots of p_F.

## 5. Localization (sharp Gershgorin) and MSS

q_{U,V}(Y) = 2F1(-d, U; V; Y/U) = const P_d^{(alpha,beta)}(1 - 2Y/U). Monic Jacobi recurrence: b_k = (beta^2 - alpha^2)/
((2k+S)(2k+S+2)), a_k^2 = 4k(k+alpha)(k+beta)(k+S)/((2k+S)^2((2k+S)^2 - 1)), S = alpha + beta = U - d - 1. In Y = U(1-s)/2:
  diag(k) - V = [V S (d-1) + 2k(S+k+1)(U - 2V)]/((2k+S)(2k+S+2)),
  off(k)^2 = U^2 k(k+V-1)(k+U-V-d)(k+S)/((2k+S)^2((2k+S)^2 - 1)).
For 0 <= k <= d: |diag - V|/sqrt(Vd) <= sqrt(d/n) [V' s1 u + 2(s1 + (d/n + nu)u)|1 - 2V'u|]/(s1^2 sqrt V'),
off^2/(Vd) <= (1 + d/(nV'))(1 - V'u)(s1 + u d/n)/(s1^2 (s1^2 - nu^2 u^2)), u = n/U, V' = V/n, s1 = S/U. Gershgorin:
roots of q_{U,V} in V +- C_J sqrt(Vd), C_J = sup|diag - V|/sqrt(Vd) + 2 sup off/sqrt(Vd).
Certified: C_J1 <= 2.00001 ((A,B)), C_J2 <= 1.41441 ((C,D)); eta := C_J2 sqrt(d/D) <= 2.2e-6.
MSS max-root theorem for p_F and for its reversal (positive lower endpoints certified): every root, hence every critical point,
y of p_F has |y - B| <= C_loc sqrt(Bd), C_loc = C_J1 (1 + eta) + C_J2 sqrt(B/D) <= 4.0006.

## 6. Critical radius (Savva Sec. 14 with sharp constants)

T_k := y^k p_F^{(k)}(y)/p_F(y). Recurrence P3 T_{m+3} + P2 T_{m+2} + P1 T_{m+1} + P0 T_m = 0 (0 <= m <= d-2), kappa = Dy/(AC):
P3 = 1 - kappa, P2 = B + D + 2m + 1 - kappa(A + C + 3m - d + 3),
P1 = (B+m)(D+m) - kappa[(A+m)(C+m) + (m - d + 1)(A + C + 2m + 1)]  (Savva's App l. 688 prints (m - d): typo; (13.6) is right),
P0 = kappa (d - m)(A + m)(C + m).
Write P1 = D(B - y) + m(B + D) + m^2 - Dy[(1 + m/A)(1 + m/C) - 1] - kappa(m - d + 1)(A + C + 2m + 1); with |y - B| <= C_loc sqrt(Bd),
m <= d <= omega^2 n: P2/n >= 1.9781, |P1|/(n sqrt(Bd)) <= 4.0007, P0/(n^2 d) <= 2.0306, |P3| <= 1 (hulls over all cells).
R = K_r sqrt(Bd), K_r = 2.3, M = max_j |T_j|/R^j; T_0 = 1, T_1 = 0, T_{d+1} = 0. If M > 1 the maximizer k >= 2 and the recurrence
at m = k-2 gives 1 <= |P3| R/P2 + |P1|/(R P2) + |P0|/(R^2 P2) =: budget, a contradiction once budget < 1. The certificate evaluates
the budget per cell with the cell's correlated enclosures: max over cells <= 0.97277 (substituting the hulls above instead gives
<= 0.979, still < 1). So |T_k| <= R^k for all k (including k = d-1, d: T_{d+1} = 0).

## 7. Residual

E_F(z) = log gamma(n+z) - log gamma(n) - z log S'' - log m(z) (principal log Gamma on Re > 0, h on D(n)); holomorphic on
|z| < n/20, real on the reals, E_F(0..5) = 0, exp E_F(j) = c(j).
E_F^{(6)}(z) = h^{(6)}(n+z) - psi^{(5)}(n + 1/2 + z) - psi^{(5)}(A+z) - psi^{(5)}(C+z) + psi^{(5)}(B+z) + psi^{(5)}(D+z).
|psi^{(5)}(w) - 24 w^{-5}| <= 240 (Re w)^{-6}; exact pairings x^{-5} - y^{-5} = (y - x) sum_{i<=4} x^{-(i+1)} y^{-(5-i)}:
n^5 L E_F^{(6)}(z) = Q_6 (n/(n+z))^5 u_0(n)/u_0(n+z) - 24 w sum_i (b + z')^{-(i+1)} (t + z')^{-(5-i)}
  - 24 delta sum_i (D/n + z')^{-(i+1)} (1 + nu/2 + z')^{-(5-i)} - 24 a^5 eps^4 (1 + a eps z')^{-5} + err,
z' = z/n, b = B/n, err <= 720 etaG 20^6 L/n + 240 (L/n) sum (Re w/n)^{-6}. On |z| <= R0 := d + lambda R (R0/n <= 6.6e-6, so
n + z and the Cauchy circles of radius n/20 lie in D(n)), complex-ball evaluation per cell: **C6 := sup n^5 L |E_F^{(6)}| <= 2.0504**
(leading value -2 = 48 - 10 - 40).

## 8. Hermite-Genocchi, stability, assembly

Omega := {dist(z, [0,d]) <= lambda R}, convex, contains 0..5 (d >= 5), inside |z| <= R0. Six-node HG:
sup_Omega |E_F| <= (C6/(720 n^5 L)) (d + lambda R)^6 = Cfin d^3/(n^2 L), Cfin = C6 (lambda K_r sqrt(B/n) + sqrt(d/n))^6/720,
computed per cell (C6 and B/n from the same cell; C6 is largest near eps = 0 where B/n ~ 2.0006, B/n is largest near eps = 1/80):
max over cells Cfin <= 18.752. (Combining the hulls C6 <= 2.0504 and B/n <= 2.0673 would give 20.6, i.e. K = 25.3.)
Stability (Savva Lemma 16.1, radius lambda R): P/p = sum_k (Delta^k c(0)/k!) T_k; for k >= 1, |Delta^k c(0)/k!| = |c[0..k]|
<= sup_[0,k] |c^{(k)}|/k! <= e/(lambda R)^k (real Hermite-Genocchi; Cauchy on circles of radius lambda R centred on [0,d]),
e := sup_Omega |c - 1|; k = 1..5 vanish; so at critical points |P/p - 1| <= e/(lambda^5 (lambda - 1)) < 1 if e < lambda^5(lambda-1)
= 1.37332 (lambda = 1.33), which holds if sup_Omega |E_F| <= thetaE := log(1 + lambda^5(lambda - 1)) - 0.001 = 0.86329. Then P
keeps the sign of p at the d-1 critical points, at 0 and at +infinity (c(d) > 0): d distinct positive zeros; y = -SX gives d
distinct negative zeros of J^{d,n}. For 1 <= d <= 5, P = p_F exactly.
Sufficient: n^2 L >= (Cfin/thetaE) d^3, implied by n^2 log(n+2) >= K_res d^3, K_res = (1 + lamL) Cfin/thetaE.
(i) L* = 80, K_r = 2.3, lambda = 1.33, 1000 cells: K_res <= 23.0332 (certificate, unrounded balls; from the rounded constants
1.06043 * 18.752/0.86329 <= 23.035), K_N(9e24) <= 5.6413, K_N(d0) <= 0.97797: **K = 24**. Per-cell witnesses (outward-rounded
C6, B/n, Cfin, radius budget, C_loc for each of the 1000 cells): cells_main.csv; the maximum Cfin = 18.7514 is attained in the last
cell eps in [0.0124875, 0.0125] (C6 <= 1.99342, B/n <= 2.02192).
(ii) L* = 24, K_r = 4, lambda = 1.33, -self (d/n <= 1/N*): C6 <= 2.4512, Cloc <= 4.0050, budget <= 0.5459, K_res = 830.0,
**K = N*^2 log(N*+2) <= 2.756e25** (400 cells, bisected to 2210; witnesses cells_self.csv).

## 9. Cited, not re-proved

MSS max-root theorem for multiplicative finite free convolution; MMP Props 2.7(iii), 2.17; zeros of Jacobi polynomials =
eigenvalues of the truncated Jacobi matrix; Gershgorin; Polya's integral and the HollandK saddle derivation; Chasse's sector lemma
+ the verified RH height (for D_all).

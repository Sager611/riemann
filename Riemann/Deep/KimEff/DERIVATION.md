# Explicit Kim–Lee-shape hyperbolicity region for Riemann's Jensen polynomials — derivation

Certificate: certify.py (python-flint Arb balls, 256 bits). Every inequality marked [cert: NAME] is evaluated there
(NAME is the label it prints); the script exits 1 on any failure. Logs: logs/certify.log (a = 0.37338196661320752),
logs/certify_half.log (a = 1/2).

## 0. Statement

Notation. xi(s) = s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)/2; 8 xi(1/2+z) = sum_m gamma(m) z^{2m}/m!;
Theta(w) := sum_m gamma(m) w^m/m! = 8 xi(1/2+sqrt w) (Kim–Lee's Xi_0(-w) up to the factor 8);
J^{d,n}(X) = sum_{j<=d} binom(d,j) gamma(n+j) X^j = J(Theta^{(n)}; d)(X), the family of our N6 records.
Z(u) := u(pi e^{2u} - 9/4) (increasing on u > 0); u(n) := Z^{-1}(n).
a := 0.37338196661320752 (every nontrivial zero rho of zeta has |Re rho - 1/2| <= a: the Lean-checked half-plane,
repo input); everything below also holds with a = 1/2 (the classical critical strip; certify_half.log).

**Theorem.** Let n_* := ceil(Z(12)) = 998,616,377,771 [cert: n_*]. For every integer n >= n_* put
R_n := (9/10)(n - 1/2)^2/(u(n) + 1)^2. Then

(i) Theta^{(n)}(w) != 0 for every w with (Im w)^2 <= 4a^2(a^2 - Re w) and Re w >= -R_n;

(ii) J^{d,n} is hyperbolic for every d <= D_n := R_n^2/(4a^2(R_n + a^2));

(iii) in particular J^{d,n} is hyperbolic whenever n >= n_* and d <= 6.45 n^2/log^2 n (a = 0.37338...), resp.
d <= 3.59 n^2/log^2 n (a = 1/2). Consequently (Kim–Lee shape; sufficient subregions, c > 1/2) J^{d,n} is hyperbolic for
every d >= 1 and every n >= max(n_*, C_c d^c), with (c, C_c) = (0.51, 7.6890), (0.55, 0.87160), (0.6, 0.069829),
(1, 1.1853e-10) [cert: power forms, constants rounded up; a = 1/2: 10.367, 1.2031, 0.099247, 2.1294e-10].
All printed constants are rounded outward (lower bounds down, upper bounds and sufficient thresholds up).

No RH verification is used. Combined with the all-n range D_all = 16,140,857,488,599,731,839,205,676 (TREE N6):
D_n > D_all for every n >= 4.6974173e13 [cert: crossover], and for every d > D_all the region of (iii) contains
the two-thirds wedge n^2 log(n+2) >= 24 d^3 (Sec. 10). Without the all-n range the containment holds for every
d > 1.04710e8 [cert: self-contained].

## 1. Cited inputs

(I1) Riemann's integral (Titchmarsh, Theory of the Riemann zeta-function, (2.16.1)-(2.16.2), variable u -> t):
xi(1/2 + z) = 4 int_0^oo phi(t) cosh(z t) dt, phi(t) = sum_{k>=1} phi_k(t),
phi_k(t) = pi k^2 e^{5t/2}(2 pi k^2 e^{2t} - 3) exp(-pi k^2 e^{2t}).
Hence, with Mt(m) := int_0^oo t^{2m} phi(t) dt: gamma(m) = 32 m! Mt(m)/(2m)!. Checked: gamma(m)/(m! Mt(m)/(2m)!) = 32.0
to 15 digits for m = 0..6 (norm.py, logs/norm.log).
(I2) Zeros of zeta: |Re rho - 1/2| <= a (a = 0.37338196661320752 repo input; a = 1/2 classical).
(I3) Obreschkoff (Acta Math. Acad. Sci. Hungar. 12 (1961) 175-184), as Kim–Lee arXiv 2105.05386 Thm 3: if P, Q are
real polynomials, Z(P) in S(delta) := {|Im z| <= delta |z|}, Q hyperbolic and deg Q <= delta^{-2}, then P(D)Q is
hyperbolic. Corollary (Kim–Lee): Z(P) in S(delta) => J(P; d) hyperbolic for d <= delta^{-2}, since
J(P;d)(z) = z^d (P(D) z^d)(1/z). (Sharp: P = z^2 - 2z cos t + 1, Q = z^d gives discriminant 4d(1 - d sin^2 t).)
(I4) Bessel functions (DLMF): J_nu(z) = (z/2)^nu/Gamma(nu+1) prod_k (1 - z^2/j_{nu,k}^2) (10.21.15);
0F1(; nu+1; -z^2/4) = Gamma(nu+1)(z/2)^{-nu} J_nu(z) (10.16.9); 0 < j_{nu,1} < j_{nu,2} < ... real, and
j_{nu,1} > j'_{nu,1} >= nu for nu >= 0 (10.21.3). Checked numerically for nu = 4.5, 19.5, 59.5 (sanity3.py).
(I5) Hadamard factorization; Hurwitz's theorem (limits of real-rooted real polynomials of bounded degree, with
nonvanishing leading coefficient in the limit, are real-rooted).

## 2. The zeros of Theta^{(n)} lie in a parabola

(2.1) A zero of Theta is w = z^2 with xi(1/2+z) = 0, z = (Re rho - 1/2) + i Im rho =: b0 + i t, |b0| <= a.
Then Re w = b0^2 - t^2, Im w = 2 b0 t, (Im w)^2 = 4 b0^2 (b0^2 - Re w) <= 4a^2(a^2 - Re w)
(b0^2 - Re w = t^2 >= 0, b0^2 <= a^2). So Z(Theta) is contained in P := {w : (Im w)^2 <= 4a^2(a^2 - Re w)},
which is closed and convex (Re w <= a^2 - (Im w)^2/(4a^2), a concave bound).
(2.2) Theta has order 1/2, so every Theta^{(n)} has order 1/2, genus 0: Theta^{(n)}(w) = gamma(n) prod_j (1 - w/w_j),
sum 1/|w_j| < oo; gamma(n) > 0 by (I1) (phi > 0, Sec. 5.1).
(2.3) Gauss–Lucas for genus 0: if f = c prod(1 - w/w_j), sum 1/|w_j| < oo, and w0 lies outside the closed convex hull H
of the w_j, a separating unit vector e gives Re(conj(e)(w0 - w_j)) > 0 for all j, hence Re(e/(w0 - w_j)) > 0 and
Re(e f'(w0)/f(w0)) = sum_j Re(e/(w0 - w_j)) > 0 (absolutely convergent), so f'(w0) != 0. By induction
Z(Theta^{(n)}) is contained in the closed convex hull of Z(Theta), hence in P, for every n.

## 3. From a zero-free region to hyperbolicity (Obreschkoff step)

Suppose Theta^{(n)} has no zero in W_R := {w in P : Re w >= -R}, R > 2a^2. By (2.3) every zero w lies in P with
s := -Re w > R, so |Im w|/|w| <= |Im w|/s <= 2a sqrt(a^2 + s)/s <= delta := 2a sqrt(a^2 + R)/R (decreasing in s).
Let P_k be the partial products gamma(n) prod_{j in F_k}(1 - w/w_j) over finite conjugation-closed sets F_k exhausting
the zeros (with multiplicity): real polynomials with Z(P_k) in S(delta), P_k -> Theta^{(n)} locally uniformly
(absolute convergence). By (I3) J(P_k; d) is hyperbolic for d <= delta^{-2}; J(P_k;d) -> J(Theta^{(n)};d)
coefficientwise (Weierstrass), and the limit has degree d (leading coefficient gamma(n+d) > 0), so by (I5)
J^{d,n} is hyperbolic for every d <= delta^{-2} = R^2/(4a^2(R + a^2)) = D(R). (This is the proof of Kim–Lee Thm 1
with their ineffective input, Kim's 1996 zero-location theorems, replaced by the explicit Sec. 8.)

## 4. Mixture identity

Put b := n + 1/2, nu := n - 1/2, G(x) := 0F1(; b; x) = sum_k x^k/((b)_k k!), p_n(t) := t^{2n} phi(t)/Mt(n) on (0, oo).
**Theta^{(n)}(w) = gamma(n) int_0^oo p_n(t) G(t^2 w/4) dt.**
Proof: Theta(w) = 32 sum_m Mt(m) w^m/(2m)! = 32 int phi(t) cosh(t sqrt w) dt; differentiate n times termwise:
Theta^{(n)}(w) = 32 int phi(t) t^{2n} F_n(t^2 w) dt, F_n(v) = sum_k (n+k)! v^k/(k!(2n+2k)!), and
(n+k)!/(2n+2k)! = (n!/(2n)!)/(4^k (n+1/2)_k), so F_n(v) = (n!/(2n)!) G(v/4); finally gamma(n) = 32 n! Mt(n)/(2n)!.
Fubini: sum_k |coefficient| |t^2 w|^k <= F_n(t^2 |w|) <= F_n(0) e^{t^2|w|/(4n+2)}, integrable against phi(t) t^{2n}.
Checked: n = 3, w = -50 + 10i, direct Taylor sum = mixture = 1.46779162173e-5 + 2.99542734655e-6 i (sanity1.py).

## 5. Facts used

Throughout n >= 1, so b = n + 1/2 >= 3/2 and nu = n - 1/2 >= 1/2 > 0 (needed for j_{nu,1} > nu).

(5.1) phi_k(t) > 0 for t >= 0 (2 pi k^2 e^{2t} > 3). phi_1(t) = 2 pi^2 e^{9t/2 - pi e^{2t}}(1 - beta e^{-2t}),
beta := 3/(2 pi). For k >= 2, phi_k/phi_1 = k^2 (2 pi k^2 q - 3)/(2 pi q - 3) e^{-pi(k^2-1)q}, q = e^{2t} >= 1, is
decreasing in q (d/dq log = -6 pi(k^2-1)/((2pi k^2 q-3)(2pi q-3)) - pi(k^2-1) < 0). Hence for t >= 0:
phi_1 <= phi <= (1 + eps0) phi_1, eps0 := sum_{k>=2} k^2(2 pi k^2-3)/(2 pi-3) e^{-pi(k^2-1)} <= 0.0021760612
[cert: eps0; k = 2..11 summed in Arb; tail: for k >= 12, k^2(2 pi k^2 - 3)/(2 pi - 3) <= 2 k^4 and
pi(k^2 - 1) - 4 log k >= 3 k^2, so the tail is <= sum_{k>=12} 2 e^{-3k^2} <= 4 e^{-432} < e^{-300}].
(5.2) log phi_1 is concave on [0, oo): (log(2 pi e^{2t} - 3))'' < 0 and (-pi e^{2t})'' < 0.
(5.3) G(x) = prod_k (1 + x/c_k), c_k := j_{nu,k}^2/4 (I4); S1 := sum 1/c_k = 1/b, S2 := sum 1/c_k^2 = 1/(b^2(b+1))
(compare x and x^2 coefficients: G = 1 + x/b + x^2/(2b(b+1)) + .., e2 = (S1^2 - S2)/2); c_1 >= nu^2/4 (I4).
Checked numerically (sanity3.py: S2 from 400 zeros agrees to 1e-9; S1 up to the slowly converging tail).
For X < c_1 put tau := X/c_1 (< 1), tau+ := max(tau, 0), x = -X + iY. Then, with all c_k - X > 0:
(5.4) |G(x)| >= G(-X) > 0 (|1 + x/c_k| >= Re(1 + x/c_k) = 1 - X/c_k > 0).
(5.5) log|G(x)| <= log G(-X) + Y^2 S2/(2(1 - tau+)^2) (|1+x/c|^2 = (1 - X/c)^2 (1 + Y^2/(c-X)^2), log(1+y) <= y,
sum 1/(c_k - X)^2 <= S2/(1 - tau+)^2).
(5.6) X >= 0: -X S1 - X^2 S2/(2(1 - tau)) <= log G(-X) <= -X S1 - X^2 S2/2 (log(1-y) in [-y - y^2/(2(1-y)), -y - y^2/2]);
X < 0: log G(-X) >= 0.
(5.7) g(X) := sum_k 1/(c_k - X): 0 <= X < c_1: S1 <= g <= S1 + X S2/(1 - tau)
(1/(c - X) = sum_m X^m/c^{m+1}, convergent as X < c_1 <= c, and sum_k c_k^{-m-1} <= S2 c_1^{-(m-1)} for m >= 1);
X < 0: 0 < 1/(c_k - X) <= 1/c_k, so 0 < g <= S1.
(5.8) sum_k c_k/(c_k - X)^2 <= S1 + X+ S2 (2 - tau+)/(1 - tau+)^2: for 0 <= X < c_1 from c/(c-X)^2 = sum_m (m+1) X^m/c^{m+1};
for X < 0 directly c/(c - X)^2 <= 1/c.
(5.9) |G(x)| <= G(|x|) <= e^{|x|/b} ((b)_k >= b^k).
(5.4)-(5.8) were checked on 120 random points, b = 5.5, 20.5, 60.5 (sanity3.py: 0 failures).

## 6. The tilted log-density

For w = -s + iv in W_R (so -a^2 <= s <= R, v^2 <= 4a^2(a^2 + s)) and t > 0 put x(t) := t^2 w/4 = -X(t) + iY(t),
X = t^2 s/4, Y = t^2 v/4. h(t) := 2n log t + 9t/2 - pi e^{2t}; u := u(n) is the unique zero of h' (h'(u) = 0 iff
2n/u + 9/2 = 2 pi e^{2u} iff Z(u) = n). Useful identities (from 2 pi e^{2u} = 2n/u + 9/2):
  h'(t) = 2n(u - t)/(tu) + 2 pi e^{2u}(1 - e^{-2(u - t)})   (t < u),
  h'(t) = -2n(t - u)/(tu) - 2 pi e^{2u}(e^{2(t - u)} - 1)   (t > u).
On 0 < t < T_s (X(t) < c_1) define psi(t) := 2n log t + log phi_1(t) + log G(-X(t)). Then
  psi'(t) = h'(t) + 2 beta e^{-2t}/(1 - beta e^{-2t}) - (t s/2) g(X(t)),
  psi''(t) = h''(t) - 4 beta e^{-2t}/(1-beta e^{-2t})^2 - (s/2) g(X) - (t^2 s^2/4) sum_k 1/(c_k - X)^2,
so -psi'' >= 2n/t^2 + 4 pi e^{2t} - a^2/(2b) > 0 (for s < 0 the term -(s/2) g is <= a^2 S1/2): psi is strictly
concave. For 0 < t < T_s, t^{2n} phi_1(t) |G(x(t))| >= e^{psi(t)} (5.4) and t^{2n} phi(t) |G(x(t))| <=
(1 + eps0) e^{E+} e^{psi(t)} for t <= U (5.1, 5.5), where U and E+ are fixed in Sec. 7: U^2 R = thetaU nu^2 <= thetaU 4 c_1
gives tau+ <= thetaU < 1 and |Y| <= U^2 vmax/4 for t <= U, so the (5.5) correction is <= E+.

## 7. Parameters

theta = 9/10, thetaU = 49/50, delta_R = 1, Delta = 17/25, K = 4, eps_eta = 1/1000, u_* = 12; rho := sqrt(thetaU/theta)
>= 1.043498389. For n >= n_* (u = u(n) >= 12): R := theta nu^2/(u + 1)^2; u_L := u - Delta, u_R := u + 1/n;
kappa0 := 4 pi e^{2(u_L - 2 eps_eta)} - a^2/(2b); eta := sqrt(K/kappa0); U := (u + 1) rho; vmax := 2a sqrt(R + a^2)
(>= |v| on W_R); Wmax := sqrt(R^2 + vmax^2) (>= |w|); Sigma := S1 + (theta nu^2/4) S2 (2 - theta)/(1 - theta)^2;
E+ := (U^2 vmax/4)^2 S2/(2(1 - thetaU)^2). Since U^2 R = thetaU nu^2 and c_1 >= nu^2/4: X(t) <= thetaU c_1 for t <= U,
and X(t) <= theta c_1, tau <= theta for t <= u + 1.

## 8. Proof of (i): Theta^{(n)} has no zero in W_R

Fix w in W_R. By Sec. 4, Theta^{(n)}(w) Mt(n)/gamma(n) = int_0^oo t^{2n} phi(t) G(x(t)) dt =: J.

8.1 Peak bracket. (Q1) psi'(u_R) < 0: psi'(u_R) <= h'(u_R) + 2beta e^{-2u}/(1-beta) + u_R a^2 S1/2 (the tilt term
-(ts/2)g is <= 0 for s >= 0 and <= t a^2 S1/2 for s >= -a^2, 5.7), and -h'(u_R) >= 2 pi e^{2u}(e^{2/n} - 1) >=
4 pi e^{2u}/n = 4/u + 9/n. So Q1 >= (1/u)[4 - 2 beta u e^{-2u}/(1-beta) - u(u+1) a^2/(2n)] > 0 [cert: M1, value
>= 3.9999999 at u_*]. (Q2) psi'(u_L) >= 0: for s <= 0 the tilt is >= 0 and h'(u_L) > 0; for 0 < s <= R,
psi'(u_L) >= h'(u_L) - (u_L R/2) g(X_L), X_L := u_L^2 R/4 (g increasing in X), with h'(u_L) >= (2n/u)(1 - e^{-2Delta})
and (u_L R/2)(S1 + X_L S2/(1 - theta)) <= (2n/u)(theta/4 + theta^2/(16(1 - theta))) (use nu^2 <= n b, nu^4 <= n b^2(b+1),
u u_L <= (u+1)^2, u u_L^3 <= (u+1)^4). So Q2 holds for every n if
1 - e^{-2 Delta} - theta/4 - theta^2/(16(1-theta)) > 0 [cert: Q2, >= 0.012089223]. Hence psi has a unique maximum
point ut in [u_L, u_R].

8.2 Local interval. I := [u1, u2], u1 := ut - 2 eta, u2 := ut + 2 eta. Need eta <= eps_eta [cert: M3a, eta <= eta_up <= 6.8561470e-6
at u_*] and 1/n + 2 eta <= delta_R [cert: M3b], so I is inside [u_L - 2eps_eta, u + 1], where -psi'' >= kappa0 and tau <= theta.
Then psi'(u1) >= 2 eta kappa0 = 2K/eta and psi'(u2) <= -2K/eta (mean value theorem, psi'(ut) = 0).

8.3 Main term. On I, G(x(t)) != 0; its continuous argument satisfies
|d/dt arg G(x(t))| = (t |v|/2) sum_k c_k/((c_k - X)^2 + Y^2) <= (t |v|/2) Sigma (5.8; Im(w G'/G(x)) = v sum c_k/|x + c_k|^2
because s Y = v X). So |arg G(x(t)) - arg G(x(ut))| <= Phi := eta (u + 1) vmax Sigma on I, and
|int_I t^{2n} phi G| >= cos(Phi) int_I t^{2n} phi |G| >= cos(Phi) int_I e^{psi} >= 2 eta cos(Phi)(e^{psi(u1)} + e^{psi(u2)})
(psi increasing on [u1, ut], decreasing on [ut, u2]). Majorant (vmax <= 2a(sqrt R + a), sqrt R <= sqrt(theta) n/(u+1),
Sigma <= (1 + Chat)/n since 1/b <= 1/n and nu^2 <= b(b+1)): Phi <= 2a(1 + Chat) eta_up (sqrt(theta) + a(u+1)/n),
Chat := theta(2 - theta)/(4(1 - theta)^2), eta <= eta_up := sqrt(K/kappa_low), kappa_low := 4(n/u)e^{-2(Delta + 2 eps_eta)}
- a^2/(2n) <= kappa0 (4 pi e^{2u} >= 4n/u, b >= n) [cert: M4, Phi <= 1.2507254e-4].

8.4 Side tails. psi concave and below its tangents: int_0^{u1} e^{psi} <= e^{psi(u1)}/psi'(u1) <= (eta/(2K)) e^{psi(u1)},
int_{u2}^{U} e^{psi} <= (eta/(2K)) e^{psi(u2)}. With Sec. 6: |int_0^{u1}| + |int_{u2}^{U}| <=
(1 + eps0) e^{E+} (eta/(2K))(e^{psi(u1)} + e^{psi(u2)}). E+ <= a^2 rho^4 (theta(u+1)^2/n + a^2(u+1)^4/n^3)/(8(1 - thetaU)^2)
[cert: M5, E+ <= 7.8678191e-9].

8.5 Crude tail (t > U). By (5.9), |int_U^oo| <= (1 + eps0) 2 pi^2 int_U^oo e^{Hc(t)} dt, Hc := h + t^2 Wmax/(4b)
(phi_1 <= 2 pi^2 e^{9t/2 - pi e^{2t}}, |x(t)| <= t^2 Wmax/4). Hc is concave on [U, oo): Hc'' = -2n/t^2 - 4 pi e^{2t}
+ Wmax/(2b) < 0, since 4 pi e^{2t} >= 4n/u and Wmax/(2b) <= (R + vmax)/(2n) <= theta n/(2(u+1)^2) + a(sqrt(theta)/(u+1) + a/n)
[cert: Hc concave, 4n/u minus this >= 3.30213e11 at u_*, increasing in u]; with Hc'(U) < 0 (below) the tangent at U
gives <= (1 + eps0) 2 pi^2 e^{Hc(U)}/|Hc'(U)|. Lower bound:
e^{psi(u2)} >= 2 pi^2 (1 - beta e^{-2 u_L}) e^{h(u2) + LG}, LG := -X2 S1 - X2^2 S2/(2(1 - theta)), X2 := theta nu^2/4 (5.6).
Estimates (from the h' identities): h(U) - h(u) <= -pi e^{2u}(e^{2y} - 1 - 2y), y := U - u = (rho - 1)u + rho;
h(u) - h(u2) <= max(n Delta^2/(u_L u) + 2 pi e^{2u}(Delta - (1 - e^{-2Delta})/2), n e^2/u^2 + pi e^{2u}(e^{2e} - 1 - 2e)),
e := 1/n + 2 eta (h concave, u2 in [u_L, u + e]); |Hc'(U)| >= 2 pi e^{2u}(e^{2y} - 1) - U Wmax/(2b) >= Hcp_low :=
(n/u)(2(e^{2 rho} - 1) - rho theta/2) - a rho(sqrt(theta) + a) >= 1.13607e12 [cert: M6 Hcp_low] (y >= rho;
U Wmax/(2b) <= (u+1) rho (R + vmax)/(2n) <= rho theta n/(2u) + a rho(sqrt(theta) + a)). Hence
|int_U^oo| <= eta e^{psi(u2)} e^{logratio}, logratio := log(1 + eps0) - log(1 - beta e^{-2u_L}) + log(5n/K)/2 - c6 n,
using -log eta <= log(4 pi e^{2u}/K)/2 <= log(5n/K)/2 and the normalized exponent bound
[Hc(U) - h(u2) - LG]/n <= -F(u) + G_up(u) =: -c6(u), F := (e^{2y} - 1 - 2y)/u, G_up := thetaU/4 +
a rho^2 (u+1)(sqrt(theta) + a(u+1)/n)/(2n) + dh1 + dh2 + theta/4 + theta^2/(32(1-theta)), dh1, dh2 the two branches above
divided by n (with 2 pi e^{2u} = (2n/u)(1 + 9u/(4n)), pi e^{2u} >= n/u, Wmax <= R + vmax, b >= n)
[cert: M6 c6 >= 0.78582827 at u_*; logratio <= -7.84740e11].

8.6 Assembly. |J| >= eta (e^{psi(u1)} + e^{psi(u2)}) [2 cos(Phi) - (1 + eps0) e^{E+}/(2K) - e^{logratio}] > 0
[cert: FINAL, >= 1.8747279]. So Theta^{(n)}(w) != 0. (Positions: 0 < u1 < u2 <= u + 1 < U.)

## 9. Uniformity in n >= n_*

Each majorant is monotone in u = u(n) on [12, oo) (n = Z(u) is increasing, n/u = pi e^{2u} - 9/4 increasing, and
(u+1)^k/n decreasing for k <= 4): M1 increasing (u e^{-2u} and (u+1)/(pi e^{2u} - 9/4) decrease); kappa_low :=
4(n/u)e^{-2(Delta + 2 eps_eta)} - a^2/(2n) <= kappa0 increasing, so eta <= eta_up := sqrt(K/kappa_low) decreasing;
Phi_up, E_up, G_up decreasing; F increasing for u >= 1/(2(rho - 1)) = 11.4947.. < 12 [cert: rho - 1 >= 1/(2u_*)]
(d/du of (e^{2y} - 1 - 2y)/u > 0 when 2(rho-1)u >= 1); hence c6 increasing; Hcp_low increasing; logratio decreasing
in n because c6 > 1/(2n) [cert]; Phi_up <= pi/2 at u_* [cert] and decreasing, so cos(Phi_up) increases; E_up and
logratio decrease; so FINAL increases. Therefore checking at u_* = 12 proves (i) for all n >= n_*. (Hcp_low increasing
and >= 1 gives Hc'(U) < 0 and -log|Hc'(U)| <= 0 for all u >= u_*; the Hc-concavity margin 4n/u - theta n/(2(u+1)^2)
- a(sqrt(theta)/(u+1) + a/n) is increasing.)
(The grid in certify.log at u = 12..300 shows the same monotone values; it is not part of the proof.)

## 10. From (i) to (ii)-(iii); comparisons

(ii) is Sec. 3 with R = R_n. (iii), for a = 0.37338196661320752 (a = 1/2: replace 6.45 by 3.59 < 3.6): D(R) >= R/(4a^2) - 1/4; u(n) + 1 <= (1/2) log n for n >= n_* because
Z(log n/2 - 1) = (log n/2 - 1)(pi n/e^2 - 9/4) >= n iff (log n/2 - 1)(pi/e^2 - 9/(4n)) >= 1, true at n_* (>= 4.448)
and increasing [cert: H2]. So D_n >= 3.6 (n - 1/2)^2/(4a^2 log^2 n) - 1/4 >= 6.45 n^2/log^2 n
(3.6/(4a^2) >= 6.4555883 [cert: H3 constant]; slack 0.0055 n^2 - 6.46 n - log^2(n)/4 > 0 [cert: H3 slack]);
a = 1/2: 3.6 > 3.59. Power forms (c > 1/2, C > 0): n >= C d^c with n >= n_* gives d <= (n/C)^{1/c} <= 6.45 n^2/log^2 n iff
C >= f(n) := n^{1-2c}(log^2 n/6.45)^c; log f is concave in L = log n with maximum at L = 2c/(2c-1), so
C_c := f(max(n_*, e^{2c/(2c-1)})) [cert: power forms].
Two-thirds wedge (n^2 log(n+2) >= 24 d^3): for d > D_all it forces n >= 3.41725e37 (> n_*) and
d <= (n^2 log(n+2)/24)^{1/3} <= 6.45 n^2/log^2 n because log(n+2) log^6 n <= 24 (6.45)^3 n^4 (24 (6.45)^3 = 6440.067 exactly
used in Arb) for every n >= 8: the ratio n^4/(log(n+2) log^6 n) has log-derivative 4 - n/((n+2) log(n+2)) - 6/log n
>= 4 - 1/log 10 - 6/log 8 > 0 for n >= 8, and the inequality holds at n = 8 [cert: containment]; so for every d > D_all our
region contains it; at d = D_all + 1 ours needs n >= 4.69742e13 (rounded up) [cert]. Crossover with the all-n range:
D_n > D_all for n >= 4.6974173e13 (rounded up). D_all = 16,140,857,488,599,731,839,205,676.
Self-contained (no all-n range): put d0 := (n_*^2 log(n_* + 2)/24)^{1/3} <= 1.04710e8 (rounded up). If d > d0 and
n^2 log(n+2) >= 24 d^3, then n > n_* (n^2 log(n+2) is increasing), so d <= 6.45 n^2/log^2 n by the containment inequality:
for every d > d0 our region contains the two-thirds region [cert: self-contained].

## 11. Constants (DERIVATION -> certify.py)

| quantity | value | script |
|---|---|---|
| theta, thetaU, delta_R, Delta, K, eps_eta, u_* | 9/10, 49/50, 1, 17/25, 4, 1/1000, 12 | theta, thetaU, dR, Delta, K, ee, u_star |
| eps0 | <= 0.0021760612 | eps0 |
| Q2 slack | >= 0.012089223 | Q2 |
| rho | >= 1.043498389 | rho |
| n_* | 998,616,377,771 | nstar_int |
| M1 (= u Q1 lower bound) | >= 3.9999999 | M1 |
| eta_up | <= 6.8561470e-6 | eta_up |
| Phi_up | <= 1.2507254e-4 | Phi_up |
| E_up | <= 7.8678191e-9 | E_up |
| F(u_*), G_up(u_*), c6 | >= 1.5637456, <= 0.77791741, >= 0.78582827 | F, G, c6 |
| Hcp_low | >= 1.13607e12 | Hcp_low |
| Hc concavity margin | >= 3.30213e11 | hcc |
| logratio(u_*) | <= -7.84740e11 | logratio |
| FINAL | >= 1.8747279 | final |
| 4a^2; 3.6/(4a^2) | 0.55765637196778566057; >= 6.4555883 | four_a2, kappa |
| H2 | >= 4.4484549 | H2 |
| crossover with D_all | n >= 4.6974173e13 | hi_n |
| two-thirds threshold at D_all + 1 | n >= 3.41725e37 (lower bound) | lo_n |
| ours at D_all + 1 | n >= 4.69742e13 (rounded up) | hi2 |
| 24 * 6.45^3 | 6440.067 | C23 |
| d0 (self-contained comparison) | <= 1.04710e8 (rounded up) | d0 |
| power forms C_c (c = .51, .55, .6, 1) | 7.6890, 0.87160, 0.069829, 1.1853e-10 (rounded up) | Cc |

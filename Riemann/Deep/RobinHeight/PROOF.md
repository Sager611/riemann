# PROOF: height-tied ranges for Nicolas' and Robin's inequalities

Agent robinh, 2026-11-05. Companion files: REPORT.md, certify.py (analytic ranges), smallcheck.py (small range).
Notation: p runs over primes; theta, psi are Chebyshev's functions; eps(t) := theta(t) - t;
h(t) := 1/log t + 1/log^2 t,  w(t) := h(t)/t^2 = -(d/dt) (1/(t log t))  (t > 1);
P(x) := prod_{p<=x} (1-1/p)^{-1};  S_2(x) := sum_{p>x} sum_{m>=2} 1/(m p^m) in [0, 1/(2(x-1))];
n(x) := sqrt(x) log x (normalisation used by the certificate); T0 := 3*10^12; gamma = Euler's constant.

## 0. Inputs

| # | statement | source |
|---|---|---|
| H1 | every zero rho = beta + i gamma with 0 < gamma <= T0 has beta = 1/2 | Platt-Trudgian, Bull. LMS 53 (2021), arXiv 2004.09765, abstract/main theorem |
| H2 | zeta(s) != 0 for Re s > b0, b0 = 7/8 (variant b0 = 0.87338196661320752); by the functional equation also none with Re s < 1 - b0 | OpenAI 2026 (Lean-checked, sorry-free port in Riemann/Deep/QRH); variant: Bebel (Lean-checked per QuasiRiemannTracker, axioms unchecked) |
| H3 | abs(N(t) - (t/2pi) log(t/(2 pi e))) <= R(t) := 0.1038 log t + 0.2573 log log t + 9.3675 for t >= e | Hasanalizade-Shen-Wong, J. Number Theory 235 (2022), arXiv 2107.06506, Corollary 1.2 (label main-thm) |
| H4 | sum over all nontrivial zeros of 1/(rho(1-rho)) = 2 + gamma - log 4pi = 0.0461914179... | classical (Davenport, Ch. 12; also Nicolas 2012, arXiv 1202.0729, eq. (beta)) |
| H5 | psi_1(x) := int_0^x psi = x^2/2 - sum_rho x^{rho+1}/(rho(rho+1)) - x log 2pi + zeta'(-1)/zeta(-1) - sum_{r>=1} x^{1-2r}/(2r(2r-1)), x >= 1, the rho-series absolutely convergent | Ingham 1932, Thm 28 |
| H6 | abs(x - psi(x)) <= 0.94 sqrt x (11 < x <= 1e19); 0.05 sqrt x < x - theta(x) <= 1.95 sqrt x (1423 <= x <= 1e19; the lower bound for 1 <= x <= 1e19) | Buthe, Math. Comp. 87 (2018), arXiv 1511.02032, Theorem 2 (t:explicit-bounds) |
| H7 | theta(y) < y log 4 for y >= 1 | Chebyshev/Erdos, Hardy-Wright Thm 415 |
| H8 | Mertens: log P(x) - log log x -> gamma; sum_{p<=x} 1/p = log log x + B + E(x), E(x) -> 0, B = gamma - sum_p sum_{m>=2} 1/(m p^m) | classical (Rosser-Schoenfeld 1962, (4.13)-(4.15)) |
| H9 | Robin's inequality holds for 5041 <= n <= 10^(10^13.11485) | Morrill-Platt, arXiv 1809.10813, Theorem "Platt" (final largest prime 29,996,208,012,611) |

Only H1-H3 involve zeros above height 14; H2 is used only for zeros with abs(gamma) > T0.

## 1. The Mertens identity (Lemma 1)

For x >= 2: log P(x) = gamma + log log x + eps(x)/(x log x) - I(x) - S_2(x), where I(x) := int_x^inf eps(t) w(t) dt.

Proof. Partial summation: sum_{p<=x} 1/p = int_{2-}^x d theta(t)/(t log t) = log log x - log log 2 + 1/log 2 + eps(x)/(x log x) + int_2^x eps w
(eps(2-) = -2). PNT gives eps(t) = O(t/log^3 t), so int_2^inf eps w converges absolutely; comparing with H8 gives
E(x) = eps(x)/(x log x) - I(x). Finally log P(x) = sum_{p<=x} 1/p + sum_{p<=x} sum_{m>=2} 1/(m p^m) = log log x + B + E(x) + (gamma - B - S_2(x)). QED
(This is RS62 (4.13)-(4.15); also eq. (ps) of arXiv 2608.01498.)

## 2. Reductions

**Lemma 2 (Nicolas).** For x >= 3 put D(x) := log log x - log log theta(x) + eps(x)/(x log x). Then D(x) >= 0 and
log P(x) - gamma - log log theta(x) = D(x) - I(x) - S_2(x).
Hence Nicolas' inequality N_k/phi(N_k) > e^gamma log log N_k holds at every k with p_k = x as soon as -I(x) > S_2(x).

Proof. The identity is Lemma 1. Put u = eps(x)/x > -1 and v = log(1+u)/log x > -1 (theta(x) > 1). Then
log log theta(x) = log log x + log(1+v) and D = u/log x - log(1+v) >= u/log x - v = (u - log(1+u))/log x >= 0. Since N_k/phi(N_k) = P(p_k)
and log N_k = theta(p_k), the claim follows. QED

**Upper bound for D (needed for Robin only).** If abs(eps(x)) <= eta x with eta <= 1/2 and log x >= 20, then
D(x) <= eta^2/log x * (1 + (1+eta)^2/log x)   (from log(1+u) >= u - u^2 and -log(1+v) <= -v + v^2 for abs(u), abs(v) <= 1/2).

**Colossally abundant numbers.** For y > 1, k >= 1 let F(y,k) := log(1 + 1/(y + y^2 + ... + y^k))/log y; F decreases in y and in k.
For eps > 0 let x_k = x_k(eps) be the solution of F(x_k, k) = eps and N_eps := prod_p p^{a_p}, a_p := #{k >= 1 : F(p,k) >= eps}.
Since sigma(p^a)/sigma(p^{a-1}) >= p^{1+eps} iff F(p,a) >= eps, N_eps maximises sigma(n)/n^{1+eps}, and a_p >= k iff p <= x_k, so
log N_eps = sum_k theta(x_k). The colossally abundant (CA) numbers are the N_eps.

**Lemma 3 (interval lemma; Robin 1984, Briggs 2006).** Let eps_c be a value where N_eps jumps, N- := N_{eps_c + 0} < N+ := N_{eps_c}.
If Robin's inequality holds (strictly) at N- and N+, it holds for every n in [N-, N+].
Proof. Both N- and N+ maximise sigma(n)/n^{1+eps_c} (ties at eps_c), so log(sigma(n)/n) <= c + eps_c log n for all n with equality at
N- and N+. G(u) := c + eps_c u - log log u is convex for u > 1, so for log N- <= log n <= log N+,
log(sigma(n)/n) - log log log n <= G(log n) <= max(G(log N-), G(log N+)), and G(log N+-) = log(sigma(N+-)/N+-) - log log log N+-. QED
Consequently Robin holds for all n in [N_{eps_1}, N_{eps_2}] (eps_2 < eps_1) once it holds at every N_eps, eps in [eps_2, eps_1].

**Lemma 4 (Robin at N_eps).** Let x = x_1(eps) >= 10^10, x_2 = x_2(eps). Then
log(sigma(N_eps)/N_eps) - gamma - log log log N_eps <= D(x) - I(x) - S_2(x) - sum_{x_2 < p <= x} p^{-2} - L_2,
L_2 := log(1 + log(1 + theta(x_2)/theta(x))/log theta(x)).
Proof. sigma(N)/N = prod_{p<=x} (1 - p^{-a_p-1})/(1 - 1/p) = P(x) prod_{p<=x}(1 - p^{-a_p-1}); for x_2 < p <= x, a_p = 1 and
log(1 - p^{-2}) <= -p^{-2}; the other factors are < 1. log N >= theta(x) + theta(x_2), so
log log log N >= log log theta(x) + L_2. Apply Lemma 2's identity. QED
So Robin holds at N_eps if  -I(x) + D(x) < sum_{x_2<p<=x} p^{-2} + L_2 + S_2(x).

The two criteria sandwich the same quantity:  S_2(x) < -I(x)  (Nicolas)  and  -I(x) + D(x) < sum p^{-2} + L_2 + S_2  (Robin).
Under RH, -I(x) n(x) -> 2 + O(0.0462), sum p^{-2} n(x) -> sqrt 2, L_2 n(x) -> sqrt 2: margins 2 - 0.046 and 2 sqrt 2 - 2 - 0.046 = 0.78.

## 3. Explicit formula for the tail (Lemma 5)

For x >= 2, assuming H2 with b0 < 1:
-I(x) = M(x) + Z(x) + log(2pi)/(x log x) + V(x),
M(x) := int_x^inf (psi(t) - theta(t)) w(t) dt,  Z(x) := sum_rho (1/rho) F(rho),  F(rho) := int_x^inf t^{rho-2} h(t) dt,
V(x) := (1/2) int_x^inf log(1 - t^{-2}) w(t) dt in [-x^{-3}/log x, 0].

Proof. -I = int (t - theta) w = int (t - psi) w + M. Let Dd(t) := psi_1(t) - t^2/2 + t log 2pi - zeta'(-1)/zeta(-1) + q(t),
q(t) := sum_r t^{1-2r}/(2r(2r-1)). By H5, Dd(t) = -sum_rho t^{rho+1}/(rho(rho+1)), and Dd'(t) = psi(t) - t + log 2pi + (1/2) log(1 - t^{-2})
for t not a prime power (q' = (1/2) log(1 - t^{-2})). Dd is continuous and piecewise C^1, abs(Dd(t)) <= t^{1+b0} sum 1/abs(rho)^2, and
w(t) = O(1/(t^2 log t)), so Dd w -> 0 and integration by parts gives
int_x^inf (psi - t + log 2pi + (1/2) log(1-t^{-2})) w = -Dd(x) w(x) - int_x^inf Dd w'.
abs(w'(t)) = (2 + 3/L + 2/L^2)/(t^3 L) <= 11/(t^3 log t) for t >= 2 (L = log t), so sum_rho int_x^inf abs(t^{rho+1} w'(t)/(rho(rho+1))) dt <= sum 1/abs(rho)^2 * int_x^inf 11 t^{b0-2}/log t dt < inf:
we may integrate termwise (dominated convergence). For each rho, x^{rho+1} w(x) + int_x^inf t^{rho+1} w' = -(rho+1) int_x^inf t^rho w
(integration by parts, t^{rho+1} w -> 0 as beta < 1). Hence int_x^inf (psi - t) w = -Z(x) - log(2pi)/(x log x) - V(x), using int_x^inf w = 1/(x log x).
Finally abs(log(1-t^{-2}))/2 <= t^{-2} (t >= 2), so 0 >= V >= -x^{-2} int_x^inf w = -x^{-3}/log x. QED
Signs: t - psi(t) = sum_rho t^rho/rho + log 2pi + (1/2) log(1 - t^{-2}) (von Mangoldt), so the constant enters -I with + sign and the
trivial zeros give V <= 0. Remark (b0 = 1, no half-plane): the interchange is still valid by Tonelli, because
int_x^inf t^{beta+1} abs(w') <= 11 x^{beta-1}/((1-beta) log x) and sum_rho 1/((1-beta) abs(rho)^2) < inf by the classical zero-free region.

**Lemma 6 (one zero).** For 0 < beta < 1 and x > e: abs(F(rho)) <= x^{beta-1} h(x) (1 + kappa(x)/abs(rho)) / abs(rho - 1),
kappa(x) := 1 + 2 x abs(h'(x))/h(x) = 1 + 2(L+2)/(L(L+1)), L = log x.
Proof. h > 0, h' < 0, h'' > 0, h(t) -> 0, t h'(t) -> 0. Integrating by parts,
F = -x^{rho-1} h(x)/(rho-1) - G/(rho-1),  G := int_x^inf t^{rho-1} h'(t) dt = -x^rho h'(x)/rho - (1/rho) int_x^inf t^rho h''(t) dt.
Since t^{beta-1} <= x^{beta-1} on [x, inf) and t h'' > 0: abs(int t^rho h'') <= x^{beta-1} int_x^inf t h'' = x^{beta-1}(h(x) + x abs(h'(x))).
So abs(G) <= x^{beta-1}(h(x) + 2x abs(h'(x)))/abs(rho). QED

**Lemma 7 (zero sum).** Under H1-H4, for x > e:
abs(Z(x)) <= x^{-1/2} h(x) beta_0 (1 + kappa/gamma_1) + (1/2)(x^{b0-1} + x^{-b0}) h(x) (1 + kappa/T0) S_T,
beta_0 := 2 + gamma - log 4pi, gamma_1 := 14.1347 (<= least positive ordinate),
S_T := (log(T0/2pi) + 1)/(pi T0) + 2R(T0)/T0^2 + (c1+c2+c3)(2 log T0 + 1)/T0^2 >= sum_{abs(gamma)>T0} gamma^{-2}   (= 2.95941e-12).
Proof. abs(F(rho)/rho) <= x^{beta-1} h (1 + kappa/abs(rho))/abs(rho(rho-1)).
(i) abs(gamma) <= T0: beta = 1/2 (H1), abs(rho(rho-1)) = abs(rho)^2 = rho(1-rho), abs(rho) >= gamma_1; and Re 1/(rho(1-rho)) > 0 for every zero, so by H4
sum_{abs(gamma)<=T0} abs(rho)^{-2} <= beta_0 and sum abs(rho)^{-3} <= beta_0/gamma_1.
(ii) abs(gamma) > T0: 1 - b0 <= beta <= b0 (H2), abs(rho(rho-1)) >= gamma^2, abs(rho) > T0. The involution rho -> 1 - conj(rho) permutes these zeros
(with multiplicity) and fixes gamma, so sum x^{beta-1}/gamma^2 = sum (x^{beta-1} + x^{-beta})/(2 gamma^2) <= (x^{b0-1} + x^{-b0})/2 * sum gamma^{-2}
(beta -> x^{beta-1} + x^{-beta} is convex, maximal at the ends).
(iii) sum_{gamma>T} gamma^{-2} = -N(T)/T^2 + 2 int_T^inf N t^{-3} <= int_T^inf M'(t) t^{-2} dt + R(T)/T^2 + 2 int_T^inf R t^{-3},
M(t) = (t/2pi) log(t/(2 pi e)), M'(t) = log(t/2pi)/(2pi), and R(t) <= (c1+c2+c3) log t for t >= e^e. QED

**Lemma 8 (prime powers).** For 2.03e6 <= x <= 1e38, with I_a(x) := int_x^inf t^{a-2} h = a E_1((1-a)L) + e^{-(1-a)L}/L and
J_a(x) := int_x^inf t^{a-2} h log t = e^{-(1-a)L}/(1-a) + E_1((1-a)L), Y := 1e38:
M(x) >= I_{1/2}(x) - I_{1/2}(Y) - 1.95 (I_{1/4}(x) - I_{1/4}(Y)),
M(x) <= I_{1/2}(x) + I_{1/3}(x) + I_{1/4}(x) + J_{1/5}(x)/log 2 + (log 4 - 1)(I_{1/2}(Y) + I_{1/3}(Y)) + log 4 * J_{1/4}(Y)/log 2.
Proof. psi(t) - theta(t) = sum_{k>=2} theta(t^{1/k}), the terms vanishing for t^{1/k} < 2. Lower: drop t > Y and k >= 3, and use
theta(sqrt t) >= sqrt t - 1.95 t^{1/4} (H6, 1423 <= sqrt t <= 1e19). Upper: for t <= Y every t^{1/k} <= 1e19, so theta(t^{1/k}) < t^{1/k} (H6) and
psi - theta <= sqrt t + t^{1/3} + t^{1/4} + (log t/log 2) t^{1/5}; for t > Y use H7: psi - theta <= log 4 (sqrt t + t^{1/3} + (log t/log 2) t^{1/4}). QED

**Lemma 9 (Robin's positive terms).** For 1423 <= x_2 < x <= 1e19:
sum_{x_2<p<=x} p^{-2} >= E_1(log x_2) - E_1(log x) - 1.95 x^{-3/2}/log x - 3.9 x_2^{-3/2}/log x_2;
for x > 1e19 the same with x replaced by 1e19 (drop the primes above 1e19).
Proof. With g = 1/(t^2 log t): sum = int_{x_2}^x g dt + eps(x) g(x) - eps(x_2) g(x_2) + int eps abs(g'); eps(x_2) < 0 (H6), eps >= -1.95 sqrt t,
abs(g') <= 3/(t^3 log t). QED
And L_2 >= log(1 + log(1 + th2/th1)/log th1) with th2 = x_2 - 1.95 sqrt(x_2) <= theta(x_2) and th1 >= theta(x) (th1 = x for x <= 1e19 by H6;
th1 = x(1 + eta) above, Lemma 10). L_2 increases in theta(x_2) and decreases in theta(x).

**Lemma 10 (pointwise bound, second-order averaging).** Let x >= 10, 0 < delta <= 1e-3, h = delta x, and
Lam(s) := min(s, 2h - s) on [0, 2h]. Then abs(psi(x) - x) <= h + log 2pi + x^{-2} + Zs(x), where
Zs(x) := x^{1/2}(1+2 delta)^{5/2} [S1(U) + 4 delta^{-2} S3(U)] + (1/2)(x^{b0} + x^{1-b0})(1+2 delta)^{b0+2} 4 delta^{-2} S3(T0),  U := 2/delta,
S1(U) >= sum_{0<abs(gamma)<=U} 1/abs(gamma), S3(T) >= sum_{abs(gamma)>T} abs(gamma)^{-3} (explicit from H3, see certify.py).
Then abs(theta(x) - x) <= that + sqrt x + x^{1/3} + x^{1/4} + (log x/log 2) x^{1/5} =: eta(x) x  (x <= 1e38).
Proof. psi is nondecreasing and int Lam = h^2, so psi(x) h^2 <= int_x^{x+2h} psi(t) Lam(t-x) dt = Delta^2 psi_2(x) and
psi(x) h^2 >= Delta^2 psi_2(x - 2h), where psi_2 = int psi_1 and Delta^2 f(y) = f(y+2h) - 2f(y+h) + f(y). Integrating H5 termwise (uniform
convergence on compacts) and using Delta^2[y^3/6] = h^2 (y + h), Delta^2[y^2/2] = h^2, Delta^2[linear] = 0,
Delta^2 psi_2(y)/h^2 = y + h - log 2pi - h^{-2} sum_rho Delta^2[t^{rho+2}](y)/(rho(rho+1)(rho+2)) - h^{-2} int q'(t) Lam(t-y) dt,
with 0 <= -h^{-2} int q' Lam <= y^{-2}. Per zero: Delta^2[t^{rho+2}](y) = (rho+1)(rho+2) int t^rho Lam, = (rho+2) int_y^{y+h}((t+h)^{rho+1} - t^{rho+1}) dt,
and abs(Delta^2[t^{rho+2}]) <= 4(y+2h)^{beta+2}, so its contribution is <= x^beta (1+2delta)^{beta+2} min(1/abs(gamma), 2/(delta gamma^2), 4/(delta^2 abs(gamma)^3))
(y + 2h <= x(1 + 2 delta)). Split at U (abs(gamma) <= U: 1/abs(gamma); U < abs(gamma): 4/(delta^2 abs(gamma)^3), beta = 1/2 up to T0) and
pair the zeros above T0 as in Lemma 7 (U = 2/delta < T0 for the deltas used). psi - theta as in Lemma 8. QED
(The backward average's trivial-zero correction has the helpful sign; bounding it by y^{-2} is merely wasteful.)
With delta = 1e-9: eta(1e19) ~ 2.6e-8, eta(1e31) ~ 1.1e-9 (cross-check: FKS 2206.12557 table theta-med gives 2.2e-12 for log x >= 71).

## 4. Theorems

**Theorem A (Nicolas).** Under H1-H3 and H2 with b0 = 7/8: N_k/phi(N_k) > e^gamma log log N_k for every k with p_k <= X_N, X_N = 10^32.287
(certify.py; variant b0 = 0.87338196661320752 and b0 = 1 in REPORT.md).
Proof. p_k <= 9,999,991: direct interval computation (smallcheck.py, all 664,579 primorials). 1e7 <= x <= X_N: by Lemmas 2, 5, 7, 8,
-I(x) - S_2(x) >= M_low - Zup - x^{-3}/log x - 1/(2(x-1)) > 0, checked in ball arithmetic on boxes in log x covering [log 1e7, log X_N]. QED

**Theorem B (Robin).** Under the same inputs and H9: sigma(n) < e^gamma n log log n for every n with 5041 <= n <= N_{eps_R},
eps_R := F(X_R, 1), X_R = 10^31.194; since log N_{eps_R} >= theta(X_R) >= X_R(1 - 2e-9), this covers all n with
log n <= 0.999999 X_R, i.e. n <= 10^(10^30.832).
Proof. CA numbers N_eps with x_1(eps) < 1e10 satisfy log N_eps < 1.01e10 < 2.9997e13, so they and all n up to 10^(10^13.11485) are covered by H9.
For x_1(eps) in [1e10, X_R], Lemma 4 with Lemmas 5, 7, 8, 9, 10 and the D bound gives Robin at N_eps (certify.py, ball arithmetic, boxes in
log x; x_2 enclosed rigorously by bisection on F(., 2)). Lemma 3 fills the gaps between consecutive CA numbers. QED

**What limits X.** At the end of the range the only growing term is the zeros above T0:
(x^{b0-1/2}/2)(1 + 1/log x) S_T against the margin 2 - 0.05 (Nicolas) or 0.78 (Robin), i.e. x^{3/8} ~ c T0/log T0.
b0 = 1 (no half-plane) gives x^{1/2} ~ c T0/log T0 instead.


## 5. Level 2: smoothed cutoff (certify2.py)

The sharp cutoff at t = x is what makes the zeros above T0 cost 1/gamma^2. Smoothing the cutoff over a tiny window costs only a
pointwise bound for theta on that window, and the zeros then cost 1/(lam abs(gamma)^3).

**Lemma 11 (ramp).** Let x >= 1e7, 0 < lam <= 1e-3, phi(t) := min(max((log t - log x)/lam, 0), 1), omega := w phi, H(u) := 1/u + 1/u^2.
Assume abs(theta(t) - t) <= eta t on [x, x e^lam]. Then
-I(x) = M_phi + Z_phi + log(2pi) int omega + V_phi + W,
M_phi := int (psi - theta) omega in [M(x e^lam), M(x)],  Z_phi := sum_rho (1/rho) F_phi(rho), F_phi(rho) := int t^{rho-2} h(t) phi(t) dt,
0 <= int omega <= 1/(x log x),  V_phi in [-x^{-3}/log x, 0],  abs(W) <= eta H(log x) lam/2,  W := int_x^{x e^lam} (t - theta) w (1 - phi).
For every zero:
(B1) abs(F_phi(rho)) <= x^{beta-1} h(x) [(1 + kappa/abs(rho))/abs(rho-1) + lam/2];
(B2) abs(F_phi(rho)) <= x^{beta-1} (2H(L)/lam + 3 abs(H'(L)))/abs(rho-1)^2,  L = log x.
Proof. The decomposition is the identity 1_{t>=x} = phi + 1_{[x, x e^lam]}(1 - phi); in u = log t, t w dt = H(u) du, so
int_x^{x e^lam} t w (1-phi) dt = int_L^{L+lam} H(u)(1 - (u-L)/lam) du <= H(L) lam/2. omega is continuous, piecewise C^1, omega(x) = 0,
so Lemma 5's proof applies verbatim with w replaced by omega (no boundary term at x). (B1): F_phi = F - int_L^{L+lam} e^{su} H (1 - Phi) du
with s = rho - 1, and Lemma 6. (B2): with G(u) := H(u) Phi((u-L)/lam), G(L) = 0, two integrations by parts on (L, L+lam) and (L+lam, inf) give
int_L^inf e^{su} G du = -(1/s)[ (1/(s lam))(e^{s(L+lam)} H(L+lam) - e^{sL} H(L)) - (1/s) int_L^inf e^{su} G'' du ],
and int_L^inf abs(G'') <= int_L^inf H'' + 2 abs(H'(L)) = 3 abs(H'(L)) (H'' > 0, abs(H') decreasing; on the window G'' = H''(u-L)/lam + 2H'/lam). QED
Zero sum: abs(Z_phi) <= x^{-1/2} h(x)[beta_0(1 + kappa/gamma_1) + (lam/2) S1(T0)] + (1/2)(x^{b0-1} + x^{-b0}) H(L)(2/lam + 3 abs(H')/H) S3(T0)
((B1) below T0, (B2) above, pairing as in Lemma 7; S1(T0) = 119.99 >= sum_{0<abs(gamma)<=T0} 1/abs(gamma), S3(T0) = 4.844e-25).
Normalised, the last two error terms are about x^{b0-1/2} S3/lam (zeros) and sqrt(x) eta lam/2 (window); lam is chosen per box near
sqrt(2 x^{b0-1} S3/eta).

**Lemma 10_k (order-k one-sided averages).** For k >= 2, 0 < delta <= 1e-3/k, h = delta x, B := (1_{[0,h]})^{*k} (mass h^k, mean kh/2),
psi_k the k-fold integral of psi: h^k psi(y) <= Delta_h^k psi_k(y) <= h^k psi(y + kh). From H5 integrated termwise,
Delta^k psi_k(y)/h^k = y + kh/2 - log 2pi - h^{-k} sum_rho Delta^k[t^{rho+k}](y)/(rho(rho+1)...(rho+k)) + [0, y^{-2}],
and each zero contributes at most x^beta (1+k delta)^{beta+k} min(1/abs(gamma), (2/delta)^k/abs(gamma)^{k+1}). Hence
abs(psi(x) - x) <= k delta x/2 + log 2pi + x^{-2} + x^{1/2}(1+k delta)^{1/2+k}[S1(U) + (2/delta)^k S_{k+1}(U)]
                    + (1/2)(x^{b0} + x^{1-b0})(1+k delta)^{b0+k}(2/delta)^k S_{k+1}(T0),   U = 2/delta,
S_m(T) := 2[(log(T/2pi)/((m-1)T^{m-1}) + 1/((m-1)^2 T^{m-1}))/(2pi) + R(T)/T^m + (c1+c2+c3)(m log T + 1)/(m T^m)] >= sum_{abs(gamma)>T} abs(gamma)^{-m}.
(Zeros with T0 < abs(gamma) <= U, if any, are majorised by the last term, so U may exceed T0.) For x <= 1e38, psi - theta <=
sqrt x + x^{1/3} + x^{1/4} + (log x/log 2) x^{1/5}; for 1e38 < x <= 1e76 (and whenever the window [x, x e^lam] reaches beyond 1e38) the same times (1 + eta_1), eta_1 := eta(1e19) = 2.73e-8
(order 2, delta = 1e-9; eta is decreasing in x for fixed k, delta). certify2.py takes the best of k in {2,...,5} and delta on a grid
10^{j/4}, and Buthe's 1.95/sqrt x for x <= 1e19. Typical values (b0 = 7/8): eta(1e30) = 5.6e-11, eta(1e38) = 4.2e-11.

**Prime powers up to 1e76.** For 1e19 <= y <= 1e38, abs(theta(y) - y) <= eta_1 y. Lemma 8 is extended by a second zone
1e38 < t <= 1e76 with theta(sqrt t) >= (1 - eta_1) sqrt t and psi - theta <= (1 + eta_1)(sqrt t + t^{1/3} + t^{1/4} + (log t/log 2) t^{1/5}),
and theta(y) < y log 4 beyond. Robin's prime-square sum gets a second piece on (1e19, min(x, 1e38)] with error
eta_1 (1/(b log b) + 4/(a log a)), and theta(x_2) >= x_2(1 - eta_1) when x_2 > 1e19.

**Theorem A2 (Nicolas, Level 2).** Under H1-H3 and H2 with b0 = 7/8: Nicolas' inequality holds for every k with p_k <= 1.04*10^40
(b0 = 0.87338196661320752: 1.31*10^40; no half-plane, b0 = 1: 1.42*10^34). Proof: smallcheck below 1e7, and for 1e7 <= x <= X:
-I(x) - S_2(x) >= M_low(x e^lam) - abs(Z_phi) - abs(W) - x^{-3}/log x - 1/(2(x-1)) > 0 (certify2.py).

**Theorem B2 (Robin, Level 2).** Same inputs and H9: Robin's inequality holds for 5041 <= n <= N_{eps_R}, eps_R = F(X_R, 1), X_R = 9.37*10^38,
hence for all n <= 10^(10^38.6095) (b0 = 0.87338196661320752: X_R = 1.19*10^39, n <= 10^(10^38.7133); b0 = 1: X_R = 2.18*10^33,
n <= 10^(10^32.9762)). Proof: Lemma 4 with -I(x) <= M_up(x) + abs(Z_phi) + log(2pi)/(x log x) + abs(W), D(x) <= eta^2/log x (1 + (1+eta)^2/log x)
with eta from Lemma 10_k, Robin's positive terms as extended above, Lemma 3, and H9 below 1e10.

**What limits X at Level 2.** The window term sqrt(x) eta lam/2 against the zeros term x^{b0-1/2} S3/lam: their optimum is
2 sqrt(x^{b0} S3 eta/2), so X ~ (c/(S3 eta))^{1/b0}. eta (~4e-11, the bias k delta/2 of Lemma 10_k) is now the bottleneck; a
Chirre-Helfgott-quality bound (eta ~ pi/T0) would add roughly a factor 50.

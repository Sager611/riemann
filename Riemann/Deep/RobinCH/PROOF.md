# PROOF (RobinCH): Chirre-Helfgott pointwise bound, k-box ramps, sharper D

Agent robinch, 2026-11-05. Base: Riemann/Deep/RobinHeight/PROOF.md (inputs H1-H9, Lemmas 1-9 and 11, Theorems A2/B2;
notation as there: L = log x, H(u) = 1/u + 1/u^2, h(t) = H(log t), w = h/t^2, n(x) = sqrt(x) log x, T0 = 3e12).
Only what changes is proved here. Lemmas 10 and 10_k (our own pointwise bounds) are no longer used anywhere.

## 0. New input

H10 (Chirre-Helfgott, arXiv 2512.15709v1, 17 Dec 2025, preprint, 52 pp.; Corollary cor:psiexpl, TeX l.402-406):
"For any $x\geq 1$, $|\psi(x)-x|\leq \frac{\pi}{3\cdot 10^{12}}\cdot x + 113.67 \sqrt{x}$."
- Source of the corollary (TeX l.4981-4991): for x >= 3e12 + 3 it is Corollary cor:psi (= Proposition prop:sagaro, sigma = 0,
  TeX l.4868-4880) with T = 3e12 + 1 + pi/3, using 1 <= (pi/T) coth(pi/T) <= 1 + pi^2/(3T^2) and
  pi^2/(3T^2) + pi/(T-1) <= pi/(T - 1 - pi/3) = pi/(3e12); below 3e12 + 3 their brute-force Lemma lem:pernic.
- Hypothesis of Prop. prop:sagaro: "the Riemann hypothesis holds up to height T >= 10^7", defined (TeX l.387-388) as
  "all non-trivial zeros of zeta(s) with 0 < Im s <= T have real part 1/2"; range x > max(T, 10^9). No zero-free region,
  no zero-density input. H1 (Platt-Trudgian: RH to 3,000,175,332,800) implies it for T = 3e12 + 1 + pi/3.
- C_T = (1/2pi) log^2(T/2pi) - (1/6pi) log(T/2pi) = 113.66888... < 113.67 at T = 3e12 + 1 + pi/3 (rechecked in arb, check_all2.py (0)).
- psi(x) = sum_{n<=x} Lambda(n) for all real x (their Lemma lem:pernic quotes extrema "at x = 2^-"), so it holds at every real x.
We use H10 only for x >= 1e19 e^{-1e-3} > 3e12 + 3 (Buthe's H6 is used whenever a box and its window lie below 1e19), i.e. only
the part that rests on Prop. prop:sagaro.

## 1. Lemma 12 (pointwise theta; replaces Lemmas 10 and 10_k)

For 3e12 + 3 <= x <= 1e76: abs(theta(x) - x) <= eta_CH(x) x, where
eta_CH(x) := pi/(3e12) + 113.67/sqrt(x) + P(x)/x,
P(x) := (1 + eta_1 [x > 1e38]) (sqrt x + x^{1/3} + x^{1/4} + (log x/log 2) x^{1/5}),  eta_1 := eta_CH(1e19) <= 3.62632e-8.
For each fixed choice of the factor (1 or 1 + eta_1) the formula eta_CH is decreasing in x; at x = 1e38 the factor switches
and eta_CH jumps up (by about 3.6e-27), so the factor-(1 + eta_1) formula, which is >= both, is used for any box or window that
reaches beyond 1e38.
Proof. psi - theta = sum_{k>=2} theta(x^{1/k}), terms with x^{1/k} < 2 vanish, so at most log x/log 2 terms have k >= 5.
For x <= 1e38 every x^{1/k} (k >= 2) is <= 1e19, where theta(y) < y (H6): 0 <= psi - theta <= P(x) (no eta_1).
Then theta - x <= psi - x <= pi x/3e12 + 113.67 sqrt x (H10) and x - theta = (x - psi) + (psi - theta) <= pi x/3e12 + 113.67 sqrt x + P(x).
This proves the claim on [3e12 + 3, 1e38]; by monotonicity abs(theta(y) - y) <= eta_1 y on [1e19, 1e38]. For 1e38 < x <= 1e76 each
x^{1/k} <= 1e38, so theta(x^{1/k}) <= (1 + eta_1) x^{1/k}, giving P(x) with the factor (1 + eta_1); that formula also bounds
abs(theta(x) - x)/x for x <= 1e38 (it is larger). Monotonicity of each fixed-factor formula: pi/3e12 is constant,
x^{-1/2}, x^{-2/3}, x^{-3/4} decrease, and (log x) x^{-4/5} decreases for x > e^{5/4}. QED
Use: eta := min(1.95/sqrt x [H6; only if the box and its window lie below 1e19], eta_CH(x)) wherever RobinHeight used
Lemma 10_k (window term, D(x), theta(x) <= x(1 + eta)); eta_1 (theta on [1e19, 1e38]: the zone 1e38 < t <= 1e76 of the
prime-power integral, Robin's prime-square sum on (1e19, 1e38], theta(x_2) >= x_2(1 - eta_1)) is eta_CH(1e19) instead of
Lemma 10's 2.70e-8 (b0 = 0.87338...; 2.73e-8 for 7/8). A box or window that reaches beyond 1e38 uses the factor-(1 + eta_1) formula throughout; with the
factor fixed, eta_CH decreases, so its value at the left end of a window bounds it on the whole window (certify3.best_eta3 at L - c).

## 2. Lemma 11_k (k-box ramp, one-sided or centred)

Let k >= 1, 0 < lam <= 1e-3, a := lam/k, V := U_1 + ... + U_k with U_i i.i.d. uniform on [0, a], Phi its CDF (continuous,
piecewise polynomial, C^{k-1}); E V = lam/2, Var V = lam^2/(12k). Let c in {0, lam/2}, L' := L - c, x' := e^{L'},
phi(t) := Phi(log t - L'), omega := w phi. (k = 1, c = 0 is the linear ramp of Lemma 11.) Then for x >= 1e7:
-I(x) = M_phi + Z_phi + log(2pi) int omega + V_phi + W,   W := int_0^inf (t - theta(t)) w(t) (1_{t>=x} - phi(t)) dt,
M_phi := int (psi - theta) omega in [M(x' e^lam), M(x')],  0 <= int omega <= 1/(x' L'),  V_phi in [-x'^{-3}/L', 0],
Z_phi := sum_rho F_phi(rho)/rho,  F_phi(rho) := int t^{rho-2} h(t) phi(t) dt = int_{L'}^inf e^{su} H(u) Phi(u - L') du  (s := rho - 1),
abs(W) <= eta_win H(L') m_c,   m_c := int abs(1_{u>=L} - Phi(u - L')) du = E abs(V - c),   eta_win >= sup over [x', x' e^lam] of abs(theta(t)/t - 1),
m_0 = lam/2,  m_{lam/2} = (lam/k) MAD_k <= sd V = lam/sqrt(12k),
MAD_k := E abs(S - k/2) = (2/(k+1)!) sum_{0<=j<=k/2} (-1)^j C(k,j) (k/2 - j)^{k+1}  (S = sum of k uniforms on [0,1]; exact rational),
and for every zero rho = beta + i gamma (0 < beta < 1):
(B1')  abs(F_phi(rho) - F(rho)) <= x'^{beta-1} H(L') m_c,   F(rho) := int_x^inf t^{rho-2} h dt;
(B2_k) abs(F_phi(rho)) <= x'^{beta-1} [ H(L') (2k/(lam abs(s)))^k/abs(s) + abs(H'(L')) lam/(2 abs(s)) + 2 abs(H'(L'))/abs(s)^2 ].
Proof. (i) Decomposition: 1_{t>=x} = phi + (1_{t>=x} - phi) for every t. omega is continuous and piecewise C^1, vanishes on
(0, x'] and omega' = O(t^{-3}) on [x', inf), so the proof of Lemma 5 applies verbatim with w replaced by omega (no boundary
term at x', termwise integration by dominated convergence, the per-zero integration by parts uses beta < 1). This gives the
identity for int (t - theta) omega; the rest is W by definition.
(ii) 0 <= phi <= 1, phi = 0 below x', phi = 1 above x' e^lam, psi - theta >= 0, w > 0: M(x' e^lam) <= M_phi <= M(x');
int omega <= int_{x'}^inf w = 1/(x' L'); abs(log(1 - t^{-2}))/2 <= t^{-2} <= x'^{-2} for t >= x' >= 2 gives V_phi.
(iii) W: t w(t) dt = H(u) du (u = log t), 1_{t>=x} - phi vanishes outside [x', x' e^lam], abs(t - theta) <= eta_win t there and
H is decreasing, so abs(W) <= eta_win H(L') int abs(1_{u>=L} - Phi(u - L')) du. With v = u - L':
int abs(1_{v>=c} - Phi(v)) dv = int_0^c Phi + int_c^inf (1 - Phi) = E (c - V)_+ + E (V - c)_+ = E abs(V - c);
for c = 0 this is E V = lam/2. For c = lam/2 = E V, V = (lam/k) S and E abs(S - k/2) = 2 E (k/2 - S)_+ (S is symmetric about k/2)
= 2 int_0^{k/2} F_S, where F_S(s) = (1/k!) sum_j (-1)^j C(k,j) (s - j)_+^k (Irwin-Hall); integrating termwise gives MAD_k,
and MAD_k <= sd S = sqrt(k/12) (Cauchy-Schwarz; asserted in the code).
(iv) (B1'): F_phi - F = int e^{su} H(u)(Phi(u - L') - 1_{u>=L}) du, and abs(e^{su}) H(u) = e^{(beta-1)u} H(u) <= x'^{beta-1} H(L') for u >= L'.
(v) (B2_k): Phi(u - L') = E 1[V <= u - L']; by Fubini (absolutely convergent since beta < 1) F_phi(rho) = E F(rho; L' + V),
F(rho; y) := int_y^inf e^{su} H(u) du. Two integrations by parts (e^{su} H^{(j)}(u) -> 0 as u -> inf):
F(rho; y) = -e^{sy} H(y)/s + e^{sy} H'(y)/s^2 + s^{-2} int_y^inf e^{su} H''(u) du.
First term: E[e^{s(L'+V)} H(L'+V)] = e^{sL'} ( H(L') E e^{sV} + E[e^{sV}(H(L'+V) - H(L'))] ).
E e^{sV} = ((e^{sa} - 1)/(sa))^k and abs(e^{sa} - 1) <= 1 + e^{(beta-1)a} <= 2 (Re s < 0), so abs(E e^{sV}) <= (2k/(lam abs(s)))^k;
abs(e^{sV}) <= 1 and abs(H(L'+v) - H(L')) <= v abs(H'(L')) (H' < 0 < H''), E V = lam/2.
Second and third terms: abs(H'(L'+v)) <= abs(H'(L')), and int_y^inf abs(e^{su}) H'' <= x'^{beta-1} int_{L'}^inf H'' = x'^{beta-1} abs(H'(L')).
Collecting, with abs(e^{sL'}) = x'^{beta-1}, gives (B2_k). QED

**Zero sum.** Low zeros (abs(gamma) <= T0, beta = 1/2 by H1): by Lemma 6 and (B1'),
abs(F_phi(rho)/rho) <= x^{-1/2} h(x)(1 + kappa/abs(rho))/abs(rho)^2 + x'^{-1/2} H(L') m_c/abs(rho); summed with H4 and S1(T0) >= sum 1/abs(gamma).
High zeros (abs(gamma) > T0, 1 - b0 <= beta <= b0 by H2): abs(rho), abs(s) >= abs(gamma); each bound has the form
x'^{beta-1} g(abs(gamma)), so the pairing rho <-> 1 - conj(rho) of Lemma 7 (ii) gives
abs(sum_high) <= (x'^{b0-1} + x'^{-b0})/2 * [ H(L')(2k/lam)^k S_{k+2}(T0) + abs(H'(L'))(lam/2 S_2(T0) + 2 S_3(T0)) ],
S_m(T) >= sum_{abs(gamma)>T} abs(gamma)^{-m} from H3 exactly as in RobinHeight (certify2.Sm; valid for every m >= 2).
For b0 = 1 (no half-plane) the same holds with x'^{0} + x'^{-1}, beta < 1 being classical.

## 3. Lemma 13 (sharper bound for D)

If abs(theta(x)/x - 1) <= eta <= 1/2 and L = log x >= 20, then with D(x) as in Lemma 2,
D(x) <= (eta^2/(2L)) [ 1 + 2 eta/(3(1 - eta)) + (1 + eta)^2/L + 2 eta (1 + eta)^3/(3 L^2 (1 - eta(1 + eta)/L)) ].
Proof. u := theta(x)/x - 1, v := log(1+u)/L, D = u/L - log(1 + v). For abs(z) <= 1/2 the series give
log(1 + z) >= z - z^2/2 - abs(z)^3/(3(1 - abs(z))) and abs(log(1 + z)) <= abs(z) + z^2. Hence abs(v) <= eta(1 + eta)/L <= 1/2,
u/L - v = (u - log(1+u))/L <= (u^2/2 + abs(u)^3/(3(1 - abs(u))))/L, and -log(1 + v) <= -v + v^2/2 + abs(v)^3/(3(1 - abs(v))).
Add and insert the bounds for abs(u), abs(v). QED  (RobinHeight's bound eta^2/L (1 + (1+eta)^2/L) is asymptotically twice as
large as eta -> 0; e.g. 1.25 times at eta = 1/2, L = 20.)

## 4. Theorems (certified by check_all2.py)

Level 2-CH (brief): Theorems A2/B2 of RobinHeight with eta, eta_1 from Lemma 12; everything else identical (linear one-sided
ramp, Lemma 11's (B2), RobinHeight's D bound).
Level 3: Lemma 11_k with c = lam/2 (centred), (k, lam) chosen per box, D by Lemma 13. The Nicolas margin is
n(x)[M_low(x' e^lam) - abs(Z_phi) - abs(W) - x'^{-3}/L' - 1/(2(x-1))] > 0; the Robin margin is Lemma 4 with
-I(x) <= M_up(x') + abs(Z_phi) + log(2pi)/(x' L') + abs(W) and D from Lemma 13, theta(x) <= x(1 + eta) with eta at x' (eta decreasing).
Ranges: see REPORT.md.

## 5. What limits X now

Level 2-CH: window sqrt(x) eta lam/2 against zeros x^{b0-1/2} S3 (2/lam): X ~ (c/(S3 eta))^{1/b0} with eta = pi/3e12.
Level 3: the zeros above T0 are damped by (2k/(lam T0))^k, which needs lam T0 of order 2k (k ~ 10-30); the window then costs
about sqrt(x) (pi/T0) lam/sqrt(12k), i.e. sqrt(x) times a few 1e-24, and Robin's D term costs sqrt(x) eta^2/2 = 5.5e-25 sqrt(x).
The half-plane matters little at this level: the limit is pointwise knowledge of theta(x) at relative precision pi/T0.

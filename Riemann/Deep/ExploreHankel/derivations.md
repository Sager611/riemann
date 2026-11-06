# explore13 derivations (Wiener-Hopf / Hankel / Nehari-AAK / H-infinity -> RH tree)

Notation. s = 1/2 + i tau on the critical line; P = Cauchy(0, 1/2) probability measure, dP = dtau / (2 pi |s|^2).
H^2_+ / H^2_- = Hardy spaces of Re s > 1/2 / Re s < 1/2 = Mellin images of L^2(0,1) / L^2(1,inf).
Mellin: M f(s) = int_0^inf f(t) t^{s-1} dt, unitary up to 1/sqrt(2 pi) on the critical line.

## 1. Baez-Duarte in Mellin form (N5)
M[ {1/(k t)} ](s) = -k^{-s} zeta(s)/s (0 < Re s < 1), M[chi_(0,1]](s) = 1/s, so
  d_N^2 = min_{C in D_N} int |1 - zeta(s) C(s)|^2 dP,   D_N = span{k^{-s}: k <= N}.
Gram: <rho_{1/h}, rho_{1/k}> = Vasyunin (coprime) and G(dh, dk) = G(h, k)/d; b_k = (log k + 1 - gamma)/k; ||chi|| = 1.
Rational dilations: <rho_{a/b}, rho_{c/d}> = a c Vas(bc, ad) (homogeneity G(lambda alpha, lambda beta) = lambda G(alpha, beta));
<chi, rho_theta> = theta (log(1/theta) + 1 - gamma).

## 2. The only Hankel operator of the dilation system is rank one
Multiplicative convolution (rho * g)(t) = int g(theta) {theta/t} dtheta/theta has Mellin multiplier -zeta(s)/s.
Hankel part L^2(0,1) -> L^2(1,inf): for theta < 1 < t, {theta/t} = theta/t, so (rho*g)(t) = t^{-1} int_0^1 g(theta) dtheta:
rank one. Kronecker: zeta(s)/s = 1/(s-1) + (H^2_+ function); the single pole at s = 1 lies across the line.
(This is why Nyman's functions rho_theta - theta rho_1 vanish on (1, inf).)  The reverse Hankel L^2(1,inf) -> L^2(0,1),
kernel {e^{u+v}} e^{-(u+v)/2} in log variables, depends LINEARLY on the symbol zeta(s)/s, which is analytic in Re s > 1/2:
zeros of zeta are not singularities of the symbol, so its singular values carry no zero information (AAK: s_n = distance in
L^inf(line) to functions meromorphic in Re s < 1/2 with <= n poles; the obstruction is the growth of chi(s) zeta(1-s)/s to the left,
an FE/Gamma-factor property shared by Davenport-Heilbronn).

## 3. Beurling floor = norm of a Hankel operator with inner symbol
Phi(s) := zeta(s)(s-1)/s^2 in H^2_+, Phi = I O (inner-outer; I = Blaschke product of zeros with Re rho > 1/2, no singular part).
zeta(s)/s = Phi(s) s/(s-1), and s/(s-1) = conj(b_1) on the line, b_1 = (s-1)/s (unimodular there).
Continuous dilations theta in [lambda, 1] -> lambda -> 0: closure of {Phi A} = I H^2_+ (Beurling-Lax for the semigroup e^{-x s}).
  D_inf^2 = lim D(lambda)^2 = dist( (s-1)/s^2 , I H^2_+ )^2 = || P_{K_I} ((s-1)/s^2) ||^2 = || H_{conj I} ((s-1)/s^2) ||^2,
with H_{conj I} f = P_-(conj(I) f).  H_{conj I}^* H_{conj I} = P_{K_I} (model space): a partial isometry, every nonzero
singular value = 1, rank = number of off-line zeros (Kronecker).  So AAK's singular values give only the COUNT of off-line
zeros; RH <=> H_{conj I} = 0 is Beurling's criterion restated.
Residue formula (finitely many zeros rho_j, B = prod (s - rho_j)/(s - (1 - conj rho_j))):
  P_-(g/B) = sum_j c_j/(s - rho_j), c_j = Res_{rho_j} g/B,  ||.||^2 = sum_{j,l} c_j conj(c_l)/(rho_j + conj(rho_l) - 1)
(the kernel: (1/2pi) int dtau /((s - rho_j) conj(s - rho_l)) = 1/(rho_j + conj rho_l - 1)).
One zero: D_inf^2 = |rho - 1|^2 (2 sigma - 1)/|rho|^4 ~ (2 sigma - 1)/|rho|^2 (zeta form, target (s-1)/s^2);
for an entire f (Davenport-Heilbronn) the target is 1/s and D_inf^2 = (2 sigma - 1)/|rho|^2 exactly for one zero.
Monotone: adding zeros shrinks I H^2, so any finite subset of off-line zeros gives a lower bound for the floor.
Natural numbers (Baez-Duarte) are a subset of the continuous dilations: d_N >= D(log N); every relaxation to a richer
(linearly ordered) dictionary bounds d_N from BELOW only.  The RH direction (d_N -> 0) is an upper bound.
On-line zeros are boundary zeros: they belong to the OUTER factor O; their effect is the finite-memory rate
D(lambda)^2 log(1/lambda) >= sum m^2/|rho|^2 (BBLS, Burnol), a boundary-zero (Fisher-Hartwig-type) effect of a truncated
Wiener-Hopf problem, invisible to Hankel operators.

## 4. Sup-norm (AAK-native) version
e_N := min_{C in D_N} sup_tau |1 - zeta C| W, W = 1/|s|.  At an on-line zero rho: |1 - zeta(rho) C(rho)| W = 1/|rho|.
Hence e_N >= 1/|rho_1| = 0.0707101 for all N.  An off-line zero rho_0 adds an interior Pick constraint S(rho_0) = 1 for
S = 1 - zeta C, i.e. e >= |W(rho_0)| (W outer) = 1/|rho_0| < 1/|rho_1|: never binding.  General outer W:
|W(rho_0)| <= Poisson average of |W| over a window of width ~ (Re rho_0 - 1/2) at height Im rho_0, and on-line zeros there
are 2 pi/log T apart, so for |W| with modulus of continuity omega the off-line constraint is dominated up to omega(gap) plus
the Poisson tail.  Only weights oscillating on the zero-gap scale (built from the known zeros) could see off-line zeros:
tautological.  Same for Davenport-Heilbronn: e_N >= 1/|rho_1^DH| = 0.1953644 (rho_1^DH = 1/2 + 5.094160 i).

## 5. Davenport-Heilbronn analogue of Baez-Duarte
f(s) = sum a_n n^{-s}, a = (1, kappa, -kappa, -1, 0) mod 5, A(x) = sum_{n<=x} a_n = S(floor(x) mod 5), S = (0, 1, 1+kappa, 1, 0).
M[A(1/(k t))](s) = k^{-s} f(s)/s on Re s > 0 (A bounded; A(1/(kt)) = 0 for t > 1/k, so everything lives in L^2(0,1) = H^2_+).
Gram: int_0^inf A(x/h) A(x/k) dx/x^2, integrand periodic with period P = 5 lcm(h,k); a piece [x_i, x_{i+1}) with value c_i gives
  sum_{q>=0} c_i (1/(x_i + qP) - 1/(x_{i+1} + qP)) = (c_i/P)(psi(x_{i+1}/P) - psi(x_i/P))   (exact, digamma).
b_k = beta/k, beta = sum_{r=1}^4 S(r)(psi((r+1)/5) - psi(r/5))/5.
Constant: C_f = sum_rho 1/(rho(1-rho)) = 2 xi_f'/xi_f(1), xi_f(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s):
  C_f = log(5/pi) - gamma + 2 f'(1)/f(1) = 0.1952726 (zeta: 2 + gamma - log 4 pi = 0.0461914).

## 6. Wiener-Hopf index theory on vertical lines (N7, N13)
Douglas: T_phi Fredholm => phi invertible in L^inf.  sigma > 1: zeta(sigma + it) = exp(sum Lambda(n)/(n^s log n)) is invertible
in the Wiener almost-periodic algebra with mean motion 0: index 0 (no zeros; nothing to count).  1/2 < sigma <= 1:
zeta(sigma + it) is unbounded above and below (Bohr-Landau; Voronin universality for sigma < 1; |zeta(1+it)| << 1/loglog t i.o.),
so no Toeplitz/Wiener-Hopf Fredholm theory and no half-plane index.  sigma = 1/2: zeros on the contour.  Only bounded contours
(argument principle: Backlund/Turing) count zeros.  DH: identical (and its zeros with Re > 1 are counted by the argument principle).

## 7. Weil positivity (N1) as a truncated Wiener-Hopf operator
Q(f) = W(f * f~) = sum_gamma |f^(gamma)|^2 under RH: the convolution operator on L^2(0, L) whose symbol is the pure point measure
sum_gamma delta_gamma.  No spectral (Wiener-Hopf) factorization exists for a singular symbol (Szegő-Kolmogorov-Krein need
int log w' /(1 + tau^2) > -inf); window positivity is Krein's extension problem (explore8) or a frame bound (R368).
Several L-functions: Weil forms decouple over characters, the matrix symbol is diagonal, partial indices are the scalar ones.

## 8. Control dictionary (H-infinity / H2)
d_N = H2 tracking cost of the plant zeta (critical line = j omega axis) for the reference 1/s with a finite-impulse-response
controller of memory log N (delays log k).  RH <=> the 'zeta plant' is minimum phase on Re s > 1/2.  Off-line zeros = right
half-plane zeros: the cheap-control floor is the Blaschke-residue formula of section 3 (Qiu-Davison 1993 type).  On-line zeros =
j omega-axis zeros: they cost ~ sum m^2/(|rho|^2 log N) at memory log N (Burnol).  Nevanlinna-Pick (Pick functions -Xi'/Xi)
= Lagarias' Re xi'/xi > 0 criterion (R20).  Kronecker for zeta'/zeta on Re s = sigma_0: rank of the right-side Hankel =
1 + #{rho : Re rho > sigma_0}: Landau/Littlewood restated (N7, N10).

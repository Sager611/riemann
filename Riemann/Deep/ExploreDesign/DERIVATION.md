# ExploreDesign DERIVATION: guard-band exchange for the lattice design of Weil's form (explore27, 11 Nov 2026)

Notation as in ../WeilZeta/DERIVATION.md: c > 1, L = ln c, I = [-L/2, L/2], F(t) = int_I f e^{itx}, g = f * conj f(-.), supp g in [-L, L],
g^(t) = int g e^{itx} = |F(t)|^2, ||f||^2 = (1/2pi) int |F|^2; pole(f) = 2|<f,cosh(x/2)>|^2 - 2|<f,sinh(x/2)>|^2;
w_n = Lambda(n) n^{-1/2}; b_n = 2n + 1/2; Q(f) - pole(f) = <K, g>, K = K_arch + K_prime, K_prime = -sum_{n<c} w_n (delta_{log n} + delta_{-log n}),
K_arch = C delta_0 - sum_n e^{-b_n |x|} (regularised at 0; its Fourier transform is phi(t) = -log pi + Re psi(1/4 + it/2)).

## 1. The symbol at an arbitrary frequency
mu(t) := int_{-L}^{L} K(x) e^{-itx} dx = phi(t) + 2 sum_{n>=0} c^{-b_n} [b_n cos(tL) - t sin(tL)]/(b_n^2 + t^2) - 2 sum_{n<c} w_n cos(t log n).
Proof: int_{-L}^{L} = int_R - int_{|x|>L}; int_L^inf e^{-bx} cos(tx) dx = e^{-bL}[b cos tL - t sin tL]/(b^2 + t^2); e^{-b_n L} = c^{-b_n}.
At t = pi k/L (sin = 0, cos = (-1)^k) this is Chekmachev's mu(k) (WeilZeta DERIVATION sec. 3). wsym.py evaluates it in float64;
checked against mpmath at k = 1, 3, 101, 3689, 3693, 10263762037 (c = 47): agreement 1e-15 relative.

## 2. Theorem G (lattice identity with a guard band)
Let P >= 2L, t_k = 2 pi k/P, and kappa a finite real even measure (or L^1 function) supported in the guard band G_P = {L <= |x| <= P/2}.
Put w_k = mu(t_k) + kappa^(t_k), kappa^(t) = int kappa e^{-itx}. Then for every f in L^2(I), in (-inf, +inf]:
  (LATG)  Q(f) = pole(f) + (1/P) sum_{k in Z} w_k |F(t_k)|^2.
Proof. For f in C_c^inf(int I), g is smooth with supp g in (-L, L). Let Kt = K 1_{(-L,L)} + kappa on [-P/2, P/2], extended P-periodically.
Since g = 0 on G_P, <Kt, g> = <K, g> = Q(f) - pole(f). On R/PZ, Kt = (1/P) sum_k w_k e^{2 pi i k x/P} in D'(R/PZ) with
w_k = int_{-P/2}^{P/2} Kt e^{-2 pi i k x/P} dx = mu(t_k) + kappa^(t_k); pairing with the smooth periodisation of g gives
(1/P) sum_k w_k int g e^{i t_k x} dx = (1/P) sum_k w_k |F(t_k)|^2. The extension to L^2(I) is Chekmachev's Thm 3.1 argument verbatim:
w_k >= -C uniformly (kappa^ bounded, mu bounded below), sum_k |F(t_k)|^2 = P ||f||^2 (I fits in one period), monotone/Fatou on the
positive part. QED.  For P = 2L and kappa = 0 this is (LAT).

## 3. Lemma R (rigidity at P = 2L; why the vertex-exchange step needs P > 2L)
A weight perturbation eta_k keeps (LATG) exact for all f in L^2(I) iff sum_k eta_k |F(t_k)|^2 = 0 for all f, i.e. the P-periodic
distribution with Fourier coefficients eta_k vanishes on (-L, L). For P = 2L the complement of (-L, L) in a period is the single point
{+-L}; distributions supported there that annihilate every g with g(+-L) = 0 (all that L^2 f give) are a delta_{+-L}, so eta_k = 2a(-1)^k.
So Chekmachev's lattice weights are canonical up to the alternating sequence (the even/odd threshold split t1 = 1 / 0 is exactly this
freedom). Every forced (negative) node of the P = 2L lattice must enter the matrix. For P > 2L the exchange directions are
eta_k = kappa^(t_k) with kappa in L^1(G_P): band-pass sequences, normalised frequency in [2 pi L/P, pi] (alternating patterns).
In design language: (LATG) is a signed design on the lattice with the moment constraint 'Fourier transform = Weil's distribution on
(-L, L)'; kappa^ is the Boehning vertex-exchange direction, which exists only when the lattice oversamples (P > 2L).

## 4. Theorem C (certificate on P = mL; m = 2 is Chekmachev's Thm 5.10)
For integer m >= 2 the sublattices Lambda_r = {t_k : k = r mod m} = (2 pi/L)(Z + r/m) are orthogonal bases of L^2(I):
sum_{k = r mod m} |F(t_k)|^2 = L ||f||^2. Fix tau_r >= 0 (r = 1..m-1), deff = (1/m) sum_r tau_r. Let K (indices not in mZ) and
E (indices j of t_{mj} = 2 pi j/L) be finite with w_k >= tau_{k mod m} for k not in K, k not in mZ, and w_{mj} >= 0 for j not in E.
Then for eps < deff:  Q(f) - eps ||f||^2 >= sum_j dl_j |<f, phi_j>|^2 + (1/P) sum_{k in K} (w_k - tau_{k mod m}) |F(t_k)|^2 + pole(f),
dl_j = w_{mj}/m + deff - eps (>= deff - eps for j not in E), phi_j = L^{-1/2} e^{2 pi i j x/L}.
Proof: split (LATG) by sublattices; for r >= 1 bound w_k by tau_r off K and use sum_{Lambda_r} |F|^2 = L||f||^2; on Lambda_0,
|F(2 pi j/L)|^2 = L |<f, phi_j>|^2 and ||f||^2 = sum_j |<f,phi_j>|^2. QED.
The right side is P0 + (rank ones) - (rank ones) +- pole as in Chekmachev's sec. 5, with P0 diagonal in phi_j and closed-form Gram entries
<phi_j, e^{i t_k x}>_{L^2(I)} = L^{-1/2} 2 sin(pi(j - k/m))/((2 pi/L)(j - k/m)) (Cauchy structure), so his Woodbury bound (Lemma 5.3),
Cholesky residual test (Lemma 5.5) and the scan carry over unchanged. m = 2, tau_1 = 1 gives deff = 1/2 and his rule exactly.

## 5. Patches (the exchange step) and why they cost nothing at the near-null level
kappa = sum_d sum_i hat_i(|x|) [a_{d,i} cos(t_d x) + b_{d,i} cos(t_d x + (pi/2) sgn x)], hat_i = hat of half-width D centred at
x_i = L + i (P/2 - L)/(M + 1) (i = 1..M), D = (P/2 - L)/(M + 1): real, even, supported in G_P. Exactly:
kappa^(t) = sum_{d,i} D sinc^2((t - t_d) D/2) [a cos((t - t_d) x_i) + b sin((t - t_d) x_i)] + (mirror at t + t_d: D sinc^2((t + t_d)D/2)
[a cos((t + t_d)x_i) - b sin((t + t_d)x_i)]), sinc u = sin u/u. All closed forms (Arb-ready).
Loss: Q = minorant + sum_{k not in K} (w_k - tau) |F(t_k)|^2 / P (>= 0, dropped). With every node up to the low-block height in K, the
patch tails there are part of the exact low-block weights (no loss); the only loss is the patched high nodes' surplus times the
near-null mode's content above the low block, which is the same quantity Chekmachev already drops (his certificate closes with K2 =
5c^2/3, so that content is below the margin). Patches therefore remove forced high nodes from the matrix without touching the
RH-sensitive low block.

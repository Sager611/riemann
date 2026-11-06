# explore22 derivations (electrochemistry -> RH tree)

Notation: F is xi, mv (xi with zeros 6, 7 replaced by 1/2 +- delta + i 39.2524, delta = 0.15 unless stated) or dh
(Davenport-Heilbronn Lambda); all satisfy F(s) = F(1-s) and are real on Re s = 1/2. z_rho = -i(rho - 1/2), w_rho = z_rho^2,
x_rho = 1/w_rho. Code: [code/](code/); outputs: [results/](results/).

## D1. The passive xi-network (T1)
- Hadamard: F'/F(1/2 + w) = sum_rho 1/(w - i z_rho). If every z_rho is real (RH for F): = sum_{gamma>0} 2w/(w^2 + gamma^2),
  Foster's first canonical form: parallel LC tanks in series, every capacitance 1/2, inductances 2/gamma_n^2. Lossless,
  positive real (PR) on Re w > 0. Conversely a zero with beta > 1/2 is a pole of F'/F in Re w > 0, impossible for PR.
  So "the xi-network is passive" <=> RH: Lagarias (Acta Arith. 89, 1999), <=> |xi(sigma+it)| increasing in sigma > 1/2
  (Sondow-Dumitrescu 2010) <=> E_h(z) = xi(1/2 + h - iz) Hermite-Biehler for every h > 0 (de Branges spaces H(E_h)).
- Shifted network: Re F'/F(s) = sum_rho (sigma - beta)/|s - rho|^2 > 0 whenever sigma > sup beta. Unconditionally the
  xi-network shifted by h > T_* - 1/2 = 0.37338 (Lean-checked half-plane) is passive (lossy); for DH only h > sup beta_DH - 1/2,
  and DH has zeros in sigma > 1 (Davenport-Heilbronn 1936), so h_DH > 1/2. The network sees exactly the zero-free half-plane (N7).
- Euler-product pieces: (1 - p^-s)^-1 = (1 + coth(s log p/2))/2, coth(sT) = input impedance of an open-circuited lossless line
  of delay T. Z_P(s) = sum_{p<=P} (log p/2) coth(s log p/2) is PR (Foster) on Re s > 0, and on Re s > 1
  -zeta'/zeta(s) = Z_P(s) - theta(P)/2 + (primes > P). Hence
  xi'/xi(s) = theta(P)/2 - Z_P(s) + (1/2)psi(s/2) - (1/2)log pi + 1/s + 1/(s-1) + (primes > P):
  the lossless Euler network enters with the WRONG sign (-Z_P, Re < 0 on Re s > 0), balanced by a resistor theta(P)/2 -> infinity
  (the prime number theorem main term), plus non-passive archimedean and pole elements (Re psi(1/4) < 0; Re 1/(s-1) < 0 for
  Re s < 1). Since Re Z_P(sigma+it) <= Z_P(sigma) = theta(P)/2 + sum_{p<=P} log p/(p^sigma - 1), the balance closes only for
  sigma > 1 (trivial region). Passivity of the parts does not transfer; the sum's passivity is RH (the sign problem of Weil's
  prime terms, R379; the network form of explore12's quantum graphs).
- Test (results/t1_pr.jsonl, t in [1, 200], step 0.1): G(eps, t) = Re F'/F(1/2 + eps + it).
  xi: G > 0 on eps = 0.02, 0.1, 0.2, 0.3 (min 2 eps sigma_1 at t = 1, sigma_1 = sum 1/gamma^2 = 0.0231050).
  mv: G < 0 only on t in [39.2, 39.4] for eps = 0.02, 0.1 (beta - 1/2 = 0.15); G > 0 for eps = 0.2, 0.3.
  dh: G < 0 only at the heights of the known off-line zeros with beta - 1/2 > eps: eps = 0.02: 85.5-85.9, 114.1-114.3, 166.5,
  176.5-176.9; eps = 0.1: 85.5-85.9, 114.1-114.2, 176.6-176.8; eps = 0.2: 85.5-85.9, 176.7; eps = 0.3: 85.65-85.75 (min -114.8).
  So the passivity test on a line = detection of zeros to its right, at their own heights: zero verification.

## D2. Cauer RC ladder of the zero network (T2)
- Z(p) = sum_rho 1/(w_rho + p) = sum_k (-1)^k sigma_{k+1} p^k, sigma_k = sum_rho w_rho^-k. If all w_rho > 0 this is an RC
  (Stieltjes) impedance: Foster RC form with time constants 1/gamma^2; Cauer's second form from the low-frequency expansion
  is the S-fraction c_1/(1 + c_2 p/(1 + c_3 p/(...))) with every c_k > 0 (element values are products of c_k's).
- c_k via Hankel LU pivots: p0_n, p1_n pivots of [m_{i+j}], [m_{i+j+1}], m_k = sigma_{k+1}; c_1 = m_0, c_{2n} = p1_n/p0_n,
  c_{2n+1} = p0_{n+1}/p1_n. All c_k > 0 <=> both Hankel sequences positive <=> Stieltjes moments <=> (Grommer 1914) all
  w_rho > 0 <=> RH for F. Same criterion: Polson 2026 (Thorin series, quoted in 2611.02944), Matiyasevich 0709.0028; R296, R307.
- Gram reading (P(x) = sum a_i x^i): sum a_i a_j m_{i+j} = sum_rho x_rho P(x_rho)^2 = sum_rho g(z_rho)^2, g(z) = z^-1 P(z^-2): Weil-type
  positivity sum_rho Fh(z_rho)^2 >= 0 on the span of z^{-2i-1} (rational, non-compactly supported test family); finite
  sections are the Gram-matrix certificate method with a low-height basis.
- Unconditional data (no zeros used): xi: Phi(u) = sum_m (2 pi^2 m^4 e^{9u} - 3 pi m^2 e^{5u}) exp(-pi m^2 e^{4u}), even,
  int_0^inf Phi cos(zu) du has zeros at 2 gamma; M_2n = int_0^inf Phi u^2n du by the trapezoid rule (strip |Im u| < pi/8, shift
  pi/16, h from the target digits; checked by hfac 0.7-0.8). dh: Lambda(1/2 + iz) = 2 int_0^inf Psi(v) cos(zv) dv,
  Psi(v) = 2 e^{3v/2} theta(e^{2v}), theta(x) = sum a_n n e^{-pi n^2 x/5}; Psi even (checked 1e-16 at v = 0.1, 0.3, 0.7;
  2 int Psi = Lambda(1/2) = 1.43319216877500844 and Lambda(1/2 + 3i), (1/2 + 7.5i) reproduced to 20 digits). Power sums by
  Newton's identities; mv = xi's sums with gamma_6, gamma_7 replaced by gbar -+ i delta. Check: sigma_1(xi) =
  0.0231049931154189707889338 (= sum 1/gamma^2), sigma_2 = 3.7172599285e-5.
- Results (results/t2_*.json): xi: c_1..c_220 > 0 (dps 650). dh: negative exactly c_149, c_151, c_218, c_220, all other
  c_1..c_232 > 0 (dps 650 and 850, h factors 0.8 and 0.7 agree to 3 decimals in log10 c_k; dps 450 drifts from c_185 on). DH sign changes on the line
  (results/t2_count.json): 43 below 85.6993, 61 below 114.1633, so the off-line pairs are atoms 44-45 and 64-65 of the
  zero measure; detection at Hankel order 75 and 109 (1.7x the rank). mv: first negative c_26, 28, 29, 32, 36, 39 for
  delta = 0.3, 0.15, 0.05, 0.01, 0.001, 0.0001 (signature '-+-'; about +3.5 per decade of delta, unbounded as delta -> 0, since
  at delta = 0 the pair is a real double zero and every c_k > 0).
- Cost: precision about 5 digits per Hankel order (dps 450 fails near order 92 for DH; dps 650 holds to order 110). An off-line pair of rank k needs
  order >= 1.7k: for heights near 3e12 (k ~ 1e13) about 1e14 digits; and positivity of c_1..c_K excludes no off-line zero at
  any height (delta small enough passes). Strictly weaker than zero verification.

## D3. Kramers-Kronig / Bode (T3)
- KK holds for every function analytic and bounded-type in the right half-plane; the gain-phase (Bode) relation holds iff the
  function is minimum phase (outer). For F(1/2 + a + w): minimum phase <=> no zeros with beta > 1/2 + a. a = 0: Beurling's outer
  criterion (N5, explore13); Bode's sensitivity integral is Jensen/Littlewood: Balazard-Saias-Yor 1999,
  (1/2 pi) int log|zeta(1/2+it)| dt/(1/4+t^2) = sum_{beta>1/2} log|rho/(1-rho)|.
- Each off-line pair adds 2 log|rho/(1-rho)| = (2 beta - 1)/|rho|^2 (results/t3_bode.json): DH 8.40e-5 (85.70), 2.31e-5
  (114.16), 5.37e-6 (166.48), 1.44e-5 (176.70), total 1.27e-4 (explore13's floor 2.54e-4 counts the 8 zeros of the
  quartets); mv 1.95e-4; a pair at height 3e12 with beta - 1/2 = 1e-3: 2.2e-28. A KK/minimum-phase test on line data must resolve
  these numbers: far worse than locating the zero.

## D4. Diffusion, Warburg, Biane-Pitman-Yor, Marcus-Hush (T4)
- Finite-length Warburg (blocking end): Z = coth(sqrt x)/sqrt x = 1/x + sum_n 2/(x + n^2 pi^2): DRT = atoms at tau_n =
  1/(n pi)^2; current response K(t) = 1 + 2 sum e^{-n^2 pi^2 t} = theta_3(0, e^{-pi^2 t}); the electrochemists' short-time
  (erfc/image) series and long-time (eigenmode) series are Jacobi's imaginary transformation; Mellin of K - 1 is
  2 pi^{-2w} Gamma(w) zeta(2w) (spectral zeta of the cell; Riemann's second proof). Transmissive cell: tanh sqrt x/sqrt x,
  rates (n - 1/2)^2 pi^2; membrane permeation (exit flux per inlet step, Devanathan-Stachurski / Daynes time-lag cell)
  sqrt x/sinh sqrt x = prod (1 + x/(n^2 pi^2))^-1; BPY (math/9912170): 2 xi is the Mellin transform of the scaled range of a
  Brownian bridge, whose square has Laplace transform (sqrt(2 lambda)/sinh sqrt(2 lambda))^2 (squared permeation transfer).
- Adversarial: DH has the same structure with a positive kernel: for v >= 0, x = e^{2v} >= 1, |a_n| <= 1 gives
  theta(x) >= e^{-pi x/5}(1 - 2e^{-3 pi/5} - 3e^{-8 pi/5} - ...) >= 0.67 e^{-pi x/5} > 0, and Psi is even, so Psi > 0 on R
  (numerically min on [0, 3] is 1.5e-108 > 0). Lambda_DH(1/2 + iz)/Lambda_DH(1/2) is the characteristic function of the
  probability density Psi/int Psi, yet has off-line zeros. Positive relaxation kernels are DH-blind (R10.2, R342).
  2D analogue: the spectral zeta of a periodic diffusion cell (flat torus R^2/(Z + 5^{-1/2} Z), dual lattice Z + sqrt5 Z) is
  (4 pi^2)^{-s} times the Epstein zeta of x^2 + 5y^2 (class number 2), which has infinitely many zeros in sigma > 1
  (Davenport-Heilbronn 1936): physical diffusion operators do not force their spectral zeta onto a line.
- Marcus-Hush-Chidsey: k(eta) is a Gaussian (variance ~ lambda k_B T) smoothing of the Fermi function, i.e. the forward heat flow,
  the de Bruijn-Newman direction t < 0, which by Rodgers-Tao (Lambda >= 0) creates non-real zeros of H_t; the complete Fermi-Dirac
  integral at zero chemical potential is Gamma(s) eta(s) = Gamma(s)(1 - 2^{1-s}) zeta(s): zeta values, no zero statement.
  Butler-Volmer i ~ e^{a f eta} - e^{-(1-a) f eta} has zeros only on Re eta = 0 (two-term exponential sums); multistep rate laws
  are exponential polynomials with Polya-Dickson strips, like Dirichlet polynomials: no handle on Xi.

## D5. Positive relaxation spectra on Nyman-Beurling (T5)
- DRT inversion imposes G >= 0. Transfer to N5: Baez-Duarte with c_k >= 0: d+_N^2 = min_{c>=0} ||chi - sum c_k {1/(kt)}||^2.
  On t > 1, {1/(kt)} = 1/(kt): error S/t, S = sum c_k/k, contributing S^2; on (1/2, 1): {1/t} = 1/t - 1 and {1/(kt)} = 1/(kt)
  (k >= 2): error 1 + c_1 - S/t. With A = 1 + c_1: d^2 >= 2S^2 - 2 ln2 A S + A^2/2 >= (A^2/2)(1 - ln^2 2), so c_1 >= 0 alone
  gives d^2 >= 0.259773 for every N; the natural coefficients are c_k = -mu(k) (sum_k mu(k)[x/k] = 1), so c_1 = -1, A ~ 0.
- NNLS (results/t5_posnb.txt; Vasyunin Gram as explore13): zeta d+^2 = 0.858, 0.506, 0.2947, 0.2733, 0.27105, 0.271034
  (N = 1, 2, 5, 10, 20, >= 40; support 3-9, 11, 13, 17, 19, 23, 31) against free d^2 = 0.00935 at N = 160. DH: d+^2 plateaus
  at 0.066867 from N = 20 (free 0.0347 at N = 80). Sign variants for zeta (N = 10, 40, 160): c_1 free, c_k >= 0 (k >= 2):
  0.028132, 0.027225, 0.027225 (stalls); Mobius pattern c_k (-mu(k)) >= 0: 0.025883, 0.013511, 0.009659; c_1 <= 0 only: equal to
  free 0.023772, 0.012715, 0.009346. Positivity stalls both, RH-true and RH-false alike; the working sign pattern is Mobius's, so
  the stall measures the sign structure of Weil's prime kernel, not the zeros.

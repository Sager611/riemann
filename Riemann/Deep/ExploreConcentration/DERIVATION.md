# explore44 DERIVATION (functional inequalities / concentration -> RH tree)

## 1. Choice (15:34-16:25 UTC, measured before choosing)
Candidates measured: N9 (Bohr lift, Herbst/LSI/Weissler-Bayart), N8 (linear statistics; Theorem U' Lemma M), N3/N5 (heat flow), N12 (LinComb3 7.4 non-concentration).
Chosen: N9 + N7 interface. Tool: tensorized log-Sobolev on T^inf (circle LSI constant 1, Ent(g^2) <= 2 E|g'|^2) + Herbst, with Talagrand T2 by
Lipschitz push-forward, and Weissler-Bayart hypercontractivity. Reason: it is the only one of the five tools that sees the Euler product
(through tensorization), so it is the only one that can separate zeta from Davenport-Heilbronn (DH); test whether that separation reaches zeros.

## 2. Statements
S1 (T^inf, proved). sigma > 1/2, Lambda_sigma(theta) := log zeta_lift(sigma, theta) = -sum_p log(1 - p^{-sigma} e^{i theta_p}).
 d/dtheta_p Lambda = i p^{-sigma} e^{i theta_p}/(1 - p^{-sigma} e^{i theta_p}), so |grad Re Lambda|, |grad Im Lambda| <= L(sigma),
 L(sigma)^2 := sum_p p^{-2 sigma} (1 - p^{-sigma})^{-2} (finite iff sigma > 1/2).
 LSI tensorizes (constant 1 on T^inf, dimension-free) => Herbst: E exp(lam (F - EF)) <= exp(lam^2 L^2/2) for F = log|zeta_lift| or arg zeta_lift,
 P(|F - EF| >= r) <= 2 exp(-r^2/(2 L^2)); E F = 0 (log|1 - a e^{i theta}| has mean 0 for |a| < 1).
 T2 (Otto-Villani + push-forward by an L-Lipschitz map): the Bohr-Jessen law mu_sigma of log zeta(sigma + it) on C satisfies
 W2(nu, mu_sigma)^2 <= 2 L(sigma)^2 H(nu | mu_sigma) for every nu.
S1' (DH, adversarial). F = c1 L(s,chi) + conj(c1) L(s,chi-bar), chi mod 5 complex, c1 = (1 - i k)/2, k = 0.2840 (a(n) = 1, k, -k, -1, 0).
 Its lift c1 e^{G1} + conj(c1) e^{G2} vanishes on T^inf for sigma in (1/2, sigma_max) (positive zero density in strips): log|F_lift| is not Lipschitz,
 P(log|F_lift| < -r) ~ c e^{-2r} (codimension-2 zero set), so no Gaussian lower tail: S1 separates zeta from DH.
S2 (the node quantity; Littlewood). For F = 1 + sum_{n>=2} a_n n^{-s} (zeta, DH), sigma > 1/2:
 (2 pi/T) sum_{beta > sigma, T < gamma <= 2T} (beta - sigma) = E_T log|F(sigma + it)| + O(log T/T),  E_T = mean over t in [T, 2T].
 On T^inf: E log|zeta_lift| = 0; E log|F_lift| = phi_F(sigma) (Jessen function), phi_DH > 0 where DH has zeros.
 So the zero sum is exactly the MEAN SHIFT E_T log|F| - E_inf log|F_lift| (zeta) of a non-Lipschitz function; S1 controls fluctuations
 about means, never this shift. That is the claimed barrier; tested below.
S3 (transfer cost, quantified). (a) Tail route: if P_T(log|zeta(sigma+it)| < -r) <= exp(-r^2/(2L^2)) for r <= r*, a zero rho with
 0 < beta - sigma = d gives log|zeta| <= log(2d) + log M_rho on |t - gamma| <= d (M_rho = max|zeta'| near rho), so
 #{rho : beta - sigma < e^{-r}/(2 M)} <= T e^{r - r^2/(2L^2)} M; integrating sigma over [sigma0, 1] at r = 2L^2:
 N(sigma0 + e^{-2L^2}, T) <= T M e^{2L^2}/2 = O_sigma(T): Bohr-Landau 1914 strength, no saving.
 (b) Mean route = Jensen: E_T log|zeta M| <= (1/2) log E_T|zeta M|^2 loses exactly Var(log|zeta M|) (lognormal), and replacing Jensen
 by Herbst needs E_T of |zeta M|^{2k}, i.e. mean values of zeta^k M^k: the Ingham/Montgomery density inputs; no gain from S1.
 (c) Weissler-Bayart: ||F(. + s')||_{H^q} <= ||F||_{H^2} iff 2^{-s'} <= sqrt(2/q): s'(q) = (1/2) log2(q/2) >= 1/2 for q >= 4, i.e. the shift
 leaves the critical strip for every moment beyond the 4th (prime 2 is the bottleneck); removing p <= P gives s' = log(q/2)/(2 log P)
 at cost prod_{p<=P}(1 - p^{-sigma})^{-1}; on [T,2T] the inequality for zeta^k is Titchmarsh 13.2's mean-value theorem (Lindelof in mean).

## 3. Tests so far (code and logs in .scratch/explore44/, peaks <= 100 MB, single thread)
Correction to S1: for F = log|zeta_lift| the sharp per-prime constant is a/(1 - a^2) (a = p^{-sigma}): L_re^2 = sum p^{-2s}/(1 - p^{-2s})^2;
 a/(1 - a) is the constant for arg. c1 (smallest case, one prime p = 2): max|f'| = a/(1 - a^2) to 1e-10; E f = 0; E e^{lam f} = 2F1(lam/2, lam/2; 1; a^2)
 (hyp2f1 = quadrature to 1e-6); Herbst bound holds and is loose (sigma = 0.6, lam = 2: 1.77 vs 15.3); LSI Ent(g^2) <= 2E|g'|^2 holds (0.17 <= 0.32 at lam 1).
 L_re(sigma)^2 (primes to 1e7 + tail): 3.735 (0.55), 2.776 (0.6), 1.778 (0.7), 1.238 (0.8), 0.681 (1.0); Var log|zeta_lift| = 0.797 (0.6),
 0.526 (0.7), 0.386 (0.8): Herbst's variance proxy is 3.2-3.5 x the true variance (factor 2 from sin^2, rest from p = 2, 3).
c2 (lift Monte Carlo, primes <= 3e4 + Gaussian tail; 2e4 zeta, 4e4 DH samples): E log|zeta_lift| = -0.009 +- 0.006, 0.001 +- 0.005, -0.005 +- 0.004
 (sigma 0.6, 0.7, 0.8; exact value 0); DH Jessen function phi(sigma) = E log|F_lift| = 0.0442 +- 0.0035, 0.0091 +- 0.0028, 0.0030 +- 0.0022.
 Lower tails: zeta P(log|.| < -2) = 0.0072 (0.6), 0 at -3 (super-Gaussian); DH P(< -3) = 1.28e-3, P(< -4) = 1.7e-4 (ratio 7.5 = e^{2.0}: e^{-2r}).
c3 ([1000, 2000], step 0.01, Euler-Maclaurin, checked against mpmath to 4e-13): zeta mean log|zeta(sigma+it)| = -0.0008 (all three sigma; RH => O(log T/T)),
 var 0.782 / 0.528 / 0.389 (lift 0.796 / 0.534 / 0.390), P(< -1) 0.137 (lift 0.135), P(< -2) 0.0064 (lift 0.0072); DH mean 0.0291 / 0.0087 / 0.0006
 (lift 0.044 / 0.009 / 0.003), P(< -3) 1.6e-3 (lift 1.3e-3), P(< -4) 2.4e-4 (lift 1.7e-4). Bohr-Jessen already accurate at T = 1000 for both.
 Herbst bound P(log|.| < -r) <= exp(-r^2/(2 L^2)), L^2 = 2.776 at 0.6: r = 3: 0.198 vs zeta 5e-5 and DH 1.6e-3: BOTH satisfy it at every observed r;
 a Gaussian bound with proxy L^2 is beaten by e^{-2r} only beyond r = 4 L^2 = 11.1 (probability e^{-22}): the Euler-product information
 that LSI/Herbst carries is invisible at accessible levels; what separates zeta from DH in the data is the lift's MEAN and exact tail shape.
c4 (Littlewood on [sigma, 3] x [10, 420], DH, the 9 known off-line zeros of explore23): LHS 8.07649 vs 2 pi sum(beta - 0.6) = 8.07572;
 3.80769 vs 3.80731 at sigma = 0.7. So the mean shift E log|F(sigma+it)| (0.0211 over [10, 420]) is the zero sum (0.0197) plus horizontal
 terms O(1/T): S2 confirmed numerically; for zeta the same mean is 0 to 1e-3.

## 4. Transport form of Littlewood's lemma (the one place the tools enter an RH-relevant inequality)
Prop TL. sigma > 1/2, mu_T = law of zeta(sigma + it), t uniform on [T, 2T]; mu_inf = law of zeta_lift(sigma) (Bohr-Jessen); R > 0. Then
 (2 pi/T) sum_{beta > sigma, T < gamma <= 2T} (beta - sigma) <= e^R W1(mu_T, mu_inf) + A(R) + O(log T/T),
 A(R) := E log+(e^{-R}/|zeta_lift|) <= int_R^inf exp(-r^2/(2 L_re(sigma)^2)) dr = sqrt(pi L^2/2) erfc(R/sqrt(2 L^2)).
Proof. g_R(z) = log max(|z|, e^{-R}) is e^R-Lipschitz on C and log|z| <= g_R. Littlewood: zero sum = E_T log|zeta| + O(log T/T)
 <= E_T g_R = E_inf g_R + (E_T - E_inf) g_R <= [E_inf log|zeta_lift| + A(R)] + e^R W1 (Kantorovich-Rubinstein); E_inf log|zeta_lift| = 0;
 A(R) = int_R^inf P(log|zeta_lift| < -r) dr, bounded by Herbst (S1). QED. Optimising R: zero sum <= W1^{1 - o(1)} as W1 -> 0.
Where the Euler product enters: only through A(R) (Herbst needs the lift's Lipschitz bound = tensorization over primes).
 For DH the same argument gives zero sum <= phi_DH(sigma) + e^R W1 + A_DH(R) with A_DH(R) ~ c e^{-2R} (lift vanishes): no information.
Converse (so TL is a re-encoding, not a gain): g(z) = (e0 - |z|)+ is 1-Lipschitz; a zero at distance d < e0/(4M) right of the line
 puts mu_T-mass >= 2d/T into B(0, e0/2) while mu_inf(B(0, e0)) = o(e0^A) (super-Gaussian lower tail): W1 >= c (2 pi/T) sum over zeros
 within e0/(4M) of the line. So W1-rates of Bohr-Jessen convergence and zero density near sigma are the same problem up to w^{o(1)}.
Known W1-type rates: Lamzouri-Lester-Radziwill (1402.6682) discrepancy << (log T)^{-sigma} for log zeta, proved using zero-density
 estimates; inserted into TL: zero sum << (log T)^{-sigma + o(1)}, a log-power saving; Carlson (1920) already gives T^{4 sigma(1-sigma)-1+eps}.
Numbers (c6, T = 1000): A(R) zeta: lift 1.6e-3 (R = 2), 0 (R >= 3) = [T,2T] values; Herbst bound 0.48, 0.15, 0.034 (R = 2, 3, 4, sigma 0.6).
 DH: lift 4.2e-3, 6.0e-4, 5.7e-5 vs [T,2T] 4.8e-3, 8.0e-4, 1.3e-4 (e^{-2R}). 1-D W1 of the laws of log|.| ([T,2T] vs lift) = 0.027 (zeta, 0.6),
 0.038 (DH): the Bohr-Jessen error at T = 1000 equals DH's entire zero sum (0.029): at accessible T no W1 argument separates them.

## 5. N8 side (brief's question on Lemma M of Theorem U')
DPC's R_T = sum_{n < T^theta} c_n n^{-s}, c_n = mu(n) n^{-1/2} f(log n/L) prod rho_eta: a chaos of degree <= k on T^inf. Weissler (analytic):
 ||P_r g||_q <= ||g||_2 for r <= sqrt(2/q) => E_inf|R_T|^q <= (sum_j (q/2)^{j/2} ||R_j||_2)^q <= (k+1)^{q/2} (q/2)^{kq/2} ||R_T||_2^q, EVERY q, EVERY theta.
On [T,2T]: q = 2m even and m theta < 1 => E_T|R_T|^{2m} = E_inf|R_T|^{2m}(1 + O(T^{m theta - 1})) (MV); for theta > 1/2 only q = 2.
 General polynomials of length T^theta, theta > 1/2, break it (c5: resonant D, 4th moment [T,2T]/Haar = 0.99, 0.99, 1.12, 1.64, 3.35, 12.0 at
 theta = 0.39, 0.50, 0.61, 0.70, 0.79, 0.91, T = 2000); Lemma Z + zeros give q < kappa/theta (Lemma M); RH gives all q.
So the T^inf inequality holds for every theta and q; the whole gap 7/11 -> 1/2 is the transfer, which only zero information (kappa,
 6e' hairpins: 2 theta < kappa independently of Lemma M) can pay for. No functional inequality on T^inf moves Lemma M.

## 6. More tests (16:08-16:16 UTC)
c3 at T = 4000 (sigma 0.6, step 0.02, 2e5 points, peak 289 MB: chunks must shrink for larger T): zeta mean -0.00025, var 0.790 (lift 0.796),
 P(< -2) 0.0067; DH mean 0.0335 (T = 1000: 0.0291; lift phi = 0.0442), P(< -3) 6.7e-4, P(< -4) 1.0e-4. DH's zero sum per unit height grows
 toward the Bohr-Jessen value as T grows: finite-T value distributions see the primes up to ~T only (the Kronecker transfer cost).
c7 rotation invariance (the precise reason T^inf inequalities cannot bound [T,2T] moments past length T^{1/2}): D_{t0}(t) = sum_{n<=N}
 n^{-1/2-i(t-t0)} has the same T^inf law for every t0 (z_p -> p^{i t0} z_p preserves Haar), hence the same LSI/Herbst/Weissler bounds; at
 T = 2000 the [T,2T] 4th moments are, as ratios to the Haar value 43.5 / 84.4 / 147.8 / 226.9 / 334.6 / 530.5 (N = 20 ... 1000):
 t0 = 0: 1.01, 1.02, 1.00, 0.95, 0.84, 0.67;  t0 = 1.5T: 0.99, 0.99, 1.12, 1.64, 3.35, 12.0. At theta = 0.91 the two differ by 17.9x.
 So any moment bound that depends only on the T^inf law of the coefficients is >= the resonant value ~ N^2/(T log N), useless for
 theta > 1/2; Lemma Z (explicit formula, not rotation-invariant) is what bounds R_T, and only zero density pays for it.
c8 Prop TL at scale (T = 1000, n = 2500-point exact assignment): W1(mu_T, mu_inf) on C = 0.198 (zeta, 0.6) with same-law floors 0.179
 (lift-lift) and 0.112 (T-T); 0.112 at 0.7 (floors 0.111, 0.088); DH 0.150 / 0.069. TL bound min_R e^R W1 + A_Herbst(R) = 1.64 (0.6), 0.94 (0.7)
 for zeta, whose true zero sum is 0; Jensen (Bohr-Landau) gives 1/2 log zeta(2 sigma) = 0.861 (0.6), 0.567 (0.7) [corrected from 0.458, mpmath]; DH's actual sum 0.029.
 TL is weaker than Jensen at accessible T, and its input (a W1 rate) is itself proved through zero-density estimates (LLR).
S3(a) corrected constant: zeros with 0 < beta - sigma <= e^{-r}/(4M) give |zeta| < e^{-r} on |t - gamma| <= e^{-r}/(4M), so
 #{...} <= 2 M T e^{r} G(r); integrating sigma over [sigma0, 1]: N(sigma0 + e^{-r}/(4M), T) <= 8 M^2 (1 - sigma0) T e^{2r} G(r); r = 2L^2:
 <= 8 M^2 (1 - sigma0) e^{2 L^2} T = 825 M^2 T at sigma0 = 0.6 (e^{5.55} = 257.8); beats N(T) ~ (T/2 pi) log T only if log T > 5180 M^2.

## 7. Reviewer (verdict-only, items inline, 16:18 UTC) and repairs
Verdicts: S1 Lipschitz/Herbst CORRECT; Prop TL CORRECT; converse GAP; tail route FLAW; Weissler-Bayart CORRECT (sqrt(k+1) removable);
 rotation invariance CORRECT but only excludes constant-times-Haar comparisons.
Repairs. (i) Converse: with M a common bound for |zeta'| on the discs, intervals clipped to [T, 2T] and m_loc the maximal overlap (with
 multiplicity, <= local zero count): W1 >= (e0^2/(8M)) N_near/(T m_loc) - e0 mu_inf(B(0, e0)), N_near = #{rho : |beta - sigma| <= e0/(8M)}.
(ii) Tail route (was wrong: r = 2L^2 maximises e^{2r - r^2/2L^2}). Correct reading: a Gaussian lower tail on [T, 2T] up to level r*, for all
 sigma >= sigma0, gives N(sigma0 + e^{-r*}/(4M), T) <= 8 M^2 (1 - sigma0) T exp(2 r* - r*^2/(2L^2)), strong once r* > 4L^2, and zero-freeness
 if r* = inf; conversely each zero within e^{-r}/(4M) of the line adds e^{-r}/(2M T) to the tail. So the [T,2T] tail at level r* IS a
 density statement near sigma0; the T^inf tail does not transfer for free: Selberg's approximation log zeta = P_X + small holds off an
 exceptional set whose measure is the density input, and k matched moments of P_X (X^k <= T) give only polynomial tails (CkL^2/r^2)^k.
 Either way the output is the input up to e^{2r}: circular, no gain. (M <= T^{c(sigma0)} by convexity.)
(iii) Weissler constant: R = P_r Q, Q = sum_j r^{-j} R_j, r = sqrt(2/q): ||R||_q <= ||Q||_2 <= (q/2)^{k/2} ||R||_2.
(iv) Rotation invariance excludes law-only bounds of size C x Haar for q >= 2/theta (resonance X^{q/2}/(log X)^C >> T, also for DPC's R_T:
 rotating by z_p -> -p^{i t0} z_p turns c_n into |c_n| n^{i t0}, peak sum|c_n| ~ X^{1/2}/log X); for 2 < q < 2/theta a law-only bound is
 not excluded: Montgomery's large-values conjecture (R << N^{2+eps} V^{-2}) would give it up to T^eps (dyadic V: MV below N^{1/2+eps}, LVC above; total T N^{q/2} + N^q, fine iff q <= 2/theta), and no T^inf inequality proves that (fractional [T,2T] moments
 are not Haar averages). Decisive for N8 anyway: Theorem U''s 6e' (hairpins, e^{2 theta a} against e^{-kappa a}) needs 2 theta < kappa by
 itself, so even Lemma M for every q leaves 7/11.

## 8. The rate TL needs, measured (c9, sigma = 0.6, 2e4 random t per T, lift 6e4 samples, primes <= 2e4 + Gaussian tail)
1-D W1(law of log|F(sigma+it)| on [T,2T], lift law), same-size lift-vs-lift floor 0.011 (zeta), 0.011 (DH):
 T = 250, 500, 1000, 2000, 4000: zeta 0.063, 0.038, 0.028, 0.030, 0.022; DH 0.067, 0.055, 0.035, 0.037, 0.030.
 Sample means: zeta 0.010, 0.004, 0.003, 0.010, -0.006 (+- 0.006; true zero sum 0); DH 0.025, 0.024, 0.037, 0.030, 0.036 (lift 0.045).
 Over T <= 4000 the W1 error (0.02-0.06) is as large as DH's entire zero sum (0.025-0.037): no W1 argument separates zeta from DH there.

## 9. Second reviewer (framing, 16:27 UTC): INCOMPLETE -> scoped
A: tensorization also keeps prime independence/chaos structure (pseudo-moments), not only small-ball: headline narrowed to "T^inf information".
B: TL and its converse are one-sided, not an equivalence: rephrased ("a W1 rate is at least as strong as a near-line density bound").
C: multiscale prime-block routes (Soundararajan/Harper under RH, HRS, Arguin-Bailey) not examined: listed as out of scope in REPORT.
D: headline scoped to the five transfers tested; LVC phrased as "would supply", not "is needed".

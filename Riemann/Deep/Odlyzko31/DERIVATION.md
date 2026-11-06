# odlyzko31: derivation

All statements below are proved here or cited; every numerical constant is computed in the scripts (names in brackets) and can be found there with rg.

## 0. Notation and inputs
- K a number field, n = [K:Q] = r1 + 2 r2, d_K its discriminant, rd_K = |d_K|^{1/n} (root discriminant).
- Strip hypothesis H(w): every nontrivial zero rho = beta + i gamma of zeta_K satisfies |beta - 1/2| <= w.
- H(15/32) holds for every K: Cureton-Thorner (Hexagon 2610.00200v2, Thm 1.1; Lean development github.com/AxiomMath/DedekindZetaZeroFree) prove zeta_K(s) != 0 for Re s > 31/32; the functional equation Lambda_K(s) = Lambda_K(1 - s) maps a zero rho to 1 - rho, so also beta >= 1/32.
- H(w_ab), w_ab = 0.3733819666132076 [certify.py: '3733819666132076/10000000000000000'], holds for every abelian K: zeta_K(s) = prod_{chi in X(K)} L(s, chi*) over the n primitive characters of the character group X(K) (conductor-discriminant formula); Bebel's Lean-checked theorem (QuasiRiemannTracker, 31 Oct / 1 Nov 2026) gives L(s, chi) != 0 for Re s > 0.87338196661320752 for every Dirichlet character; for primitive chi, Lambda(s, chi) = eps(chi) Lambda(1 - s, conj chi), so if rho is a nontrivial zero of L(s, chi) then 1 - rho is one of L(s, conj chi), hence 1 - beta <= 0.87338196661320752. Since 0.87338196661320752 - 1/2 < w_ab, H(w_ab) holds.
- Conditional variants for every K: Liu's B_2 = (777 - sqrt 1081)/852 = 0.87338196661320745... (Hexagon 2611.00082, all finite-order Hecke L over every number field, unrefereed) gives H(w_ab) for every K, because zeta_K is itself a Hecke L-function; Liu's arXiv 2610.26835 (B_flat = 0.874886647040, unrefereed) gives H(0.374886647040) [certify.py: '374886647040/1000000000000'].
- H(1/2) holds trivially (no zeros on Re s = 0, 1); H(0) is GRH for zeta_K.

## 1. The explicit formula (Weil; Poitou's form)
Let F be even, real, continuous, compactly supported (support in [-X, X]), piecewise C^infinity and C^2, F(0) = 1. Put Phi(z) = int_R F(x) e^{(z - 1/2)x} dx (entire; Phi(c + it) << (1 + |t|)^{-2} uniformly for c in compacts).

Completed zeta: Lambda_K(s) = |d_K|^{s/2} Gamma_R(s)^{r1} Gamma_C(s)^{r2} zeta_K(s), Gamma_R(s) = pi^{-s/2} Gamma(s/2), Gamma_C(s) = 2 (2 pi)^{-s} Gamma(s); Lambda_K is meromorphic of order 1 with simple poles at 0, 1 and Lambda_K(s) = Lambda_K(1 - s).

(a) Contour step. For c > 1 the residue theorem on the rectangle [1 - c, c] x [-T, T] (T -> infinity along heights avoiding zeros; the horizontal integrals vanish by the standard bound Lambda'/Lambda(s) << log(|d_K| (|t| + 2)^n) off zeros and the decay of Phi) gives
  sum_rho Phi(rho) - Phi(0) - Phi(1) = (1/2 pi i) [int_(c) - int_(1-c)] Phi(z) (Lambda'/Lambda)(z) dz.
Since F is even, Phi(1 - z) = Phi(z), and Lambda'/Lambda(1 - z) = -Lambda'/Lambda(z); substituting z -> 1 - z in the left-line integral gives
  sum_rho Phi(rho) = Phi(0) + Phi(1) + (2/(2 pi i)) int_(c) Phi(z) (Lambda'/Lambda)(z) dz.   (1.1)
(b) Pieces of Lambda'/Lambda = (1/2) log|d_K| + r1 Gamma_R'/Gamma_R + r2 Gamma_C'/Gamma_C + zeta_K'/zeta_K:
  - constant: (2/(2 pi i)) int_(c) Phi(z) dz (1/2) log|d_K| = F(0) log|d_K| (Fourier inversion: Phi(1/2 + it) = Fhat(t), (1/2 pi) int Fhat = F(0); the line may be moved to Re z = 1/2 since Phi is entire and decays);
  - zeta_K'/zeta_K(z) = -sum_a Lambda_K(a) N a^{-z} (Re z = c > 1, absolutely convergent): (1/(2 pi i)) int_(c) Phi(z) e^{-u z} dz = e^{-u/2} F(u), so this piece is -2 sum_a Lambda_K(a) N a^{-1/2} F(log N a) = -2 sum_{p, m} (log N p) N p^{-m/2} F(m log N p);
  - Gamma pieces: moved to Re z = 1/2 (Gamma_R, Gamma_C are holomorphic on Re z > 0): (1/pi) int_R Fhat(t) [r1 Re (Gamma_R'/Gamma_R)(1/2 + it) + r2 Re (Gamma_C'/Gamma_C)(1/2 + it)] dt (the imaginary parts are odd in t and Fhat is even).
(c) Gamma pieces in x-space. With psi(z) = -gamma + int_0^infinity (e^{-x} - e^{-zx})/(1 - e^{-x}) dx (Re z > 0) and (1/pi) int_R Fhat(t) cos(tx) dt = 2 F(x):
  - (1/pi) int Fhat(t) Re psi(1/2 + it) dt = 2 [psi(1/2) + I_c(F)], I_c(F) = int_0^infinity (1 - F(x))/(2 sinh(x/2)) dx, because e^{-x/2}/(1 - e^{-x}) = 1/(2 sinh(x/2));
  - (1/pi) int Fhat(t) Re psi(1/4 + it/2) dt = 2 [psi(1/4) + I_c(F) + I_r(F)], I_r(F) = int_0^infinity (1 - F(x))/(2 cosh(x/2)) dx (substitute u = 2x in the psi integral; 2 e^{-x/2}/(1 - e^{-2x}) = 1/(2 sinh(x/2)) + 1/(2 cosh(x/2)));
  - Gamma_C'/Gamma_C(z) = -log 2 pi + psi(z), Gamma_R'/Gamma_R(z) = -(1/2) log pi + (1/2) psi(z/2), psi(1/2) = -gamma - 2 log 2, psi(1/4) = -gamma - pi/2 - 3 log 2.
(d) Collecting (1.1)-(c) (Poitou 1977; Odlyzko's survey, Sem. Theor. Nombres Bordeaux 2 (1990), eq. for log d):
  log|d_K| = n (gamma + log 8 pi - I_c(F)) + r1 (pi/2 - I_r(F)) - E(F) + sum_rho Phi(rho) + 2 sum_{p, m} (log N p) N p^{-m/2} F(m log N p),   (P)
  E(F) = Phi(0) + Phi(1) = 4 int_0^infinity F(x) cosh(x/2) dx.
Checks of (P), all in checks.py (logs/checks.log): (i) K = Q (n = r1 = 1, d = 1) with F = H(x/2) sech(15x/32): Gamma terms 2.9575358243, -E = -5.4052137685, 300 Riemann zeros 0.0000004638, primes 2.4476774803, residual -4.4e-14 (a second kernel: residual -3.8e-9, tail bound 8.6e-9); (ii) the t-integrals in (b) computed directly with mpmath's digamma agree with the x-space forms in (c) to 1e-13 and 1e-10 for two kernels; (iii) (P) with Odlyzko's kernels reproduces all 322 (b, A, B, E) rows of his 1976 Tables 3 (GRH) and 4 (unconditional): our A, B exceed his "lower estimates" by 3e-5 to 1e-3 relative.

## 2. Positivity under H(w)
Lemma A (cosh ratio). For |a| < w: int_R cosh(ax)/cosh(wx) e^{itx} dx = (2 pi / w) cos(pi a/(2w)) cosh(pi t/(2w)) / (cosh(pi t/w) + cos(pi a/w)) > 0. For |a| = w the ratio is 1.
Proof. int_R e^{zx}/cosh(x) dx = pi / cos(pi z/2) for |Re z| < 1 (residues). With alpha = pi a/(2w), tau = pi t/(2w): (1/w)(pi/2)[1/cos(alpha + i tau) + 1/cos(-alpha + i tau)] = (pi/w) 2 cos alpha cosh tau / (2 (cosh^2 tau - sin^2 alpha)) and cosh^2 tau - sin^2 alpha = (cosh 2 tau + cos 2 alpha)/2. Positivity: cos alpha > 0 and cosh 2 tau + cos 2 alpha > 0 for |alpha| < pi/2. QED
Lemma B (kernel shapes). G(x) = (1 - x/2) cos(pi x/2) + sin(pi x/2)/pi and H(x) = (2 - x)(1 + cos(pi x)/2)/3 + sin(pi x)/(2 pi) on [0, 2] (Odlyzko 1976, extended evenly, 0 for |x| > 2) are the normalised autocorrelations int h(t) h(t + x) dt / int h^2 of h_G(t) = cos(pi t/2) and h_H(t) = 1 + cos(pi t) on [-1, 1] (direct integration; normalisers 1 and 3). Hence 0 <= S <= 1, S(0) = 1, and Fhat_S = |hhat|^2/||h||^2 >= 0.
Proposition C. Let F(x) = S(x/b) sech(w' x), S in {G, H}, b > 0, w' >= w. Under H(w) every zero contributes Re Phi(rho) >= 0, and the prime terms are >= 0. Hence
  log|d_K| >= r1 log A + 2 r2 log B - E,  log B = gamma + log 8 pi - I_c(F),  log A = log B + pi/2 - I_r(F).   (T)
Proof. F >= 0 gives the primes. For the zeros: the multiset of zeros is closed under rho -> conj(rho), so sum_rho Phi(rho) = sum_rho Re Phi(rho) (absolutely convergent: Phi(beta + i gamma) << (1 + |gamma|)^{-2} uniformly in beta, and #{gamma in [T, T+1]} << log(|d_K| T^n)). With a = beta - 1/2 and F even, Re Phi(rho) = int F(x) e^{ax} cos(gamma x) dx = int F(x) cosh(ax) cos(gamma x) dx = (Fourier transform of x -> S(x/b) * cosh(ax)/cosh(w'x)) at gamma = (1/2 pi) (Fhat_{S(./b)} * FT[cosh(a.)/cosh(w'.)])(gamma) >= 0 by Lemmas A, B since |a| <= w <= w'. QED
Remarks. (1) In every optimisation run the best w' was w itself (64 runs, logs/opt_wp.log), so all tabulated kernels have w' = w. Optimising the shape over nonnegative autocorrelations h = q^2 (q a 3-term cosine series; kernopt.py, logs/kernopt.log) gains at most 0.1% (n = 8) and <= 0.002% for n >= 1000 over H, so the tabulated values are within 0.1% of the optimum of this kernel class. (2) Both A >= B and monotonicity: 0 <= F <= 1 gives I_r(F) <= int_0^infinity dx/(2 cosh(x/2)) = pi/2, so log A >= log B; therefore the formal r1 = 0 value (n log B - E)/n is a lower bound for every signature, and for a fixed kernel (log B - E/n) and (log A - E/n) increase with n: each tabulated value holds for all larger degrees.

## 3. Best possible constants of the method (asymptotic optimality)
Call F admissible for H(w) if F is even, F(0) = 1, F >= 0, F cosh(x/2) integrable, and x -> F(x) cosh(ax) has nonnegative Fourier transform for all |a| <= w (the condition under which (P) yields (T) whatever the zeros under H(w) are: a zero at 1/2 + a + i gamma with |a| <= w and any gamma is not excluded, so nonnegativity at every gamma is needed).
Proposition D. For admissible F: F(x) <= sech(wx) for all x; hence I_c(F) >= I_c(w) := I_c(sech(w.)), I_r(F) >= I_r(w), E(F) > 0, and the bound (T) is < A_r(w)^{r1/n} A_c(w)^{2 r2/n} for every n, with
  A_c(w) = exp(gamma + log 8 pi - I_c(w)),  A_r(w) = A_c(w) exp(pi/2 - I_r(w)).
Conversely F_b = H(x/b) sech(wx) is admissible and (T) with F_b tends to these constants (n -> infinity first, then b -> infinity): 1 - F_b = (1 - sech(wx)) + sech(wx)(1 - H(x/b)) with 0 <= 1 - H(t) <= C t^2 on [0, 2] (C = sup (1 - H(t))/t^2 < infinity since H(t) = 1 - (pi^2/6) t^2 + O(t^4)) and 1 - H <= 1, so the extra integrand is <= min(C x/b^2, 1/(2 sinh(x/2))) (resp. 1/(2 cosh(x/2))), which tends to 0 pointwise and is dominated by min(C x, 1/(2 sinh(x/2))) for b >= 1; hence I_c(F_b) -> I_c(w), I_r(F_b) -> I_r(w), while E(F_b)/n -> 0 for each fixed b. So rd_K >= (A_r(w) - o(1))^{r1/n} (A_c(w) - o(1))^{2 r2/n}, and no admissible kernel gives more.
Proof of F <= sech(wx): phi(x) = F(x) cosh(wx) is continuous, integrable, with nonnegative Fourier transform, hence positive definite (Bochner), so |phi(x)| <= phi(0) = 1. QED
Values (certify.py, Arb enclosures, logs/certify.log):
  w = 1/2:   I_c = log 2, I_r = pi/2 - 1; A_c = 4 pi e^gamma = 22.381616, A_r = 4 pi e^{1+gamma} = 60.839540 (Odlyzko-Poitou);
  w = 15/32: I_c = 0.6491685445, I_r = 0.5387691951; A_c = 23.387894, A_r = 65.643966;
  w = w_ab:  I_c = 0.5068009132, I_r = 0.4310016305; A_c = 26.966253, A_r = 84.299917;
  w = 0.374886647040 (Liu arXiv): A_c = 26.903236, A_r = 83.949939;
  w = 0 (GRH): A_c = 8 pi e^gamma = 44.763232, A_r = 8 pi e^{gamma + pi/2} = 215.332516.

## 4. Why not the s-family suggested in the brief
- The kernel F_s(x) = e^{-(s - 1/2)|x|} turns (P) into the Hadamard identity at real s:
  log|d_K| = r1 (log pi - psi(s/2)) + 2 r2 (log 2 pi - psi(s)) - 2/s - 2/(s - 1) + 2 sum_rho Re 1/(s - rho) + 2 (-zeta_K'/zeta_K)(s).
- For s > 1 both sums are >= 0 (Stark); positive mixtures int F_s dmu(s), mu >= 0 on s > 1, give a probability average of the single-s bounds, hence nothing beyond the best single s. Its limits are 2 pi e^gamma = 11.19 (all signatures) and 4 pi e^gamma = 22.38 (totally real) (checks.py s_only: n = 10^4 gives 10.792 and 21.408), far below the kernel method (22.050, 59.750 classical).
- For 31/32 < s <= 1 the identity holds but -zeta_K'/zeta_K(s) is no longer a convergent Dirichlet series: -zeta_K'/zeta_K(s) = 1/(s - 1) - gamma_K + O(s - 1) is negative near s = 1 and cancels the "positive pole term" 2/(1 - s); nothing on the prime side can be dropped. In (P) this shows up as E(F) = +infinity for any kernel with mu-mass on s <= 1 (F decays slower than e^{-|x|/2}). So the strip cannot be used through s < 1; it enters only through the positivity condition of Prop. C/D, whose linear program has the explicit asymptotic optimum sech(wx).

## 5. Abelian fields: the elementary conductor count
For abelian K, |d_K| = prod_{chi in X(K)} f_chi with n distinct primitive characters (all even if K is real; abelian fields are totally real or totally complex). The number of primitive characters of conductor q is phi*(q) = sum_{d | q} mu(q/d) phi(d), of which (phi*(q) + mu(q) + [2 | q] mu(q/2))/2 are even. Hence rd_K >= (product of the n smallest conductors, with these multiplicities)^{1/n} [tables.py: elementary]. It is sharp for n = 2 (3 and 5) and grows like sqrt n; among the tabulated degrees it beats the B-bound from n = 220 (all) and n = 1332 (totally real) (TABLES.md, column elem).

## 6. Certification (certify.py)
- Thresholds quoted in REPORT.md (rd > 22 for n >= 1581, totally real rd > 60 for n >= 1223, abelian/B_2 totally real rd > 81 for n >= 35448) are certified as Arb comparisons (THRESHOLD lines).
- Kernels F = S(x/b) sech(wx), S in {G, H}, b a decimal with <= 5 significant digits (kernels.json, chosen by the float optimiser make_kernels.py; the choice needs no proof).
- Arb (python-flint 0.9.0, 128-bit): I_c = int_delta^{2b} (1 - F)/(2 sinh(x/2)) dx + [0, M delta] - log tanh(b/2), delta = 2^{-24}, M = sup_{[0, delta]} |F'| (ball evaluation of F' = S'(x/b)/b sech(wx) - S(x/b) w tanh(wx) sech(wx); |1 - F(x)| <= M x and 2 sinh(x/2) >= x); I_r = int_0^{2b} (1 - F)/(2 cosh(x/2)) dx + pi/2 - 2 arctan(tanh(b/2)); E = 4 int_0^{2b} F cosh(x/2) dx; the tails use int_X^infinity dx/(2 sinh(x/2)) = -log tanh(X/4) and int_X^infinity dx/(2 cosh(x/2)) = pi/2 - 2 arctan(tanh(X/4)) with X = 2b, where F = 0. acb.integral is rigorous for these meromorphic integrands (poles at x = 2 pi i k and x = i pi (k + 1/2)/w are off [delta, 2b]).
- Each table entry: the ball v = (r1 log A + (n - r1) log B - E)/n; reported c = k/1000 is the largest with exp(v) > c verified by Arb comparison (floor3).
- Crossovers: exp(v) > 4 pi e^gamma resp. 4 pi e^{1 + gamma} checked as Arb comparisons (PASS lines).
- Limits: kernel sech(wx) on [delta, 200], [0, delta] piece in [0, w^2 delta^2/4] (1 - sech y <= y^2/2), tails in [0, -log tanh 50] resp. [0, pi/2 - 2 arctan tanh 50].

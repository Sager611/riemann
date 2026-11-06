# DERIVATION: deflated mixed-precision certificate for Weil positivity on a window (weildefl, 2026-11-06)

Window convention (TREE.md N1, same as Chekmachev Zenodo 23848597 and ../WeilZeta): c > 1, L = ln c, W = [-(ln c)/2, (ln c)/2],
  Q(f) = 2 <f, cosh(x/2)>^2 - 2 <f, sinh(x/2)>^2 + (1/2 pi) int (Re psi(1/4 + it/2) - log pi) |F(t)|^2 dt - 2 sum_{n<c} Lambda(n) n^{-1/2} Re (f * f~)(log n),
F(t) = int f e^{itx}; f in L^2(W). Statement certified: Q(f) >= eps ||f||^2 for all f in L^2(W) of the stated parity; both parities give all f
(Q is real, symmetric and commutes with x -> -x; complex f = u + iv splits). Notation below: |.| Euclidean norm, ||.||_F Frobenius, A >= B PSD order.

## 1. The finite reduction (unchanged from ../WeilZeta; proofs there)
(A1) Lattice identity and minorant (../WeilZeta/DERIVATION.md sec. 3-5; Chekmachev Thm 3.1, 5.10, Def. 6.1; scan RIGOR.md):
  Q(f) - eps ||f||^2 >= sum_j dl_j <f, phi_j>^2 + sum_{k in K} w_k <f, e_k>^2 +- 2 <f, v>^2 =: (P - N)(f),
  valid when mu(k) >= t1 = 1 for odd k not in K and mu(2j) >= 0 for j not in E (proved by the rigorous scan up to kscan and the monotone bound
  beyond; cert_zeta.rscan, unchanged). P = P0 + U D U^T (positive columns: e_k with w_k > 0, and +2 v v^T for even f), N = B B^T
  (columns sqrt(deff - dl_j) phi_j for dl_j <= 0, sqrt(-w_k) e_k for w_k <= 0, sqrt 2 v for odd f), deff = 1/2 - eps > 0, P0 diagonal in phi_j.
(A1') The high-frequency bound, re-derived here independently of ../WeilZeta (it is what replaces "Galerkin + tail bound" in full space):
  mu(k) = h(r_k) + 2(-1)^k tau(r_k) - 2 sum_{m<c} w_m cos(r_k ln m), h(r) = -log pi + Re psi(1/4 + i r/2), tau(r) = sum_{n>=0} c^{-b_n} b_n/(b_n^2 + r^2) >= 0.
  Gauss: Re psi(s + it) = -gamma + sum_n [1/(n+1) - (n + s)/((n + s)^2 + t^2)], each term increasing in t^2, so h is increasing on r >= 0; each term
  of tau decreases in r. Hence for k' >= k >= 1: mu(k') >= h(r_k) - 2 tau(r_k) - 2 sum |w_m| =: m(k) (code: SymZeta.monotone, which also adds the
  truncation remainder trem twice, a safe overcount). kscan := first k with m(k + 1) >= t1 (odd) or >= 0 (even), found by bisection in Arb
  (cert_zeta.scan_end); every node in (K2, kscan] is checked by the proved C scan (RIGOR.md). So mu >= t1 on all odd k not in K and mu >= 0
  on all even 2j not in E: the minorant drops only nonnegative terms, on the whole infinite-dimensional space. c = 16: kscan = 313,858 odd,
  115,462 even; c = 29: 88,378,021 and 32,512,457.
(A2) Woodbury: P > 0, so P - N >= 0 <=> T := I - B^T P^{-1} B >= 0, and B^T P^{-1} B = X - Y^T S^{-1} Y with
  X = B^T P0^{-1} B (m x m), Y = U^T P0^{-1} B (nU x m), S = D^{-1} + U^T P0^{-1} U (nU x nU) > 0. All entries are closed forms (Cauchy sums
  over E, ../WeilZeta/DERIVATION.md sec. 4-5); cert.py build() is cert_zeta.py main() up to and including this assembly, line for line.
  Proof of "<=>": for P > 0, x^T (P - B B^T) x >= 0 for all x  <=>  |B^T P^{-1/2} y|^2 <= |y|^2 for all y  <=>  B^T P^{-1} B <= I.
(A3) S >= D^{-1}, because U^T P0^{-1} U >= 0; hence ||S^{-1}||_2 <= max_i D_ii =: maxD (D = positive weights w_k, and 2 for the pole).

What ../WeilZeta does next: forms W_up = X - Y^T Z - Z^T Y + Z^T S Z (any Z, Chekmachev Lemma 5.3) and proves I - W_up >= nu > 0 by a Cholesky
of the whole m x m matrix at the full working precision (512 bits at c = 16). Here that step is replaced by sec. 2-5.

## 2. Two elementary facts
(F1) Woodbury lower bound (Chekmachev Lemma 5.3): for S > 0 and ANY matrices G, W of matching size,
  G^T S^{-1} G >= G^T W + W^T G - W^T S W,  since  (W - S^{-1} G)^T S (W - S^{-1} G) >= 0 expands to exactly the difference.
(F2) Shifted-Cholesky test (Chekmachev Lemma 5.5): for symmetric R, s0 > 0 and ANY lower-triangular L, put E = R - L L^T - s0 I.
  If max_i sum_j |E_ij| <= s0 - nu then lambda_min(R) >= s0 - ||E||_2 >= s0 - ||E||_inf >= nu  (||E||_2 <= ||E||_inf for symmetric E).
  L comes from a floating factorisation; only E is enclosed in Arb.

## 3. Deflation lemma (the new ingredient)
Lemma. Let T be real symmetric m x m, Yt any real m x k matrix of rank k, P_perp = I - Yt (Yt^T Yt)^{-1} Yt^T, delta > 0. If
  (i)  z^T T z >= delta |z|^2 for every z with Yt^T z = 0, and
  (ii) Sk := Yt^T T Yt - (P_perp T Yt)^T (P_perp T Yt)/delta >= 0,
then T >= 0. More generally (ii) may use any Lk <= Yt^T T Yt and any Pk with a^T Pk a >= |P_perp T Yt a|^2.
Proof. x = Yt a + z with z = P_perp x (so Yt^T z = 0). x^T T x = a^T Yt^T T Yt a + 2 z^T T Yt a + z^T T z, and z^T T Yt a = z^T (P_perp T Yt a)
(P_perp z = z, P_perp symmetric). By Cauchy-Schwarz and 2|u||w| <= |u|^2/delta + delta |w|^2:
  x^T T x >= a^T Yt^T T Yt a - |P_perp T Yt a|^2/delta + (z^T T z - delta |z|^2) >= a^T Sk a >= 0.  []
The lemma needs no basis of the complement and no exactness of Yt: Yt is whatever matrix the code holds (exact dyadic entries).
If Yt spans the invariant subspace of the k smallest eigenvalues of T exactly, P_perp T Yt = 0 and Sk = Yt^T T Yt: no loss. The loss is
quadratic in the residual |P_perp T Yt|, which is why Yt needs only about half the bits of the final margin (sec. 5).

## 4. Certifying (i) and (ii)
(i), variant WD_COMPL=arb, at LOW precision p_lo = 96 bits (the default is (i') below). By (F1) with G = Y, W = Ws (Ws = float64 solve of S Ws = Y, any matrix is allowed):
  T >= T_lo := I - X + Y^T Ws + Ws^T Y - Ws^T S Ws.
Put R_aug := T_lo + kappa Yt Yt^T - delta I with kappa = 1. If R_aug >= 0 then for z with Yt^T z = 0: z^T T z >= z^T T_lo z = z^T R_aug z + delta |z|^2
>= delta |z|^2, i.e. (i). R_aug >= 0 is proved by (F2): numpy Cholesky of mid(R_aug) - s0 I (s0 = lambda/2, then lambda/4, lambda/10, lambda = smallest
float64 eigenvalue of mid(R_aug)), nu = s0/2, residual enclosed in Arb at 96 bits. kappa Yt Yt^T lifts the deflated directions to about kappa, so
R_aug is well conditioned (its smallest eigenvalue is about delta = (first non-deflated float eigenvalue of T)/2, 0.03-0.3 in all runs).
(i'), default implementation (compl_double in cert.py): the same R_aug in binary64 with a rigorous error bound, so no O(m^3) Arb work.
  Assumption: IEEE binary64, round to nearest; numpy matmul = classical dgemm, so |fl(AB) - AB| <= gamma_n |A||B| entrywise (Higham Thm 3.5,
  any summation order, FMA allowed); gamma_n = n u/(1 - n u), u = 2^-53. WD_COMPL=arb runs the assumption-free 96-bit Arb version (c = 5-25: both pass).
  Inputs: Xd = fl(mid X) etc. with eX >= max_ij |X_ij - Xd_ij| (Arb: ball radius plus rounding), same eY, eS, and eYt for Yt -> Ytd; Ws = Wsd exact.
  Computed: A1 = fl(Yd^T Ws), B1 = fl(Sd Ws), A2 = fl(Ws^T B1), C = fl(Ytd Ytd^T), Rt = fl(I - Xd + A1 + A1^T - A2 + C - delta I), Rs = fl((Rt + Rt^T)/2).
  ||Rs - R_aug||_2 <= Etot = 1.01 (E_in + E_mm + E_sum), with (||.||_2 <= ||.||_F throughout; F-norms computed in binary64, inflated by 1 + 1e-12, rounded up)
    E_in  = m eX + 2 sqrt(nU m) eY ||Ws||_F + ||Ws||_F^2 nU eS + 2 (||Ytd||_F + sqrt(mk) eYt) sqrt(mk) eYt     (input errors; ||M||_F <= sqrt(#entries) max|M_ij|),
    E_mm  = 2 gamma_nU ||Yd||_F ||Ws||_F + gamma_nU ||Sd||_F ||Ws||_F^2 + gamma_nU ||Ws||_F ||B1||_F + gamma_k ||Ytd||_F^2   (products),
    E_sum = 8u (sqrt m + ||Xd||_F + 2||A1||_F + ||A2||_F + ||C||_F + delta sqrt m) + u ||Rt||_F       (entrywise sums of 7 terms, symmetrisation;
            (Rt + Rt^T)/2 - R_aug = (D + D^T)/2 with D = Rt - R_aug, so symmetrising does not increase the error).
  Then (F2) in binary64: Lc = cholesky(Rs - s0 I), Ed = fl(Rs - s0 I - fl(Lc Lc^T)), with ||(Rs - s0 I - Lc Lc^T) - Ed||_F <= gamma_m ||Lc||_F^2 +
  3u(||Rs||_F + s0 sqrt m + ||fl(Lc Lc^T)||_F) =: errE. lambda_min(R_aug) >= lambda_min(Rs) - Etot >= s0 - 1.01(||Ed||_F + errE) - Etot >= nu > 0
  (Lc Lc^T >= 0). Observed: Etot = 4e-11-5e-11 at c = 16, 1.5e-8 at c = 29, against s0 - nu = 0.03-0.008.
(ii), at HIGH precision p_hi, k columns only. With G := Y Yt (nU x k), any W (nU x k), Rk := G - S W:
  (ii-a) Yt^T T Yt = Yt^T Yt - Yt^T X Yt + G^T S^{-1} G >= Lk := Yt^T Yt - Yt^T X Yt + G^T W + W^T G - W^T S W      (F1).
  (ii-b) P_perp T Yt = P_perp(Yt - X Yt + Y^T S^{-1} G) = P_perp E1 + P_perp Y^T S^{-1} Rk,  E1 := Y^T W - X Yt   (P_perp Yt = 0, S^{-1} G = W + S^{-1} Rk),
         |P_perp Y^T S^{-1} Rk a| <= ||Y||_2 ||S^{-1}||_2 ||Rk||_2 |a| <= rho |a|,  rho := ||Y||_F maxD ||Rk||_F   (A3),
         so |P_perp T Yt a|^2 <= 2 |P_perp E1 a|^2 + 2 rho^2 |a|^2 =: a^T Pk a   ((u + w)^2 <= 2u^2 + 2w^2).
  Sk := Lk - Pk/delta (k x k, Arb at p_hi; (Yt^T Yt)^{-1} enclosed by Arb's inverse). Sk >= nu2 > 0 by (F2) with an mpmath Cholesky of
  mid(Sk) - s0 I (s0 = (smallest eigenvalue of mid(Sk), by inverse iteration)/2, /4, /10; nu2 = s0/2) and the Arb residual. (F2) with an
  absolute residual is robust for the graded block (entries 1e-2 ... 1e-124); a plain ball Cholesky of the raw block is not (c = 20 even:
  last pivot [0 +- 1e-68] at 384 bits), because the float64 eigenvectors mix the cluster below 1e-15 at random.
Then T >= 0 by the lemma, P - N >= 0 by (A2), and Q(f) >= eps ||f||^2 by (A1).

## 5. How Yt, W are produced (no rigor needed; only their quality matters)
- float64: Xd, Yd, Sd = midpoints; Td = I - Xd + Yd^T Sd^{-1} Yd; eigh(Td). k := #{eigenvalues < tau}, tau = 1e-2 (WD_KTAU); delta := ev[k]/2.
  The cluster below about 1e-15 is noise in float64, but its invariant SUBSPACE is well determined (gap: ev[k] = 0.06-0.25 in all runs).
- Refinement at p_mid = min(p_hi, p_hi/2 + 64) bits (correction equation, inverse-iteration type): TY = T Yt (with S^{-1} G by float64-preconditioned
  iterative refinement, warm-started), H = (Yt^T Yt)^{-1} Yt^T TY, PR = P_perp (TY - Yt H); solve T_perp E - E H = PR in float64 with
  K = Td + V_k V_k^T (Cholesky) by E <- K^{-1}(PR + E H), 10 sweeps, project E orthogonal to Yt, Yt <- mid(Yt - E). Each step gains 13-15 digits
  (logs: refine list). Stop when ||PR||_F < 2^-(p_hi/2 + 16) (then 2 ||PR||^2/delta is far below the margin, which is >= 2^-(p_hi - 64) by the
  choice of p_hi); S-solves stop at max(10^-6 target, 10^-17 x the previous ||PR||_F) (first step: 10^-13 x 10^-17), i.e. only as accurate as
  that step can use. Final W refined once more at p_hi to 10^-6 target.
- Bits: p_hi ~ log2(1/margin) + 64 (c = 16: 320; the margin 7.29e-77 = 2^-253). The entries of X, Y, S are O(1)-O(10^3) and Yt^T T Yt is a
  cancellation down to the margin, so (ii-a) needs absolute accuracy far below the margin: p_hi. (ii-b) and the refinement need only the square root:
  p_mid. (i) needs a few digits: binary64 with the rounding bound of (i') (default) or 96-bit Arb.

## 5b. Lean mode and checkpoints (engineering only)
- WD_LEAN=1: X, Y, S are never stored; MatGen evaluates the same closed forms (gram) in row blocks of 48 rows for every product, and one pass
  gives the binary64 copies with their error bounds and an upper bound of ||Y||_F. Same numbers, O(n) high-precision memory.
- WD_CKPT=<file>: the Cauchy sums F, F2 are saved as decimal strings with radius (python-flint arb.str(radius=True, more=True): the printed
  ball contains the computed ball) and re-read with arb(str); a resumed run therefore uses enclosures of the same quantities.

## 6. Constants (each appears in cert.py; rg-checked, logs/rg_constants.log)
t1 = 1 (WZ_T1); K2 = max(201, ceil(5 c^2/3)) made odd; NBO = 32, NBEK = 1, NBEJ = 3 (Chekmachev Def. 6.1 neighbour rule); scan flag margin 1e-6 and
E < margin/10 (cert_zeta.rscan); tau = 1e-2; delta = ev[k]/2; kappa = 1; p_lo = 96; p_mid = min(p_hi, p_hi//2 + 64); target = 2^-(p_hi//2 + 16);
S-solve tolerance max(target*1e-6, npr_prev*1e-17) with npr_prev = 1e-13 at the first step; 10 sweeps for E; factor 2 in Pk (theta = 1);
s0 factors 0.5, 0.25, 0.1 and nu = s0/2 in all (F2) tests; mpmath dps = 0.3 p_hi + 10; binary64 bound factors 8u, 3u, u, 1.01, 1 + 1e-12, gamma_n;
lean row blocks bs = 48; checkpoint strings with 0.302 p_hi + 15 digits.

## 7. Why the tail-integral ladder of explore17 cannot carry the full-space certificate (the reviewer's GAP, quantified)
Setting of the brief: V = span{Phi^{(j)}|_W}, B and Q_V from the tails t_j = Phi^{(j)} 1_{|u|>L/2} via Q(Phi^{(j)}|_W, h) = -Q_R(t_j, h)
(Q_R(Phi^{(j)}, h) = 0 for every h because Phi^ = Xi vanishes at every zero), Q_{V-perp} by Galerkin plus a high-frequency bound.
(O1) The coupling functions are rough. As a function on W, beta_j(z) := Q(Phi^{(j)}|_W, delta_z) = -Q_R(t_j, delta_z) contains
  +w_n [t_j(z + log n) + t_j(z - log n)] for every n >= 2, which jumps by w_n Phi^{(j)}(L/2) at z = +-(L/2 - log n) (inside W for n < c), and
  the archimedean part (1/2) int_{|y|>L/2} t_j(y) k(|y - z|) dy, k(x) = e^{x/2}/sinh x ~ 1/x, which is ~ (Phi^{(j)}(L/2)/2) log(1/d) at
  distance d from +-L/2. So the Fourier coefficients of beta_j decay like Phi^{(j)}(L/2)/omega and ||beta_j above frequency Omega||^2 ~
  (sum of squared jumps) / Omega. The ladder vector itself (extended by 0) jumps by Phi^{(j)}(L/2) at +-L/2.
(O2) The exact reduction Q >= 0 <=> [Q_{V-perp} > 0 and S_V = Q_V - B Q_{V-perp}^{-1} B^T >= 0] requires B Q_{V-perp}^{-1} B^T to absolute
  accuracy below lambda_min. Its part carried by frequencies above Omega is about ||beta_high||^2/log(Omega), which is NOT small against lambda_min
  (table in REPORT.md: e.g. c = 16, Omega = 1000: 1e-38 against 3.7e-75): it must be computed to relative accuracy ~ e^{-3.7c}, not bounded.
  A Galerkin cut-off plus a coercivity/tail bound gives only a bound, so S_V's lower bound goes negative unless Omega ~ e^{3.7c}
  (numbers in REPORT.md, from obstruct.py).
(O3) In the only full-space framework available (the lattice minorant (A1), which already needs the bad nodes up to kscan = 3.1e5 at c = 16), the
  minorant drops R(f) = (1/L) sum_{odd k not in K} (mu(k) - 1) <f, e_k>^2 + sum_{j not in E} (dl_j - deff) <f, phi_j>^2. For the raw ladder vector
  <v, e_k> ~ 2 v(L/2)/r_k, so R(v, v) ~ (2 v(L/2)^2/pi)(log(r/2pi) - 1)/r at r = r_{K2}: 1.2e-37 at c = 16, of which the part beyond kscan
  (4.7e-40) involves mu(k) known only to O(1) relative accuracy (REPORT.md): M(v, v) = Q(v, v) - eps ||v||^2 - R(v, v) is uncomputable at the
  needed accuracy 1e-77. The true near-null vectors are fine
  (R(psi) < lambda: that is why the minorant works), so the ladder must first be corrected at the ~e^{-5c} level, which is the high-precision step.
(O4) Boundary-adapted ladders P(D)Phi|_W vanishing to order p at +-L/2 (tested in e1_adapted.py, float64 + Arb jets; logs/e1_adapted.log) cut the jump but their
  spans become nearly dependent (Gram condition 6e17 at c = 10, p = 4) and some unit vectors carry boundary derivatives 1.9e6 (p = 4) and 3.5e12 (p = 8);
  the resulting worst-case drop estimate R ~ 5e-13 (p = 4), 7e-23 (p = 8) at K2 = 427 is still far above lambda_min(10) = 1.9e-43 (U(10)).
  (A graded analysis direction by direction was not done.)
Conclusion: the deflation idea is sound and pays only in the exact finite reduction (sec. 3-5), with the near-null subspace taken from T itself
(float64 eigensolve + refinement). The Phi-ladder explains that subspace: k(c) = 4, 9, 15, 19, 24, 28 (even) at c = 5, 10, 16, 20, 25, 29 for tau = 1e-2,
close to explore17's k ~ c, and at c = 16 the cosines of the principal angles between span{B^T Phi^(2i)|_W : i < 19} (B^T psi is the T-vector of a
near-null function psi, sec. 1) and span(Yt) are all >= 0.84 (ladder.py, logs/ladder.log). But its tail identity cannot replace the high-precision
information.
(O5) Scan growth (why no linear algebra reaches the record range by this route): kscan = 2 L e^{2 sigma + 1} reproduces the code's kscan exactly
(sigma = sum_{n<c} Lambda(n)/sqrt n; h(r) ~ log(r/2pi)); 5.1e10 at c = 47, 2.5e13 at c = 71, 1.2e16 at c = 100.

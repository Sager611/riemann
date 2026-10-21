# PROOF.md: Weil positivity on [-log 4, log 4], checked step by step (2026-10-20)

**Theorem (computer-assisted, Arb ball arithmetic).** For every complex C^2 function f supported in I = [-L, L], L = log 4 (mu = e^{2L} = 16):
Q(f) >= 10^-76 ||f||^2.
**Not a record.** Elston (Zenodo 23174403, 2026-10-06) already claims a certified larger window, half-width 1.495 (section 8).
Code: `bc.py` in this folder. Logs from the 2026-10-20 end-to-end run: `run2/` (schur, grams, band, cert, consistency check, `cert_16_summary.json`).

## 0. Conventions
F(t) = int f(x) e^{itx} dx, so (1/2pi) int |F|^2 = ||f||^2. a(t) = Re psi(1/4 + it/2) - log pi. S_s f(x) = 1_I(x-s) f(x-s).
C_p = sum_{n<16} Lambda(n) n^{-1/2} (S_{log n} + S_{-log n}), with n = 2,3,4,5,7,8,9,11,13.
Q(f) = 2|int f cosh(x/2)|^2 - 2|int f sinh(x/2)|^2 + (1/2pi) int a |F|^2 - <f, C_p f>. This is Weil's explicit formula for g = f * f~, in the same form as Mori VII eq. (Q) and Zhu.
It is checked against zeta zeros to 1e-28 (`../Weil/cert/check_ef.py`) and against an independent x-space matrix (`cons16.py`, section 7).
Every operator below commutes with x -> -x. So Q splits into even and odd parts with no cross term, and each part is certified on its own (700 modes each). The matrices are real symmetric, so positive semidefinite over R implies the same over C.

## 1. The chain of inequalities
- (A) ||C_p|| <= s = 3.6421 (section 2).
- (B) gamma = 9/2 and M = gamma - C_p, so m := gamma - s <= M <= b := gamma + s = 8.1421. Take Omega = 640 and c = Omega L = 887.23. Then C := a(Omega) - gamma = 0.12359101 > 0. The function a is even and increasing on t >= 0 (each term of the psi series decreases), so a - gamma >= C for |t| >= Omega.
  Let B be the band projection compressed to I, <f,Bf> = (1/2pi) int_{|t|<=Omega} |F|^2, so 0 <= B <= I. Then T_tail(f) := (1/2pi) int_{|t|>Omega} (a - gamma)|F|^2 >= C <f,(I-B)f>.
  Cauchy-Schwarz for the form <.,(I-B).> on each parity gives C <f,(I-B)f> >= u_0 |<h_0,f>|^2 + u_1 |<h_1,f>|^2. Here h_j = (I-B)v_j, v_0 = (2L)^{-1/2}, v_1 = (3/2L)^{1/2} x/L, and u_j E_j <= C with E_j >= <v_j,(I-B)v_j>.
  The bounds E_0 = (c+1)/(pi c^2) and E_1 = (3c^2+6c+2)/(pi c^3) come from integration by parts:
  <v_0,(I-B)v_0> = (2/(pi L)) int_Omega^inf sin^2(Lt)/t^2 dt <= (1/pi)(1/c + 1/c^2), and
  <v_1,(I-B)v_1> = (6/pi) int_c^inf (sin s - s cos s)^2/s^4 ds <= (6/pi)(1/(2c) + 1/c^2 + 1/(3c^3)).
  With u = (343, 114): u_0 E_0 = 0.12319641 <= C and u_1 E_1 = 0.12297563 <= C. Hence Q(f) >= <f,(M+V)f>, where V = K + U, K = poles + band part with symbol (a - gamma) 1_{|t|<=Omega}, and U = sum_j u_j |h_j><h_j|. No prolate defect is needed.
- (C) If M_tau + V >= 0, with M_tau := M - tau and tau = 10^-76, then Q(f) >= tau ||f||^2. The steps below are for M_tau, which satisfies m_t := gamma - s - tau = 0.8579 <= M_tau <= b - tau.
- (D) Legendre basis v_j = ((2j+1)/2L)^{1/2} P_j(x/L), with E: C^N -> L^2(I) and P = EE*, N = 1400. Put A = E*M_tau E, B2 = E*M_tau^2 E, J = E*VE, G = E*M_tau^{-1} E and D = B2 - A^2 = E*M_tau(I-P)M_tau E >= 0. Block elimination (section 3): if G^{-1} + J >= alpha I + beta D, then M_tau + V >= 0.
- (E) For every matrix X, W(X) = X + X* - X*AX + m_t^{-1}(I - AX - X*A + X*B2X) satisfies W - G = R_X*(m_t^{-1} - M_tau^{-1})R_X >= 0, where R_X = E - M_tau E X. So W >= G > 0 and W^{-1} <= G^{-1}. The code uses X = midpoint of (B2 - m_t A)^{-1}(A - m_t I).
- (F) Let Z = J_quad - eps_J I - alpha I - beta D (section 4). Then G^{-1} + J - alpha I - beta D >= W^{-1} + Z =: T. Also T >= 0 if and only if H := W + WZW = W T W >= 0.
- (G) H >= 0 is proved by interval Cholesky (section 5).

## 2. Gap 1: the Schur test now checks every cell
**Lemma.** Suppose w > 0 and (C_p w)(x) = sum_s c_s 1_I(x-s) w(x-s) <= s w(x) for almost every x. Then ||C_p|| <= s.
**Proof.** Start from |g(x)||f(x-s)| <= (1/2)(|g(x)|^2 w(x-s)/w(x) + |f(x-s)|^2 w(x)/w(x-s)). Multiply by c_s 1_I(x) 1_I(x-s), sum over the symmetric set of shifts (c_s = c_{-s}), integrate, and substitute y = x - s in the second sum. This gives |<g,C_p f>| <= (s/2)(||g||^2 + ||f||^2).
**Check (`bc.py schur`, which every cert run now repeats from the weight file).** The weight w is a step function on M = 32000 equal cells of u = x/L. For each cell i, R_i = sum over the shifts +-t_n of c_n max{w_j : cell j meets (cell i +- t_n) and [-1,1]} / w_i. This bounds the essential supremum of C_p w / w on cell i. The index ranges are widened by 1e-6 of a cell, and the shifts t = 1/2, 1, 3/2 (n = 2, 4, 8) are exact rationals.
**Result.** All 32000 cells are checked; the old code checked only i >= M/2. max_i R_i <= 3.642073270 <= s, and the maximum sits in cell 10506, which is in the left half. The weight is also exactly even (w_i = w_{M-1-i} as doubles), so the old half check was right, but the proof no longer depends on that.

## 3. Gap 2: Legendre truncation (e, n, h, alpha, beta)
**Chebyshev bound (Trefethen, ATAP Thm 8.2).** Let g be analytic with |g| <= M_g on the Bernstein ellipse E_rho of [-L,L]. Then ||(I-P)g||_{L^2(I)} <= (2L)^{1/2} 2 M_g rho^{1-N}/(rho-1).
On E_rho, |Im z| <= L(rho - 1/rho)/2, so for |t| <= Omega:
||(I-P)e^{it.}|| <= r_N := (8L rho^2/(rho-1)^2 e^{c(rho-1/rho)} rho^{-2N})^{1/2} = 3.07e-154 (rho = 3).
Let R = V - EJE*, so that PRP = 0, (I-P)RP = (I-P)VP and (I-P)R(I-P) = (I-P)V(I-P). Bound each part of V separately:
- **Band part K_b.** For ||f|| = 1, K_b f extends to an entire function with |K_b f(z)| <= (1/2pi) int_{|t|<=Omega} |a-gamma||F_f| e^{Omega|Im z|} <= K_0 (2L)^{1/2} e^{Omega|Im z|}, where K_0 = (Omega/pi) max(|a(0)-gamma|, |a(Omega)-gamma|) and |F_f| <= ||f||_1 <= (2L)^{1/2}. So ||(I-P)K_b P|| <= K_0 (2L)^{1/2} r_N = 1.03e-150.
  If z is orthogonal to all polynomials of degree < N, then |F_z(t)| <= ||(I-P)e^{-it.}|| ||z|| <= r_N ||z||, so |<z,K_b z>| <= K_0 r_N^2 ||z||^2 = 1.89e-304 ||z||^2.
- **U.** (I-P)h_j = -(I-P)Bv_j, and |Bv_j(z)| <= (Omega/pi)^{1/2} e^{Omega|Im z|} by Cauchy-Schwarz and Plancherel. Also ||h_j|| <= E_j^{1/2}, because 0 <= I-B <= I. Hence ||(I-P)UP|| <= (Omega/pi)^{1/2} r_N sum_j u_j E_j^{1/2} = 4.48e-152, (I-P)U(I-P) >= 0, and ||(I-P)U(I-P)|| <= (Omega/pi) r_N^2 (u_0 + u_1).
- **Poles.** The Taylor remainder gives ||(I-P)e^{+-x/2}|| <= q_N = (2L)^{1/2} e^{L/2} (L/2)^N/N!, and ||cosh(x/2)||, ||sinh(x/2)|| <= H_p = (2L)^{1/2} cosh(L/2). So ||(I-P)K_p P|| <= 4 H_p q_N = 1.15e-4020, and the norm of (I-P)K_p(I-P) is at most 4 q_N^2.
- **Sums (Mori VII Lemma 3.2).** e = 1.07e-150, n = 1.89e-304 < m_t, h = 8.94e-303. All are evaluated in Arb and printed by every cert run.
**Block elimination, with proof.** Write x = Ey + z with z in Ran(I-P). Put C' = (I-P)M_tau E and M22 = (I-P)M_tau(I-P) >= m_t, and split z = z* + w with z* = -M22^{-1} C' y.
Then <x,M_tau x> = <y,G^{-1}y> + <w,M22 w> (Schur complement), and ||z*||^2 <= <y,Dy>/m_t^2.
Also <x,Vx> >= <y,Jy> - 2e||y||(||z*|| + ||w||) - n||z*||^2 - 2h||z*|| ||w|| - n||w||^2.
Minimise over ||w|| (this uses m_t > n) and use 2e||y|| ||z*|| <= e(||y||^2 + ||z*||^2). The result is <x,(M_tau+V)x> >= <y,(G^{-1} + J - alpha - beta D)y>, with
alpha = e + 2e^2/(m_t - n) = 1.07e-150 and beta = (e+n)/m_t^2 + 2h^2/(m_t^2(m_t - n)) = 1.46e-150.
**Size of the correction.** alpha I + beta D is subtracted exactly inside Z. Its norm is at most alpha + beta (b - tau)^2 = 9.76e-149, which is 1e-72 of tau. (D <= ||M_tau||^2 <= (b - tau)^2.)

## 4. Gap 3: quadrature errors, their size and where they enter
- **A0 = E*C_p E and B0 = E*C_p^2 E have no quadrature error.** On each piece of [0,1] between the breakpoints |1 - t_n|, the integrands are polynomials of degree <= 2N-2, so (N+1)-point Gauss is exact. Nodes, weights and Legendre values are Arb balls, and the largest radius is <= 1.02e-113. The two pieces of length exactly 0 (repeated rational breakpoints 0 and 1/2) are now skipped only when the length is the exact ball 0. Before, they were skipped whenever the length was below 1e-40.
- **The pole vector has a closed form.** Av_j = (2(2j+1)L)^{1/2} (pi/L)^{1/2} I_{j+1/2}(L/2). E_j, C, K_0, r_N and q_N are closed forms too.
- **Band matrices K (weight a - gamma) and Bb (weight 1) on [0, Omega].** These use composite Gauss: 17 panels, 2376 nodes. The error bound for one normalised entry has three parts:
  - Gauss error: the sum over panels of h(64/15) rho^{-2n}/(rho^2 - 1) (2L/pi) e^{2LB} max|g| on the ellipse boundary (ATAP Thm 19.3). This uses |j_k(z)| <= e^{|Im z|}/(2k+1)^{1/2}. The maximum of |a_ext - gamma| comes from ball evaluation on arcs, and the ellipses avoid the digamma poles +-i(2k+1/2).
  - Rounding nodes to points: 2h x (max node radius) x sup|F'|, with sup|F'| <= (2L/pi)(9 + 2L sup|a - gamma|), because |a'| <= psi'(1/4)/2 < 9.
  - Flushing Bessel values below 2^-400.
  - Totals: errK = 5.59e-105 + 2.65e-109 + 2.00e-113 = 5.59e-105, and errB = 6.32e-105.
- **How entry errors become operator errors.** ||Delta||_2 <= ||Delta||_F <= d max|Delta_ik|. The only place these errors enter the proof is J: ||J_true - J_quad||_2 <= eps_J = d errK + u_j(2 e_v + e_v^2), with e_v = d errB = 4.42e-102. The U vector e_0 - Bb e_0 has true norm <= 1.
  eps_J = 3.92e-102 + 3.03e-99 = 3.04e-99 (even) and 1.01e-99 (odd). Z uses J_quad - eps_J I <= J_true.

## 5. Gap 4: from the interval Cholesky to a norm bound and a margin
H_true is the exact matrix inside the ball matrix H, and H_mid is its midpoint. Take a floating-point Cholesky factor L of H_mid - sigma I, and enclose Res = H_mid - sigma I - L L^T in Arb. Then, for every real x,
x^T H_true x = sigma |x|^2 + |L^T x|^2 + x^T(Res + dH)x >= (sigma - ||Res + dH||_2)|x|^2.
Here ||Res + dH||_2 <= ||Res||_F + ||dH||_F <= (sum |Res_ik|^2)^{1/2} + (sum rad(H_ik)^2)^{1/2} =: eps_chol, using upper bounds of each entry, and |dH_ik| <= rad H_ik. This is a matrix 2-norm bound, not an entrywise one.
So lambda_min(H_true) >= sigma - eps_chol. Since T = W^{-1} H W^{-1}, lambda_min(T) >= (sigma - eps_chol)/||W||_2^2, and ||W||_2 <= min(||W||_inf, ||W||_F) because W is symmetric.
The chain then closes: G^{-1} + J_true - alpha - beta D >= T > 0, so M_tau + V >= 0 by section 3, so Q(f) >= <f,(M+V)f> >= tau ||f||^2.

| parity | sigma | ||Res||_F | ||rad H||_F | lambda_min(H) >= | ||W||_2 <= | lambda_min(T) >= |
|---|---|---|---|---|---|---|
| even | 2^-256 = 8.64e-78 | 8.04e-81 | 2.75e-110 | 8.628e-78 | 3.039 | 9.341e-79 |
| odd | 2^-240 = 5.66e-73 | 1.89e-82 | 2.39e-110 | 5.660e-73 | 2.986 | 6.346e-74 |

## 6. Gap 5: end-to-end re-run (2026-10-20, `run2/run16.sh`, one process at a time, nice 10)
| step | wall | CPU (user) | peak RSS |
|---|---|---|---|
| schur | 1 s | 1 s | small |
| grams | 1478 s | 977 s | 493 MB |
| band | 318 s | 217 s | 438 MB |
| cert | 323 s | 194 s | 784 MB |
The cert step re-verifies H1 (Schur) and checks H2 to H7, exiting non-zero if any fails. It prints every quantity above and writes `cert_16_summary.json`.
Inputs (sha256, first 16 hex digits): grams 31c8af4e83e8a636, band 2f730902e2bc4bb2, weight 37d59ebe68f3ccab. Software: python-flint 0.9.0, Python 3.14.7.
The margins are identical to the 2026-10-20 cert log in this folder (`cert_mu16_tau1e-76.log`). Smoke test of the new code at mu = 5: certified at >= 1e-20 (`run2/cert_5.log`).
Reproduce:
`python bc.py schur 16 w_16_32000.txt; python bc.py grams 16 1400 384; python bc.py band 16 1400 640 9/2 384; BC_WFILE=w_16_32000.txt BC_SUMMARY=s.json python bc.py cert 16 grams_16_1400_384.bin band_16_1400_640_9o2_384.bin 36421/10000 1/2 1/10^76 343 114 384 - 1400`
In the last command, write 1/10^76 out as `1/1` followed by 76 zeros.

## 7. Independent checks (not part of the proof)
- `cons16.py` (`run2/cons16.log`) compares the certificate matrices with the x-space Galerkin matrix of Q from `../Weil/cert/weil.py`, which uses separate code and the x-space form of the archimedean term, on the first 40 Legendre modes.
  - lambda_min(Q_40) = 8.0977e-30 against lambda_min((M+V)_40) = 7.9980e-30 (even), and 3.7642e-27 against 3.7479e-27 (odd).
  - lambda_min(Q_40 - (M+V)_40) = 7.3e-76 and 1.8e-75, both >= 0, as Q >= M + V requires.
  - Terms of order 1 cancel down to 1e-30 in both codes, so the normalisations agree.
- Earlier validation (REPORT.md): Mori VI's parameters (mu = 10) reproduce his margin, and the mu = 3 and mu = 5 checks against weil.py pass.

## 8. Novelty: the record does not stand
- **Elston**, "Weil positivity on finite windows across L-functions: certified windows, a pre-registered decay law, ...", Zenodo 23174403 (concept doi 10.5281/zenodo.23163863), 2026-10-06; code at github.com/zacharyelston/weil-windows.
  It uses the same functional in Zhu's normalisation, with f supported in [-a, a] and n < e^{2a}. The logs `data/connes/grid_certificate_a1495_{even,odd}.run.txt` (runs of 2026-10-02) read "RESULT CERTIFIED ... a = 1.495: Q(f) >= 3.50e-96 ||f||^2" (even) and 8.26e-92 (odd).
  Its README says "2.38 audited; 2.6 and 2.99 regression-tested". Half-width 1.495 is larger than log 4 = 1.3863, and the claim is two weeks older than ours. We have not checked Elston's certificate.
- **Mori.** VIII (half-width 1/2 log 15, 2026-10-10) is still Mori's largest. His later deposits, up to 2026-10-14, are upper bounds and quasi-RH results.
- **Others.** Liu reaches full widths 2 and 2.125, Zhu 1.6. The other Zenodo hits are smaller windows or other problems: Ursachi 0.805, Senke 77/64, Ibbotson (finite taper families only), Umihara "a = 3" (a finite Gram matrix under external inputs), Kastoris (the semilocal operator).
- **Searches run on 2026-10-20.** Zenodo API: "Weil positivity", Mori's records, and related queries. GitHub API: moriryota, luciferyu666, and a repository search for "weil positivity" (33 repositories). The arXiv API returned HTTP 429 three times, so it was not re-checked today.
- The REPORT's earlier statement "nothing beyond 1/2 log 15" missed Elston.

## 9. Trust assumptions
- Arb / python-flint 0.9.0 encloses correctly.
- The analytic lemmas above are derived here and re-derived from Mori VII and Liu, but are not peer reviewed.
- The explicit formula is standard.
- The Python code is not checked by a proof assistant.

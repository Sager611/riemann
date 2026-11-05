# DERIVATION.md -- the out-of-core certificate (weilooc, 2026-11-05)

Same theorem and proof chain as ../PROOF.md (sections 1-5) and cert4.py.  This file proves the four changes made for
the memory-lean pipeline (oocert.py): (1) the Legendre recurrence used for the Gram matrices, (2) the error budget of
fixed-point Gram matrices, (3) the division-free form of W and H, (4) exact fixed-point products and the final
inequality.  Notation as in PROOF.md: M_tau = gamma - tau - C_p, m = gamma - s - tau, A = E*M_tau E, B2 = E*M_tau^2 E.

## 1. Legendre values by a norm-one recurrence
Let x = cos(theta) in [-1, 1], z = e^{i theta}, c = x, s = sin(theta) >= 0, and D_j = P_j(x) - z P_{j-1}(x) (P_{-1} = 0).
From (j+1) P_{j+1} = (2j+1) x P_j - j P_{j-1} and x = (z + 1/z)/2:
  (j+1) D_{j+1} = [(j+1/2) zbar - z/2] P_j - j P_{j-1} = lambda_j D_j + (1/2)(1 - z^2) P_{j-1},  lambda_j = (j+1/2) zbar - z/2.
Since (1/2)(1 - z^2) = -i z s and Im D_j = -s P_{j-1}:  (j+1) D_{j+1} = lambda_j D_j + i z Im(D_j).
Writing D = a + i b (a, b real), lambda_j = j c - i (j+1) s and i(lambda_j + z) = i((j+1) c - i j s) give
  (a_{j+1}, b_{j+1}) = T_j (a_j, b_j),  T_j = diag(j/(j+1), 1) [[c, s], [-s, c]],  a_0 = 1, b_0 = 0,
and P_{j+1} = a_{j+1} + c P_j (real part of D_{j+1}).  ||T_j||_2 = 1 exactly (a contraction times a rotation).
Check: a_1 = 0, b_1 = -s, P_1 = c; a_2 = -s^2/2, P_2 = c^2 - s^2/2 = (3c^2 - 1)/2.
**Error bound (gvals.c).**  Fixed point with FV fractional bits; c = u is the exact dyadic node, s^ = sqrt(1-u^2) rounded,
|s^ - s| <= delta_s <= 2^-FV.  Each step rounds a once after the product (|.| <= 2^-FV-1), then once after the
multiplication by j/(j+1) (|.| <= 2^-FV-1), b once (<= 2^-FV-1): the rounding vector has 2-norm <= 1.12 2^-FV.
With e_j = (a^_j - a_j, b^_j - b_j): e_{j+1} = T_j e_j + (T^_j - T_j) d^_j + r_j, ||T^_j - T_j|| = delta_s,
||d^_j|| <= |D_j| + ||e_j|| <= 2 + ||e_j|| (|P_j| <= 1 on [-1, 1]).  So ||e_{j+1}|| <= (1 + delta_s)||e_j|| + 2 delta_s +
1.12 2^-FV, hence ||e_j|| <= 3.2 j 2^-FV (j delta_s <= 2^-400).  For P: eps_{j+1} = e_{j+1,a} + c eps_j + rho'_j with
|rho'_j| <= 2^-FV-1, |c| <= 1, eps_0 = 0, so |eps_j| <= sum_{i<=j} (3.2 i + 0.5) 2^-FV <= 2.1 (j+1)^2 2^-FV.
Node displacement: |P_j(u^) - P_j(u)| <= P_j'(1) |u^ - u| = j(j+1)/2 |u^ - u| (Markov; |P_j'| <= P_j'(1) on [-1, 1]).
Hence eta_P = 2.1 (N+1)^2 2^-FV + N(N+1)/2 delta_u for every j < N, and for C_p P_j(u) = sum_a c_a P_j(u + t_a):
eta_C = c_sum (2.1 (N+1)^2 2^-FV + N(N+1)/2 delta_y) + n_active max|c^_a - c_a| + 2^-FV-1, delta_y = delta_u + max|t^_a - t_a|.
(The old grams_c needed ball arithmetic at prec + 1.28 N + 64 bits because a ball recurrence wraps by up to 1 + sqrt 2
per step; here no ball is propagated, the a-priori bound replaces it.)  Validated: Gram matrices agree with grams_c within
the old radii at mu = 5 (N = 140) and L = 31/20 (N = 800).

## 2. Gram matrices from scaled values and their error
Per parity, nodes u_q (bound-mode plan of grams_c, pieces of [0, 1]) with weights w_q; sigma_qj = sqrt(w_q (2j+1)).
gvals writes Phat_qi = round_F(sigma^_qj P^_j(u_q)), Chat_qi = round_F(sigma^_qj C^_j(u_q)) (j = 2i + p), so the exact
integer products Phat^T Chat and Chat^T Chat are the normalised quadrature sums of A0 = E*C_pE and B0 = E*C_p^2E.
For X, Y in {P, C} with |X| <= M_X (M_P = 1, M_C = c_sum), |Xhat_qi - sigma_qi X_i(u_q)| <= sigma_qi eta_X + rho_X with
rho_X = 2^-F-1 + delta_sigma (M_X + eta_X):
  |sum_q Xhat_qi Yhat_qk - sum_q sigma_qi sigma_qk X_i Y_k| <= s_ik W1 [eta_X (M_Y + eta_Y) + M_X eta_Y]
        + sqrt(2j_i+1) SW (eta_X + M_X) rho_Y + sqrt(2j_k+1) SW rho_X (M_Y + eta_Y) + nq rho_X rho_Y,
s_ik = sqrt((2j_i+1)(2j_k+1)), W1 >= sum_q w_q, SW >= sum_q sqrt(w_q).  Frobenius norms: ||[s_ik]||_F = S1 = sum_i (2j_i+1),
||[sqrt(2j_i+1)]||_F = sqrt(d S1).  Quadrature (grams_c bound mode, per raw entry eA, eB per piece): S1 sum(eA).
Final rounding of the product: d 2^-F-1.  E_A, E_B (oocert.gv_budget) are the sums, bounds of ||. ||_2 <= ||. ||_F.

## 3. Division-free W and H
With dyadic gamma, s, tau, m = gamma - s - tau is dyadic.  A = (gamma-tau) I - A0s (A0s = sym A0),
B2 = (gamma-tau)^2 I - 2(gamma-tau) A0s + B0.  Then G1 := m I - A = A0s - s I and
G2 := B2 - m A = (gamma-tau) s I - (gamma-tau+s) A0s + B0, and for every X
  W' := m W(X) = I + (mI - A) X + X^T (mI - A) + X^T (B2 - mA) X = I + G1 X + X^T G1 + X^T G2 X,
the W(X) of PROOF.md (E), multiplied by m.  X = -G2^-1 G1 (the minimiser, = cert4's X) is computed by iterative
refinement (any X is valid).  H' := m^2 H = m^2 (W + W Z W) = m W' + W' Z W'.  Since H = W T W (T = W^-1 + Z):
  lambda_min(T) >= lambda_min(H)/||W||^2 = lambda_min(H')/||W'||^2   (m cancels).

## 4. Exact products, symmetric storage, Z and the final inequality
fxmm: every stored entry is round(V/2^F) of an EXACT integer V (multimodular, 45-48 primes < 2^21, dgemm exact since every
partial sum < 512 p^2 < 2^51, explicit CRT with a range check from the inputs' top digits): one rounding <= 2^-F-1 per entry.
Residue reduction in one float step is exact for |R| < 2^51 + 2^22: the computed quotient R fl(1/p) is within 2^-21.9 < 1/p
of R/p, so it has the right floor unless p | R (then the representative is 0 or p; the explicit CRT accepts [0, p]).
The CRT quotient floor(sum_i w_i/p_i) is exact because V/P lies in [2^-27, 1/4] and the float sum errs by < 2^-41.
Symmetric outputs store the upper blocks and mirror them.  If computed = S + N (S symmetric, N not), the mirrored matrix
equals S + mirror(N), ||mirror(N)||_2 <= ||mirror(N)||_F <= sqrt(2) ||N||_F; rounding residues are bounded in Frobenius.
- G1, G2: ||G1_true - G1~|| <= E_A + d rho, ||G2_true - G2~|| <= (gamma-tau+s) E_A + E_B + d rho (rho = 2^-F-1).
- W': computed = W'_true - (dG1 X + X^T dG1) - X^T dG2 X + r1 + r1^T + X^T r2 (dG symmetric):
  e_W' <= 2(||X|| E_G1 + d rho) + d rho + sqrt2 (||X|| (||X|| E_G2 + d rho) + d rho)   (looser than needed; valid).
- Z: Kc = sum_{g>0} Bhat^T Bhat - sum_{g<0} Bhat^T Bhat (exact), R2 = sg Av~ Av~^T + u Eh~ Eh~^T (exact),
  ||J_true - Kc - R2||_2 <= epsJ = d errK + E_valK + E_wg + |sg|(2||Av|| dAv + dAv^2) + u (2 ev + ev^2), ev = sqrt(d)(errB + rad b0 + rho)
  (||Eh_true|| <= 1, PROOF.md section 4).  Z~ := round(Kc + R2 - shift I) with shift >= epsJ + alpha + beta (b-tau)^2 + 2 d rho,
  so J_true - alpha I - beta D >= Z~ (D <= (b - tau)^2).  Z~ is stored exactly: it IS the Z of PROOF.md (F).
- H': ||H'_true - H'~||_2 <= e_H' = m e_W' + sqrt2 (e_W' ||Z|| (2||W'~|| + e_W') + ||W'~|| d rho + d rho),
  H'_true := m W'_true + W'_true Z~ W'_true.
- Cholesky: L (any data) with H'~ - sigma' I = L L^T + Res, Res computed exactly by CRT (no rounding), ||Res||_F bounded
  from float64 conversions with relative error 2^-45.  lambda_min(H'_true) >= sigma' - ||Res||_F - e_H', and
  lambda_min(W^-1 + Z) >= (sigma' - ||Res||_F - e_H') / (||W'~||_2 + e_W')^2 > 0 closes the chain of PROOF.md.

## 5. Band
bandplan.py = cert4.band_lean's plan (panels, ellipses, node counts, errK = Gauss + node rounding, errB likewise; no flush;
per-entry target 2^-band_target_bits; Gauss nodes and weights at 64 bits above the entry precision).
bvals.c: Bessel values by cert2.sbessel's algorithm in C (balls at F + x_max + 128 bits); Bhat_qi = round_F(
sqrt((2L/pi)(2j+1) w^_q |g^_q|) s_j(t_q)) with w^, g^ the ball midpoints: entry error maxerr (ball radius + 2^-F-1).
E_valK = maxerr (2 sigsum sqrt(d S1) + d nq maxerr), sigsum = sqrt((2L/pi) nq Omega (supg + rad g)) >= sum_q sqrt((2L/pi) w^ |g^|);
E_wg = (2L/pi) S1 nq (rad w (supg + rad g) + w_max rad g) covers w g -> w^ g^ (|s_j| <= 1).  b0 in ball arithmetic + errB per entry.

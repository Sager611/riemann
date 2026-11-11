# DERIVATION (explore28): Chekmachev's matrix stage as an N-body problem

Notation as in ../WeilZeta/DERIVATION.md sec. 4-5. Columns c_a: U = {e_k : k in K+} (+ v, even f), B = {phi_j : dl_j <= 0} u {e_k : k in K-} (+ v, odd f).
P0 = deff Id + sum_{j in G} (dl_j - deff) phi_j phi_j^T, G = {j in E : dl_j > 0}; weights D_a (U: w_k, pole 2), s_a (B: deff - dl_j, -w_k, pole 2).

## Lemma P (particle form)
Points y and charges (u, w): e_k: y = k^2, u = beta_k, w = beta_k F(k^2); v: y = y_v = -L^2/(4 pi^2), u = beta_v, w = beta_v F(y_v) + kappa/deff
(kappa = (L/pi) cosh(L/4) even, (L/pi)^2 cosh(L/4) odd); phi_j: y = 4 j^2, u = 0, w = A_j/deff. Here F(y) = sum_{j in G} gam_j/(y - 4j^2),
gam_j = (1/dl_j - 1/deff) A_j^2. Then for a != b
  G_ab := <c_a, P0^{-1} c_b> = (u_a w_b - w_a u_b)/(y_a - y_b),
and G_aa = |c_a|^2/deff + beta_a^2 F2(y_a) (F2 = sum gam_j/(y - 4j^2)^2; phi_j: 1/deff).
Proof. P0^{-1} = Id/deff + sum_G (1/dl_j - 1/deff) phi_j phi_j^T. <phi_j, c> = A_j beta_c/(y_c - 4j^2) for c = e_k, v; partial fractions
1/((y-x)(y'-x)) = [1/(y'-x) - 1/(y-x)]/(y - y') give the F-part beta_a beta_b (F(y_b) - F(y_a))/(y_a - y_b). The Id/deff part: <e_k, e_k'> = 0
(k != k'), <e_k, v> = beta_k kappa/(y_k - y_v) (r_k^2 + 1/4 = (pi/L)^2 (k^2 - y_v)), <phi_j, c> = A_j beta_c/(y_c - 4j^2), <phi_j, phi_j'> = 0.
Each case equals (u_a w_b - w_a u_b)/(y_a - y_b) with the charges above (phi_j is the residue of e at y = 4j^2: u = 0). Checked numerically:
the structured LDL reproduces WeilZeta's dense margins to 10-15 digits (REPORT sec. 4).
So D_y G - G D_y = u w^T - w u^T: displacement rank 2. The same structure (H_qr = (u_q v_r - v_q u_r)/(lambda_q - lambda_r)) is the core of
Wiktorowicz-Wiktorowicz's pipeline (../WiktAudit/REPORT.md, rows 3, 5, 10); here it is shown for Chekmachev's minorant matrix.

## Lemma I (inertia form)
J'(nu) = [[S, Y], [Y^T, X - (1 - nu) Id]], S = D^{-1} + U^T P0^{-1} U, Y = U^T P0^{-1} B, X = B^T P0^{-1} B (B with columns sqrt(s_a) c_a).
Haynsworth: In J' = In S + In(X - (1-nu) Id - Y^T S^{-1} Y) = (nU, 0, 0) + In(-(Id - W - nu Id)), W = X - Y^T S^{-1} Y. Hence
  Id - W > nu Id  <=>  In J'(nu) = (nU, nB, 0).
Congruence by diag(Id, diag(1/sqrt s_a)) gives J''(nu) = G + diag(delta), delta_U = 1/D_a, delta_B = -(1 - nu)/s_a (no square roots);
then by diag(2^{e_a}) (exact Jacobi scaling) J_s. Both preserve inertia and the Cauchy-like form (charges scale by 2^{e_a}).

## Structured LDL^T (Schur/GKO recursion with the diagonal carried)
Pivot p: G_ap = (u_a w_p - w_a u_p)/(y_a - y_p), g_a = G_ap/d_p, u_a -= g_a u_p, w_a -= g_a w_p, d_a -= g_a G_ap (a > p).
Exactness: (y_a - y_b) G'_ab = R_ab - g_b R_ap - g_a R_pb = u'_a w'_b - w'_a u'_b (uses G_pb = G_bp). Memory 5n numbers; time ~ 4.5 n^2 operations.
Stabilisers (exact, inertia-preserving): symmetric pivoting (largest |d| within U, then within B); balancing the pair (u, w) by an exact dyadic
shear w -= x u and u *= 2^e, w *= 2^-e (det 1, R unchanged; Gu 1998: GKO loses accuracy when the generator pair is ill-conditioned).

## Lemma B (backward-error certificate; "point mode")
Run the recursion from exact points: before each step all remaining charges and diagonals are exact dyadics; the step is done in Arb, so the
radii bound the local error. With H_a = sum_{b != a} 1/|y_a - y_b| (all points; exact integer differences) every local error matrix E obeys
  ||E||_2 <= ||E||_inf <= max_a H_a (ru_a (W + RW) + |u_a| RW + rw_a (U + RU) + |w_a| RU) + rd_a
(U, W, RU, RW maxima over the remaining indices; E_aa = rd_a). Let eta = input error + sum of all step and shear errors.
(1) A_0 := M_0 + sum_p (F_p + E_p), F_p (shear rounding, support [p, n), zero diagonal) and E_p (step rounding, support (p, n)), has exactly the
computed pivots: eliminating p from M_p + F_p + (trailing terms) gives M_{p+1} + F_{p+1} + (trailing terms). So In A_0 = signs of the pivots.
(2) ||A_0 - J_s|| <= eta, so J_s - eta Id <= A_0 and #pos(J_s - eta Id) <= #pos(A_0).
(3) If #pos(A_0) = nU: sigma := min_U delta_a 2^{2 e_a} <= lambda_min(S_s) (the Gram part is PSD), S_s - eta > 0, and Haynsworth gives
Schur_B(J_s - eta Id) <= 0. With 0 <= (S_s - eta)^{-1} - S_s^{-1} <= eta/(sigma(sigma - eta)) and Schur_B(J_s) = -D (Id - W - nu Id) D,
D_a = 2^{e_a}/sqrt(s_a):
  lambda_min(Id - W) >= nu - kappa eta,   kappa = (1 + ||Y_s||^2/(sigma(sigma - eta))) max_B s_a 2^{-2 e_a},
with ||Y_s|| <= Frobenius norm of the entrywise bound (|u_a| + ru_a)(|w_b| + rw_b) + (|w_a| + rw_a)(|u_b| + ru_b) over |y_a - y_b|.
(4) If #pos(A_0) >= nU + 1: the same argument with J_s + eta Id gives lambda_min(Id - W) < nu + kappa eta.
So two runs give a rigorous two-sided enclosure. All bounds are O(n) per step (O(n^2) once for H and ||Y_s||).

## Lemma F (exact-integer Cauchy sums)
For k^2 and 4j^2 (integers < 2^68) D = y - x_j is an exact integer, |D| >= 1. With Gq_j = gam_j 2^P rounded (error <= 1 + rad_j 2^P units),
q = tdiv(Gq_j, D), q2 = tdiv(q, D): |F 2^P - sum q| <= sum_j (rad_j 2^P/|D_j| + 2), |F2 2^P - sum q2| <= sum_j (rad_j 2^P/|D_j| + 3). Added as radii.

## Ball mode and wrapping
The same recursion directly in Arb balls (inputs as balls) is rigorous too, but the radii are re-amplified at every step (dependency/wrapping):
digits lost grow from 8 (c = 10) to about 314 (c = 29 odd), so the needed precision grows with c. Lemma B replaces it (17 to 33 digits lost).

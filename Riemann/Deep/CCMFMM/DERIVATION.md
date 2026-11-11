# DERIVATION (ccmfmm): CCM's simple even ground state in the particle form

Parents: ../WeilSimple (criterion, deflated certificate, chain; its DERIVATION.md secs. 3-5) and ../ExploreFMM (particle form, inertia, backward error;
its DERIVATION.md). Notation: c = lambda^2, L = ln c, half-width a = log lambda = L/2; columns c_a, P0, G, deff, dl_j, w_k, s_a, beta, A_j, y_v,
F, F2 as in ../ExploreFMM/DERIVATION.md and ../WeilZeta/DERIVATION.md secs. 4-5; g = sum_{l <= NG} g_l phi_l (even, exact dyadic g_l), M > 0.

## 0. What is certified (../WeilSimple/DERIVATION.md sec. 3)
(a) U = Q(g)/||g||^2 >= lambda_1^even (g in E_NG; Arb from the closed-form block, gstate.py = cert_simple.py's ground_vec + rayleigh).
(b) odd sector at eps = L_odd: Q(f) - L_odd ||f||^2 >= (P - N)(f) >= 0 for odd f, decided by In J''(nu) (../ExploreFMM Lemmas P, I, B; unchanged).
(c) defl sector at eps = L_2: Q(f) - L_2 ||f||^2 + M <f, g>^2 >= (P + M g g^T - N)(f) >= 0 for even f (new here in the particle form).
U < min(L_odd, L_2) => spec(A_lambda) meets (-inf, min(L_odd, L_2)) in one simple eigenvalue with an even eigenvector (Courant-Fischer).

## Lemma P' (the deflation column)
P0 phi_l = p_l phi_l with p_l = dl_l (l in G) and p_l = deff (l not in G), since P0 = deff Id + sum_G (dl_j - deff) phi_j phi_j^T and
{0, ..., NG} lies in E (NG <= K2). With Phi(y) = sum_{l <= NG} g_l A_l/(p_l (y - 4 l^2)):
  r_b := <g, P0^{-1} c_b> = sum_l (g_l/p_l) <phi_l, c_b> = beta_b Phi(y_b)       (c_b = e_k: y_b = k^2; c_b = v: y_b = y_v),
  r_b = g_j/deff (c_b = phi_j, dl_j <= 0, j <= NG), r_b = 0 (j > NG),     <g, P0^{-1} g> = sum_l g_l^2/p_l.
Proof: <phi_l, e_k> = A_l beta_k/(k^2 - 4l^2), <phi_l, v> = A_l beta_v/(y_v - 4l^2), <phi_l, phi_j> = delta_lj (../WeilZeta sec. 4).
These are the inner products of cert_simple.py's gram() (sector defl) term by term. g is one more U column with weight M:
J'' gets one dense row/column: J''_gb = r_b (no 1/sqrt(s_b): the congruence of Lemma I keeps G unscaled), J''_gg = <g, P0^{-1} g> + 1/M.

## Lemma I' (inertia)
S' = D'^{-1} + U'^T P0^{-1} U' > 0 (D' = D plus M; Gram part PSD), so Lemma I holds with n_U + 1 positive indices:
  Id - B^T (P + M g g^T)^{-1} B > nu Id  <=>  In J''(nu) = (n_U + 1, n_B, 0).
Then P + M g g^T - N >= 0 (nu >= 0), and (c) holds; complex f split into real and imaginary parts (../WeilSimple sec. 3).

## Lemma R (g first: displacement rank 4)
Eliminate g first (pivot d_g = J''_gg >= 1/M > 0). The Schur complement on the particles is M = J''_pp - r r^T/d_g. With D_y = diag(y):
  D_y G - G D_y = u w^T - w u^T (Lemma P),   D_y (r r^T) - (r r^T) D_y = (D_y r) r^T - r (D_y r)^T,
so D_y M - M D_y = (u w^T - w u^T) + (u2 w2^T - w2 u2^T) with u2 = r, w2 = D_y r/d_g, i.e. for a != b
  M_ab = [(u_a w_b - w_a u_b) + (u2_a w2_b - w2_a u2_b)]/(y_a - y_b),   M_aa = d_a - r_a^2/d_g.
(Check: (u2_a w2_b - w2_a u2_b)/(y_a - y_b) = r_a r_b (y_b - y_a)/(d_g (y_a - y_b)) = -r_a r_b/d_g.) The Schur/GKO step of ../ExploreFMM
(pivot p, g_a = M_ap/d_p, every generator x -= g_a x_p, d_a -= g_a M_ap) is exact for any number of generator pairs: with R = sum of pairs,
R'_ab = R_ab - g_b R_ap - g_a R_pb (the g_a g_b R_pp term is 0) and g_b R_ap + g_a R_pb = (M_ap M_pb/d_p)(y_a - y_b) (M symmetric).
Memory 7n numbers; the stabilisers (pivoting within U, then B; exact dyadic shear and 2^e scaling) act on each pair separately (det 1 on
the pair, so its u w^T - w u^T is unchanged). Exact 2^{e_g} scaling of g (J''_gg -> 2^{2 e_g} J''_gg, r_a -> 2^{e_a + e_g} r_a) as for particles.

## Lemma B' (backward error with the border and two pairs)
As Lemma B, from exact points at every step. Changes:
- Input: row a of the input error gets + rad(r_a); row g is sum_a rad(r_a) + rad(J''_gg). eta_in = max(errnorm + max_a rad r_a, sum_a rad r_a + rad J''_gg).
- The g step is done in Arb from the recentred (exact) border; its recentring error E_g lives in the particle block (d_a and w2 radii) and is
  bounded by errnorm. Every later step: errnorm with both pairs, ||E||_inf <= max_a H_a sum_pairs (ru_a (W + RW) + |u_a| RW + rw_a (U + RU)
  + |w_a| RU) + rd_a (the bound of Lemma B applied to each pair; the pair-2 entries have the same 1/(y_a - y_b) kernel).
- A_0 := M_0 + E_g + sum_p (F_p + E_p) (E_g, F_p, E_p zero in row/column g) has exactly the computed pivots: eliminating g from A_0 gives
  Schur(M_0) + E_g + ... = M_1 + ..., and the rest is Lemma B (1).
- kappa: g is a U index, so sigma = min(min_U delta_a 2^{2 e_a}, delta_g 2^{2 e_g}) <= lambda_min(S_s) (S_s >= diag(delta_U) scaled since the
  Gram part is PSD; delta_g = 1/M); ||Y_s||_F^2 gets the g row sum_{b in B} (|r_b| + rad r_b)^2 (scaled entries).
Hence #pos = n_U + 1 => lambda_min(Id - W') >= nu - kappa eta; #pos >= n_U + 2 => lambda_min(Id - W') < nu + kappa eta.
Point domain (reviewer's gap, closed): H_a needs lower bounds of |y_a - y_b|. Integer pairs use exact fmpz differences (rounded toward 0).
The only non-integer point is the pole y_v = -L^2/(4 pi^2) in (-0.35, -0.22) for c in [20, 40]; every other point is an integer >= 0, so
|y_b - y_v| >= 0.22 and its binary64 evaluation |fl(y_b) - fl(mid y_v)| has relative error <= 3 * 2^-53 (y_b exact below 2^53, else
relative error 2^-53) plus absolute 2^-55: inside slt.c's margin (1 - 1e-12) x - 1e-12. (Two close huge non-integers, where double
rounding can overstate a distance, cannot occur.)
slt.c: SLT_BRD=<prep_ccm.py's .brd>; lines marked CCM. With no border the code path is ../ExploreFMM's (regression: c = 10 even, logs/val).

## Levels, windows, chain
Integer windows: L_2 = (e1 e2)^{1/2}, L_odd = (e1 o1)^{1/2} from approximate eigenvalues of the truncated blocks at N = NG (orientation only;
any level works if its certificate passes; c > 30: e1 := U, o1 by inverse iteration, L_2 := L_odd, no full eigendecomposition). Each certificate runs the LDL at a ladder nu = eps 10^k; the certified bound is the largest passing nu
(with nu > kappa eta), the first failing nu brackets lambda_min(Id - W) from above. Chain: L(c_{i+1}) = 2.00001 U(c_i) for both sectors at
c_{i+1}; monotonicity in c of every min-max value of each parity (../WeilSimple sec. 5) covers [c_i, c_{i+1}].
Scan premise at every window: rscan.c (../WeilZeta/RIGOR.md, proved rounding bound; flagged nodes re-checked in Arb), K2 = max(401, ceil(5c^2/3)).

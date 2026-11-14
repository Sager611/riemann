# explore48 derivation (appends of <= 60 lines)

## 1. Setup (23:10 UTC)
Unimodular rank-2 lattice L_z = (1/sqrt y)(Z + zZ), z = x + iy reduced, y >= sqrt(3)/2. Shortest vector lambda_1 = y^{-1/2};
aspect ratio lambda_2/lambda_1 = |z| ~ y. E*(z,s) = pi^{-s} Gamma(s) (1/2) sum' y^s |mz+n|^{-2s} = E*(z,1-s).
Chowla-Selberg = Poisson summation along the shortest vector (1-D smoothing of each line sum):
  E*(z,s) = a_y(s) + 4 sqrt(y) sum_{k>=1} k^{s-1/2} sigma_{1-2s}(k) K_{s-1/2}(2 pi k y) cos(2 pi k x),
  a_y(s) = zeta*(2s) y^s + zeta*(2s-1) y^{1-s}  (constant term; Lagarias-Suzuki 2006: for y >= 1 all its zeros lie on Re s = 1/2,
  except one real pair when y > 4 pi e^{-gamma} = 7.0555; matches BST 2110.09368's 1/Delta_c* = 1/0.141733 = 7.0555).
The k-th Bessel term is the k-th dual frequency of the Gaussian-smoothed line sum: the lattice-crypto 'smoothing error'.
Binary forms: Lambda_Q(s) = (sqrt|D|/2pi)^s Gamma(s) Z_Q(s) = 2 E*(z,s), z = (-b + i sqrt|D|)/(2a). Checked (t1_validate.py):
CS vs zeta*L products for D = -20, -15 agree to 1e-15; Potter-Titchmarsh zero of x^2+5y^2 = 0.93296970 + 15.66824953 i
(t/(2 pi y) = 1.115); explore21's x^2+xy+4y^2 zero 0.80001099 + 12.03859863 i (t/(2 pi y) = 0.989).

## 2. On the line
s = 1/2 + it: E* = 2 sqrt(y) |zeta*(1+2it)| F(t), F = cos theta_a + beta,
  theta_a = arg(zeta*(1+2it) y^{it}),  beta = 2 sqrt(cosh pi t) sum_k C_k(t) K_{it}(2 pi k y) cos(2 pi k x) / |zeta(1+2it)|,
  C_k(t) = sum_{d|k} cos(t log(k/d^2)) = k^{it} sigma_{-2it}(k) (real, |C_k| <= d(k)).
(|zeta*(1+2it)| = |zeta(1+2it)| / sqrt(cosh pi t).)

## 3. Lemma K (Gaussian tail of the dual frequencies; contour shift, rigorous)
For 0 <= nu < u: e^{pi nu/2} |K_{i nu}(u)| <= K_0(sqrt(u^2 - nu^2)) exp(nu arccos(nu/u)) <= sqrt(pi/(2w)) e^{-u g(nu/u)},
  w = sqrt(u^2-nu^2), g(a) = sqrt(1-a^2) - a arccos a (g(0)=1, g decreasing, g(1-e) ~ (2 sqrt2/3) e^{3/2}).
Proof: K_{i nu}(u) = (1/2) int_R exp(-u cosh v + i nu v) dv; shift v -> v + i alpha, sin alpha = nu/u (integrand entire, decays in
|Im v| < pi/2); |integrand| = e^{-nu alpha} e^{-u cos(alpha) cosh v}; int = 2 K_0(u cos alpha). Then pi/2 - alpha = arccos(nu/u).
So for t < 2 pi y: |beta(t)| <= sqrt2 (1+e^{-2pi t})^{1/2} |zeta(1+2it)|^{-1} sum_k d(k) sqrt(pi/(2 w_k)) e^{-2 pi k y g(t/(2 pi k y))}:
exponentially small until t = 2 pi y - O(y^{1/3}) (Airy window). Smoothing height of L_z: T_s(z) = 2 pi y = 2 pi / lambda_1^2.

## 4. Counting lemma (Rouche + sign changes) -> 'RH below the smoothing height'
If on [t1,t2] (t1,t2 extrema of cos theta_a) |beta| < 1 at every extremum of cos theta_a, F has >= one sign change per
half-period of theta_a, i.e. >= #zeros of a_y in the box. If moreover |E* - a_y| < |a_y| on the boundary of the box
[1/2-d, 1/2+d'] x [t1,t2] (true when the Lemma K bound is < 1 there), Rouche gives equality of counts, so every zero of E* in the
box is on the line and simple. Hence: RH(L_z) holds for 0 < |t| < T_1(z) := first t with |beta| >= 1 (sharp version), and
unconditionally for |t| <= 2 pi y - C y^{1/3} (log log y)^{2/3} via Lemma K and |zeta(1+2it)|^{-1} << log t.

## 5. Exact counting (count.py; validated)
N(T) = 1 + (1/pi)[Im lgamma(sR+iT) - T log pi + T log y + Arg Et(sR+iT) + D_hor arg G + pi(1/2-sR)/2] - N_real/2,
Et = y^{-s}E (|Et-1| <= Et(3)-1 <= 0.07 on Re s = 3 for the tested z), G = e^{-i pi s/2}E*, N_real = 2 iff y > 7.0555.
Derivation: argument principle for s(s-1)E*(z,s) on [1-sR, sR] x [-T, T], fourfold symmetry; arg s(s-1) contributes pi.
Off-line count = N(T) - N0(T), N0 = sign changes of F. Integers to 1e-4 at every checkpoint; D=-163 (h=1) gives 0 to 97.2.
Pitfall: N0 misses close on-line pairs (zeta x L(chi_D) has nearly coincident zeros) -> spurious '2'; fixed by local refinement.

## 6. Off-line onset data (t4_onset.py, logs/t4_*.log); T_s = 2 pi y = 2 pi / lambda_1^2
First off-line zero rho_1, ratio Im(rho_1)/T_s:
- h=2 principal forms (16 discriminants, y = 1.94..10.33): 0.982 .. 1.166 (D-15 0.989, D-24 0.982, D-232 0.997, D-427 1.034).
- D-20 fat form 2x^2+2xy+3y^2 (y = 1.118): 2.0794+7.6289i, 1.086 (Re > 1: DH zero).
- rectangular z = iy: y = 2.5..30: 0.954 .. 1.229; y = 40: 0.5354+522.2719i, 2.078 (no off-line zero in [0, 2 T_s]).
- Re(rho_1) - 1/2 shrinks with y: 0.30 (y=2.4) .. 0.03 (y=25); about 0.8/y.
- h=1 controls D=-163, -67, -43: |beta| reaches 2.7-3.5 above T_s, yet no off-line zero (pending refined rerun for -67, -43).

## 7. Refined reruns (local N0 refinement, logs/t4_h1r.log, t4_redo.log)
- h=1 controls (Euler products zeta*L(chi_D)): D=-163 (y=6.38), -67 (4.09), -43 (3.28): earlier '2' were missed close on-line
  pairs; after refinement 0 off-line zeros up to 3.25, 3.39, 3.49 T_s, although |beta| reaches 2.7, 3.5, 3.2 above T_s.
- D=-403 (h=2, y=10.04): spurious bracket at 1.03 T_s removed; true first off-line zero 0.575046+129.867571i = 2.059 T_s (k = 2).
- So the law is 'just above k T_s' with k = 1 or 2, not always k = 1.

## 8. Birth law: off-line zeros are born at the Airy extrema of the dual frequencies (analysis.py, logs/analysis.log)
Near t = u_k = 2 pi k y: e^{pi t/2} K_{it}(u_k) ~ pi (2/t)^{1/3} Ai(-(2/u_k)^{1/3} (t - u_k)) (turning point). Extrema of |Ai| at
the zeros a'_j of Ai': |a'_j| = 1.0188, 3.2482, 4.8201. Define xi = (Im rho_1 - k T_s)(2/(k T_s))^{1/3}, k = floor(ratio + 0.06).
- 44 first off-line zeros: k = 1 (40 cases, ratio 0.954..1.322), k = 2 (D-403, rect 15, 30, 40); none for h = 1.
- y >= 15: xi = 1.222, 1.016, 0.968, 0.943 (rhombic 15..30), 1.119, 0.877 (rect 20, 25) vs |a'_1| = 1.019;
  rect 30 (k=2) 0.835 ~ |a'_1|; rect 40 (k=2) 3.109 ~ |a'_2| = 3.248; rect 15 (k=2) 4.734 ~ |a'_3| = 4.820.
- small y (T_s < 40): xi spread -0.36..1.56 (hump far above threshold -> births earlier, even before the onset).
- rect y = 10, all births to 4.35 T_s (t6_dist.py): pairs born in T/T_s in (1.00,1.10] x2, (1.95,2.29] x4, (2.60,2.65],
  (3.05,3.10], (4.06,4.35] x3; NONE in (1.10,1.95), (2.65,3.05), (3.10,4.05): births quantized just above k T_s.

## 9. Every birth verified (t7_verify.py: exact count per interval, fine on-line scan, Newton), logs/t7_*.log
- rect y=10 (T_s = 62.83), 11 pairs, Im/T_s: 1.029 1.066 | 1.971 2.069 2.120 2.290 | 2.616 | 3.065 | 4.106 4.122 4.299.
- D-427 (h=2, T_s = 64.92), 13 pairs: 1.034 1.305 1.427 | 2.079 2.177 | 3.048 3.086 3.102 3.121 | 3.804 3.947 | 4.040 4.284.
- 20 of 24 births lie in [k - 0.06, k + 0.30] T_s (36% of the line): binomial tail P(>= 20 | 24, 0.36) = 2.7e-6.
- D-163 (h=1): the one jump (3.456-3.492 T_s) is 3 zeros, 3 on the line, 0 off: a missed close pair. GRH-consistent to 4 T_s.

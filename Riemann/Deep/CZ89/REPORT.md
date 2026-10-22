NO REFUTATION: Theorem 1 of Cobeli-Zaharescu (arXiv 2610.23467), liminf N_{s or 0}(T)/N(T) >= 0.888931 > 8/9, holds. Its one computer-checked step, Lemma 5.2 (S(g) >= 11/250 for all g >= 0), is proved here by an independent rigorous Arb branch-and-bound. But the paper's own computer evidence for that lemma does not hold up: the true margin is 6.4e-5, not about 3e-3, and the reported 252,035-box run cannot be reproduced with the tests the appendix describes.

## Quantity that must exceed 11/250
S(g) = sum_{0<=i<j<=6} w_ij K_*(y_j - y_i)^2 + (1/200) sum g_i, with y_k = g_1+...+g_k, weights from Table 1, and K_* the Fourier transform of f_* in (5.1).
eps = 11/250 gives alpha = 114 eps/120 = 209/5000, then tau = 2 + sqrt(2 - 2 alpha), then (6.2). Nothing else is numerical.

## Proof chain re-derived (every step checks)
- Lemma 4.1 (Hilbert-Schmidt square inequality), n_+(Q) <= b, (4.1) b <= (N-n)/2 - B/4, and tr phi(P) = n + R(H) give (3.1); A_tau = 1 + alpha.
- (5.9) clip, (5.10) block cut via Jensen/Peierls, Lemma 5.3: 114 windows, each diagonal sums to 2, each gap counted at most 6 times, each window in at most 114 blocks. This gives (5.11) with d = 57/2000. Cap (5.12): 5.016 < 5.0625 < (tau-1)^2 = 5.685.
- (6.1) is Lamzouri Lemma 3.1 = BGSTB Lemma 5 (checked in Lamzouri's TeX). It is applied to Q and Q''; for the smoothed f_delta both are supported strictly inside (-1,1). The weight-removal identity checks. Smoothing costs eps_delta = eps - 24 ||f_delta - f_*||_1, and the limit gives (6.2).
- Arb values: C = 1.33133947994646 (closed form agrees with 40-digit quadrature); U = 0.888931016871589; U - 8/9 = 4.2128e-5; f_* >= 0.79249 on I; K_*(1/2) = 0.675873162 > 2/3; K_*(0.743) = 0.380298365 > 0.38029.

## The true landscape of S
- Global minimum 0.0440644948 at g* = (1.08371, 1.89392, 1.04136, 1.02488, 1.02891, 1.02556) and at its reversal (sum 7.098). Arb: S(g*) - 11/250 = 6.4495e-5.
- Next basin 0.0440648987. At least 19 basins lie below 0.046 (188k multi-start runs).
- A commented-out line in the paper's TeX gives "Minimum S = 0.0469769, g = (1.877, 1.152, 1.879, 1.130, 1.825, 1.136), sum 9.000000, margin 2.98e-3". The published S at that g is 0.0489185, and on the plane sum = 9 the published S gets down to 0.0465045. So the authors' numerics were not evaluating the published S.

## Certificate from this audit (cert/, logs/)
- Arb tables (python-flint 0.9.0, 96-bit): enclosures of q = K_*^2 and q' at x = m 2^-13 (cross-check grid 2^-15) on [0,9], plus per-cell lower bounds for q and upper bounds for |q''|. The sinc singularities are handled by a series with an explicit tail bound.
- Branch-and-bound on {g >= 0, sum <= 8.8}; above sum 8.8 the linear term alone gives 0.044. Dyadic boxes are tightened by the sum constraint. Two lower-bound tests: L1 takes the range-minimum of q for each pair; L2 is a first-order Taylor bound at the box centre using a local |q''| bound on each pair's range. Every float sum carries an explicit 64u error term.
- Target 0.044: COMPLETE. 2,562,599 boxes on the 2^-13 grid and 2,821,419 on 2^-15; 0 unresolved; about 5 s each.
- Target 0.04406: COMPLETE (3,997,755 boxes). Negative control at target 0.04407: FAILS at the second basin, as it should.
- 2,913 random boxes were spot-checked: no computed lower bound exceeded the sampled minimum of S.
- Result: min over g >= 0 of S lies in [0.04406, 0.0440644948], so Lemma 5.2 is true.

## The paper's Appendix A run
- I reran with exactly the appendix's tests (L1 with exact infima, L2 with the global bound |q''| < 40). L2 never accepts a box. The cube of radius 0.002 around g* alone needs 9,882,657 boxes. On Omega' the run is still unfinished after 30,004,135 boxes. The paper reports 252,035 boxes in total.
- Most likely their code computed a different S (compare the commented-out minimum above). This is a gap in the published proof, and the certificate above closes it.

## Side remark (routine)
With eps = 0.04406 (certified), alpha = 0.041857 and the constant rises to 0.888948. This is not a discovery.


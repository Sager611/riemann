# RIGOR: the rigorous node scan (rscan.c) and its a-priori error bound (weilchi, 2026-11-05)

## What must be proved
The minorant (DERIVATION.md sec. 5) needs: every odd k in (K2, kscan] with mu(k) < t1 lies in K, and every even k = 2j in (2J, kscan_e] with
mu(2j) < 0 lies in E. (Beyond kscan, resp. kscan_e, the monotone bound proved in DERIVATION.md sec. 5 covers all nodes; it is evaluated in Arb.)
rscan.c computes, for EVERY node k of these ranges, a double mu_c(k) with the proved ONE-SIDED bound mu(k) >= mu_c(k) - E, and prints k
whenever mu_c(k) < thr + margin, where E < margin/10 (checked in cert_chi.py, function rscan). No absolute bound |mu - mu_c| <= E is
claimed: the block bounds H, T (below) bias mu_c downward by an amount that is not bounded here. The no-miss argument needs only the
one-sided bound. If mu(k) < thr, then mu_c(k) <= mu(k) + E < thr + margin, so k is printed. Printed nodes are re-evaluated in Arb (96 bits)
and selected if the Arb lower endpoint is below thr. A printed node with mu(k) < thr is always kept, because its Arb ball contains mu(k), so
the lower endpoint is <= mu(k) < thr. Selecting extra nodes is always allowed.
No libm function is called in the loop: only integer arithmetic and IEEE-754 binary64 +, -, * (correctly rounded, round to nearest,
compiled with -ffp-contract=off). u = 2^-53.

## The quantity computed
mu(k) = h(r_k) + 2(-1)^k tail(r_k) - 2 sum_n w_n cos(r_k log n), r_k = pi k/L, h(r) = log(q/pi) + Re psi(1/4 + a/2 + ir/2),
tail(r) = sum_m c^{-b_m} b_m/(b_m^2 + r^2) (+ remainder <= trem).
- Blocks. The k-range is cut into blocks [k0, k1] (same parity, k1 <= k0 (1 + 0.001)). For each block Arb gives H <= h(r_{k0}) and
  T >= 2 tail(r_{k0}) + 2 trem, outward-rounded to doubles (nextafter). Since h increases and tail decreases in r > 0 (DERIVATION.md
  sec. 5), h(r_k) + 2(-1)^k tail(r_k) >= H - T for every k in the block.
- Phases. alpha_n = log n/(2L), so r_k log n = 2 pi k alpha_n. B_n = round(2^64 alpha_n) (mpmath, 60 digits); Arb gives
  e_n >= |2^64 alpha_n - B_n| (e_n <= 0.4648 in all runs, asserted <= 0.5000001). phi = k B_n mod 2^64 is exact (uint64 wrap-around),
  and theta_used = 2 pi phi/2^64 satisfies |theta_used - 2 pi k alpha_n (mod 2 pi)| <= 2 pi k e_n 2^-64.
- Table. phi = A 2^48 + b with A < 2^16, b < 2^48. C_A, S_A are cos, sin(2 pi A/2^16) from Arb (128 bits) rounded to nearest:
  |C_A - cos| <= 2^-54 + 2^-120 <= 6e-17, same for S_A.
- Low part. beta = 2 pi b 2^-64 < 2 pi 2^-16 = 9.6e-5. K = fl(2 pi 2^-64) (relative error <= u); (double) b is exact (b < 2^53), so
  beta_c = fl(b K) has |beta_c - beta| <= beta (2u + u^2) <= 2.2e-20.
  cb = fl(1 - fl(0.5 fl(beta_c^2))): |cb - (1 - beta_c^2/2)| <= u + u beta_c^2/2; Taylor: 0 <= cos x - (1 - x^2/2) <= x^4/24 <= (2 pi 2^-16)^4/24 = 3.5204e-18 <= 3.54e-18;
  so |cb - cos beta| <= 1.1103e-16 + 3.54e-18 + 2.2e-20 = 1.1459e-16 <= 1.15e-16 (the later budgets 1.8e-16, 4.1e-16 <= 1e-15 are unchanged).
  sb = fl(beta_c fl(1 - fl(t fl(1/6)))), t = fl(beta_c^2): |sb - beta_c(1 - beta_c^2/6)| <= 2.3 u beta_c <= 2.5e-20; Taylor remainder
  <= x^5/120 <= 7e-23; so |sb - sin beta| <= 2.5e-20 + 7e-23 + 2.2e-20 <= 5e-20.
- Combination. cs = fl(fl(C_A cb) - fl(S_A sb)), and cos(theta_used) = cos A' cos beta - sin A' sin beta (A' = 2 pi A/2^16):
  |C_A cb - cos A' cos beta| <= 6e-17 (1 + 1e-15) + 1.15e-16 <= 1.8e-16; |S_A sb - sin A' sin beta| <= 6e-17 * 1e-4 + 5e-20 <= 1e-19;
  three roundings <= u (1 + 1e-15) + u 1e-4 + u (1 + 1e-4) <= 2.23e-16. Total |cs - cos theta_used| <= 4.1e-16 <= 1e-15.
- Sum. s = sum_n fl-products W_n cs_n, summed left to right: |s - sum W_n cs_n| <= gamma_N sum |W_n cs_n| with gamma_N = N u/(1 - N u)
  (N atoms; each product adds its rounding inside gamma_N's count), and |cs_n| <= 1 + 1e-15. W_n = fl(mid w_n), and Arb gives
  dW >= sum |W_n - w_n|.
- Hence (an absolute bound for the prime sum only), with sigma_W = sum |W_n|:
  |s - sum_n w_n cos(r_k log n)| <= sigma_W (1+u)(gamma_N (1+u) + u + 1e-15) + dW + sigma_W 2 pi k e_max 2^-64 =: E_s(k).
- Final. mu_c = fl(fl(H - T) - 2 s) (2 s exact): rounding <= 3u (|H| + |T| + 2 sigma_W (1.0001)).
  So mu(k) >= H - T - 2 sum w_n cos(...) >= mu_c(k) - E with E = 2 E_s(kmax) + 3u(max|H| + max|T| + 2.0002 sigma_W), multiplied by 1.01
  for the roundings in evaluating E itself. The first inequality is one-sided: mu(k) - (H - T - 2 sum w_n cos) >= 0 can be as large as the
  variation of h and tail over the block. So only mu(k) >= mu_c(k) - E is proved, and only that is used.
- Re-check after the correction above. The corrected Taylor constant leaves every stated budget unchanged, and the code's E uses 1e-15 per
  cosine, so E and the node sets are unchanged.
Observed: E = 2.1e-11 (chi_-3, c = 30), 1.1e-10 (chi_-3, c = 35), 2.6e-9 (chi_-4, c = 47, kmax = 8.1e8); the flag margin is 1e-6
and the code asserts E < 1e-7. The rigorous scans select exactly the node sets of the earlier float scans (224/63, 721/31, 794/253). The dominant term is the phase error
2 pi k e_n 2^-64 (k <= kmax).

## Assumptions
IEEE-754 binary64 arithmetic with round-to-nearest for +, -, * (the C compiler emits no fused or reassociated operations:
-O2 -ffp-contract=off, no -ffast-math); Arb (python-flint 0.9.0) encloses correctly; mpmath only chooses B_n (its error is bounded in Arb).

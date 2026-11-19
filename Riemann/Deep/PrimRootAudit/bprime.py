#!/usr/bin/env python3
"""bprime.py (ppraudit): Theorem B' for p > 10^7, omega <= 8: own least kappa at L_m = log 10^7 (own sieve, own kernel),
256-bit margin at its 3-digit ceiling, c = (kappa/L_m^2)^(1/3) vs the claimed 0.814."""
import sys
sys.argv = [sys.argv[0]]; sys.dont_write_bytecode = True
import indep_kappa as IK
from fractions import Fraction as Fr
M = IK.M
Lm = M.log(10 ** 7); cmax = (0, 0)
for w in range(1, 9):
    y = IK.tab[w][0]; Tn, C, _ = IK.rosser(w, y); T = M.mpf(Tn) / IK.ONE
    M.mp.prec = 64; lo_, hi_ = M.mpf(1), M.mpf(10) ** 6
    for _ in range(50):
        mid = M.sqrt(lo_ * hi_); mm = IK.margin(mid, Lm, T, C)
        if mm is not None and mm > 0: hi_ = mid
        else: lo_ = mid
    M.mp.prec = 256
    e = int(M.floor(M.log10(hi_))) - 2; k = M.ceil(hi_ / M.mpf(10) ** e) * M.mpf(10) ** e
    mg = IK.margin(k, Lm, T, C); c = M.cbrt(k / Lm ** 2)
    if c > cmax[0]: cmax = (c, w)
    print("w = %d: kappa(L >= log 1e7) = %s, margin %s, c = %s" % (w, M.nstr(k, 4), M.nstr(mg, 5), M.nstr(c, 6)))
print("max c over w <= 8 at L_m = log 1e7: %s (w = %d); claimed 0.814" % (M.nstr(cmax[0], 6), cmax[1]))

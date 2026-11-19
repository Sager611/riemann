#!/usr/bin/env python3
"""sizes.py (ppraudit): recompute the page's Sizes / DERIVATION sec. 7 table from logs/kappa_table.txt."""
import math
IB = 1 / (1 - 0.8733819666132076)
kap, Lw = {}, {}
for line in open("../PrimRootNum/logs/kappa_table.txt"):
    if not line.startswith("#"):
        f = line.split(); kap[int(f[0])] = float(f[6]); Lw[int(f[0])] = float(f[5])
for N in (10, 20, 100, 1000, 10000):
    L = N * math.log(10)
    wmax = max(w for w in kap if Lw[w] <= L)
    A = max(IB * math.log10(kap[w] * L) for w in kap if w <= wmax)
    print("p ~ 10^%d: max omega %d; worst A 10^%.1f; B' 10^%.1f; D 10^%.1f; w=2 10^%.1f; w=5 10^%.1f" % (
        N, wmax, A, 3 * IB * math.log10(0.814 * L), 3.4 * IB * math.log10(1.956 * L), IB * math.log10(kap[2] * L), IB * math.log10(kap[5] * L)))

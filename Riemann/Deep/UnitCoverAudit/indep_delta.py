#!/usr/bin/env python3
"""indep_delta.py (unitaudit) -- independent recomputation of UnitCover Lemma D (DERIVATION sec 2, 2.3) with mpmath
(prec from env INDEP_PREC, default 256), using the ExplProg auditor's own quadrature implementation
(Riemann/Deep/ExplProgAudit/indep.py: kernel constants, envelopes ks(Y), zero sums zco/Zof, majorant zmaj, Sum 1/gamma^2),
NOT cert.py.  Modes: time | cells (all 664 cells, worst) | select (extremes + random cells + tail + table rows + interior
and exponent-7.8978 checks)."""
import sys, os, time, random, importlib.util
sys.dont_write_bytecode = True          # no __pycache__ in UnitCover/ or ExplProgAudit/
t0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ind", os.path.join(HERE, "..", "ExplProgAudit", "indep.py"))
ind = importlib.util.module_from_spec(spec); spec.loader.exec_module(ind)
mp = ind.mp
B, b, LOG2, LOG9 = ind.B, ind.b, ind.LOG2, ind.LOG9
RS = mp.mpf(101624) / 100000
def delta(cp, L1, Nm1, L0=None, X=None, Lz=None):
    """delta bound; y = (cp L0)^(1/b) unless X (= log y) given; k's at Y = X - 7; Z at conductor e^L1; N - 1 = Nm1"""
    L0 = L1 if L0 is None else L0
    X = mp.log(cp * L0) / b if X is None else X
    Y = X - ind.ELL
    assert Y >= 10
    (n0, k1, k2), (m0, k1L, k2L) = ind.ks(Y)
    Z = ind.Zof(ind.zco(n0, k1, k2), n0, L1); ZL = ind.Zof(ind.zco(m0, k1L, k2L), m0, L1)
    z1 = ind.SIG * k2 / 2; z1L = ind.SIG * k2L / 2
    y = mp.exp(X); yB = mp.exp(B * X); yh = mp.exp(X / 2)
    triv = ind.W0 / (1 - mp.exp(-2 * Y))
    PP = ind.WINF * RS * (mp.exp(X / 2) + mp.exp(X / 3) * X / LOG2)
    num = (yB * Z + yh * ZL) / 2 + triv + PP + 2 * ind.WINF * (Nm1 * X + L1)
    den = y * ind.W1 - (yB * z1 + yh * z1L) - triv - PP
    assert den > 0
    return num / den
def tail(cp, d0, Lm, C1=4, C2=200):
    X = mp.log(cp * Lm) / b; Y = X - ind.ELL
    (n0, k1, k2), (m0, k1L, k2L) = ind.ks(Y)
    co = ind.zco(n0, k1, k2); coL = ind.zco(m0, k1L, k2L)
    yB = mp.exp(B * X); xhB = mp.exp((mp.mpf(1) / 2 - B) * X); x3B = mp.exp((mp.mpf(1) / 3 - B) * X)
    z1 = ind.SIG * k2 / 2; z1L = ind.SIG * k2L / 2
    lhs = ind.zmaj(co, n0, Lm) / 2 + xhB * ind.zmaj(coL, m0, Lm) / 2
    lhs += (1 + d0) * (ind.W0 / ((1 - mp.exp(-2 * Y)) * yB) + ind.WINF * RS * (xhB + x3B * X / LOG2)) / Lm
    lhs += 2 * ind.WINF * ((C1 + C2 / Lm) * X + 1) / yB + d0 * (z1 + xhB * z1L) / Lm
    rhs = d0 * cp * ind.W1
    return (rhs - lhs) / rhs
def cells():
    out = []; L0 = LOG9
    while L0 < 60:
        L1 = L0 + max(mp.mpf(1) / 100, L0 / 200); out.append((L0, L1)); L0 = L1
    return out
mode = sys.argv[1] if len(sys.argv) > 1 else "time"
CS = cells()
print("[%s] prec %d; %d cells, last L1 = %s; W1 = %s, W0 = %s, winf = %s, SIG = %s" % (mode, mp.mp.prec, len(CS), mp.nstr(CS[-1][1], 8),
      mp.nstr(ind.W1, 10), mp.nstr(ind.W0, 10), mp.nstr(ind.WINF, 9), mp.nstr(ind.SIG, 9)), flush=True)
if mode == "time":
    t1 = time.time(); d = delta(820, CS[0][1], 4 * CS[0][1] + 200, CS[0][0]); print("first cell delta %s in %.2f s" % (mp.nstr(d, 8), time.time() - t1))
if mode == "cells":
    sel = [(820, mp.mpf(1) / 1000), (820000, mp.mpf(1) / 10**6)]
    if len(sys.argv) > 2: sel = [x for x in sel if x[0] == int(sys.argv[2])]
    i0, i1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (0, len(CS))
    print("cells %d..%d" % (i0, i1 - 1))
    for (cp, d0) in sel:
        worst = (0, None)
        for (L0, L1) in CS[i0:i1]:
            d = delta(cp, L1, 4 * L1 + 200, L0)
            if d > worst[0]: worst = (d, L0)
        print("c'=%d: worst cell delta %s at L0 = %s (claim 0.00098366 / 9.8258e-7); <= d0: %s  [%.0f s]" % (cp, mp.nstr(worst[0], 8), mp.nstr(worst[1], 8), worst[0] <= d0, time.time() - t0), flush=True)
if mode == "select":
    rng = random.Random(7)
    pick = [0, 1, 2, len(CS) - 1] + rng.sample(range(3, len(CS) - 1), 6)
    for (cp, d0) in [(820, mp.mpf(1) / 1000), (820000, mp.mpf(1) / 10**6)]:
        for i in pick:
            L0, L1 = CS[i]; d = delta(cp, L1, 4 * L1 + 200, L0)
            # interior: true L in the cell, y = (cp L)^E with k's at true Y, Z at L, N-1 = 4L+200; and exponent 7.8978
            Lm = (L0 + L1) / 2; di = delta(cp, Lm, 4 * Lm + 200)
            dx = delta(cp, Lm, 4 * Lm + 200, X=mp.mpf('7.8978') * mp.log(cp * Lm))
            print("c'=%d cell %d L0=%s: cell bound %s; interior %s; exp 7.8978 %s; ok %s" % (cp, i, mp.nstr(L0, 6), mp.nstr(d, 8), mp.nstr(di, 8), mp.nstr(dx, 8), d <= d0 and di <= d and dx <= di), flush=True)
        m = tail(cp, d0, CS[-1][1]); print("c'=%d tail from L = %s: margin %s (claim 0.3443 / 0.3452)" % (cp, mp.nstr(CS[-1][1], 7), mp.nstr(m, 6)), flush=True)
        for L in [mp.mpf(61), mp.mpf(100), mp.mpf(1000), mp.mpf(10**4)]:
            print("   direct at L = %s: delta %s" % (mp.nstr(L, 6), mp.nstr(delta(cp, L, 4 * L + 200), 6)), flush=True)
    # table rows: (log10 Q, c', N max of the row, d0)
    for (e, cp, N, d0) in [(3, 750, 33, 1e-3), (6, 670, 52, 1e-3), (12, 600, 97, 1e-3), (30, 530, 231, 1e-3), (100, 490, 745, 1e-3),
                           (1000, 460, 7336, 1e-3), (12, 600000, 82, 1e-6), (100, 490000, 497, 1e-6), (1000, 460000, 4787, 1e-6)]:
        L = e * mp.log(10); d = delta(cp, L, N - 1)
        lgy = mp.log(cp * L) / b / mp.log(10)
        print("table Q=10^%d c'=%d N=%d: delta %s <= %g: %s; log10 y = %s" % (e, cp, N, mp.nstr(d, 6), d0, d <= d0, mp.nstr(lgy, 6)), flush=True)
print("[%s] %.1f s" % (mode, time.time() - t0))

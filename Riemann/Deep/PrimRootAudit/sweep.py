#!/usr/bin/env python3
"""sweep.py (ppraudit): every row w = 1..10^4 of PrimRootNum/logs/kappa_table.txt: consecutive coverage, own L_w = log(P_w + 1),
own margin (4.4) at the table kappa (own kernel, Lemma Z quadrature, sup scan; T_w from the table lowered by 1e-8 relative,
C_w from the table), c_w^3 L_w^2 >= kappa_w, Theorem D need c = 1.956, Corollary E (table kappa) per w > 120."""
import sys, math, time
sys.argv_ = list(sys.argv); sys.argv = [sys.argv[0]]
import indep_kappa as IK
M = IK.M; M.mp.prec = 53
def supZL(co, N0, Lm):
    best = M.mpf(0); L = Lm
    for k in range(90):             # L = Lm e^{k/3}, k < 90 (up to Lm e^30); the limit za is added
        best = max(best, IK.Z(co, N0, max(L, IK.L9)) / L); L = L * M.e ** (M.mpf(1) / 3)
    return max(best, co[0])
IK.supZL = supZL
def clear_quad_caches():           # mpmath caches transformed nodes per interval: memory grows with distinct intervals
    for v in list(vars(M.mp).values()):
        tc = getattr(v, "transformed_cache", None)
        if isinstance(tc, dict):
            for key in [k for k in tc if not (k[0] == -1 and k[1] == 1)]: del tc[key]   # keep standard [-1, 1] nodes
W0, W1_ = (int(sys.argv_[1]), int(sys.argv_[2])) if len(sys.argv_) > 2 else (1, 10000)
STEP = int(sys.argv_[3]) if len(sys.argv_) > 3 else 1          # sample every STEP-th row (and the last 20 rows)
t0 = time.time(); rows = sorted(IK.tab); assert rows == list(range(1, 10001)), "table rows not 1..10000"
Pw = 1; worst = None; bad = []; badc = []; badD = []; badE = []; cmax = (0, 0); short = 0.0
for q in IK.PR[:W0 - 1]: Pw *= q
for w in range(W0, W1_ + 1):
    y, T, C, Lmt, kap, cw = IK.tab[w]
    Pw *= IK.PR[w - 1]
    if (w - W0) % STEP and w <= W1_ - 20: continue
    Lm = M.log(Pw + 1) if Pw + 1 >= 3 else M.log(3)
    if abs(Lm - M.mpf(Lmt)) > 2e-6: bad.append(("Lm", w))
    k = M.mpf(kap.numerator) / kap.denominator
    mg = IK.margin(k, Lm, M.mpf(T) * (1 - M.mpf(10) ** -8), C)
    if mg is None or mg <= 0: bad.append(("margin", w, mg))
    rel = mg / (M.mpf(T) * IK.W1 * k)
    if worst is None or rel < worst[0]: worst = (rel, w)
    if M.mpf(cw) ** 3 * Lm ** 2 < k:
        badc.append(w); short = max(short, float((k - M.mpf(cw) ** 3 * Lm ** 2) / k))
    clear_quad_caches()
    if cw > cmax[0]: cmax = (cw, w)
    if M.mpf("1.956") ** M.mpf(3.4) * Lm ** M.mpf(2.4) < k: badD.append(w)
    if w > 120:
        okp = M.log(k * Lm) * IK.IB < Lm and M.log(k * Lm) * IK.IB + M.mpf(1) / 1000 < Lm / 2
        if not okp: badE.append(w)
    if w % (250 if STEP == 1 else 1000) == 0: print("  w = %d done (%.0f s): min rel margin so far %s at w = %d; failures %d/%d/%d/%d" % (w, time.time() - t0, M.nstr(worst[0], 5), worst[1], len(bad), len(badc), len(badD), len(badE)), flush=True)
print("chunk w = %d..%d: printed-c_w shortfall rows %d, max relative shortfall (kappa - c_w^3 L^2)/kappa = %.3g" % (W0, W1_, len(badc), short))
print("rows 1..10000 consecutive; Lm/margin failures: %s" % (bad[:10] if bad else "none"))
print("min relative margin (4.4)/(T W1 kappa) = %s at w = %d" % (M.nstr(worst[0], 6), worst[1]))
print("c_w^3 L_w^2 >= kappa failures: %s; max c_w = %s at w = %d; Theorem D (1.956) failures: %s; Corollary E (w > 120) failures: %s" % (badc[:10] or "none", cmax[0], cmax[1], badD[:10] or "none", badE[:10] or "none"))
print("c_w <= 0.474 for w > 100: %s; c_10000 = %s" % (all(IK.tab[w][5] <= 0.474 for w in rows if w > 100), IK.tab[10000][5]))
print("RSS-limited chunk done")
print("done %.0f s" % (time.time() - t0))

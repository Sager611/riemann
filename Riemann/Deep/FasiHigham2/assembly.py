
# assembly.py : tau-sweep verification of  c_n > 0  (for nu >= nu0, n >= 8) on tau in [tau_a, tau_b]
import json, math, sys, time, pickle, os
import fh2lib as L
from fh2lib import arb, acb, ab, ball, polar, up, S, Psie, PI, lg_box
def lo(x): return arb(x.lower())
def hi(x): return arb(x.upper())
YJ = arb('0.85'); RMAX = 1.6
RG = {0: 0.0, 1: 0.03, 2: 0.05, 3: 0.07, 4: 0.08}; RG5 = 0.116
def r_g(l): return RG.get(l, RG5)
# ---------- LG grid ----------
RB = [0.0] + [0.004 * (1.04 ** k) for k in range(0, 160)]
RB = [x for x in RB if x < 1.75] + [1.75]
NA = 16
AB = [float(PI / 2) * j / NA for j in range(NA + 1)]
LGCACHE = {}
def lg_grid(nu0, path):
    global LGCACHE
    if os.path.exists(path):
        LGCACHE = pickle.load(open(path, 'rb')); return
    t0 = time.time()
    for i in range(len(RB) - 1):
        for j in range(NA):
            r0, r1 = RB[i], RB[i + 1]
            if r1 >= 1.0 and AB[j + 1] > 1.45: continue     # near/beyond turning point: not needed
            res = lg_box(arb(r0), arb(r1), arb(AB[j]), arb(AB[j + 1]) if j + 1 < NA else PI / 2, nu0)
            if res is not None:
                LGCACHE[(i, j)] = (float(res[2].upper()) * (1 + 1e-9), float(res[3].upper()) * (1 + 1e-9))
    pickle.dump(LGCACHE, open(path, 'wb'))
    print('lg grid built', len(LGCACHE), '%.0fs' % (time.time() - t0), flush=True)
def lam_region(r0, r1, a0, a1):
    # max (lamI, lamK) over grid cells covering the polar region
    # conservative cover: cells [RB[i],RB[i+1]] x [AB[j],AB[j+1]] covering [r0-1e-12, r1+1e-12] x [a0-1e-12, a1+1e-12]
    a0 = min(max(a0, 0.0), AB[NA]); a1 = min(max(a1, 0.0), AB[NA])
    if a0 > a1: a0, a1 = a1, a0
    if r1 + 1e-12 > RB[-1]:
        if r1 > RB[-1] + 1e-9: return None
        r1 = RB[-1] - 1e-12
    r0 = max(r0 - 1e-12, 0.0); r1 = r1 + 1e-12; a0 = max(a0 - 1e-12, 0.0); a1 = min(a1 + 1e-12, AB[NA])
    i0 = max(i for i in range(len(RB) - 1) if RB[i] <= r0)
    i1 = min(i for i in range(len(RB) - 1) if RB[i + 1] >= r1)
    j0 = max(j for j in range(NA) if AB[j] <= a0)
    j1 = min(j for j in range(NA) if AB[j + 1] >= a1)
    lI = lK = 0.0
    for i in range(i0, i1 + 1):
        for j in range(j0, j1 + 1):
            v = LGCACHE.get((i, j))
            if v is None: return None
            lI = max(lI, v[0]); lK = max(lK, v[1])
    return arb(lI), arb(lK)
# ---------- saddle / rates ----------
def saddle(a, t):
    t = ab(t)
    if not (t > 0): return arb(0)
    if a == 1: return t * (t + 2) / (2 * (t + 1))
    return (-(t + a) + a * ((t + a) ** 2 - a * a + 1).sqrt()) / (a * a - 1)
def Phi(a, t):
    t = ab(t)
    if not (t > 0): return Psie(acb(0)).real * a     # A_a(0) = a Psie(0)
    r = saddle(a, t)
    return r + a * Psie(acb(r)).real - t * r.log()
def Aa(a, z): return z + a * Psie(z)
TAB = None
def good_factor(l, rlo_, rhi_, nu0):
    """log of (1-e)(1+r^2)^{1/4} e^{-a lw0} e^{-l/(12 nu0)} / sqrt(2 pi nu0 a bh) ; None if not good"""
    a = 2 * l + 1
    if rhi_ > RMAX or rlo_ < r_g(l): return None
    if l == 0 and rhi_ <= 0.15:
        rows = TAB['0small']
    else:
        rows = TAB[str(l)] if l <= 4 else TAB['>=5']
    if l == 0 and rlo_ < 0.15 < rhi_:
        rows = TAB['0small'] + TAB['0']
    sel = [x for x in rows if x[1] >= rlo_ and x[0] <= rhi_]
    # the selected boxes must cover [rlo_, rhi_]
    sel.sort(key=lambda x: x[0])
    if not sel or sel[0][0] > rlo_ or sel[-1][1] < rhi_: return None
    for u, v in zip(sel, sel[1:]):
        if v[0] > u[1]: return None
    if not sel or any(x[2] is None for x in sel): return None
    e = max(x[2] for x in sel)
    if e >= 0.9: return None
    # beta_hi / nu : general rows store per-unit-a bound (beta <= nu a bh), small rows store beta/nu bound
    bnu = max(x[3] * (a if x[5] != 'small' else 1) for x in sel)
    bh = arb(bnu)
    lam = lam_region(0.8 * rlo_, min(1.25 * rhi_ + 1e-6, 1.74), 0, 0.001)
    if lam is None: return None
    lw0 = l * lam[0] + (l + 1) * lam[1]
    val = (1 - arb(e)).log() - lw0 - arb(l) / (12 * nu0) - (2 * PI * nu0 * bh).log() / 2
    return lo(val)
# ---------- crude bounds: return list of (rate_per_nu (upper), c (upper), len) ----------
def crude_pieces(l, tlo, thi, rc, nu0):
    a = 2 * l + 1; out = []
    rc = arb(rc)
    if rc <= YJ * arb('0.999'):
        phic = PI / 2
    else:
        phic = (YJ / rc).asin()
    lam = lam_region(float(rc) * 0.999, float(rc) * 1.001, 0, float(phic) if float(phic) < 1.5707 else 1.5708)
    if lam is None: return None
    lW = l * lam[0] + (l + 1) * lam[1]
    Arc = Aa(a, acb(rc)).real
    lr = rc.log()
    tt = ab(tlo) if lr > 0 else ab(thi)
    out.append((hi(Arc - tt * lr), hi((1 + rc * rc).log() / 4 + lW), hi(phic / PI)))
    if rc <= YJ * arb('0.999'):
        lam2 = lam_region(float(rc), float(YJ), 1.5707, 1.5708)
        if lam2 is None: return None
        lW2 = l * lam2[0] + (l + 1) * lam2[1]
        R = a * Psie(acb(0, rc)).real - ab(thi) * lr
        out.append((hi(R), hi(lW2 + (YJ / rc).log().log()), hi(1 / PI)))
    else:
        xc = (rc * rc - YJ * YJ).sqrt()
        Kx = 30
        for k in range(Kx):
            x0 = xc * k / Kx; x1 = xc * (k + 1) / Kx
            z = acb(ball(x0, x1), YJ)
            lz = abs(z).log()
            R1 = Aa(a, z).real - ab(tlo) * lz; R2 = Aa(a, z).real - ab(thi) * lz
            R = hi(R1) if hi(R1) > hi(R2) else hi(R2)
            angl = float(lo((YJ / x1).atan())); angh = float(hi((YJ / x0).atan())) if x0 > 0 else 1.5708
            lam3 = lam_region(float(lo(abs(z))), float(hi(abs(z))), angl, min(angh, 1.5708))
            if lam3 is None: return None
            lW3 = l * lam3[0] + (l + 1) * lam3[1]
            c = hi((abs(1 + z * z)).log() / 4 - abs(z).log()) + lW3
            out.append((R, c, hi((x1 - x0) / PI)))
    return out
def neg_pieces(LN, tau_b, nu0):
    # sum over l >= LN of |J_l| (N <= -1), imaginary-axis contour
    lam = lam_region(0.0, float(YJ), 1.5707, 1.5708)
    if lam is None: return None
    lI, lK = lam
    kap = (2 * LN - ab(tau_b))
    if kap < 0: return None
    big = bool(nu0 * kap >= 1)       # |N| = nu(2LN - tau) >= 1 always; use kap only if nu0*kap >= 1
    # X(y) <= X(yJ) = exp(2 nu Re eta(i yJ) + lI + lK)
    reta = Psie(acb(0, YJ)).real + YJ.log()
    lnX = 2 * nu0 * reta + lI + lK
    if not (lnX < 0): return None
    c0 = LN * (lI + lK) + lK - (1 - lnX.exp()).log()
    out = []
    ys = [arb(0)] + [YJ * k / 40 for k in range(1, 41)]
    for k in range(len(ys) - 1):
        y0, y1 = ys[k], ys[k + 1]
        if k == 0:
            # int_0^{y1} y^{|N|-1} dy = y1^{|N|}/|N| <= y1^{nu kap} (big) or <= y1 ; Re Psie(iy) <= Psie(0)
            if big:
                R = (2 * LN + 1) * Psie(acb(0)).real + kap * y1.log(); c = c0
            else:
                R = (2 * LN + 1) * Psie(acb(0)).real; c = c0 + y1.log()
            out.append((hi(R), hi(c), arb(1) / PI))
        else:
            yb = ball(y0, y1)
            if big:
                R = (2 * LN + 1) * Psie(acb(0, yb)).real + kap * yb.log(); c = c0 - y0.log()
            else:
                R = (2 * LN + 1) * Psie(acb(0, yb)).real; c = c0
            out.append((hi(R), hi(c), hi((y1 - y0) / PI)))
    return out
def IL_pieces(tau_b, nu0):
    out = []
    K = 24
    for k in range(K):
        f0 = PI / 2 + PI / 2 * k / K; f1 = PI / 2 + PI / 2 * (k + 1) / K
        z = polar(YJ, YJ, f0, f1)
        R = (z + Psie(z)).real - ab(tau_b) * YJ.log()
        lam = lam_region(float(YJ) * 0.999, float(YJ) * 1.001, float(PI - f1), float(PI - f0))
        if lam is None: return None
        eK = lam[1].exp() - 1      # |h_K| <= e^{lamK}-1 ... (lamK = -log(1-eps) => eps = 1-e^{-lamK})
        epsK = 1 - (-lam[1]).exp()
        c = hi(abs(1 + z * z).log() / 4 - (1 - epsK).log())
        out.append((hi(R), c, hi((f1 - f0) / PI)))   # (1/2pi)*2 (symmetric halves)
    return out
def check_box(tau_a, tau_b, nu0, verbose=False):
    tau_a, tau_b = float(tau_a), float(tau_b)
    Lmax = int(math.floor(tau_b / 2 + 1e-12))
    best = None; bads = []
    for l in range(0, Lmax + 1):
        a = 2 * l + 1
        tlo_a = arb(tau_a) - 2 * l; thi_a = arb(tau_b) - 2 * l      # exact (arb)
        if thi_a < 0: continue
        if not (tlo_a > 0): tlo_a = arb(0)
        tlo = tlo_a; thi = thi_a
        rlo_ = max(float(lo(saddle(a, tlo_a))) - 1e-12, 0.0) if tlo_a > 0 else 0.0
        rhi_ = float(hi(saddle(a, thi_a))) + 1e-12
        gf = None
        if rlo_ >= r_g(l) and rhi_ <= RMAX and not (l >= 1 and not (tlo > 0)):
            gf = good_factor(l, rlo_, rhi_, nu0)
        if gf is not None:
            ph = min(lo(Phi(a, tlo)), lo(Phi(a, thi)))
            val = nu0 * ph + gf
            if best is None or val > best[0]: best = (val, l, ph, gf)
        else:
            if rhi_ > RMAX: rc = RMAX
            else:
                tm = (tlo + thi) / 2
                rc = max(float(saddle(a, tm).mid()), 0.01) if tm > 0 else 0.01
            pcs = crude_pieces(l, tlo, thi, rc, nu0)
            if pcs is None: return dict(ok=False, why='crude %d' % l)
            bads += [(p, 'l%d' % l) for p in pcs]
    LN = int(math.floor(tau_a / 2 + 1e-12)) + 1
    npcs = neg_pieces(LN, tau_b, nu0)
    if npcs is None: return dict(ok=False, why='neg')
    bads += [(p, 'neg') for p in npcs]
    ipcs = IL_pieces(tau_b, nu0)
    if ipcs is None: return dict(ok=False, why='IL')
    bads += [(p, 'IL') for p in ipcs]
    if best is None: return dict(ok=False, why='no good')
    val, lb, ph, gf = best
    # total bad / best  at nu0, with monotonicity requirement nu0*delta >= 1
    tot = arb(0); mind = None; worst = None
    for ((R, c, ln), tag) in bads:
        d = ph - R
        if mind is None or d < mind: mind = d; worst = tag
        tot += ln * (nu0 * (R - ph) + c).exp()
    ratio = tot / gf.exp()
    ok = bool(ratio < arb('0.5')) and bool(nu0 * mind > 1)
    return dict(ok=ok, ratio=float(ratio.upper()) if ratio.is_finite() else float('inf'), lbest=lb, mind=float(mind.lower()), worst=worst)
def load(nu0, tabfile, lgfile):
    global TAB
    TAB = json.load(open(tabfile))
    lg_grid(nu0, lgfile)
if __name__ == '__main__':
    nu0 = arb(sys.argv[1]); load(nu0, sys.argv[2], sys.argv[3])
    ta, tb, w = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
    t0 = time.time(); tau = ta; nfail = 0; worst_ratio = 0; mind_all = 1e9
    edges = []
    x = ta
    while x < tb - 1e-12:
        nx = min(tb, x + w)
        ev = math.floor(x / 2 + 1e-12) * 2 + 2
        if x < ev < nx: nx = ev
        edges.append((x, nx)); x = nx
    for (x, nx) in edges:
        res = check_box(x, nx, nu0)
        if not res['ok']:
            nfail += 1
            if nfail <= 20: print('FAIL box', round(x, 5), round(nx, 5), res, flush=True)
        else:
            worst_ratio = max(worst_ratio, res['ratio']); mind_all = min(mind_all, res['mind'])
    print('tau in [%g,%g]: boxes %d fails %d worst ratio %.3e min delta %.4f  %.0fs' % (ta, tb, len(edges), nfail, worst_ratio, mind_all, time.time() - t0), flush=True)

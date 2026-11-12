#!/usr/bin/env python3
"""uaudit N3: Step 5c on actual zeros (and with planted off-line pairs).
S(t) = pi^-1 sum_{gamma in [T/2,5T/2]} log|(t+h/2-gamma)/(t-h/2-gamma)| = H[F~](t) exactly (on-line zeros),
F~(u) = #{gamma in [T/2,5T/2] : |u-gamma| < h/2}.
(1) Re V_exact(t) - S(t) at random t in [T,2T] via mpmath log|xi|: the far-zero + gamma part, claimed O(1).
(2) int_T^{2T} |S|^m / T, int |Im V|^m / T, ||F~||_{L^m(R)}, and the M. Riesz ratio ||S||_{L^m[T,2T]}/||F~||_{L^m(R)}
    against Pichorides' norm cot(pi/(2m)).
(3) planted: 5% of the zeros replaced by off-line pairs at depth d: S_pl = H[F~_pl] with Poisson-smeared F~_pl.
"""
import sys, json, time, math
import numpy as np
from scipy.integrate import quad

ZF = sys.argv[2] if len(sys.argv) > 2 else "../ExploreHelio/results/zeta_zeros.npy"
Z = np.sort(np.load(ZF))
T = 2000.0
L = math.log(T)
PHI = 0.7
H = 2 * math.pi * PHI / L
ZW = Z[(Z >= T / 2) & (Z <= 2.5 * T)]
rng = np.random.default_rng(7)

def S_on(t, g):
    return float(np.sum(np.log(np.abs((t + H / 2 - g) / (t - H / 2 - g))))) / math.pi

def S_pl(t, g_on, g_off, d):
    v = np.sum(np.log(np.abs((t + H / 2 - g_on) / (t - H / 2 - g_on))))
    v += np.sum(np.log(((t + H / 2 - g_off) ** 2 + d * d) / ((t - H / 2 - g_off) ** 2 + d * d)))  # pair: 2 * (1/2) log
    return float(v) / math.pi

def moment_S(fun, sing, a, b, m):
    pts = np.sort(np.concatenate([[a, b], sing[(sing > a) & (sing < b)]]))
    tot = 0.0
    for i in range(len(pts) - 1):
        lo, hi = pts[i], pts[i + 1]
        if hi - lo < 1e-14:
            continue
        tot += quad(lambda t: abs(fun(t)) ** m, lo, hi, limit=50, epsrel=1e-5, epsabs=1e-10)[0]
    return tot

def count_moment(events_plus, events_minus, a, b, m, weight=None):
    """int_a^b c(t)^m, c = #{gamma : |t-gamma|<h/2}: c jumps +1 at gamma-h/2, -1 at gamma+h/2."""
    ev = np.concatenate([np.stack([events_plus, np.ones_like(events_plus)], 1), np.stack([events_minus, -np.ones_like(events_minus)], 1)])
    ev = ev[np.argsort(ev[:, 0], kind="stable")]
    c, tot, prev = 0, 0.0, a
    # initial count at a
    c = int(np.sum((events_plus <= a) & (events_minus > a)))
    for (x, s) in ev:
        if x <= a:
            continue
        if x >= b:
            break
        tot += (x - prev) * c ** m
        c += int(s)
        prev = x
    tot += (b - prev) * c ** m
    return tot

def smeared_moment(g_on, g_off, d, a, b, m):
    def F(u):
        v = np.sum(np.abs(u - g_on) < H / 2)
        v += (2 / math.pi) * np.sum(np.arctan((u + H / 2 - g_off) / d) - np.arctan((u - H / 2 - g_off) / d))
        return float(v)
    sing = np.concatenate([g_on - H / 2, g_on + H / 2])
    return moment_S(F, sing, a, b, m)

if __name__ == "__main__":
    t0 = time.time()
    out = dict(T=T, L=L, h=H, nW=int(len(ZW)))
    # (1) exact check
    import mpmath as mp
    mp.mp.dps = 15
    def logxi_abs(tau):
        s = mp.mpc(0.5, tau)
        return float(mp.log(abs(s * (s - 1) / 2)) - 0.25 * mp.log(mp.pi) + mp.re(mp.loggamma(s / 2)) + mp.log(abs(mp.zeta(s))))
    diffs = []
    for t in rng.uniform(T, 2 * T, 120):
        ex = (logxi_abs(t + H / 2) - logxi_abs(t - H / 2)) / math.pi
        diffs.append(ex - S_on(t, ZW))
    diffs = np.array(diffs)
    out["far_part"] = dict(mean=float(diffs.mean()), std=float(diffs.std()), maxabs=float(np.abs(diffs).max()))
    print("Re V_exact - S: mean %.4f std %.4f max|.| %.4f  (%.0fs)" % (diffs.mean(), diffs.std(), np.abs(diffs).max(), time.time() - t0), flush=True)
    # (2) moments on actual zeros
    sing = np.concatenate([ZW - H / 2, ZW + H / 2])
    rows = []
    for m in (2, 4, 6):
        Sm = moment_S(lambda t: S_on(t, ZW), sing, T, 2 * T, m)
        Im = count_moment(ZW - H / 2, ZW + H / 2, T, 2 * T, m)
        Fm = count_moment(ZW - H / 2, ZW + H / 2, T / 2 - 1, 2.5 * T + 1, m)
        ratio = (Sm ** (1 / m)) / (Fm ** (1 / m))
        row = dict(m=m, ReV_m_over_T=Sm / T, ImV_m_over_T=Im / T, Ft_norm_m_over_T=Fm / T, riesz_ratio=ratio, pichorides=1 / math.tan(math.pi / (2 * m)))
        rows.append(row)
        print("m=%d  int|S|^m/T=%.4f  int|ImV|^m/T=%.4f  ||F~||_m^m/T=%.4f  ratio=%.3f <= cot(pi/2m)=%.3f  (%.0fs)" % (
            m, Sm / T, Im / T, Fm / T, ratio, row["pichorides"], time.time() - t0), flush=True)
    out["actual"] = rows
    # (3) planted off-line pairs
    pl = []
    for d in (1 / L, 3 / L, 0.45):
        mask = rng.random(len(ZW)) < 0.05
        g_on, g_off = ZW[~mask], ZW[mask]
        sing = np.concatenate([g_on - H / 2, g_on + H / 2])
        for m in (2, 4):
            Sm = moment_S(lambda t: S_pl(t, g_on, g_off, d), sing, T, 2 * T, m)
            Fm = smeared_moment(g_on, g_off, d, T / 2 - 1, 2.5 * T + 1, m)
            # sharp count with pairs twice
            ep = np.concatenate([g_on - H / 2, g_off - H / 2, g_off - H / 2]); em = np.concatenate([g_on + H / 2, g_off + H / 2, g_off + H / 2])
            Im = count_moment(np.sort(ep), np.sort(em), T, 2 * T, m)
            row = dict(d=d, m=m, ReV_m_over_T=Sm / T, smeared_m_over_T=Fm / T, ImV_sharp_m_over_T=Im / T,
                       riesz_ratio=(Sm ** (1 / m)) / (Fm ** (1 / m)))
            pl.append(row)
            print("planted d=%.3f m=%d  int|ReV|^m/T=%.4f  ||F~_smeared||^m/T=%.4f  int|ImV|^m/T=%.4f  ratio=%.3f  (%.0fs)" % (
                d, m, Sm / T, Fm / T, Im / T, row["riesz_ratio"], time.time() - t0), flush=True)
    out["planted"] = pl
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "vmoments.json", "w"), indent=1)
    print("time %.1f s" % (time.time() - t0))

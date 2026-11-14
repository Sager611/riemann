# Lever equivalents (BeyondB2 sigma, lambda), full-system entry check, and gain-multiplier table.
import numpy as np, sys
from mpmath import mp, mpf, sqrt, findroot
from system import Rd_native, R_line, G0, G1, WIN, B2, D2, R2
from infimum import Bstar, golden, mk_line, DOM
from ceiling import curve, mk_curve, best
mp.dps = 30
def Bsig(s):
    r = sqrt(964 * s**2 - 1156 * s + 1081)
    return (103528 * s**2 + 7316 * s * r + 51710 * s - 3823 * r - 51257) / (12 * (8720 * s**2 + 640 * s * r + 4912 * s - 356 * r - 5179))
def sig_for(B): return findroot(lambda s: Bsig(s) - B, mpf('0.0003'))
# BeyondB2 fmodel entries with a general count function Rf (sv = 0)
P1 = 1124661 / 250000; LST = 0.21
def entries(Bv, bv, Rf, n=241):
    dd, xx = np.meshgrid(np.linspace(11 / 50, 13 / 20, n), np.linspace(6 / 25, 0.5, n))
    l = 11 / 3 - 4 * Bv; ly = (1 + bv - l) / 2; h = (1 + bv + 3 * l) / 2; T = 2 * ly - h
    C0 = -0.25 + bv / 6 + 5 * l / 4; th0 = dd * xx / P1
    Rdd = Rf(dd, xx)
    E = lambda R: C0 + dd / 2 + dd * xx * (l + (h / 2 - ly) / P1) - 3 * h * (1 - R) / 4
    S = lambda R: (3 * R + 1) / 4 + th0 / 2 - 17 / 50
    N = lambda R: (3 * R + 1) / 4 - th0 / 2
    Dsh = lambda R: R + (0.5 - th0) * 3 / 40 - N(R)
    JN = 5 / 2 - 5 / P1
    GN = dd * xx * JN / 10 * T; GL = th0 * (0.70 * h - ly); GR = (dd / 2 - th0) * (ly - 0.53 * h)
    ent = {"G_N-E (active)": GN - E(Rdd), "G_L-E": GL - E(Rdd), "G_R-E": GR - E(Rdd), "S(R_d)-1/5": S(Rdd) - 0.2,
           "1/250-E": 1 / 250 - E(Rdd), "(p1-2dx)L-R_d": (P1 - 2 * th0 * P1) * LST - Rdd, "-E-hDsh": -E(Rdd) - h * Dsh(Rdd)}
    return {k: float(np.min(v)) for k, v in ent.items()}
if __name__ == '__main__':
    thinf = 0.8733466805; th = 2183367 / 2500000
    for nm, Bv in [("PR137 theta", th), ("PR137 line infimum", thinf), ("ceiling (continuum)", 0.8733458652)]:
        s = sig_for(mpf(Bv))
        print("%-22s B=%.10f  B_2-B (=lambda) = %.4e  sigma_eq = %.6e  (BeyondB2: 0.873 needs sigma 0.003518, lambda 0.000382)" % (nm, Bv, float(B2) - Bv, float(s)))
    # full system entries at PR137 point and at the ceiling point
    Rl = mk_line(float(G0), float(G1))
    bh = 2 * 0.81447 - 1 - 3 * (11 / 3 - 4 * th)
    print("entries at PR137 (theta, h=0.81447):", {k: "%.2e" % v for k, v in entries(th, bh, Rl).items()})
    Gf, lo, hi = curve(LEFT=0.70, STEP=1 / 400000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0)
    Rc = mk_curve(Gf, lo, hi, (0, 1, 0, 1))
    B, bo, h, d, x = best(Rc)
    print("entries at ceiling (B=%.10f, b=%.5f):" % (B, bo), {k: "%.2e" % v for k, v in entries(B + 1e-9, bo, Rc).items()})
    sys.stdout.flush()
    # gain multiplier: count paid by mu * z_inf(r)
    for mu in (2, 5, 10, 12):
        Rm = mk_curve(lambda w, mu=mu: mu * Gf(w), lo, hi, (0, 1, 0, 1))
        B, bo, h, d, x = best(Rm, 0.06, 0.20)
        e = entries(B + 1e-9, bo, Rm)
        print("mu=%-3d B*=%.7f b=%.4f contact (%.4f,%.3f) min other entry %.2e (%s)" % (mu, B, bo, d, x, min(v for k, v in e.items() if k != "G_N-E (active)"), min((v, k) for k, v in e.items() if k != "G_N-E (active)")[1]))
        sys.stdout.flush()

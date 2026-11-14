# Ceiling: the scalar system with the count paid by the exact gain curve of the anchor chain (instead of the line).
import numpy as np, sys
from system import Rd_native, R_line, G0, G1, WIN, Broot, B2
from infimum import Bstar, golden, mk_line, DOM
from chain_limit import run
def curve(**kw):
    out, _ = run(**kw)
    rs = np.array([o[0] for o in out][::-1] + [out[0][1]])           # increasing r
    zs = np.array([o[4] for o in out][::-1] + [out[0][2]])
    return (lambda w: np.interp(w, rs, zs, left=np.nan, right=0.0)), rs[0], rs[-1]
def mk_curve(Gf, lo, hi, win=WIN):
    def Rf(dd, xx):
        V = 5 / 6 - dd
        f = lambda c: (1 - dd * (c + 2 * xx * ((1 - c) / 2 + Gf(c)))) - (1 - dd + V * (1.5 * c - 1))
        a = np.full_like(dd, lo); b = np.full_like(dd, hi); fa = f(a); fb = f(b); ok = (fa > 0) != (fb > 0)
        for _ in range(60):
            m = (a + b) / 2; fm = f(m); s = (fm > 0) == (fa > 0)
            a = np.where(s, m, a); fa = np.where(s, fm, fa); b = np.where(s, b, m)
        c = (a + b) / 2; v = 1 - dd + V * (1.5 * c - 1)
        inwin = (dd >= win[0]) & (dd <= win[1]) & (xx >= win[2]) & (xx <= win[3]) & np.isfinite(v) & ok
        return np.where(inwin, np.minimum(v, Rd_native(dd, xx)), Rd_native(dd, xx))
    return Rf
def best(Rf, blo=0.08, bhi=0.14):
    bo = golden(lambda b: Bstar(b, Rf, 801, 41)[0], blo, bhi, 1e-8)
    B, d, x = Bstar(bo, Rf, 2401, 81)
    l = 11 / 3 - 4 * B
    return B, bo, (1 + bo + 3 * l) / 2, d, x
if __name__ == '__main__':
    rows = []
    for nm, kw in [("PR137 chain curve (pads as PR137)", dict(LEFT=0.7335)),
                   ("continuum, zero pads (STEP 1/400000)", dict(LEFT=0.70, STEP=1 / 400000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0)),
                   ("continuum, zero pads, TOP=7/9", dict(TOP=7 / 9 - 1e-12, LEFT=0.70, STEP=1 / 400000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0))]:
        Gf, lo, hi = curve(**kw)
        for wn, win in [("PR137 window", WIN), ("window = all", (0, 1, 0, 1))]:
            B, bo, h, d, x = best(mk_curve(Gf, lo, hi, win))
            print("%-40s %-13s B*=%.10f  b=%.5f h=%.5f contact (%.4f, %.3f)  G(0.7524)=%.7f  curve left end %.4f" % (nm, wn, B, bo, h, d, x, float(Gf(0.7524)), lo))
            sys.stdout.flush()
    B, bo, h, d, x = best(mk_line(float(G0), float(G1)))
    print("PR137 line (reference)                   B*=%.10f  b=%.5f h=%.5f" % (B, bo, h))

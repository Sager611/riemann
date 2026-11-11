
# T2 analysis: structured pseudospectrum of L(s, chi_13) inside its FE class.
# F_c = cos c * Lambda(s, chi_12) + sin c * G(s), G = DH-type block from the order-3 pair (labels 3, 9).
# On the line: Z_c(t) = cos c * Zq(t) + sin c * ZG(t). Zeros of Z_c  <=>  r(t) = Zq/ZG = -tan c.
# Off-line pairs are born (EP) at local minima of r with r* < 0 when tan c passes |r*|.
import sys, json, numpy as np
d = np.load(sys.argv[1]); ts = d["ts"]; Z = d["Z"]
Zq, ZG, ZH = Z[0], Z[1] + Z[2], Z[3] + Z[4]
from math import lgamma, log, pi
import cmath
def theta(t):  # arg gamma(1/2+it), gamma(s) = (13/pi)^{s/2} Gamma(s/2), continuous
    import mpmath
    return float(t/2*log(13/pi) + mpmath.im(mpmath.loggamma(0.25 + 0.5j*t)))
def nsc(z):   # cumulative sign-change count
    s = np.signbit(z).astype(np.int8)
    return np.concatenate([[0], np.cumsum(np.abs(np.diff(s)))])
T = ts[-1]
print("grid: %d points, t in [%.3f, %.3f], h = %.4f" % (len(ts), ts[0], T, ts[1]-ts[0]))
th = theta(T)/pi
for name, z in (("chi_12 (quad)", Zq), ("block A (3,9)", ZG), ("block B (4,10)", ZH)):
    print("%-16s sign changes to T: %d ; theta(T)/pi = %.3f" % (name, nsc(z)[-1], th))
# smooth comparison at many heights (Turing-like mean of S)
idx = np.searchsorted(ts, np.arange(50, T, 50))
thv = np.array([theta(t)/pi for t in ts[idx]])
for name, z in (("quad", Zq), ("A", ZG), ("B", ZH)):
    D = nsc(z)[idx] - thv
    print("D(t) = N_signchanges - theta/pi at t=50..: %s: mean %.2f, values %s" % (name, D.mean(), np.round(D,1).tolist()[:16]))
# local extrema of r = Zq/ZG between consecutive zeros of ZG (parabolic refinement)
r = Zq / ZG
crit = []
for i in range(1, len(ts)-1):
    if ZG[i-1]*ZG[i] <= 0 or ZG[i]*ZG[i+1] <= 0: continue
    a, b, c_ = r[i-1], r[i], r[i+1]
    if (b < a and b < c_) or (b > a and b > c_):
        den = a - 2*b + c_
        x = 0.5*(a - c_)/den if den != 0 else 0.0
        tv = ts[i] + x*(ts[1]-ts[0]); rv = b - 0.25*(a - c_)*x
        crit.append((float(tv), float(rv), "min" if b < a else "max"))
negmin = [(t, -v) for (t, v, k) in crit if k == "min" and v < 0]
negmax = [(t, -v) for (t, v, k) in crit if k == "max" and v < 0]
print("critical points of r: %d; negative minima (EP births for c>0): %d; negative maxima (re-entries): %d" % (len(crit), len(negmin), len(negmax)))
# record curve: lowest birth height for each threshold
negmin.sort()
rec = []; best = 1e9
for t, a in negmin:
    if a < best:
        best = a; rec.append((t, a, float(np.arctan(a))))
print("record EP births (height t*, |r*| = tan c*, c*):")
for t, a, c in rec: print("  t* = %8.3f  tan c* = %.5f  c* = %.5f" % (t, a, c))
json.dump({"negmin": negmin, "negmax": negmax, "rec": rec}, open("t2_crit.json", "w"))
# per c: sign-change deficit relative to c = 0 and wrong-sign extrema of Z_c
N0 = nsc(Zq)
print("c, tan c, deficit (N_line(0,T)-N_line(c,T))/2 at T=200,400,600,%d, #wrong-sign extrema of Z_c, first wrong-sign extremum" % int(T))
for c in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, np.pi/2]:
    zc = np.cos(c)*Zq + np.sin(c)*ZG
    Nc = nsc(zc)
    out = []
    for TT in (200, 400, 600, T):
        k = np.searchsorted(ts, TT) - 1
        out.append((N0[k] - Nc[k])/2)
    ws = []
    for i in range(1, len(ts)-1):
        if (zc[i] > zc[i-1] and zc[i] > zc[i+1] and zc[i] < 0) or (zc[i] < zc[i-1] and zc[i] < zc[i+1] and zc[i] > 0):
            ws.append(float(ts[i]))
    print("c=%.3f tan=%.4f deficits %s  wrong-sign extrema %d  first %s" % (c, np.tan(c), out, len(ws), ws[:4]))
    json.dump(ws, open("t2_ws_%.3f.json" % c, "w"))

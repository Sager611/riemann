#!/usr/bin/env python3
# dbnhead prep.py: coefficient enclosures + exact uint64 phase data of row i's finite head (all 4095 labels n <= 4096) for bb.c
# head(t,x) = sum_n c_n(t,x) cos(tau log n),  c_n = (1/2) d_t(n) n^{-1/2-tL/4} F_i(n),  tau = x/2 + pi t/8,  L = log(x/4pi)  (U.tex l.1044-1077, 1595-1604)
# d_t(n) >= 0 nondecreasing in t (Lemma src:coeff) -> c_n in sign(F) [0.5 d_t0 n^{-1/2-t1 Lhi/4}|F|, 0.5 d_t1 n^{-1/2-t0 Llo/4}|F|] over the row's time cell and x in [X, Xhi]
# phase: tau log n = (x'/2) log n with x' = x + pi t/4 in [X + pi t0/4, Xhi + pi t1/4]
# usage: prep.py ARC row u gfrac out.dat [u0] [tight]   (window x in [X e^u0, X e^u], u0 = 0 default; thr = (1 - gfrac) * mx_pp, mx_pp = frozen prime-power torus max as in ceilaudit5)
import sys, math, json
import numpy as np, mpmath as mp
sys.path.insert(0, 'Riemann/Deep/DBNHigh')
import model as M
mp.mp.dps = 50
arc, ri, u, gfrac, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
K = 4096
row = M.load_rows(arc)[ri]; t0, t1 = row['t0'], row['t1']
u0 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
Xm = mp.mpf(6000000185827); Xhim = Xm*mp.exp(mp.mpf(u)); Xlom = Xm*mp.exp(mp.mpf(u0))
Llo = float(mp.log(Xlom/(4*mp.pi))); Lhi = float(mp.log(Xhim/(4*mp.pi)))
d0, ln = M.dt_coeffs(t0, N=K); d1, _ = M.dt_coeffs(t1, N=K)
n = np.arange(2, K + 1); lnn = ln[n]
Fv = M.v_row(row, lnn)
hi_mag = 0.5*d1[n]*np.exp(-(0.5 + t0*Llo/4)*lnn)*np.abs(Fv)*(1 + 1e-11) + 1e-15
lo_mag = np.maximum(0.0, 0.5*d0[n]*np.exp(-(0.5 + t1*Lhi/4)*lnn)*np.abs(Fv)*(1 - 1e-11) - 1e-15)
TIGHT = int(sys.argv[7]) if len(sys.argv) > 7 else 0
if TIGHT:   # primes: d_t(p) = log p b_t(p) exactly, so c_p = 0.5 log p p^-1/2 exp(t log p (log p - L)/4)|F|, decreasing in t and L (log p < L)
    isp = np.array([all(int(m) % q for q in range(2, int(m**0.5) + 1)) for m in n])
    e_hi = np.exp(t0*lnn*(lnn - Llo)/4); e_lo = np.exp(t1*lnn*(lnn - Lhi)/4)
    hi_mag = np.where(isp, 0.5*lnn*np.exp(-0.5*lnn)*e_hi*np.abs(Fv)*(1 + 1e-11) + 1e-15, hi_mag)
    lo_mag = np.where(isp, np.maximum(0.0, 0.5*lnn*np.exp(-0.5*lnn)*e_lo*np.abs(Fv)*(1 - 1e-11) - 1e-15), lo_mag)
sg = np.sign(Fv)
lo = np.where(sg > 0, lo_mag, -hi_mag); hi = np.where(sg > 0, hi_mag, -lo_mag)
c_fr = 0.5*d1[n]*np.exp(-(0.5 + t0*M.c0)*lnn)*Fv          # ceilaudit5 frozen (d_{t+}, c0)
c_act = 0.5*d0[n]*np.exp(-(0.5 + t0*Llo/4)*lnn)*Fv        # actual at t = t0, x = X (witness only)
def factor(m):
    f = {}; p = 2
    while p*p <= m:
        while m % p == 0: f[p] = f.get(p, 0) + 1; m //= p
        p += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f
groups = {}   # base -> list of (k, lo, hi, cfr, cact)
for i, m in enumerate(n):
    f = factor(int(m))
    if len(f) == 1:
        (p, e), = f.items(); groups.setdefault(p, []).append((e, lo[i], hi[i], c_fr[i], c_act[i]))
    else:
        groups[int(m)] = [(1, lo[i], hi[i], c_fr[i], c_act[i])]
th = np.linspace(0, 2*math.pi, 400001)
def gmax(hs, idx):
    f = sum(h[idx]*np.cos(h[0]*th) for h in hs); return float(f.max())
mx_pp = sum(gmax(hs, 3) for b, hs in groups.items() if len(factor(b)) == 1)
comp = float(sum(abs(hs[0][3]) for b, hs in groups.items() if len(factor(b)) > 1))
tri = float(np.abs(c_fr).sum())
thr = (1 - gfrac)*mx_pp
D = 4*math.pi/math.log(2)          # double; anchors x_k = xs + k D exactly
M1, M2, M3 = 32, 8, 8
xs = Xlom + mp.pi*t0/4; xe = Xhim + mp.pi*t1/4
Kan = int(mp.ceil((xe - xs)/D)) + 1
Dm = mp.mpf(D); h1 = Dm/M1; h2 = h1/M2; h3 = h2/M3
two64 = mp.mpf(2)**64
def u64(v): return int(mp.floor(mp.frac(v)*two64)) % (1 << 64)
def u64r(v): return int(mp.nint(mp.frac(v)*two64)) % (1 << 64)
order = sorted(groups, key=lambda b: -sum(max(abs(h[1]), abs(h[2])) for h in groups[b]))
order.remove(2); order = [2] + order
with open(out, 'w') as fo:
    fo.write('%d %d %.17g %d %d %d %.17g %.17g %.17g %.17g\n' % (len(order), Kan, D, M1, M2, M3, thr, mx_pp, tri, comp))
    for b in order:
        w = mp.log(b)/(4*mp.pi)        # turns per unit x
        hs = groups[b]
        fo.write('%d %d %d %d %d %d %d %.17g\n' % (b, len(hs), u64(xs*w), u64r(Dm*w), u64r(h1*w), u64r(h2*w), u64r(h3*w), math.log(b)/2))
        for h in hs: fo.write('%d %.17g %.17g %.17g\n' % (h[0], h[1], h[2], h[4]))
print(json.dumps(dict(row=ri, t0=t0, t1=t1, u0=u0, u=u, Xhi=float(Xhim), K=Kan, groups=len(order), mx_pp=mx_pp, comp=comp, tri=tri, thr=thr, gfrac=gfrac,
                      hi_sum=float(np.abs(hi).sum()), lo_rel_min=float(np.min(lo_mag/np.maximum(hi_mag, 1e-300))))))

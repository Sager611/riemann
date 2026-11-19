"""padhi2 follow-up, Thm 5.2: independent float check of Padhi's margin M2 (l.577) at box corners.
W(z) = int_0^D e^{-zu} w(u) du by mpmath quadrature of BTY Definition 1 (independent of his closed form in rep_cert.W);
kappa_m and c_k (hence a_0, a_1, a) parsed from the BTY TeX; constants E2, E3, c10, G from his bty2 (as his certify uses them,
eta0 = ceil(1e4 e_max)/1e4, sigma0 = floor(1e5 sigma_lo)/1e5). Compares with his ball margin (rep_cert.margin) at the same point
and finds, per corner, the distance d_fail where M2 first drops to 0 (falsification target: d_fail < d* would refute the table).
usage: m2corner.py lk,k [lk,k ...]"""
import os, sys, json, math, re
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
SCR = os.path.join(ROOT, '.scratch', 'padhi2')
REP = os.path.join(SCR, 'repo', 'ramanujans-inequality-main', 'code', 'rep'); os.chdir(REP); sys.path.insert(0, REP)
import mpmath as mp
mp.mp.dps = 20
tex = open(os.path.join(SCR, 'bty', 'Bellotti_Trudgian_Yang_v1.tex')).read()
# --- BTY constants from the TeX
th = mp.mpf('1.1338'); tt = mp.tan(th); ct = 1/tt; s2 = 1/mp.cos(th)**2; D = 2*th*ct
def w(u): return s2*(s2*(th*ct - u/2)*mp.cos(u*tt) + 2*th*ct - u + mp.sin(2*th - u*tt)/mp.sin(2*th) - 2*(1 + mp.sin(th - u*tt)/mp.sin(th)))
w0 = w(0)
seg = tex[tex.index('P(x) := \\frac{1}{14912370}'):][:3000]          # the c_k block of BTY Section 4 (TeX l.757-780)
ck = [int(v) for i, v in sorted(((int(i), v) for i, v in re.findall(r'c_\{?(\d+)\}?\s*=\s*(-?\d+)', seg)), key=lambda p: p[0])]
assert len(ck) == 17, ck
NORM = sum(c*c for c in ck)
ak = [mp.mpf(sum(ck[j]*ck[j+k] for j in range(17-k))*(1 if k == 0 else 2))/14912370 for k in range(17)]
a0, a1, a = ak[0], ak[1], sum(ak[1:])
print('BTY TeX: c_k =', ck, ' sum c^2 =', NORM, ' a0 = %s a1 = %s a = %s' % (mp.nstr(a0, 12), mp.nstr(a1, 12), mp.nstr(a, 12)))
print('   a1 vs 865534/497079: %s ; a vs 2919857/828465: %s ; w0 = %s' % (mp.nstr(a1 - mp.mpf(865534)/497079, 3), mp.nstr(a - mp.mpf(2919857)/828465, 3), mp.nstr(w0, 12)))
import rep_cert as RC
RC.lamt = {float(k): v for k, v in json.load(open('lamstar2c.json')).items()}
km = [mp.mpf(1), mp.mpf(-851)/859, mp.mpf(780)/859, mp.mpf(-525)/859, mp.mpf(171)/859, mp.mpf(28)/859, mp.mpf(-29)/859]
print('   kappa_m (his) =', [mp.nstr(x, 8) for x in km], ' sum =', mp.nstr(sum(km), 10), '= 433/859:', abs(sum(km) - mp.mpf(433)/859) < 1e-15)
import numpy as np
_x, _wt = np.polynomial.legendre.leggauss(400)                      # Gauss-Legendre on [0, D], float64 (independent of his closed form)
_Df = float(D); _u = (_x + 1)*_Df/2; _wt = _wt*_Df/2
_wu = np.array([float(w(mp.mpf(float(x)))) for x in _u])
def W(z):
    z = complex(z); return complex(np.sum(_wt*_wu*np.exp(-z*_u)))
# sanity: quadrature W vs mpmath quad at two points
for zz in (mp.mpc(0.5, -3), mp.mpc(70, -8)):
    q = mp.quad(lambda u: mp.exp(-zz*u)*w(u), [D*i/24 for i in range(25)])
    print('   W check at', complex(zz), ': GL400', W(zz), ' mp.quad', complex(q), flush=True)
KM = [float(x) for x in km]
def F(z, eta, sig):
    z, eta, sig = complex(z), float(eta), float(sig)
    return sum(k*W((z + (2*sig - 1)*m)/eta) for m, k in enumerate(KM))
import bty2
def consts(eta0c, sig0):
    E2, E3 = bty2.errconsts(eta0c, sig0); c10 = bty2.lemma10(eta0c, sig0)[0]; G = bty2.gamma_const(sig0)
    return [float(x) for x in (E2, E3, c10, G)]
A0f, A1f, Af, W0f = float(a0), float(a1), float(a), float(w0)
def M2(Lt, sig, e0, e1, d, C):
    E2, E3, c10, G = C; Lt, sig, e0, e1, d = float(Lt), float(sig), float(e0), float(e1), float(d)
    lhs = e0*W0f*((Af*(433/859)/2)*Lt - G) - c10*e0*W0f + A0f*F(sig - 1, e0, sig).real + 1e-10*e0 + 2e-7 + E2*(e0*Lt)*e0**2 + E3*e0**3
    r1 = A1f*(F(sig - 1 + e0, e0, sig) + F(sig - e0, e0, sig)).real
    r2 = A1f*(F(complex(sig - 1 + e1, -d), e0, sig) + F(complex(sig - e1, -d), e0, sig)).real
    return r1 + r2 - lhs, r1 - lhs, r2
DA = json.load(open('dstar2c_all.json'))

# --- table consistency: dstar2c_all(lk, k) <= min over l <= lk of the certified value inherited from the next certified k above
DC = json.load(open('dcert_29_30_31_32_33_34_35_36_37.json'))
def inherited(l, k):
    ks_ = sorted(int(c.split(',')[1]) for c in DC if int(c.split(',')[0]) == l)
    if not ks_ or k < ks_[0]: return None            # k below the zero-free threshold of l: no zero, imposes nothing
    above = [c for c in ks_ if c >= k]
    return DC['%d,%d' % (l, min(above))] if above else 0.0
bad = []; n = 0
for key, v in DA.items():
    lk, k = [int(x) for x in key.split(',')]
    vals = [inherited(l, k) for l in range(29, lk + 1)]; vals = [x for x in vals if x is not None]
    need = min(vals) if vals else 0.0; n += 1
    if v > need + 1e-12: bad.append((key, v, need))
    # zero-free threshold of l: certified ks start at kmin(l) = floor(lambda*(l)/l*1e4) + 1
for l in range(29, 38):
    ks_ = sorted(int(c.split(',')[1]) for c in DC if int(c.split(',')[0]) == l)
    kmin = math.floor(RC.lam_at(l)/l*1e4) + 1
    if ks_[0] != kmin: bad.append(('kmin', l, ks_[0], kmin))
print('table check: %d entries; entries above min_l(certified, inherited) or kmin mismatches: %s' % (n, bad[:8] if bad else 'none'))
from flint import arb, acb
for a_ in ([] if '--nocorner' in sys.argv else [x for x in sys.argv[1:] if ',' in x]):
    lk, k = [int(x) for x in a_.split(',')]; dst = DA['%d,%d' % (lk, k)]
    emin = RC.lam_at(lk)/lk; emax = k*1e-4
    A = RC.lam_at(math.log(16*math.exp(lk) + 1e10)); sig = 1 - A/math.log(16*math.exp(lk) + 1e10)
    sig_lo = 1 - A/math.log(16*math.exp(max(lk - 1, 28.7)) + 1e10)
    C = consts(math.ceil(emax*1e4)/1e4, math.floor(sig_lo*1e5)/1e5)
    print('cell (%d,%d): d* = %.6f, e_min = %.6f, e_max = %.4f, A = %.6f, sigma(top) = %.7f, consts E2 %.2f E3 %.1f c10 %.5f G %.5f'
          % (lk, k, dst, emin, emax, A, sig, *[float(x) for x in C]), flush=True)
    # his ball margin at the top corner for comparison
    Cb = RC.consts(math.ceil(emax*1e4)/1e4, math.floor(sig_lo*1e5)/1e5)
    for (e0, e1) in [(emax, emax), (emin, emax), (emax, emin), (emin, emin)]:
        m, m1, r2 = M2(mp.mpf(lk), mp.mpf(sig), mp.mpf(e0), mp.mpf(e1), mp.mpf(dst), C)
        hb = RC.margin(arb(lk), arb(sig), arb(e0), arb(e1), arb(dst), Cb)
        # d_fail: smallest d in (0, 1] with M2 <= 0 for this ordering (scan then bisect); None if M2 > 0 up to d = 1
        prev, dfail = None, None
        for dd in [dst*j/8 for j in range(1, 65)]:
            if M2(mp.mpf(lk), mp.mpf(sig), mp.mpf(e0), mp.mpf(e1), mp.mpf(dd), C)[0] <= 0: dfail = dd; break
            prev = dd
        if dfail is not None:
            lo_, hi_ = prev or 0.0, dfail
            for _ in range(25):
                mid = (lo_ + hi_)/2
                if M2(mp.mpf(lk), mp.mpf(sig), mp.mpf(e0), mp.mpf(e1), mp.mpf(mid), C)[0] <= 0: hi_ = mid
                else: lo_ = mid
            dfail = hi_
        print('  log t = %d, (eta0, eta1) = (%.6f, %.6f): M2(d*) = %.6e [M1 = %.6e, pair-2 term %.6e]; his ball margin %s; d_fail = %s'
              % (lk, e0, e1, float(m), float(m1), float(r2), hb.str(5, radius=False), ('%.6f (= %.3f d*)' % (dfail, dfail/dst)) if dfail else '> 8 d*'), flush=True)

# --- interior grid (both orderings), independent float: min over the grid of max(M2(order 1), M2(order 2))
if '--grid' in sys.argv:
    for a_ in [x for x in sys.argv[1:] if ',' in x]:
        lk, k = [int(x) for x in a_.split(',')]; dst = DA['%d,%d' % (lk, k)]
        emin = RC.lam_at(lk)/lk; emax = k*1e-4; A = RC.lam_at(math.log(16*math.exp(lk) + 1e10))
        l0 = max(lk - 1, 28.7); sig_lo = 1 - A/math.log(16*math.exp(l0) + 1e10)
        C = consts(math.ceil(emax*1e4)/1e4, math.floor(sig_lo*1e5)/1e5)
        best = None; n = 0
        for Lt in np.linspace(l0, lk, 5):
            sig = 1 - A/math.log(16*math.exp(Lt) + 1e10)
            for e0 in np.linspace(emin, emax, 7):
                for e1 in np.linspace(emin, emax, 7):
                    for dd in np.linspace(0, dst, 11):
                        m = max(M2(Lt, sig, e0, e1, dd, C)[0], M2(Lt, sig, e1, e0, dd, C)[0]); n += 1
                        if best is None or m < best[0]: best = (m, Lt, e0, e1, dd)
        print('grid (%d,%d): %d points, min over grid of max(M2 both orderings) = %.4e at log t %.3f, eta0 %.6f, eta1 %.6f, d %.5f'
              % (lk, k, n, *best), flush=True)

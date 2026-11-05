# C5: Zenodo 23272725 - independent check of the 44 new kernel triples (r, rho, lambda) of Prop. "New certificate family":
#   2 K_y(d,b) >= lambda (y/h) J_{h,r}(d,b),  0<y<=rho h, 0<=b<=y, d real   (scale h=1)
# (a) adversarial minimisation of the relative defect in the ORIGINAL variables (not the author's polynomial);
# (b) the polynomial reduction: J_{1,r} = 2G/D4 and 2K_y = 4y(D+v-B)/Q at random points, D4>0, so P>=0 <=> inequality;
# (c) self caps s_r(w) <= M on 0<=w<=rho^2; (d) the 7/8 seed and startup inequality.
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 kernels116_check.py GEOMETRY_CERTIFICATES.json HISTORY.json.gz
import signal; signal.alarm(1200)  # hard time limit
import json, sys, math, gzip
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
mp.mp.dps = 40
recs = json.load(open(sys.argv[1]))['records']
rng = np.random.default_rng(1)
def K(h, d, b):
    return (h - b)/(d*d + (h - b)**2) + (h + b)/(d*d + (h + b)**2)
def J(r, d, b):
    return 0.5*(K(1.0, d + r, b) + K(1.0, d - r, b))
def reldef(lam, r, d, y, b):
    k2 = 2*K(y, d, b)
    return (k2 - lam*y*J(r, d, b))/k2
worst = []
ident_err = 0.0
for rec in recs:
    r = float(F(rec['r'])); rho = float(F(rec['rho'])); lam = float(F(rec['lambda']))
    # samples
    n = 60000
    y = rho*rng.random(n)**0.5; y = np.maximum(y, 1e-9)
    b = y*rng.random(n); b[: n//6] = y[: n//6]*(1 - 1e-9*rng.random(n//6))
    d = np.exp(rng.uniform(math.log(1e-7), math.log(1e4), n)); d[: n//10] = 0.0
    v = reldef(lam, r, d, y, b)
    k = int(np.argmin(v))
    best = (float(v[k]), float(d[k]), float(y[k]), float(b[k]))
    # refine with Nelder-Mead in (log d, y/rho, b/y) on box
    def fobj(x):
        dd = math.exp(x[0]); yy = rho*min(max(x[1], 1e-9), 1.0); bb = yy*min(max(x[2], 0.0), 1.0 - 1e-12)
        return reldef(lam, r, dd, yy, bb)
    starts = [[math.log(max(best[1], 1e-7)), best[2]/rho, best[3]/best[2]]] + [[math.log(x), 1.0, s] for x in (1e-3, 0.3, 1.0, 3.0, 30.0) for s in (0.0, 0.5, 0.999)]
    fb = (1e9, None)
    for s0 in starts:
        res = minimize(fobj, s0, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-15, 'maxiter': 3000})
        if res.fun < fb[0]: fb = (res.fun, res.x)
    x = fb[1]; dd = math.exp(x[0]); yy = rho*min(max(x[1], 1e-9), 1.0); bb = yy*min(max(x[2], 0.0), 1.0 - 1e-12)
    # high precision recheck
    Km = lambda h, d_, b_: (h - b_)/(d_*d_ + (h - b_)**2) + (h + b_)/(d_*d_ + (h + b_)**2)
    D_, Y_, B_ = mp.mpf(dd), mp.mpf(yy), mp.mpf(bb)
    k2 = 2*Km(Y_, D_, B_); Jm = (Km(1, D_ + mp.mpf(F(rec['r']).numerator)/F(rec['r']).denominator, B_) + Km(1, D_ - mp.mpf(F(rec['r']).numerator)/F(rec['r']).denominator, B_))/2
    relmp = (k2 - mp.mpf(F(rec['lambda']).numerator)/F(rec['lambda']).denominator*Y_*Jm)/k2
    # far-field limit 1 - lambda/2
    worst.append(dict(r=rec['r'], rho=rec['rho'], lam=rec['lambda'], sample_min=best[0], refined_min=float(relmp), at=(dd, yy, bb), farfield=1 - lam/2))
    # (b) identity check of the polynomial reduction at random points
    for _ in range(200):
        dd2 = float(np.exp(rng.uniform(-5, 4))); yy2 = rho*float(rng.random()); bb2 = yy2*float(rng.random())
        D, vv, B, R = dd2**2, yy2**2, bb2**2, r*r
        C = D + R + 1 + B; A = C*C + 4*D*R - 4*B; D4 = A*A - 16*D*R*C*C; G = A*(C - 2*B) - 8*D*R*C
        Q = (D + vv - B)**2 + 4*D*B
        e1 = abs(J(r, dd2, bb2) - 2*G/D4)/abs(J(r, dd2, bb2)); e2 = abs(2*K(yy2, dd2, bb2) - 4*yy2*(D + vv - B)/Q)/abs(2*K(yy2, dd2, bb2))
        assert D4 > 0 and Q > 0
        ident_err = max(ident_err, e1, e2)
mins = sorted(worst, key=lambda w: w['refined_min'])
for w in mins[:6]: print(json.dumps(w))
print('MIN_REL_DEFECT_44', mins[0]['refined_min'], 'NEGATIVE', [ (w['r'], w['rho']) for w in worst if w['refined_min'] < 0 or w['sample_min'] < 0])
print('POLY_REDUCTION_MAX_REL_ERR', ident_err)
# (c) self caps
for r, M in ((F(23, 40), F(1503067131, 10**9)), (F(3, 5), F(25, 17))):
    rr = float(r); ws = np.linspace(0, 0.45**2, 200001)
    s = 2*(1 + rr*rr - ws)/((1 + rr*rr - ws)**2 + 4*rr*rr*ws)
    # exact check via the quadratic: M w^2 + (M(2r^2-2)+2) w + M(1+r^2)^2 - 2(1+r^2) >= 0 on [0, 0.2025]
    a2 = M; a1 = M*(2*r*r - 2) + 2; a0 = M*(1 + r*r)**2 - 2*(1 + r*r)
    wv = -a1/(2*a2); cand = [F(0), F(81, 400)] + ([wv] if 0 <= wv <= F(81, 400) else [])
    qmin = min(a2*w*w + a1*w + a0 for w in cand)
    print('SELF_CAP', str(r), 'max_s_grid', float(s.max()), 'M', float(M), 'exact_quadratic_min', float(qmin), 'OK' if qmin >= 0 and s.max() <= float(M) + 1e-15 else 'FAIL')
# (d) seed / startup from the main primary history
H = json.loads(gzip.decompress(open(sys.argv[2], 'rb').read()))
q0 = F(H['initial_q']); delay = F(H['initial_delay']); c0 = F(H['records'][0]['rate'])
print('SEED', str(q0), 'q0-9/16', float(q0 - F(9, 16)), 'startup', float(q0 - (2 + c0*q0)*delay - F(9, 16)), 'OK' if q0 - (2 + c0*q0)*delay > F(9, 16) else 'FAIL')

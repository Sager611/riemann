# Independent spot check of the 67 kernel inequalities of Zenodo 23752103 (Theorem src:kernel):
#   R_i(d, sqrt q) <= 4/(d^2+4q) for d real, 0<q<=cap_i;  O_i(q)=R_i(0,sqrt q) nondecreasing on [0,cap_i]
# Written from the TeX formulas only (not from the author's code).
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 kernel_check.py PLAN.json
import signal; signal.alarm(1500)  # hard time limit
import json, sys, math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
mp.mp.dps = 40
plan = json.load(open(sys.argv[1]))
rows = plan['rows']
rng = np.random.default_rng(12345)

def parse(r):
    nodes = [(float(F(n['height'])), float(F(n['shift'])), float(F(n['weight']))) for n in r['nodes']]
    nodesQ = [(F(n['height']), F(n['shift']), F(n['weight'])) for n in r['nodes']]
    jj = r.get('joint_jet_overlay')
    alpha = float(F(jj['alpha'])) if jj else 1.0
    cents = [(float(F(c['p'])), [float(F(b)) for b in c['b']]) for c in jj['centres']] if jj else []
    centsQ = [(F(c['p']), [F(b) for b in c['b']]) for c in jj['centres']] if jj else []
    alphaQ = F(jj['alpha']) if jj else F(1)
    cap = float(F(r['q_interval'][1])); lo = float(F(r['q_interval'][0]))
    return nodes, alpha, cents, cap, lo, nodesQ, alphaQ, centsQ, F(r['q_interval'][1])

def Kh(h, d, y):
    return (h - y)/(d*d + (h - y)**2) + (h + y)/(d*d + (h + y)**2)

def R(nodes, alpha, cents, d, y):
    s = 0.0
    for h, a, w in nodes:
        s += w*0.5*(Kh(h, d + a, y) + Kh(h, d - a, y))
    s *= alpha
    for p, b in cents:
        z1 = (p - y) - 1j*d; z2 = (p + y) - 1j*d
        acc = 0.0
        for r_, br in enumerate(b):
            acc += br*(p**r_)*((-1)**r_)*((z1**(-r_-1)) + (z2**(-r_-1))).real
        s += 2.0*acc/p
    return s

def Rmp(nodesQ, alphaQ, centsQ, d, y):
    d = mp.mpf(d); y = mp.mpf(y)
    def K(h, dd):
        h = mp.mpf(h.numerator)/h.denominator
        return (h - y)/(dd*dd + (h - y)**2) + (h + y)/(dd*dd + (h + y)**2)
    s = mp.mpf(0)
    for h, a, w in nodesQ:
        A = mp.mpf(a.numerator)/a.denominator; W = mp.mpf(w.numerator)/w.denominator
        s += W*(K(h, d + A) + K(h, d - A))/2
    s *= mp.mpf(alphaQ.numerator)/alphaQ.denominator
    for p, b in centsQ:
        P = mp.mpf(p.numerator)/p.denominator
        z1 = mp.mpc(P - y, -d); z2 = mp.mpc(P + y, -d)
        acc = mp.mpf(0)
        for r_, br in enumerate(b):
            B = mp.mpf(br.numerator)/br.denominator
            acc += B*P**r_*(-1)**r_*mp.re(z1**(-r_-1) + z2**(-r_-1))
        s += 2*acc/P
    return s

out = []
worst_all = (1e9, None)
for r in rows:
    i = r['row']
    nodes, alpha, cents, cap, lo, nodesQ, alphaQ, centsQ, capQ = parse(r)
    hmin = min([h for h, a, w in nodes] + [p for p, b in cents])
    assert hmin*hmin > cap, ('pole', i)
    # random samples: q in (0,cap], d >= 0 (R even in d)
    nq, nd = 160, 400
    qs = np.concatenate([cap*np.logspace(-8, 0, nq//2), cap*rng.random(nq//2)])
    ds = np.concatenate([[0.0], np.logspace(-6, 4, nd - 1)])
    best = (1e9, 0, 0)
    for q in qs:
        y = math.sqrt(q)
        Rv = R(nodes, alpha, cents, ds, y)
        tgt = 4.0/(ds*ds + 4*q)
        rel = 1.0 - Rv/tgt
        k = int(np.argmin(rel))
        if rel[k] < best[0]:
            best = (float(rel[k]), float(ds[k]), float(q))
    # local refinement from best sample: minimise relative defect over (log d, q)
    def f(v):
        d = math.exp(v[0]); q = min(max(v[1], 1e-14), cap); y = math.sqrt(q)
        return 1.0 - R(nodes, alpha, cents, d, y)*(d*d + 4*q)/4.0
    x0 = [math.log(max(best[1], 1e-6)), best[2]]
    res = minimize(f, x0, method='Nelder-Mead', options={'xatol':1e-10, 'fatol':1e-14, 'maxiter':4000})
    dmin = math.exp(res.x[0]); qmin = min(max(res.x[1], 1e-14), cap)
    relmp = 1 - Rmp(nodesQ, alphaQ, centsQ, dmin, mp.sqrt(qmin))*(mp.mpf(dmin)**2 + 4*mp.mpf(qmin))/4
    # far field coefficient lim d^2 R / 4
    dbig = 1e6
    ff = float(Rmp(nodesQ, alphaQ, centsQ, dbig, mp.sqrt(cap/2))*mp.mpf(dbig)**2/4)
    # own pair monotonicity on a grid
    qg = np.linspace(0.0, cap, 2001)
    Og = R(nodes, alpha, cents, np.zeros_like(qg), np.sqrt(qg))
    mono = float(np.min(np.diff(Og)))
    rec = dict(row=i, cap=cap, nodes=len(nodes), jets=len(cents), sample_min_rel=best[0], at_d=best[1], at_q=best[2],
               refined_min_rel=float(relmp), refined_d=dmin, refined_q=qmin, farfield_ratio=ff, own_min_step=mono,
               own_at_cap=float(Og[-1]))
    out.append(rec)
    if float(relmp) < worst_all[0]: worst_all = (float(relmp), i)
    print(json.dumps(rec), flush=True)
print('WORST_REL_DEFECT', worst_all, 'NEGATIVE_ROWS', [o['row'] for o in out if o['refined_min_rel'] < 0 or o['sample_min_rel'] < 0], 'NONMONOTONE_ROWS', [o['row'] for o in out if o['own_min_step'] < -1e-12])


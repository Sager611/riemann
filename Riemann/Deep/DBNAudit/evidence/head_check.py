# C3: independent recomputation of the finite arithmetic head of Zenodo 23752103 (App. B (ii), App. D):
#   H(t_ev) = 1/2 sum_{n=2}^{4096} d_{t+}(n) n^{-1/2 - t_ev c0} |F(n)|,  F(n) = sum_j w_j n^{-h_j/2} cos(a_j log n/2)
# and an adversarial search for the torus supremum of Re P_t(chi) = sum_n sgn(F(n)) a_n(t) Re chi_n
# (chi completely multiplicative, |chi_p|=1), compared with the charged head H - C (C = arithmetic credits).
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 head_check.py PLAN.json DIRECT.json rows...
import signal; signal.alarm(1500)  # hard time limit
import json, sys, math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
mp.mp.dps = 30
plan = {r['row']: r for r in json.load(open(sys.argv[1]))['rows']}
direct = {r['row']: r for r in json.load(open(sys.argv[2]))['rows']}
rows = [int(x) for x in sys.argv[3:]]
K = 4096; N0 = 690988
c0 = float(mp.log(N0)/2)
# factorizations
spf = list(range(K + 1))
for p in range(2, int(K**0.5) + 1):
    if spf[p] == p:
        for m in range(p*p, K + 1, p):
            if spf[m] == m: spf[m] = p
def fac(n):
    out = {}
    while n > 1:
        p = spf[n]; out[p] = out.get(p, 0) + 1; n //= p
    return out
primes = [p for p in range(2, K + 1) if spf[p] == p]
pidx = {p: i for i, p in enumerate(primes)}
E = np.zeros((K + 1, len(primes)))   # exponent matrix
for n in range(2, K + 1):
    for p, e in fac(n).items(): E[n, pidx[p]] = e
divs = [[] for _ in range(K + 1)]
for d in range(2, K + 1):
    for m in range(2*d, K + 1, d): divs[m].append(d)
logn = [mp.mpf(0)] + [mp.log(n) for n in range(1, K + 1)]
cache = {}
def dcoef(t):
    key = str(t)
    if key in cache: return cache[key]
    t = mp.mpf(t.numerator)/t.denominator
    b = [mp.mpf(0)] + [mp.exp(t*logn[n]**2/4) for n in range(1, K + 1)]
    d = [mp.mpf(0)]*(K + 1)
    for n in range(2, K + 1):
        d[n] = logn[n]*b[n] - mp.fsum(d[r]*b[n//r] for r in divs[n])
    arr = np.array([float(x) for x in d])
    cache[key] = arr
    return arr
ns = np.arange(K + 1, dtype=float); ns[0] = 1.0
lns = np.log(ns)
def Fvec(r):
    Fv = np.zeros(K + 1)
    for nd in r['nodes']:
        h = float(F(nd['height'])); a = float(F(nd['shift'])); w = float(F(nd['weight']))
        Fv += w*ns**(-h/2)*np.cos(a*lns/2)
    Fv[0] = Fv[1] = 0.0
    return Fv
def sup_torus(a_signed, starts=int(__import__('os').environ.get('STARTS','24')), seed=0):
    # maximise G(theta) = sum_n a_signed[n] cos(E[n].theta)
    idx = np.nonzero(a_signed)[0]
    A = a_signed[idx]; Ee = E[idx]
    def negG(th):
        ph = Ee @ th
        return -np.sum(A*np.cos(ph)), (Ee.T @ (A*np.sin(ph)))
    rng = np.random.default_rng(seed)
    gauge = np.array([0.0 if a_signed[p] >= 0 else math.pi for p in primes])
    best = None
    for s in range(starts):
        th0 = gauge.copy() if s == 0 else gauge + rng.normal(0, 0.6 if s < starts//2 else 2.0, len(primes))
        res = minimize(negG, th0, jac=True, method='L-BFGS-B', options={'maxiter': int(__import__('os').environ.get('MAXIT','3000')), 'gtol': 1e-12})
        if best is None or res.fun < best.fun: best = res
    return -best.fun, float(np.sum(np.abs(A)))
out = []
for i in rows:
    r = plan[i]; dr = direct[i]; hc = dr['head_construction']
    tl, tr = F(r['time'][0]), F(r['time'][1])
    Fv = Fvec(r); sgn = np.sign(Fv)
    dR = dcoef(tr)
    Hm = 0.5*np.sum(dR[2:]*ns[2:]**(-0.5 - float(tl)*c0)*np.abs(Fv[2:]))
    Hp = 0.5*np.sum(dR[2:]*ns[2:]**(-0.5 - float(tr)*c0)*np.abs(Fv[2:]))
    if 'raw_head_bounds' in hc:
        rawm, rawp = [float(F(x)) for x in hc['raw_head_bounds']]
        Cfree = float(F(hc['composite_credit']))
        Cclock = float(F(hc['old_angular_credit'])) + float(F(hc['additional_clock_credit']))
        mode = 'retained'
    else:
        rawm, rawp = [float(F(x)) for x in hc['raw_cap_sums']]
        Cfree = float(F(hc['fixed_credit'])) + float(F(hc['coordinate_credit']))
        Cclock = float(F(hc['orbit_credit']))
        mode = 'shared'
    res = dict(row=i, mode=mode, t=[str(tl), str(tr)], my_raw_first=Hm, their_raw_first=rawm, my_raw_last=Hp, their_raw_last=rawp,
               credit_phase_free=Cfree, credit_relation=Cclock)
    for lab, tt in (('t-', tl), ('t+', tr)):
        dt = dcoef(tt)
        a = 0.5*dt*ns**(-0.5 - float(tt)*c0)*np.abs(Fv); a[:2] = 0
        sup, absum = sup_torus(sgn*a)
        H = rawm if lab == 't-' else rawp
        res[lab] = dict(abs_sum=absum, torus_sup=sup, charged_after_free=H - Cfree, charged_after_all=H - Cfree - Cclock,
                        slack_vs_free=H - Cfree - sup, slack_vs_all=H - Cfree - Cclock - sup)
    out.append(res)
    print(json.dumps(res), flush=True)
bad = [o['row'] for o in out if o['t-']['slack_vs_free'] < 0 or o['t+']['slack_vs_free'] < 0]
rel = [o['row'] for o in out if o['t-']['slack_vs_all'] < 0 or o['t+']['slack_vs_all'] < 0]
print('ROWS_WHERE_TORUS_SUP_EXCEEDS_PHASE_FREE_CHARGE', bad)
print('ROWS_WHERE_FLOOR_NEEDS_RELATION_CREDIT', rel)

# C3b: falsification test of the common-clock (angular) credit of the retained rows of Zenodo 23752103.
# Their chain: head(theta) <= H - cc - sum_p c_p (1-cos theta_p)  and  sum_p c_p(1-cos theta_p) >= d  on {|k_j.theta| >= rho_j}.
# So for every gauged principal phase vector theta with |k_j.theta| >= rho_j (all j) the TRUE head must be <= H - cc - d.
# We maximise the true head (t = t-, actual heat coefficients d_{t-}) over that set for every sign pattern (SLSQP),
# relation primes in [-pi,pi], other primes at their unconstrained optimum.  A value > H - cc - d would falsify the credit.
# usage: uv run --no-project --with numpy --with scipy python3 clock_check.py PLAN.json DIRECT.json JOINT_ANGULAR.json rows...
import signal; signal.alarm(1500)
import json, sys, math, itertools
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
plan = {r['row']: r for r in json.load(open(sys.argv[1]))['rows']}
direct = {r['row']: r for r in json.load(open(sys.argv[2]))['rows']}
ja = json.load(open(sys.argv[3]))
rows = [int(x) for x in sys.argv[4:]]
primes31 = ja['prime_list']; cat = ja['relation_catalogue']; recs = {r['row']: r for r in ja['records']}
K = 4096; c0 = 0.5*math.log(690988)
spf = list(range(K + 1))
for p in range(2, int(K**0.5) + 1):
    if spf[p] == p:
        for m in range(p*p, K + 1, p):
            if spf[m] == m: spf[m] = p
primes = [p for p in range(2, K + 1) if spf[p] == p]; pidx = {p: i for i, p in enumerate(primes)}
E = np.zeros((K + 1, len(primes)))
for n in range(2, K + 1):
    m = n
    while m > 1:
        p = spf[m]; E[n, pidx[p]] += 1; m //= p
divs = [[] for _ in range(K + 1)]
for d in range(2, K + 1):
    for m in range(2*d, K + 1, d): divs[m].append(d)
ns = np.arange(K + 1, dtype=float); ns[0] = 1; lns = np.log(ns)
def dcoef(t):
    b = np.exp(t*lns**2/4); d = np.zeros(K + 1)
    for n in range(2, K + 1):
        d[n] = lns[n]*b[n] - sum(d[r]*b[n//r] for r in divs[n])
    return d
out = []
for i in rows:
    r = plan[i]; dr = direct[i]; hc = dr['head_construction']; jr = recs[i]
    tl = float(F(r['time'][0]))
    Fv = np.zeros(K + 1)
    for nd in r['nodes']:
        Fv += float(F(nd['weight']))*ns**(-float(F(nd['height']))/2)*np.cos(float(F(nd['shift']))*lns/2)
    Fv[:2] = 0
    s = np.sign(Fv)
    gp = np.array([s[p] for p in primes])           # gauge g_p = sign F(p)
    gn = np.ones(K + 1)
    for n in range(2, K + 1): gn[n] = np.prod(gp**E[n])
    a = 0.5*dcoef(tl)*ns**(-0.5 - tl*c0)*np.abs(Fv); a[:2] = 0
    A = s*gn*a                                        # head(theta) = sum_n A_n cos(E_n.theta)
    idx = np.nonzero(A)[0]; Ai = A[idx]; Ei = E[idx]
    def head(th): return float(np.sum(Ai*np.cos(Ei @ th)))
    def neg(th):
        ph = Ei @ th; return -np.sum(Ai*np.cos(ph)), Ei.T @ (Ai*np.sin(ph))
    # unconstrained sup (start at gauge optimum theta=0)
    res0 = minimize(neg, np.zeros(len(primes)), jac=True, method='L-BFGS-B', options={'maxiter': 2000})
    th_free = res0.x; sup_free = -res0.fun
    Hm = float(F(hc['raw_head_bounds'][0])); cc = float(F(hc['composite_credit'])); d0 = float(F(jr['credit']))
    dtot = d0 + float(F(hc['additional_clock_credit']))
    rel = jr['relations']; rho = [float(F(ang['rho_lower'])) for ang in jr['angles']]
    Kmat = np.zeros((len(rel), len(primes)))
    for j, ident in enumerate(rel):
        for p, k in zip(primes31, cat[ident]['exponents']): Kmat[j, pidx[p]] = k
    R = [pidx[p] for p in primes31]
    best = (-1e9, None, None)
    for signs in itertools.product((1, -1), repeat=len(rel) - 1):
        sg = np.array((1,) + signs)
        cons = [{'type': 'ineq', 'fun': (lambda th, j=j: sg[j]*(Kmat[j] @ th) - rho[j]), 'jac': (lambda th, j=j: sg[j]*Kmat[j])} for j in range(len(rel))]
        bounds = [(-math.pi, math.pi) if k in R else (None, None) for k in range(len(primes))]
        # start: free optimum, then nudge relation primes toward feasibility along the minimum-capacity direction
        x0 = th_free.copy()
        for j in range(len(rel)):
            viol = rho[j] - sg[j]*(Kmat[j] @ x0)
            if viol > 0: x0 += sg[j]*Kmat[j]*viol/np.dot(Kmat[j], Kmat[j])
        x0[R] = np.clip(x0[R], -math.pi, math.pi)
        res = minimize(neg, x0, jac=True, method='SLSQP', bounds=bounds, constraints=cons, options={'maxiter': 400, 'ftol': 1e-13})
        th = res.x
        ok = all(sg[j]*(Kmat[j] @ th) >= rho[j] - 1e-9 for j in range(len(rel))) and np.all(np.abs(th[R]) <= math.pi + 1e-9)
        if ok:
            # use the principal-angle condition |k.theta| >= rho (sign-pattern-free) as a final check
            hv = head(th)
            if hv > best[0]: best = (hv, signs, float(min(abs(Kmat[j] @ th) - rho[j] for j in range(len(rel)))))
    o = dict(row=i, t=tl, relations=len(rel), sup_free=sup_free, raw_H=Hm, cc=cc, d_old=d0, d_total=dtot,
             bound_old=Hm - cc - d0, bound_total=Hm - cc - dtot, constrained_max=best[0], min_constraint_slack=best[2],
             margin_old=Hm - cc - d0 - best[0], margin_total=Hm - cc - dtot - best[0])
    out.append(o); print(json.dumps(o), flush=True)
print('CLOCK_FALSIFIED_ROWS', [o['row'] for o in out if o['margin_old'] < 0])

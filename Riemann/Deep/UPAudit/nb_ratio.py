# N-B: Lemma Z's pointwise bound on real zeta zeros (RH holds there, so W = W_{-y} = On):
# R+ = |M(1/2+y+it)| / (1+On)^k,  R- = |M(1/2-y+it)| / (X^y (1+On)^k), t on a grid in [T, 2T] (2T <= 19900), A = 12.
# Also sup |M(1/2-y+it)|/X^y (trivial bound X^{1/2}·#terms-ish) to see whether the X^y factor is right.
import numpy as np, sys, json
from math import log
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit')
from na_identity import Case, ZEROS
def coeffs(c):
    ns, cs = [], []
    def rec(start, lg, coef, mu):
        ns.append(lg); cs.append(mu*complex(c.f(lg/c.L))*coef)
        for i in range(start, len(c.P)):
            q = c.P[i]; l2 = lg+log(q)
            if l2 >= c.lX: break
            rec(i+1, l2, coef*c.rho(log(q)/c.L), -mu)
    rec(0, 0.0, 1.0, 1); return np.array(ns), np.array(cs)
def On(c, t, A=12):
    lo = np.searchsorted(ZEROS, t-4.0); hi = np.searchsorted(ZEROS, t+4.0)
    out = np.empty(len(t))
    for i in range(len(t)):
        x = c.L*np.abs(ZEROS[lo[i]:hi[i]]-t[i]); out[i] = np.sum(np.minimum(1.0, np.maximum(x, 1e-300)**(-A)))
    return out
res = []
cfgs = [(0.56, 0.10, 0.66), (0.60, 0.15, 0.70)]
Ts = [1e3, 3e3, 9.9e3] if len(sys.argv) < 2 else [1e3]
for (th, eta, b) in cfgs:
    for T in Ts:
        c = Case(T, th, eta, b); lg, cs = coeffs(c); X = T**th; k = int(th/eta)
        t = np.arange(T, 2*T, 0.05 if len(sys.argv) < 2 else 0.5); on = On(c, t)
        row = dict(T=T, th=th, eta=eta, k=k, terms=len(cs), jmax_real=c.jmax)
        for y in (0.0, 0.25, 0.5):
            Mp = np.empty(len(t)); Mm = np.empty(len(t))
            for i in range(0, len(t), 2000):
                E = np.exp(-1j*np.outer(t[i:i+2000], lg))
                Mp[i:i+2000] = np.abs(E @ (cs*np.exp(-(0.5+y)*lg))); Mm[i:i+2000] = np.abs(E @ (cs*np.exp(-(0.5-y)*lg)))
            row[f'y{y}'] = dict(supRp=float(np.max(Mp/(1+on)**k)), supRm=float(np.max(Mm/(X**y*(1+on)**k))),
                                supRp1=float(np.max(Mp/(1+on))), supRm1=float(np.max(Mm/(X**y*(1+on)))),
                                supMm_over_Xy=float(np.max(Mm)/X**y), q99Mm_over_Xy=float(np.quantile(Mm, 0.99)/X**y),
                                trivial_over_Xy=float(np.sum(np.abs(cs)*np.exp(-(0.5-y)*lg))/X**y))
        res.append(row); print(json.dumps(row), flush=True)
json.dump(res, open('/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit/logs/nb_ratio.json', 'w'), indent=1)

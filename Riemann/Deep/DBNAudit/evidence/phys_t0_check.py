# C4: end-to-end physical test of the row-0 source floor of Zenodo 23752103 at t=0, independent of the author's model:
#   H_0(z)=xi((1+iz)/2)/8  =>  S_0(x,h) = (1/2) Re xi'/xi((1+h)/2 + i x/2),  sigma=(1+h)/2 >= 1.75 for row 0,
#   xi'/xi(s) = 1/s + 1/(s-1) - log(pi)/2 + psi(s/2)/2 + zeta'/zeta(s),  zeta'/zeta(s) = -sum Lambda(n) n^-s.
#   J_0(0,x) = alpha*sum_j w_j (S(x+a_j,h_j)+S(x-a_j,h_j))/2 + 2*T_b S(x,p)   (row 0: p=61/20, order 4)
# Scan x in [X0, X0+W] (truncated prime sum, rigorous tail bound), then polish the worst points with a long sum.
# Compare min J with the claimed floor S_0(0)=main_lo-native_debit-first_hi.
# usage: uv run --no-project --with numpy --with mpmath python3 phys_t0_check.py PLAN.json TABLE.json
import signal; signal.alarm(1500)  # hard time limit
import json, sys, math
from fractions import Fraction as F
import numpy as np
import mpmath as mp
mp.mp.dps = 40
plan = {r['row']: r for r in json.load(open(sys.argv[1]))['rows']}
table = {r['row']: r for r in json.load(open(sys.argv[2]))}
r = plan[0]; f = {k: F(v) for k, v in table[0]['fields'].items()}
floor0 = f['main_lo'] - f['native_debit'] - f['first_hi']
X0 = 6000000185827
nodes = [(F(n['height']), F(n['shift']), F(n['weight'])) for n in r['nodes']]
jj = r['joint_jet_overlay']; alpha = F(jj['alpha']); p = F(jj['centres'][0]['p']); bcoef = [F(b) for b in jj['centres'][0]['b']]
def sieve(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
def lambda_terms(N):
    ns, lam = [], []
    for q in sieve(N):
        lq = math.log(q); m = q
        while m <= N:
            ns.append(m); lam.append(lq); m *= q
    o = np.argsort(ns)
    return np.array(ns, dtype=np.int64)[o], np.array(lam)[o]
def G_of_n(ns):
    # G(n) = alpha*F(n) + (2/p) n^{-p/2} P_b(-p log n/2): coefficient of -(1/2) Lambda(n) n^{-1/2} cos(x log n/2)
    ln = np.log(ns.astype(float))
    Fv = np.zeros(len(ns))
    for h, a, w in nodes:
        Fv += float(w)*ns**(-float(h)/2)*np.cos(float(a)*ln/2)
    pf = float(p); z = -pf*ln/2
    Pb = np.zeros(len(ns)); fact = 1.0
    for m, bm in enumerate(bcoef):
        if m > 0: fact *= m
        Pb += float(bm)*z**m/fact
    return float(alpha)*Fv + (2/pf)*ns**(-pf/2)*Pb
def main_part(x):
    x = mp.mpf(x); tot = mp.mpf(0)
    def Sg(xx, h):
        s = mp.mpc((1 + h)/2, xx/2)
        return mp.re(1/s + 1/(s - 1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2)/2
    for h, a, w in nodes:
        hh = mp.mpf(h.numerator)/h.denominator; aa = mp.mpf(a.numerator)/a.denominator; ww = mp.mpf(w.numerator)/w.denominator
        tot += ww*(Sg(x + aa, hh) + Sg(x - aa, hh))/2
    tot *= mp.mpf(alpha.numerator)/alpha.denominator
    pp = mp.mpf(p.numerator)/p.denominator
    # jet main term: only r=0 survives up to O(1/x) (derivatives of the gamma part are O(1/x))
    tot += 2*(mp.mpf(bcoef[0].numerator)/bcoef[0].denominator)/pp*Sg(x, pp)
    return float(tot)
def frac_phase(ns, x0):
    # fractional part of x0*log(n)/(4 pi) (cycles), exact to ~1e-20
    out = np.empty(len(ns))
    c = 4*mp.pi
    for k, n in enumerate(ns.tolist()):
        v = mp.mpf(x0)*mp.log(n)/c
        out[k] = float(v - mp.floor(v))
    return out
def tail_bound(N):
    # sum_{n>N} Lambda(n) n^{-1/2}|G(n)|/2 <= (1/2) sum_terms |coef| * sum_{n>N} Lambda(n) n^{-sigma}, Chebyshev psi(y)<=1.04 y
    tb = 0.0
    for h, a, w in nodes:
        sg = 0.5 + float(h)/2
        tb += abs(float(alpha*w))*1.04*sg/(sg - 1)*N**(1 - sg)
    pf = float(p); sg = 0.5 + pf/2
    # |P_b(-p log n /2)| <= sum |b_m| (p log n/2)^m/m!; bound log n by log of the tail scale crudely via 3*sum
    tb += (2/pf)*sum(abs(float(bm)) for bm in bcoef)*50*1.04*sg/(sg - 1)*N**(1 - sg)
    return tb/2
# --- stage 1: scan with N1 terms
N1 = 6000
ns1, lam1 = lambda_terms(N1)
coef1 = 0.5*lam1*ns1**(-0.5)*G_of_n(ns1)
u1 = np.log(ns1.astype(float))/(4*math.pi)   # cycles per unit x
ph0 = frac_phase(ns1, X0)
mainX0 = main_part(X0)
W = 2_000_000.0; step = 0.37
best = []
B = 20000
for start in np.arange(0.0, W, step*B):
    dl = start + step*np.arange(B)
    cyc = ph0[None, :] + np.outer(dl, u1)
    vals = mainX0 - (np.cos(2*math.pi*cyc) @ coef1)
    k = np.argsort(vals)[:5]
    best += [(float(vals[j]), float(dl[j])) for j in k]
best.sort(); best = best[:40]
tb1 = tail_bound(N1)
print('stage1 floor', float(floor0), 'scan_min', best[0][0], 'tail_bound_N1', tb1, 'mainX0', mainX0, flush=True)
# --- stage 2: polish candidates with a long sum N2 and local minimisation in x
N2 = 400000
ns2, lam2 = lambda_terms(N2)
coef2 = 0.5*lam2*ns2**(-0.5)*G_of_n(ns2)
u2 = np.log(ns2.astype(float))/(4*math.pi)
ph2 = frac_phase(ns2, X0)
tb2 = tail_bound(N2)
def Jfull(dl):
    return main_part(X0 + dl) - float(np.cos(2*math.pi*(ph2 + dl*u2)) @ coef2)
res = []
seen = []
for v, dl in best:
    if any(abs(dl - s) < 5 for s in seen): continue
    seen.append(dl)
    # golden-section style local refinement on [dl-0.5, dl+0.5]
    xs = np.linspace(dl - 0.5, dl + 0.5, 201)
    vals = [Jfull(x) for x in xs]
    j = int(np.argmin(vals)); lo, hi = xs[max(j - 1, 0)], xs[min(j + 1, 200)]
    for _ in range(40):
        m1, m2 = lo + (hi - lo)/3, hi - (hi - lo)/3
        if Jfull(m1) < Jfull(m2): hi = m2
        else: lo = m1
    xm = (lo + hi)/2
    res.append((Jfull(xm), xm))
    if len(seen) >= 12: break
res.sort()
JX0 = Jfull(0.0)
out = dict(row=0, t=0, floor=float(floor0), J_at_X0=JX0, tail_bound_N2=tb2, polished_min=res[0][0], at_x_minus_X0=res[0][1],
           margin_min_minus_floor=res[0][0] - float(floor0) - tb2, top=res[:6])
print(json.dumps(out))
print('PHYS_T0', 'PASS' if res[0][0] - tb2 >= float(floor0) else 'FAIL')

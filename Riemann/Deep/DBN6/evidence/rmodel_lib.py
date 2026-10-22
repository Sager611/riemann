# Non-rigorous model of Raymond's H(7/8) chain (Zenodo 23272725) from the two-provider source in Planat's
# arXiv 2610.18214 ancillary files, and of the gain from adding a vertical-derivative term (HB Lemma 10 idea,
# translated as in Planat's pilot) with an ASSUMED physical bound G >= mean d_h S(x +- rh, h, t).
import json, numpy as np, mpmath as mp, sys
from fractions import Fraction as Fr
import os
base = os.environ.get('DBN6_SRC', '../pl/anc/certificate/inputs/original/')  # Planat arXiv 2610.18214 ancillary: anc/certificate/inputs/original/
src = {k: json.load(open(base + 'SOURCE_%s.json' % k)) for k in ('RETUNED', 'OLD_ENHANCED')}
f = lambda s: float(Fr(s))
kern = [(f(k['r']), f(k['rho']), f(k['lambda']), f(k['self_cap'])) for k in src['RETUNED']['kernels']]
def rkey(r):
    fr = Fr(r).limit_denominator(1000); return 'r_%d_%d' % (fr.numerator, fr.denominator)
cells = []
for i in range(608):
    opts = []
    for k in src:
        row = src[k]['rows'][i]; h = f(src[k]['plan']['rows'][i]['height'])
        fl = {kk: f(v) for kk, v in row.items() if kk.startswith('r_')}
        opts.append((h, fl))
    cells.append(opts)
L4 = float(mp.log(mp.mpf(6000000185827)/(4*mp.pi)))/4
def K(h, d, b): return (h-b)/(d*d+(h-b)**2) + (h+b)/(d*d+(h+b)**2)
def dK(h, d, b): return sum((d*d - c*c)/(d*d + c*c)**2 for c in (h-b, h+b))
dg = np.concatenate([np.linspace(0, 3, 401), np.linspace(3, 60, 200), [1e3, 1e4]])
def lam_table(r, rho, mus):
    LHS=[]; Jv=[]; dJv=[]
    for y in rho*np.linspace(0.02, 1, 40):
        for bb in np.linspace(0, 1, 31):
            b = bb*y; d = dg if bb < 1 else dg[dg > 1e-9]
            LHS.append(2*K(y, d, b)/y); Jv.append(0.5*(K(1, d+r, b)+K(1, d-r, b))); dJv.append(0.5*(dK(1, d+r, b)+dK(1, d-r, b)))
    LHS = np.concatenate(LHS); Jv = np.concatenate(Jv); dJv = np.concatenate(dJv)
    return [np.min((LHS + mu*dJv)/Jv) for mu in mus]
def own(r, rho):
    w = np.linspace(0, rho*rho, 400); y = np.sqrt(w)
    M = np.max(K(1, r, y)); g = np.min(0.5*(dK(1, r, y)+dK(1, -r, y)))
    return M, g
rs = sorted({rr for rr, _, _, _ in kern})
rhos = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
mus = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
dtab = []
for r in rs:
    for rho in rhos:
        lams = lam_table(r, rho, mus); M, g = own(r, rho)
        for mu, lam in zip(mus, lams):
            dtab.append((r, rho, 0.997*lam, mu, M*1.001, g*1.001))
def psi2(s): return float(mp.zeta(s, derivative=2)/mp.zeta(s) - (mp.zeta(s, derivative=1)/mp.zeta(s))**2)
gcache = {}
def Gmodel(t, h, fac, add):
    key = (round(t, 9), round(h, 9))
    if key not in gcache: gcache[key] = 0.25*psi2((1+h)/2 + L4*t)
    return fac*gcache[key] + add
def run(deriv, fac=1.0, add=0.0, sub=32, q0=0.5625+1e-6, qe=1e-8):
    q = q0; dt = 1/4000/sub
    for i in range(608):
        t0 = i/4000
        for j in range(sub):
            best = 0
            for h, fl in cells[i]:
                if deriv:
                    G = Gmodel(t0, h, fac, add)
                    for r, rho, lam, mu, M, g in dtab:
                        if q > (rho*h)**2 or rkey(r) not in fl: continue
                        c = (2*lam/h)*(fl[rkey(r)] - M/h) - 2*mu*(G - g/h**2)
                        best = max(best, c)
                else:
                    for r, rho, lam, M in kern:
                        if q > (rho*h)**2 or rkey(r) not in fl: continue
                        best = max(best, (2*lam/h)*(fl[rkey(r)] - M/h))
            c = best; assert c > 0
            qn = (q + 2/c)*np.exp(-c*dt) - 2/c
            t = t0 + j*dt
            if qn <= qe:
                tau = np.log(1 + c*(q-qe)/(2+c*qe))/c
                return t + tau + qe/2
            q = qn
    return None

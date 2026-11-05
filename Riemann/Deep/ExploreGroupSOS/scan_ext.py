# scan_ext.py: count negative lattice weights with and without the separable prime-power extension (c real).
import sys, math, numpy as np
from scipy.optimize import linprog
from weil_lat import atoms, m_sym, hplus, Ttail, pole
def beta_p(p, c, extra, ngrid=4096):
    ks_fixed = [k for k in range(1, 200) if p**k < c]
    k0 = (ks_fixed[-1] + 1) if ks_fixed else 1
    ks_free = list(range(k0, k0 + extra))
    th = np.linspace(0, math.pi, ngrid)
    fixed = -2 * sum((math.log(p) / p**(k / 2)) * np.cos(k * th) for k in ks_fixed)
    nf = len(ks_free)
    A = np.hstack([-np.ones((ngrid, 1)), -2 * np.cos(np.outer(th, ks_free))])
    res = linprog(np.concatenate([[1.0], np.zeros(nf)]), A_ub=A, b_ub=fixed, bounds=[(None, None)] * (nf + 1), method='highs')
    return float(res.x[0]), [(k, float(e)) for k, e in zip(ks_free, res.x[1:])]
def run(c, extra, Pfac=None):
    L = math.log(c); at = atoms(c)
    primes = sorted({p for p in range(2, int(c) + 1) if all(p % q for q in range(2, int(p**0.5) + 1)) and p < c})
    ext = []; btot = 0; Pneed = 2 * L
    for p in primes:
        b, co = beta_p(p, c, extra); btot += b
        for k, e in co:
            ext.append((k * math.log(p), e)); Pneed = max(Pneed, 2 * k * math.log(p))
    out = {}
    for label, P, E in [('P=2L no ext', 2 * L, []), (f'P={Pneed:.2f} no ext', Pneed, []), (f'P={Pneed:.2f} sep ext', Pneed, ext)]:
        sig = 2 * sum(w for _, w in at) if not E else btot
        Eb = 0.0
        # envelope crossing: hplus(r) - |T| - |pole| - sig > 0
        r = 20.0
        while True:
            env = hplus(r) - 2 * sum(math.exp(-(2*n+0.5)*L) * ((2*n+0.5) + r) / ((2*n+0.5)**2 + r*r) for n in range(40)) - 4 * math.exp(L/2) / (r - 0.5) - sig
            if env > 0: break
            r *= 1.02
        K = int(r * P / (2 * math.pi)) + 2
        neg = 0; negs = []; minw = 1e9
        for k0 in range(0, K + 1, 200000):
            k = np.arange(k0, min(K + 1, k0 + 200000)); rr = 2 * math.pi * k / P
            w = m_sym(rr, L, at)
            for x, e in E: w = w + 2 * e * np.cos(rr * x)
            neg += int(np.sum(w < 0)); minw = min(minw, float(w.min()))
            negs += list(rr[w < 0][:5000])
        negs = np.array(negs)
        lowcut = 2 * math.pi * 3 * c   # 'low' frequency band ~ the near-null band T* ~ 2 pi e c
        print(f'c={c:g} {label}: envelope constant {sig:.3f}, nodes checked {K}, negative weights {neg} (r<{lowcut:.0f}: {int(np.sum(negs < lowcut))}, beyond: {int(np.sum(negs >= lowcut))}; largest negative r = {negs.max() if len(negs) else 0:.1f}); min weight {minw:.3f}', flush=True)
c = float(sys.argv[1]); extra = int(sys.argv[2])
run(c, extra)

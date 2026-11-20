#!/usr/bin/env python3
# ceilaudit5 torus.py: is Raymond's finite-head triangle bound attained on the torus of prime phases?
# head(x) = (1/2) sum_{n<=4096} d_t(n) n^{-1/2-t c0} v(log n) cos(<e(n), theta>), theta_p = (x/2) log p mod 2pi (Kronecker: dense)
# worst case = max over the torus of Q(theta) = sum_n beta_n v_n cos(<e(n),theta>); triangle = sum |beta_n v_n|.
import sys, math, json
import numpy as np
sys.path.insert(0, 'Riemann/Deep/DBNHigh')
import model as M
K = 4096
def factor(n):
    f = {}; m = n; p = 2
    while p*p <= m:
        while m % p == 0: f[p] = f.get(p, 0) + 1; m //= p
        p += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f
FAC = [None, None] + [factor(n) for n in range(2, K + 1)]
def coeffs(row, t, tcoef):
    d, ln = M.dt_coeffs(tcoef, N=K)
    n = np.arange(2, K + 1)
    beta = 0.5*d[n]*np.exp(-(0.5 + t*M.c0)*ln[n])
    v = M.v_row(row, ln[n])
    return n, beta*v
def analyse(c, label):
    n = np.arange(2, K + 1)
    tri = float(np.sum(np.abs(c))); q0 = float(np.sum(c))
    pp = {}; comp = 0.0; compsigned = 0.0
    for i, m in enumerate(n):
        f = FAC[m]
        if len(f) == 1:
            (p, k), = f.items(); pp.setdefault(p, {})[k] = c[i]
        else:
            comp += abs(c[i]); compsigned += c[i]
    th = np.linspace(0, 2*math.pi, 200001)
    sep = 0.0; sepmin = 0.0; gaps = []; best = {}
    for p, ks in pp.items():
        if len(ks) == 1:
            a = list(ks.values())[0]; sep += abs(a); sepmin += -abs(a); best[p] = 0.0 if a > 0 else math.pi; continue
        f = sum(a*np.cos(k*th) for k, a in ks.items())
        j = int(np.argmax(f)); mx = float(f[j]); mn = float(np.min(f))
        sep += mx; sepmin += mn; best[p] = float(th[j])
        g = sum(abs(a) for a in ks.values()) - mx
        if g > 1e-6: gaps.append((p, round(g, 6), round(sum(abs(a) for a in ks.values()) + mn, 6)))
    upper = sep + comp
    # lower bound: evaluate full Q at the separable maximiser, then coordinate ascent
    ths = dict(best)
    def Q(ths):
        return float(sum(c[i]*math.cos(sum(e*ths[p] for p, e in FAC[m].items())) for i, m in enumerate(n)))
    lo = Q(ths)
    print('%s: triangle %.6f  Q(0) %.6f  torus-max upper (sep pp max + composite triangle) %.6f  lower (Q at sep argmax) %.6f  composite triangle %.6f signed %.6f' % (label, tri, q0, upper, lo, comp, compsigned))
    print('   gap triangle - upper = %.6f (%.3f%%);  separable min side: -min = %.6f' % (tri - upper, 100*(tri - upper)/tri, -(sepmin - comp)))
    print('   primes with sign-conflict gaps (p, gap at max, gap at min):', gaps[:20])
    return dict(tri=tri, q0=q0, upper=upper, lo=lo, comp=comp)
if __name__ == '__main__':
    rows = M.load_rows(sys.argv[1])
    out = {}
    for i in [int(x) for x in sys.argv[2].split(',')]:
        row = rows[i]
        for t in (row['t0'], row['t1']):
            n, c = coeffs(row, t, row['t1'])
            out['%d_%g' % (i, t)] = analyse(c, 'row %d t %.5f (file first_hi %.6f last_hi %.6f)' % (i, t, row['A'], row['B']))
    json.dump(out, open(sys.argv[3], 'w'), indent=1)

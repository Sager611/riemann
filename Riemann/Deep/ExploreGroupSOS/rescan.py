import sys, math, numpy as np
from weil_lat import atoms, m_sym, hplus
d = np.load(sys.argv[1]); sig = float(sys.argv[2]); c = float(sys.argv[3])
F, e, lg, P, L = d['F'], d['e'], d['lg'], float(d['P']), float(d['L'])
at = atoms(c)
r = 20.0
while True:
    envr = hplus(r) - 2 * sum(math.exp(-(2*n+0.5)*L) * ((2*n+0.5) + r) / ((2*n+0.5)**2 + r*r) for n in range(40)) - 4 * math.exp(L/2) / (r - 0.5) - sig
    if envr > 0: break
    r *= 1.01
K = int(r * P / (2 * math.pi)) + 2
k = np.arange(0, K + 1); rr = 2 * math.pi * k / P
w = m_sym(rr, L, at)
for a, ee in zip(F, e): w = w + 2 * ee * np.cos(rr * float(np.dot(a, lg)))
neg = w < 0
print(f'c={c:g} P={P:.4f}: certified envelope {sig}: scan to r={r:.1f}, nodes {K}, negative {int(neg.sum())}, last neg r {rr[neg].max():.1f}, '
      f'weights < 1: {int((w < 1).sum())}; beyond r_K every weight >= envelope margin > 0')

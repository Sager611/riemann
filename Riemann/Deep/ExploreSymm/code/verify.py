# independent check of a competitor: coefficients via np.poly, S on every line, roots via np.roots, Jensen gaps
import numpy as np, sys
d = np.load(sys.argv[1]); u = d['u']; v = d['v']; Pz = float(d['Pz']); xs = d['xs']; hl = d['hl']; M = int(d['M'])
wj = np.exp(u + 1j*v)
c = np.poly(1/wj)[::-1]            # prod (1 - w/w_j) = prod(-1/w_j)(w - w_j)... use roots 1/w_j of the reversed poly
c = np.poly(wj)[::-1]; c = c/c[0]   # p(w) = c0 + c1 w + ... with c0 = 1
rts = np.roots(c[::-1]); nin = int(np.sum(np.abs(rts) <= 1 + 1e-9))
worst = -1e9; rows = []
for x, h in zip(xs, hl):
    r = np.exp(-x/Pz); S = np.sum(np.abs(c)**2 * r**(2*np.arange(len(c))))
    nf = 4096; w = r*np.exp(2j*np.pi*np.arange(nf)/nf); pv = np.polyval(c[::-1], w)
    m = np.mean(np.log(np.abs(pv))); hc = 0.5*np.log(S)
    worst = max(worst, hc - h); rows.append((x, h, hc, m, hc - m))
print(f"{sys.argv[1]}: K={len(c)-1} roots |w|<=1: {nin} (target M={M}); nu = {nin/Pz:.4f}; max (h_comp - h_target) = {worst:+.2e}")
for (x, h, hc, m, J) in rows:
    if x in (-12, -6, -3, -2, -1.4, -1.1, -1, -0.5, 0, 0.5, 2) or abs(x + 1.24618) < 1e-6:
        print(f"  x={x:7.3f}  h_target={h:.4f}  h_comp={hc:.4f}  m_comp={m:.4f}  J_comp={J:.4f}")

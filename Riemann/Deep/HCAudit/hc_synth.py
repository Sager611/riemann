# hcaudit synthetic falsification of (iii') with the cert768 data (EPS = 1/1000): (A) the exact finite inequality
#   S = sum_{i,j} f^(x_i - x_j) <= n f^(0) + 2 f^(g_{n-1}) + 2[(n-2)c + A(g_{n-1}) - A(g_1) + beta(x_{n-1} - x_1 - (n-2))]
# on adversarial hard-core sequences (w = 1); (B) max over cyclic gap sequences of E_beta = avg[f^(g_i) + f^(g_i+g_{i+1}) - beta(g_i-1)]
# (<= c by (H) summed cyclically), coordinate ascent with restarts. Usage: EPS=1/1000 python3 hc_synth.py cert768.json -0.3478616532
import sys, json, os
from fractions import Fraction as Fr
import numpy as np
from flint import fmpq_poly, fmpq, arb, arb_poly, ctx
d = json.load(open(sys.argv[1])); c = float(sys.argv[2]); rng = np.random.default_rng(7)
lam = float(Fr(d['lam'])); D = d['D']; s = float(Fr(d['s'])); U = float(Fr(d['U'])); dg = float(Fr(d['dg']))
Bq = [Fr(b) for b in d['B']]; Bq[0] += Fr(os.environ.get('EPS', '0')); Bf = np.array([float(b) for b in Bq])
A = np.array([float(Fr(a)) for a in d['A']]); beta = float(Fr(d['beta'])); nodes = lam + dg*np.arange(len(A))
k = np.sqrt(2*np.pi)
def fh(u):
    v = k*s*np.asarray(u, float); h0, h1 = np.ones_like(v), 2*v; tot = Bf[0]*h0; hs = [h0, h1]
    for m in range(1, 2*D):
        hs.append(2*v*hs[-1] - 2*m*hs[-2])
    for n in range(1, D+1): tot = tot + Bf[n]*(-1)**n*hs[2*n]
    return s*tot*np.exp(-v*v/2)
# float vs Arb sanity at a few points
ctx.prec = 200; q2 = lambda x: fmpq(x.numerator, x.denominator)
H = [fmpq_poly([1]), fmpq_poly([0, 2])]
for n in range(1, 2*D): H.append(fmpq_poly([0, 2])*H[n] - 2*n*H[n-1])
Qp = sum((q2(Bq[n]*(-1)**n)*H[2*n] for n in range(D+1)), fmpq_poly([0])); QA = arb_poly([arb(Qp[i]) for i in range(Qp.degree()+1)])
fa = lambda u: float((arb(s)*QA(arb(k*s*u))*(-(arb(k*s*u)**2)/2).exp()).mid())
print('float-vs-Arb f^ max abs err on 9 points in [0,8]: %.2e' % max(abs(fh(u) - fa(u)) for u in np.linspace(0, 8, 9)))
Af = lambda x: np.interp(x, nodes, A, right=A[-1]); f0 = float(fh(0.0))
def pair_sum(x):   # all ordered pairs, f^ ~ 0 beyond 14
    S = len(x)*f0; kk = 1
    while kk < len(x) and (x[kk:] - x[:-kk]).min() < 14:
        S += 2*fh(x[kk:] - x[:-kk]).sum(); kk += 1
    return S
def bound(x):
    g = np.diff(x); n = len(x)
    return n*f0 + 2*fh(g[-1]) + 2*((n-2)*c + Af(g[-1]) - Af(g[0]) + beta*(x[n-2] - x[0] - (n-2)))
fams = {}
fams['Tonks mean 1'] = [lam + rng.exponential(1 - lam, 3000) for _ in range(40)]
fams['all gaps lambda'] = [np.full(500, lam)]
fams['runs r*lambda + long'] = [np.tile(np.r_[np.full(r, lam), r*(1 - lam) + 1], 400 // (r+1)) for r in range(1, 9)]
fams['alt lambda, 2-lambda'] = [np.tile([lam, 2 - lam], 300)]
fams['near A nodes'] = [nodes[rng.integers(0, len(A), 2000)] + rng.choice([-1e-9, 1e-9], 2000).clip(0) for _ in range(20)]
fams['argmax pair (lam,1.11709)'] = [np.tile([lam, 1.117091], 300)]
fams['uniform [lam, 2-lam]'] = [rng.uniform(lam, 2 - lam, 3000) for _ in range(20)]
fams['boundary cluster'] = [np.r_[np.full(50, lam), rng.uniform(lam, 3, 500), np.full(50, lam)] for _ in range(20)]
for name, gs in fams.items():
    worst = -1e18; per = []
    for gg in gs:
        x = np.r_[0.0, np.cumsum(gg)]; S, Bd = pair_sum(x), bound(x); worst = max(worst, S - Bd); per.append(S/len(x))
    print('(A) %-26s %3d seqs: max (S - bound) = %.4f (<= 0 required: %s); S/n in [%.4f, %.4f]' % (name, len(gs), worst, worst <= 1e-9, min(per), max(per)))
# (B) coordinate ascent on cyclic sequences
grid = np.arange(lam, 5.0, 0.0005)
def Ebeta(g): return np.mean(fh(g) + fh(g + np.roll(g, -1)) - beta*(g - 1))
best_all = -9; best_g = None
for n in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 30]:
    bn = -9
    for rest in range(12 if n > 1 else 1):
        g = rng.uniform(lam, 2.5, n) if rest else np.full(n, 1.0)
        for sweep in range(60):
            old = Ebeta(g)
            for i in rng.permutation(n):
                gp, gn = g[(i-1) % n], g[(i+1) % n]
                if n == 1: val = fh(grid) + fh(2*grid) - beta*(grid - 1)
                elif n == 2: val = fh(grid) + fh(grid + gn) + fh(gp + grid) - beta*(grid - 1)
                else: val = fh(grid) + fh(grid + gn) + fh(gp + grid) - beta*(grid - 1)
                g[i] = grid[val.argmax()]
            if Ebeta(g) <= old + 1e-13: break
        e = Ebeta(g)
        if e > bn: bn = e; bg = g.copy()
    if bn > best_all: best_all = bn; best_g = bg.copy()
    print('(B) n = %2d: max E_beta = %.7f (c_cert %.7f, slack %.2e); its mean gap %.4f, gaps %s' % (n, bn, c, c - bn, bg.mean(), np.round(np.sort(bg)[:8], 4)))
x = np.r_[0.0, np.cumsum(np.tile(best_g, max(1, 3000 // len(best_g))))]
print('(B) best cyclic E_beta overall %.7f <= c_cert: %s ; full pair energy per point of its periodic extension %.5f vs f^(0) + 2c + 2 beta (mean-1) = %.5f' % (
      best_all, best_all <= c, pair_sum(x)/len(x), f0 + 2*c + 2*beta*(best_g.mean() - 1)))

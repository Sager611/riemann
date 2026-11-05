# lp8.py -- cutting-plane LP for 8-point local data (critical line), window of arXiv 2610.08965 (second window)
# usage: python lp8.py [delta] [rounds] [seed] [mode]   (single-threaded; floating point; proposes data only)
import json, sys, time
import numpy as np
from scipy.optimize import linprog
t0 = time.time()
delta = float(sys.argv[1]) if len(sys.argv) > 1 else 1.3e-5
rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 30
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
mode = sys.argv[4] if len(sys.argv) > 4 else 'free'      # 'free' rewards, or 'k3' single reward at site 3
rng = np.random.default_rng(seed)
d = json.load(open('src/anc/certificate_line/inputs.json'))
a = np.array(d['window_coefficients'], float) / 1e6
theta = float(sys.argv[7]) if len(sys.argv) > 7 else 0.85
C0 = 6639737 / 5000000; beta = 13707941 / 2**28
import os
if os.environ.get('WINDOW'):   # trial window: comma-separated a_j (floats); C(f) by the closed form, beta numerically
    a = np.array([float(v) for v in os.environ['WINDOW'].split(',')])
    jj = np.arange(1, len(a) + 1)
    C0 = 4 / 3 + 0.5 * (a**2).sum() + (((-1.0) ** jj * a / jj**2 - a**2 / (4 * jj**2)).sum()) / np.pi**2 + 1e-9
    print('trial window: C(f) = %.10f, 1-sum|a| = %.6f' % (C0, 1 - np.abs(a).sum()))
NP = int(sys.argv[6]) if len(sys.argv) > 6 else 8
NG = NP - 1
pairs = [(i, j) for i in range(NP) for j in range(i + 1, NP)]
PI = np.array([p[0] for p in pairs]); PJ = np.array([p[1] for p in pairs])
knots = 0.75 + np.arange(14) / 8
def Kex(t):
    s = np.sinc(t)
    for j, aj in enumerate(a, 1):
        s = s + aj / 2 * (np.sinc(t - j) + np.sinc(t + j))
    return s
H = 2e-5; tg = np.arange(0, 60, H); K2g = Kex(tg) ** 2
if os.environ.get('WINDOW'):
    beta = K2g[tg <= theta].min() * (1 - 1e-6)
    uu = np.linspace(0, 0.5, 50001); ff = 1 + sum(aj * np.cos(2 * np.pi * j * uu) for j, aj in enumerate(a, 1))
    g1 = (np.sinc(theta * (uu - 0.5)) + np.sinc(theta * (uu + 0.5))) ** 2
    gam = (ff / g1).max(); Lam = 2 * gam * (1 + np.sinc(theta)) / theta
    print('trial window: beta %.6f, min f %.4f, gamma0 needed %.5f, Lambda %.5f (< 2 needed)' % (beta, ff.min(), gam, Lam))
def K2(t): return np.interp(t, tg, K2g)
def hats(g):  # g: (m,) -> (m,12) values of the interior hat functions k=1..12
    u = (g[:, None] - knots[None, 1:13]) * 8.0
    return np.maximum(0.0, 1.0 - np.abs(u))
NV = 1 + NP + len(pairs) + NG * 12
iA, iK, iW, iB = 0, 1, 1 + NP, 1 + NP + len(pairs)
def features(G):
    m = G.shape[0]
    y = np.concatenate([np.zeros((m, 1)), np.cumsum(G, axis=1)], axis=1)
    D = y[:, PJ] - y[:, PI]
    F = np.zeros((m, NV))
    F[:, iA] = G.sum(axis=1)
    F[:, iK:iK + NP] = -1.0
    F[:, iW:iW + len(pairs)] = K2(D)
    for r in range(NG):
        F[:, iB + 12 * r:iB + 12 * (r + 1)] = hats(G[:, r])
    return F
def slack(G, z):
    m = G.shape[0]
    y = np.concatenate([np.zeros((m, 1)), np.cumsum(G, axis=1)], axis=1)
    v = z[iA] * G.sum(axis=1) + K2(y[:, PJ] - y[:, PI]) @ z[iW:iW + len(pairs)] - z[iK:iK + NP].sum()
    for r in range(NG):
        kv = np.concatenate([[0.0], z[iB + 12 * r:iB + 12 * (r + 1)], [0.0]])
        v = v + np.interp(G[:, r], knots, kv, left=0.0, right=0.0)
    return v
MAXIT = 150
def refine(G, z, minstep=1e-5):
    G = G.copy(); m = G.shape[0]
    val = slack(G, z); step = np.full(m, 0.05)
    # move set: gap moves, point moves (g_{r-1}+h, g_r-h), and random directions
    E = [np.eye(NG)[r] for r in range(NG)]
    for r in range(1, NG):
        e = np.zeros(NG); e[r - 1] = 1; e[r] = -1; E.append(e / np.sqrt(2))
    E = np.array(E)
    it = 0
    while (step > minstep).any() and it < MAXIT:
        it += 1
        act = step > minstep
        improved = np.zeros(m, bool)
        dirs = np.concatenate([E, -E, rng.normal(size=(6, NG))])
        dirs[-6:] /= np.linalg.norm(dirs[-6:], axis=1)[:, None]
        for e in dirs:
            Gn = np.maximum(theta, G + step[:, None] * e[None, :])
            v = slack(Gn, z)
            b = act & (v < val - 1e-14)
            G[b] = Gn[b]; val[b] = v[b]; improved |= b
        step[act & ~improved] /= 2
    return G, val
def sample(m):
    G = np.empty((m, NG))
    k = m // 5
    cand = np.concatenate([knots[1:13], [1.04, 1.97, 1.04, 1.97, 2.95]])
    G[4 * k:] = cand[rng.integers(0, len(cand), size=(m - 4 * k, NG))]
    centers = np.array([1.04, 1.97, 2.95])
    idx = rng.choice(3, size=(2 * k, NG), p=[0.55, 0.4, 0.05])
    G[:2 * k] = centers[idx] + rng.normal(0, 0.03, size=(2 * k, NG))
    G[2 * k:3 * k] = theta + rng.uniform(0, 2.4, size=(k, NG))
    G[3 * k:4 * k] = theta + rng.exponential(0.5, size=(k, NG))
    return np.maximum(G, theta)
def paper_z():  # the paper's 8-point data, embedded in frames of NP >= 8 points
    z = np.zeros(NV); D = 1e6
    z[iA] = d['alpha'] / D
    z[iK:iK + 8] = np.array(d['h'][0]) / D
    p8 = [(i, j) for i in range(8) for j in range(i + 1, 8)]
    for q, (i, j) in enumerate(p8):
        z[iW + pairs.index((i, j))] = d['weights'][0][q] / D
    for r in range(7):
        z[iB + 12 * r:iB + 12 * (r + 1)] = np.array(d['cob_knots'][r][0][1:13]) / D
    return z
def sval(z):
    x = z[iK:iK + NP].sum(); return (2 - C0 - NG * z[iA]) / (1 - x)
def oracle(z, m=60000, keep=700, extra=None):
    G = sample(m)
    if extra is not None: G = np.concatenate([G, extra])
    v = slack(G, z)
    idx = np.argsort(v)[:keep]
    Gr, vr = refine(G[idx], z)
    return Gr, vr
if mode == 'checkpaper':
    z = paper_z()
    Gr, vr = oracle(z, m=100000, keep=1000)
    o = np.argsort(vr)
    print('paper data: s = %.7f, min slack %.3e' % (sval(z), vr[o[0]]))
    for k in o[:5]: print('  slack %.3e gaps %s' % (vr[k], np.round(Gr[k], 4)))
    print('time %.0fs' % (time.time() - t0)); sys.exit()
# LP pieces
Aub_budget = np.zeros((NG, NV))
for p, (i, j) in enumerate(pairs): Aub_budget[j - i - 1, iW + p] = 1.0
Aeq = np.zeros((12, NV))
for r in range(NG):
    for k in range(12): Aeq[k, iB + 12 * r + k] = 1.0
xrow = np.zeros(NV); xrow[iK:iK + NP] = 1.0
BB = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0015
bounds = [(0, 0.01)] + [(0, 0.03)] * NP + [(0, 2)] * len(pairs) + [(-BB, BB)] * (NG * 12)
if mode == 'k3':
    bounds[iK:iK + NP] = [(0, 0)] * NP; bounds[iK + 3] = (0, 0.03)
G0 = sample(20000)
G0 = np.concatenate([G0, oracle(paper_z())[0]])
Fs = features(G0); Gall = G0.copy()
s0 = sval(paper_z()); z = None
hist = []
for rd in range(rounds):
    for dk in range(1 if rd > 0 else 3):  # Dinkelbach on the current cut set
        c = np.zeros(NV); c[iA] = float(NG); c[iK:iK + NP] = -s0
        Aub = np.concatenate([-Fs, Aub_budget, xrow[None, :]])
        bub = np.concatenate([-np.full(Fs.shape[0], delta), np.full(NG, 2.0), [beta]])
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=np.zeros(12), bounds=bounds, method='highs')
        if res.status != 0: print('LP status', res.status, res.message); sys.exit(1)
        z = res.x; s1 = sval(z)
        if abs(s1 - s0) < 1e-10: s0 = s1; break
        s0 = s1
    Gr, vr = oracle(z, extra=G0[np.argsort(Fs[:G0.shape[0]] @ z)[:500]] if rd < 0 else None)
    viol = vr < delta - 2e-9
    hist.append((rd, s0, vr.min(), int(viol.sum()), Fs.shape[0]))
    print('round %d  s_LP %.7f  x %.6f  alpha %.6e  min slack %.3e  new cuts %d  rows %d  t %.0fs' % (rd, s0, z[iK:iK+NP].sum(), z[iA], vr.min(), viol.sum(), Fs.shape[0], time.time() - t0), flush=True)
    if not viol.any(): break
    o = np.argsort(vr)[:3]
    for k in o: print('   worst %.3e gaps %s' % (vr[k], np.round(Gr[k], 3)), flush=True)
    keeprow = (Fs @ z) < 0.004
    Fs = Fs[keeprow]; Gall = Gall[keeprow]
    newG = Gr[viol]
    # dedupe roughly
    _, ui = np.unique(np.round(newG, 3), axis=0, return_index=True)
    newG = newG[ui]
    Fs = np.concatenate([Fs, features(newG)]); Gall = np.concatenate([Gall, newG])
    G0 = newG
c = np.zeros(NV); c[iA] = float(NG); c[iK:iK + NP] = -s0
res = linprog(c, A_ub=np.concatenate([-Fs, Aub_budget, xrow[None, :]]), b_ub=np.concatenate([-np.full(Fs.shape[0], delta), np.full(NG, 2.0), [beta]]), A_eq=Aeq, b_eq=np.zeros(12), bounds=bounds, method='highs')
z = res.x; s0 = sval(z)
print('final re-solve on %d cuts: s_LP %.7f' % (Fs.shape[0], s0))
np.save('lp_z_%dpt_%s_%g_th%g.npy' % (NP, mode, delta, theta), z)
duals = -res.ineqlin.marginals[:Fs.shape[0]]
np.savez('lp_cuts_%dpt_%s_%g_th%g.npz' % (NP, mode, delta, theta), G=Gall, duals=duals, z=z, s=s0)
print('binding cuts with dual > 1e-9: %d' % (duals > 1e-9).sum())
Gr, vr = oracle(z, m=100000, keep=1500)
near = vr < delta + 2e-5
print('near-binding frames (slack < delta+2e-5): %d; smallest gap among them %.4f; largest %.4f' % (near.sum(), Gr[near].min() if near.any() else -1, Gr[near].max() if near.any() else -1))
print('min K^2 on [0,theta] (float) %.6f vs x %.6f' % ((Kex(np.linspace(0, theta, 20001))**2).min(), z[iK:iK+NP].sum()))
print('final s_LP %.7f (delta %g, mode %s), rewards %s' % (s0, delta, mode, np.round(z[iK:iK + NP] * 1e6, 1)))
print('W*1e6', np.round(z[iW:iW + len(pairs)] * 1e6).astype(int).tolist())
print('time %.0fs' % (time.time() - t0))

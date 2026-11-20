#!/usr/bin/env python3
"""npaframe CHANGE vs NPAFals/colgen_simple.py: set specs T<n0> (pairs at all levels <= N, words of length 3-4 only at levels <= n0),
Po (pairs <= N-1, open band as in the AF finite programme), To<n0> (Po + length 3-4 words at levels <= n0); final duals y,
feature list and law support saved to an .npz (argv[7]).
npafals CHANGE vs ExploreNPA/colgen_open.py: objective min E[#simple]/N (not max E[sum m^2]/N); sets add L4oo (length>=3 words only
|lam| <= N-2) and L4c; patterns = doubles + at most one triple for every N; reports minSimple and the heuristic lower bound;
pricing warm-starts half the Adam batch from the previous best positions, adds 24 columns per pattern, stops at 850 s wall clock;
3000 CUE(N) samples as initial columns (jobs launched after 23:00 UTC); vectorised gradient (jobs after 23:05 UTC).
Column generation, continuous positions, RH world (on-line points only), CUE(N) Diaconis-Shahshahani constraints E p_lam conj p_mu = d z_lam.
usage: colgen_simple.py N SET [iters] [B] [seed] [steps]"""
import sys, math, time, resource
import numpy as np
from scipy.optimize import linprog

def zl(lam):
    z = 1
    for j in set(lam):
        a = lam.count(j); z *= j**a * math.factorial(a)
    return z

def parts(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in parts(n-k, k): yield (k,) + r

def keep(S, N, n, lam, mu):
    L = len(lam) + len(mu)
    if S == 'P': return L <= 2
    if S == 'P3': return L <= 3
    if S == 'L4': return L <= 4
    if S == 'P3o': return L <= 3 and not (L >= 3 and n == N)
    if S == 'L4o': return L <= 4 and not (L >= 3 and n == N)
    if S == 'L4oo': return L <= 4 and not (L >= 3 and n >= N - 1)
    if S == 'L4c': return L <= 4 and len(lam) < 3 and len(mu) < 3
    if S == 'Po': return L <= 2 and n <= N - 1
    if S.startswith('To'): return (L <= 2 and n <= N - 1) or (3 <= L <= 4 and n <= int(S[2:]))
    if S.startswith('T'): return L <= 2 or (3 <= L <= 4 and n <= int(S[1:]))
    raise ValueError(S)

def featlist(N, S):
    F = []
    for n in range(1, N+1):
        P = list(parts(n))
        for i, lam in enumerate(P):
            for mu in P[i:]:
                if not keep(S, N, n, lam, mu): continue
                s = math.sqrt(zl(lam)*zl(mu))
                F.append((lam, mu, (zl(lam) if lam == mu else 0.0)/s, s))
    return F

class Feats:
    def __init__(self, N, S):
        self.N = N; self.F = featlist(N, S); nf = len(self.F)
        W = max(len(l)+len(m) for l, m, _, _ in self.F)
        self.idx = np.zeros((nf, W), dtype=int); self.cj = np.zeros((nf, W), dtype=bool)
        for f, (l, m, t, s) in enumerate(self.F):
            sl = [(k, False) for k in l] + [(k, True) for k in m]
            for j, (k, c) in enumerate(sl): self.idx[f, j] = k; self.cj[f, j] = c
        self.t = np.array([t for _, _, t, _ in self.F]); self.sc = np.array([s for _, _, _, s in self.F]); self.W = W
    def pvals(self, th, m):
        k = np.arange(1, self.N+1)
        E = np.exp(1j*th[:, :, None]*k[None, None, :])
        P = np.einsum('i,bik->bk', m, E)
        return np.concatenate([np.ones((th.shape[0], 1)), P], axis=1), E
    def q(self, Pe):
        Q = Pe[:, self.idx]
        return np.where(self.cj[None], np.conj(Q), Q)
    def values(self, th, m):
        Pe, _ = self.pvals(th, m); Q = self.q(Pe)
        return (np.prod(Q, axis=2).real) / self.sc[None]
    def g_grad(self, th, m, y):   # npafals: vectorised (prefix/suffix exclusive products, scatter by matmul); same values
        Pe, E = self.pvals(th, m); Q = self.q(Pe); B = th.shape[0]; W = self.W
        if not hasattr(self, 'Sel'):
            self.Sel = np.zeros((len(self.F)*W, self.N+1)); self.Sel[np.arange(len(self.F)*W), self.idx.reshape(-1)] = 1.0
        pre = np.ones_like(Q); suf = np.ones_like(Q)
        for s in range(1, W): pre[:, :, s] = pre[:, :, s-1]*Q[:, :, s-1]
        for s in range(W-2, -1, -1): suf[:, :, s] = suf[:, :, s+1]*Q[:, :, s+1]
        g = ((pre[:, :, W-1]*Q[:, :, W-1]).real / self.sc[None]) @ y
        ex = pre*suf; coef = y / self.sc
        c = 0.5*coef[None, :, None]*np.where(self.cj[None], np.conj(ex), ex)
        G = c.reshape(B, -1) @ self.Sel
        G = G[:, 1:]; k = np.arange(1, self.N+1)
        dth = -2.0*m[None, :]*np.einsum('bk,bik->bi', G*k[None, :], E).imag
        return g, dth

def adam(fe, m, y, B, steps, rng, lr=0.05, warm=None):
    th = rng.uniform(0, 2*np.pi, (B, len(m)))
    if warm is not None and len(warm):
        h = min(B // 2, len(warm)); th[:h] = warm[:h] + rng.normal(0, 0.05, (h, len(m)))
    th[:, 0] = 0.0
    mo = np.zeros_like(th); v = np.zeros_like(th)
    for t in range(1, steps+1):
        g, d = fe.g_grad(th, m, y); d[:, 0] = 0.0
        mo = 0.9*mo + 0.1*d; v = 0.999*v + 0.001*d*d
        th = th - lr*(mo/(1-0.9**t))/(np.sqrt(v/(1-0.999**t)) + 1e-12)
        if t == steps//2: lr *= 0.2
    g, _ = fe.g_grad(th, m, y)
    return th, g

def patterns(N):
    out = [[2]*d + [1]*(N-2*d) for d in range(0, N//2 + 1)]
    out += [[3] + [2]*d + [1]*(N-3-2*d) for d in range(0, (N-3)//2 + 1)]
    return [np.array(p, dtype=float) for p in out]

def main():
    N = int(sys.argv[1]); S = sys.argv[2]
    iters = max(150, int(sys.argv[3]) if len(sys.argv) > 3 else 150)   # npafals: cap 150 iterations (850 s wall clock)
    B = int(sys.argv[4]) if len(sys.argv) > 4 else 384
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    steps = int(sys.argv[6]) if len(sys.argv) > 6 else 300
    PRUNE = 7000; import gc; gc.collect()   # npafals: 4000 lost the CUE columns at N = 9 (stalled); 7000 from 23:38 UTC
    import os; WALL = float(os.environ.get('WALL', '850'))   # npaframe: wall clock via env WALL
    TOL = 2e-5   # npafals: stop when max reduced cost found < 2e-5 (jobs after 23:16 UTC; lb shows the precision)
    rng = np.random.default_rng(seed); fe = Feats(N, S); nf = len(fe.F); pats = patterns(N)
    simple = lambda m: float((m == 1).sum())/N
    print(f"N={N} S={S} features={nf} patterns={len(pats)} B={B} steps={steps}", flush=True)
    cols, objs = [], []
    for m in pats:
        th = rng.uniform(0, 2*np.pi, (200, len(m))); th[:, 0] = 0
        V = fe.values(th, m); cols.append(V); objs += [-simple(m)]*len(V)
    # npafals: 3000 CUE(N) samples (Haar eigenangles, all simple) as extra initial columns -> the master is feasible at once
    Z = (rng.normal(size=(3000, N, N)) + 1j*rng.normal(size=(3000, N, N)))/np.sqrt(2)
    Qm, Rm = np.linalg.qr(Z); dg = np.diagonal(Rm, axis1=1, axis2=2); Qm = Qm*(dg/np.abs(dg))[:, None, :]
    thc = np.angle(np.linalg.eigvals(Qm)); cols.append(fe.values(thc, np.ones(N))); objs += [-1.0]*3000
    A = np.vstack(cols).T; obj = np.array(objs); BIG = 20.0; t0 = time.time(); rmax = np.inf; r = None; WARM = {}
    for it in range(iters):
        nc = A.shape[1]
        Aeq = np.vstack([np.hstack([A, np.eye(nf), -np.eye(nf)]), np.hstack([np.ones((1, nc)), np.zeros((1, 2*nf))])])
        beq = np.concatenate([fe.t, [1.0]]); c = np.concatenate([-obj, BIG*np.ones(2*nf)])
        r = linprog(c, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
        lam = r.eqlin.marginals; y = -lam[:nf]; y0 = -lam[nf]
        val = float(obj @ r.x[:nc]); slack = float(r.x[nc:].sum())
        rmax = -np.inf; newc = []
        for pi, m in enumerate(pats):
            th, g = adam(fe, m, y, B, steps, rng, warm=WARM.get(pi))
            o = -simple(m); red = o - y0 - g
            order = np.argsort(-red)[:24]; WARM[pi] = th[np.argsort(-red)[:B // 2]].copy()
            for j in order:
                if red[j] > 1e-9: newc.append((m, th[j].copy(), o))
            rmax = max(rmax, float(red[order[0]]))
        print(f"it={it:3d} minSimple={-val:.6f} slack={slack:.2e} maxred={rmax:.3e} lb~{-(val+max(rmax,0)):.6f} cols={nc} t={time.time()-t0:.0f}s", flush=True)
        if rmax < TOL and slack < 1e-9: break
        if time.time() - t0 > WALL: print(f"wall-clock stop ({WALL} s)", flush=True); break
        if not newc: break
        Vn = np.array([fe.values(th[None], m)[0] for m, th, o in newc]).T
        A = np.hstack([A, Vn]); obj = np.concatenate([obj, [o for _, _, o in newc]])
        if A.shape[1] > PRUNE:   # npafals: prune threshold (was 6000)
            kp = np.unique(np.concatenate([np.where(r.x[:nc] > 1e-12)[0], np.arange(nc, A.shape[1])]))
            A = A[:, kp]; obj = obj[kp]
    if len(sys.argv) > 7:   # npaframe: save duals, features, the law (weights, patterns, positions of the support)
        np.savez(sys.argv[7], y=y, y0=y0, feats=np.array([str((l, m)) for l, m, _, _ in fe.F]), t=fe.t, sc=fe.sc,
                 val=val, rmax=rmax)
    print(f"FINAL N={N} S={S} minSimple={-val:.6f} lb~{-(val+max(rmax,0)):.6f} slack={slack:.1e} features={nf} "
          f"maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.0f}MB t={time.time()-t0:.0f}s", flush=True)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Column generation, continuous positions, for the finite-N CUE multiplicity LP (explore50).
max E[sum m^2]/N over laws of N-multisets on the circle with Haar (Diaconis-Shahshahani) moments for a constraint set:
 P (E|p_k|^2 = k, k<=N), P3 (+ E p_a p_b conj p_{a+b} = 0, a+b<=N), L4 (+ all len(lam)+len(mu)=4, |lam|<=N).
Pricing: batch Adam over positions for each multiplicity pattern (s singles, d doubles[, t triples]).
Bounds: master value (feasible law => lower bound on the max); master + max reduced cost found (heuristic upper bound).
usage: colgen.py N SET [iters] [B] [seed]"""
import sys, math, time
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

def featlist(N, S):
    F = []  # (ups, downs, target, scale)
    for n in range(1, N+1):
        P = list(parts(n))
        for i, lam in enumerate(P):
            for mu in P[i:]:
                L = len(lam) + len(mu)
                if S == 'P' and L > 2: continue
                if S == 'P3' and L > 3: continue
                if S == 'L4' and L > 4: continue
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
    def pvals(self, th, m):   # th (B,n), m (n,) -> Pext (B, N+1) with Pext[:,0]=1
        k = np.arange(1, self.N+1)
        E = np.exp(1j*th[:, :, None]*k[None, None, :])          # (B,n,K)
        P = np.einsum('i,bik->bk', m, E)
        return np.concatenate([np.ones((th.shape[0], 1)), P], axis=1), E
    def q(self, Pe):
        Q = Pe[:, self.idx]                                      # (B,F,W)
        return np.where(self.cj[None], np.conj(Q), Q)
    def values(self, th, m):
        Pe, _ = self.pvals(th, m); Q = self.q(Pe)
        return (np.prod(Q, axis=2).real) / self.sc[None]
    def g_grad(self, th, m, y):   # g = sum_f y_f feat_f ; returns g (B,), dg/dth (B,n)
        Pe, E = self.pvals(th, m); Q = self.q(Pe); B = th.shape[0]
        g = (np.prod(Q, axis=2).real / self.sc[None]) @ y
        G = np.zeros((B, self.N+1), dtype=complex); coef = y / self.sc
        for s in range(self.W):
            Ls = np.prod(np.delete(Q, s, axis=2), axis=2)       # (B,F)
            c = 0.5*coef[None]*np.where(self.cj[None, :, s], np.conj(Ls), Ls)
            np.add.at(G.T, self.idx[:, s], c.T)
        G = G[:, 1:]; k = np.arange(1, self.N+1)
        dth = -2.0*m[None, :]*np.einsum('bk,bik->bi', G*k[None, :], E).imag
        return g, dth

def adam(fe, m, y, B, steps, rng, th=None, lr=0.05):
    n = len(m)
    if th is None: th = rng.uniform(0, 2*np.pi, (B, n))
    th[:, 0] = 0.0
    mo = np.zeros_like(th); v = np.zeros_like(th)
    for t in range(1, steps+1):
        g, d = fe.g_grad(th, m, y); d[:, 0] = 0.0
        mo = 0.9*mo + 0.1*d; v = 0.999*v + 0.001*d*d
        th = th - lr*(mo/(1-0.9**t))/(np.sqrt(v/(1-0.999**t)) + 1e-12)
        if t == steps//2: lr *= 0.2
    g, _ = fe.g_grad(th, m, y)
    return th, g

def patterns(N, triples):
    out = []
    for d in range(0, N//2 + 1):
        out.append([2]*d + [1]*(N-2*d))
    if triples:
        for d in range(0, (N-3)//2 + 1):
            out.append([3] + [2]*d + [1]*(N-3-2*d))
    return [np.array(p, dtype=float) for p in out]

def main():
    N = int(sys.argv[1]); S = sys.argv[2]
    iters = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    B = int(sys.argv[4]) if len(sys.argv) > 4 else 512
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    rng = np.random.default_rng(seed); fe = Feats(N, S); nf = len(fe.F)
    pats = patterns(N, triples=(N <= 8))
    print(f"N={N} S={S} features={nf} patterns={len(pats)} B={B}", flush=True)
    cols, objs, desc = [], [], []
    for m in pats:                                   # initial random columns
        th = rng.uniform(0, 2*np.pi, (200, len(m))); th[:, 0] = 0
        V = fe.values(th, m); cols.append(V); objs += [float((m**2).sum())/N]*len(V); desc += [(m, x) for x in th]
    A = np.vstack(cols).T; obj = np.array(objs); BIG = 20.0; t0 = time.time()
    best_hist = []
    for it in range(iters):
        nc = A.shape[1]
        Aeq = np.vstack([np.hstack([A, np.eye(nf), -np.eye(nf)]), np.hstack([np.ones((1, nc)), np.zeros((1, 2*nf))])])
        beq = np.concatenate([fe.t, [1.0]]); c = np.concatenate([-obj, BIG*np.ones(2*nf)])
        r = linprog(c, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
        lam = r.eqlin.marginals; y = -lam[:nf]; y0 = -lam[nf]
        val = float(obj @ r.x[:nc]); slack = float(r.x[nc:].sum())
        # pricing
        rmax = -np.inf; newc = []
        for m in pats:
            th, g = adam(fe, m, y, B, 300, rng)
            o = float((m**2).sum())/N; red = o - y0 - g
            order = np.argsort(-red)[:12]
            for j in order:
                if red[j] > 1e-9: newc.append((m, th[j].copy(), o))
            rmax = max(rmax, float(red[order[0]]))
        best_hist.append((val, rmax))
        print(f"it={it:3d} val={val:.6f} slack={slack:.2e} maxred={rmax:.3e} ub~{val+max(rmax,0):.6f} cols={nc} t={time.time()-t0:.0f}s", flush=True)
        if rmax < 1e-7 and slack < 1e-9: break
        if not newc: break
        Vn = np.array([fe.values(th[None], m)[0] for m, th, o in newc]).T
        A = np.hstack([A, Vn]); obj = np.concatenate([obj, [o for _, _, o in newc]]); desc += [(m, th) for m, th, o in newc]
        if A.shape[1] > 6000:                        # drop non-basic old columns
            keep = np.where(r.x[:nc] > 1e-12)[0]
            keep = np.unique(np.concatenate([keep, np.arange(nc, A.shape[1])]))
            A = A[:, keep]; obj = obj[keep]; desc = [desc[i] for i in keep]
    # report support of final master
    nc = A.shape[1]
    Aeq = np.vstack([np.hstack([A, np.eye(nf), -np.eye(nf)]), np.hstack([np.ones((1, nc)), np.zeros((1, 2*nf))])])
    r = linprog(np.concatenate([-obj, BIG*np.ones(2*nf)]), A_eq=Aeq, b_eq=np.concatenate([fe.t, [1.0]]), bounds=(0, None), method='highs')
    w = r.x[:nc]; sup = np.where(w > 1e-7)[0]
    print(f"FINAL N={N} S={S} val={float(obj@w):.6f} slack={float(r.x[nc:].sum()):.1e} support={len(sup)}")
    ds = sum(w[i] for i in sup if desc[i][0].max() >= 2)
    print(f"  P(config has a multiple point)={ds:.4f}; mass by #doubles: " +
          ", ".join(f"d={d}:{sum(w[i] for i in sup if (desc[i][0]==2).sum()==d):.4f}" for d in range(N//2+1)))

if __name__ == '__main__':
    main()

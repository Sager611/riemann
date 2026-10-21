"""Numerics only (not part of the proof): positive step weight w on a uniform grid of [-1,1] (u = x/L) whose
weighted-Schur ratio for the prime operator C_p is close to its norm.  Output: w_<tag>.txt (one float per cell).
The rigorous bound is computed from this w by bc.py (schur)."""
import numpy as np, math, sys
def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0
def setup(tag):
    L = 1.5 if tag == 'e3' else 0.5*math.log(float(tag))
    mu = math.exp(2*L)
    ns = [n for n in range(2, int(math.ceil(mu)) + 1) if n < mu and vm(n) > 0]
    return L, [(math.log(n)/L, vm(n)/math.sqrt(n)) for n in ns]
def run(tag, M, iters=5000):
    L, sh = setup(tag)
    h = 2.0/M; edges = -1 + h*np.arange(M+1)
    ops = []
    for t, c in sh:
        for sg in (1, -1):
            lo = edges[:-1] + sg*t; hi = edges[1:] + sg*t
            j0 = np.clip(np.floor((lo + 1)/h - 1e-9).astype(int), 0, M-1)
            j1 = np.clip(np.ceil((hi + 1)/h + 1e-9).astype(int) - 1, 0, M-1)
            valid = (hi > -1) & (lo < 1)
            ops.append((c, j0, j1, valid))
    w = np.ones(M); best = (np.inf, w); stall = 0; it = 0
    for it in range(iters):
        T = np.zeros(M)
        for c, j0, j1, valid in ops:
            m = np.maximum(w[j0], w[j1])
            for d in (1, 2):
                jj = np.minimum(j0 + d, j1); m = np.maximum(m, w[jj])
            T += np.where(valid, c*m, 0.0)
        ln = (T/w).max()
        if ln < best[0] - 1e-15: best = (ln, w.copy()); stall = 0
        else: stall = stall + 1 if it else 0
        if it > 300 and stall > 50: break
        w = T/T.max(); w = np.maximum(w, 1e-12)
        w = 0.5*(w + w[::-1])
    return best[1], best[0], it
if __name__ == "__main__":
    tag, M = sys.argv[1], int(sys.argv[2])
    w, lam, it = run(tag, M)
    np.savetxt(f"w_{tag}_{M}.txt", w, fmt="%.17g")
    print(tag, M, "numeric Schur ratio", lam, "iters", it)

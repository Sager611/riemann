
import numpy as np, time, sys
import explore as E
rng = np.random.default_rng(int(sys.argv[1]))
nb = int(sys.argv[2]); M = int(sys.argv[3]); tmax = float(sys.argv[4]); thr = float(sys.argv[5])
t0 = time.time()
mins = []   # list of (S, g)
def sample(M):
    out = []
    while sum(len(o) for o in out) < M:
        g = rng.uniform(0.743, 5.085, size=(4*M,6))
        g = g[g.sum(axis=1) <= 8.8]
        out.append(g)
    return np.concatenate(out)[:M]
tot = 0
for batch in range(nb):
    g0 = sample(M)
    if batch % 2 == 1:
        g0 = rng.uniform(0.0, 2.5, size=(M,6))
    g, S = E.optimize(g0, iters=800, lo=0.0, hi=9.0)
    tot += M
    sel = np.where(S < thr)[0]
    for i in sel:
        gi = g[i]
        found = False
        for m in mins:
            if np.abs(m[1]-gi).max() < 5e-3:
                found = True
                if S[i] < m[0]: m[0] = S[i]; m[1] = gi.copy()
                m[2] += 1
                break
        if not found:
            mins.append([S[i], gi.copy(), 1])
    if time.time()-t0 > tmax: break
mins.sort(key=lambda m: m[0])
print("starts", tot, "time", round(time.time()-t0,1), "distinct minima below", thr, ":", len(mins))
for m in mins[:40]:
    print("%.9f" % m[0], np.round(m[1],4), "sum=%.4f" % m[1].sum(), "hits", m[2])


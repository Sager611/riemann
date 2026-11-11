# explore31 s6: random designs of the Baez-Duarte density at larger N with a fine pool (no pool-coarseness artefact):
# theta = 1 pinned, N-1 dilations Q/k with k uniform in (Q, QN] (theta density prop. 1/theta^2 on [1/N, 1)), Gram on the fly (Vasyunin).
# usage: uv run --no-project --with numpy --with scipy python3 s6_random.py N Q REPS
import sys, math, time
import numpy as np
from scipy.linalg import cho_factor, cho_solve
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
from s3_hybrid import vas  # noqa: E402  (s3_hybrid's main loop is guarded below)
EG = 0.57721566490153286061
N, Q, REPS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(310 + N)
t0 = time.time(); out = []
for r in range(REPS):
    idx = np.concatenate(([Q], np.sort(rng.choice(np.arange(Q + 1, Q * N + 1), size=N - 1, replace=False))))
    n = len(idx); Gs = np.empty((n, n))
    for a in range(n):
        for b in range(a, n):
            Gs[a, b] = Gs[b, a] = vas(int(idx[a]), int(idx[b]))
    t = (np.log(idx / Q) + 1 - EG) / idx
    cf = cho_factor(Gs, lower=True, check_finite=False)
    v = 1 - Q * float(t @ cho_solve(cf, t, check_finite=False))
    out.append(v)
    check(250, f'rep {r}')
    print(f'N={N} Q={Q} rep {r}: d2={v:.6f}  t={time.time()-t0:.0f}s', flush=True)
out = np.array(out)
print(f'N={N} Q={Q}: random+rho1 mean {out.mean():.6f} sd {out.std():.6f} min {out.min():.6f} max {out.max():.6f}; footprint/resident MB {mem()}')

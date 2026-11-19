from fintest_lib import *
t0 = time.time()
K0 = kern(np.sqrt(2))
for xi in (0.3, 0.2+0.7j, 1.1+2.3j, 2.5j):
    print('K(%s): closed %s quad %s' % (xi, np.round(K0(xi), 10), np.round(quadK(np.sqrt(2), xi), 10)))
def gaps(n, tau):
    g = []
    for _ in range(n):
        r = rng.random()
        if r < 0.25: g.append(rng.uniform(0, 2*tau))
        elif r < 0.4: g.append(rng.uniform(tau, 3*tau))
        elif r < 0.7: g.append(rng.uniform(0.7, 1.3))
        else: g.append(rng.exponential(1.0) + 1e-6)
    return np.cumsum(g)
def rand_config(tau, lattice=False):
    n = int(rng.integers(1, 41))
    xs = (np.arange(n)*rng.uniform(0.8, 1.4) + rng.normal(0, 0.02, n)) if lattice else gaps(n, tau)
    pts = []
    for x in xs:
        m = 1 if rng.random() > 0.08 else int(rng.integers(2, 4))
        pts.append((complex(x, 0), m))
    for _ in range(int(rng.integers(0, 7))):
        x0 = xs[int(rng.integers(0, n))]
        r = rng.random()
        d = tau*(1 if r < 0.2 else (-1 if r < 0.4 else rng.uniform(-1, 1))) if rng.random() < 0.85 else rng.uniform(-2, 2)
        y = float(np.exp(rng.uniform(np.log(0.005), np.log(4.0))))
        pts.append((complex(x0 + d, y), 1 if rng.random() < 0.85 else 2))
    return pts
worst = {}
cases = [(c, w2) for c in (0.001, 0.00476, 0.01, 0.05, 0.1, 0.15, 1/6, 0.2, 0.25) for w2 in (2.0, 3.75)] + [(0.05, 1.5), (0.25, 1.5)]
ntot = 0
for c, w2 in cases:
    K, lam, H = setup(c, np.sqrt(w2))
    wF0 = wF = np.inf
    for it in range(400):
        pts = rand_config(c, lattice=(it % 3 == 0))
        F0, F, N = check(pts, c, K, lam, H); ntot += 1
        # local search: perturb to decrease min(F0, F)
        best = min(F0, F)
        for _ in range(15):
            q = list(pts); j = int(rng.integers(0, len(q))); z, m = q[j]
            z2 = complex(z.real + rng.normal(0, 0.05), max(z.imag + rng.normal(0, 0.1), 0.001) if z.imag > 0 else 0.0)
            q[j] = (z2, m)
            a, b_, _ = check(q, c, K, lam, H); ntot += 1
            if min(a, b_) < best: pts, best, F0, F = q, min(a, b_), a, b_
        wF0, wF = min(wF0, F0), min(wF, F)
    worst[(c, w2)] = (wF0, wF, lam)
    print('c=%.5f w^2=%.2f lam=%.3e: min margin (F0) %.6f, (F) %.6f' % (c, w2, lam, wF0, wF)); sys.stdout.flush()
viol = [k for k, v in worst.items() if min(v[0], v[1]) < -1e-9]
print('multisets checked: %d; violations: %s; elapsed %.1f s' % (ntot, viol if viol else 'none', time.time() - t0))

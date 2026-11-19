# fintest2.py: as fintest.py, restricted to non-trivial multisets (N >= 3 and at least one isolated point with only complex witnesses),
# plus a targeted family: lattice of spacing s with one point x0 having a complex pair at x0 + delta + i y, y near rho_max. Margins printed in %.3e.
from fintest_lib import *
import fintest_lib as FL
t0 = time.time()
def info(pts, tau):
    Z = []
    for z, m in pts:
        Z.append((z, m))
        if z.imag > 0: Z.append((z.conjugate(), m))
    R1 = [z.real for z, m in Z if z.imag == 0 and m == 1]
    R2 = [z.real for z, m in Z if z.imag == 0 and m > 1]
    S = [z for z, m in Z if z.imag > 0]
    I = [x for k, x in enumerate(R1) if all(abs(x - y) > 2*tau for j, y in enumerate(R1) if j != k)]
    n2s = sum(1 for x in I if not any(abs(x - y) <= tau for y in R2) and any(abs(x - z.real) <= tau for z in S))
    return n2s
rng2 = np.random.default_rng(7)
viol = []
for c, w2 in ((0.001, 2.0), (0.00476, 2.0), (0.05, 2.0), (0.1, 2.0), (1/6, 2.0), (1/6, 3.75), (0.25, 2.0)):
    K, lam, H = setup(c, np.sqrt(w2))
    # rho_max: first y with H(y) >= 3 lam
    ys = np.linspace(0.01, 8, 8000); hv = np.array([H(y) for y in ys]) - 3*lam
    rmax = ys[np.argmax(hv >= 0)] if (hv >= 0).any() else 8.0
    best = (np.inf, None)
    for it in range(1500):
        n = int(rng2.integers(3, 30)); s = rng2.uniform(0.6, 1.6)
        xs = np.arange(n)*s + rng2.normal(0, 0.01, n)
        pts = [(complex(x, 0), 1) for x in xs]
        for _ in range(int(rng2.integers(1, 4))):
            x0 = xs[int(rng2.integers(0, n))]; d = rng2.uniform(-c, c)
            y = float(rmax*np.exp(rng2.normal(0, 0.3))) if rng2.random() < 0.7 else float(np.exp(rng2.uniform(np.log(0.01), np.log(4))))
            pts.append((complex(x0 + d, y), 1))
        if rng2.random() < 0.3:   # add a cluster pair and a double point
            x1 = xs[-1] + 2.0; pts += [(complex(x1, 0), 1), (complex(x1 + rng2.uniform(0, 2*c), 0), 1), (complex(x1 + 3.0, 0), 2)]
        if info(pts, c) == 0: continue
        F0, F, N = check(pts, c, K, lam, H)
        m = min(F0, F)
        if m < best[0]: best = (m, N)
        if m < -1e-9: viol.append((c, w2, m, N))
    print('c=%.5f w^2=%.2f lam=%.3e rho_max~%.3f: min margin over non-trivial multisets %.3e (N=%s)' % (c, w2, lam, rmax, best[0], best[1])); sys.stdout.flush()
print('violations:', viol if viol else 'none', '; elapsed %.1f s' % (time.time() - t0))

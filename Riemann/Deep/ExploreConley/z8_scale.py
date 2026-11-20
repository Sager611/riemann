# z8_scale.py — explore66 bounded-scale run (P4, P5) on the first 20000 zeros of zeta (Arb zeta_zeros, rigorous isolation).
# For each consecutive pair k: g_k, u_k = (g/2)^2 (z = 2 gamma units), far tidal constant S_k = sum_{j not in pair, both signs} 1/(c_k - x_j)^2
# (exact over the 20000 listed zeros, density tail beyond), forward rate gdot_k = 4/g - 2 g Sigma_k (Sigma_k = sum 1/((x_{k+1}-x_j)(x_k-x_j))),
# P5 exit time in the frozen far field t_c(k) = -(1/(4S)) log(1/(1 - 2 S u)) when 2 S u < 1 (else the frozen model never collides).
import time, os, sys, numpy as np
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
f = LOG + 'zeros_1_%d.npy' % N
T0 = time.time()
if not os.path.exists(f):
    from flint import acb, ctx
    ctx.prec = 64
    zs = []
    for n0 in range(1, N + 1, 1000):
        zz = acb.zeta_zeros(n0, min(1000, N + 1 - n0))
        zs += [float(z.imag.mid()) for z in zz]
        print('zeros to', n0 + len(zz) - 1, round(time.time() - T0, 1), 's', flush=True)
    np.save(f, np.array(zs))
g = np.load(f); x = 2*g
print('N =', len(x), ' max gamma =', g[-1], flush=True)
W = 600
rho = lambda X: np.log(X/(4*np.pi))/(4*np.pi)
res = []
for k in range(W, len(x) - W - 2):
    a, b = x[k], x[k+1]; c = (a + b)/2; gk = b - a; u = (gk/2)**2
    lo, hi = k - W, k + W + 2
    near = np.concatenate([x[lo:k], x[k+2:hi]])
    S = np.sum(1/(c - near)**2) + np.sum(1/(c + x[:hi])**2)     # listed zeros near + all mirrored -x_j up to hi
    Sig = np.sum(1/((b - near)*(a - near))) + np.sum(1/((b + x[:hi])*(a + x[:hi])))
    # zeros below lo (positive side) and above hi: density tails  int rho(y)/(c-y)^2 dy  ~ rho/dist
    S += rho(x[lo])/(c - x[lo]) + rho(x[hi])/(x[hi] - c) if lo > 0 else rho(x[hi])/(x[hi] - c)
    Sig += rho(x[lo])/(c - x[lo]) + rho(x[hi])/(x[hi] - c) if lo > 0 else rho(x[hi])/(x[hi] - c)
    tc = -np.log(1/(1 - 2*S*u))/(4*S) if 2*S*u < 1 else np.nan
    res.append((k + 1, g[k], gk, u, S, 4/gk - 2*gk*Sig, tc, gk*rho(c)))
r = np.array(res)
np.save(LOG + 'z8_scale.npy', r)
ok = ~np.isnan(r[:, 6])
print('pairs analysed: %d (zeros %d..%d); frozen-model collisions (2Su < 1): %d = %.2f%%' % (len(r), int(r[0, 0]), int(r[-1, 0]) + 1, ok.sum(), 100*ok.mean()))
order = np.argsort(-np.where(ok, r[:, 6], -np.inf))
print('largest exit times t_c(k) (closest to 0):')
for i in order[:8]:
    print('  zeros %d-%d  gamma=%.4f  normalised gap=%.5f  u=%.3e  S=%.4f  t_c=%.6e' % (r[i, 0], r[i, 0] + 1, r[i, 1], r[i, 7], r[i, 3], r[i, 4], r[i, 6]))
for frac in (0.25, 0.5, 1.0):
    m = int(frac*len(r)); sub = r[:m]; okm = ~np.isnan(sub[:, 6])
    print('first %5.0f%% of pairs (gamma <= %.0f): max t_c = %.4e; min normalised gap = %.5f' % (100*frac, sub[-1, 1], np.nanmax(sub[:, 6]), sub[:, 7].min()))
km = np.argmin(r[:, 2]); print('minimal z-gap pair: zeros %d-%d, g=%.5f, gdot=%.4f (P4: >= 0 expected only for a global minimum)' % (r[km, 0], r[km, 0] + 1, r[km, 2], r[km, 5]))
print('#gdot<0: %d of %d; time %.1f s' % ((r[:, 5] < 0).sum(), len(r), time.time() - T0))

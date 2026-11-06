# t2_gridcode.py (explore21): the Euler-product phases t*log p mod 2pi read as a grid code (module p, period 2pi/log p).
# a) capacity: first alias tau_N(delta) > 0 with ||tau log p/2pi|| <= delta for the first N primes (exact interval sieve);
#    grid-code theory (Dirichlet/Kronecker) predicts log tau_N ~ (N-1) log(1/(2 delta)) + C.
# b) random-torus (Bohr-Jessen) model vs t-averages: small-value probabilities P(|F(sigma+it)| < eps), F = zeta and the
#    Davenport-Heilbronn function DH = sum a(n) n^-s, a = (1, kappa, -kappa, -1, 0) mod 5 (no Euler product).
# Usage: t2_gridcode.py a | b sigma T0 T1 step
import sys, math, time
import numpy as np

def primes(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

def first_alias(N, delta, kmax=4*10**7, chunk=2*10**6):
    P = primes(200)[:N]; lam = 2*math.pi/np.log(P.astype(float))
    for k0 in range(1, kmax, chunk):
        k = np.arange(k0, min(k0+chunk, kmax), dtype=float)
        a = (k - delta)*lam[0]; b = (k + delta)*lam[0]
        for j in range(1, N):
            L = lam[j]; na, nb = [], []
            m0 = np.floor(a/L)
            for dm in (0, 1, 2):
                m = m0 + dm
                lo = np.maximum(a, (m - delta)*L); hi = np.minimum(b, (m + delta)*L)
                ok = lo <= hi
                na.append(lo[ok]); nb.append(hi[ok])
            a = np.concatenate(na); b = np.concatenate(nb)
            o = np.argsort(a); a, b = a[o], b[o]
            if a.size == 0: break
        if a.size: return a[0], b[0]
    return None

def part_a():
    t0 = time.time()
    print('N  p_N   delta  first alias tau_N   log(tau_N)  (N-1)log(1/2delta)  diff')
    for delta in (0.1, 0.05):
        for N in range(2, 13):
            r = first_alias(N, delta)
            if r is None: print(N, 'none below kmax'); break
            tau = 0.5*(r[0] + r[1]); P = primes(200)[:N]
            print('%2d %4d  %.2f  %16.3f  %9.3f  %9.3f  %7.3f' % (N, P[-1], delta, tau, math.log(tau), (N-1)*math.log(1/(2*delta)), math.log(tau) - (N-1)*math.log(1/(2*delta))), flush=True)
            if time.time() - t0 > 600: break

def part_b(sigma, T0, T1, step):
    from flint import acb, arb, ctx
    ctx.prec = 64
    kap = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
    A = [1.0, kap, -kap, -1.0]
    t0 = time.time()
    ts = np.arange(T0, T1, step)
    zv = np.empty(ts.size); dv = np.empty(ts.size)
    five = acb(5)
    for i, t in enumerate(ts):
        s = acb(sigma, t)
        zz = s.zeta(); zv[i] = abs(complex(float(zz.real.mid()), float(zz.imag.mid())))
        tot = acb(0)
        for r in range(1, 5): tot += A[r-1]*s.zeta(acb(r)/5)
        tot *= five**(-s)
        dv[i] = abs(complex(float(tot.real.mid()), float(tot.imag.mid())))
    # sanity: DH at sigma=2 vs direct Dirichlet sum
    s = acb(2, 3); tot = acb(0)
    for r in range(1, 5): tot += A[r-1]*s.zeta(acb(r)/5)
    tot *= five**(-s)
    n = np.arange(1, 200001); an = np.array([0.0, 1.0, kap, -kap, -1.0])[n % 5]
    direct = (an*n**(-2.0 - 3j)).sum()
    print('DH check s=2+3i: hurwitz %.10f%+.10fi  direct %.10f%+.10fi' % (float(tot.real.mid()), float(tot.imag.mid()), direct.real, direct.imag))
    print('t-average sigma=%.2f t in [%g,%g) step %g, %d points, %.0fs' % (sigma, T0, T1, step, ts.size, time.time()-t0), flush=True)
    # random model: theta_p iid uniform, primes up to 10^5
    P = primes(100000).astype(float); M = 60000; rng = np.random.default_rng(1)
    u = P**(-sigma)
    chi = {1: 1, 2: 1j, 3: -1j, 4: -1}
    cp = np.array([chi.get(int(p) % 5, 0) for p in P])
    rz = np.empty(M); rd = np.empty(M)
    for c0 in range(0, M, 200):
        th = rng.uniform(0, 2*np.pi, (200, P.size)); e = np.exp(1j*th); del th
        lz = -np.log(1 - u*e).sum(1)
        l1 = -np.log(1 - cp*u*e).sum(1); l2 = -np.log(1 - np.conj(cp)*u*e).sum(1)
        rz[c0:c0+200] = np.abs(np.exp(lz))
        rd[c0:c0+200] = np.abs(0.5*(1 - 1j*kap)*np.exp(l1) + 0.5*(1 + 1j*kap)*np.exp(l2))
    print('random model: %d samples, primes < 1e5, %.0fs' % (M, time.time()-t0))
    print('eps     P(|zeta|<eps) t-avg  random   |  P(|DH|<eps) t-avg  random   | DH ratio to eps^2')
    for eps in (1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01):
        a, b, c, d = (zv < eps).mean(), (rz < eps).mean(), (dv < eps).mean(), (rd < eps).mean()
        print('%-6g  %.2e        %.2e  |  %.2e        %.2e  | %.3f %.3f' % (eps, a, b, c, d, c/eps**2, d/eps**2))
    print('min |zeta| on grid %.4f, min |DH| on grid %.5f; quantiles log|zeta| t-avg sd %.3f random sd %.3f' % (zv.min(), dv.min(), np.log(zv).std(), np.log(rz).std()))

if __name__ == '__main__':
    if sys.argv[1] == 'a': part_a()
    else: part_b(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))

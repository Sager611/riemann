# prolate.py: exact rungs delta_n(lambda) = 1 - |chi_n(lambda)|, n even, of ACLG 2611.13129 (App. A, eqs. tridiag l. 1411-1414,
# chicomp l. 1426-1428): Bouwkamp tridiagonal matrix in orthonormal Legendre polynomials (even k), gamma = 2 pi lambda^2,
# chi_n = lambda sqrt2 d_0^(n) / phi_n(0). Double-precision eigh start, then Rayleigh-quotient iteration in mpmath on the
# tridiagonal matrix. Usage: prolate.py c [nmax] [dps]  (c = lambda^2). Prints JSON. Also importable: rungs(c, nmax, dps).
import sys, json, mpmath as mp, numpy as np

def tri(c, dps):
    mp.mp.dps = dps
    lam = mp.sqrt(c); g = 2*mp.pi*c; g2 = g*g
    K = int(2*float(g) + 60); K += K % 2
    ks = list(range(0, K + 1, 2)); m = len(ks)
    dg = [k*(k + 1) + g2*mp.mpf(2*k*k + 2*k - 1)/((2*k - 1)*(2*k + 3)) for k in ks]
    of = [g2*mp.mpf((k + 1)*(k + 2))/((2*k + 3)*mp.sqrt(mp.mpf((2*k + 1)*(2*k + 5)))) for k in ks[:-1]]
    return lam, ks, dg, of

def tsolve(dg, of, mu, x):            # (T - mu) y = x, Thomas algorithm
    m = len(dg); cp = [0]*m; dp = [0]*m
    b = dg[0] - mu; cp[0] = of[0]/b if m > 1 else 0; dp[0] = x[0]/b
    for i in range(1, m):
        b = dg[i] - mu - of[i - 1]*cp[i - 1]
        cp[i] = of[i]/b if i < m - 1 else 0
        dp[i] = (x[i] - of[i - 1]*dp[i - 1])/b
    y = [0]*m; y[-1] = dp[-1]
    for i in range(m - 2, -1, -1): y[i] = dp[i] - cp[i]*y[i + 1]
    return y

def tmul(dg, of, x):
    m = len(dg)
    return [dg[i]*x[i] + (of[i - 1]*x[i - 1] if i > 0 else 0) + (of[i]*x[i + 1] if i < m - 1 else 0) for i in range(m)]

def rungs(c, nmax=None, dps=None, iters=4):
    c = mp.mpf(c)
    if dps is None: dps = int(30 + 1.15*float(4*mp.pi*c/mp.log(10)))
    lam, ks, dg, of = tri(c, dps); m = len(ks)
    Tn = np.diag([float(v) for v in dg]) + np.diag([float(v) for v in of], 1) + np.diag([float(v) for v in of], -1)
    w, V = np.linalg.eigh(Tn)
    P0 = []
    for k in ks:
        r = mp.mpf(1)
        for j in range(2, k + 1, 2): r *= mp.mpf(j - 1)/j
        P0.append(mp.sqrt(mp.mpf(2*k + 1)/2)*(-1)**(k//2)*r)
    nmax = nmax if nmax is not None else int(4*float(c)) + 12
    out = {}
    for idx in range(min(m, nmax//2 + 1)):
        n = 2*idx
        x = [mp.mpf(float(v)) for v in V[:, idx]]
        nr = mp.sqrt(mp.fsum(v*v for v in x)); x = [v/nr for v in x]
        for _ in range(iters):
            mu = mp.fsum(a*b for a, b in zip(x, tmul(dg, of, x)))
            y = tsolve(dg, of, mu, x)
            nr = mp.sqrt(mp.fsum(v*v for v in y)); x = [v/nr for v in y]
        Tx = tmul(dg, of, x); mu = mp.fsum(a*b for a, b in zip(x, Tx))
        res = mp.sqrt(mp.fsum((a - mu*b)**2 for a, b in zip(Tx, x)))
        phi0 = mp.fsum(x[r]*P0[r] for r in range(m))
        chi = lam*mp.sqrt(2)*x[0]/phi0
        delta = 1 - abs(chi)
        out[n] = {'chi_sign': int(mp.sign(chi)), 'delta': delta, 'res': res, 'tail': abs(x[-1])}
    return out, dps, ks[-1]

if __name__ == '__main__':
    c = sys.argv[1]; nmax = int(sys.argv[2]) if len(sys.argv) > 2 else None
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else None
    out, dps, K = rungs(mp.mpf(c), nmax, dps)
    print(json.dumps({'c': c, 'dps': dps, 'K': K, 'rungs': {n: {'chi_sign': v['chi_sign'], 'delta': mp.nstr(v['delta'], 20),
          'mlog': mp.nstr(-mp.log(v['delta']), 12) if v['delta'] > 0 else None, 'res': mp.nstr(v['res'], 3), 'tail': mp.nstr(v['tail'], 3)}
          for n, v in out.items()}}))

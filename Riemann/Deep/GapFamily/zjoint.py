"""gapfamily: order-J criterion (ExploreShape Criterion J) for Bui-Hall's f_k, amplifier P of degree d (BH: d = 1).
Joint loop as ExploreShape joint.run: Wirtinger optimum (V, P) -> SDP with cutting planes (sdpcut.sdp_cut) for the form ->
re-optimise V and P for that form -> repeat.  h_J > 1 means 'criterion J proves Lambda^(k) > lam' (float).
usage: python3 zjoint.py k lam d J1,J2,.. [nit] [eta|mid]"""
import sys, json, time, os
from math import comb
import numpy as np, scipy.linalg as sl
from scipy.optimize import minimize_scalar
import zk
from sdpcut import sdp_cut
HERE = os.path.dirname(os.path.abspath(__file__))


def lmax(A, B):
    ev, U = sl.eigh(A, B)
    return ev[-1], U[:, -1]


def shifted(M, t, lam, n):
    """lam^n * form of (S - (v0 + t))^n from forms of (S - v0)^i"""
    return lam ** n * sum(comb(n, i) * M[i] * (-t) ** (n - i) for i in range(n + 1))


def run(k, lam, d, Js, nit=3, eta='mid', nodes=None, v0=1.237, verbose=True):
    eta = lam / 4 if eta == 'mid' else float(eta)
    nodes = nodes or (8, 8, 6, 8, 6, 12 + k, 12 + k)
    t0 = time.time()
    M, herm, imag = zk.forms(k, lam, eta, v0, d, 12, nodes=nodes)
    out = {'k': k, 'lam': lam, 'd': d, 'eta': eta, 'nodes': nodes, 'herm': herm, 'imag': imag, 'tforms': time.time() - t0}
    fw = lambda t: -lmax(M[0], shifted(M, t, lam, 2))[0]
    r = minimize_scalar(fw, bounds=(-0.15, 0.15), method='bounded', options={'xatol': 1e-9})
    tw = r.x; out['h1'] = -r.fun; out['V1'] = v0 + tw
    _, y1 = lmax(M[0], shifted(M, tw, lam, 2))
    mu = [y1 @ shifted(M, tw, 1.0, n) @ y1 / (y1 @ M[0] @ y1) for n in range(13)]
    out['mu'] = mu; out['kurt'] = mu[4] / mu[2] ** 2; out['m6s'] = mu[6] / mu[2] ** 3
    if verbose:
        print('   moments about V (Wirtinger P): mu1 %.3g, mu2 %.6g, kurtosis %.5f, m6/m2^3 %.4f, m8/m2^4 %.4f' % (mu[1], mu[2], out['kurt'], out['m6s'], mu[8] / mu[2] ** 4), flush=True)
    if verbose:
        print('k=%d lam=%.4f d=%d: h1 = %.6f (v = %.5f)  [forms %.0fs, herm %.1e imag %.1e]' % (k, lam, d, out['h1'], v0 + tw - 1, out['tforms'], herm, imag), flush=True)
    for J in Js:
        t = tw; _, y = lmax(M[0], shifted(M, t, lam, 2))
        hJ = None
        for it in range(nit):
            n0 = y @ M[0] @ y
            m = [y @ shifted(M, t, lam, 2 * j) @ y / n0 for j in range(J + 1)]
            val, a, ph, mnv, rds = sdp_cut(m, J, rounds=10)
            def fa(tt):
                R = sum(a[j] * shifted(M, tt, lam, 2 * j) for j in range(1, J + 1))
                return a[0] * lmax(M[0], R)[0]
            rr = minimize_scalar(fa, bounds=(t - 0.01, t + 0.01), method='bounded', options={'xatol': 1e-9})
            t = rr.x; hJ = -rr.fun
            _, y = lmax(M[0], sum(a[j] * shifted(M, t, lam, 2 * j) for j in range(1, J + 1)))
            if verbose:
                print('   J=%d it %d: SDP %.6g (cut rounds %d, minvalid %.2g) -> h_J = %.6f (v = %.5f)' % (J, it, val, rds, mnv, hJ, v0 + t - 1), flush=True)
        out['h%d' % J] = hJ; out['a%d' % J] = [float(x) for x in a]; out['phi%d' % J] = {str(kk): vv for kk, vv in ph.items()}
        out['V%d' % J] = v0 + t; out['P%d' % J] = [float(x) for x in y / y[0]]; out['minvalid%d' % J] = mnv
    return out


if __name__ == '__main__':
    k = int(sys.argv[1]); lam = float(sys.argv[2]); d = int(sys.argv[3]); Js = [int(x) for x in sys.argv[4].split(',')]
    nit = int(sys.argv[5]) if len(sys.argv) > 5 else 3; eta = sys.argv[6] if len(sys.argv) > 6 else 'mid'
    out = run(k, lam, d, Js, nit, eta)
    fn = os.path.join(HERE, 'logs', 'zjoint.jsonl')
    with open(fn, 'a') as f:
        f.write(json.dumps(out) + '\n')

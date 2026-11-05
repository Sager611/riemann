# t1b_branches.py -- explore3 transfer 1, diagnostics:
# (a) Eisenstein matrix element F(b) = <(1-L_b)^{-1} 1, delta_0> versus 2 zeta(2b-1)/zeta(2b);
# (b) eigenvalue branches through 1 at b = rho/2 and through +-1 at cusp forms, along horizontal lines;
# (c) the branch through rho_1/2 followed along Re b = 1/4.
import sys, time
import numpy as np
from flint import acb, acb_mat, ctx
sys.path.insert(0, '.')
from t1_flint import build_M

def build_acb(b, K, prec):
    ctx.prec = prec
    b = acb(b.real, b.imag)
    two = acb(2)
    Z = [(2*b + m).zeta(two) for m in range(2*K - 1)]
    H = []
    for j in range(K):
        a = 2*b + j
        p = acb(1); row = []
        for n in range(K):
            row.append(p * Z[j + n]); p = p * (a + n) / (n + 1)
        H.append(row)
    M = acb_mat(K, K)
    for n in range(K):
        d = [H[j][n] for j in range(K)]
        for k in range(K):
            M[n, k] = d[0] if (n + k) % 2 == 0 else -d[0]
            d = [d[j] - d[j + 1] for j in range(len(d) - 1)]
    return M

def F_eis(b, K=90, prec=400):
    M = build_acb(b, K, prec)
    I = acb_mat(K, K)
    for i in range(K): I[i, i] = acb(1)
    rhs = acb_mat(K, 1)
    rhs[0, 0] = acb(1)
    x = (I - M).solve(rhs)
    s = acb(0)
    for k in range(K):
        s = s + x[k, 0] if k % 2 == 0 else s - x[k, 0]
    ctx.prec = prec
    bb = acb(b.real, b.imag)
    ref = 2 * (2*bb - 1).zeta() / (2*bb).zeta()
    return complex(s), complex(ref)

def track(points, K, start, target):
    # follow the eigenvalue closest to 'target' at points[start] along the list in both directions
    ev = {}
    for i, b in enumerate(points):
        ev[i] = np.linalg.eigvals(build_M(b, K, prec=64 + 3*K))
    lam = [None]*len(points)
    lam[start] = ev[start][np.argmin(abs(ev[start] - target))]
    for i in range(start + 1, len(points)):
        lam[i] = ev[i][np.argmin(abs(ev[i] - lam[i-1]))]
    for i in range(start - 1, -1, -1):
        lam[i] = ev[i][np.argmin(abs(ev[i] - lam[i+1]))]
    return lam, ev

if __name__ == '__main__':
    t0 = time.time()
    g = [14.134725141734693790, 21.022039638771554993, 25.010857580145688763, 30.424876125859513210]
    print('(a) Eisenstein matrix element F(b) vs 2zeta(2b-1)/zeta(2b), K=90')
    tests = [complex(1.5, 2), complex(0.8, 5), complex(0.6, 9), complex(0.3, 5), complex(0.25, 9),
             complex(0.5, 13.779751351890738 + 1e-6), complex(0.5, 9.533695261353557), complex(0.5, 17.73856338),
             complex(0.25, g[0]/2 + 0.01), complex(0.25, g[0]/2 + 0.001)]
    for b in tests:
        f, ref = F_eis(b)
        print('  b=%-22s F=%-40s ref=%-40s |diff|=%.2e' % (b, '%.10f%+.10fi' % (f.real, f.imag), '%.10f%+.10fi' % (ref.real, ref.imag), abs(f - ref)))
    print('(b) |lambda| along horizontal lines through eigenvalue-one points (K=100)')
    sig = [0.05 + 0.025*i for i in range(29)]
    lines = [('zeta rho1/2', g[0]/2, 0.25, 1.0), ('zeta rho2/2', g[1]/2, 0.25, 1.0), ('zeta rho3/2', g[2]/2, 0.25, 1.0),
             ('even cusp 13.7798', 13.779751351890738, 0.5, 1.0), ('odd cusp 9.5337', 9.533695261353557, 0.5, -1.0),
             ('odd cusp 12.1730', 12.173008324679, 0.5, -1.0)]
    for name, t, s0, tgt in lines:
        pts = [complex(s, t) for s in sig]
        start = int(round((s0 - 0.05) / 0.025))
        lam, ev = track(pts, 100, start, tgt)
        print('  %-18s' % name + ' '.join('%.2f:%.3f' % (s, abs(l)) for s, l in zip(sig, lam) if abs(round((s - 0.05)/0.1) - (s-0.05)/0.1) < 1e-9 or abs(s - s0) < 1e-9))
        rr = [abs(e).max() for e in [ev[i] for i in range(len(sig))]]
        print('  %-18s' % ' spectral radius' + ' '.join('%.2f:%.3f' % (s, r) for s, r in zip(sig, rr) if abs(round((s - 0.05)/0.1) - (s-0.05)/0.1) < 1e-9 or abs(s - s0) < 1e-9))
    print('(c) branch through rho1/2 followed up the line Re b = 1/4 (K=100)')
    ts = [g[0]/2 + 0.05*i for i in range(0, 70)]
    pts = [complex(0.25, t) for t in ts]
    lam, ev = track(pts, 100, 0, 1.0)
    for i in range(0, 70, 5):
        print('  t=%.3f lambda=%.4f%+.4fi |lambda|=%.4f  r=%.3f' % (ts[i], lam[i].real, lam[i].imag, abs(lam[i]), abs(ev[i]).max()))
    print('time %.1fs' % (time.time() - t0))

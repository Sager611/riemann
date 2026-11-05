# t1_flint.py -- explore3 transfer 1: Mayer operator via python-flint (Arb), fast version of t1_mayer.py
# Same matrix: M[n,k] = (-1)^(n+k) sum_j C(k,j) (-1)^j (2b+j)_n/n! zeta(2b+j+n, 2)  (basis (z-1)^k, disc |z-1|<3/2)
# Run: uv run --no-project --with python-flint --with numpy python t1_flint.py <mode>
import sys, time
import numpy as np
from flint import acb, arb, ctx

def build_M(b, K, prec=256):
    ctx.prec = prec
    b = acb(b.real, b.imag) if isinstance(b, complex) else acb(b)
    two = acb(2)
    Z = [(2*b + m).zeta(two) for m in range(2*K - 1)]
    M = np.zeros((K, K), dtype=complex)
    # h[j][n] = (2b+j)_n / n! * Z[j+n]
    H = []
    for j in range(K):
        a = 2*b + j
        p = acb(1)
        row = []
        for n in range(K):
            row.append(p * Z[j + n])
            p = p * (a + n) / (n + 1)
        H.append(row)
    for n in range(K):
        d = [H[j][n] for j in range(K)]
        for k in range(K):
            v = d[0]
            c = complex(float(v.real.mid()), float(v.imag.mid()))
            M[n, k] = c if (n + k) % 2 == 0 else -c
            d = [d[j] - d[j + 1] for j in range(len(d) - 1)]
            if not d: break
    return M

def eigs(b, K, prec=256):
    ev = np.linalg.eigvals(build_M(b, K, prec))
    return ev[np.argsort(-abs(ev))]

if __name__ == '__main__':
    mode = sys.argv[1]
    t0 = time.time()
    if mode == 'sanity':
        ev = eigs(complex(1, 0), 40)
        print('b=1 K=40:', ev[:3])
        M = build_M(complex(0.25, 14.134725141734693790/2), 50)
        print('det(1-M) at rho1/2, K=50:', abs(np.linalg.det(np.eye(50) - M)), ' det(1+M):', abs(np.linalg.det(np.eye(50) + M)))
    elif mode == 'conv':
        for (sg, tt) in ((1.0, 20.0), (0.5, 20.0), (0.25, 20.0), (1.0, 30.0), (0.5, 30.0), (0.25, 30.0), (0.5, 40.0)):
            row = []
            for K in (60, 80, 100, 130):
                t1 = time.time()
                ev = eigs(complex(sg, tt), K, prec=64 + 3*K)
                row.append('K=%d r=%.8f (%.1fs)' % (K, abs(ev[0]), time.time() - t1))
            print('b=%.2f+%.1fi: ' % (sg, tt) + ' | '.join(row), flush=True)
    print('total %.1fs' % (time.time() - t0))

def Kfor(t):
    return 100 if t <= 20 else (150 if t <= 30 else 190)

if __name__ == '__main__' and sys.argv[1] == 'grid':
    t0 = time.time()
    sigmas = [round(0.20 + 0.05*i, 2) for i in range(19)]
    ts = [0.5*i for i in range(1, 65)]
    out = open('t1_grid.tsv', 'w')
    out.write('sigma\tt\tK\tr1\tr2\tr3\tl1_re\tl1_im\n')
    for t in ts:
        K = Kfor(t)
        for s in sigmas:
            ev = eigs(complex(s, t), K, prec=64 + 3*K)
            out.write('%.2f\t%.1f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n' % (s, t, K, abs(ev[0]), abs(ev[1]), abs(ev[2]), ev[0].real, ev[0].imag))
        out.flush()
    # convergence spot checks
    for (s, t) in ((0.25, 32.0), (0.5, 32.0), (1.0, 32.0), (0.5, 25.0)):
        a = abs(eigs(complex(s, t), Kfor(t), prec=64 + 3*Kfor(t))[0]); b = abs(eigs(complex(s, t), 220, prec=64 + 660)[0])
        print('conv check b=%.2f+%.1fi: K=%d r=%.6f, K=220 r=%.6f' % (s, t, Kfor(t), a, b))
    print('grid done %.1fs' % (time.time() - t0))

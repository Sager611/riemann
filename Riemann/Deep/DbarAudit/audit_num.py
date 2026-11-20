# audit_num.py (dbaraudit, 20 Nov 2026): falsification of Lemma B (endpoint ratio, gamma^2 gain) and Lemma D
# (almost orthogonality of groups) in the rescaled variable x = u/ell in [0,1]; frequencies in units 1/ell, so a gap s0 = L^{-3/4}
# becomes D = s0*ell ~ L^{1/4}.  Window psi = Gevrey-2 bump (e^{-1/y} construction), psi = 1 on [1/4,3/4], chi = psi^2.
# Group spaces: Newton (Hermite) basis of exponentials at nodes (repeated nodes = jets), built in mpmath, then double QR.
import sys, numpy as np, mpmath as mp, time
mp.mp.dps = 40
t0 = time.time()
def f_(y): return np.where(y > 0, np.exp(-1.0 / np.where(y > 0, y, 1.0)), 0.0)
def S(y):
    y = np.clip(y, 0.0, 1.0); return f_(y) / (f_(y) + f_(1.0 - y))
def psi(x): return S(4 * x) * S(4 * (1 - x))
def quad(P):
    """composite Gauss-Legendre: P panels (P multiple of 8, so 1/8, 1/4, 3/4, 7/8 are panel ends), 48 nodes each.
    (leggauss(n) builds an n x n companion matrix: n = 4800 took 520 MB and the first D run was killed by memcap.)"""
    P = int(8 * np.ceil(P / 8)); g, w = np.polynomial.legendre.leggauss(48)
    a = np.arange(P) / P; h = 1.0 / P
    X = (a[:, None] + h / 2 + h / 2 * g[None, :]).ravel(); Wt = np.tile(h / 2 * w, P)
    return X, Wt
def newton(nodes, x, deriv=False):
    """Newton basis nu_k(v) = (e^{i . v})[d_1..d_k] at v = x - 1/2 (and d/dx of it); nodes sorted, repeats adjacent."""
    nodes = sorted(nodes); n = len(nodes); out = np.zeros((len(x), n), complex); dout = np.zeros_like(out)
    for q, xv in enumerate(x):
        v = mp.mpf(xv) - mp.mpf(1) / 2
        for which in ((0, 1) if deriv else (0,)):
            T = []
            for i in range(n):
                d = mp.mpf(nodes[i]); e = mp.expj(d * v)
                T.append(e if which == 0 else 1j * d * e)
            col = [T[0]]
            for j in range(1, n):
                newT = []
                for i in range(n - j):
                    if nodes[i + j] == nodes[i]:
                        d = mp.mpf(nodes[i]); iv = 1j * v
                        if which == 0: val = iv ** j * mp.expj(d * v) / mp.factorial(j)
                        else: val = 1j * (d * iv ** j + j * iv ** (j - 1)) * mp.expj(d * v) / mp.factorial(j)
                    else:
                        val = (T[i + 1] - T[i]) / (mp.mpf(nodes[i + j]) - mp.mpf(nodes[i]))
                    newT.append(val)
                T = newT; col.append(T[0])
            arr = np.array([complex(c) for c in col])
            if which == 0: out[q] = arr
            else: dout[q] = arr
    return (out, dout) if deriv else out
def config(kind, M, D, seed=1):
    if kind == 'jet': return [0.0] * M
    if kind == 'chain': return [j * 0.99 * D for j in range(M)]
    if kind == 'tight': return [j * 0.01 * D for j in range(M)]
    rng = np.random.default_rng(seed); return list(np.concatenate([[0], np.cumsum(rng.uniform(0, D, M - 1))]))
mode = sys.argv[1] if len(sys.argv) > 1 else 'B'
if mode == 'B':
    # N_B = sup |R(0)|^2 / int |R|^2 chi  (gamma*ell -> inf limit of Lemma B's ratio); R_B(gl) with R + R'/(i gl)
    for kind in ('jet', 'tight', 'chain', 'rand'):
        for D in (2.1, 4.2, 10.0, 31.6):
            if kind == 'jet' and D != 2.1: continue
            row = []
            for M in (1, 2, 3, 4, 6, 8):
                vals = []
                for nmid in (16, 32):
                    x, w = quad(nmid); wc = w * psi(x) ** 2
                    nodes = config(kind, M, D)
                    A, dA = newton(nodes, x, deriv=True)
                    b0 = newton(nodes, np.array([0.0]))[0]
                    res = []
                    for gl in (np.inf, 280.0):
                        Bm = A if gl == np.inf else A + dA / (1j * gl)
                        Q, R = np.linalg.qr(np.sqrt(wc)[:, None] * Bm)
                        y = np.linalg.solve(R.conj().T, b0.conj())
                        res.append(float(np.vdot(y, y).real))
                    vals.append(res)
                hc = abs(vals[1][0] - vals[0][0]) / vals[1][0]
                row.append(f'M={M}: N_B={vals[1][0]:.4e} R_B(280)/N_B={vals[1][1]/vals[1][0]:.4f} N_B^(1/M)={vals[1][0]**(1/M):.2f} h={hc:.0e}')
            print(f'B {kind:5s} D={D:5.1f} | ' + ' | '.join(row), flush=True)
if mode == 'D':
    # 2 groups at gap exactly D: cos of the principal angle in L^2(chi); 3 groups: sup sum_g||q_g||^2/||q||^2 = 1/sigma_min^2
    for kind in ('jet', 'chain'):
        for M in (1, 2, 4, 8):
            row = []
            for D in (4.2, 10.0, 31.6, 100.0):
                g1 = config(kind, M, D); width = max(g1) - min(g1)
                om = 3 * (width + D) + 10
                vals = []
                for P in (max(16, om / 20), max(32, om / 10)):
                    x, w = quad(P); sw = np.sqrt(w * psi(x) ** 2)
                    Qs = [np.linalg.qr(sw[:, None] * newton([z + k * (width + D) for z in g1], x))[0] for k in range(3)]
                    cos = np.linalg.svd(Qs[0].conj().T @ Qs[1], compute_uv=False)[0]
                    s3 = np.linalg.svd(np.concatenate(Qs, axis=1), compute_uv=False)[-1]
                    vals.append((cos, 1 / s3 ** 2))
                row.append(f'D={D}: cos={vals[1][0]:.3e} K3={vals[1][1]:.4e} (h {abs(vals[1][0]-vals[0][0]):.0e})')
            print(f'D {kind:5s} M={M} | ' + ' | '.join(row), flush=True)
if mode == 'D5':
    for kind in ('jet', 'chain'):
        for M in (2, 4, 6):
            row = []
            for D in (2.1, 4.2, 10.0, 31.6):
                nmid = int(max(320, 30 * M * D)); x, w = quad(nmid); sw = np.sqrt(w * psi(x) ** 2)
                g1 = config(kind, M, D); width = max(g1) - min(g1)
                Q = np.concatenate([np.linalg.qr(sw[:, None] * newton([z + k * (width + D) for z in g1], x))[0] for k in range(5)], axis=1)
                s = np.linalg.svd(Q, compute_uv=False)
                row.append(f'D={D}: sup sum||q_g||^2/||q||^2 = {1/s[-1]**2:.6f}')
            print(f'D5 {kind:5s} M={M} | ' + ' | '.join(row), flush=True)
print('time', round(time.time() - t0, 1), 's')


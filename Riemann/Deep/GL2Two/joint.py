"""Joint alternating optimisation of the two-piece functional c = c1 + 2 c12 + c2 (float, untrusted) at fixed R.
Variables: P1 = x + x(1-x) sum_j p_j T_j(2x-1) (j < J), g = P2'' = sum_k g_k s T_{k-1}(2s-1) (k <= K),
Q = 1 + sum_n h_n (T_{2n-1}(1-2y) - 1) (n <= MQ).  For fixed Q, c is a quadratic in (p, g): exact solve.
For fixed (p, g), c = 1/2 + v^T H v, v = (1, h): exact solve.  Each step lowers c.
Start: Degree2Mol one-piece optimum (P1 = sinh, Q, R from params.json or onepiece.best)."""
import sys, time, json, os
from fractions import Fraction
import numpy as np
from numpy.polynomial import chebyshev as C
from tp import *
from onepiece import basis as qbasis, best

def cheb_derivs(z, kmax, mmax=4):
    """T_k^{(m)}(z) for k <= kmax, m <= mmax via the differentiated recurrence (z-derivatives)."""
    T = np.zeros((mmax + 1, kmax + 1) + z.shape)
    T[0, 0] = 1; T[0, 1] = z
    if mmax >= 1: T[1, 1] = 1
    for k in range(1, kmax):
        for m in range(mmax + 1):
            T[m, k + 1] = 2 * z * T[m, k] - T[m, k - 1] + (2 * m * T[m - 1, k] if m else 0)
    return T

def rho_derivs(y, MQ, mmax=4):
    """rho_a^{(m)}(y), a = 0..MQ: rho_0 = 1, rho_n = T_{2n-1}(1-2y) - 1 (y-derivatives). shape (mmax+1, MQ+1)+y.shape"""
    T = cheb_derivs(1 - 2 * y, 2 * MQ - 1, mmax)
    out = np.zeros((mmax + 1, MQ + 1) + y.shape)
    out[0, 0] = 1
    for n in range(1, MQ + 1):
        for m in range(mmax + 1):
            out[m, n] = (-2.0) ** m * T[m, 2 * n - 1] - (1 if m == 0 else 0)
    return out

class P1Basis:
    def __init__(self, J): self.J = J
    def ev(self, x, der=0):
        """(J+1, ...) : row 0 = affine part x (or derivatives), rows 1..J = x(1-x)T_{j}(2x-1)."""
        x = np.asarray(x, float); z = 2 * x - 1
        base = [x, np.ones_like(x), np.zeros_like(x)][der]
        rows = [base]
        for j in range(self.J):
            e = np.zeros(j + 1); e[-1] = 1
            T0 = C.chebval(z, e)
            T1 = 2 * C.chebval(z, C.chebder(e)) if j >= 1 else np.zeros_like(z)
            T2 = 4 * C.chebval(z, C.chebder(e, 2)) if j >= 2 else np.zeros_like(z)
            q0 = x * (1 - x); q1 = 1 - 2 * x; q2 = -2.0
            rows.append([q0 * T0, q1 * T0 + q0 * T1, q2 * T0 + 2 * q1 * T1 + q0 * T2][der])
        return np.array(rows)

def fit_P1(Pfun, J):
    x, w = gl01(80); Bp = P1Basis(J); Phi = Bp.ev(x)
    r = Pfun(x) - Phi[0]
    G = (Phi[1:] * w) @ Phi[1:].T; b = (Phi[1:] * w) @ r
    return np.linalg.solve(G, b)

class PEval:
    def __init__(self, Bp, p): self.Bp, self.p = Bp, np.concatenate([[1.0], p])
    def __call__(self, x, k=0): return np.tensordot(self.p, self.Bp.ev(np.asarray(x, float), k), axes=(0, 0))

def H12_Q(Pf, R, th1, th2, Bg, g, MQ, n=40, chunk=8):
    """matrix K with c12 = v^T K v (symmetrised), for fixed P1 (callable) and g; chunked over the first variable."""
    p_, wp = gl01(n); q_, wq = gl01(n); u_, wu = gl01(n)
    K = np.zeros((MQ + 1, MQ + 1))
    for i0 in range(0, n, chunk):
        PP, QQ, UU = [a.ravel() for a in np.meshgrid(p_[i0:i0 + chunk], q_, u_, indexing='ij')]
        W = (wp[i0:i0 + chunk, None, None] * wq[None, :, None] * wu[None, None, :]).ravel() * (1 - PP)
        a = PP; b = (1 - PP) * QQ; s = (1 - PP) * (1 - QQ) * UU
        al = a * UU * th2; be = 1 - b * UU * th2
        ra = rho_derivs(al, MQ, 1); rb = rho_derivs(be, MQ, 1)
        A0 = ra[0]; A1 = -th1 * (R * ra[0] + ra[1]); B0 = rb[0]; B1 = th1 * (R * rb[0] + rb[1])
        z = 1 - (1 - UU) * th2 / th1
        gs = g @ Bg.ev(s, None)
        base = 4 * (th2 / th1) ** 2 * np.exp(R) * W * UU ** 2 * (1 - UU) * np.exp(R * UU * th2 * (a - b)) * gs
        f1 = base * Pf(z); f2 = base * Pf(z, 1); f3 = base * Pf(z, 2)
        K += (A1 * f1) @ B1.T + (A1 * f2) @ B0.T + (A0 * f2) @ B1.T + (A0 * f3) @ B0.T
        memguard()
    return (K + K.T) / 2

def H2_Q(R, th2, Bg, g, MQ, nt=120, nr=24, nu=24, tchunk=30):
    """matrix H with c2 = v^T H v for fixed g."""
    t, wt = gl01(nt); r, wr = gl01(nr); u, wu = gl01(nu)
    H = np.zeros((MQ + 1, MQ + 1))
    fac = [1, 1, 2]
    for ir in range(nr):
        rr = r[ir]
        # X jets in x for (r, u) and Y jets in y for (r, v) with fixed g
        c_ = 1 - u; s0 = c_ * rr
        f0 = g @ Bg.ev(s0, None, 0); f1 = g @ Bg.ev(s0, None, 1); f2 = g @ Bg.ev(s0, None, 2)
        Xj = np.stack([rr * f0, f0 + c_ * rr * f1, (2 * c_ * f1 + c_ * c_ * rr * f2) / 2])     # (3, nu)
        for t0 in range(0, nt, tchunk):
            tt = t[t0:t0 + tchunk]; wtt = wt[t0:t0 + tchunk]
            T, U, V = np.meshgrid(tt, u, u, indexing='ij')
            S0 = -rr * (U + V); kap = th2 * R * (2 * T - 1)
            e0 = np.exp(2 * R * T + kap * S0); cx = kap * (1 - U); cy = kap * (1 - V)
            ex = [np.ones_like(cx), cx, cx * cx / 2]; ey = [np.ones_like(cy), cy, cy * cy / 2]
            E = np.zeros(T.shape + (3, 3))
            for i in range(3):
                for j in range(3): E[..., i, j] = e0 * ex[i] * ey[j]
            Lin = np.zeros(T.shape + (3, 3))
            Lin[..., 0, 0] = 1 / th2 + S0; Lin[..., 1, 0] = 1 - U; Lin[..., 0, 1] = 1 - V
            E = jmul(Lin, E)
            XY = np.zeros(T.shape + (3, 3))
            for i in range(3):
                for j in range(3): XY[..., i, j] = Xj[i][None, :, None] * Xj[j][None, None, :]
            Z = jmul(E, XY)
            wgt = (wtt[:, None, None] * wu[None, :, None] * wu[None, None, :]) * wr[ir] * (1 - rr) ** 4 * (8.0 / 3.0)
            al1 = th2 * U * rr + T * (1 + th2 * S0); be1 = th2 * U + T * th2 * (1 - U); ga1 = -th2 + T * th2 * (1 - V)
            al2 = th2 * V * rr + T * (1 + th2 * S0); be2 = -th2 + T * th2 * (1 - U); ga2 = th2 * V + T * th2 * (1 - V)
            R1 = rho_derivs(al1.ravel(), MQ); R2 = rho_derivs(al2.ravel(), MQ)
            b1, g1, b2, g2 = be1.ravel(), ga1.ravel(), be2.ravel(), ga2.ravel()
            Zf = Z.reshape(-1, 3, 3) * wgt.ravel()[:, None, None]
            for i1 in range(3):
                for j1 in range(3):
                    D1 = R1[i1 + j1] * (b1 ** i1 * g1 ** j1 / (fac[i1] * fac[j1]))
                    for i2 in range(3 - i1):
                        for j2 in range(3 - j1):
                            D2 = R2[i2 + j2] * (b2 ** i2 * g2 ** j2 / (fac[i2] * fac[j2]))
                            H += (D1 * Zf[:, 2 - i1 - i2, 2 - j1 - j2]) @ D2.T
            memguard()
    return (H + H.T) / 2

def c1_parts(Pf, R, th1, MQ):
    G0, G1, end = qbasis(R, MQ)
    x, w = gl01(80); Cp = np.sum(w * Pf(x, 1) ** 2); Bp = np.sum(w * Pf(x) ** 2)
    return G0, G1, end, Cp, Bp

def run(name, R, th1, th2, lam, h, MQ, J=12, K=8, iters=6, res=(120, 24, 24), two=True):
    t0 = time.time()
    Bp = P1Basis(J); Bg = GBasis(K)
    P1s = P_sinh(lam); p = fit_P1(P1s, J); Pf = PEval(Bp, p)
    v = np.concatenate([[1.0], h])
    G0, G1, end = qbasis(R, MQ)
    x, w = gl01(80)
    Phi0, Phi1 = Bp.ev(x, 0), Bp.ev(x, 1)
    GP0 = (Phi0 * w) @ Phi0.T; GP1 = (Phi1 * w) @ Phi1.T          # Gram of (affine, basis) for int P^2, int P'^2
    def c1_of(v, Pf):
        W0 = v @ G0 @ v; W1 = v @ G1 @ v
        Cp = np.sum(w * Pf(x, 1) ** 2); Bq = np.sum(w * Pf(x) ** 2)
        return 0.5 * (1 + (end @ v) ** 2) + (Cp / th1) * W0 + th1 * Bq * W1
    cA = c1_of(v, Pf); print(f"== {name}: R={R} th1={th1:.9f} th2={th2:.9f} MQ={MQ} J={J} K={K}; start c1={cA:.10e} kappa1={1-np.log(cA)/R:.12f}", flush=True)
    g = np.zeros(K); hist = []
    for it in range(iters):
        # (p, g) step for fixed Q
        Q = QCheb.from_h(v[1:])
        W0 = v @ G0 @ v; W1 = v @ G1 @ v
        A = (W0 / th1) * GP1 + th1 * W1 * GP0                       # quadratic form of c1 in (1, p)
        # c12 is linear in P1: L(P1) = sum_i P1coef_i L(pi_i); compute L for each basis element
        Ls = []
        for i in range(J + 1 if two else 0):
            coef = np.zeros(J + 1); coef[i] = 1
            class Pi:
                def __init__(s, c): s.c = c
                def __call__(s, xx, k=0): return np.tensordot(s.c, Bp.ev(np.asarray(xx, float), k), axes=(0, 0))
            Ls.append(c12_vec(Pi(coef), Q, R, th1, th2, Bg))
        Ls = np.array(Ls) if two else np.zeros((J + 1, K))           # (J+1, K)
        M = c2_mat(Q, R, th2, Bg, *res)[0] if two else np.eye(K)
        # unknowns y = (p, g); c = const + 2 b^T y + y^T S y
        n1 = J; Smat = np.zeros((n1 + K, n1 + K)); bvec = np.zeros(n1 + K)
        Smat[:n1, :n1] = A[1:, 1:]; bvec[:n1] = A[1:, 0]
        Smat[:n1, n1:] = Ls[1:]; Smat[n1:, :n1] = Ls[1:].T; bvec[n1:] = Ls[0]
        Smat[n1:, n1:] = M
        y = -np.linalg.solve(Smat, bvec); p, g = y[:n1], y[n1:]
        Pf = PEval(Bp, p)
        cpg = A[0, 0] + 0.5 * (1 + (end @ v) ** 2) + 2 * bvec @ y + y @ Smat @ y
        # Q step for fixed (p, g)
        Cp = np.sum(w * Pf(x, 1) ** 2); Bq = np.sum(w * Pf(x) ** 2)
        H1 = (Cp / th1) * G0 + th1 * Bq * G1 + 0.5 * np.outer(end, end)
        H12 = H12_Q(Pf, R, th1, th2, Bg, g, MQ) if two else 0 * G0
        H2 = H2_Q(R, th2, Bg, g, MQ, *res) if two else 0 * G0
        Ht = H1 + 2 * H12 + H2
        chk12 = v @ H12 @ v - g @ (Ls[0] + p @ Ls[1:]); chk2 = v @ H2 @ v - g @ M @ g
        hnew = np.linalg.solve(Ht[1:, 1:], -Ht[1:, 0]); v = np.concatenate([[1.0], hnew])
        cq = 0.5 + v @ Ht @ v
        c1v = c1_of(v, Pf)
        k = 1 - np.log(cq) / R
        hist.append(k)
        print(f"  it {it}: after (p,g) c={cpg:.12e}  after Q c={cq:.12e}  kappa2={k:.13f}  c1-part={c1v:.10e} (kappa from c1 alone {1-np.log(c1v)/R:.13f})"
              f"  |g|={np.abs(g).max():.2e}  chk c12 {chk12:.1e} c2 {chk2:.1e}  [{time.time()-t0:.0f}s]", flush=True)
    return v, p, g, hist

if __name__ == "__main__":
    D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Degree2Mol', 'params.json')))
    which = sys.argv[1]
    if which == 'at':          # joint.py at xi1 xi2 R MQ iters : start from the one-piece optimum at this R
        th1, th2, R = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]); MQ = int(sys.argv[5]); iters = int(sys.argv[6])
        c, h, lam = best(R, th1, MQ)
        run(f'at R={R}', R, th1, th2, lam, h, MQ, iters=iters); sys.exit()
    if which == 'ms':          # multistart: joint.py ms GRC|uncond krzQ|rand<seed>|Rq<R> iters  (only the start Q matters: the first step is the exact (p,g) solve)
        setting, kind, iters = sys.argv[2], sys.argv[3], int(sys.argv[4]); MQ = 30
        prm = D['RSC'] if setting == 'GRC' else D['Bernard_uncond']; th1 = float(Fraction(prm['xi'])); th2 = th1 if setting == 'GRC' else 25 / 298 - 1e-7
        R = float(Fraction(prm['R'])); lam = float(Fraction(prm['lam'])); hs = np.array([float(x) for x in prm['h']])
        if kind == 'krzQ':
            qz = np.array([.498939, 1.53685, 0, -2.7925, 0, 2.77524, 0, -1.01853] if setting == 'GRC' else [.499386, 1.58992, 0, -2.99061, 0, 3.01825, 0, -1.11694])
            cz = C.poly2cheb(qz / qz.sum()); h0 = np.zeros(MQ); h0[:4] = cz[1:8:2]
        elif kind.startswith('rand'):
            rng = np.random.default_rng(int(kind[4:])); h0 = hs + 0.3 * rng.standard_normal(MQ) / np.arange(1, MQ + 1) ** 2
        elif kind.startswith('Rq'):
            _, h0, _ = best(float(kind[2:]), th1, MQ)
        run(f'multistart {setting} {kind}', R, th1, th2, lam, h0, MQ, iters=iters, res=(100, 20, 20)); sys.exit()
    MQ = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    iters = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    two = not (len(sys.argv) > 4 and sys.argv[4] == 'onepiece')
    if which == 'GRC':
        prm = D['RSC']; th1 = float(Fraction(prm['xi'])); th2 = th1
    elif which == 'uncond':
        prm = D['Bernard_uncond']; th1 = float(Fraction(prm['xi'])); th2 = 25 / 298 - 1e-7
    elif which == 'tier2':
        prm = D['AT21_uncond']; th1 = float(Fraction(prm['xi'])); th2 = 25 / 298 - 1e-7
    R = float(Fraction(prm['R'])); lam = float(Fraction(prm['lam'])); h = np.array([float(s) for s in prm['h']])
    if MQ != len(h):   # re-optimise the one-piece Q at this degree first
        c, h, lam = best(R, th1, MQ)
    run(which + ('' if two else ' [control: one piece, same code]'), R, th1, th2, lam, h, MQ, iters=iters, two=two)

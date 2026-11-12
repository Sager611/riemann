# explore42 s15: joint lower bound for the continuous NB distance D(1/N)^2:
#   D^2 >= || P_{R + span{y_rho}} chi_1 ||^2,  R = P_(0,1) V^{-1} L^2(0,lam) (so W_lam^perp = L^2(1,inf) + R), y_rho = P_(0,1) Y_rho (Burnol/explore38 zero vectors).
# Pieces: Burnol Gram H, c (explore38 s3 formulas); P_(1,inf) Y_rho = J P_(0,1) C h_rho, h_rho = (1-M')^3 [1_(0,lam) V psi_rho] (exact log-polynomial algebra,
# cosine transform by even moments); R via prolate h-basis on H3 (s14). If completeness holds (W ⊖ C = span of W's kernels at zeros) this equals D^2.
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 s15_joint.py K N1,N2,...
import sys, math, time, resource
import numpy as np
import mpmath as mp
from scipy.special import sici
mp.mp.dps = 25
t0 = time.time()
K = int(sys.argv[1]); NS = [int(x) for x in sys.argv[2].split(',')]
OUT = '/Users/adriansager/playground/riemann/.scratch/explore42/'
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
J = 24
zc = np.load(OUT + 'zeros_cache.npz')
g = zc['g']; th = zc['th'][:K]; th1all = zc['th1'][:K]
gt = g[K:]; th1t = 0.5 * np.log(gt / (2 * np.pi)) + 1.0 / (48 * gt ** 2); Ttop = g[-1]
gam = np.concatenate([g[:K], -g[:K]]); rho = 0.5 + 1j * gam
ph = np.concatenate([th, -th])
V = (rho / (1 - rho)) ** 3 * np.exp(2j * ph)
dphi = np.concatenate([2 * th1all, 2 * th1all]) + 3.0 / np.abs(rho) ** 2
jj = np.arange(1, J + 1)
fact = np.array([math.factorial(2 * j + 1) for j in jj], dtype=float)
A = (4 * jj[None, :] * ((-1.0) ** jj)[None, :] * ((2 * jj / (2 * jj + 1.0)) ** 2)[None, :] / (fact[None, :] * (rho[:, None] + 2 * jj[None, :])))
bj = 2 * ((-1.0) ** jj) * (2 * jj / (2 * jj + 1.0)) / fact
tp = (2 * math.pi) ** (2 * jj)
n = len(gam)
def tail_diag(L):
    s = 2 * np.sum((1.0 / (0.25 + gt ** 2)) / (L + 2 * th1t + 3.0 / (0.25 + gt ** 2)))
    u = np.linspace(0, 60, 60001); t = Ttop * np.exp(u); lg = np.log(t / (2 * np.pi))
    s += float(np.trapezoid(2 * (lg / (2 * np.pi)) / (t * t * (L + lg)) * t, u))
    return s
def burnol(N):   # explore38 s3_scan.burnol (H_ij = (Y_j, Y_i), c_i = (chi_1, Y_i))
    lam = 1.0 / N; L = math.log(N)
    S = jj[:, None] + jj[None, :]
    Wjk = (2 * math.pi) ** (2 * S) * lam ** (2 * S + 1) / (2 * S + 1)
    P1 = lam ** (2 * jj[None, :] + 1 - rho[:, None]) / (2 * jj[None, :] + 1 - rho[:, None])
    Q1 = (np.conj(A) * tp[None, :]); A2 = A * tp[None, :]
    P2 = lam ** (2 * jj[:, None] + 1 - np.conj(rho)[None, :]) / (2 * jj[:, None] + 1 - np.conj(rho)[None, :])
    AW = A @ Wjk; Vc = np.conj(V)
    H = np.empty((n, n), dtype=np.complex128)
    for r0 in range(0, n, 128):
        r1 = min(n, r0 + 128)
        dg = gam[None, :] - gam[r0:r1, None]
        with np.errstate(divide='ignore', invalid='ignore'):
            blk = (1 - V[r0:r1, None] * Vc[None, :] * np.exp(-1j * L * dg)) / (1j * dg)
        idx = np.arange(r0, r1)
        blk[idx - r0, idx] = L + dphi[r0:r1]
        blk -= V[r0:r1, None] * (P1[r0:r1] @ Q1.T) + Vc[None, :] * (A2[r0:r1] @ P2) + AW[r0:r1] @ np.conj(A).T
        H[r0:r1, :] = np.conj(blk)
    c = (rho - 1) / rho ** 2
    corr = (bj * tp)[None, :] * (Vc[:, None] * lam ** (2 * jj[None, :] + 1 - np.conj(rho)[:, None]) / (2 * jj[None, :] + 1 - np.conj(rho)[:, None]))
    c = c - corr.sum(axis=1) - np.conj(A) @ (Wjk @ bj)
    return H, c
NMOM = 45
def h_moments(N):
    # even moments M_n = int_0^lam t^{2n} h_rho(t) dt of h_rho = (1-M')^3 g_rho, g_rho = V t^{-rho} + sum_j A_j (2pi)^{2j} t^{2j} on (0,lam)
    lam = 1.0 / N
    nn = np.arange(NMOM)[:, None]
    kf = np.array([1.0, 1.0, 2.0, 6.0])
    Mom = np.zeros((NMOM, n), dtype=np.complex128)
    def add(alpha, p):
        # alpha: (n,) coeffs of t^p (p array (n,) or scalar); apply (1-M')^3 on span{t^p, l^0..l^3}
        al = alpha.astype(np.complex128); be = np.zeros((4, n), dtype=np.complex128)
        lp = np.exp(p * math.log(lam)) if np.ndim(p) else lam ** p
        for _ in range(3):
            nb = be.copy()
            nb[0] = be[0] - al * lp / p
            for k in range(1, 4):
                nb[k] = be[k] - be[k - 1] / k
            be = nb; al = al * (1 + 1.0 / p)
        e = 2 * nn + 1
        Mom[:] += al[None, :] * np.exp((e + p) * math.log(lam)) / (e + p)
        for k in range(4):
            Mom[:] += be[k][None, :] * (lam ** e * kf[k] / e ** (k + 1))
    add(V, -rho)
    for j in range(1, J + 1):
        add(A[:, j - 1] * (2 * math.pi) ** (2 * j), 2.0 * j)
    return Mom
J1 = float(mp.quadosc(lambda u: mp.cos(u) * mp.log(u) / u, [1, mp.inf], omega=1))
def qf(v):
    a = 2 * math.pi * v
    r = float(mp.quad(lambda u: (mp.cos(u) - 1) * mp.log(u) / u, [a, 1]))
    return -2 * ((1 + math.log(a)) * sici(a)[1] + J1 - math.log(a) ** 2 / 2 + r)
def Tk(v, w):
    d = v[:, None] - w[None, :]; s = v[:, None] + w[None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        t1 = np.where(np.abs(d) < 1e-14, 2.0, np.sin(2 * np.pi * d) / (np.pi * d))
    return t1 + np.sin(2 * np.pi * s) / (np.pi * s)
for N in NS:
    lam = 1.0 / N; L = math.log(N)
    H, c = burnol(N)
    Mom = h_moments(N)
    # u-quadrature on (0,1)
    xu, wu = np.polynomial.legendre.leggauss(60); u = (xu + 1) / 2; wu = wu / 2
    nn = np.arange(NMOM)
    coef = np.array([(-1.0) ** k / math.factorial(2 * k) for k in nn])
    Pw = (2 * np.pi * u[:, None]) ** (2 * nn[None, :]) * coef[None, :]       # (Mu, NMOM)
    Ch = 2 * (Pw @ Mom)                                                       # (Mu, n): C h_rho on u-nodes
    B = float(np.real(np.vdot(c, np.linalg.solve(H, c))))
    # H <- Myy = H - H1, H1[i,j] = int Ch_j conj(Ch_i) du  (in place, row blocks)
    Chc = np.conj(Ch)
    for r0 in range(0, n, 256):
        r1 = min(n, r0 + 256)
        H[r0:r1, :] -= (Chc[:, r0:r1] * wu[:, None]).T @ Ch
    h1max = float(np.abs(np.einsum('ki,ki,k->i', Chc, Ch, wu)).max())
    BW = float(np.real(np.vdot(c, np.linalg.solve(H, c))))
    # R part (prolate h-basis on H3)
    Mv = 60
    xv, wxv = np.polynomial.legendre.leggauss(Mv); xv = (xv + 1) / 2; wxv = wxv / 2
    pp = 3; v = lam * xv ** pp; wv = wxv * lam * pp * xv ** (pp - 1); sw = np.sqrt(wv)
    qv = np.array([qf(t) for t in v])
    Kz = sw[:, None] * Tk(v, v) * sw[None, :]
    Cc = np.stack([sw * np.log(lam / v) ** k for k in range(3)], axis=1)
    Qm, _ = np.linalg.qr(Cc, mode='complete'); Z = Qm[:, 3:]
    MRR = Z.T @ (np.eye(Mv) - Kz) @ Z
    tauR = -(Z.T @ (sw * qv))
    rmat = 2 * np.cos(2 * np.pi * v[:, None] * u[None, :]) @ (wu[:, None] * Ch)
    MRy = Z.T @ (sw[:, None] * rmat)
    Fv = float(tauR @ np.linalg.solve(MRR, tauR))
    # Schur complement on the small R block: joint = F' + ...; solve G via block elimination without forming G
    X = np.linalg.solve(MRR, MRy)                    # (Mr, n)
    H -= np.conj(MRy).T @ X                          # Schur complement S = Myy - MyR MRR^{-1} MRy
    cc = c - np.conj(MRy).T @ np.linalg.solve(MRR, tauR)
    val = Fv + float(np.real(np.vdot(cc, np.linalg.solve(H, cc))))
    ev = np.array([0.0])
    td = tail_diag(L)
    pass
    print(f'N={N:5d} K={K} F={Fv:.10f} B_K={B:.10f} BW_K(y only)={BW:.10f} joint_K={val:.10f} tail={td:.4e} | joint+tail={val+td:.10f} (x L {L*(val+td):.8f}) | B+eps-style F+B={Fv+B+td:.10f} | maxH1diag={h1max:.2e} {time.time()-t0:.1f}s', flush=True)
print(f'maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')

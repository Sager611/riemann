# explore38 s2: Burnol's zero vectors Y^lambda_rho (math/0103058, Thm-Def + Thm theo) at FINITE lambda = 1/N:
# exact Gram (Y_rho, Y_rho') and scalar products (chi_1, Y_rho); B_K(N) = c^H H^{-1} c <= D(1/N)^2 <= d_N^2 (unconditional lower bound).
# (Y_w,Y_w') = (psi_w,psi_w') - int_0^lambda Vpsi_w conj(Vpsi_w'), Vpsi_w(t) = V(w) t^{-w} + F_w(t), F_w = sum_{j>=1} a_j(w) (2 pi t)^{2j},
# a_j(w) = 4j(-1)^j (2j/(2j+1))^2/((2j+1)!(w+2j)),
# V(w) = (w/(1-w))^3 / chi(w); on the line V(1/2+it) = ((1/2+it)/(1/2-it))^3 e^{2 i theta(t)}.
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 s2_burnol.py K N1,N2,...   (K zero pairs)
import sys, math, time, resource
import numpy as np
import mpmath as mp
mp.mp.dps = 30
t0 = time.time()
K = int(sys.argv[1]); NS = [int(x) for x in sys.argv[2].split(',')]
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
J = 24
def zeros(K):
    zf = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')
    g = zf[:K].copy()
    nex = min(K, 300)
    for n in range(1, nex + 1):  # file zeros are 3e-6 off at the bottom; exact low zeros from mpmath
        g[n - 1] = float(mp.zetazero(n).imag)
    return g
def prep(K):
    g = zeros(K)
    th = np.array([float(mp.siegeltheta(t)) for t in g])
    th1 = np.array([float(mp.siegeltheta(t, derivative=1)) for t in g])
    gam = np.concatenate([g, -g])
    rho = 0.5 + 1j * gam
    ph = np.concatenate([th, -th])
    V = (rho / (1 - rho)) ** 3 * np.exp(2j * ph)
    dphi = np.concatenate([2 * th1, 2 * th1]) + 3.0 / np.abs(rho) ** 2   # Phi'(gamma) - L = Re V'/V(rho)
    jj = np.arange(1, J + 1)
    fact = np.array([math.factorial(2 * j + 1) for j in jj], dtype=float)
    # a_j(w): phi_w = U psi_w = k(t) + w int_1^inf u^{w-1} k(tu) du (the boundary term k(t) = sin(2 pi t)/(pi t) included),
    # phi_w - U(w) t^{-w} = sum_{j>=1} 2(-1)^j (2pi t)^{2j} (2j/(w+2j))/(2j+1)!  (Burnol's series, k = 0); F_w = (1-M)^2 of it
    A = (4 * jj[None, :] * ((-1.0) ** jj)[None, :] * ((2 * jj / (2 * jj + 1.0)) ** 2)[None, :]
         / (fact[None, :] * (rho[:, None] + 2 * jj[None, :])))          # a_j(rho)
    bj = 2 * ((-1.0) ** jj) * (2 * jj / (2 * jj + 1.0)) / fact
    return gam, rho, V, dphi, A, jj, bj
def burnol(N, gam, rho, V, dphi, A, jj, bj):
    lam = 1.0 / N; L = math.log(N)
    n = len(gam)
    tp = (2 * math.pi) ** (2 * jj)
    # main term
    dg = gam[None, :] - gam[:, None]          # gamma_b - gamma_a
    with np.errstate(divide='ignore', invalid='ignore'):
        M = (1 - V[:, None] * np.conj(V)[None, :] * np.exp(-1j * L * dg)) / (1j * dg)
    M[np.arange(n), np.arange(n)] = L + dphi
    # F terms
    P1 = lam ** (2 * jj[None, :] + 1 - rho[:, None]) / (2 * jj[None, :] + 1 - rho[:, None])    # [a, j]
    I1 = P1 @ (np.conj(A) * tp[None, :]).T                                                   # [a, b]
    I2 = (A * tp[None, :]) @ (lam ** (2 * jj[:, None] + 1 - np.conj(rho)[None, :]) / (2 * jj[:, None] + 1 - np.conj(rho)[None, :]))
    S = jj[:, None] + jj[None, :]
    Wjk = (2 * math.pi) ** (2 * S) * lam ** (2 * S + 1) / (2 * S + 1)
    I3 = A @ Wjk @ np.conj(A).T
    M -= V[:, None] * I1 + np.conj(V)[None, :] * I2 + I3
    del I1, I2, I3, P1
    c = (rho - 1) / rho ** 2
    corr = (bj * tp)[None, :] * (np.conj(V)[:, None] * lam ** (2 * jj[None, :] + 1 - np.conj(rho)[:, None]) / (2 * jj[None, :] + 1 - np.conj(rho)[:, None]))
    c2 = np.conj(A) @ (Wjk @ (bj * tp / tp))  # sum_k conj(a_k) (2pi)^{2k} sum_j b_j (2pi)^{2j} lam^{2j+2k+1}/(2j+2k+1)
    c = c - corr.sum(axis=1) - c2
    H = np.conj(M)          # H_ij = (Y_j, Y_i) = M_ji = conj(M_ij)
    del M
    x = np.linalg.solve(H, c)
    B = float(np.real(np.vdot(c, x)))
    Bd = float(np.sum(np.abs(c) ** 2 / np.real(np.diag(H))))
    ev = np.linalg.eigvalsh(H) if n <= 800 else None
    return B, Bd, (ev.min() if ev is not None else float('nan')), float(np.sum(np.abs(c) ** 2))
gam, rho, V, dphi, A, jj, bj = prep(K)
print(f'K={K} pairs, gamma_K={gam[K-1]:.3f}, C_K=sum|c|^2(limit)={float(np.sum(1/np.abs(rho)**2)):.10f}, C={C:.10f}, prep {time.time()-t0:.1f}s', flush=True)
res = []
for N in NS:
    B, Bd, emin, cc = burnol(N, gam, rho, V, dphi, A, jj, bj)
    L = math.log(N)
    res.append((N, B, Bd))
    print(f'N={N:6d} L={L:.4f} B={B:.12e} BlogN={B*L:.8f} Bdiag*L={Bd*L:.8f} L^2B-C_K L={L*L*B - cc*L:+.5f} eigmin(H)={emin:.4f}', flush=True)
np.save(OUT + f'burnol_K{K}.npy', np.array(res))
print(f'done {time.time()-t0:.1f}s maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')

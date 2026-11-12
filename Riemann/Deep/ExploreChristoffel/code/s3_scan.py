# explore38 s3: memory-lean Burnol/Christoffel lower bound B_K(N) for many N; diagonal tail over zeros > K.
# Same formulas as s2_burnol.py (verified: H positive definite for N >= 2). Row-block construction keeps temporaries small.
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 s3_scan.py K Nspec OUTTAG
#   Nspec: comma list, or 'all:NMAX' (every N in [2, NMAX]) or 'geo:a:b:m' (m log-spaced N in [a,b])
import sys, math, time, resource
import numpy as np
import mpmath as mp
mp.mp.dps = 30
t0 = time.time()
K = int(sys.argv[1]); spec = sys.argv[2]; TAG = sys.argv[3]
if spec.startswith('all:'):
    NS = list(range(2, int(spec[4:]) + 1))
elif spec.startswith('geo:'):
    _, a, b, m = spec.split(':'); NS = sorted(set(int(round(x)) for x in np.geomspace(float(a), float(b), int(m))))
else:
    NS = [int(x) for x in spec.split(',')]
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
J = 24
zc = np.load(OUT + 'zeros_cache.npz')   # from s0_cache.py
g = zc['g']; th = zc['th'][:K]; th1all = zc['th1'][:K]
# tail (zeros K+1 .. 22491 from file, then density integral beyond)
gt = g[K:]
th1t = 0.5 * np.log(gt / (2 * np.pi)) + 1.0 / (48 * gt ** 2)
Ttop = g[-1]
gam = np.concatenate([g[:K], -g[:K]])
rho = 0.5 + 1j * gam
ph = np.concatenate([th, -th])
V = (rho / (1 - rho)) ** 3 * np.exp(2j * ph)
dphi = np.concatenate([2 * th1all, 2 * th1all]) + 3.0 / np.abs(rho) ** 2
jj = np.arange(1, J + 1)
fact = np.array([math.factorial(2 * j + 1) for j in jj], dtype=float)
A = (4 * jj[None, :] * ((-1.0) ** jj)[None, :] * ((2 * jj / (2 * jj + 1.0)) ** 2)[None, :]
     / (fact[None, :] * (rho[:, None] + 2 * jj[None, :])))
bj = 2 * ((-1.0) ** jj) * (2 * jj / (2 * jj + 1.0)) / fact
tp = (2 * math.pi) ** (2 * jj)
n = len(gam)
CK = float(np.sum(1 / np.abs(rho) ** 2))
def tail_diag(L):
    s = 2 * np.sum((1.0 / (0.25 + gt ** 2)) / (L + 2 * th1t + 3.0 / (0.25 + gt ** 2)))
    # beyond Ttop: density (1/2pi) log(t/2pi), weight 1/t^2, Christoffel 1/(L + log(t/2pi))
    # substitute t = Ttop e^u: integrand dt = Ttop e^u du
    u = np.linspace(0, 60, 60001); t = Ttop * np.exp(u); lg = np.log(t / (2 * np.pi))
    f = 2 * (lg / (2 * np.pi)) / (t * t * (L + lg)) * t
    s += float(np.trapezoid(f, u))
    return s
def tail_C():
    s = 2 * np.sum(1.0 / (0.25 + gt ** 2))
    u = np.linspace(0, 60, 60001); t = Ttop * np.exp(u); lg = np.log(t / (2 * np.pi))
    s += float(np.trapezoid(2 * (lg / (2 * np.pi)) / t, u))
    return s
TC = tail_C()
H = np.empty((n, n), dtype=np.complex128)
def burnol(N):
    lam = 1.0 / N; L = math.log(N)
    S = jj[:, None] + jj[None, :]
    Wjk = (2 * math.pi) ** (2 * S) * lam ** (2 * S + 1) / (2 * S + 1)
    P1 = lam ** (2 * jj[None, :] + 1 - rho[:, None]) / (2 * jj[None, :] + 1 - rho[:, None])
    Q1 = (np.conj(A) * tp[None, :])                     # [b, j]
    A2 = A * tp[None, :]                                # [a, j]
    P2 = lam ** (2 * jj[:, None] + 1 - np.conj(rho)[None, :]) / (2 * jj[:, None] + 1 - np.conj(rho)[None, :])   # [j, b]
    AW = A @ Wjk
    Vc = np.conj(V)
    for r0 in range(0, n, 128):
        r1 = min(n, r0 + 128)
        dg = gam[None, :] - gam[r0:r1, None]
        with np.errstate(divide='ignore', invalid='ignore'):
            blk = (1 - V[r0:r1, None] * Vc[None, :] * np.exp(-1j * L * dg)) / (1j * dg)
        idx = np.arange(r0, r1)
        blk[idx - r0, idx] = L + dphi[r0:r1]
        blk -= V[r0:r1, None] * (P1[r0:r1] @ Q1.T) + Vc[None, :] * (A2[r0:r1] @ P2) + AW[r0:r1] @ np.conj(A).T
        H[r0:r1, :] = np.conj(blk)                      # H_ij = (Y_j, Y_i) = conj(M_ij)
    c = (rho - 1) / rho ** 2
    corr = (bj * tp)[None, :] * (Vc[:, None] * lam ** (2 * jj[None, :] + 1 - np.conj(rho)[:, None]) / (2 * jj[None, :] + 1 - np.conj(rho)[:, None]))
    c = c - corr.sum(axis=1) - np.conj(A) @ (Wjk @ bj)
    Bd = float(np.sum(np.abs(c) ** 2 / np.real(np.diag(H))))
    x = np.linalg.solve(H, c)
    B = float(np.real(np.vdot(c, x)))
    return B, Bd, L
rows = []
for i, N in enumerate(NS):
    B, Bd, L = burnol(N)
    td = tail_diag(L)
    rows.append((float(N), L, B, Bd, td))
    if i < 3 or N in (10, 40, 100, 320, 640, 1280) or i == len(NS) - 1 or (len(NS) < 40):
        print(f'N={N:7d} L={L:.4f} B_K={B:.10e} Bdiag_K={Bd:.10e} tail={td:.4e} (B_K+tail)L={(B+td)*L:.8f} L^2(B_K+tail)-CL={L*L*(B+td)-C*L:+.5f}', flush=True)
rows = np.array(rows)
np.save(OUT + f'scan_{TAG}.npy', rows)
print(f'K={K} gamma_K={g[K-1]:.2f} C_K={CK:.10f} tailC={TC:.6e} C_K+tailC={CK+TC:.10f} C={C:.10f}')
print(f'done {len(NS)} N in {time.time()-t0:.1f}s maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')

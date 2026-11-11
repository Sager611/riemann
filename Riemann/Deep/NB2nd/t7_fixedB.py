# nb2nd t7: fixed (N-independent) head corrections. For B*_K = argmax of the N-free head problem (support 1..K), check
# L^2 I(V_N + B*_K/L) - C L - [J - Sigma_K] - R(N) -> 0 (claim: = -2<E, zeta B>-type coupling = O(N^-eps)), and the
# second-order values of these constructions vs the optimum.
import math, numpy as np, mpmath as mp
O = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
Z = np.load(O + 'G640.npz'); G = Z['G']; b = Z['b']; mu = Z['mu'].astype(float); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
mp.mp.dps = 20
a1 = float(mp.zeta(0, derivative=2) / mp.zeta(0) - (mp.zeta(0, derivative=1) / mp.zeta(0)) ** 2)
z = [0.0] + [float(mp.zeta(2 * k + 1)) for k in range(1, 200)]
def ell(n):
    s0 = -(math.log(2 * math.pi) * (1 + math.log(n) - EG) + a1) / n
    if n == 1: tail = sum((z[k] - 1) / (2 * k * (2 * k + 1)) for k in range(1, 199)) + (1 - math.log(2))
    else:
        x = 1.0 / n; tail = sum(z[k] * x ** (2 * k + 1) / (2 * k * (2 * k + 1)) for k in range(1, 60))
    return s0 - tail
ks = np.arange(1, NMAX + 1); lth = np.array([ell(int(n)) for n in ks])
A2 = np.load(O + 't2_allN.npz'); NsA = A2['N']; V1 = A2['V1']; RA = A2['R']
for K in [5, 20, 60]:
    bs = np.linalg.solve(G[:K, :K], lth[:K]); SigK = float(lth[:K] @ bs)
    out = []
    for N in range(max(100, 2 * K), NMAX + 1):
        L = math.log(N); c = -mu[:N] * (1 - np.log(ks[:N]) / L); c[:K] -= bs / L   # A = V_N + B/L  <=>  c -> c - beta/L
        I = 1 - 2 * c @ b[:N] + c @ G[:N, :N] @ c
        i = N - 2; out.append((N, L * L * I - C * L, L * L * I - C * L - (J - SigK) - RA[i], V1[i]))
    o = np.array(out)
    for lo, hi in [(100, 250), (250, 450), (450, 640)]:
        m = (o[:, 0] >= lo) & (o[:, 0] <= hi)
        if m.any():
            print(f'K={K} (J-Sigma_K={J-SigK:.4f}) N in [{lo},{hi}]: construction L^2I-CL mean {o[m,1].mean():+.4f} sd {o[m,1].std():.4f}; '
                  f'residual vs J(B)+R(N): mean {o[m,2].mean():+.4f} sd {o[m,2].std():.4f}; optimum mean {o[m,3].mean():+.4f}')

# nb2nd t3: tail jump/tilt family. A_j = sum_{n<=N} mu(n)(1 - (1 - j/L) log n/L) n^-s  (profile P(x) = 1 - (1-j/L)x, P(1) = j/L).
# Prediction (heuristic, derivation.md): L^2 I(A_j) - C L = J + D(N) - j C + j^2 X(N) + O(1/L), X(N) = ||1 - zeta M_N||^2 (sharp Moebius cut),
# and L <R_N, R_N^sharp> -> C/2. Plus head correction B_K/L (beta* from the N-free head problem): construction value vs optimum.
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
Lc = np.linalg.cholesky(G); y = np.linalg.solve(Lc, b); d2 = 1 - np.cumsum(y * y)
def ip(c1, c2, N):  # <F_c1, F_c2>, F_c = 1 + zeta sum c_k k^-s  (NB convention c)
    return 1 - b[:N] @ c1 - b[:N] @ c2 + c1 @ G[:N, :N] @ c2
rows = []
Kh = 60; bh = np.linalg.solve(G[:Kh, :Kh], lth[:Kh]); SigK = float(lth[:Kh] @ bh)
for N in range(100, NMAX + 1):
    L = math.log(N); lk = np.log(ks[:N])
    cr = -mu[:N] * (1 - lk / L); csh = -mu[:N]
    I = ip(cr, cr, N); X = ip(csh, csh, N); RS = ip(cr, csh, N)
    # A_j: c_j = -mu (1 - (1 - j/L) lk/L) = cr - mu j lk/L^2 ; I(c_j) = I + 2 j I1 + j^2 I2
    dlt = -mu[:N] * lk / L ** 2
    I1 = (G[:N, :N] @ cr - b[:N]) @ dlt; I2 = dlt @ G[:N, :N] @ dlt
    jopt = -I1 / I2; Ij = I - I1 * I1 / I2
    # add head correction: A = A_j + B/L, B = sum beta*_n n^-s  -> c += -beta*/L on n <= Kh; optimise jointly (j, scale t of beta*) exactly
    hb = np.zeros(N); hb[:Kh] = -bh / L
    M = np.array([[I2, dlt @ G[:N, :N] @ hb], [dlt @ G[:N, :N] @ hb, hb @ G[:N, :N] @ hb]])
    g = np.array([I1, (G[:N, :N] @ cr - b[:N]) @ hb])
    Ijh = I - g @ np.linalg.solve(M, g)
    rows.append((N, L, L * L * I - C * L - J, X, L * RS, L * L * 2 * I1, L * L * I2, jopt, L * L * Ij - C * L, L * L * Ijh - C * L, L * L * d2[N - 1] - C * L))
R = np.array(rows); np.save(O + 't3_rows.npy', R)
cols = ['N', 'L', 'R=L2I-CL-J', 'X(N)', 'L<R,Rsh>', 'L2*2I1', 'L2*I2', 'j_opt', 'tilt:L2I-CL', 'tilt+head:L2I-CL', 'opt:L2d2-CL']
for lo, hi in [(100, 200), (200, 400), (400, 640)]:
    m = (R[:, 0] >= lo) & (R[:, 0] <= hi)
    print(f'N in [{lo},{hi}]: ' + '; '.join(f'{c} {R[m, i].mean():+.5f}({R[m, i].std():.4f})' for i, c in enumerate(cols) if i >= 2))
print(f'C = {C:.6f}, C/2 = {C/2:.6f}; head K={Kh}: Sigma_K = {SigK:.6f}, J - Sigma_K = {J - SigK:.5f}')
for N in [160, 320, 640]:
    r = R[R[:, 0] == N][0]; print('N=%d: ' % N + '; '.join(f'{c} {r[i]:+.5f}' for i, c in enumerate(cols) if i >= 2))

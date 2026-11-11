# nb2nd t4: how much of the tail gain is the single jump/tilt direction? For each N in [100,640]:
# S_N (full), S_h (head support 1..K free), S_ht (head + tilt direction delta = -mu log n/L^2), S_hte (head + tilt + end support [N e^-0.5, N]).
import math, numpy as np
O = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
Z = np.load(O + 'G640.npz'); G = Z['G']; b = Z['b']; mu = Z['mu'].astype(float); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
ks = np.arange(1, NMAX + 1); K = 60
LcF = np.linalg.cholesky(G)  # Cholesky of G_N = leading block of Cholesky of G_640
from scipy.linalg import solve_triangular
def gain(idx, GN, r, extra=None):  # max over span{e_i : i in idx} (+ column 'extra') of 2 x.r - x.G x = g^T M^-1 g
    cols_g = [r[idx]]; M = GN[np.ix_(idx, idx)]
    if extra is not None:
        ge = extra @ r; Ge = GN[:, idx].T @ extra; ee = extra @ GN @ extra
        M = np.block([[M, Ge[:, None]], [Ge[None, :], np.array([[ee]])]]); g = np.concatenate([r[idx], [ge]])
    else:
        g = r[idx]
    return float(g @ np.linalg.solve(M, g))
rows = []
for N in range(100, NMAX + 1, 3):
    L = math.log(N); lk = np.log(ks[:N]); GN = G[:N, :N]
    cr = -mu[:N] * (1 - lk / L); r = GN @ cr - b[:N]
    y = solve_triangular(LcF[:N, :N], r, lower=True); S = L * L * float(y @ y)
    H = np.arange(K); dl = -mu[:N] * lk / L ** 2
    lo = int(math.ceil(N * math.exp(-0.5))); E = np.arange(lo - 1, N); HE = np.concatenate([H, E])
    Sh = L * L * gain(H, GN, r); Sht = L * L * gain(H, GN, r, dl)
    Shte = L * L * gain(HE, GN, r, dl); She = L * L * gain(HE, GN, r)
    rows.append((N, S, Sh, Sht, She, Shte))
R = np.array(rows); np.save(O + 't4_rows.npy', R)
for lo_, hi_ in [(100, 250), (250, 450), (450, 640)]:
    m = (R[:, 0] >= lo_) & (R[:, 0] <= hi_)
    S, Sh, Sht, She, Shte = (R[m, i] for i in range(1, 6))
    print(f'N in [{lo_},{hi_}]: tail gain S-Sh mean {np.mean(S-Sh):.4f} sd {np.std(S-Sh):.4f}; tilt only Sht-Sh {np.mean(Sht-Sh):.4f} ({np.std(Sht-Sh):.4f}); '
          f'end[Ne^-.5,N] She-Sh {np.mean(She-Sh):.4f}; tilt+end Shte-Sh {np.mean(Shte-Sh):.4f}; fraction tilt {np.mean((Sht-Sh)/(S-Sh)):.2f}, tilt+end {np.mean((Shte-Sh)/(S-Sh)):.2f}')

# explore61 s6: zero formula along the hybrid optimizer.  Q = dist(P_W1 h, V_N)^2 vs (1/L) sum_rho |M[P_W1perp T e](rho)|^2,
# R = |P_W1perp T e|^2, Z = sum_rho |M[g](rho)|^2/|g|^2 (zero visibility of the roughness g = P_perp Te).  usage: s6_optzero.py npzfile NZ MMAX
import sys, os, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memcap import check, mem
from s1_split import w1_perp_gram
t0 = time.time()
f, NZ, MMAX = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
z = np.load(f); head = [int(k) for k in z['head']]; e = z['e']
N = int(f.split('_N')[1].split('_')[0]); Qg = int(f.split('_Q')[1].split('_')[0]); L = math.log(N)
gam = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')[:NZ]
s = 0.5 + 1j * gam
def Dsum(s):
    out = np.zeros(len(s), complex)
    for m0 in range(1, MMAX + 1, 10000):
        m = np.arange(m0, min(m0 + 10000, MMAX + 1), dtype=np.float64); l, r = 1 / (m + 1), 1 / m
        av = m * (m + 1) * np.log1p(1 / m)
        for i0 in range(0, len(s), 16):
            ss = s[i0:i0 + 16, None]
            out[i0:i0 + 16] += np.sum((r ** (ss - 1) - l ** (ss - 1)) / (ss - 1) - av * (r ** ss - l ** ss) / ss, axis=1)
    return out + 1 / (1 - s)
Ds = Dsum(s)
Om, B, Cp = w1_perp_gram(head, Qg, MMAX)
Ge = np.zeros(len(s), complex)
for j, k in enumerate(head):
    n = np.arange(1, (MMAX * Qg) // k + 2, dtype=np.int64); kn = k * n; m = kn // Qg
    ok = (kn % Qg != 0) & (m <= MMAX); m = m[ok].astype(np.float64); x = Qg / kn[ok].astype(np.float64)
    l, r = 1 / (m + 1), 1 / m; q = (x - l) / (r - l); lx, ll, lr = np.log(x), np.log(l), np.log(r)
    for i0 in range(0, len(s), 16):
        ss = s[i0:i0 + 16, None]
        Ge[i0:i0 + 16] += -e[j] * np.sum((np.exp(ss * lx) - np.exp(ss * ll)) - q * (np.exp(ss * lr) - np.exp(ss * ll)), axis=1) / s[i0:i0 + 16]
    check(280, f'head {j}')
beta = float(e @ B) / Cp
Ge -= beta * Ds
g2 = 2 * np.abs(Ge) ** 2; cs = np.cumsum(g2)
A = np.mean(np.abs(Ge[NZ // 2:]) ** 2 * gam[NZ // 2:] ** 2); T = gam[-1]
tail = 2 * A * (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi * T)
Zs = cs[-1] + tail
R = float(e @ Om @ e) - beta * float(e @ B)
P, Rz = float(z['P']), float(z['R']); Qd = P - Rz
print(f'N={N} Q={Qg}: R(s1)={Rz:.4e} R(here, M={MMAX})={R:.4e} | measured L^3 Q={L**3*Qd:.4f}  predicted L^2 sum|G_e|^2={L*L*Zs:.4f} '
      f'(partial 10/100/{NZ}: {L*L*cs[9]:.4f} {L*L*cs[min(99,NZ-1)]:.4f} {L*L*cs[-1]:.4f}) ratio {L**3*Qd/(L*L*Zs):.3f} | Z = sum|G|^2/R = {Zs/R:.3f}, c_e = LQ/R = {L*Qd/Rz:.3f} | t={time.time()-t0:.0f}s mem {mem()}')

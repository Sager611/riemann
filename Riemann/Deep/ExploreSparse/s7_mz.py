# explore61 s7: is the zero visibility Z of the optimizer's roughness g = P_W1perp(Te) an MZ average?  Compare
#   S_z = 2 sum_{0<gamma<=T} |G(1/2+i gamma)|^2   with   S_I = 2 int_0^T |G(1/2+it)|^2 theta'(t)/pi dt  (smooth zero density),
# and |g|^2 = (1/pi) int_0^inf |G|^2 dt (Plancherel, G(conj s) = conj G(s)).  usage: s7_mz.py npz T dt MMAX
import sys, os, math, time
import numpy as np
from scipy.special import digamma
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memcap import check, mem
from s1_split import w1_perp_gram
t0 = time.time()
f, T, dt, MMAX = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
z = np.load(f); head = [int(k) for k in z['head']]; e = z['e']; Qg = int(f.split('_Q')[1].split('_')[0])
gam = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy'); gam = gam[gam <= T]
tt = np.arange(dt / 2, T, dt)
def Gfun(s):
    out = np.zeros(len(s), complex)
    for m0 in range(1, MMAX + 1, 5000):
        m = np.arange(m0, min(m0 + 5000, MMAX + 1), dtype=np.float64); l, r = 1 / (m + 1), 1 / m; av = m * (m + 1) * np.log1p(1 / m)
        for i0 in range(0, len(s), 16):
            ss = s[i0:i0 + 16, None]
            out[i0:i0 + 16] += -beta * np.sum((r ** (ss - 1) - l ** (ss - 1)) / (ss - 1) - av * (r ** ss - l ** ss) / ss, axis=1)
    out += -beta / (1 - s)
    for j, k in enumerate(head):
        n = np.arange(1, (MMAX * Qg) // k + 2, dtype=np.int64); kn = k * n; m = kn // Qg
        ok = (kn % Qg != 0) & (m <= MMAX); m = m[ok].astype(np.float64); x = Qg / kn[ok].astype(np.float64)
        l, r = 1 / (m + 1), 1 / m; q = (x - l) / (r - l); lx, ll, lr = np.log(x), np.log(l), np.log(r)
        for i0 in range(0, len(s), 16):
            ss = s[i0:i0 + 16, None]
            out[i0:i0 + 16] += -e[j] * np.sum((np.exp(ss * lx) - np.exp(ss * ll)) - q * (np.exp(ss * lr) - np.exp(ss * ll)), axis=1) / s[i0:i0 + 16]
    return out
Om, B, Cp = w1_perp_gram(head, Qg, MMAX)
beta = float(e @ B) / Cp
R = float(e @ Om @ e) - beta * float(e @ B)
Gz = Gfun(0.5 + 1j * gam); check(280, 'Gz')
Gt = Gfun(0.5 + 1j * tt); check(280, 'Gt')
thp = np.real(digamma(0.25 + 0.5j * tt)) / 2 - math.log(math.pi) / 2  # theta'(t)
Sz = 2 * np.sum(np.abs(Gz) ** 2); SI = 2 * np.sum(np.abs(Gt) ** 2 * np.maximum(thp, 0) / math.pi) * dt
g2 = np.sum(np.abs(Gt) ** 2) * dt / math.pi
for Tc in (50, 100, 200, 500, T):
    a = 2 * np.sum(np.abs(Gz[gam <= Tc]) ** 2); b = 2 * np.sum((np.abs(Gt) ** 2 * np.maximum(thp, 0) / math.pi)[tt <= Tc]) * dt
    c = np.sum((np.abs(Gt) ** 2)[tt <= Tc]) * dt / math.pi
    print(f'T<={Tc:6.0f}: S_z={a:.4e} S_I={b:.4e} S_z/S_I={a/b:.3f} | share of |g|^2 below T: {c/R:.3f}', flush=True)
print(f'R={R:.4e}  Plancherel |g|^2 (t<{T}) = {g2:.4e}  Z = S_z/R = {Sz/R:.3f}  MZ-average Z_I = S_I/R = {SI/R:.3f}  t={time.time()-t0:.0f}s mem {mem()}')

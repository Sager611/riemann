# explore61 s5: zero formula for the in-W1 rate.  Test  kappa(theta) := lim L dist(P_W1 rho_theta, V_N)^2  against
#   K0(theta) := sum_rho |M[P_W1perp rho_theta](rho)|^2   (M = Mellin transform  int_0^inf f(x) x^{s-1} dx; for chi: sum 1/|rho|^2 = C),
# using M[P_W1 rho_theta](rho) = M[rho_theta](rho) - M[P_perp rho_theta](rho) = -M[P_perp rho_theta](rho) (zeta(rho) = 0).
# P_perp rho_theta = (on I_m) -sum_jumps (1[x<x_j] - q_j) - beta (1/x - avg_m 1/x);  (on (1,inf)) -beta/x;  beta = B/C' (s1_split).
# Zeros: Riemann/Deep/ExploreHelio/results/zeta_zeros.npy (22491 ordinates < 20000).
# usage: uv run --no-project --with numpy --with scipy python3 s5_zero.py NZ MMAX
import sys, os, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memcap import check, mem
from s1_split import w1_perp_gram
t0 = time.time()
NZ = int(sys.argv[1]); MMAX = int(sys.argv[2])
gam = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')[:NZ]
s = 0.5 + 1j * gam
C = 2 + 0.57721566490153286061 - math.log(4 * math.pi)
print(f'check: 2 sum_{{n<={NZ}}} 1/|rho|^2 = {2*np.sum(1/np.abs(s)**2):.6f} (C = {C:.6f})', flush=True)
# D(s) = sum_m int_{I_m} (1/x - avg_m) x^{s-1} dx  + 1/(1-s)   (the beta-direction), m <= MMAX
def Dsum(s):
    out = np.zeros(len(s), complex)
    for m0 in range(1, MMAX + 1, 20000):
        m = np.arange(m0, min(m0 + 20000, MMAX + 1), dtype=np.float64)
        l, r = 1 / (m + 1), 1 / m
        av = m * (m + 1) * np.log1p(1 / m)
        for i0 in range(0, len(s), 16):
            ss = s[i0:i0 + 16, None]
            out[i0:i0 + 16] += np.sum((r ** (ss - 1) - l ** (ss - 1)) / (ss - 1) - av * (r ** ss - l ** ss) / ss, axis=1)
    return out + 1 / (1 - s)
Ds = Dsum(s)
check(280, 'Ds')
print(f'D(s) done t={time.time()-t0:.0f}s', flush=True)
TH = [(2, 3), (3, 4), (3, 5), (4, 5), (5, 6), (4, 7), (5, 7), (5, 8), (2, 5), (3, 7), (3, 8), (2, 7), (2, 9), (3, 10)]
for a, b in TH:
    Q, k = a, b
    n = np.arange(1, (MMAX * Q) // k + 2, dtype=np.int64)
    kn = k * n; m = kn // Q
    ok = (kn % Q != 0) & (m <= MMAX)
    m = m[ok].astype(np.float64); x = Q / kn[ok].astype(np.float64)
    l, r = 1 / (m + 1), 1 / m
    q = (x - l) / (r - l)
    lx, ll, lr = np.log(x), np.log(l), np.log(r)
    Js = np.zeros(len(s), complex)
    for i0 in range(0, len(s), 32):
        ss = s[i0:i0 + 32, None]
        Js[i0:i0 + 32] = -np.sum((np.exp(ss * lx) - np.exp(ss * ll)) - q * (np.exp(ss * lr) - np.exp(ss * ll)), axis=1) / s[i0:i0 + 32]
    Om, B, Cp = w1_perp_gram([k], Q, MMAX)
    beta = float(B[0]) / Cp
    om2 = float(Om[0, 0]) - beta * float(B[0])
    G = Js - beta * Ds
    g2 = 2 * np.abs(G) ** 2  # rho and conj(rho)
    cs = np.cumsum(g2)
    # tail beyond gamma_NZ: |G|^2 ~ A/gamma^2 on average, zero density log(t/2pi)/(2pi):  2 int_T^inf A/t^2 log(t/2pi)/(2pi) dt
    A = np.mean(np.abs(G[NZ // 2:]) ** 2 * gam[NZ // 2:] ** 2)
    T = gam[-1]; tail = 2 * A * (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi * T)
    print(f'{a}/{b}: om^2={om2:.6f} beta={beta:+.5f}  K0 partial: n<=10 {cs[9]:.5f}, 100 {cs[min(99,NZ-1)]:.5f}, {NZ} {cs[-1]:.5f}, +tail {tail:.5f} -> {cs[-1]+tail:.5f}'
          f'  |G(rho1)| {abs(G[0]):.4f}  t={time.time()-t0:.0f}s', flush=True)
    check(280, f'{a}/{b}')
print(f'done {time.time()-t0:.0f}s mem {mem()}')

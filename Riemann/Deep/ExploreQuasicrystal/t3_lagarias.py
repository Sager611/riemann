# t3_lagarias.py -- explore5 T-C: Hermite-Biehler / Lee-Yang threshold of the differenced functions
#   A_h(s) = F(s+h) + F(s-h)  (Taylor 1945 for h = 1/2; Lagarias math/0601653: all zeros on Re s = 1/2 for |h| >= 1/2)
# Claim tested: A_h has an off-line zero near an off-line zero rho0 of F exactly when h < Re(rho0) - 1/2
# (HB threshold = zero-free half-plane). F = completed Davenport-Heilbronn function (rho0 = 0.80851718246 + 85.69934848538 i);
# control F = xi. L(s,chi) mod 5 via Hurwitz zeta by Euler-Maclaurin in numpy (checked against mpmath below).
# run: OMP_NUM_THREADS=1 uv run --no-project --with mpmath --with numpy --with scipy python3 t3_lagarias.py
import numpy as np, math, time, sys
import mpmath as mp
from scipy.special import loggamma
KAP = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
C1 = (1 - 1j * KAP) / 2; C2 = (1 + 1j * KAP) / 2
CH = [0, 1, 1j, -1j, -1]
B2 = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330, 854513/138, -236364091/2730]
def hurwitz(s, a, N=60, J=12):
    s = np.asarray(s, dtype=complex)
    k = np.arange(N)[:, None] + a
    acc = np.exp(-s[None, :] * np.log(k)).sum(axis=0)
    Na = N + a; lN = math.log(Na)
    acc += np.exp((1 - s) * lN) / (s - 1) + 0.5 * np.exp(-s * lN)
    poch = s.copy(); term_pow = np.exp((-s - 1) * lN)   # (N+a)^{-s-1}
    fact = 2.0
    for j in range(1, J + 1):
        acc += B2[j - 1] / fact * poch * term_pow
        poch = poch * (s + 2 * j - 1) * (s + 2 * j); term_pow = term_pow / (Na * Na)
        fact *= (2 * j + 1) * (2 * j + 2)
    return acc
def L(s, chi):
    s = np.asarray(s, dtype=complex)
    return np.exp(-s * math.log(5)) * sum(chi[a] * hurwitz(s, a / 5) for a in range(1, 5))
def logPhi(s):
    s = np.asarray(s, dtype=complex)
    f = C1 * L(s, CH) + C2 * L(s, [np.conj(c) for c in CH])
    return (s / 2) * math.log(5 / math.pi) + loggamma((s + 1) / 2) + np.log(f + 0j)
def logxi(s):
    s = np.asarray(s, dtype=complex)
    z = hurwitz(s, 1.0)
    return np.log(s * (s - 1) / 2) - (s / 2) * math.log(math.pi) + loggamma(s / 2) + np.log(z + 0j)
# checks against mpmath
mp.mp.dps = 20
z0 = 0.80851718246 + 85.69934848538j
for s in [0.3 + 10j, 0.81 + 85.7j, 1.2 + 91.0j]:
    ref = complex(C1 * mp.dirichlet(s, CH) + C2 * mp.dirichlet(s, [0, 1, -1j, 1j, -1]))
    got = complex(C1 * L(np.array([s]), CH)[0] + C2 * L(np.array([s]), [0, 1, -1j, 1j, -1])[0])
    print("check f_DH(%s): rel err %.1e" % (s, abs(got - ref) / abs(ref)))
print("check |f_DH(rho0)| =", abs(C1 * L(np.array([z0]), CH)[0] + C2 * L(np.array([z0]), [0, 1, -1j, 1j, -1])[0]))
print("check zeta(0.6+88i): rel err %.1e" % (abs(hurwitz(np.array([0.6 + 88j]), 1.0)[0] - complex(mp.zeta(0.6 + 88j))) / abs(complex(mp.zeta(0.6 + 88j)))))
def count(logF, h, s1, s2, T1, T2, step=0.002):
    n = int(round((s2 - s1) / step)); m = int(round((T2 - T1) / step))
    bot = s1 + (s2 - s1) * np.arange(n) / n + 1j * T1
    rig = s2 + 1j * (T1 + (T2 - T1) * np.arange(m) / m)
    top = s2 - (s2 - s1) * np.arange(n) / n + 1j * T2
    lef = s1 + 1j * (T2 - (T2 - T1) * np.arange(m) / m)
    z = np.concatenate([bot, rig, top, lef, [bot[0]]])
    a = logF(z + h); b = logF(z - h); mm = np.maximum(a.real, b.real)
    v = np.exp(a - mm) + np.exp(b - mm)
    ph = np.unwrap(np.angle(v))
    return (ph[-1] - ph[0]) / (2 * math.pi), np.max(np.abs(np.diff(ph))), np.min(np.abs(v))
t0 = time.time()
print("\nDH: zeros of A_h in [1/2+0.003, 1.25] x [80, 92]  (rho0 - 1/2 = 0.30852)")
for h in [0.05, 0.15, 0.25, 0.28, 0.29, 0.30, 0.305, 0.3085, 0.31, 0.32, 0.35, 0.45]:
    c, mx, mn = count(logPhi, h, 0.503, 1.25, 80.0, 92.0)
    print("  h=%.4f  zeros of A_h with Re s > 1/2: %6.3f   (max phase step %.3f, min scaled|A_h| %.2e)   [Phi(s+h) vanishes at Re %.4f]" % (h, c, mx, mn, 0.80851718246 - h))
    sys.stdout.flush()
print("\ncontrol xi: zeros of A_h in [1/2+0.003, 1.25] x [80, 92]")
for h in [0.05, 0.15, 0.25]:
    c, mx, mn = count(logxi, h, 0.503, 1.25, 80.0, 92.0)
    print("  h=%.4f  zeros of A_h with Re s > 1/2: %6.3f   (max phase step %.3f, min scaled|A_h| %.2e)" % (h, c, mx, mn))
print("elapsed %.0f s" % (time.time() - t0))

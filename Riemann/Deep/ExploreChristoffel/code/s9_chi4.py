# explore38 s9: control with an Euler product: L(s, chi_{-4}) (GRH true numerically). NB Gram (digamma, period 4), on-line zeros
# by sign changes of Z_chi, Burnol/Christoffel bound (generic bound() from s7), deficit L^2 (d_N^2 - B) vs L.
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 s9_chi4.py NMAX TZ KON
import sys, math, time, json, resource
import numpy as np
import mpmath as mp
from math import gcd
from scipy.special import digamma
t0 = time.time()
NMAX = int(sys.argv[1]); TZ = float(sys.argv[2]); KON = int(sys.argv[3])
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
src = open(OUT + 's7_dh.py').read()
exec(src[src.index('def bound('):src.index('# ---------------- zeta validation')])
q = 4
S = [0.0, 1.0, 1.0, 0.0]   # A(y) = sum_{n<=y} chi_{-4}(n) = S[floor(y) mod 4]
def g_coprime(h, k):
    P = q * h * k
    br = np.union1d(np.arange(0, P + 1, h), np.arange(0, P + 1, k)).astype(np.float64)
    a, b = br[:-1], br[1:]
    c = np.array(S)[(a // h).astype(np.int64) % q] * np.array(S)[(a // k).astype(np.int64) % q]
    m = c != 0
    return float(np.sum(c[m] / P * (digamma(b[m] / P) - digamma(a[m] / P))))
import os
d2file = OUT + f'chi4_d2_{NMAX}.npy'
if os.path.exists(d2file):
    d2 = np.load(d2file)
    print('loaded', d2file)
G = None if os.path.exists(d2file) else np.zeros((NMAX, NMAX)); cache = {}
if G is not None:
    for h in range(1, NMAX + 1):
        for k in range(h, NMAX + 1):
            g_ = gcd(h, k); a_, b_ = h // g_, k // g_
            if (a_, b_) not in cache:
                cache[(a_, b_)] = g_coprime(a_, b_)
            G[h - 1, k - 1] = G[k - 1, h - 1] = cache[(a_, b_)] / g_
    del cache
    beta = sum(S[r] * (digamma((r + 1) / q) - digamma(r / q)) / q for r in range(1, q))
    bvec = np.array([beta / k for k in range(1, NMAX + 1)])
    Lc = np.linalg.cholesky(G); y = np.linalg.solve(Lc, bvec); d2 = 1 - np.cumsum(y * y)
    del G, Lc
    np.save(d2file, d2)
print(f'chi_-4 Gram N<={NMAX} {time.time()-t0:.1f}s: d2(10)={d2[9]:.6f} d2(10)logN={d2[9]*math.log(10):.5f} d2(N)logN={d2[-1]*math.log(NMAX):.5f}', flush=True)
# zeros
mp.mp.dps = 20
chi = [0, 1, 0, -1]
from scipy.special import loggamma as lgam
from scipy.optimize import brentq
from math import comb
B2 = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330]
def hurwitz(s, a, M):   # Euler-Maclaurin, vectorized over s (numpy complex array)
    n = np.arange(M)[:, None] + a
    tot = np.sum(np.exp(-s[None, :] * np.log(n)), axis=0)
    x = M + a
    tot = tot + x ** (1 - s) / (s - 1) + 0.5 * x ** (-s)
    fac = s.copy(); pw = x ** (-s - 1)
    for k, b in enumerate(B2):
        tot = tot + b / math.factorial(2 * k + 2) * fac * pw
        fac = fac * (s + 2 * k + 1) * (s + 2 * k + 2); pw = pw / (x * x)
    return tot
def Zvec(t):
    t = np.atleast_1d(np.asarray(t, dtype=float)); s = 0.5 + 1j * t
    Lv = np.exp(-s * math.log(4)) * (hurwitz(s, 0.25, 600) - hurwitz(s, 0.75, 600))
    th = (t / 2) * math.log(q / math.pi) + np.imag(lgam(0.75 + 0.5j * t))
    return np.real(np.exp(1j * th) * Lv)
def Zc(t):
    th = (t / 2) * mp.log(mp.mpf(q) / mp.pi) + mp.im(mp.loggamma(mp.mpc(0.75, t / 2)))
    return mp.re(mp.exp(1j * th) * mp.dirichlet(mp.mpc(0.5, t), chi))
for tt in [10.0, 300.0, 850.0]:
    print(f'Z check t={tt}: EM {Zvec(tt)[0]:.12f} mpmath {float(Zc(tt)):.12f}', flush=True)
ts = np.arange(0.05, TZ, 0.02)
vals = np.concatenate([Zvec(ts[i:i + 2000]) for i in range(0, len(ts), 2000)])
idx = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
zs = np.array([brentq(lambda u: Zvec(u)[0], ts[i], ts[i + 1], xtol=1e-13) for i in idx])
print('first zeros vs mpmath refine:', [float(mp.findroot(Zc, mp.mpf(z))) - z for z in zs[:3]], flush=True)
NT = (TZ / (2 * math.pi)) * math.log(q * TZ / (2 * math.pi * math.e)) - 1 / 8
im0 = float(mp.im(mp.dirichlet(mp.mpc(0.5, 10.0), chi) * mp.exp(1j * ((5.0) * mp.log(mp.mpf(q) / mp.pi) + mp.im(mp.loggamma(mp.mpc(0.75, 5.0)))))))
print(f'zeros of L(s,chi_-4) in (0,{TZ}): {len(zs)} by sign changes (smooth count {NT:.1f}); first {zs[:4]}; Im Z check {im0:.2e}; {time.time()-t0:.1f}s', flush=True)
g1 = zs[:KON]
rho = np.concatenate([0.5 + 1j * g1, 0.5 - 1j * g1])
def Vf(w):
    w = mp.mpc(w)
    return complex((w / (1 - w)) * mp.power(mp.mpf(q) / mp.pi, w - 0.5) * mp.gamma((1 + w) / 2) / mp.gamma((2 - w) / 2))
Vv = np.array([Vf(r) for r in rho])
dph = np.array([float(1 / (0.25 + t * t) + math.log(q / math.pi) + mp.re(mp.digamma(mp.mpc(0.75, t / 2)))) for t in np.concatenate([g1, -g1])])
kk = np.arange(0, 20); pd = -(2.0 * kk + 1)
Rd = np.array([float((p_ / (1 - p_)) * mp.power(mp.mpf(q) / mp.pi, p_ - 0.5) * 2 * (-1) ** int(k_) / (mp.factorial(int(k_)) * mp.gamma((2 - p_) / 2))) for k_, p_ in zip(kk, pd)])
td = Rd / pd
gt = zs[KON:]
def tail(L):
    s = 2 * float(np.sum((1 / (0.25 + gt ** 2)) / (L + np.log(q * gt / (2 * math.pi)))))
    u = np.linspace(0, 60, 60001); t = zs[-1] * np.exp(u); lg = np.log(q * t / (2 * math.pi))
    return s + float(np.trapezoid(2 * (lg / (2 * math.pi)) / (t * t * (L + lg)) * t, u))
Con = float(np.sum(1 / np.abs(rho) ** 2)) + 2 * float(np.sum(1 / (0.25 + gt ** 2)))
rows = []
for N in range(5, NMAX + 1):
    L = math.log(N)
    if N in (5, 100, NMAX):
        B, emin = bound(N, rho, Vv, dph, np.ones(len(rho), bool), pd, Rd, 1 / (1 - np.conj(rho)), td, ret_parts=True)
    else:
        B, emin = bound(N, rho, Vv, dph, np.ones(len(rho), bool), pd, Rd, 1 / (1 - np.conj(rho)), td), float('nan')
    tl = tail(L)
    rows.append((N, L, d2[N - 1], B + tl, B + tl))
    if N in (5, 10, 20, 40, 100, 200, 400, 600, 800, 1000) or N == NMAX:
        print(f'chi4 N={N:5d} d2={d2[N-1]:.6e} d2L={d2[N-1]*L:.5f} B+tail={B+tl:.6e} ratio={(B+tl)/d2[N-1]:.4f} L^2(d2-B)={L*L*(d2[N-1]-B-tl):+.5f} eigmin={emin:.3f}', flush=True)
np.save(OUT + 'chi4_compare.npy', np.array(rows))
print(f'sum_on 1/|rho|^2 (zeros < {TZ}) = {Con:.6f}; done {time.time()-t0:.1f}s maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')

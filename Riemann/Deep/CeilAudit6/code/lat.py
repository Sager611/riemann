# ceilaudit6 N12 test. (a) quaternary diagonal forms x^2+a2 y^2+a3 z^2+a4 w^2 (1<=a2<=a3<=a4<=AMAX, a2 a3 a4 a square: trivial
# character): does E_Q(s) have zeta(s)zeta(s-1) as factor, i.e. r_Q(m k) = f(m) sigma_1(k) for k coprime to 2 a2 a3 a4 (all n <= M)?
# T_s = 2 pi sqrt(a2) (Chowla-Selberg along e1; tail K_{s-1/2}(2 pi sqrt(a2) |y|)). Zeta zeros are OFF Re s = 1 for rank 4.
# (b) rank 3: is zeta(2s-1) an analytic factor of E_{Z^3}? (Shimura: formally E = zeta(2s)zeta(2s-1) * sum_t r3(t) t^-s / L(2s,chi_-4t)).
#     Evaluate E_{Z^3} at s = (1+rho_1)/2 = 3/4 + i gamma_1/2 (on the rank-3 critical line) by the theta integral.
import sys, math, numpy as np, mpmath as mp
M = int(sys.argv[1]) if len(sys.argv) > 1 else 3000; AMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20
def r1(a):  # coefficients of theta(a z) up to M
    v = np.zeros(M + 1, dtype=np.int64); x = 0
    while a * x * x <= M:
        v[a * x * x] += 1 if x == 0 else 2; x += 1
    return v
def conv(u, v): return np.convolve(u, v)[:M + 1]
sig = np.zeros(M + 1, dtype=np.int64)
for d in range(1, M + 1): sig[d::d] += d
base = r1(1); cache = {}
def isq(n): r = int(math.isqrt(n)); return r * r == n
res = []
for a2 in range(1, AMAX + 1):
    for a3 in range(a2, AMAX + 1):
        for a4 in range(a3, AMAX + 1):
            if not isq(a2 * a3 * a4): continue
            r = conv(conv(conv(base, r1(a2)), r1(a3)), r1(a4))
            Nl = 2 * a2 * a3 * a4
            ks = [k for k in range(1, M + 1) if math.gcd(k, Nl) == 1]
            ratio = r[ks] / sig[ks]
            zf = bool(np.all(np.abs(ratio - ratio[0]) < 1e-9))
            # multiplicativity check for m = 2, 4, 3 (if coprime part ok)
            Ts = 2 * math.pi * math.sqrt(a2)
            res.append((a2, a3, a4, zf, Ts, float(ratio[0]), float(np.std(ratio))))
zfs = [x for x in res if x[3]]
print(f'(a) {len(res)} forms with square a2 a3 a4 <= {AMAX}^3, n <= {M}; zeta(s)zeta(s-1)-factor (r = c sigma_1 on k coprime to level): {len(zfs)}')
for x in zfs: print(f'   ZETA-FACTOR a=(1,{x[0]},{x[1]},{x[2]}) T_s=2pi sqrt(a2)={x[4]:.2f} c={x[5]:.4f}')
big = [x for x in res if x[4] > 14.134725]
print(f'   forms with T_s > gamma_1 = 14.1347: {len(big)}, of which zeta-factor: {sum(1 for x in big if x[3])}; min relative sd of r/sigma among them: {min((x[6]/abs(x[5]) for x in big), default=float("nan")):.3f}')
# (b)
mp.mp.dps = 30
def th3m1(t):  # theta(t)^3 - 1, theta(t) = sum_n exp(-pi n^2 t)
    th = 1 + 2 * mp.nsum(lambda n: mp.exp(-mp.pi * n * n * t), [1, mp.inf]); return th ** 3 - 1
def E3(s):
    I = mp.quad(lambda t: th3m1(t) * (t ** (s - 1) + t ** (mp.mpf(3) / 2 - s - 1)), [1, 4, 16, mp.inf])
    return (I - 1 / s - 1 / (mp.mpf(3) / 2 - s)) * mp.pi ** s / mp.gamma(s)
# check: E3(2) vs direct partial sum is slow; check against known E_{Z^3}(2) ~ 10.1949 (Madelung-type constant check below)
g1 = mp.zetazero(1).imag
s0 = mp.mpf(3) / 4 + 1j * g1 / 2
v = E3(s0); z = mp.zeta(2 * s0 - 1)
print(f'(b) E_Z3(3/4 + i gamma_1/2) = {mp.nstr(v, 10)}  |E|={mp.nstr(abs(v), 6)};  zeta(2s-1) there = {mp.nstr(z, 5)}')
for t in [6.5, 7.0, 7.5]:
    print(f'    |E_Z3(3/4 + {t} i)| = {mp.nstr(abs(E3(mp.mpf(3)/4 + 1j*t)), 6)}')
v2 = E3(mp.mpf(2)); brute = mp.mpf(sum(int(r) * mp.mpf(n) ** -2 for n, r in enumerate(conv(conv(base, base), base)) if n > 0))
print(f'    sanity: E_Z3(2) theta-integral {mp.nstr(v2, 10)} vs partial sum n<={M} {mp.nstr(brute, 10)} (tail ~ 4 pi/sqrt(M) = {4*math.pi/math.sqrt(M):.3f})')

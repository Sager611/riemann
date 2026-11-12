# Exact second moment S(D) = sum_{chi != 1} L(1/2,chi)^2 over class group characters of fundamental D < -4 (w = 2),
# via S = h sum_C zeta(1/2,C)^2 - zeta_K(1/2)^2, zeta(1/2,C) = Z_Q(1/2)/2, Chowla-Selberg at s = 1/2:
# Z_Q(1/2) = a^{-1/2} [2 gamma + 2 log(sqrt|D|/(8 pi a)) + 8 sum_n d(n) cos(pi n b/a) K_0(pi n sqrt|D|/a)].
import numpy as np, math, sys, json, time, resource
from scipy.special import k0
def fund(D):
    n = -D
    if D % 4 == 1 % 4 and D % 4 == 1:
        pass
    best = 1
    for g in range(1, math.isqrt(n) + 1):
        if n % (g*g) == 0 and (D // (g*g)) % 4 in (0, 1): best = g
    return best == 1
def forms(D):
    n = -D; A = []; B = []
    for a in range(1, math.isqrt(n // 3) + 1):
        b = np.arange(-a + 1, a + 1, dtype=np.int64)
        num = b*b - D
        ok = num % (4*a) == 0
        b = b[ok]; c = num[ok] // (4*a)
        ok = (c >= a) & ~((b < 0) & (c == a))
        b = b[ok]; c = c[ok]
        g = np.gcd(np.gcd(a, np.abs(b)), c)
        b = b[g == 1]
        A.append(np.full(len(b), a, dtype=np.int64)); B.append(b)
    return np.concatenate(A), np.concatenate(B)
dn = np.array([0] + [sum(1 for d in range(1, n+1) if n % d == 0) for n in range(1, 13)], dtype=float)
def S_of(D):
    a, b = forms(D); h = len(a); n = -D; r = math.sqrt(n)
    z = 2*np.euler_gamma + 2*np.log(r/(8*np.pi*a))
    for m in range(1, 13):
        z = z + 8*dn[m]*np.cos(np.pi*m*b/a)*k0(np.pi*m*r/a)
    z = z/np.sqrt(a)/2.0
    return h, float(h*(z*z).sum() - z.sum()**2), float((1/a).sum())
if __name__ == '__main__':
    out = sys.argv[1]; Ds = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else None
    t0 = time.time(); rows = []
    for D in Ds:
        h, S, ia = S_of(D); n = -D
        rows.append(dict(D=D, h=h, S=S, R=S/math.sqrt(n), L1=math.pi*h/math.sqrt(n), sum_inv_a=ia))
        if time.time() - t0 > 1500: print('TIME CAP'); break
    json.dump(rows, open(out, 'w'))
    for x in rows[:6]: print(x)
    print('n', len(rows), 'time', round(time.time()-t0, 1), 'maxrss MB', round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20, 1))

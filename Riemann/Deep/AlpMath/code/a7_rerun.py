# A7 rerun: ClassMean (Chowla-Selberg) with exact real parts of tau (first run used a float -b/(2a): error floor 8e-20)
import sys, math, mpmath as mp
from sympy import primerange
mp.mp.dps = 100
def eta(tau):
    q = mp.exp(2j * mp.pi * tau); return mp.exp(2j * mp.pi * tau / 24) * mp.qp(q)
def kron_m(p, t):  # Kronecker (-p / t) for p = 3 mod 4 prime, t >= 1: equals Legendre (t/p)
    r = pow(t % p, (p - 1) // 2, p); return 0 if r == 0 else (1 if r == 1 else -1)
err = 0; n = 0
for p in [p for p in primerange(3, 200) if p % 4 == 3]:
    F = []; a = 1
    while 3 * a * a <= p:
        for b in range(-a + 1, a + 1):
            if (b * b + p) % (4 * a): continue
            c = (b * b + p) // (4 * a)
            if c < a or (c == a and b < 0) or math.gcd(math.gcd(a, abs(b)), c) != 1: continue
            F.append((a, b, c))
        a += 1
    mean = mp.fsum(mp.log(4 * mp.pi * (mp.sqrt(p) / (2 * a)) * abs(eta(mp.mpc(mp.mpf(-b) / (2 * a), mp.sqrt(p) / (2 * a)))) ** 4) for (a, b, c) in F) / len(F)
    s = 1 + mp.mpf(10) ** -30
    chi = [kron_m(p, t) for t in range(1, p + 1)]
    L = p ** (-s) * mp.fsum(chi[t - 1] * mp.zeta(s, mp.mpf(t) / p) for t in range(1, p + 1))
    Ld = -mp.log(p) * L + p ** (-s) * mp.fsum(chi[t - 1] * mp.zeta(s, mp.mpf(t) / p, 1) for t in range(1, p + 1))
    rhs = mp.euler + mp.log(2 * mp.pi / mp.sqrt(p)) - Ld / L
    e = abs(mean - rhs); err = max(err, e); n += 1
    print(p, len(F), mp.nstr(e, 3)); sys.stdout.flush()
print(f"A7 ClassMean rerun: {'PASS' if err < 1e-25 else 'FAIL'} {n} primes p = 3 mod 4, p < 200, dps 100, s-1 = 1e-30, max |err| {mp.nstr(err, 3)}")

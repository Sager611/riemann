# R323 (node N10): a prime in (x, x + x^(7/8)] for every real x >= 1 (routine corollary; see REPORT.md)
import mpmath as mp
mp.mp.dps = 40
A, B, C = mp.mpf('0.04934'), mp.mpf('0.02636425'), mp.mpf('1.4262')
def E(z):  # Fiori-Jaskari 2609.23222 Thm 2 (given Re rho <= 7/8), valid for z > e^3
    return A*z**mp.mpf('0.75')*mp.log(z)**2 + B*z**mp.mpf('0.875')
def r(L):  # [E(x) + E(x+h) + C sqrt(x+h)] / h with x = e^L, h = x^(7/8)
    x = mp.e**L; h = x**mp.mpf('0.875')
    return (E(x) + E(x+h) + C*mp.sqrt(x+h))/h
lo, hi = mp.mpf(16), mp.mpf(60)
for _ in range(120):
    m = (lo+hi)/2
    if r(m) < 1: hi = m
    else: lo = m
print('r = 1 at log x =', mp.nstr(hi, 15), ' x =', mp.nstr(mp.e**hi, 10))
Lc = mp.mpf('41.56')
print('certified cutoff: r(41.56) =', mp.nstr(r(Lc), 15), '< 1;  e^41.56 =', mp.nstr(mp.e**Lc, 10), '< 4e18')
print('limit of r as x -> oo: 2B =', mp.nstr(2*B, 10))
# small x: worst case x = p_n; need p_{n+1} - p_n <= p_n^(7/8) for p_n < 4187 (1476^(8/7) = 4186.13); checked to 10^7
N = 10**7
s = bytearray([1])*(N+1); s[0] = s[1] = 0
for i in range(2, int(N**.5)+1):
    if s[i]: s[i*i::i] = bytearray(len(range(i*i, N+1, i)))
ps = [i for i in range(N+1) if s[i]]
bad = [(p, q) for p, q in zip(ps, ps[1:]) if q-p > p**0.875]
print('primes below 1e7:', len(ps), ' gaps with p_{n+1}-p_n > p_n^(7/8):', len(bad))
print('1476^(8/7) =', 1476**(8/7))

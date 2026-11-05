"""Z.-W. Sun, arXiv 2607.07638 (v3, 21 Jul 2026), Conjecture 4.1(i): for an odd prime p with (p/3) = 1 or (p/13) = 1,
sum_{k=0}^{p-1} D(k) u_k == 0 (mod p^2), D = Domb numbers, u_0 = 0, u_1 = 1, u_{n+1} = 11 u_n - u_{n-1}.
Fails at p = 3 ((3/13) = 1): 1*0 + 4*1 + 28*11 = 312 == 6 (mod 9). No other failure for p < 400 (also the mod p^3 and (ii) parts).
Run: python3 verify_2607.py (stdlib only)."""
from math import comb
def D(n): return sum(comb(n, k) ** 2 * comb(2 * k, k) * comb(2 * (n - k), n - k) for k in range(n + 1))
def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
def primes(n): return [q for q in range(3, n) if all(q % d for d in range(2, int(q ** 0.5) + 1))]
N = 400; Dv = [D(n) for n in range(N)]; u = [0, 1]
for n in range(1, N): u.append(11 * u[n] - u[n - 1])
assert (Dv[:3], u[:3]) == ([1, 4, 28], [0, 1, 11]) and leg(3, 13) == 1
s3 = sum(Dv[k] * u[k] for k in range(3)); assert s3 == 312 and s3 % 9 == 6
bad = [p for p in primes(N) if (leg(p, 3) == 1 or leg(p, 13) == 1) and sum(Dv[k] * u[k] for k in range(p)) % (p * p)]
assert bad == [3]
print("Conj. (i) fails exactly at p = 3 among odd primes p < 400: sum = 312 == 6 (mod 9)")


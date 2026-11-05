# Numerical sanity check of Lemmas 2-3 (single weight, B = 1/2 bounds; GRH holds numerically for these small conductors).
# Run: uv run --no-project --with mpmath --with numpy python3 -u sanity.py   (~8 min; output: sanity_output.txt)

# Sanity check of the smoothed explicit formula (EF1),(EF2) on small conductors, numerically (floats/mpmath).
import numpy as np, mpmath as mp, math
mp.mp.dps = 30
def kron(D, n):
    # Kronecker symbol (D/n) for fundamental discriminant D, n>=1
    return int(mp.mpf(0)) if False else _kron(D, n)
def _kron(a, n):
    # Kronecker symbol (a/n), n>0
    if n == 0: return 1 if abs(a) == 1 else 0
    res = 1
    # factor out 2s
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0: return 0
        if a % 8 in (3, 5): res = -res
    # Jacobi (a/n) with n odd
    a = a % n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: res = -res
        a = a % n
    return res if n == 1 else 0

def test(D, L, h, Bz):
    q = abs(D); a = 1 if D < 0 else 0
    chi = [_kron(D, n) for n in range(q)]
    A = -mp.fsum(chi[r]*mp.digamma(mp.mpf(r)/q) for r in range(1, q) if chi[r])
    C = mp.fsum(chi[r]*mp.stieltjes(1, mp.mpf(r)/q) for r in range(1, q) if chi[r])
    Lv = A/q
    logL = float(mp.log(Lv)); LdL = float(-mp.log(q) - C/A)
    if q <= 12:
        Lv2 = mp.dirichlet(1, chi); dLv2 = mp.dirichlet(1, chi, 1)
        print('  check digamma/stieltjes vs mp.dirichlet:', float(Lv-Lv2), float(LdL - dLv2/Lv2), flush=True)
    N = int(math.exp(L*(1+h))) + 1
    # Lambda(n) via sieve
    lam = np.zeros(N+1)
    isp = np.ones(N+1, dtype=bool); isp[:2] = False
    for p in range(2, int(N**0.5)+1):
        if isp[p]: isp[p*p::p] = False
    primes = np.nonzero(isp)[0]
    P0 = 0.0; P1 = 0.0
    for p in primes:
        lp = math.log(p); pk = p; k = 1
        c = chi[p % q]
        while pk <= N:
            u = k*lp/L
            v = 1.0 if u <= 1 else (max(0.0, (1+h-u)/h))
            if v > 0:
                ck = c**k
                P0 += lp*ck*v/(pk*k*lp)
                P1 += lp*ck*v/pk
            pk *= p; k += 1
    b = 1 - Bz
    e1 = math.exp(-b*L); e2 = math.exp(-(1+h)*b*L)
    k0 = (e2/(1+h) + e1)/(h*L*L*b); k1 = (e2 + e1)/(h*L*b)
    psi = float(mp.digamma((1+a)/2.0))
    Z = LdL + 0.5*math.log(q/math.pi) + 0.5*psi
    print(flush=True); print(f"D={D} L={L} h={h}: logL={logL:.6f} P0={P0:.6f} diff={logL-P0:+.2e} bound(B={Bz})={k0*Z:.2e} | -L'/L={-LdL:.6f} P1={P1:.6f} diff={-LdL-P1:+.2e} bound={k1*Z:.2e} | Z={Z:.4f}")

for D in [5, -4, -3, 8, -7, 12, -23, 13]:
    test(D, 10.0, 0.5, 0.5)

"""C2: cubic basin of the Hermite polynomial in finite free cumulant coordinates.
p_eps = exp(-y^2/2 + eps y^3)(D) x^d, i.e. sum_j c_j (d)_j x^{d-j} with sum c_j y^j = exp(-y^2/2 + eps y^3) mod y^{d+1}.
r3(d) = sup{eps > 0 : p_eps real-rooted} by bisection over exact rationals (flint exact root isolation).
Necessary condition (centred root moments, Hankel 3x3): eps^2 <= 1/(9(d-2)), i.e. r3(d) <= 1/(3 sqrt(d-2)).
Also the line through the actual cumulant vector of xi's J*^{d,n} (all rho_k, k <= d): t* = sup{t: p(t*rho) real}."""
from fractions import Fraction as Fr
import flint, math, sys
def coeffs_from_logsym(l, d):   # l: dict k->Fraction (log-symbol), returns c_0..c_d of exp(sum l_k y^k)
    c = [Fr(0)]*(d+1); c[0] = Fr(1)
    # exp via recurrence: j c_j = sum_{k=1}^{j} k l_k c_{j-k}
    for j in range(1, d+1):
        s = Fr(0)
        for k in range(1, j+1):
            if k in l: s += k*l[k]*c[j-k]
        c[j] = s/j
    return c
def is_real_rooted(c, d):
    co = [Fr(0)]*(d+1)          # ascending powers of x: coefficient of x^{d-j} is c_j (d)_j
    for j in range(d+1):
        co[d-j] = c[j]*math.perm(d, j)
    den = 1
    for v in co: den = den*v.denominator//math.gcd(den, v.denominator)
    p = flint.fmpz_poly([int(v*den) for v in co])
    return sum(m for r, m in p.complex_roots() if r.imag == 0) == d
def r3(d, steps=40):
    lo, hi = Fr(0), Fr(1)
    for _ in range(steps):
        mid = (lo+hi)/2
        if is_real_rooted(coeffs_from_logsym({2: Fr(-1, 2), 3: mid}, d), d): lo = mid
        else: hi = mid
    return lo, hi
if __name__ == '__main__':
    D = int(sys.argv[1])
    print(' d   r3(d) in [lo,hi]          Hankel bound 1/(3sqrt(d-2))   r3*sqrt(d)')
    for d in list(range(3, 13)) + list(range(15, D+1, 5)):
        lo, hi = r3(d)
        hb = 1/(3*math.sqrt(d-2)) if d > 2 else float('inf')
        print('%3d  %.10f %.10f   %.10f   %.6f' % (d, float(lo), float(hi), hb, float(lo)*math.sqrt(d)))

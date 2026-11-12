# explore47 N1 measurement: Toeplitz multiplier σ_c of Weil's form on [-a,a], a = ln(c)/2: Q(f) = (1/2π)∫|f̂|²σ_c, so λ_min ≥ inf σ_c
# (plain symbol / anti-Wick-type bound). σ_c = Reψ(1/4+iξ/2) − log π + FT[2cosh(s/2)1_{|s|≤2a}](ξ) − Σ_{n<c} 2Λ(n)n^{-1/2}cos(ξ log n).
import math, numpy as np
from scipy.special import psi
def lam(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0
xi = np.linspace(0, 400, 400001)
arch = np.real(psi(0.25 + 0.5j * xi)) - math.log(math.pi)
print('arch < 0 for ξ <', xi[np.nonzero(arch >= 0)[0][0]])
for c in [1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 47.0]:
    a = math.log(c) / 2
    pole = 4 * (0.5 * math.sinh(a) * np.cos(2 * a * xi) + xi * math.cosh(a) * np.sin(2 * a * xi)) / (xi ** 2 + 0.25)
    pr = np.zeros_like(xi)
    for n in range(2, int(math.ceil(c))):
        L = lam(n)
        if L and n < c: pr += 2 * L / math.sqrt(n) * np.cos(xi * math.log(n))
    s = arch + pole - pr; neg = s < 0
    meas = 2 * np.sum(neg) * (xi[1] - xi[0])          # both signs of ξ
    cells = 2 * a * meas / (2 * math.pi)              # phase-space area / 2π
    print(f'c={c:5.1f} inf σ_c={s.min():8.3f} at ξ={xi[np.argmin(s)]:7.3f}  |{{σ_c<0}}|={meas:7.3f}  cells=(2a·|·|)/2π={cells:6.2f}  max|primes|={2*sum(lam(n)/math.sqrt(n) for n in range(2,int(math.ceil(c))) if n<c):.2f}')

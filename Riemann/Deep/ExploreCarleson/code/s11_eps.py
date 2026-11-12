# explore42: eps_N = ||(1-Q_lam) V chi_1||^2 = int_0^lam |V chi_1|^2, V chi_1 = sum_j b_j (2 pi t)^{2j} (explore38 DERIVATION s2)
import math, sys
import mpmath as mp
mp.mp.dps = 40
def eps(N, J=60):
    lam = mp.mpf(1) / N
    b = [2 * (-1) ** j * (mp.mpf(2 * j) / (2 * j + 1)) / mp.factorial(2 * j + 1) for j in range(1, J + 1)]
    s = mp.mpf(0)
    for j in range(1, J + 1):
        for k in range(1, J + 1):
            e = 2 * j + 2 * k
            s += b[j - 1] * b[k - 1] * (2 * mp.pi) ** e * lam ** (e + 1) / (e + 1)
    return s
for N in [2, 3, 5, 10, 20, 40, 50, 100, 1280]:
    e = eps(N)
    print(f'N={N:5d} eps={mp.nstr(e, 12)} eps*logN={mp.nstr(e * mp.log(N), 12)}')

# bvirraudit: exact size of the term the discovery page's Cor 3 wording omits (roots of unity of ORDER <= N instead of DEGREE < N).
# Littlewood model f = sum_{k=0}^n a_k x^k, a_k iid +-1. Phi_4 | f  <=>  f(i) = 0  <=>  two independent sign sums vanish.
import math
def p0(m):  # P(sum of m iid signs = 0)
    return math.comb(m, m // 2) / 2 ** m if m % 2 == 0 else 0.0
for n in [15, 99, 999, 9999, 99999]:
    if n % 4 != 3: continue
    E = n // 2 + 1; O = (n + 1) // 2     # number of even / odd indices in 0..n
    P4 = p0(E) * p0(O)
    print("n=%6d  P(Phi_4 | f) = %.4e   4/(pi n) = %.4e   n^(-3/2) = %.4e   ratio = %.1f" % (n, P4, 4 / (math.pi * n), n ** -1.5, P4 / n ** -1.5))

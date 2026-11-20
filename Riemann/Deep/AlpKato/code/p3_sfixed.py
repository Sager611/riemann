# P3 (alpkato): the hypothesis hS : S = Above L l of LData is load-bearing for LemmaK.
# Localization (L totally complex, i = e+1 >= 2, S' = S_l u T, T finite primes not over l):
#   0 -> H^2(O_L[1/l], Z_l(i)) -> H^2(O_{L,S'}, Z_l(i)) -> (+)_{v in T} Z_l/(Nv^(i-1) - 1) -> 0   (H^3 = 0: no real places)
# so log #H^2(G_{L,S'}, Z_l(i)) >= sum_{v in T} v_l(Nv^e - 1) log l, unbounded in T, while kBound depends only on (e, n, D).
import sympy as sp, math
def vl(x, l):
    c = 0
    while x % l == 0: x //= l; c += 1
    return c
def norms_Qi(X):      # primes of Z[i] not over 2: p = 1 mod 4 gives two primes of norm p; p = 3 mod 4 one prime of norm p^2
    for p in sp.primerange(3, X+1):
        if p % 4 == 1: yield p; yield p
        elif p*p <= X: yield p*p
def norms_Qz3(X):     # primes of Z[zeta_3] not over 3: p = 1 mod 3 two of norm p; p = 2 mod 3 one of norm p^2
    for p in sp.primerange(2, X+1):
        if p == 3: continue
        if p % 3 == 1: yield p; yield p
        elif p*p <= X: yield p*p
for name, gen, l in [("L = Q(i), l = 2", norms_Qi, 2), ("L = Q(zeta_3), l = 3", norms_Qz3, 3)]:
    for e in [1, 2]:
        row = []
        for X in [10**2, 10**3, 10**4, 10**5, 10**6]:
            T = list(gen(X)); s = sum(vl(N**e - 1, l) for N in T)
            row.append("X=%g: #T=%d, sum v_l=%d" % (X, len(T), s))
        print("%s, e = %d (Z_l(%d)):  " % (name, e, e+1) + "; ".join(row))
print("With S = S_l: H^2(Z[i,1/2], Z_2(2)) = K_2(Z[i]) (x) Z_2 = 0 and H^2(Z[zeta_3,1/3], Z_3(2)) = K_2(Z[zeta_3]) (x) Z_3 = 0 (Tate 1976), consistent with LemmaK.")
print("kBound(24, e=1, n=2, D=4) = exp(4^24 * log 8) = exp(%.4g): with S arbitrary, log #H^2 >= (sum v_l) log l exceeds it once #T is of order exp(%.4g); finite but unbounded in T, so the uniform statement would be false" % (4**24*math.log(8), 4**24*math.log(8)))

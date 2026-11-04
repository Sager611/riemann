"""Saxl's statement g(rho_m, rho_m, lam) > 0 for every lam |- N_m, computed from
characters (own Murnaghan-Nakayama, mn.py), exact integer arithmetic.
usage: python3 saxl_kron.py MMAX"""
import sys
from math import factorial
import mn

def kron_staircase(m):
    N = m*(m+1)//2
    rho = tuple(range(m, 0, -1))
    P = list(mn.partitions(N))
    acc = {lam: 0 for lam in P}
    for mu in P:
        mn._mn.cache_clear()
        w = mn.class_size(mu) * mn.chi(rho, mu)**2
        if w == 0: continue
        for lam in P:
            acc[lam] += w * mn.chi(lam, mu)
    g = {}
    for lam in P:
        q, r = divmod(acc[lam], factorial(N))
        assert r == 0
        g[lam] = q
    return N, g

if __name__ == "__main__":
    M = int(sys.argv[1])
    for m in range(1, M+1):
        N, g = kron_staircase(m)
        zero = [l for l, v in g.items() if v <= 0]
        mn_ = min(g.values()); mx = max(g.values())
        # consistency: sum_lam g * dim(lam) = dim(rho)^2
        rho = tuple(range(m, 0, -1))
        assert sum(v*mn.dim(l) for l, v in g.items()) == mn.dim(rho)**2
        print(f"m={m} N={N} #lam={len(g)} min g={mn_} max g={mx} zero-or-negative={zero}")

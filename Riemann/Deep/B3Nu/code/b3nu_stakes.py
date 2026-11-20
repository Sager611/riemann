# b3nu: stakes of further dual improvements with the frozen detector (P, Q, R = 44/5): theta at which kappa(nu) = k*(theta),
# for fixed nu, and for the moment range nu = theta/2 (limit of Lemma C's hypotheses) and nu = (3theta-1)/4 + g(1-theta) for gain factors g.
import sys
src = open(__file__.replace('b3nu_stakes.py', 'b3c_arb_copy.py')).read().split('# admissibility')[0].replace('ctx.prec = 4000', 'ctx.prec = 3000')
exec(src)
s2 = arb(2).sqrt(); cth = lambda th: 2 - th/2 - (th/s2).cot()/s2; kst = lambda th: 1 - 2/(2 - cth(th)); R = fmpq(44, 5)
def solve(fn, lo=fmpq(40, 100), hi=fmpq(55, 100)):
    if not (fn(lo) < 0 and fn(hi) > 0): return None
    for _ in range(26):
        mid = (lo + hi)/2
        if fn(mid) > 0: hi = mid
        else: lo = mid
    return float(hi)
for nu in [fmpq(13, 100), fmpq(15, 100), fmpq(20, 100), fmpq(1, 4)]:
    k = kappa(nu, R)[0]; print(f"fixed nu = {float(nu):.3f}: kappa = {k.str(8)}, theta_det = {solve(lambda th: k - kst(arb(th)))}", flush=True)
print("nu = theta/2:", solve(lambda th: kappa(th/2, R)[0] - kst(arb(th))), flush=True)
for g in [fmpq(1, 100), fmpq(1, 20), fmpq(1, 10)]:
    print(f"nu = (3th-1)/4 + {float(g)}(1-th):", solve(lambda th: kappa((3*th - 1)/4 + g*(1 - th), R)[0] - kst(arb(th))), flush=True)

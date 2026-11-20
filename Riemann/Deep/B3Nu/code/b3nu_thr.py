# b3nu: simple-zero threshold with the frozen SL/EXP-010 detector (construction = B3Cert/code/b3c_arb.py lines up to kappa, executed verbatim
# from the read-only copy b3c_arb_copy.py), Wang's c(theta), k*(theta), and the Lemma C moment range nu < nu1(theta).
# usage: python3 b3nu_thr.py <params_exp010.json> <prec bits> ; float scan + bisection for theta1: kappa(nu1(theta), R) = k*(theta).
import sys
src = open(__file__.replace('b3nu_thr.py', 'b3c_arb_copy.py')).read().split('# admissibility')[0]
prec = int(sys.argv[2]); src = src.replace('ctx.prec = 4000', f'ctx.prec = {prec}')
exec(src)   # defines fmpq, arb, P, Qp, kappa(nu, R), ...
s2 = arb(2).sqrt(); cth = lambda th: 2 - th/2 - (th/s2).cot()/s2; kst = lambda th: 1 - 2/(2 - cth(th))
h = lambda th, kk: (3 + kk - ((1 - kk)*(9 - kk - 8*cth(th))).sqrt())/4
R = fmpq(44, 5)
nu1 = lambda th, k, l: (th*(1 + 2*l) - (k + l))/(2*(1 + l - k))
k0, _ = kappa(fmpq(121, 1000), R); k1, _ = kappa(fmpq(122, 1000), R)
print("kappa(0.121) =", k0.str(15), " dkappa/dnu ~", ((k1 - k0)*1000).str(6), " dk*/dtheta ~", ((kst(arb(fmpq(4948,10000))) - kst(arb(fmpq(4947,10000))))*10000).str(6), flush=True)
for name, (k, l) in [("(0,1) = Cor B'", (fmpq(0), fmpq(1))), ("A^3B(0,1) = (1/30,13/15)", (fmpq(1, 30), fmpq(13, 15))), ("A^2(Bourgain) = (13/414,359/414)", (fmpq(13, 414), fmpq(359, 414))), ("A^3BA^2BA^2B(0,1) = (13/502,445/502)", (fmpq(13, 502), fmpq(445, 502)))]:
    lo, hi = fmpq(48, 100), fmpq(50, 100)
    f = lambda th: kappa(nu1(th, k, l), R)[0] - kst(arb(th))
    assert f(lo) < 0 and f(hi) > 0
    for _ in range(30):
        mid = (lo + hi)/2
        if f(mid) > 0: hi = mid
        else: lo = mid
    print(f"{name}: theta1 in [{float(lo):.9f}, {float(hi):.9f}], nu1(theta1) = {float(nu1(hi, k, l)):.6f}, kappa = {kappa(nu1(hi, k, l), R)[0].str(10)}", flush=True)

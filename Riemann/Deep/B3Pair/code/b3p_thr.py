# b3pair (copy of B3Nu/code/b3nu_thr.py; only the pair list changed): simple-zero threshold with the frozen SL/EXP-010 detector (construction = B3Cert/code/b3c_arb.py lines up to kappa, executed verbatim
# from the read-only copy b3c_arb_copy.py), Wang's c(theta), k*(theta), and the Lemma C moment range nu < nu1(theta).
# usage: python3 b3p_thr.py <params_exp010.json> <prec bits> ; float scan + bisection for theta1: kappa(nu1(theta), R) = k*(theta).
import sys
src = open(__file__.replace('b3p_thr.py', 'b3c_arb_copy.py')).read().split('# admissibility')[0]
prec = int(sys.argv[2]); src = src.replace('ctx.prec = 4000', f'ctx.prec = {prec}')
exec(src)   # defines fmpq, arb, P, Qp, kappa(nu, R), ...
s2 = arb(2).sqrt(); cth = lambda th: 2 - th/2 - (th/s2).cot()/s2; kst = lambda th: 1 - 2/(2 - cth(th))
h = lambda th, kk: (3 + kk - ((1 - kk)*(9 - kk - 8*cth(th))).sqrt())/4
R = fmpq(44, 5)
nu1 = lambda th, k, l: (th*(1 + 2*l) - (k + l))/(2*(1 + l - k))
k0, _ = kappa(fmpq(121, 1000), R); k1, _ = kappa(fmpq(122, 1000), R)
print("kappa(0.121) =", k0.str(15), " dkappa/dnu ~", ((k1 - k0)*1000).str(6), " dk*/dtheta ~", ((kst(arb(fmpq(4948,10000))) - kst(arb(fmpq(4947,10000))))*10000).str(6), flush=True)
def Aq(k, l): return (k/(2*k + 2), (k + l + 1)/(2*k + 2))
e6 = fmpq(1, 10**6)
PAIRS = [("control A^3BA^2BA^2B(0,1) = (13/502,445/502)", (fmpq(13, 502), fmpq(445, 502))),
         ("A^2(Bourgain) exact = (13/414,359/414)", Aq(*Aq(fmpq(13, 84), fmpq(55, 84)))),
         ("A^2(Bourgain + 1e-6)", Aq(*Aq(fmpq(13, 84) + e6, fmpq(55, 84) + e6))),
         ("A(TTY 89/1282,997/1282) = (89/2742,1184/1371)", Aq(fmpq(89, 1282), fmpq(997, 1282))),
         ("A(TTY 89/1282 + 1e-6, 997/1282 + 1e-6)", Aq(fmpq(89, 1282) + e6, fmpq(997, 1282) + e6)),
         ("A(TTY 652397/9713986,7599781/9713986)", Aq(fmpq(652397, 9713986), fmpq(7599781, 9713986))),
         ("A^2(Huxley 32/205,269/410)", Aq(*Aq(fmpq(32, 205), fmpq(269, 410)))),
         ("A^2(Huxley + 1e-6)", Aq(*Aq(fmpq(32, 205) + e6, fmpq(269, 410) + e6))),
         ("stake: EPC A(0,1/2) = (0,3/4)", Aq(fmpq(0), fmpq(1, 2))),
         ("stake: EPC (0,1/2) on F_u directly", (fmpq(0), fmpq(1, 2))),
         ("R515 check: outer A^2(Bourgain) + 1e-6", (fmpq(13, 414) + e6, fmpq(359, 414) + e6))]
if len(sys.argv) > 3: PAIRS = [q for q in PAIRS if q[0].startswith(sys.argv[3])]
for name, (k, l) in PAIRS: print(name, '=', k, l, flush=True)
for name, (k, l) in PAIRS:
    lo, hi = fmpq(48, 100), fmpq(50, 100)
    f = lambda th: kappa(nu1(th, k, l), R)[0] - kst(arb(th))
    assert f(lo) < 0 and f(hi) > 0
    for _ in range(30):
        mid = (lo + hi)/2
        if f(mid) > 0: hi = mid
        else: lo = mid
    print(f"{name}: theta1 in [{float(lo):.9f}, {float(hi):.9f}], nu1(theta1) = {float(nu1(hi, k, l)):.6f}, kappa = {kappa(nu1(hi, k, l), R)[0].str(10)}", flush=True)

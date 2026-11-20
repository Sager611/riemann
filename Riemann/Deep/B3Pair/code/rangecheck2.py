# b3pair range check, row R16b (E1 l.67: Thm 5.1 needs kappa < 1 and c(theta) < 1): c(theta) = 2 - theta/2 - cot(theta/sqrt2)/sqrt2 (b3c_arb_copy/b3nu_thr 'cth'), Arb.
from flint import arb, fmpq, ctx
ctx.prec = 200
s2 = arb(2).sqrt(); c = lambda th: 2 - th/2 - (th/s2).cot()/s2
for th in (fmpq(4939, 10000), fmpq(1, 2), fmpq(9, 10), fmpq(1)):
    print('theta =', th, ': c(theta) =', c(arb(th)).str(12), '; 1 - c =', (1 - c(arb(th))).str(12), '(>0:', (1 - c(arb(th))) > 0, ')')
print("c' = cot^2(theta/sqrt2)/2 > 0 on (0,1) (B3F l.177), so c(theta) < c(1) < 1 on (theta_c, 1); kappa(nu0, 44/5) = 0.0861980736963448913832780779143 +/- 2.23e-32 < 1 (logs/b3p_arb_tty_4000.log)")

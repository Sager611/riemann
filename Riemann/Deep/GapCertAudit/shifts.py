"""gapcertaudit: what the two shifts in cert4.py do (Lemma 2 quantities, exact rationals, mpmath 40 digits).
Forms: 'cert' (certified), 'noshift' (a0 - 2e-5), 'sdp' (a0 - 2e-5, phi12 - 2e-6 = the SDP witness after rationalisation).
Prints p(1), n1..n3, and S_odd(g), S_even(g) at g = pi and at small g (sum over k <= KMAX plus the tail bound/estimate)."""
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40
A = [Fr('-5899179108889603/6134092019850000'), Fr(1), Fr('14789027537/948575124080'), Fr('-4037122085/176350902122'),
     Fr('1982999731/771179571628')]
PH12 = Fr('3691572375143311/104390384155500000')
FORMS = {'cert': (A, PH12), 'noshift': ([A[0] - Fr(2, 100000)] + A[1:], PH12),
         'sdp': ([A[0] - Fr(2, 100000)] + A[1:], PH12 - Fr(2, 10 ** 6))}
for name, (a, p12) in FORMS.items():
    ph32 = -a[3]; c = a[3] + ph32; e = 2 * a[3] + ph32
    n = [4 * a[0] * c, 4 * (a[0] * a[4] + a[1] * c) - p12 ** 2, 4 * (a[1] * a[4] + a[2] * c) + 2 * p12 * e,
         4 * (a[2] * a[4] + a[3] * c) - (e * e - 4 * p12 * a[4])]
    am = [mp.mpf(x.numerator) / x.denominator for x in a]; nm = [mp.mpf(x.numerator) / x.denominator for x in n]
    rho = lambda x: (nm[0] + x * (nm[1] + x * (nm[2] + x * nm[3]))) / (4 * (am[0] + x * (am[1] + x * (am[2] + x * (am[3] + x * am[4])))))
    out = []
    for g, K in ((mp.pi, 4000), (mp.mpf('0.01'), 4000), (mp.mpf('1e-3'), 20000), (mp.mpf('1e-4'), 20000)):
        so = mp.fsum(rho((k * mp.pi / g) ** 2) for k in range(1, K + 1, 2))
        se = mp.fsum(rho((k * mp.pi / g) ** 2) for k in range(2, K + 1, 2))
        tail = nm[3] / (4 * am[4]) * (g / mp.pi) ** 2 / (2 * K)     # leading tail of either parity (rho ~ n3/(4 a4 x))
        out.append('g=%s: S_odd %s S_even %s (tail ~%s)' % (mp.nstr(g, 4), mp.nstr(so + tail, 6), mp.nstr(se + tail, 6), mp.nstr(tail, 2)))
    print(name, 'p(1) =', mp.nstr(sum(am), 10), 'n1..n3 =', [mp.nstr(x, 6) for x in nm[1:]])
    for s in out:
        print('   ', s)

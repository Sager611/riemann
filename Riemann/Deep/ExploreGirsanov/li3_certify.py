# li3_certify.py (explore56): Arb certificate that lambda_3 < 0 for an explicit positive symmetric law (so the
# change-of-measure positivity of BPY 2001 / arXiv 2411.08863 stops at lambda_2), plus a law matching zeta's first two
# Esscher cumulants (kt_1 = lambda_1, kt_2 = lambda_2 - 2 lambda_1) with lambda_3 < 0.
# M(w) = p cosh(a w) + (1-p) cosh(b w); F(s) = M(s-1/2)/M(1/2) = E[Y^s] for a positive Y, F(s) = F(1-s), F(0) = F(1) = 1.
# lambda_n = sum_m C(n,m) L_m/(m-1)!, L_m = (log M)^{(m)}(1/2).
from flint import arb, ctx
import mpmath as mp
ctx.prec = 200
def Ls(p, a, b):
    w = arb(1)/2
    M  = p*(a*w).cosh() + (1-p)*(b*w).cosh()
    M1 = p*a*(a*w).sinh() + (1-p)*b*(b*w).sinh()
    M2 = p*a**2*(a*w).cosh() + (1-p)*b**2*(b*w).cosh()
    M3 = p*a**3*(a*w).sinh() + (1-p)*b**3*(b*w).sinh()
    l1 = M1/M; l2 = M2/M - l1**2; l3 = M3/M - 3*M2*M1/M**2 + 2*l1**3
    return l1, l2, l3
def lams(p, a, b):
    l1, l2, l3 = Ls(p, a, b)
    return l1, 2*l1 + l2, 3*l1 + 3*l2 + l3/2
p, a, b = arb(999)/1000, arb(1)/2, arb(18)
L = lams(p, a, b)
print('law A: p=999/1000, a=1/2, b=18: lambda_1 =', L[0].str(12), ' lambda_2 =', L[1].str(12), ' lambda_3 =', L[2].str(12))
assert L[2] < 0 and L[0] > 0 and L[1] > 0
# zeros of M nearest the real axis (s = 1/2 + w): where the negative lambda_3 comes from
mp.mp.dps = 30
Mm = lambda w: mp.mpf('0.999')*mp.cosh(w/2) + mp.mpf('0.001')*mp.cosh(18*w)
zs = set()
for x0 in [0.3, 0.4, 0.5, 0.6]:
    for y0 in [0.05, 0.1, 0.2]:
        try:
            z = mp.findroot(Mm, mp.mpc(x0, y0)); zs.add((round(float(z.real), 6), round(float(abs(z.imag)), 6)))
        except Exception: pass
print('law A zeros w (Re, |Im|) near w = 1/2:', sorted(zs)[:6], ' i.e. s = 1/2 + w')
# law B: match zeta's kt_1, kt_2 at fixed b, solve for (p, a); then certify lambda_3 < 0
kt1, kt2 = mp.mpf('0.0230957089661210338'), mp.mpf('0.0461543172958046')
def Lm(p, a, b):
    M = lambda w: p*mp.cosh(a*w) + (1-p)*mp.cosh(b*w)
    return [mp.diff(lambda w: mp.log(M(w)), mp.mpf(1)/2, k) for k in (1, 2, 3)]
for bb in [6, 9, 12, 18]:
    try:
        sol = mp.findroot(lambda pp, aa: [Lm(pp, aa, bb)[0] - kt1, Lm(pp, aa, bb)[1] - kt2], (mp.mpf('0.9999'), mp.mpf('0.3')))
        pp, aa = sol[0], sol[1]
        l1, l2, l3 = Lm(pp, aa, bb)
        print('law B b=%d: p=%s a=%s kt1=%s kt2=%s lambda_3=%s' % (bb, mp.nstr(pp, 12), mp.nstr(aa, 8), mp.nstr(l1, 8), mp.nstr(l2, 8), mp.nstr(3*l1 + 3*l2 + l3/2, 8)))
        if 0 < pp < 1 and aa > 0:
            P, A = arb(mp.nstr(pp, 25)), arb(mp.nstr(aa, 25))
            LB = lams(P, A, arb(bb))
            print('   Arb: lambda_1 =', LB[0].str(10), ' lambda_2 =', LB[1].str(10), ' lambda_3 =', LB[2].str(10), ' certified <0:', bool(LB[2] < 0))
    except Exception as e:
        print('law B b=%d: no solution (%s)' % (bb, str(e)[:60]))

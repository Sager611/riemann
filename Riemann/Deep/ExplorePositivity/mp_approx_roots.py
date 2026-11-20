# mp_approx_roots.py — zeros of the MP approximants  Xi_N(t) = sum_{n<=N} (-1)^n c_2n f_2n(t/2)
# (polynomial of degree 2N in x = t/2), all roots isolated by Arb (python-flint acb_poly.roots).
# usage: mp_approx_roots.py <c2n file | toy:atom:u0 | toy:box:a | toy:gauss:s> <N list comma> <prec bits>
import sys, mpmath as mp
from fractions import Fraction as Fr
from flint import arb, acb, acb_poly, fmpq, ctx
from mp_core import f_coeffs
src, Ns, prec = sys.argv[1], [int(v) for v in sys.argv[2].split(',')], int(sys.argv[3])
ctx.prec = prec
NMAX = max(Ns)
mp.mp.dps = int(prec * 0.30) + 10
if src.startswith('toy:'):
    _, kind, par = src.split(':'); par = mp.mpf(par)
    if kind == 'atom':   # Phi = delta(u - u0) + delta(u + u0)  ->  cos(t u0) (LP)
        g = lambda u: 1
        cs = [mp.cosh(par)**(-1.5) * mp.tanh(par)**(2*n) for n in range(NMAX + 1)]
    else:
        if kind == 'box':   # Phi = 1 on [-a, a]  ->  2 sin(a t)/t (LP, Polya)
            w = lambda u: 1; lim = [0, par]
        else:               # gauss: Phi = exp(-u^2/(2 s^2)) -> no zeros (LP)
            w = lambda u: mp.exp(-u**2/(2*par**2)); lim = [0, 2*par, 5*par, 12*par]
        cs = [mp.quad(lambda u: w(u) * mp.cosh(u)**(-1.5) * mp.tanh(u)**(2*n), lim) for n in range(NMAX + 1)]
    cerr = mp.mpf(10)**(-(mp.mp.dps - 8))
else:
    cs = [mp.mpf(l.split()[1]) for l in open(src) if l.strip()][:NMAX + 1]
    cerr = mp.mpf(10)**(-int(sys.argv[4]) if len(sys.argv) > 4 else -44)   # cross-check level of mp_c2n.py (rel diff < 5e-45; 'hi' files: argv[4])
F = f_coeffs(2 * NMAX)
def ball(v, relerr):
    m = arb(mp.nstr(v, mp.mp.dps)); return arb(m.mid(), abs(m.mid()) * arb(mp.nstr(relerr, 5)) + m.rad())
for N in Ns:
    coef = [arb(0)] * (2 * N + 1)
    for n in range(N + 1):
        cb = ball(cs[n], cerr) * (1 if n % 2 == 0 else -1)
        for j, q in enumerate(F[2 * n]):
            if q != 0:
                coef[j] += cb * arb(fmpq(q.numerator, q.denominator))
    P = acb_poly([acb(c) for c in coef])
    try:
        rts = P.roots(maxprec=4 * prec)
    except Exception as e:
        print('N=%d  roots FAILED (%s)' % (N, str(e)[:80])); continue
    real = [r for r in rts if r.imag.contains(0)]
    nonreal = [r for r in rts if not r.imag.contains(0)]
    tr = sorted(float(abs(r.real.mid())) * 2 for r in real)   # t = 2x
    nr = sorted(nonreal, key=lambda r: float(abs(r.real.mid())))
    info = ''
    if nonreal:
        r0 = nr[0]; info = ' first nonreal t=%.4f%+.4fi' % (2*float(r0.real.mid()), 2*float(abs(r0.imag.mid())))
        rm = min(nonreal, key=lambda r: float(abs(r.imag.mid())))
        info += ' min|Im t|=%.3e at |Re t|=%.2f (resolution edge 4.51 sqrt N = %.1f)' % (2*float(abs(rm.imag.mid())), 2*float(abs(rm.real.mid())), 4.51*N**0.5)
    print('N=%3d deg=%3d real=%3d nonreal=%3d  smallest real t>0: %s  max real t: %s%s' % (
        N, 2*N, len(real), len(nonreal), ['%.5f' % v for v in tr if v > 0][:3:2], '%.2f' % tr[-1] if tr else '-', info), flush=True)

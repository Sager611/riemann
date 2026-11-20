# damped_mp.py — ceilaudit10 T2a: damped Meixner–Pollaczek (lam = 3/4) truncations of Xi, sum_{n<=N} mu_n (-1)^n c_2n f_2n(t/2):
# Fejer, Jackson, Jensen (N)_n/N^n and Gaussian dampings vs none; Arb root isolation (python-flint), c_2n from explore63 logs/c2n_zeta.txt.
import sys, mpmath as mp
from fractions import Fraction as Fr
from flint import arb, acb, acb_poly, ctx
sys.path.insert(0, 'Riemann/Deep/ExplorePositivity')
from mp_core import f_coeffs
mp.mp.dps = 60; ctx.prec = 800
c = {}
for line in open('Riemann/Deep/ExplorePositivity/logs/c2n_zeta.txt'):
    k, v = line.split()[:2]; c[int(k)] = mp.mpf(v)
NS = [int(v) for v in sys.argv[1].split(',')]
F = f_coeffs(2*max(NS))
def damp(kind, n, N):
    if kind == 'none': return mp.mpf(1)
    if kind == 'fejer': return 1 - mp.mpf(n)/(N+1)
    if kind == 'jackson':
        a = mp.pi/(N+1); return ((N-n+1)*mp.cos(n*a) + mp.sin(n*a)*mp.cot(a))/(N+1)
    if kind == 'jensen': return mp.rf(N-n+1, n)/mp.mpf(N)**n
    if kind == 'gauss': return mp.exp(-mp.mpf(n)**2/N)
for N in NS:
    for kind in ('none', 'fejer', 'jackson', 'jensen', 'gauss'):
        co = [mp.mpf(0)]*(2*N+1)
        for n in range(N+1):
            w = damp(kind, n, N) * (-1)**n * c[2*n]
            for j, v in enumerate(F[2*n]):
                if v != 0: co[j] += w * mp.mpf(v.numerator)/v.denominator / mp.mpf(2)**j
        rts = acb_poly([acb(arb(mp.nstr(v, 55), mp.nstr(abs(v)*mp.mpf(10)**-52, 3))) for v in co]).roots(maxprec=6000)
        nonreal = [r for r in rts if not r.imag.contains(0)]
        real = sorted(float(r.real.mid()) for r in rts if r.imag.contains(0) and float(r.real.mid()) > 0)
        mim = min((abs(float(r.imag.mid())) for r in nonreal), default=None)
        foot = min((abs(float(r.real.mid())) for r in nonreal), default=None)
        print('N=%d %-8s nonreal %d/%d  min|Im|=%s  foot|Re|=%s  real>0 <=30: %s' % (N, kind, len(nonreal), 2*N,
              None if mim is None else round(mim,3), None if foot is None else round(foot,2), [round(x,4) for x in real if x <= 30]), flush=True)

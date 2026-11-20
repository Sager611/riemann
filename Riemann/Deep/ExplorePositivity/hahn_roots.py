# hahn_roots.py — zeros of the continuous Hahn approximants  sum_{n<=N} (-1)^n d_2n g_2n(t/2)  (Arb root isolation)
import sys, mpmath as mp
from fractions import Fraction as Fr
from flint import arb, acb, acb_poly, fmpq, ctx
src, Ns, prec = sys.argv[1], [int(v) for v in sys.argv[2].split(',')], int(sys.argv[3])
ctx.prec = prec; mp.mp.dps = int(prec*0.3) + 10
NMAX = max(Ns)
ds = [mp.mpf(l.split()[1]) for l in open(src) if l.strip()][:NMAX + 1]
G = [[Fr(1)], [Fr(0), Fr(8, 3)]]
for n in range(1, 2*NMAX):
    c = [Fr(0)] * (n + 2)
    for k, v in enumerate(G[n]): c[k+1] += 8*v
    for k, v in enumerate(G[n-1]): c[k] -= (2*n + 1)*v
    G.append([v / (2*n + 3) for v in c])
cerr = mp.mpf(10)**(-38)
for N in Ns:
    coef = [arb(0)] * (2*N + 1)
    for n in range(N + 1):
        m = arb(mp.nstr(ds[n], 50)); cb = arb(m.mid(), abs(m.mid()) * arb(mp.nstr(cerr, 5))) * (1 if n % 2 == 0 else -1)
        for j, q in enumerate(G[2*n]):
            if q != 0: coef[j] += cb * arb(fmpq(q.numerator, q.denominator))
    try:
        rts = acb_poly([acb(c) for c in coef]).roots(maxprec=4*prec)
    except Exception as e:
        print('N=%d roots FAILED (%s)' % (N, str(e)[:60])); continue
    real = [r for r in rts if r.imag.contains(0)]; nonreal = [r for r in rts if not r.imag.contains(0)]
    tr = sorted(2*float(abs(r.real.mid())) for r in real)
    info = ''
    if nonreal:
        rm = min(nonreal, key=lambda r: float(abs(r.imag.mid())))
        info = ' min|Im t|=%.3e at |Re t|=%.2f (edge 2.96 N^(2/3) = %.1f)' % (2*float(abs(rm.imag.mid())), 2*float(abs(rm.real.mid())), 2.96*N**(2/3))
    print('N=%3d deg=%3d real=%3d nonreal=%3d first real t>0: %s max real t: %s%s' % (N, 2*N, len(real), len(nonreal), ['%.5f' % v for v in tr[::2][:4]], '%.1f' % tr[-1] if tr else '-', info), flush=True)

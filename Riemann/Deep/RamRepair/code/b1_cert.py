# b1_cert.py -- certify (B1): on B = {0.6 <= sigma <= 1.3, |t| <= 450}: |1/zeta| <= 14.57, and |zeta| <= 11.89 on B minus D(1, 0.1).
# 1/zeta is analytic on a neighbourhood of B (no zeros with sigma > 1/2 below height 3e12; the pole of zeta is a zero of 1/zeta), so
# max |1/zeta| over B is attained on the boundary (maximum modulus); zeta is analytic on B minus D(1,0.1), so max|zeta| is attained on
# the rectangle boundary or the circle.  |zeta(conj s)| = |zeta(s)|: edges with t >= 0 suffice.  Each edge is covered by balls.
import math, time
from flint import arb, acb, ctx
ctx.prec = 64
def edge(fixed_sigma, t_lo, t_hi, n):          # vertical edge sigma = fixed, t in [t_lo, t_hi]
    lo = None; hi = None
    for i in range(n):
        a = t_lo + (t_hi - t_lo)*i/n; b = t_lo + (t_hi - t_lo)*(i + 1)/n
        Z = acb(arb(fixed_sigma), arb((a + b)/2, (b - a)/2*(1 + 1e-12))).zeta()
        al, au = abs(Z).lower(), abs(Z).upper()
        lo = al if lo is None or al < lo else lo; hi = au if hi is None or au > hi else hi
    return lo, hi
def hedge(t, s_lo, s_hi, n):                   # horizontal edge
    lo = None; hi = None
    for i in range(n):
        a = s_lo + (s_hi - s_lo)*i/n; b = s_lo + (s_hi - s_lo)*(i + 1)/n
        Z = acb(arb((a + b)/2, (b - a)/2*(1 + 1e-12)), arb(t)).zeta()
        al, au = abs(Z).lower(), abs(Z).upper()
        lo = al if lo is None or al < lo else lo; hi = au if hi is None or au > hi else hi
    return lo, hi
t0 = time.time(); res = {}
res['left'] = edge(0.6, 0.0, 450.0, 225000); print('left sigma=0.6', res['left'], '%.0fs' % (time.time() - t0), flush=True)
res['right'] = edge(1.3, 0.0, 450.0, 45000); print('right sigma=1.3', res['right'], '%.0fs' % (time.time() - t0), flush=True)
res['top'] = hedge(450.0, 0.6, 1.3, 2000); print('top t=450', res['top'], '%.0fs' % (time.time() - t0), flush=True)
hi = None; n = 4000
for i in range(n):
    th = 2*arb.pi()*(arb(i) + arb('0.5'))/n
    cx, cy = 1 + arb('0.1')*th.cos(), arb('0.1')*th.sin(); r = float((arb('0.1')*arb.pi()/n).upper())*1.01
    Z = acb(arb(cx.mid(), r + float(cx.rad())), arb(cy.mid(), r + float(cy.rad()))).zeta()
    au = abs(Z).upper(); hi = au if hi is None or au > hi else hi
print('circle |s-1| = 0.1 max', hi, '%.0fs' % (time.time() - t0), flush=True)
mn = min(res['left'][0], res['right'][0], res['top'][0]); mx = max(res['left'][1], res['right'][1], res['top'][1], hi)
print('RESULT  min|zeta| on boundary >= %s  =>  max|1/zeta| <= %s  (<= 14.57: %s)' % (mn, (1/arb(mn)).str(8), 1/arb(mn) < arb('14.57')))
print('RESULT  max|zeta| on boundary of B minus disc <= %s  (<= 11.89: %s)' % (mx, arb(mx) < arb('11.89')))

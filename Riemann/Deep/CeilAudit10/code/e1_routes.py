# e1_routes.py — ceilaudit10: re-run e1rate's own exponent bookkeeping (Riemann/Deep/E1Rate/exponents.py, imported read-only) with
# route (a): no separation charge on the DI type II terms (phase kept; esep -> etr), type I by Lemma I unchanged;
# route (b): route (a) and type I pieces removed (eI -> -9).  Output: nu_K(theta) per route vs E1 (3t-1)/4 and 4t/7.
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import exponents as X
orig_eII, orig_eI = X.eII, X.eI
def eII_kept(th, pv, a, b, n):
    etr, esep, x0 = X.bounds(th, pv, a+b, n); return etr + X.gDI(a, b, n, pv)
ths = [float(v) for v in sys.argv[1].split(',')]
mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
for th in ths:
    if mode == 'a':
        X.eII = eII_kept; ra = X.nuK(th); X.eII = orig_eII
        w = (X.eII, X.eI); X.eII = eII_kept; worst = X.E(th, ra + 0.002); X.eII = orig_eII
        print('theta=%.4f  E1=(3t-1)/4=%.5f  route(a) nu=%.5f  (2/7)(3t-1)=%.5f  (2/3)(2t-1)=%.5f  binding block at nu+0.002: %s' % (
            th, (3*th-1)/4, ra, 2*(3*th-1)/7, 2*(2*th-1)/3, worst), flush=True)
        continue
    X.eII, X.eI = orig_eII, orig_eI; base = X.nuK(th)
    X.eII = eII_kept; ra = X.nuK(th)
    X.eI = lambda th_, pv, kq, a, n: -9.0; rb = X.nuK(th)
    X.eII, X.eI = orig_eII, orig_eI
    print('theta=%.4f  E1=(3t-1)/4=%.5f  nu_K(e1rate)=%.5f  route(a) nu=%.5f [(2/3)(2t-1)=%.5f]  route(b) nu=%.5f [4t/7=%.5f]' % (
        th, (3*th-1)/4, base, ra, 2*(2*th-1)/3, rb, 4*th/7), flush=True)

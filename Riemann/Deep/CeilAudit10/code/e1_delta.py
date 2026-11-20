# e1_delta.py — ceilaudit10: partial phase keeping. Type II (DI) charged X0^delta instead of e1rate's X0^(1/2) (delta = 1/2: e1rate; 0: route (a)).
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import exponents as X
th = float(sys.argv[1])
for d in [float(v) for v in sys.argv[2].split(',')]:
    def eII_d(th_, pv, a, b, n, d=d):
        etr, esep, x0 = X.bounds(th_, pv, a+b, n); return etr + d*max(0.0, x0) + X.gDI(a, b, n, pv)
    X.eII = eII_d
    print('theta=%.3f delta=%.3f nu_K=%.5f  (E1 %.5f, route(a) %.5f)' % (th, d, X.nuK(th), (3*th-1)/4, 2*(3*th-1)/7), flush=True)

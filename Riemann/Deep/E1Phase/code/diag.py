# e1phase diag.py — which block / piece / DI term binds route (a) (type II charged X0^delta) at theta; e1rate exponents.py imported read-only.
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import exponents as X
th = float(sys.argv[1]); d = float(sys.argv[2])
def eII_d(th_, pv, a, b, n):
    etr, esep, x0 = X.bounds(th_, pv, a+b, n); return etr + d*max(0.0, x0) + X.gDI(a, b, n, pv)
X.eII = eII_d
nu = X.nuK(th); print('theta=%.3f delta=%.3f nu*=%.5f E1=%.5f' % (th, d, nu, (3*th-1)/4))
nu2 = nu + 0.003
w, arg = X.E(th, nu2); print('at nu*+0.003: worst exponent %.5f at (P,Q,n)=%s' % (w, arg))
P, Q, n = arg
us = [nu2*j/40 for j in range(0, 21)]
def terms(a, b, n, v):
    m = max(n, a)
    return {'diag(VB)^1/2': -(a+b+n+v)/2+(v+b)/2, 'spec': -(a+b+n+v)/2 + m/4 + (v+b+m+max(v, 2*a))/4, 'NA2B2': -(a+b+n+v)/2 + m/4 + (n+2*a+2*b)/4}
for (pv, kq) in ((P, Q), (Q, P)):
    etr = X.bounds(th, pv, kq, n)[0]
    best = (9, None)
    for u in us:
        if kq <= u: worst = (etr, 'kq<=u trivial', None)
        else:
            worst = (-9, None, None)
            for i in range(17):
                a = min(2*u, kq)*i/16
                c = {'etr': etr, 'eI': X.eI(th, pv, kq, a, n), 'eII(a in)': X.eII(th, pv, a, kq-a, n), 'eII(b in)': X.eII(th, pv, kq-a, a, n)}
                k = min(c, key=c.get)
                if c[k] > worst[0]: worst = (c[k], 'piece2 a=%.4f %s' % (a, k), (a, kq-a))
            if kq-u > u:
                for i in range(17):
                    a = u+(kq-2*u)*i/16
                    c = {'etr': etr, 'eII(a in)': X.eII(th, pv, a, kq-a, n), 'eII(b in)': X.eII(th, pv, kq-a, a, n)}
                    k = min(c, key=c.get)
                    if c[k] > worst[0]: worst = (c[k], 'piece3 a=%.4f %s' % (a, k), (a, kq-a))
        if worst[0] < best[0]: best = (worst[0], (round(u, 4),) + worst[1:])
    print('modulus pv=%.4f kq=%.4f n=%.4f x0=%.4f etr=%.4f best-cut %.5f %s' % (pv, kq, n, X.bounds(th, pv, kq, n)[2], etr, best[0], best[1]))
    if best[1][2]:
        a, b = best[1][2]
        for (aa, bb) in ((a, b), (b, a)):
            print('   DI terms (A=%.4f,B=%.4f,N=%.4f,V=%.4f): %s' % (aa, bb, n, pv, {k: round(v, 4) for k, v in terms(aa, bb, n, pv).items()}))

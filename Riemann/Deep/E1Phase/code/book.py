# e1phase book.py — e1rate's exponent bookkeeping (exponents.py imported read-only) with ONE extra, known bound per block:
# the Bombieri–Iwaniec double large sieve with the Bessel phase kept as a bilinear phase (DERIVATION in REPORT §4):
#   relative error eDLS = etr + max((max(0,x0) - P - Q)/2, -n/2)   [arbitrary coefficients; needs the spacing count B_a << T^eps (R + R^2/(N X0))]
# mode 'dls': min(etr, eDLS) only;  mode 'all': min(etr, eDLS, e1rate's DI/Weil blocks Bdec (separated, delta = 1/2)).
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import exponents as X
def eDLS(th, P, Q, n):
    etr, esep, x0 = X.bounds(th, P, Q, n)
    return etr + max((max(0.0, x0) - P - Q)/2, -n/2)
def E2(th, nu, mode, G=16):
    us = [nu*j/40 for j in range(0, 21)]; worst = -9.0; arg = None
    for i in range(G+1):
        for j in range(G+1):
            P = nu*i/G; Q = nu*j/G; n0 = P+Q+1-2*th
            if n0 < 0: continue
            for n in (n0, 0.75*n0, 0.5*n0, 0.25*n0, 0.0):
                etr = X.bounds(th, P, Q, n)[0]
                b = min(etr, eDLS(th, P, Q, n))
                if mode == 'all' and b >= -0.01: b = min(b, X.Bdec(th, P, Q, n, us), X.Bdec(th, Q, P, n, us))
                if b > worst: worst, arg = b, (round(P, 4), round(Q, 4), round(n, 4))
    return worst, arg
def nu_of(th, mode, lo=0.0, hi=0.75):
    for _ in range(17):
        mid = (lo+hi)/2
        if E2(th, mid, mode)[0] < 0: lo = mid
        else: hi = mid
    return lo
mode = sys.argv[1]
for th in [float(v) for v in sys.argv[2].split(',')]:
    v = nu_of(th, mode)
    print('mode=%s theta=%.4f nu=%.5f  min(2t-1,t/2)=%.5f  E1=%.5f  e1rate nu_K=%.5f  ceilaudit route(a) (2/7)(3t-1)=%.5f  worst@nu+0.002: %s' % (
        mode, th, v, min(2*th-1, th/2), (3*th-1)/4, min(0.5, max((3*th-1)/4, 4*(2*th-1)/7)), 2*(3*th-1)/7, E2(th, v+0.002, mode)), flush=True)

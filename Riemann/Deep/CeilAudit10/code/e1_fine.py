# e1_fine.py — ceilaudit10 falsification of route (a): finer block grids (G = 24 for P, Q; 13 values of n in [0, n0]; 31 Vaughan cuts; K = 24 pieces)
# evaluate the worst block exponent E at nu just below / above (2/7)(3t-1) with type II phase kept (eII -> etr + gDI).
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, 'Riemann/Deep/E1Rate')
import exponents as X
def eII_kept(th, pv, a, b, n):
    etr, esep, x0 = X.bounds(th, pv, a+b, n); return etr + X.gDI(a, b, n, pv)
X.eII = eII_kept
def Efine(th, nu, G=24, NN=12, NU=30, K=24):
    us = [nu*j/(2*NU) for j in range(0, NU+1)]; worst = -9.0; arg = None
    for i in range(G+1):
        for j in range(G+1):
            P = nu*i/G; Q = nu*j/G; n0 = P+Q+1-2*th
            if n0 < 0: continue
            for k in range(NN+1):
                n = n0*k/NN
                etr = X.bounds(th, P, Q, n)[0]
                b = min(etr, X.Bdec(th, P, Q, n, us, K=K), X.Bdec(th, Q, P, n, us, K=K))
                if b > worst: worst, arg = b, (round(P,4), round(Q,4), round(n,4))
    return worst, arg
th = float(sys.argv[1])
for nu in [float(v) for v in sys.argv[2].split(',')]:
    print('theta=%.3f nu=%.5f fine-grid worst block exponent %s (negative = admissible); (2/7)(3t-1)=%.5f' % (th, nu, Efine(th, nu), 2*(3*th-1)/7), flush=True)

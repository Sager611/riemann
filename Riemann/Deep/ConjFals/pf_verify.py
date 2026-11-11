# High-precision re-evaluation (mpmath, 2 precisions) of det[K(x_i - y_j)] for the de Bruijn-Newman kernel K(u)=Phi(|u|)
import sys, json
from mpmath import mp, mpf, exp, pi, matrix, det, nstr, iv
def K(u, nmax=30):
    u = abs(u)
    return sum((2*pi**2*n**4*exp(9*u) - 3*pi*n**2*exp(5*u))*exp(-pi*n**2*exp(4*u)) for n in range(1, nmax+1))
def D(X, Y, dps):
    mp.dps = dps
    X = [mpf(x) for x in X]; Y = [mpf(y) for y in Y]
    M = matrix(len(X), len(Y))
    for i in range(len(X)):
        for j in range(len(Y)): M[i, j] = K(X[i] - Y[j])
    return det(M), M
X = json.loads(sys.argv[1]); Y = json.loads(sys.argv[2])
for dps in (40, 80):
    d, M = D(X, Y, dps)
    print('dps', dps, 'det =', nstr(d, 20), ' diag prod =', nstr(M[0,0]*M[1,1]*M[2,2]*(M[3,3] if len(X)>3 else 1), 6))

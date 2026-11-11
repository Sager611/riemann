# analysis 3: fit the log coefficient b in the smoothed count residual on large E, for each lambda
# R_s(E) = (Theta(E) - Phi(E))/pi - 1/4  = mean of N_+(E) - (E/pi)(log(2E/(pi lam^2)) - 1)   (per parity)
# For lam = sqrt2 and t = 2E this is the mean of N(t) - (t/2pi)(log(t/2pi) - 1) for the Dirac operator 2D.
import math, sys
D = sys.argv[1] if len(sys.argv) > 1 else '.'
def lstsq(X, y):
    # normal equations, small dense
    n = len(X[0]); A = [[sum(r[i]*r[j] for r in X) for j in range(n)] for i in range(n)]
    b = [sum(r[i]*yy for r, yy in zip(X, y)) for i in range(n)]
    # gaussian elimination
    M = [A[i] + [b[i]] for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda k: abs(M[k][i])); M[i], M[p] = M[p], M[i]
        for k in range(n):
            if k != i:
                f = M[k][i]/M[i][i]
                M[k] = [a - f*c for a, c in zip(M[k], M[i])]
    coef = [M[i][n]/M[i][i] for i in range(n)]
    res = [yy - sum(c*x for c, x in zip(coef, r)) for r, yy in zip(X, y)]
    dof = len(y) - n; s2 = sum(r*r for r in res)/dof if dof > 0 else float('nan')
    # covariance diag via inverse of A
    inv = [[float(i == j) for j in range(n)] for i in range(n)]
    Ac = [row[:] for row in A]
    for i in range(n):
        p = max(range(i, n), key=lambda k: abs(Ac[k][i])); Ac[i], Ac[p] = Ac[p], Ac[i]; inv[i], inv[p] = inv[p], inv[i]
        piv = Ac[i][i]; Ac[i] = [a/piv for a in Ac[i]]; inv[i] = [a/piv for a in inv[i]]
        for k in range(n):
            if k != i:
                f = Ac[k][i]; Ac[k] = [a - f*c for a, c in zip(Ac[k], Ac[i])]; inv[k] = [a - f*c for a, c in zip(inv[k], inv[i])]
    se = [math.sqrt(s2*inv[i][i]) for i in range(n)]
    return coef, se, math.sqrt(sum(r*r for r in res)/len(res))
for lab, fn, lam in (('sqrt2', 'lg1.txt', math.sqrt(2)), ('1', 'lg_1.txt', 1.0), ('0.5', 'lg_0.5.txt', 0.5), ('2', 'lg_2.txt', 2.0)):
    kap = 2*math.pi*lam**2
    rows = [tuple(map(float, l.split()[:2])) for l in open(f'{D}/{fn}')]
    pts = [(E, (T - (E*math.log(4*E/kap) - E + math.pi/4))/math.pi - 0.25) for E, T in rows if E >= 999]
    L = [math.log(E) for E, _ in pts]; y = [r for _, r in pts]
    # free b with (alpha log E + beta)/E corrections
    X = [[1, l, l/math.exp(l), 1/math.exp(l)] for l in L]
    (c, b, al, be), se, rms = lstsq(X, y)
    # B26 model: b = 0
    X0 = [[1, l/math.exp(l), 1/math.exp(l)] for l in L]
    (c0, al0, be0), se0, rms0 = lstsq(X0, y)
    # CM22 model: b = -1/(2pi) fixed  (in t = 2E the log coefficient is the same)
    yc = [yy + l/(2*math.pi) for yy, l in zip(y, L)]
    (c1, al1, be1), se1, rms1 = lstsq(X0, yc)
    print(f'lambda={lab:5s} E in [1e3,1e5] ({len(pts)} pts): free fit b = {b: .2e} +- {se[1]:.1e} (CM22: -0.159155), c = {c: .6f};'
          f' b=0 fit: c = {c0: .6f} rms {rms0:.1e}; CM22-forced fit rms {rms1:.1e}')
    print(f'      R_s(1e3) = {pts[0][1]: .6f}, R_s(1e4) = {[r for E, r in pts if abs(E-1e4) < 1][0]: .6f}, R_s(1e5) = {pts[-1][1]: .6f}   (B26 limit -1/4)')

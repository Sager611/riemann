# Local PF_r test of the de Bruijn-Newman kernel K(u) = Phi(|u|) (Michalowski 2602.20313): the coalescent limit of every
# r x r minor det[K(x_i - y_j)] has the sign of H_r(u) = det[(-1)^j K^{(i+j)}(u)]_{i,j<r} (u = x - y). Scan H_3, H_4 (and H_5 as control).
import sys
from mpmath import mp, mpf, exp, pi, matrix, det, nstr, taylor, factorial
mp.dps = 60
def K(u):
    return sum((2*pi**2*n**4*exp(9*u) - 3*pi*n**2*exp(5*u))*exp(-pi*n**2*exp(4*u)) for n in range(1, 25))
def H(u, r, der):
    M = matrix(r, r)
    for i in range(r):
        for j in range(r): M[i, j] = (-1)**j*der[i+j]
    return det(M)
a, b, steps = mpf(sys.argv[1]), mpf(sys.argv[2]), int(sys.argv[3])
prev = {}
for k in range(steps+1):
    u = a + (b-a)*k/steps
    c = taylor(K, u, 8)                     # Taylor coefficients c_m = K^{(m)}(u)/m!
    der = [c[m]*factorial(m) for m in range(9)]
    h = {r: H(u, r, der) for r in (2, 3, 4, 5)}
    sc = {r: h[r]/der[0]**r for r in h}      # scale-free
    for r in h:
        if r in prev and prev[r]*sc[r] < 0: print('SIGN CHANGE r=%d between u=%s and %s' % (r, nstr(u-(b-a)/steps, 8), nstr(u, 8)))
        prev[r] = sc[r]
    if k % max(1, steps//20) == 0:
        print('u=%s  H2/K^2=%s H3/K^3=%s H4/K^4=%s H5/K^5=%s' % (nstr(u,5), nstr(sc[2],6), nstr(sc[3],6), nstr(sc[4],6), nstr(sc[5],6)))

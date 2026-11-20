"""Finite free cumulants of Jensen polynomials (x-family, GORZ), C1 and C2.
Reversed J*^{d,n}(x) = sum_j C(d,j) gamma(n+j) x^{d-j} = P(D) x^d, P(y) = sum_j gamma(n+j)/gamma(n) y^j/j! (mod y^{d+1}).
l_k(n) = [y^k] log P (additive under (+)_d, independent of d for k <= d).  Normalised rho_k = l_k/(-2 l_2)^{k/2}:
after an affine change of x, J* = exp(-D^2/2 + sum_{k>=3} rho_k D^k) x^d (Hermite He_d has rho_k = 0).
Original J^{d,n}(X) = sum_i C(d,i) gamma(n+d-i) X^{d-i}: symbol Q(y) = sum_i gamma(n+d-i)/gamma(n+d) y^i/i! (depends on d)."""
import json, sys, mpmath as mp
mp.mp.dps = 120
G = json.load(open('logs/gamma.json'))
seq = {k: [mp.mpf(x) for x in v] for k, v in G.items()}
def logcoef(p, K):           # p[0] = 1; returns l[1..K] of log(sum p_j y^j)
    l = [mp.mpf(0)]*(K+1)
    for k in range(1, K+1):
        s = k*p[k] - mp.fsum(i*l[i]*p[k-i] for i in range(1, k))
        l[k] = s/k
    return l
def rho_rev(g, n, K):
    p = [g[n+j]/g[n]/mp.factorial(j) for j in range(K+1)]
    l = logcoef(p, K); lam = mp.sqrt(-2*l[2])
    return [None, None, mp.mpf(-0.5)] + [l[k]/lam**k for k in range(3, K+1)], l
def rho_orig(g, n, d, K):
    p = [g[n+d-i]/g[n+d]/mp.factorial(i) for i in range(K+1)]
    l = logcoef(p, K); lam = mp.sqrt(-2*l[2])
    return [None, None, mp.mpf(-0.5)] + [l[k]/lam**k for k in range(3, K+1)], l
if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'C1':   # d = 4 heat/affine invariant I = l_3^4/l_4^3 along the de Bruijn-Newman time t
        print('n   I_rev(t=0)  I_rev(.05)  I_rev(.1)  I_rev(.2) | I_orig(t=0) I_orig(.05) I_orig(.1) I_orig(.2)   [d=4]')
        for n in (0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000):
            row = []
            for mode in ('rev', 'orig'):
                for t in ('0', '0.05', '0.1', '0.2'):
                    g = seq['xi_t' + t]
                    _, l = rho_rev(g, n, 4) if mode == 'rev' else rho_orig(g, n, 4, 4)
                    row.append(mp.nstr(l[3]**4/l[4]**3, 8))
            print(n, ' '.join(row))
    if what == 'C2':
        for name in ('xi_t0', 'dh_t0'):
            g = seq[name]
            print(name, ' n  rho3  rho3*sqrt(n)  rho4*n  rho5*n^1.5  l2')
            for n in (0, 1, 2, 5, 10, 20, 50, 100, 200, 400, 800, 1150):
                r, l = rho_rev(g, n, 6)
                print(' %5d %s %s %s %s %s' % (n, mp.nstr(r[3], 8), mp.nstr(r[3]*mp.sqrt(max(n,1)), 8), mp.nstr(r[4]*max(n,1), 8), mp.nstr(r[5]*mp.mpf(max(n,1))**1.5, 8), mp.nstr(l[2], 6)))

"""Taylor data of Riemann's xi and of the Davenport-Heilbronn (DH) function, x-variable (GORZ):
  8 xi(1/2+sqrt x) = sum_k gamma(k) x^k/k!,  gamma(k) = 16 k! b_k/(2k)!,  b_k = int_0^inf Phi(u) u^{2k} du,
  Phi(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})   (xi(1/2+z) = 2 int_0^inf Phi cosh(zu) du).
DH: f = sum a_n n^-s, a = (1, k, -k, -1, 0) mod 5, k = (sqrt(10-2sqrt5)-2)/(sqrt5-1); completed
  Lam_f(1/2+z) = int psi(u) e^{uz} du, psi(u) = 2 e^{3u/2} sum_n n a_n exp(-pi n^2 e^{2u}/5) (even in u),
  gamma_DH(k) = 2 k! B_k/(2k)!, B_k = int_0^inf psi u^{2k} du (overall constants are irrelevant for real-rootedness).
Heat flow (de Bruijn-Newman, Xi_t = int Phi e^{t u^2} cos(zu)): b_k(t) = int Phi e^{t u^2} u^{2k} du.
Gauss-Legendre on [0,U] in panels; accuracy = difference between two node counts (printed).
usage: python3 gamma.py K dps"""
import sys, json, mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
K = int(sys.argv[1]); mp.mp.dps = int(sys.argv[2])
def Phi(u):
    s = mp.mpf(0)
    for n in range(1, 40):
        e = mp.exp(-mp.pi*n*n*mp.exp(4*u))
        t = (2*mp.pi**2*n**4*mp.exp(9*u) - 3*mp.pi*n*n*mp.exp(5*u))*e
        s += t
        if n > 2 and abs(t) < mp.mpf(10)**(-mp.mp.dps-30)*abs(s): break
    return s
kap = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
aDH = [0, 1, kap, -kap, -1]
def psi(u):
    s = mp.mpf(0)
    for n in range(1, 400):
        a = aDH[n % 5]
        if a == 0: continue
        t = n*a*mp.exp(-mp.pi*n*n*mp.exp(2*u)/5)
        s += t
        if n > 5 and abs(t) < mp.mpf(10)**(-mp.mp.dps-30)*abs(s): break
    return 2*mp.exp(1.5*u)*s
def nodes(U, panels, deg):
    gl = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
    out = []
    h = mp.mpf(U)/panels
    for p in range(panels):
        a = p*h
        for x, w in gl:
            out.append((a + h*(x+1)/2, w*h/2))
    return out
def moments(f, U, panels, deg, K, tlist=(0,)):
    nd = nodes(U, panels, deg)
    vals = [(u, w*f(u)) for u, w in nd]
    res = {}
    for t in tlist:
        cur = [wf*mp.exp(t*u*u) for u, wf in vals]
        u2 = [u*u for u, wf in vals]
        m = []
        for k in range(K+1):
            m.append(mp.fsum(cur))
            cur = [c*v for c, v in zip(cur, u2)]
        res[t] = m
    return res
if __name__ == '__main__':
    print('DH evenness psi(0.3)-psi(-0.3) =', mp.nstr(psi(mp.mpf('0.3'))-psi(mp.mpf('-0.3')), 5), ' psi(0.3) =', mp.nstr(psi(mp.mpf('0.3')), 10))
    tl = (0, mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('0.2'))
    out = {}
    for name, f, U, P in (('xi', Phi, 2.5, 20), ('dh', psi, 5.5, 44)):
        A = moments(f, U, P, 6, K, tl if name == 'xi' else (0,))
        B = moments(f, U, P, 7, K, (0,))
        err = max(abs(A[0][k]/B[0][k]-1) for k in range(K+1))
        print(name, 'max rel diff deg6 vs deg7 over k<=K:', mp.nstr(err, 3))
        for t, m in A.items():
            c = 16 if name == 'xi' else 2
            g = [c*mp.factorial(k)*m[k]/mp.factorial(2*k) for k in range(K+1)]
            out['%s_t%s' % (name, mp.nstr(t, 3))] = [mp.nstr(x, mp.mp.dps-10) for x in g]
        neg = [k for k in range(K+1) if A[0][k] < 0]
        print(name, 'negative moments at k =', neg[:20])
    json.dump(out, open('logs/gamma.json', 'w'))
    print('gamma_xi(0..3) =', out['xi_t0'][:4] if False else [x[:20] for x in out['xi_t0'][:4]])
    print('gamma_dh(0..3) =', [x[:20] for x in out['dh_t0'][:4]])

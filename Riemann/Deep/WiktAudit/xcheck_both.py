# Independent x-space Galerkin matrix of the Weil form (eq. (2.1) of rStation #60 / our PROOF.md form) in BOTH parity
# bases used by Wiktorowicz (even: cos(q pi x/a); odd: sin((q+1/2) pi x/a), extended to (-a,a)), at a = 23029/10000.
# Archimedean term by Weil's regularisation  base*g(0) + int_0^L 2varpi(u) (g(0) - Re g(u)) du  (no Green identity,
# no response family, no displacement formulas); autocorrelations in closed form; primes from the autocorrelation;
# poles from closed-form moments. Compared with W-W source.py rows: H = (a/2) E A E, E_qq = (-1)^q.
import sys, time
from flint import arb, arb_mat, ctx
KOD = sys.argv[1]; IDX = [int(s) for s in sys.argv[2].split(',')]; NB = len(IDX); PAN = int(sys.argv[3]) if len(sys.argv) > 3 else 40; ctx.prec = 256
a = arb(23029)/10000; L = 2*a; pi = arb.pi()
def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return arb(p).log() if m == 1 else None
atoms = [(arb(n).log(), vm(n)/arb(n).sqrt()) for n in range(2, 101) if vm(n) is not None]
assert len(atoms) == 35
def basis(par, q):
    k = arb(q) + arb(par)/2; kap = pi*k/a
    nq = (1/(2*a)).sqrt() if (par == 0 and q == 0) else (1/a).sqrt()
    return kap, nq
def Icos(al, be, x0, x1):          # int_{x0}^{x1} cos(al x + be) dx
    if al == 0: return (x1-x0)*be.cos()
    return ((al*x1+be).sin() - (al*x0+be).sin())/al
def G1(par, q, r, u):               # int_{-a}^{a-u} f_q(x+u) f_r(x) dx, 0 <= u <= 2a
    kq, nq = basis(par, q); kr, nr = basis(par, r)
    s = 1 if par == 0 else -1
    x0, x1 = -a, a-u
    t1 = Icos(kq-kr, kq*u, x0, x1) if q != r else (x1-x0)*(kq*u).cos()
    t2 = Icos(kq+kr, kq*u, x0, x1)
    return nq*nr*(t1 + s*t2)/2
def build(par, panels=None, nodes=60):
    panels = panels or PAN
    gl = [arb.legendre_p_root(nodes, i, weight=True) for i in range(nodes)]
    A = arb_mat(NB, NB); KL = arb(0)
    acc = [[arb(0)]*NB for _ in range(NB)]
    for p in range(panels):
        a0 = L*p/panels; a1 = L*(p+1)/panels; h = (a1-a0)/2; m0 = (a1+a0)/2
        for (x, w) in gl:
            u = m0 + h*x; ww = w*h
            den = -((-2*u).expm1()); two = 2*(-u/2).exp()/den
            KL += ww*(2*(-2*u).exp()/den - two)
            for i in range(NB):
                for j in range(i, NB):
                    q, r = IDX[i], IDX[j]
                    sym = (G1(par, q, r, u) + G1(par, r, q, u))/2
                    acc[i][j] += ww*two*((1 if q == r else 0) - sym)
    base = -arb.const_euler() - pi.log() - (-((-2*L).expm1())).log() + KL
    for i in range(NB):
        for j in range(i, NB):
            q, r = IDX[i], IDX[j]
            v = acc[i][j] + (base if q == r else 0)
            for d, wn in atoms:
                v -= wn*(G1(par, q, r, d) + G1(par, r, q, d))
            kq, nq = basis(par, q); kr, nr = basis(par, r)
            sq = 1 if q % 2 == 0 else -1; sr = 1 if r % 2 == 0 else -1
            if par == 0:   # 2 <f,cosh(x/2)>^2 ; <f_q,cosh> = n_q sinh(a/2) (-1)^q/(1/4+kap^2)
                v += 2*(nq*sq*(a/2).sinh()/(arb(1)/4+kq*kq))*(nr*sr*(a/2).sinh()/(arb(1)/4+kr*kr))
            else:          # -2 <f,sinh(x/2)>^2 ; <f_q,sinh> = n_q cosh(a/2) (-1)^q/(1/4+kap^2)
                v -= 2*(nq*sq*(a/2).cosh()/(arb(1)/4+kq*kq))*(nr*sr*(a/2).cosh()/(arb(1)/4+kr*kr))
            A[i, j] = v; A[j, i] = v
    return A
sys.path.insert(0, KOD)
from source import Source, matrix
for par in (0, 1):
    t0 = time.time()
    A = build(par)
    src = Source('2.3029', par, 256)
    H = matrix([src.row(q) for q in IDX])
    worst = 0.0; wq = None
    for i in range(NB):
        for j in range(NB):
            q, r = IDX[i], IDX[j]
            mine = a/2*(1 if (q+r) % 2 == 0 else -1)*A[i, j]
            rel = float((abs(mine-H[i, j])/abs(H[i, j])).mid())
            if rel > worst: worst = rel; wq = (q, r, H[i, j].str(15, radius=False), mine.str(15, radius=False))
    print('indices %s parity %d, %dx%d block: max relative difference %.3e at %s  (%.1f s)' % (IDX, par, NB, NB, worst, wq, time.time()-t0), flush=True)

# test_hp.py: 50-digit check of the closed-form entries of weil_matrix.py (digamma / Hurwitz forms, c^{-b_k} sums) against
# mpmath quadrature of the x-space formula A = -(log pi + gamma)<f1,f2> + int_0^inf [2e^{-x}<f1,f2> - 2e^{x/2}C12(x)]/(e^x - e^{-x}) dx
# (C12 by quadrature of int b_i(y) b_j(y - x) dy), plus pole and prime parts by quadrature. Usage: test_hp.py
import sys, mpmath as mp
from flint import arb, ctx, fmpq
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weil_matrix import WeilMat
mp.mp.dps = 50; ctx.prec = 256
def f2m(x): return mp.mpf(x.mid().str(70, radius=False))
for cq, pairs in [(fmpq(5), [('even', 0, 0), ('even', 0, 3), ('even', 2, 2), ('even', 1, 4), ('odd', 1, 1), ('odd', 2, 3)]),
                  (fmpq(25, 8), [('even', 0, 2), ('even', 3, 3), ('odd', 1, 3), ('odd', 4, 4)])]:
    W = WeilMat(cq); L = mp.log(mp.mpf(int(cq.p))/int(cq.q))
    blocks = {s: W.block(5, s) for s in ('even', 'odd')}
    for sector, i, j in pairs:
        T, idx = blocks[sector]; a, b = idx.index(i), idx.index(j)
        Ni = 1/mp.sqrt(L) if i == 0 else mp.sqrt(2/L); Nj = 1/mp.sqrt(L) if j == 0 else mp.sqrt(2/L)
        wi, wj = 2*mp.pi*i/L, 2*mp.pi*j/L
        if sector == 'even': bi = lambda y: Ni*mp.cos(wi*y); bj = lambda y: Nj*mp.cos(wj*y); v = lambda y: mp.cosh(y/2); ps = 2
        else: bi = lambda y: Ni*mp.sin(wi*y); bj = lambda y: Nj*mp.sin(wj*y); v = lambda y: mp.sinh(y/2); ps = -2
        C = lambda x: mp.quad(lambda y: bi(y)*bj(y - x), [x - L/2, 0, L/2] if x < L/2 else [x - L/2, L/2])
        dij = 1 if i == j else 0
        arch = -(mp.log(mp.pi) + mp.euler)*dij
        arch += mp.quad(lambda x: (2*mp.exp(-x)*dij - 2*mp.exp(x/2)*C(x))/(mp.exp(x) - mp.exp(-x)), [0, L/4, L/2, 3*L/4, L])
        arch += mp.quad(lambda x: 2*mp.exp(-x)*dij/(mp.exp(x) - mp.exp(-x)), [L, mp.inf])
        pole = ps*mp.quad(lambda y: bi(y)*v(y), [-L/2, 0, L/2])*mp.quad(lambda y: bj(y)*v(y), [-L/2, 0, L/2])
        pr = 0
        for n, ln, wt in W.atoms: pr += mp.log(n if all(n % p for p in range(2, int(n**0.5) + 1)) else [p for p in range(2, n) if n % p == 0][0])/mp.sqrt(n)*C(mp.log(n))
        direct = pole + arch - 2*pr
        cf = f2m(T[a, b])
        print(str(cq), sector, i, j, mp.nstr(direct, 30), 'diff', mp.nstr(direct - cf, 3), 'radius', T[a, b].rad().str(3, radius=False), flush=True)

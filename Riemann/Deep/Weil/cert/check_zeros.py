# pyright: reportMissingImports=false
# python-flint lives in a local venv (see README.md), not the system interpreter.

import mpmath as mp
from flint import arb_mat, ctx
import weil
ctx.prec = 200; mp.mp.dps = 40
for Ls in ['1', '1.5']:
    L = weil.parseL(Ls); a = L/2; N = 30
    M = weil.weil_matrix(L, N)
    u = arb_mat(N, 1)
    for (x, w) in weil.gauss(N + 20):
        v = weil.leg_on(N, x); f = (1 - x*x)**6
        for j in range(N): u[j, 0] += f*v[j]*w*a.sqrt()
    qm = (u.transpose()*M*u)[0, 0]
    am = mp.mpf(Ls)/2
    fh = lambda r: 2*mp.quad(lambda y: (1-(y/am)**2)**6*mp.cos(r*y), [0, am])
    zs = 0
    for k in range(1, 151):
        zs += 2*fh(mp.im(mp.zetazero(k)))**2
    print(Ls, 'matrix', qm.mid().str(15, radius=False), 'zeros(150)', mp.nstr(zs, 15))

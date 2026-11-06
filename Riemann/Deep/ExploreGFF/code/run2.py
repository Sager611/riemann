# run2.py (explore17, T2): Beurling-Deny (Markov) structure of W - pole; Perron-Frobenius check; sign changes of the true ground state.
import os, sys, time, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'WeilSimple'))
from gw import *
t0 = time.time()
def evalf(B, v, x):   # evaluate sum v_j b_j(x)
    if B.sector == 'even': return sum(v[a]*B.Nj[a]*np.cos(B.w[a]*x) for a in range(len(v)))
    return sum(v[a]*math.sqrt(2/B.L)*np.sin(B.w[a]*x) for a in range(len(v)))
def sign_changes(y, tol):
    s = np.sign(np.where(abs(y) < tol, 0, y)); s = s[s != 0]
    return int(np.sum(s[1:] != s[:-1]))
print('# float64 part (N = 40): even sector; columns: lam0(Q) lam1(Q) | sign changes of Q ground state on (0, L/2) | lam0(A=Q-pole), sign changes of A ground state, min/max ratio')
N = 40
for c in (1.5, 2.0, 2.5, 3.0, 4.0):
    for sec in ('even', 'odd'):
        P = weil_parts(c, N, sec, 600); B = P['B']
        Q = P['LK'] + P['V'] + P['POLE'] + P['PR']; A = Q - P['POLE']
        e, V = np.linalg.eigh(Q); ea, Va = np.linalg.eigh(A)
        x = np.linspace(1e-6, B.L/2 - 1e-6, 4001)
        g = evalf(B, V[:, 0], x); ga = evalf(B, Va[:, 0], x)
        if ga[len(ga)//2] < 0: ga = -ga
        print('  c=%-4g %s lam0(Q)=%.3e lam1(Q)=%.3e  sc(Q-gs)=%d | lam0(A)=%+.4f sc(A-gs)=%d min/max(A-gs)=%+.3f' % (
            c, sec, e[0], e[1], sign_changes(g, 1e-9*abs(g).max()), ea[0], sign_changes(ga, 1e-9*abs(ga).max()), ga.min()/ga.max()), flush=True)
print('time', round(time.time() - t0, 1), flush=True)
# high precision part (flint WeilMat closed forms): even ground state of Q by inverse iteration
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat
ctx.prec = 600
for cq, NN in ((5, 40), (10, 50), (20, 60)):
    Wm = WeilMat(fmpq(cq))
    for sec in ('even', 'odd'):
        T, idx = Wm.block(NN, sec); n = len(idx)
        v = arb_mat(n, 1); 
        for a in range(n): v[a, 0] = arb(1)
        for it in range(25):
            v = T.solve(v)
            nrm = sum((v[a, 0]**2 for a in range(n)), arb(0)).sqrt()
            v = v*(1/nrm)
        Tv = T*v; lam = sum((v[a, 0]*Tv[a, 0] for a in range(n)), arb(0))
        vf = np.array([float(v[a, 0].mid()) for a in range(n)])
        B = Basis(float(cq), NN, sec)
        x = np.linspace(1e-6, B.L/2 - 1e-6, 8001); g = evalf(B, vf, x)
        # residual check
        res = max(abs(float((Tv[a, 0] - lam*v[a, 0]).mid())) for a in range(n))
        print('  c=%d N=%d %s Rayleigh=%s resid=%.1e sign changes on (0,L/2)=%d' % (cq, NN, sec, lam.str(5, radius=False), res, sign_changes(g, 1e-6*abs(g).max())), flush=True)
print('time', round(time.time() - t0, 1))

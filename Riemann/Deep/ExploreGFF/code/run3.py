# run3.py (explore17, T2b): Sturm/Gantmacher-Krein oscillation test: sign changes of the k-th eigenfunction of the Weil form Q_W
# (all eigenvalues of both sectors merged and sorted), float64 for small c, mpmath (from WeilMat Arb closed forms) for larger c.
import os, sys, time, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'WeilSimple'))
from gw import *
t0 = time.time()
def evalf(L, sector, v, x):
    idx = np.arange(0, len(v)) if sector == 'even' else np.arange(1, len(v) + 1)
    w = 2*math.pi*idx/L
    if sector == 'even':
        Nj = np.where(idx == 0, 1/math.sqrt(L), math.sqrt(2/L)); return sum(v[a]*Nj[a]*np.cos(w[a]*x) for a in range(len(v)))
    return sum(v[a]*math.sqrt(2/L)*np.sin(w[a]*x) for a in range(len(v)))
def sc(y):
    tol = 1e-7*abs(y).max(); s = np.sign(np.where(abs(y) < tol, 0, y)); s = s[s != 0]; return int(np.sum(s[1:] != s[:-1]))
def report(c, L, eig):   # eig: list of (lam, sector, vec)
    eig = sorted(eig, key=lambda e: e[0]); x = np.linspace(-L/2 + 1e-9, L/2 - 1e-9, 20001)
    out = []
    for k, (lam, sec, v) in enumerate(eig[:8]):
        y = evalf(L, sec, v, x); out.append('%d:%s%.1e/%d' % (k, sec[0], lam, sc(y)))
    print('  c=%-5g ' % c + '  '.join(out), flush=True)
print('# k:sector lambda_k / sign changes of the k-th eigenfunction on the window (Sturm: k)')
for N in (40, 60):
    print(' N =', N)
    for c in (1.5, 2.0, 3.0, 4.0):
        eig = []
        for sec in ('even', 'odd'):
            P = weil_parts(c, N, sec, 700); Q = P['LK'] + P['V'] + P['POLE'] + P['PR']
            e, V = np.linalg.eigh(Q)
            for k in range(6): eig.append((e[k], sec, V[:, k]))
        report(c, math.log(c), eig)
print('time', round(time.time() - t0, 1), flush=True)
import mpmath as mp
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat
for cq, NN, dps in ((5, 40, 60), (10, 40, 90), (20, 40, 160)):
    ctx.prec = int(dps*3.4) + 64; mp.mp.dps = dps
    Wm = WeilMat(fmpq(cq)); eig = []
    for sec in ('even', 'odd'):
        T, idx = Wm.block(NN, sec); n = len(idx)
        M = mp.matrix(n, n)
        for a in range(n):
            for b in range(n): M[a, b] = mp.mpf(T[a, b].mid().str(dps + 10, radius=False))
        E, Vm = mp.eigsy(M)
        order = sorted(range(n), key=lambda k: E[k])
        for k in order[:6]: eig.append((float(E[k]), sec, np.array([float(Vm[a, k]) for a in range(n)])))
    report(cq, math.log(cq), eig)
    print('   time', round(time.time() - t0, 1), flush=True)

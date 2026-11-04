# solve_system.py (was t4): independent solution of the coupled interpolation system (float64), comparison with the table,
# and recomputation of the finite matrix norms quoted in the paper.
import sys, os, json, re, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fastnp import *
import os as _os
_ROOT = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))
SRC = _os.environ.get('CIRC2D_SRC', _os.path.join(_ROOT, '.scratch', 'circ2d', 'src'))
SCR = _os.environ.get('CIRC2D_SCRATCH', _os.path.join(_ROOT, '.scratch', 'circ2d'))
tab = read_table(SRC + '/build/data/coefficients.tex')
t0 = time.time()
Q = Quad(gl=48, sub=4)
NM = 240
ns = nodes(NM)
N = len(ns)
QD = np.array([node_QD(n) for n in ns]); Qn = QD[:, 0]; Dn = QD[:, 1]
Dens = np.zeros((2 * N, len(Q.t)), complex)
for k, n in enumerate(ns):
    cd, dd = Q.basis(n); Dens[k] = cd; Dens[N + k] = dd
Ms = np.array(ns, float)
Eval = np.exp(1j * PI * np.outer(Ms, Q.z)) * (Q.w * Q.lam)
Eder = Eval * (1j * PI * Q.z)
Kv = (Eval @ Dens.T).real; Kd = (Eder @ Dens.T).real
S = np.vstack([-Kv / Qn[:, None], (Dn[:, None] * Kv - Kd) / Qn[:, None]])   # rows c..., d...; cols c..., d...
def gvec(i):
    (c0, d0), C = FIXED[3 - i]
    cd, dd = Q.basis(0)
    f = c0 * cd + d0 * dd
    K0 = Q.K(f, C, Ms, 0); K1 = Q.K(f, C, Ms, 1)
    g = np.concatenate([-K0 / Qn, (Dn * K0 - K1) / Qn])
    if i == 1:
        g[N + ns.index(1)] -= 1.0
    return g
g1, g2 = gvec(1), gvec(2)
I = np.eye(2 * N)
Big = np.block([[I, -S], [-S, I]])
x = np.linalg.solve(Big, np.concatenate([g1, g2]))
x1, x2 = x[:2 * N], x[2 * N:]
print('solve done', round(time.time() - t0, 1), 's; cond(I-S)=%.3g cond(I+S)=%.3g' % (np.linalg.cond(I - S), np.linalg.cond(I + S)))
# compare with table
tot = [0.0, 0.0]; mx = [0.0, 0.0]; tail = [0.0, 0.0]
for k, n in enumerate(ns):
    for i, xx in ((0, x1), (1, x2)):
        c, d = xx[k], xx[N + k]
        if n <= 100:
            ct, dt = tab[n][2 * i], tab[n][2 * i + 1]
            e = abs(c - ct) + abs(d - dt)
            tot[i] += e; mx[i] = max(mx[i], abs(c - ct), abs(d - dt))
        else:
            tail[i] += abs(c) + abs(d)
print('sum_{n<=100} |exact-table| : list1 %.3e list2 %.3e   (paper bound for all n: 3e-5)' % tuple(tot))
print('max entry |exact-table|     : list1 %.3e list2 %.3e' % tuple(mx))
print('exact coefficient mass n>100: list1 %.3e list2 %.3e' % tuple(tail))
print('total error incl. tail      : list1 %.3e list2 %.3e' % (tot[0] + tail[0], tot[1] + tail[1]))
for n in (1, 3, 4, 7, 100, 103, 108, 109, 120, 121):
    k = ns.index(n)
    print('  n=%3d exact c1=% .12e d1=% .12e c2=% .12e d2=% .12e' % (n, x1[k], x1[N + k], x2[k], x2[N + k]))
json.dump({'nodes': ns, 'c1': list(x1[:N]), 'd1': list(x1[N:]), 'c2': list(x2[:N]), 'd2': list(x2[N:])},
          open(SCR + '/exact_coeffs.json', 'w'))
# --- operator norms (induced l1 = max column abs sum) ---
def idx(nodeset):
    k = [ns.index(n) for n in nodeset]
    return k + [N + q for q in k]
J = [1, 3, 4]; E = [n for n in ns if 7 <= n <= 100]; E16 = [n for n in E if n <= 16]
def blk(rows, cols): return S[np.ix_(idx(rows), idx(cols))]
def n1(M): return np.abs(M).sum(axis=0).max()
D = blk(J, J)
for eta in (1, -1):
    V = np.linalg.inv(np.eye(6) - eta * D)
    print('eta=%+d ||(I-eta D)^-1||=%.6f  ||(I-eta D)^-1 S_JE||=%.6f' % (eta, n1(V), n1(V @ blk(J, E))))
print('||S_{E cap[7,16],J}|| = %.17f  (paper encloses .08688748182289641591..592)' % n1(blk(E16, J)))
for xcut in (0, 4, 16, 100):
    rows = [n for n in ns if n > xcut]
    print('||S_{>%d,*}|| (cols<=%d) = %.6g' % (xcut, NM, n1(S[np.ix_(idx(rows), list(range(2 * N)))])))
def readW(fn):
    rows = []
    for line in open(SRC + '/build/data/' + fn):
        v = re.findall(r'-?\d+', line)
        if len(v) == 6: rows.append([int(a) * 1e-4 for a in v])
    return np.array(rows)
for eta, fn in ((1, 'W_plus.tex'), (-1, 'W_minus.tex')):
    W = readW(fn)
    Ae = np.eye(6) - eta * D
    U = blk(J, E)
    U79 = blk(J, [7, 9]); U1216 = blk(J, [12, 13, 15, 16]); U19 = blk(J, [n for n in E if n >= 19])
    print('eta=%+d ||W||=%.6f ||I-W A||=%.20f ||W U79||=%.4f ||W U12:16||=%.4f ||W U19:100||=%.4f' % (
        eta, n1(W), n1(np.eye(6) - W @ Ae), n1(W @ U79), n1(W @ U1216), n1(W @ U19)))


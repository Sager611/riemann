# ew_scan.py (explore67): Galerkin spectra of Weil's form Q = QW_lambda (c = lambda^2) and of its Markov part B = Q - pole,
# both parity sectors, on a grid of windows c; floating point at working precision (orientation only, not rigorous).
# B^ev = Q^ev - 2|cosh(x/2)><cosh(x/2)|, B^odd = Q^odd + 2|sinh(x/2)><sinh(x/2)| (pole: +2<f,cosh>^2 even, -2<f,sinh>^2 odd).
# Usage: ew_scan.py N c1 c2 ...  (c decimal or p/q). Output: one JSON line per c.
import sys, os, json, time
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weil_matrix import WeilMat
def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    if '.' in s:
        a, b = s.split('.'); return fmpq(int(a + b), 10**len(b))
    return fmpq(int(s))
def low(T, k):
    E = acb_mat(T).eig(algorithm='approx')
    ev = sorted([e.real.mid() for e in E], key=lambda z: z)
    return ev[:k]
N = int(sys.argv[1])
for cs in sys.argv[2:]:
    t0 = time.time(); cq = parse_c(cs); cf = int(cq.p)/int(cq.q)
    ctx.prec = max(192, int(18*cf) + 160)
    W = WeilMat(cq); out = {'c': cs, 'N': N, 'prec': ctx.prec}
    for sec, ps in (('even', 2), ('odd', -2)):
        T, idx = W.block(N, sec)
        pv = [W.pole_vec(j, sec) for j in idx]; n = len(idx)
        B = arb_mat(T)
        for a in range(n):
            for b in range(n): B[a, b] = T[a, b] - ps*pv[a]*pv[b]
        out['Q' + sec[0]] = [x.str(8, radius=False) for x in low(T, 4)]
        out['B' + sec[0]] = [x.str(8, radius=False) for x in low(B, 4)]
    out['sec'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)

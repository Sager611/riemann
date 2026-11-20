# ew_bisect.py (explore67): Galerkin location of the parity crossing c_x of the second level of B = Q - pole (e2ev(B) = e1odd(B)),
# floating point at working precision (orientation only; the certified bracket is 5/2 < c_x < 11/4, logs/certs.jsonl). Usage: ew_bisect.py N
import sys, os, json
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weil_matrix import WeilMat
ctx.prec = 256; N = int(sys.argv[1])
def lowB(cq, sec, k):
    W = WeilMat(cq); T, idx = W.block(N, sec); n = len(idx); pv = [W.pole_vec(j, sec) for j in idx]; ps = 2 if sec == 'even' else -2
    B = arb_mat(T)
    for a in range(n):
        for b in range(n): B[a, b] = T[a, b] - ps*pv[a]*pv[b]
    return sorted([float(e.real.mid()) for e in acb_mat(B).eig(algorithm='approx')])[:k]
def f(cq):
    e = lowB(cq, 'even', 2); o = lowB(cq, 'odd', 1)
    return e[1] - o[0], e[1], o[0]
lo, hi = fmpq(5, 2), fmpq(11, 4)
flo, fhi = f(lo)[0], f(hi)[0]
for it in range(14):
    mid = (lo + hi)/2; fm = f(mid)
    if (fm[0] > 0) == (flo > 0): lo, flo = mid, fm[0]
    else: hi, fhi = mid, fm[0]
print(json.dumps({'N': N, 'c_x_lo': float(lo), 'c_x_hi': float(hi), 'level_at_hi': f(hi)[1:]}))

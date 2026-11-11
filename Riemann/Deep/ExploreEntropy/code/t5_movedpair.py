# t5_movedpair.py (explore32): parity of the ground state for xi-like forms with a zero pair moved off the line (explore19's model).
# For each (k, delta): bisection for the first window c where min(lambda_even, lambda_odd) < 0; report which sector crosses first and
# the ground-state parity on a grid of positive windows below the threshold.
import sys, json, math, time
import mpmath as mp
from flint import fmpq
import movedpair as M
mp.mp.prec = 200
K = int(sys.argv[1]); deltas = [float(x) for x in sys.argv[2].split(',')]; clo = float(sys.argv[3]); chi = float(sys.argv[4])
zeros = [str(mp.mp.zetazero(i).imag) for i in range(1, K + 2)]
def both(c):
    cq = fmpq(int(round(c*256)), 256)
    le, l0e, N = M.lam_min(cq, K, d, zeros, 'even')
    lo, l0o, _ = M.lam_min(cq, K, d, zeros, 'odd')
    return float(cq), le, lo, l0e, l0o, N
for d in deltas:
    t0 = time.time()
    a, b = clo, chi
    A = both(a); B = both(b)
    if not (min(A[1], A[2]) > 0 and min(B[1], B[2]) < 0):
        print(json.dumps({'k': K, 'delta': d, 'note': 'no bracket', 'lo': A[:5], 'hi': B[:5]}), flush=True); continue
    for _ in range(9):
        m = (a + b)/2; Mv = both(m)
        if min(Mv[1], Mv[2]) > 0: a, A = m, Mv
        else: b, B = m, Mv
    first = 'even' if B[1] < 0 and (B[2] >= 0 or B[1] < B[2]) else 'odd'
    grid = []
    for t in [0.6, 0.75, 0.85, 0.92, 0.97, 0.995]:
        c = 1 + (a - 1)*t
        G = both(c)
        grid.append({'c': round(G[0], 4), 'even': '%.3e' % G[1], 'odd': '%.3e' % G[2], 'ground': 'even' if G[1] < G[2] else 'odd', 'zeta_e': '%.3e' % G[3], 'zeta_o': '%.3e' % G[4]})
    print(json.dumps({'k': K, 'delta': d, 'threshold': [round(A[0], 5), round(B[0], 5)], 'at_lo': ['%.3e' % A[1], '%.3e' % A[2]], 'at_hi': ['%.3e' % B[1], '%.3e' % B[2]], 'first_negative_sector': first, 'grid_below': grid, 's': round(time.time() - t0, 1)}), flush=True)

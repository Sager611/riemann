# li_range2.py (explore56): universal Li range of the change-of-measure method given zeta's first two Esscher cumulants.
# Class: positive symmetric laws <-> probability nu on [0,inf) (atoms v_i, weights w_i): Q-moments (Q = size-biased law of V=log Y)
#   m_k = sum w_i v_i^k            (k even)
#   m_k = sum w_i v_i^k tanh(v_i/2) (k odd)
# kt_m = cumulants of Q; lambda_n = sum_m C(n,m) kt_m/(m-1)!.  Constraint: kt_1, kt_2 = zeta's values (lambda_1, lambda_2 - 2 lambda_1).
# For 3 atoms the weights solve a 3x3 linear system; scan atoms, keep positive weights, minimize lambda_n.  Float scan, Arb check after.
import math, itertools, sys
from fractions import Fraction
C1 = 0.0230957089661210338; C2 = 0.0461543172958046  # zeta: kt_1, kt_2
def moments(ws, vs, K):
    return [sum(w*(v**k)*(1.0 if k % 2 == 0 else math.tanh(v/2)) for w, v in zip(ws, vs)) for k in range(K+1)]
def cumulants(m, K):
    k = [0.0]*(K+1)
    for n in range(1, K+1):
        k[n] = m[n] - sum(math.comb(n-1, j-1)*k[j]*m[n-j] for j in range(1, n))
    return k
def lam(n, kt):
    return sum(math.comb(n, m)*kt[m]/math.factorial(m-1) for m in range(1, n+1))
def solve3(vs):
    # rows: sum w = 1, sum w v tanh(v/2) = C1, sum w v^2 = C2 + C1^2
    A = [[1.0, 1.0, 1.0], [v*math.tanh(v/2) for v in vs], [v*v for v in vs]]
    b = [1.0, C1, C2 + C1*C1]
    det = lambda M: (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    D = det(A)
    if abs(D) < 1e-300: return None
    ws = []
    for j in range(3):
        Mj = [row[:] for row in A]
        for i in range(3): Mj[i][j] = b[i]
        ws.append(det(Mj)/D)
    return ws
NMAX = 14
best = {n: (float('inf'), None) for n in range(1, NMAX+1)}
grid1 = [i*0.01 for i in range(0, 41)]
grid2 = [0.5*i for i in range(1, 121)]
cnt = 0
for v1 in grid1:
    for v2 in grid2:
        for v3 in grid2:
            if v3 <= v2: continue
            ws = solve3([v1, v2, v3])
            if ws is None or min(ws) < 0: continue
            cnt += 1
            m = moments(ws, [v1, v2, v3], NMAX); kt = cumulants(m, NMAX)
            for n in range(1, NMAX+1):
                L = lam(n, kt)
                if L < best[n][0]: best[n] = (L, (v1, v2, v3, ws))
print('feasible 3-atom laws scanned:', cnt)
for n in range(1, NMAX+1):
    L, arg = best[n]
    print('n=%2d  min lambda_n = %+.6e  at atoms %s weights %s' % (n, L, tuple(round(x, 3) for x in arg[:3]), ['%.3e' % w for w in arg[3]]))

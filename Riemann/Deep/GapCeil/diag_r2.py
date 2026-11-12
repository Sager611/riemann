import math, sys
from flint import arb
import gcert as GC, gcert_r as GR
kap = GC.dec(sys.argv[1]); N = 72; r = 2
gam = int(float(2 * GC.A_(kap) * arb.pi().mid() + 3)) + 1
D = max(50, int(2.5 * gam / 4 + 34))
GC.set_prec(128 + int(D * math.log2(max(gam / 4, 2))) + 32)
ops = GR.OpsR(kap, N, D, gam, r)
for hp in [('a', 0), ('a', 1), ('b', 0), ('c', 0)]:
    A = ops.A(*hp).tolist()
    col = [(abs(float(A[w][0].mid())), w) for w in range(N, ops.Nw)]
    top = sorted(col)[-3:]
    vv = [abs(float(A[k][0].mid())) for k in range(0, 6)]
    print(hp, 'D=%d Nw=%d' % (D, ops.Nw), 'col0 top W entries', ['%.2e@%d' % t for t in top], ' low entries', ['%.2e' % x for x in vv], flush=True)
Q = ops.QT.tolist()
print('max |QT| row Nw-1: %.3e ; row N: %.3e' % (max(abs(float(x.mid())) for x in Q[-1]), max(abs(float(x.mid())) for x in Q[N])))

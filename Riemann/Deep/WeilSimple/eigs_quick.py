# eigs_quick.py: approximate smallest eigenvalues of the even / odd blocks of QW_lambda^N (floating point at the working
# precision; orientation only, nothing rigorous). Usage: eigs_quick.py c N prec
import sys, time
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weil_matrix import WeilMat
def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    return fmpq(int(s))
cq = parse_c(sys.argv[1]); N = int(sys.argv[2]); ctx.prec = int(sys.argv[3]); t0 = time.time()
W = WeilMat(cq); out = {}
for sector in ('even', 'odd'):
    T, idx = W.block(N, sector)
    E = acb_mat(T).eig(algorithm="approx")
    ev = sorted([e.real.mid() for e in E], key=lambda z: float(z) if abs(float(z)) > 0 else 0.0)
    ev = sorted(ev, key=lambda z: z.mid())
    out[sector] = ev
print('c', str(cq), 'N', N, 'prec', ctx.prec, 'even', [e.str(6, radius=False) for e in out['even'][:4]], 'odd', [e.str(6, radius=False) for e in out['odd'][:3]], round(time.time() - t0, 1), 's', flush=True)

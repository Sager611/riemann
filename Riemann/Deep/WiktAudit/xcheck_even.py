# Independent check of the W-W source formulas (source.py: H = (a/2) E A E with the Green-identity displacement form)
# against our ChekAudit x-space Galerkin code (galerkin_x.py: Weil's regularisation + closed-form autocorrelations),
# even sector, at the certificate window a = 23029/10000 (L = 2a, all 35 prime powers m < e^{2a}).
import sys, re, time
from flint import arb, ctx
KOD = sys.argv[1]; GX = sys.argv[2]; NB = int(sys.argv[3]) if len(sys.argv) > 3 else 8
ctx.prec = 256
src = open(GX).read()
src = src.replace("L = arb(c).log(); pi = arb.pi()", "L = c; pi = arb.pi()")
src = src.replace("while n < c:", "while n < L.exp():")
src = src.split('if __name__ == "__main__":')[0]
ns = {}; exec(src, ns)
a = arb(23029)/10000
t0 = time.time()
M = ns['build'](2*a, NB-1, 40, 60)        # physical form A on basis N_j cos(pi j x/a), j=0..NB-1
sys.path.insert(0, KOD)
from source import Source, matrix
S = Source('2.3029', 0, 256)
H = matrix(S.rows(NB))
worst = arb(0); worst_rel = 0.0
for q in range(NB):
    for r in range(NB):
        mine = a/2*(1 if (q+r) % 2 == 0 else -1)*M[q, r]
        d = abs(mine - H[q, r])
        rel = float((d/abs(H[q, r])).mid())
        if rel > worst_rel: worst_rel = rel
        if q <= 3 and r <= 3: print(q, r, 'WW', H[q, r].str(18, radius=False), 'ours', mine.str(18, radius=False), 'rel %.2e' % rel)
print('even sector, %dx%d block: max relative difference %.3e  (%.1f s)' % (NB, NB, worst_rel, time.time()-t0))

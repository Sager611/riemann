# Small-scale end-to-end test of the W-W pipeline algebra at the real window a = 2.3029:
# their core.py (congruence + certified solves), pipeline.full_postcore (rational transfer (5.1)-(5.6), moment branch
# (6.2)-(6.3) for q >= 1024), pipeline.assemble (distant convolutions (7.2)-(7.4), binary64 output) are compared with
# DENSE direct Schur complements computed here from the (independently verified) source matrix.
import sys, json, time, struct
from flint import arb, arb_mat, ctx
KOD = sys.argv[1]; par = int(sys.argv[2]); K, N, M = 24, 512, 4096
sys.path.insert(0, KOD)
from source import Source, matrix
from core import certify
from pipeline import full_postcore, assemble
import tail
t0 = time.time()
core = certify('2.3029', par, K, 1000)
assert core['status'] == 'POTWIERDZONE', core
print('core K=%d gap %s eta %s (%.1fs)' % (K, core['gap'][:30], core['eta'][:30], time.time()-t0), flush=True)
U, V, D, _ = full_postcore(core, M)
print('postcore done (%.1fs)' % (time.time()-t0), flush=True)
# dense direct at 1000 bits
ctx.prec = 1000
src = Source('2.3029', par, 1000)
rowsK = src.rows(K); M0 = matrix(rowsK)
def col(q):   # (A_iq)_{i<K}
    x, u, v, d = src.row(q)
    return arb_mat(K, 1, [(ri[1]*v - ri[2]*u)/(ri[0]-x) for ri in rowsK]), (x, u, v, d)
uP = arb_mat(K, 1, [r[1] for r in rowsK]); vP = arb_mat(K, 1, [r[2] for r in rowsK])
aS = M0.solve(uP); bS = M0.solve(vP)
def post_direct(q):
    c, (x, u, v, d) = col(q)
    Uq = u - (c.transpose()*aS)[0, 0]; Vq = v - (c.transpose()*bS)[0, 0]
    Sqq = d - (c.transpose()*M0.solve(c))[0, 0]
    return Uq, Vq, Sqq
worst = 0.0
for q in (K, K+1, 100, 511, 512, 1023, 1024, 1025, 2047, 3000, M-1):
    Uq, Vq, Sqq = post_direct(q)
    i = q-K
    rel = max(float((abs(Uq-U[i])/abs(Uq)).mid()), float((abs(Vq-V[i])/abs(Vq)).mid()), float((abs(Sqq-D[i])/abs(Sqq)).mid()))
    worst = max(worst, rel)
    print('q=%5d U %s V %s Sqq %s  rel %.2e  (theirs contains direct: %s)' % (q, Uq.str(12, radius=False), Vq.str(12, radius=False), Sqq.str(12, radius=False), rel,
          U[i].overlaps(Uq) and V[i].overlaps(Vq) and D[i].overlaps(Sqq)), flush=True)
print('post-core generators/diagonal: worst rel %.2e' % worst, flush=True)
# assembly with a fake (positive) distant diagonal t_q = dfree(q) + 10 ; the algebra does not depend on t's values
ctx.prec = 192
far = {'N': N, 'M': M, 'rho': arb(0), 'reserve': arb(-10), 'vM': arb(1), 'aM': arb(1), 'bM': arb(1), 'J2M': arb(1), 'J4M': arb(1), 'tM': arb(1)}
out = '/tmp/wiktaudit_test_%d.bin' % par
rec = assemble(core, far, U, V, D, out)
raw = open(out, 'rb').read(); n = N-K
assert len(raw) == n*40
cols = struct.unpack('<%dd' % (5*n), raw)
row = lambda i: cols[5*i:5*i+5]
ctx.prec = 1000
srcd = Source('2.3029', par, 192)
tq = {q: tail.diagonal_free(srcd, q) + 10 for q in range(N, M)}
def S_direct(i, q):   # entry of the post-core Schur complement S at (i,q), i != q, dense
    ci, (xi, ui, vi, di) = col(i); cq, (xq, uq, vq, dq) = col(q)
    return (ui*vq - vi*uq)/(xi - xq) - (ci.transpose()*M0.solve(cq))[0, 0]
pairs = [(K, K), (K, K+1), (K+5, 300), (200, 201), (N-1, N-1), (K, N-1), (300, 300)]
cache = {}
for (i, j) in pairs:
    for x in (i, j):
        if x not in cache:
            cx, _ = col(x); sol = M0.solve(cx)
            cache[x] = (cx, sol)
ws = 0.0
for (i, j) in pairs:
    # direct: Smat_ij - sum_{N<=q<M} S_iq S_jq / t_q
    if i == j:
        _, (x, u, v, d) = col(i); base = d - (cache[i][0].transpose()*cache[i][1])[0, 0]
    else:
        base = S_direct(i, j)
    tot = arb(0)
    xi = src.row(i); xj = src.row(j)
    for q in range(N, M):
        cq, (xq, uq, vq, dq) = col(q)
        siq = (xi[1]*vq - xi[2]*uq)/(xi[0]-xq) - (cache[i][1].transpose()*cq)[0, 0]
        sjq = (xj[1]*vq - xj[2]*uq)/(xj[0]-xq) - (cache[j][1].transpose()*cq)[0, 0]
        tot += siq*sjq/tq[q]
    direct = base - tot
    ri, rj = row(i-K), row(j-K)
    if i == j: theirs = ri[0]
    else:
        lam = lambda q: (q + par/2)**2
        theirs = ((ri[1]*rj[2] - ri[2]*rj[1]) + (ri[3]*rj[4] - ri[4]*rj[3]))/(lam(i) - lam(j))
    rel = abs(float(direct.mid()) - theirs)/abs(float(direct.mid()))
    ws = max(ws, rel)
    print('S_M(%d,%d): direct %s  binary-reconstructed %.15g  rel %.2e' % (i, j, direct.str(15, radius=False), theirs, rel), flush=True)
print('assembled matrix vs dense: worst rel %.2e (binary64 storage); total %.1fs' % (ws, time.time()-t0))

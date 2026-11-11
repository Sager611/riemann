# spec_small.py: low spectrum of the Galerkin truncation QW_lambda^N (both iota-sectors) at c = lambda^2, from the closed-form
# Arb matrix of ../../WeilSimple/weil_matrix.py (same space E_N as ACLG App. D). Eigenvalues: acb_mat.eig(approx) at prec bits
# (midpoints; floating point, orientation only; min-max: each is an upper bound of the true eigenvalue of the same index).
# Usage: spec_small.py c N prec [k]   -> JSON line {c, N, prec, even: [...k], odd: [...k]}
import sys, json, time, os
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'WeilSimple'))
from weil_matrix import WeilMat

def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    return fmpq(int(s))

if __name__ == '__main__':
    cq = parse_c(sys.argv[1]); N = int(sys.argv[2]); ctx.prec = int(sys.argv[3]); k = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    t0 = time.time(); W = WeilMat(cq); out = {'c': sys.argv[1], 'N': N, 'prec': ctx.prec}
    for sector in ('even', 'odd'):
        T, idx = W.block(N, sector)
        E = acb_mat(T).eig(algorithm="approx")
        ev = sorted([e.real.mid() for e in E], key=lambda z: z)
        out[sector] = [e.str(30, radius=False) for e in ev[:k]]
    out['sec'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)

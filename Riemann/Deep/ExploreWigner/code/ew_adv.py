# ew_adv.py (explore67): adversarial check of odd-sector interlacing E1odd <= E2ev (and parity alternation) for the Weil-type forms of
# zeta, L(s, chi_5) and Davenport-Heilbronn (explore19's wmat.py: gamma factor, conductor 5, atoms; DH and chi5 have no pole and signed
# prime weights, so no Markov part). Galerkin, floating point at working precision (orientation only). Usage: ew_adv.py N kind c1 c2 ...
import sys, os, json, time
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wmat import Form
def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    return fmpq(int(s))
N = int(sys.argv[1]); kind = sys.argv[2]
for cs in sys.argv[3:]:
    t0 = time.time(); cq = parse_c(cs); cf = int(cq.p)/int(cq.q); ctx.prec = max(192, int(18*cf) + 160)
    F = Form(cq, kind); out = {'kind': kind, 'c': cs, 'N': N, 'prec': ctx.prec}
    for sec in ('even', 'odd'):
        T = F.block(N, sec); T = T[0] if isinstance(T, tuple) else T
        ev = sorted([e.real.mid() for e in acb_mat(T).eig(algorithm='approx')])
        out[sec] = [x.str(6, radius=False) for x in ev[:3]]
    e, o = [float(x) for x in out['even']], [float(x) for x in out['odd']]
    out['alt'] = bool(e[0] < o[0] < e[1] < o[1] < e[2]); out['interlace'] = bool(o[0] <= e[1]); out['sec'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)

# o1fix.py (ccmfmm): approximate lowest odd eigenvalue o1 of the truncated odd block at N (inverse iteration + Arb Rayleigh quotient,
# an upper bound for lambda_1^odd), used only to place L_odd and the chain steps. (The NO = 120 values of the first run_ccm.py rows
# 21..25 were unconverged: 10x to 1e4x too high at c = 22..25; NO = NG since.)  Usage: o1fix.py <out.jsonl> N bits c1 c2 ...
import sys, os, json, time
from flint import arb, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weil_matrix import WeilMat
from gstate import ground_vec, rayleigh, parse_c
out = sys.argv[1]; N = int(sys.argv[2]); ctx.prec = int(sys.argv[3])
for cs in sys.argv[4:]:
    t0 = time.time(); W = WeilMat(parse_c(cs)); T, _ = W.block(N, 'odd')
    g = ground_vec(T, N); num, gg = rayleigh(T, g); o1 = num/gg
    rec = dict(c=cs, N=N, prec=ctx.prec, o1=o1.mid().str(12, radius=False), o1_ball=o1.str(10), sec=round(time.time() - t0, 1))
    open(out, 'a').write(json.dumps(rec) + '\n'); print(json.dumps(rec), flush=True)

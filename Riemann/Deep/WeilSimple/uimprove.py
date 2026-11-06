# uimprove.py: rigorous upper bounds U(c, N) = Q(g)/||g||^2 (g = inverse-iteration ground vector of the even block of
# QW_lambda^N, exact dyadics; Arb matrix from weil_matrix.py). Usage: uimprove.py out.jsonl c:N:bits [c:N:bits ...]
import sys, os, json, time
from flint import arb, ctx, fmpq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from weil_matrix import WeilMat
from cert_simple import ground_vec, rayleigh, parse_c
out = sys.argv[1]
for spec in sys.argv[2:]:
    cs, N, bits = spec.split(':'); N = int(N); ctx.prec = int(bits); t0 = time.time()
    cq = parse_c(cs); W = WeilMat(cq); T, _ = W.block(N, 'even')
    g = ground_vec(T, N + 1); num, gg = rayleigh(T, g); U = num/gg
    rec = dict(c=cs, N=N, bits=int(bits), U_up=arb(U.upper()).str(15, radius=False), U=U.str(10, radius=True), sec=round(time.time() - t0, 1))
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)

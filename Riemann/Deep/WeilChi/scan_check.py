# scan_check.py: compare the C pre-scan (cscan.c) with the cached numpy pre-scan (both followed by the same Arb confirmation)
import sys, os, json, time
from flint import ctx, arb
import cert_chi as C
key, cs, K2 = sys.argv[1], sys.argv[2], int(sys.argv[3])
ctx.prec = 128
cq = C.parse_c(cs); S = C.SymChi(key, cq); J = K2
ks = C.scan_end(S, K2, 1); kse = C.scan_end(S, 2*J, 0)
t0 = time.time()
fo, fe = C.cscan(S, K2, J, ks, kse, os.environ['WEILCHI_CSCAN'])
t1 = time.time()
ctx.prec = 96; S96 = C.SymChi(key, cq)
selo = [k for k in fo if S96.mu(k).lower() < 1]
sele = [j for j in fe if S96.mu(2*j).lower() < 0]
cache = json.load(open(os.path.join(os.environ['WEILCHI_SCAN'], f"scan_{key}_{cs.replace('/', 'o')}_{K2}.json")))
print(f"{key} c={cs} K2={K2} ks={ks} kse={kse} cscan {t1-t0:.0f}s flagged {len(fo)}/{len(fe)} (numpy {cache['flagged']}); selected {len(selo)}/{len(sele)} (numpy {len(cache['selo'])}/{len(cache['sele'])}); identical: {selo == cache['selo'] and sele == cache['sele']}", flush=True)

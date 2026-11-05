# mu_dist.py: distribution of mu(k) over the selected odd nodes (scan cache) -- how many would stay selected for a smaller t1
import sys, json, os
from flint import ctx
import cert_chi as C
key, cs, K2 = sys.argv[1], sys.argv[2], int(sys.argv[3])
ctx.prec = 96
S = C.SymChi(key, C.parse_c(cs))
d = json.load(open(os.path.join(os.environ['WEILCHI_SCAN'], f"scan_{key}_{cs}_{K2}.json")))
vals = sorted(float(S.mu(k).mid()) for k in d['selo'])
for t in (1, 0.75, 0.5, 0.25, 0.1, 0.0, -1.0):
    print(f"t1={t}: {sum(v < t for v in vals)} of {len(vals)} selected odd nodes have mu < t1")
vals2 = sorted(float(S.mu(2*j).mid()) for j in d['sele'])
print('even selected: min mu', vals2[0], ' # < -1:', sum(v < -1 for v in vals2), 'of', len(vals2))

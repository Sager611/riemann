import json, sys
from fractions import Fraction as F
src, dst = sys.argv[1], sys.argv[2]
KEYS = ['Z', 'psi', 'psi2', 'd0', 'd1', 'd2', 'r0', 'r1', 'r2']
rows = [json.loads(l) for l in open(src)]
with open(dst, 'w') as f:
    f.write("# i,t_i,d_i,1e10*q_i,1e10*q_{i+1},1e8*C_i,1e8*G_i,1e6*lambda_i,1e6*mu_i," + ",".join("1e12*U_" + k for k in KEYS) + "\n")
    for r in rows:
        q = F(r['q']); qn = F(r['qnext']); C = F(r['C']); G = F(r['G']); la = F(r['lam']); mu = F(r['mu'])
        def I(x, den):
            y = x * den; assert y.denominator == 1, (x, den); return str(y.numerator)
        f.write(",".join([str(r['i']), r['t'], r['d'], I(q, 10**10), I(qn, 10**10), I(C, 10**8), I(G, 10**8), I(la, 10**6), I(mu, 10**6)] + [str(r['U'][k]) for k in KEYS]) + "\n")
print("wrote", len(rows), "rows")

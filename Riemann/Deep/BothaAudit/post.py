# post.py (bothaaudit): eigenvalues of deig.c's Rayleigh-Ritz matrices (mpmath, 60 digits); appends one JSON record.
# Usage: post.py out.jsonl c N prec even.out odd.out [extra-json]
import sys, json
from mpmath import mp, mpf, matrix, eigsy
mp.dps = 60
def parse(fn):
    info = {}; H = None; its = []
    for l in open(fn):
        w = l.split()
        if not w: continue
        if w[0] == 'n': info = dict(n=int(w[1]), neg=int(w[3]), minpiv_log10=float(w[5]), fill_s=float(w[7]), ldl_s=float(w[9]))
        elif w[0] == 'it': its.append([float(z) for z in w[2:]])
        elif w[0] == 'H': H = [mpf(z) for z in w[1:]]
        elif w[0] == 'sec': info['sec'] = float(w[1])
    k = int(round(len(H) ** 0.5)); M = matrix(k, k)
    for a in range(k):
        for b in range(k): M[a, b] = (H[a*k + b] + H[b*k + a])/2
    E, Q = eigsy(M); ev = sorted([E[i] for i in range(k)])
    info['ev'] = ev; info['conv'] = [round(r[0] - its[-1][0], 10) for r in its]
    return info
out, c, N, prec = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
e = parse(sys.argv[5]); o = parse(sys.argv[6])
f = lambda z: mp.nstr(z, 12)
rec = dict(c=c, N=N, prec=prec, e1=f(e['ev'][0]), e2=f(e['ev'][1]), e3=f(e['ev'][2]) if len(e['ev']) > 2 else None, o1=f(o['ev'][0]),
           o2=f(o['ev'][1]) if len(o['ev']) > 1 else None, r_o1e1=float(o['ev'][0]/e['ev'][0]), r_e2e1=float(e['ev'][1]/e['ev'][0]),
           neg=[e['neg'], o['neg']], conv_e=e['conv'], conv_o=o['conv'], sec=[e['sec'], o['sec']], n=[e['n'], o['n']])
if len(sys.argv) > 7: rec.update(json.loads(sys.argv[7]))
open(out, 'a').write(json.dumps(rec) + '\n'); print(json.dumps(rec), flush=True)

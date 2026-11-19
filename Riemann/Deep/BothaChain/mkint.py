# mkint.py (bothachain): step heuristics logs/chain_int.jsonl = bothaaudit's rows c = 30..35 (unchanged) + rows 36..56:
# o1 := 0.6 * (float-scan odd Galerkin value Uo at Botha's windows, ../BothaAudit/logs/points.jsonl, log10-interpolated), bits := run_ccm.bitsrule
# of the interpolated scan e1 with factor BITSF (env; default 1.6 = run_ccm.bitsrule) for rows >= 36; every row gets e2 := 0.6 * scan U2
# (chain_pair.py's test for a shared defl certificate). They choose step sizes and precision only; no rigour.
import json, math, os
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(H, '..', 'BothaAudit', 'logs')
def lg(s): m, _, e = s.partition('e'); return math.log10(float(m)) + int(e)
BF = float(os.environ.get('BITSF', '1.6'))
def bitsrule(le1): return max(192, 64*math.ceil((BF*(-le1)*math.log2(10) + 64)/64))
def e10(x): return '%.6fe%d' % (10**(x - math.floor(x)), math.floor(x))
P = [json.loads(l) for l in open(os.path.join(A, 'points.jsonl'))]
xs = [float(d['c']) for d in P]
def interp(x, key):
    ys = [lg(d[key]) for d in P]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]: t = (x - xs[i])/(xs[i + 1] - xs[i]); return (1 - t)*ys[i] + t*ys[i + 1]
    i = 0 if x < xs[0] else len(xs) - 2   # c = 30 and c = 55, 56: linear extrapolation from the nearest two windows
    return ys[i] + (x - xs[i])*(ys[i + 1] - ys[i])/(xs[i + 1] - xs[i])
rows = [json.loads(l) for l in open(os.path.join(A, 'chain_int.jsonl'))]
for d in rows: d['approx']['e2'] = e10(interp(int(d['c']), 'U2_up') + math.log10(0.6))
for n in range(36, 57):
    lo = interp(n, 'Uo_up') + math.log10(0.6); le = interp(n, 'U_up')
    rows.append(dict(c=str(n), bits=bitsrule(le), approx=dict(o1=e10(lo), e2=e10(interp(n, 'U2_up') + math.log10(0.6))), src='0.6*scan Uo, U2; BITSF %g' % BF))
with open(os.path.join(H, 'logs', 'chain_int.jsonl'), 'w') as f:
    for d in rows: f.write(json.dumps(d) + '\n')
print(len(rows), rows[-1])

# a401172.py: certified intervals L <= kappa(n) <= U and a(n) = floor(-log10 kappa(n)) (OEIS A401172) for n = 2..20 from the logs.
# kappa(n) = inf Q(f)/||f||^2 over f supported in (-log(n)/2, log(n)/2) = lambda_1^even(n), since lambda_1^odd >= L_odd > U (REPORT sec. 4).
# Lower bounds: plain even certificates (logs/bracket.jsonl, logs/abracket.jsonl); upper bounds: Rayleigh quotients (logs/int.jsonl U_up,
# logs/uimprove.jsonl). Usage: python3 a401172.py   (prints the markdown table)
import json, os
import mpmath as mp
mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
def rd(f): return [json.loads(l) for l in open(os.path.join(HERE, 'logs', f))] if os.path.exists(os.path.join(HERE, 'logs', f)) else []
low, up = {}, {}
def setlow(n, v, src):
    if n not in low or mp.mpf(v) > mp.mpf(low[n][0]): low[n] = (v, src)
def setup(n, v, src):
    if n not in up or mp.mpf(v) < mp.mpf(up[n][0]): up[n] = (v, src)
for r in rd('bracket.jsonl'):
    if r.get('L1'): setlow(int(r['c']), r['L1'], 'even cert, K2 = max(401, 5n^2/3)')
for r in rd('abracket.jsonl'):
    if r.get('rigorous'): setlow(int(r['c']), r['eps'], 'even cert, K2 = %d' % r['K2'])
for r in rd('int.jsonl'): setup(int(r['c']), r['U_up'], 'g in E_%d' % r['NG'])
for r in rd('uimprove.jsonl'):
    if '/' not in r['c']: setup(int(r['c']), r['U_up'], 'g in E_%d' % r['N'])
oeis = {2: 2, 3: 7, 4: 12, 5: 17, 6: 22, 7: 27, 8: 32, 9: 37}
print('| n | certified interval for kappa(n) | -log10 kappa(n) in | a(n) | OEIS data | sources (lower; upper) |')
print('|---|---|---|---|---|---|')
for n in range(2, 21):
    L, ls = low[n]; U, us = up[n]
    lo, hi = -mp.log10(mp.mpf(U)), -mp.log10(mp.mpf(L))
    a = int(mp.floor(lo)) if int(mp.floor(lo)) == int(mp.floor(hi)) else None
    astr = str(a) if a is not None else '%d or %d (open)' % (int(mp.floor(lo)), int(mp.floor(hi)))
    print(f"| {n} | [{mp.nstr(mp.mpf(L), 6)}, {mp.nstr(mp.mpf(U), 6)}] | [{mp.nstr(lo, 9)}, {mp.nstr(hi, 9)}] | {astr} | {oeis.get(n, '-')} | {ls}; {us} |")

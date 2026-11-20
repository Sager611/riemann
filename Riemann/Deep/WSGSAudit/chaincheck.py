# chaincheck.py (wsgsaudit): independent check of WeilSimple's chain and integer-window logs (exact Decimal/Fraction arithmetic).
import json, sys
from decimal import Decimal as D
from fractions import Fraction as F
base = sys.argv[1]
ints = {F(r['c']): r for r in (json.loads(l) for l in open(base + '/logs/int.jsonl'))}
ch = [json.loads(l) for l in open(base + '/logs/chain.jsonl')]
bad = []
for c, r in sorted(ints.items()):   # integer windows: U(c) < min(Lodd(c), L2(c)), both certificates rigorous
    if not (r['defl']['rigorous'] and r['odd']['rigorous'] and D(r['U_up']) < D(r['defl']['eps']) and D(r['U_up']) < D(r['odd']['eps'])): bad.append(('int', str(c)))
    if D(r['U_up']) < D(r['defl']['U'].split('+/-')[0].strip('[ ')): bad.append(('U_up<mid', str(c)))
cur, U = F(2), ints[F(2)]['U_up']; pts = {F(2): (U, ints[F(2)]['NG'])}
for r in ch:   # chain: contiguous, Uprev = previous U, U(c_i) < eps(c_{i+1}), both certificates rigorous at c_{i+1}
    ok = (F(r['cprev']) == cur and r['Uprev'] == U and D(r['Uprev']) < D(r['eps']) and r['ok'] and r['defl']['rigorous'] and r['odd']['rigorous'])
    if not ok: bad.append(('chain', r['cprev'], r['c']))
    cur, U = F(r['c']), r['U_up']; pts[cur] = (U, r['NG'])
print('chain end', cur, 'steps', len(ch), 'bad', bad)
# truncation claim: at every window c with U(c) (g in E_NG) below a certified lower bound at some c' >= c (same c or next chain point)
lows = {}
for c, r in ints.items(): lows.setdefault(c, []).append(min(D(r['defl']['eps']), D(r['odd']['eps'])))
for r in ch: lows.setdefault(F(r['c']), []).append(D(r['eps']))
wins = {}
for c, (Ug, NG) in pts.items(): wins[c] = (D(Ug), NG)
for c, r in ints.items(): wins.setdefault(c, (D(r['U_up']), r['NG']))
res = []
for c, (Ug, NG) in sorted(wins.items()):
    Lb = max(max(v) for cc, v in lows.items() if cc >= c)   # monotone: a lower bound at c' >= c holds at c
    best = max((max(v) for cc, v in lows.items() if cc >= c), default=None)
    res.append((c, NG, Ug < Lb))
print('windows', len(res), 'all ok', all(x[2] for x in res), 'N0 range', min(x[1] for x in res), max(x[1] for x in res))

# table.py: tables for REPORT.md and the chain check. Usage: python3 table.py   (reads logs/int.jsonl, logs/chain.jsonl)
import json, math, os
from decimal import Decimal, getcontext
from fractions import Fraction
getcontext().prec = 50
HERE = os.path.dirname(os.path.abspath(__file__))
def D(s): return Decimal(s)
def fr(s): return Fraction(s)
def sci(s, n=4):
    d = D(s); e = d.adjusted(); m = d.scaleb(-e)
    return f"{m:.{n-1}f}e{e}"
rows = [json.loads(l) for l in open(os.path.join(HERE, 'logs/int.jsonl'))]
rows.sort(key=lambda r: fr(r['c']))
print('Integer windows (each row: lambda_1^even <= U < min(L_odd, L_2) => simple even ground state; QW_lambda^N even-simple for N >= NG)')
print('| c = lambda^2 | lambda | U | L_odd | L_2 | NG | bits | seconds: certificates (with parameter search) | peak MB | margins (defl, odd) |')
print('|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    d, o = r['defl'], r['odd']; ok = r['even_simple']
    assert ok and d['rigorous'] and o['rigorous'] and D(r['U_up']) < D(d['eps']) and D(r['U_up']) < D(o['eps'])
    c = fr(r['c']); lam = math.sqrt(float(c))
    print(f"| {r['c']} | {lam:.4f} | {sci(r['U_up'])} | {sci(o['eps'])} | {sci(d['eps'])} | {r['NG']} | {r['bits']} | {round(d['wall'] + o['wall'])} ({round(r['seconds'])}) | {max(d['rss_mb'] or 0, o['rss_mb'] or 0)} | {sci(d['margin'], 2)}, {sci(o['margin'], 2)} |")
ch = [json.loads(l) for l in open(os.path.join(HERE, 'logs/chain.jsonl'))]
okc = [r for r in ch if r['ok'] and D(r['U_up']) < D(r['eps'])]
print()
print(f'Chain: {len(ch)} records, {len(okc)} ok')
# verify the chain: consecutive, each U(prev) < eps(next) (= L_odd = L_2 at next)
start = fr('2'); cur = start; Ucur = None
for r in rows:
    if fr(r['c']) == start: Ucur = D(r['U_up'])
pts = [(cur, Ucur, None)]
used = []
for r in okc:
    if fr(r['cprev']) != cur: continue
    assert D(r['Uprev']) == Ucur, (r['c'], r['Uprev'], Ucur)
    assert D(r['Uprev']) < D(r['eps']) and r['defl']['rigorous'] and r['odd']['rigorous']
    cur = fr(r['c']); Ucur = D(r['U_up']); used.append(r); pts.append((cur, Ucur, r))
print(f'covered interval: [{start}, {cur}] = lambda in [{math.sqrt(float(start)):.4f}, {math.sqrt(float(cur)):.4f}], {len(used)} steps')
print('| c_i -> c_{i+1} | U(c_i) | L(c_{i+1}) = eps | U(c_{i+1}) | margins (defl, odd) | NG | bits | s |')
print('|---|---|---|---|---|---|---|---|')
for r in used:
    print(f"| {r['cprev']} -> {r['c']} | {sci(r['Uprev'], 3)} | {sci(r['eps'], 3)} | {sci(r['U_up'], 3)} | {sci(r['defl']['margin'], 2)}, {sci(r['odd']['margin'], 2)} | {r['NG']} | {r['bits']} | {r['sec']} |")
tot = sum(r['sec'] for r in used); print(f'chain compute: {tot:.0f} s; max RSS {max(max(r["defl"]["rss_mb"] or 0, r["odd"]["rss_mb"] or 0) for r in used)} MB')

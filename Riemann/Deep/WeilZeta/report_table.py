# report_table.py: the curated certificate table of REPORT.md sec. 2 (one row per c and parity; latest successful run wins),
# with Chekmachev's own margin at his eps from his manifests (.scratch/weilzeta/chek, Zenodo 23847003) where our eps equals his.
import json, glob, re, math, os
from decimal import Decimal
runs = {}
for f in ['logs/val.log', 'logs/z20.log', 'logs/z25.log', 'logs/z29.log']:
    if not os.path.exists(f): continue
    tt = open(f[:-4] + '.time').read() if os.path.exists(f[:-4] + '.time') else ''
    reals = re.findall(r'([0-9.]+) real', tt); rss = re.findall(r'(\d+)\s+maximum resident set size', tt)
    js = [json.loads(l) for l in open(f) if l.startswith('{"char"')]
    # time entries include killed runs (no JSON line); align from the end of each file only if counts match
    for i, d in enumerate(js):
        key = (int(d['c']), d['sector'])
        tr = (reals[i], rss[i]) if len(reals) == len(js) else (None, None)
        runs[key] = (d, tr)
chek = '/Users/adriansager/playground/riemann/.scratch/weilzeta/chek'
print('| c | half-width a | parity of f | eps | margin (ours) | Chekmachev, same eps | bits | order / U | time (s) | peak RSS (GB) | status |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for (c, s), (d, (rr, mm)) in sorted(runs.items()):
    if not d['rigorous']: continue
    their = ''
    mf = f'{chek}/certificates__c{c}_{s}_eps__manifest.json'
    if os.path.exists(mf):
        m = json.load(open(mf)); e1 = Decimal(m['eps'])
        try: e0 = Decimal(d['eps'])
        except Exception: e0 = None
        if e0 is not None and abs(e0 - e1) <= Decimal('1e-9')*e1: their = f"{m['margin_lower']:.6e} (eps {float(e1):.10e})"
    t = f"{float(rr):.0f}" if rr else '?'; g = f"{int(mm)/1e9:.2f}" if mm else '?'
    print(f"| {c} | {math.log(c)/2:.4f} | {s} | {d['eps']} | {d['margin']} | {their} | {d['prec']} | {d['order']} / {d['nU']} | {t} | {g} | {'rigorous (Arb + proved scan)' if d.get('scan') == 'rigorous' else 'float scan budget'} |")

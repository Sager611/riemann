# table.py: collect cert_chi.py JSON lines from logs/*.log into a markdown table
import json, glob, sys
rows = []
for f in sorted(glob.glob('logs/*.log')):
    for line in open(f):
        line = line.strip()
        if not line.startswith('{'): continue
        try: d = json.loads(line)
        except Exception: continue
        rows.append((f.split('/')[-1], d))
print('| log | chi | c | half-width | sector | eps | bits | K2 | kscan | sel odd/even | K+ | K- | J- | G | order | margin | rigorous | max rad | time (s) |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
import math
for f, d in rows:
    p, q = (d['c'].split('/') + ['1'])[:2]
    hw = math.log(int(p)/int(q))/2
    print(f"| {f} | {d['char']} | {d['c']} | {hw:.4f} | {d['sector']} | {d['eps']} | {d['prec']} | {d['K2']} | {d['kscan']} | {d['sel_odd']}/{d['sel_even']} | {d['Kp']} | {d['Km']} | {d['Jm']} | {d['G']} | {d['order']} | {d['margin']} | {d['rigorous']} | {d['maxrad']} | {d['times'][-1]} |")

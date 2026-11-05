# table.py: cert_zeta.py JSON lines (logs/*.log) joined with /usr/bin/time -l output (logs/*.time, same order) -> markdown table
import json, glob, re, math, sys
rows = []
for f in sorted(glob.glob('logs/*.log')):
    tf = f[:-4] + '.time'
    try: tt = open(tf).read()
    except OSError: tt = ''
    reals = re.findall(r'([0-9.]+) real', tt); rss = re.findall(r'(\d+)\s+maximum resident set size', tt)
    js = []
    for line in open(f):
        line = line.strip()
        if line.startswith('{"char"'):
            try: js.append(json.loads(line))
            except Exception: pass
    for i, d in enumerate(js):
        rr = f"{float(reals[i]):.0f}" if i < len(reals) else '?'
        mm = f"{int(rss[i])/1e9:.2f}" if i < len(rss) else '?'
        rows.append((f.split('/')[-1], d, rr, mm))
print('| log | c | half-width | sector | eps | bits | K2 | kscan odd/even | sel odd/even | K+ / K- / J- / G | order (N) / U | margin | rigorous | resid | time (s) | peak RSS (GB) |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for f, d, rr, mm in rows:
    p, q = (d['c'].split('/') + ['1'])[:2]
    hw = math.log(int(p)/int(q))/2
    print(f"| {f} | {d['c']} | {hw:.4f} | {d['sector']} | {d['eps']} | {d['prec']} | {d['K2']} | {d['kscan']}/{d['kscan_e']} | {d['sel_odd']}/{d['sel_even']} | {d['Kp']} / {d['Km']} / {d['Jm']} / {d['G']} | {d['order']} / {d.get('nU')} | {d['margin']} | {d['rigorous']} | {d.get('resid')} | {rr} | {mm} |")

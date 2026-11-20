#!/usr/bin/env python3
# dbnrec summarize.py: markdown tables from logs (no computation): per row certified fields; per shell time and footprint.
import json, glob, re, os
L = 'Riemann/Deep/DBNRec/logs/'
print('| row | shells | u_F | T (shell-0 binding) | D | A′ (cert) | B′ (cert) | A (Raymond) | B (Raymond) | dA | dB | Raymond credit C+g | bb CPU s | peak MB (memcap / max RSS) |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for cf in sorted(glob.glob(L + 'cert_r*.json'), key=lambda p: int(re.findall(r'cert_r(\d+)', p)[0])):
    try: c = json.load(open(cf))
    except Exception: continue
    r = c['row']; sh = c['shells']; secs = 0; pk = 0; rss = 0
    for s in sh:
        t = open(L + s['tag'] + '.bb.out').read()
        m = re.findall(r'\| (\d+) s\n', t); secs += int(m[-1]) if m else 0
        m = re.findall(r'memcap: peak (\d+) MB', t); pk = max([pk] + [int(x) for x in m])
        m = re.findall(r'(\d+)\s+maximum resident set size', t); rss = max([rss] + [int(x)//2**20 for x in m])
    print('| %d | %d | %s | %.6f | %.6f | %.9f | %.9f | %.9f | %.9f | %.6f | %.6f | — | %d | %d / %d |' % (r, len(sh), c['far']['uF'], max(s['S'] for s in sh), sh[0]['D'],
          c['A_new_float'], c['B_new_float'], c['A_new_float'] + c['dA'], c['B_new_float'] + c['dB'], c['dA'], c['dB'], secs, pk, rss))
print()
print('| shell | u-range | thr | unresolved cells | S_j (Arb) | D_j | max(UB_Arb − UB_float) |')
print('|---|---|---|---|---|---|---|')
for cf in sorted(glob.glob(L + 'cert_r*.json')):
    try: c = json.load(open(cf))
    except Exception: continue
    for s in c['shells']:
        print('| %s | [%s, %s] | %.6f | %d | %.9f | %.6f | %.2e |' % (s['tag'], s['u0'], s['u1'], s['thr'], s['cells'], s['S'], s['D'], s['max_ub_arb_minus_float']))

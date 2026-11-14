# table_chain.py (bothaaudit): markdown table of logs/chain.jsonl (CCMFMM chain continued from c = 30) and the covered interval.
import json, os
H = os.path.dirname(os.path.abspath(__file__))
R = [json.loads(l) for l in open(os.path.join(H, 'logs', 'chain.jsonl'))]
print('| c_i → c_{i+1} | U(c_i) | ε = L(c_{i+1}) | U(c_{i+1}) | certified ν (defl, odd) | n (defl / odd) | NG | bits | ok | peak RSS MB (g; prep, LDL defl; odd) | s |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
end = '30'; mx = 0
for d in R:
    rd, ro = d['defl'], d['odd']
    rss = [d.get('g_rss')] + (rd.get('rss_mb') or []) + (ro.get('rss_mb') or [])
    mx = max([mx] + [x for x in rss if x])
    print('| %s → %s | %s | %s | %s | %s, %s | %s / %s | %d | %d | %s | %s; %s; %s | %s |' % (d['cprev'], d['c'], d['Uprev'], d['eps'], d['U_up'], rd.get('nu_cert'), ro.get('nu_cert'),
          rd.get('n'), ro.get('n'), d['NG'], d['bits'], 'ok' if d['ok'] else 'FAIL', d.get('g_rss'), rd.get('rss_mb'), ro.get('rss_mb'), d['sec']))
    if d['ok'] and d['cprev'] == end: end = d['c']
from fractions import Fraction
print()
print('covered: [30, %s] = [30, %.6f]; ok steps %d, failed attempts %d; max RSS %d MB' % (end, float(Fraction(end)), sum(1 for d in R if d['ok']), sum(1 for d in R if not d['ok']), mx))

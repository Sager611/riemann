# falsify.py (bothaaudit): falsification scan summary. Ratios r_odd = o1/e1 and r_l2 = e2/e1 of the truncated blocks on Botha's 31 windows in
# (30, 54.6] (logs/points.jsonl: U, Uo, U2 = Rayleigh(-Ritz) values of deig's vectors, i.e. e1, o1, e2 of T_N up to rounding), and the truncation
# study at both grid ends (logs/conv30.jsonl + logs/scan_N240.jsonl at c = 30.570758; logs/conv.jsonl at c = 54.598).
import json, os
H = os.path.dirname(os.path.abspath(__file__)); lg = lambda s: float(s)
P = [json.loads(l) for l in open(os.path.join(H, 'logs', 'points.jsonl'))]
mn_o = min((lg(d['Uo_up'])/lg(d['U_up']), d['c']) for d in P); mn_2 = min((lg(d['U2_up'])/lg(d['U_up']), d['c']) for d in P)
mx_o = max((lg(d['Uo_up'])/lg(d['U_up']), d['c']) for d in P)
print('points', len(P), 'c from', P[0]['c'], 'to', P[-1]['c'])
print('min o1/e1 %.4g at c=%s; max %.4g at c=%s; min e2/e1 %.4g at c=%s' % (mn_o[0], mn_o[1], mx_o[0], mx_o[1], mn_2[0], mn_2[1]))
print('min relative gap (o1-e1)/o1 = %.8f' % (1 - 1/mn_o[0]))
print('c | N | e1 | o1/e1 | e2/e1')
for d in P: print(d['c'], d['N'], d['U_up'], '%.4g' % (lg(d['Uo_up'])/lg(d['U_up'])), '%.4g' % (lg(d['U2_up'])/lg(d['U_up'])))
print('-- truncation study')
for fn in ('scan_N240.jsonl', 'conv30.jsonl', 'conv.jsonl'):
    for l in open(os.path.join(H, 'logs', fn)):
        d = json.loads(l)
        if fn == 'scan_N240.jsonl' and d['c'] != '30.570758': continue
        print(d['c'], d['N'], d['e1'], d['o1'], d['e2'], '%.4g %.4g' % (d['r_o1e1'], d['r_e2e1']))

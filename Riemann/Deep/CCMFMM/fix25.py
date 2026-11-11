# fix25.py (ccmfmm): the c = 25 row of the first run used L_odd from an unconverged o1 (NO = 120) and its odd certificate failed there;
# rerun the odd certificate at L_odd = (e1 o1)^{1/2} with o1 from o1fix.py (N = 200) and replace the row (old row kept in logs/int_failed.jsonl).
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_ccm import cert, lg, e10
cs = sys.argv[1]; lf = 'logs/int.jsonl'
rows = [json.loads(l) for l in open(lf)]
o1 = [json.loads(l) for l in open('logs/o1.jsonl') if json.loads(l)['c'] == cs][-1]
r = [x for x in rows if x['c'] == cs][-1]
open('logs/int_failed.jsonl', 'a').write(json.dumps(r) + '\n')
Lo = '%.6e' % 10**((lg(r['approx']['e1']) + lg(o1['o1']))/2)
ks = (-2, 0, 1, 2, 3, 4); x = lg(Lo)
od = cert(cs, 'odd', Lo, r['bits'], r['K2'], [e10(x + k) for k in ks], 'i%so2' % cs)
r['odd'] = od; r['Lodd'] = Lo; r['approx']['o1'] = o1['o1']; r['approx']['o1_note'] = 'o1fix.py N = %d' % o1['N']
r['even_simple'] = bool(r['defl'].get('rigorous') and od.get('rigorous') and lg(r['U_up']) < lg(r['L2']) and lg(r['U_up']) < lg(Lo))
with open(lf, 'w') as f:
    for x in rows: f.write(json.dumps(r if x['c'] == cs else x) + '\n')
print(json.dumps(dict(c=cs, Lodd=Lo, rig=od.get('rigorous'), nu=od.get('nu_cert'), even_simple=r['even_simple'], rss=od.get('rss_mb'), sec=od.get('sec'))))

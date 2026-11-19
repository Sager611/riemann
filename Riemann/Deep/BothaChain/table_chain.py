# table_chain.py (bothachain): markdown table of logs/chain.jsonl (bothaaudit's 3 steps from c = 30, then ours; splice rows = CCMFMM windows)
# and re-check of every step inequality from the records: U(c_i) < eps_i, both sectors rigorous (nu > kappa eta), Lemma G check, consecutive.
import json, os, math
from fractions import Fraction as F
H = os.path.dirname(os.path.abspath(__file__))
def lg(s): s = str(s).strip('[').split(' ')[0]; m, _, e = s.partition('e'); return math.log10(float(m)) + int(e)
import sys
R = [d for d in map(json.loads, open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, 'logs', 'chain.jsonl'))) if 'defl' in d]   # arg: other chain file (tests)
print('| # | c_i → c_{i+1} | U(c_i) ≥ λ₁ᵉᵛᵉⁿ on the step | ε: λ₁ᵒᵈᵈ, λ₂ᵉᵛᵉⁿ(c_{i+1}) ≥ ε | ε/U(c_i) | certified ν (defl, odd); max κη | U(c_{i+1}) | NG; n (defl/odd) | bits | peak RSS MB (g; prep, LDL defl; prep, LDL odd; scan) | s | by |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
end = '30'; mx = 0; nok = nf = 0; k = 0; bad = []
for i, d in enumerate(R):
    rd, ro = d['defl'], d['odd']; by = 'splice: CCMFMM window' if d.get('splice') else ('bothaaudit' if i < 3 else 'bothachain')
    sc = (d.get('scan2') or {}).get('rss_mb')
    rss = [d.get('g_rss')] + (rd.get('rss_mb') or []) + (ro.get('rss_mb') or []) + [sc]
    if not d.get('splice'): mx = max([mx] + [x for x in rss if x])
    ke = max([x for x in (rd.get('kappa_eta'), ro.get('kappa_eta')) if x] or [0])
    if d['ok']:
        k += 1
        dat, dlev = d.get('defl_at') or d['c'], d.get('defl_level') or d['eps']   # pair records: defl certificate at defl_at >= c, level defl_level
        chk = (lg(d['Uprev']) < lg(d['eps']) and rd.get('rigorous') and ro.get('rigorous') and float(rd['nu_cert']) > rd['kappa_eta']
               and float(ro['nu_cert']) > ro['kappa_eta'] and d.get('lemmaG', i < 3) and F(d['c']) > F(d['cprev'])
               and ro['eps'] == d['eps'] and rd['eps'] == dlev and lg(d['Uprev']) < lg(dlev) and F(dat) >= F(d['c']))
        if not chk: bad.append(d['c'])
        if d['cprev'] == end: end = d['c']
        nok += 1
    else: nf += 1
    if d.get('pair') and d.get('defl_at') and d['defl_at'] != d['c']: by += ' (defl shared: at %s, level %s)' % (d['defl_at'], d['defl_level'])
    print('| %s | %s → %s | %s | %s | %.3f | %s, %s; %.1e | %s | %s; %s/%s | %s | %s; %s; %s; %s | %s | %s |' % (k if d['ok'] else 'fail', d['cprev'], d['c'], d['Uprev'], d['eps'],
          10**(lg(d['eps']) - lg(d['Uprev'])), rd.get('nu_cert'), ro.get('nu_cert'), ke, d['U_up'], d['NG'], rd.get('n'), ro.get('n'), d['bits'],
          d.get('g_rss'), rd.get('rss_mb'), ro.get('rss_mb'), sc, d['sec'], by))
print()
print('covered by the chain: [30, %s] = [30, %.6f]; ok steps %d (bothaaudit 3), failed attempts %d; max RSS (our runs) %d MB; re-check failures: %s' % (
      end, float(F(end)), nok, nf, mx, bad or 'none'))

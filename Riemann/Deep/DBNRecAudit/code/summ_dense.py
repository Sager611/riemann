#!/usr/bin/env python3
# dbnrecaudit summ_dense.py: markdown summary of logs/dense_r*.json (falsify_dense.py) and logs/scan_r8_*off.json (scan_wrap.py).
import json, glob, os
L = 'Riemann/Deep/DBNRecAudit/logs/'
print('| row | θ | window (points) | max Q_i | G = max Q_i − s·log(x/X) | (1 − θ)A′ + θB′ | margin |')
print('|---|---|---|---|---|---|---|')
extra = []
for f in sorted(glob.glob(L + 'dense_r*.json'), key=lambda p: int(os.path.basename(p)[7:-5])):
    r = json.load(open(f))
    for x in r['runs']:
        npts = '10⁶' if x['window'] == 'cell' else '2·10⁵'
        print('| %d | %g | %s (%s) | %.6f | %.6f | %.6f | %.2e |' % (r['row'], x['theta'], 'binding cell' if x['window'] == 'cell' else 'x′_c ± 2', npts, x['Qmax'], x['G'], x['rhs'], x['margin']))
    extra.append('| %d | %.7f | %.7f | %.1e | %s | %.1e | %.1e | %.3g / %.3g |' % (r['row'], r['Qplus_cell']['max'], r['ub_arb'], r['Qplus_cell']['slack'],
                 '%.7f' % r['Qplus_wide']['max'], r['amp_upper_worst_rel'], r['amp_lower_worst_rel'], r['D_shell0_own'], r['D_shell0_cert']))
print()
print('| row | max Q⁺ on binding cell (10⁶ points) | certify Arb UB of the cell | slack | max Q⁺ on x′_c ± 2 | max (m_n − chord)/chord | max (ℓ_n − m_n)/m_n | D own / logged |')
print('|---|---|---|---|---|---|---|---|')
for e in extra: print(e)
for f in sorted(glob.glob(L + 'scan_r8_*off.json')):
    s = json.load(open(f)); best = s['full'][0]
    print('\nscan %s: K = %d anchors, %d offsets in [%.3f, %.3f]; binding anchor rank %s; best full max Q⁺ %.6f at anchor %d (+%.3f); top-3 ranks %s' % (
        os.path.basename(f), s['K'], len(s['offsets']), s['g2_high_region'][0], s['g2_high_region'][1], s['binding_rank'], best['Qplus_max'], best['k'], best['at'], s['top10'][:3]))

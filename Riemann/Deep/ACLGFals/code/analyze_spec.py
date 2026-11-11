# analyze_spec.py: Hypothesis LA / gap / interlacing / count tests (ACLG 2611.13129 l. 588-599, 698-757, 904-920) at small lambda,
# from logs/spec.jsonl (Galerkin eigenvalues, upper bounds by min-max) and exact rungs delta_n (prolate.py), plus the certified
# bounds of weilsimple (U >= eps(1); L_odd <= eps_odd(1); L_2 <= eps(2)). Output: markdown tables to stdout, JSON to logs/la.json.
import sys, os, json, math, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prolate import rungs
here = os.path.dirname(os.path.abspath(__file__))
spec = [json.loads(l) for l in open(os.path.join(here, '..', 'logs', 'spec.jsonl'))]
# weilsimple certified table (REPORT.md l. 55-73): c: (U, L_odd, L_2)
cert = {}
for l in open(os.path.join(here, '..', '..', 'WeilSimple', 'REPORT.md')):
    p = [x.strip() for x in l.split('|')]
    if len(p) > 6 and p[1].isdigit() and p[2].replace('.', '').isdigit() and 'e' in p[3]:
        cert[int(p[1])] = tuple(float(x) for x in p[3:6])
byc = {}
for r in spec: byc.setdefault(int(r['c']), {})[r['N']] = r
out = {}
rows_la, rows_gap, rows_cnt, rows_cert = [], [], [], []
for c in sorted(byc):
    lam = math.sqrt(c); L = math.log(lam)
    R, dps, K = rungs(mp.mpf(c))
    d = {n: v['delta'] for n, v in R.items()}
    Ns = sorted(byc[c]); Nb = Ns[-1]
    ev = [mp.mpf(x) for x in byc[c][Nb]['even']]; od = [mp.mpf(x) for x in byc[c][Nb]['odd']]
    ev0 = [mp.mpf(x) for x in byc[c][Ns[0]]['even']]; od0 = [mp.mpf(x) for x in byc[c][Ns[0]]['odd']]
    mins_e = [j for j in range(len(ev)) if ev[j] < 1e-2]; mins_o = [j for j in range(len(od)) if od[j] < 1e-2]
    re = [float(ev[j]/d[4*(j + 1)]) for j in mins_e if 4*(j + 1) in d]
    ro = [float(od[j]/d[4*(j + 1) + 2]) for j in mins_o if 4*(j + 1) + 2 in d]
    lg = float(mp.log(ev[1]/ev[0])); lgd = float(mp.log(d[8]/d[4])); lgf = 8*L + 8.2425
    lo = float(mp.log(od[0]/ev[0])); lod = float(mp.log(d[6]/d[4])); lof = 4*L + 4.4342
    seq = sorted([(x, 'e') for x in ev if x < 1e-2] + [(x, 'o') for x in od if x < 1e-2])
    inter = all(seq[i][1] != seq[i + 1][1] for i in range(len(seq) - 1)) and (len(seq) == 0 or seq[0][1] == 'e')
    def cnt(th):
        qe = sum(1 for x in ev if x < th); qo = sum(1 for x in od if x < th)
        u0 = sum(1 for n, x in d.items() if n % 4 == 0 and x < th); u2 = sum(1 for n, x in d.items() if n % 4 == 2 and x < th)
        return f'{qe}/{qo} | {u0}/{u2}'
    conv = float(mp.log(ev0[0]/ev[0]))
    out[c] = {'N': Nb, 'eps1': mp.nstr(ev[0], 6), 'eps2': mp.nstr(ev[1], 6), 'odd1': mp.nstr(od[0], 6), 'D1': float(-mp.log(d[4])),
              'ratios_even': re, 'ratios_odd': ro, 'loggap': [lg, lgd, lgf], 'logodd': [lo, lod, lof], 'interlace': inter,
              'r1_log': float(ev[0]/(d[4]*L)), 'conv_logeps1_N120_to_N240': conv}
    rows_la.append(f"| {c} | {lam:.4f} | {Nb} | {float(-mp.log(d[4])):.3f} | {float(-mp.log(ev[0])):.3f} | {re[0] if re else float('nan'):.3f} | "
                   f"{float(ev[0]/(d[4]*L)):.3f} | {(min(re) if re else float('nan')):.2f}-{(max(re) if re else float('nan')):.2f} ({len(re)}) | "
                   f"{(ro[0] if ro else float('nan')):.3f} | {(min(ro) if ro else float('nan')):.2f}-{(max(ro) if ro else float('nan')):.2f} ({len(ro)}) | {conv:.4f} |")
    rows_gap.append(f"| {c} | {lg:.3f} | {lgd:.3f} | {lgd - lg:.3f} | {lgf:.3f} | {lo:.3f} | {lod:.3f} | {lod - lo:.3f} | {lof:.3f} | {'yes' if inter else 'NO'} |")
    rows_cnt.append(f"| {c} | {cnt(1e-10)} | {cnt(1e-5)} | {cnt(1e-2)} |")
    if c in cert:
        U, Lo, L2 = cert[c]
        rows_cert.append(f"| {c} | {U/float(d[4]):.3f} | {U/(float(d[4])*L):.3f} | {Lo/float(d[6]):.4f} | {L2/float(d[8]):.4f} | "
                         f"{min(Lo, L2)/float(d[6]):.4f} | {float(od[0])/float(d[6]):.3f} | {float(ev[1])/float(d[8]):.3f} |")
print('| c = lambda^2 | lambda | N | D_1 = -log delta_4 | -log eps(1) | eps(1)/delta_4 | eps(1)/(delta_4 log lambda) | eps(j)/delta_4j range (#minuscule) | eps_odd(1)/delta_6 | eps_odd(j)/delta_4j+2 range (#) | log eps1(N120)/eps1(N240) |')
print('|---|---|---|---|---|---|---|---|---|---|---|'); print('\n'.join(rows_la)); print()
print('| c | log eps(2)/eps(1) | log delta_8/delta_4 | shortfall | Fuchs 8 log lam + 8.2425 | log eps_odd(1)/eps(1) | log delta_6/delta_4 | shortfall | Fuchs 4 log lam + 4.4342 | interlacing |')
print('|---|---|---|---|---|---|---|---|---|---|'); print('\n'.join(rows_gap)); print()
print('| c | < 1e-10: QW ev/odd, U 0/2 | < 1e-5 | < 1e-2 |'); print('|---|---|---|---|'); print('\n'.join(rows_cnt)); print()
print('| c | U/delta_4 (cert. upper) | U/(delta_4 log lambda) | L_odd/delta_6 (cert. lower) | L_2/delta_8 (cert. lower) | min(L_odd, L_2)/delta_6 | Galerkin eps_odd(1)/delta_6 | Galerkin eps(2)/delta_8 |')
print('|---|---|---|---|---|---|---|---|'); print('\n'.join(rows_cert))
json.dump(out, open(os.path.join(here, '..', 'logs', 'la.json'), 'w'), indent=1, default=str)

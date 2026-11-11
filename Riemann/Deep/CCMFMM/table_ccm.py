# table_ccm.py (ccmfmm): tables for REPORT.md from logs/ (int.jsonl, chain.jsonl, val/, adv/) and a re-check of every chain inequality.
import json, os, math, glob, re, sys
from fractions import Fraction
H = os.path.dirname(os.path.abspath(__file__)); LG = os.path.join(H, 'logs')
def lg(s):
    m, _, ex = str(s).partition('e'); return math.log10(abs(float(m))) + (int(ex) if ex else 0)
def sh(s, d=3):
    x = lg(s); e = math.floor(x); return '%.*fe%d' % (d, 10**(x - e), e)
def rss(r): return max(x for x in r if x is not None)
def main():
    rows = [json.loads(l) for l in open(os.path.join(LG, 'int.jsonl'))] if os.path.exists(os.path.join(LG, 'int.jsonl')) else []
    print('| c = λ² | λ | half-width log λ | U ≥ λ₁ᵉᵛᵉⁿ | L_odd | L_2 | U < min | NG | bits | n (defl / odd) | certified ν (defl, odd) | max κη | peak RSS MB (g, prep, LDL) | s (g + certificates) |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in rows:
        c = float(Fraction(r['c'])); d, o = r['defl'], r['odd']
        ok = r['even_simple']
        print(f"| {r['c']} | {math.sqrt(c):.4f} | {math.log(c)/2:.4f} | {sh(r['U_up'])} | {sh(r['Lodd'])} | {sh(r['L2'])} | {'yes' if ok else 'NO'} | {r['NG']} | {r['bits']} | {d.get('n')} / {o.get('n')} | "
              f"{sh(d['nu_cert'],2) if d.get('nu_cert') else '-'}, {sh(o['nu_cert'],2) if o.get('nu_cert') else '-'} | {max(d.get('kappa_eta') or 0, o.get('kappa_eta') or 0):.1e} | "
              f"{r['g_rss']}, {max(d['rss_mb'][0] or 0, o['rss_mb'][0] or 0)}, {max(d['rss_mb'][1] or 0, o['rss_mb'][1] or 0)} | {r['sec']:.0f} |")
    ch = os.path.join(LG, 'chain.jsonl')
    if os.path.exists(ch):
        recs = [json.loads(l) for l in open(ch)]
        print('\n| c_i → c_{i+1} | U(c_i) | L(c_{i+1}) = ε | U(c_{i+1}) | certified ν (defl, odd) | n (defl / odd) | NG | bits | ok | s |')
        print('|---|---|---|---|---|---|---|---|---|---|')
        for r in recs:
            d, o = r['defl'], r['odd']
            print(f"| {r['cprev']} → {r['c']} | {sh(r['Uprev'],2)} | {sh(r['eps'],2)} | {sh(r['U_up'],2)} | {sh(d['nu_cert'],1) if d.get('nu_cert') else '-'}, {sh(o['nu_cert'],1) if o.get('nu_cert') else '-'} | {d.get('n')} / {o.get('n')} | {r['NG']} | {r['bits']} | {'ok' if r['ok'] else 'FAIL'} | {r['sec']:.0f} |")
        # re-check: contiguous ok steps, U(c_i) < eps(c_{i+1}), both certificates rigorous (nu > kappa eta)
        good = [r for r in recs if r['ok']]
        segs = []; bad = 0; maxr = 0
        for r in good:
            assert lg(r['Uprev']) < lg(r['eps']), r
            for s in ('defl', 'odd'):
                x = r[s]; assert x['rigorous'] and float(x['nu_cert']) > x['kappa_eta']; maxr = max(maxr, max(x['rss_mb']))
            a, b = Fraction(r['cprev']), Fraction(r['c'])
            if segs and segs[-1][1] == a: segs[-1][1] = b
            else: segs.append([a, b])
        print('\nchain check: %d ok records, %d failed attempts; covered: %s; max RSS %d MB' % (len(good), len(recs) - len(good), ', '.join('[%s, %s]' % (str(a), str(b)) for a, b in segs), maxr))
    vd = os.path.join(LG, 'val')
    if os.path.isdir(vd):
        print('\n| c | sector | ε (WeilSimple) | WeilSimple margin m | ν = m(1 − 1e-9) | ν = m(1 + 1e-9) | κη | n | RSS MB (prep, LDL) | s (prep, LDL) |')
        print('|---|---|---|---|---|---|---|---|---|---|')
        for c in (10, 16, 20):
            for s, nm in (('d', 'defl'), ('o', 'odd')):
                f = os.path.join(vd, f'v{c}{s}.out')
                if not os.path.exists(f): continue
                ls = [json.loads(l) for l in open(f) if l.startswith('{')]
                p, a, b = ls[0], ls[1], ls[2]
                m = float(a['nu'])/(1 - 1e-9)
                tp = open(os.path.join(vd, f'v{c}{s}.prep.time')).read(); tl = open(os.path.join(vd, f'v{c}{s}.slt.time')).read()
                rp = int(re.search(r'(\d+)\s+maximum resident', tp).group(1))//2**20; rl = int(re.search(r'(\d+)\s+maximum resident', tl).group(1))//2**20
                sp = float(re.search(r'([\d.]+) real', tp).group(1)); sl = float(re.search(r'([\d.]+) real', tl).group(1))
                print(f"| {c} | {nm} | {p['eps']} | {m:.11e} | {'pass' if a['inertia_ok'] else 'FAIL'} | {'fails' if not b['inertia_ok'] else 'PASSES'} | {a['kappa_eta']:.1e} | {a['n']} | {rp}, {rl} | {sp:.1f}, {sl:.1f} |")
    ad = os.path.join(LG, 'adv')
    if os.path.isdir(ad):
        print('\n| c | sector | ε | ν = −1 | ν = −0.1 | ν = 0 | ν = ε | n_U (expected #pos) | verdict |')
        print('|---|---|---|---|---|---|---|---|---|')
        for f in sorted(glob.glob(os.path.join(ad, '*.slt.json'))):
            ls = [json.loads(l) for l in open(f) if l.startswith('{')]
            if not ls: continue
            tag = os.path.basename(f).split('.')[0]; c = tag[3:-1]; s = 'defl' if tag.endswith('d') else 'odd'
            cells = ['%s (%d pos)' % ('ok' if x['inertia_ok'] else 'wrong', x['U'][0] + x['B'][0]) for x in ls]
            cert = any(x['inertia_ok'] and float(x['nu']) > x['kappa_eta'] for x in ls)
            print(f"| {c} | {s} | {ls[-1]['nu']} | " + ' | '.join(cells) + f" | {ls[0]['nU']} | {'CERTIFIED (bad)' if cert else 'not certified (correct)'} |")
def check_runs():   # every point-mode LDL line: sigma > eta (Lemma B needs S_s - eta > 0), kappa > 0 finite; certified lines have nu > kappa eta
    fs = sorted(glob.glob(os.path.join(LG, 'runs', '*.slt.out')))
    nl = 0; worst = 0.0; bad = []
    for f in fs:
        for l in open(f):
            if not l.startswith('{'): continue
            d = json.loads(l); nl += 1
            if not (d['sigma'] > 2*d['eta'] and 0 < d['kappa'] < 1e300 and d['point'] == 1): bad.append((os.path.basename(f), d['nu']))
            if d['eta'] > 0: worst = max(worst, d['eta']/d['sigma'])
    if fs: print('\nLDL lines checked: %d in %d files; sigma > 2 eta, kappa finite, point mode: %s; max eta/sigma = %.1e' % (nl, len(fs), 'all' if not bad else 'FAILS ' + str(bad[:5]), worst))
main()
check_runs()

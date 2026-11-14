# run_ccm.py (ccmfmm): CCM's first missing step (simple even ground state of QW_lambda, c = lambda^2) at a list of windows, in the
# particle form (prep_ccm.py + slt.c point mode). Per window c: gstate.py (g, U = Q(g)/||g||^2 >= lambda_1^even, approximate e1, e2, o1),
# levels L_2 = (e1 e2)^{1/2}, L_odd = (e1 o1)^{1/2} as in ../WeilSimple/run_windows.py, then the defl certificate at L_2 and the odd
# certificate at L_odd: each LDL at nu in a ladder; "pass" = inertia (nU, nB) with nu > kappa eta, i.e. lambda_min(Id - W) >= nu - kappa eta > 0.
# even_simple = U < L_2 and U < L_odd and both pass (Courant-Fischer, ../WeilSimple/DERIVATION.md sec. 3).
# Usage: run_ccm.py <log.jsonl> c1 c2 ...      (one process at a time; every child under alarm and /usr/bin/time -l)
import sys, os, json, time, subprocess, math, re
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.abspath(os.environ.get('CCM_RUNS', os.path.join(HERE, '..', '..', '..', '.scratch', 'ccmfmm', 'runs')))
BIN = os.path.abspath(os.environ.get('CCM_BIN', os.path.join(HERE, '..', '..', '..', '.scratch', 'ccmfmm')))
UV = ['uv', 'run', '--no-project', '--with', 'mpmath', '--with', 'numpy', '--with', 'python-flint', 'python3']
ENV = dict(os.environ, WZ_RSCAN=os.path.join(BIN, 'rscan'), WZ_SCAN=os.path.join(BIN, 'scan'), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')

def lg(s):
    m, _, ex = str(s).partition('e'); return math.log10(abs(float(m))) + (int(ex) if ex else 0)
def e10(x): return '%.6e' % (10**(x - math.floor(x))) if False else ('%.6fe%d' % (10**(x - math.floor(x)), math.floor(x)))

def run(cmd, out, tl=1790, env=ENV):
    t0 = time.time()
    with open(out + '.out', 'w') as fo, open(out + '.time', 'w') as fe:
        p = subprocess.run(['/usr/bin/time', '-l', 'perl', '-e', 'alarm %d; exec @ARGV' % tl] + cmd, stdout=fo, stderr=fe, env=env, cwd=HERE)
    err = open(out + '.time').read(); m = re.search(r'(\d+)\s+maximum resident set size', err)
    js = [json.loads(l) for l in open(out + '.out') if l.startswith('{')]
    return dict(rc=p.returncode, js=js, rss_mb=round(int(m.group(1))/2**20) if m else None, sec=round(time.time() - t0, 1), err=err[-300:] if p.returncode else '')

def gstate(cs, NG, bits, NO, tag):
    gf = os.path.join(RUNS, tag + '.json')
    if not os.path.exists(gf):
        r = run(UV + ['gstate.py', cs, str(NG), str(bits), gf, str(NO)], os.path.join(RUNS, tag))
        if r['rc']: raise SystemExit('gstate failed ' + r['err'])
        r['js'][0]['rss_mb'] = r['rss_mb']; r['js'][0]['sec'] = r['sec']
        json.dump(r['js'][0], open(gf + '.meta', 'w'))
    return json.load(open(gf)), json.load(open(gf + '.meta'))

def cert(cs, sector, eps, bits, K2, nus, tag, gfile=None):
    b = os.path.join(RUNS, tag); env = dict(ENV)
    if gfile: env['CCM_G'] = gfile
    rp = None
    if not os.path.exists(b + '.idx'):
        rp = run(UV + ['prep_ccm.py', cs, sector, eps, str(bits), b, str(K2)], b + '.prep', env=env)
        if rp['rc']: return dict(error='prep ' + rp['err'])
    env2 = dict(ENV, SLT_POINT='1', SLT_BRD=(b + '.brd' if sector == 'defl' else ''))
    rs = run([os.path.join(BIN, 'slt'), b, str(bits), '-g', b + '.gen'] + nus, b + '.slt', env=env2)
    if rs['rc']: return dict(error='slt ' + rs['err'])
    prep = json.load(open(b + '.prep.out')) if os.path.exists(b + '.prep.out') else {}
    lad = []
    for d in rs['js']:
        ok = d['inertia_ok'] and float(d['nu']) > d['kappa_eta']
        lad.append(dict(nu=d['nu'], ok=bool(d['inertia_ok']), pos=d['U'][0] + d['B'][0], nU=d['nU'], kappa_eta=d['kappa_eta'], certified=bool(ok), t=d['time']))
    good = [x for x in lad if x['certified']]; bad = [x for x in lad if not x['ok']]
    return dict(sector=sector, eps=eps, bits=bits, K2=K2, n=prep.get('n'), nU=prep.get('nU'), order=prep.get('order'), G=prep.get('G'),
                rigorous=bool(good), nu_cert=(max(good, key=lambda x: float(x['nu']))['nu'] if good else None),
                nu_fail=(min(bad, key=lambda x: float(x['nu']))['nu'] if bad else None), kappa_eta=max(x['kappa_eta'] for x in lad) if lad else None,
                ladder=lad, rss_mb=[rp['rss_mb'] if rp else None, rs['rss_mb']], sec=[rp['sec'] if rp else None, rs['sec']], tscan=prep.get('tscan'))

def ngrule(cf):
    return 160 if cf <= 20 else (200 if cf <= 25 else (240 if cf <= 30 else (320 if cf <= 35 else 400)))
def bitsrule(le1):
    return max(192, 64*math.ceil((1.6*(-le1)*math.log2(10) + 64)/64))

def main():
    log = sys.argv[1]
    for cs in sys.argv[2:]:
        t0 = time.time(); cf = eval(cs.replace('/', '/1.0/')) if '/' in cs else float(cs)
        K2 = max(401, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
        NG = int(os.environ.get('CCM_NG', ngrule(cf)))
        guess = -5.35*cf + 11.0                              # log10 lambda_1^even (fit to ../WeilSimple), only for the precision
        bits = bitsrule(guess)
        tagc = cs.replace('/', 'o')
        if cf <= 30:
            g, gm = gstate(cs, NG, bits, NG, 'g' + tagc)   # odd block at NO = NG (NO = 120 was unconverged for c >= 22)
            e1, e2, o1 = gm['even'][0], gm['even'][1], gm['odd'][0]
        else:   # no full eigendecompositions (cost): e1 := U, o1 by inverse iteration on the odd block at N = NG (o1fix.py), L_2 := L_odd
            g, gm = gstate(cs, NG, bits, 0, 'g' + tagc)
            fo = os.path.join(RUNS, 'o1_' + tagc + '.jsonl')
            if not os.path.exists(fo):
                ro = run(UV + ['o1fix.py', fo, str(NG), str(bits), cs], os.path.join(RUNS, 'o1_' + tagc))
                if ro['rc']: raise SystemExit('o1fix failed ' + ro['err'])
                gm['o1_rss'] = ro['rss_mb']; gm['o1_sec'] = ro['sec']; gm['rss_mb'] = max(gm['rss_mb'], ro['rss_mb'] or 0)
            e1 = gm['U_up']; o1 = json.loads(open(fo).readline())['o1']; e2 = o1
        bits = bitsrule(lg(e1)); U = gm['U_up']
        L2 = '%.6e' % 10**((lg(e1) + lg(e2))/2) if (lg(e1) + lg(e2))/2 > -300 else e10((lg(e1) + lg(e2))/2)
        Lo = '%.6e' % 10**((lg(e1) + lg(o1))/2) if (lg(e1) + lg(o1))/2 > -300 else e10((lg(e1) + lg(o1))/2)
        row = dict(c=cs, NG=NG, bits=bits, K2=K2, U_up=U, U=gm['U'], approx=dict(e1=e1, e2=e2, o1=o1), g_sec=gm['sec'], g_rss=gm['rss_mb'], L2=L2, Lodd=Lo)
        ks = (-2, 0, 1, 2, 3, 4) if cf <= 25 else ((-2, 0, 2, 4) if cf <= 30 else (-1, 1))   # each nu is one O(n^2) LDL
        def ladder(eps): x = lg(eps); return [e10(x + k) for k in ks]
        row['defl'] = cert(cs, 'defl', L2, bits, K2, ladder(L2), 'i%sd' % tagc, os.path.join(RUNS, 'g' + tagc + '.json'))
        row['odd'] = cert(cs, 'odd', Lo, bits, K2, ladder(Lo), 'i%so' % tagc)
        row['even_simple'] = bool(row['defl'].get('rigorous') and row['odd'].get('rigorous') and lg(U) < lg(L2) and lg(U) < lg(Lo))
        row['sec'] = round(time.time() - t0, 1)
        with open(log, 'a') as f: f.write(json.dumps(row) + '\n')
        print(json.dumps(dict(c=cs, NG=NG, bits=bits, U=U, L2=L2, Lodd=Lo, rig=[row['defl'].get('rigorous'), row['odd'].get('rigorous')],
                              nu=[row['defl'].get('nu_cert'), row['odd'].get('nu_cert')], even_simple=row['even_simple'], sec=row['sec'],
                              rss=[row['defl'].get('rss_mb'), row['odd'].get('rss_mb'), gm['rss_mb']])), flush=True)

if __name__ == '__main__':
    main()

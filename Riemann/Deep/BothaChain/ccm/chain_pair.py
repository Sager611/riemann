# chain_pair.py (bothachain): the chain of chain_ccm.py with one deflated-even certificate per PAIR of steps (same certificates, same rules).
# Pair c_a < c_b < c_c, U_a = U(c_a) given; U_b, U_c from gstate (Lemma G); eps_a = 2.00001 U_a, eps_b = 2.00001 U_b (chain_ccm's rounding).
#   odd at c_b, level eps_a;  odd at c_c, level eps_b;  defl at c_c, level eps_a (g = gstate's vector at c_c).
# Monotonicity (every lambda_k non-increasing in c) gives, for c in [c_a, c_b]: lambda_1^even(c) <= U_a < eps_a <= lambda_1^odd(c_b) <= lambda_1^odd(c),
# lambda_2^even(c) >= lambda_2^even(c_c) >= eps_a; for c in [c_b, c_c]: lambda_1^even(c) <= U_b < eps_b <= lambda_1^odd(c_c) <= lambda_1^odd(c),
# lambda_2^even(c) >= lambda_2^even(c_c) >= eps_a > eps_b.  Each record carries its odd certificate and the defl certificate it uses (defl_at, defl
# level); table_chain.py re-checks: odd level = eps at c, defl level > U(c_i) at a point >= c.  Shared defl only if the e2 heuristic (0.6 x scan U2)
# at c_c exceeds 7.1 eps_a; if it is not tried or fails: standard defl at c_b (level eps_a) and at c_c (level eps_b).
# usage: chain_pair.py <int.jsonl> <chain.jsonl> <c1>    (resumes from the last ok record; CCM_MAXSTEPS = pairs per call)
import sys, os, json, math, time, re
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from run_ccm import gstate, cert, lg, e10, ngrule, bitsrule, RUNS, run, UV, ENV

rows = {}
STEP = Fraction(1, 24)

def lo(c, key):
    a = math.floor(c); b = a + 1; t = float(c - a)
    if b not in rows: b = a
    return (1 - t)*lg(rows[a]['approx'][key]) + t*lg(rows[b]['approx'][key])

def cstr(cn): return str(cn.numerator) if cn.denominator == 1 else '%d/%d' % (cn.numerator, cn.denominator)
def level(U): return e10(lg(2.00001*float(U)) if float(U) > 1e-300 else lg(U) + math.log10(2.00001))

def choose(c, U, c1, failed):   # chain_ccm's greedy step, skipping (c, c_next) pairs that failed before
    target = lg(U) + math.log10(4) + 0.25; h = STEP
    while c + h + STEP <= c1 and lo(c + h + STEP, 'o1') >= target: h += STEP
    while (str(c), cstr(min(c1, c + h))) in failed:
        if h == STEP: sys.exit('stuck at ' + str(c))
        h = max(STEP, h/2)
    return min(c1, c + h)

def point(cn):   # scan premise (rscan2) + gstate at cn; returns parameters, U and the Lemma G check
    n = math.ceil(cn); cf = float(cn); cs = cstr(cn)
    NG = ngrule(n); bits = rows[n]['bits']
    K2 = max(401, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
    scan2 = None
    sf = os.path.join(ENV['WZ_SCAN'], 'scan_zeta_%s_%d_rig.json' % (cs.replace('/', 'o'), K2))
    if not os.path.exists(sf):
        rs2 = run(UV + ['scan_zeta.py', cs, str(K2)], os.path.join(RUNS, 'scan' + cs.replace('/', 'o')), env=dict(ENV, PYTHONPATH=os.path.join(HERE, '..', '..', 'WeilZeta')))
        if rs2['rc'] or not os.path.exists(sf): sys.exit('scan_zeta failed ' + rs2['err'])
        scan2 = dict(rss_mb=rs2['rss_mb'], sec=rs2['sec'], **{k: rs2['js'][-1].get(k) for k in ('nodes', 'E', 'selected', 'kscan')})
    g, gm = gstate(cs, NG, bits, 0, 'g' + cs.replace('/', 'o'))
    mU = re.fullmatch(r'\[([0-9.]+e-?[0-9]+) \+/- ([0-9.]+e-?[0-9]+)\]', gm['U'])   # Lemma G: finite ball, mid > 0, radius < 1e-3 mid
    return dict(cs=cs, cf=cf, NG=NG, bits=bits, K2=K2, scan2=scan2, gm=gm, lemG=bool(mU) and lg(mU.group(2)) < lg(mU.group(1)) - 3,
                gfile=os.path.join(RUNS, 'g' + cs.replace('/', 'o') + '.json'))

def lad(eps): return [e10(lg(eps) - 1)]   # one rung nu = eps/10 (bothaaudit's rule beyond c = 30)

def odd(p, eps, tag): return cert(p['cs'], 'odd', eps, p['bits'], p['K2'], lad(eps), tag)
def defl(p, eps, tag): return cert(p['cs'], 'defl', eps, p['bits'], p['K2'], lad(eps), tag, p['gfile'])

def rec(p, cprev, Uprev, eps, ro, rd, dat, dlev, t0):
    ok = bool(ro.get('rigorous') and rd and rd.get('rigorous') and p['lemG'] and lg(p['gm']['U_up']) < lg(eps))
    return dict(c=p['cs'], cprev=cprev, Uprev=Uprev, eps=eps, NG=p['NG'], bits=p['bits'], K2=p['K2'], ok=ok, U_up=p['gm']['U_up'], U=p['gm']['U'],
                g_sec=p['gm']['sec'], g_rss=p['gm']['rss_mb'], lemmaG=p['lemG'], scan2=p['scan2'], defl=rd or {}, odd=ro, defl_at=dat, defl_level=dlev,
                pair=True, sec=round(time.time() - t0, 1))

def main():
    for l in open(sys.argv[1]):
        d = json.loads(l)
        if '/' not in d['c']: rows[int(d['c'])] = d
    out = sys.argv[2]; c1 = Fraction(sys.argv[3])
    R = [json.loads(l) for l in open(out)]
    failed = {(d['cprev'], d['c']) for d in R if not d.get('ok')}
    ok = [d for d in R if d.get('ok')]; c = Fraction(ok[-1]['c']); U = ok[-1]['U_up']
    maxp = int(os.environ.get('CCM_MAXSTEPS', '0')); npair = 0
    def emit(r):
        with open(out, 'a') as f: f.write(json.dumps(r) + '\n')
        print(json.dumps(dict(c=r['c'], cprev=r['cprev'], eps=r['eps'], U=r['U_up'], ok=r['ok'], defl_at=r['defl_at'], defl_level=r['defl_level'],
                              nu=[r['defl'].get('nu_cert'), r['odd'].get('nu_cert')], n=[r['defl'].get('n'), r['odd'].get('n')], sec=r['sec'])), flush=True)
    while c < c1 and not (maxp and npair >= maxp):
        t0 = time.time()
        cb = choose(c, U, c1, failed); pb = point(cb); ea = level(U)
        ra = odd(pb, ea, 'ch' + pb['cs'].replace('/', 'o') + 'o')
        if not (ra.get('rigorous') and pb['lemG'] and lg(pb['gm']['U_up']) < lg(ea)):
            r = rec(pb, str(c), U, ea, ra, None, None, None, t0); emit(r); failed.add((str(c), pb['cs'])); continue
        Ub = pb['gm']['U_up']; eb = level(Ub)
        if cb >= c1:   # last single step: standard defl at c_b
            rd = defl(pb, ea, 'ch' + pb['cs'].replace('/', 'o') + 'd'); emit(rec(pb, str(c), U, ea, ra, rd, pb['cs'], ea, t0))
            if not rd.get('rigorous'): sys.exit('defl failed at ' + pb['cs'])
            c, U = cb, Ub; npair += 1; continue
        t1 = time.time()
        cc = choose(cb, Ub, c1, failed); pc = point(cc)
        rb = odd(pc, eb, 'ch' + pc['cs'].replace('/', 'o') + 'o')
        if not (rb.get('rigorous') and pc['lemG'] and lg(pc['gm']['U_up']) < lg(eb)):   # second step failed: close the first one alone
            emit(rec(pc, str(cb), Ub, eb, rb, None, None, None, t1)); failed.add((str(cb), pc['cs']))
            rd = defl(pb, ea, 'ch' + pb['cs'].replace('/', 'o') + 'd'); emit(rec(pb, str(c), U, ea, ra, rd, pb['cs'], ea, t0))
            if not rd.get('rigorous'): sys.exit('defl failed at ' + pb['cs'])
            c, U = cb, Ub; npair += 1; continue
        shared = None
        if lo(cc, 'e2') >= lg(ea) + math.log10(4) + 0.25:
            shared = defl(pc, ea, 'ch' + pc['cs'].replace('/', 'o') + 'dA')
        if shared and shared.get('rigorous'):
            emit(rec(pb, str(c), U, ea, ra, shared, pc['cs'], ea, t0)); emit(rec(pc, str(cb), Ub, eb, rb, shared, pc['cs'], ea, t1))
        else:   # fallback: standard defl at both points
            rda = defl(pb, ea, 'ch' + pb['cs'].replace('/', 'o') + 'd'); rdb = defl(pc, eb, 'ch' + pc['cs'].replace('/', 'o') + 'd')
            emit(rec(pb, str(c), U, ea, ra, rda, pb['cs'], ea, t0)); emit(rec(pc, str(cb), Ub, eb, rb, rdb, pc['cs'], eb, t1))
            if not (rda.get('rigorous') and rdb.get('rigorous')): sys.exit('defl failed in pair at ' + pb['cs'])
        c, U = cc, pc['gm']['U_up']; npair += 1
    print('CHAIN DONE ->', str(c), flush=True)

if __name__ == '__main__':
    main()

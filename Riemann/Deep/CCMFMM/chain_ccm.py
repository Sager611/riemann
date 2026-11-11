# chain_ccm.py (ccmfmm): certified chain c_0 < c_1 < ... covering a continuum of windows (../WeilSimple/chain.py, particle form).
# Every lambda_k of each parity is non-increasing in c (L^2(I_c) sits in L^2(I_c'), same form; ../WeilSimple/DERIVATION.md sec. 5), so
#   lambda_1^even(c) <= U(c_i) < L(c_{i+1}) <= lambda_1^odd(c), lambda_2^even(c)   for every c in [c_i, c_{i+1}]
# once both certificates hold at c_{i+1} with eps = L(c_{i+1}) = 2.00001 U(c_i); U(c_{i+1}) = Q(g)/||g||^2 from gstate.py at c_{i+1}.
# Step sizes from the approximate odd eigenvalues of the integer rows (interpolated in log10), as in chain.py.
# Usage: chain_ccm.py <int.jsonl> <chain.jsonl> <c0> <c1>     (c0 needs U: a row of int.jsonl or CCM_U0=<U(c0)> with CCM_O0=<o1(c0)>)
import sys, os, json, math, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from run_ccm import gstate, cert, lg, e10, ngrule, bitsrule, RUNS

def main():
    rows = {}
    for l in open(sys.argv[1]):
        d = json.loads(l)
        if '/' not in d['c']: rows[int(d['c'])] = d
    out = sys.argv[2]; c0 = Fraction(sys.argv[3]); c1 = Fraction(sys.argv[4])
    fx = os.path.join(os.path.dirname(sys.argv[1]), 'o1.jsonl')   # converged o1 (o1fix.py) override the int rows
    if os.path.exists(fx):
        for l in open(fx):
            d = json.loads(l)
            if '/' not in d['c'] and int(d['c']) in rows: rows[int(d['c'])]['approx']['o1'] = d['o1']
    if os.environ.get('CCM_O0'):
        rows[int(c0)] = dict(approx=dict(o1=os.environ['CCM_O0']), U_up=os.environ['CCM_U0'])
    def lo1(c):
        a = math.floor(c); b = a + 1; t = float(c - a)
        if b not in rows: b = a
        return (1 - t)*lg(rows[a]['approx']['o1']) + t*lg(rows[b]['approx']['o1'])
    c = c0; U = rows[int(c0)]['U_up']
    if os.path.exists(out):
        done = [json.loads(l) for l in open(out)]; done = [d for d in done if d.get('ok')]
        if done: c = Fraction(done[-1]['c']); U = done[-1]['U_up']
    step = Fraction(1, 24)
    while c < c1:
        target = lg(U) + math.log10(4) + 0.25
        h = step
        while c + h + step <= c1 and lo1(c + h + step) >= target: h += step
        while True:
            cn = min(c1, c + h); n = math.ceil(cn); cf = float(cn)
            NG = ngrule(n); bits = rows[n]['bits'] if n in rows else bitsrule(lg(U) - 5.4*float(h))
            K2 = max(401, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
            eps = e10(lg(2.00001*float(U)) if float(U) > 1e-300 else lg(U) + math.log10(2.00001))
            cs = str(cn.numerator) if cn.denominator == 1 else '%d/%d' % (cn.numerator, cn.denominator)
            tag = 'ch' + cs.replace('/', 'o'); t0 = time.time()
            g, gm = gstate(cs, NG, bits, 0, 'g' + cs.replace('/', 'o'))
            x = lg(eps); lad = [e10(x + k) for k in ((-2, 0, 2) if cf <= 25 else (-1, 1))]   # fewer O(n^2) LDLs beyond 25 (time)
            rd = cert(cs, 'defl', eps, bits, K2, lad, tag + 'd', os.path.join(RUNS, 'g' + cs.replace('/', 'o') + '.json'))
            ro = cert(cs, 'odd', eps, bits, K2, lad, tag + 'o') if rd.get('rigorous') else dict(rigorous=None)
            ok = bool(rd.get('rigorous') and ro.get('rigorous'))
            rec = dict(c=cs, cprev=str(c), Uprev=U, eps=eps, NG=NG, bits=bits, K2=K2, ok=ok, U_up=gm['U_up'], U=gm['U'], g_sec=gm['sec'], g_rss=gm['rss_mb'],
                       defl=rd, odd=ro, sec=round(time.time() - t0, 1))
            with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
            print(json.dumps(dict(c=cs, cprev=str(c), eps=eps, U=rec['U_up'], ok=ok, nu=[rd.get('nu_cert'), ro.get('nu_cert')], n=[rd.get('n'), ro.get('n')],
                                  rss=[rd.get('rss_mb'), ro.get('rss_mb'), gm['rss_mb']], sec=rec['sec'])), flush=True)
            if ok and lg(rec['U_up']) < lg(eps): break
            if h == step: sys.exit('stuck at ' + str(c))
            h = max(step, h/2)
        c = cn; U = rec['U_up']
    print('CHAIN DONE', str(c0), '->', str(c), flush=True)

if __name__ == '__main__':
    main()

# chain.py: certified chain c_0 = 2 < c_1 < ... < c_m = 16 covering every window c in [2, 16] (lambda in [sqrt 2, 4]).
# For consecutive grid points: lambda_1^even(c) <= U(c_i) and lambda_1^odd(c) >= L(c_{i+1}), lambda_2^even(c) >= L(c_{i+1})
# for all c in [c_i, c_{i+1}] (eigenvalues of each parity sector are non-increasing in c: min-max on L^2(I_c) in L^2(I_c')).
# At c_{i+1} both certificates of cert_simple.py are run at eps = L(c_{i+1}) = 2 U(c_i) (rounded up), and the defl run also
# gives U(c_{i+1}). Step sizes come from interpolating the approximate spectra of the integer windows (logs/int.jsonl).
# Usage: chain.py logs/int.jsonl logs/chain.jsonl [cstart cend]
import sys, os, json, math, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from run_windows import run_cert, lg

def main():
    rows = {int(json.loads(l)['c']): json.loads(l) for l in open(sys.argv[1])}
    out = sys.argv[2]
    c0 = Fraction(sys.argv[3]) if len(sys.argv) > 3 else Fraction(2); c1 = Fraction(sys.argv[4]) if len(sys.argv) > 4 else Fraction(16)
    def lo1(c):   # interpolated log10 of the approximate smallest odd eigenvalue (from the integer rows)
        a = math.floor(c); b = min(a + 1, max(rows)); t = float(c - a)
        f = lambda n: 2*lg(rows[n]['eps_odd']) - lg(rows[n]['approx']['e1'])   # eps_odd = (e1 o1)^{1/2}
        return (1 - t)*f(a) + t*f(b)
    c = c0; U = rows[int(c0)]['U_up'] if c0.denominator == 1 else None
    assert U is not None
    done = []
    if os.path.exists(out):
        for l in open(out):
            d = json.loads(l)
            if d.get('ok'): done.append(d)
        if done: c = Fraction(done[-1]['c']); U = done[-1]['U_up']
    while c < c1:
        target = lg(U) + math.log10(4) + 0.25
        h = Fraction(1, 24)
        while c + h + Fraction(1, 24) <= c1 and lo1(c + h + Fraction(1, 24)) >= target: h += Fraction(1, 24)
        while True:
            cn = min(c1, c + h); n = math.ceil(cn); NG = rows[n]['NG']; bits = rows[n]['bits']; K2 = max(401, math.ceil(5*float(cn)**2/3))
            K2 += (K2 % 2 == 0)
            eps = '%.6e' % (2.00001*float(U))   # > 2 U(c_i)
            env = dict(os.environ, WS_NG=str(NG), WS_M='1', OMP_NUM_THREADS='1')
            cs = str(cn.numerator) if cn.denominator == 1 else '%d/%d' % (cn.numerator, cn.denominator)
            t0 = time.time()
            rd = run_cert([cs, 'defl', eps, str(bits), str(K2)], env, 3000)
            ro = run_cert([cs, 'odd', eps, str(bits), str(K2)], env, 3000) if rd.get('rigorous') else dict(rigorous=None)
            ok = bool(rd.get('rigorous') and ro.get('rigorous'))
            rec = dict(c=cs, cprev=str(c), Uprev=U, eps=eps, NG=NG, bits=bits, K2=K2, ok=ok, U_up=rd.get('U_up'), U=rd.get('U'),
                       defl=dict((k, rd.get(k)) for k in ('margin', 'rigorous', 'order', 'nU', 'resid', 'wall', 'rss_mb')),
                       odd=dict((k, ro.get(k)) for k in ('margin', 'rigorous', 'order', 'nU', 'resid', 'wall', 'rss_mb')), sec=round(time.time() - t0, 1))
            with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
            print(json.dumps(dict(c=cs, cprev=str(c), eps=eps, U=rec['U_up'], ok=ok, margins=[rec['defl']['margin'], rec['odd']['margin']], sec=rec['sec'])), flush=True)
            if ok and lg(rec['U_up']) < lg(eps): break
            if h == Fraction(1, 24): sys.exit('stuck at ' + str(c))
            h = max(Fraction(1, 24), h/2)
        c = cn; U = rec['U_up']
    print('CHAIN DONE', str(c0), '->', str(c), flush=True)

if __name__ == '__main__':
    main()

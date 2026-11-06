# run_windows.py: for each window c = lambda^2, choose the truncation NG and the levels eps_odd, eps_2 from approximate
# eigenvalues (eigs of the closed-form blocks), then run the two certificates of cert_simple.py, one process at a time:
#   defl (even f, deflated by g): prints U = Q(g)/||g||^2 (rigorous) and proves lambda_2^even >= eps_2;
#   odd: proves Q(f) >= eps_odd ||f||^2 for all odd f.
# Usage: run_windows.py <log.jsonl> c1 c2 ...   (c integer or p/q). Env as for cert_simple.py (WZ_RSCAN, WZ_SCAN).
import sys, os, json, time, subprocess, math, re
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from weil_matrix import WeilMat

def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    return fmpq(int(s))

def approx_eigs(cq, N, prec):
    ctx.prec = prec; W = WeilMat(cq); out = {}
    for sector in ('even', 'odd'):
        T, idx = W.block(N, sector)
        ev = sorted([e.real.mid() for e in acb_mat(T).eig(algorithm="approx")], key=lambda z: z.mid())
        out[sector] = [float(e.mid().str(20, radius=False)) if abs(e) > arb(2)**-1000 else 0.0 for e in ev[:3]]
        out[sector] = [mpf_str(e) for e in ev[:3]]
    return out

def mpf_str(e): return e.mid().str(15, radius=False)

def lg(s):   # log10 of a decimal string (handles tiny values)
    m, _, ex = s.partition('e'); return math.log10(abs(float(m))) + (int(ex) if ex else 0)

def run_cert(args, env, tl):
    t0 = time.time()
    p = subprocess.run(['perl', '-e', 'alarm %d; exec @ARGV' % tl, '/usr/bin/time', '-l', sys.executable, os.path.join(HERE, 'cert_simple.py')] + args,
                       capture_output=True, text=True, env=env)
    js = [l for l in p.stdout.splitlines() if l.startswith('{')]
    m = re.search(r'(\d+)\s+maximum resident set size', p.stderr)
    res = json.loads(js[-1]) if js else dict(error=p.stderr[-400:])
    res['wall'] = round(time.time() - t0, 1); res['rss_mb'] = round(int(m.group(1))/2**20) if m else None
    return res

def main():
    log = sys.argv[1]
    for cs in sys.argv[2:]:
        cq = parse_c(cs); cf = int(cq.p)/int(cq.q); t0 = time.time()
        # approximate spectrum; NG = first N with lambda_1^even(N) <= lambda_1^odd(N)/10
        # (and lambda_1^even, lambda_1^odd changed by less than a factor 1.5 from the previous N: converged enough)
        NG = None; prev = None
        for N in (20, 30, 40, 60, 80, 100, 120, 160, 200, 240):
            prec = 64*math.ceil((1.7*3.33*max(5.3*cf - 10, 10) + 64)/64)   # ~1.7 log2(1/lambda_1) + 64 bits
            ev = approx_eigs(cq, N, prec)
            e1, e2, o1 = ev['even'][0], ev['even'][1], ev['odd'][0]
            conv = prev is not None and abs(lg(e1) - lg(prev[0])) < 0.18 and abs(lg(o1) - lg(prev[1])) < 0.18
            prev = (e1, o1)
            if conv and lg(e1) <= lg(o1) - 1: NG = N; break
        if NG is None: print(json.dumps(dict(c=cs, error='no NG')), flush=True); continue
        bits = max(192, 64*math.ceil((1.6*(-lg(e1))*math.log2(10) + 64)/64))
        le1, le2, lo1 = lg(e1), lg(e2), lg(o1)
        eps_odd = '%.6e' % 10**((le1 + lo1)/2) if (le1 + lo1)/2 > -300 else None
        eps_2 = '%.6e' % 10**((le1 + le2)/2) if (le1 + le2)/2 > -300 else None
        env = dict(os.environ, WS_NG=str(NG), WS_M='1', OMP_NUM_THREADS='1')
        row = dict(c=cs, NG=NG, bits=bits, approx=dict(e1=e1, e2=e2, o1=o1), eps_odd=eps_odd, eps_2=eps_2)
        K2 = max(401, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)   # K2 = 201 (cert_zeta default) fails for odd f at c = 7, 9
        for kind, eps in (('defl', eps_2), ('odd', eps_odd)):
            tries = []
            for (e, k2) in ((eps, K2), (eps, 2*K2 + 1), ('%.6e' % (float(eps)/10), 2*K2 + 1)):
                r = run_cert([cs, kind, e, str(bits), str(k2)], env, 3000)
                tries.append(dict(eps=e, K2=k2, margin=r.get('margin'), rigorous=r.get('rigorous')))
                if r.get('rigorous'): break
            r['tries'] = tries; row[kind] = r
        U = row['defl'].get('U_up'); row['U_up'] = U
        ok = (U is not None and row['defl'].get('rigorous') and row['odd'].get('rigorous')
              and lg(U) < lg(row['defl']['eps']) and lg(U) < lg(row['odd']['eps']))
        row['even_simple'] = bool(ok); row['seconds'] = round(time.time() - t0, 1)
        with open(log, 'a') as f: f.write(json.dumps(row) + '\n')
        print(json.dumps(dict(c=cs, NG=NG, bits=bits, U=U, L2=row['defl'].get('eps'), Lodd=row['odd'].get('eps'), rig=[row['defl'].get('rigorous'), row['odd'].get('rigorous')],
                              margins=[row['defl'].get('margin'), row['odd'].get('margin')], even_simple=row['even_simple'], sec=row['seconds'],
                              rss=[row['defl'].get('rss_mb'), row['odd'].get('rss_mb')])), flush=True)

if __name__ == '__main__':
    main()

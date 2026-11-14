# point.py (bothaaudit): one window c: generators in Arb (gen.py), deig.c (float LDL + block inverse iteration) on the mid-points of
# both blocks, then RIGOROUS Galerkin upper bounds from deig's vectors with the Arb generator balls (min-max on E_N):
#   U   = y1^T T y1/y1^T y1 >= lambda_1^even,  U2 = lambda_max of the 2x2 Rayleigh-Ritz pencil on span(y1, y2) >= lambda_2^even,
#   Uo  = z1^T T z1/z1^T z1 >= lambda_1^odd.   (T = the exact Weil matrix of the orthonormal basis; ../CCMFMM/weil_matrix.py closed forms.)
# Floats e1, e2, e3, o1, o2 from deig's Rayleigh-Ritz matrices (orientation). Run the whole script under memcap.sh (it watches children). Usage: point.py out.jsonl c N prec [L_botha]
import sys, os, json, time, subprocess, re
from flint import arb, ctx, fmpz
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from scan import FastWeilMat, parse_c
from gen import gens, hx
X = os.path.abspath(os.path.join(HERE, '..', '..', '..', '.scratch', 'bothaaudit'))
MC = ['bash', os.path.join(HERE, '..', '..', '..', '.scratch', 'memcap.sh'), '290', '1700']

def readvec(fn):
    L = open(fn).read().split('\n'); n, k = map(int, L[0].split()); V = [[None]*n for _ in range(k)]
    for i in range(n):
        w = L[1 + i].split()
        for c in range(k): V[c][i] = arb(fmpz(int(w[2*c], 16))) * arb(2)**int(w[2*c + 1])
    return V

def quad(rows, ps, Vs):   # rigorous H_ab = V_a^T T V_b for exact-dyadic vectors V_a, T from generator balls (one pass over pairs)
    n = len(rows); k = len(Vs); d = [r[0] for r in rows]; p = [r[1] for r in rows]; a = [r[2] for r in rows]; u = [r[3] for r in rows]; x = [r[4] for r in rows]
    pv = [sum((p[i]*V[i] for i in range(n)), arb(0)) for V in Vs]; hv = [[a[i]*V[i] for i in range(n)] for V in Vs]
    H = [[ps*pv[s]*pv[t] + sum(((d[i] - ps*p[i]*p[i])*Vs[s][i]*Vs[t][i] for i in range(n)), arb(0)) for t in range(k)] for s in range(k)]
    prs = [(s, t) for s in range(k) for t in range(s, k)]
    for i in range(n):
        acc = {st: arb(0) for st in prs}; ui = u[i]; xi = x[i]
        for j in range(i + 1, n):
            K = (u[j] - ui)/(xi - x[j])
            for (s, t) in prs: acc[(s, t)] += (hv[s][j]*hv[t][i] + hv[t][j]*hv[s][i])*K
        for (s, t) in prs: H[s][t] += acc[(s, t)]
    for s in range(k):
        for t in range(s): H[s][t] = H[t][s]
    return H

def main():
    out, cs, N, prec = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]); Lb = sys.argv[5] if len(sys.argv) > 5 else ''
    ctx.prec = prec; t0 = time.time(); W = FastWeilMat(parse_c(cs)); tag = os.path.join(X, 'p_' + cs.replace('/', 'o') + '_' + str(N))
    rec = dict(c=cs, N=N, prec=prec, L_botha=Lb); G = {}
    for sec, k in (('even', 3), ('odd', 2)):
        rows, ps = gens(W, N, sec); G[sec] = (rows, ps)
        with open(tag + '.' + sec, 'w') as f:
            f.write('%d %d\n' % (len(rows), ps))
            for r in rows: f.write(' '.join(hx(v) for v in r) + '\n')
        t1 = time.time()
        pr = subprocess.run(['/usr/bin/time', '-l', os.path.join(X, 'deig2'), tag + '.' + sec, str(prec), str(k), '4', tag + '.' + sec + '.vec'],
                            capture_output=True, text=True)
        if pr.returncode: raise SystemExit('deig failed ' + pr.stderr[-500:])
        m = re.search(r'(\d+)\s+maximum resident set size', pr.stderr); pk = re.search(r'memcap: peak (\d+)', pr.stderr)
        rec[sec + '_deig'] = dict(out=pr.stdout.strip().split('\n')[0], rss_mb=round(int(m.group(1))/2**20) if m else None,
                                  memcap_peak=int(pk.group(1)) if pk else None, sec=round(time.time() - t1, 1))
        open(tag + '.' + sec + '.out', 'w').write(pr.stdout)
        V = readvec(tag + '.' + sec + '.vec'); t2 = time.time()
        if sec == 'even':
            y1, y2 = V[0], V[1]
            H = quad(rows, ps, [y1, y2]); h11, h12, h22 = H[0][0], H[0][1], H[1][1]
            g11 = sum((z*z for z in y1), arb(0)); g12 = sum((y1[i]*y2[i] for i in range(len(y1))), arb(0)); g22 = sum((z*z for z in y2), arb(0))
            U = h11/g11
            # largest root of det(H - t G) = 0: (g11 g22 - g12^2) t^2 - (h11 g22 + h22 g11 - 2 h12 g12) t + (h11 h22 - h12^2) = 0
            A2 = g11*g22 - g12*g12; B2 = h11*g22 + h22*g11 - 2*h12*g12; C2 = h11*h22 - h12*h12
            U2 = (B2 + (B2*B2 - 4*A2*C2).sqrt())/(2*A2)
            rec['U'] = U.str(12, radius=True); rec['U_up'] = U.upper().str(12, radius=False)
            rec['U2'] = U2.str(12, radius=True); rec['U2_up'] = U2.upper().str(12, radius=False)
        else:
            z1 = V[0]; Uo = quad(rows, ps, [z1])[0][0]/sum((q*q for q in z1), arb(0))
            rec['Uo'] = Uo.str(12, radius=True); rec['Uo_up'] = Uo.upper().str(12, radius=False)
        rec[sec + '_rq_sec'] = round(time.time() - t2, 1)
        for fn in (tag + '.' + sec, tag + '.' + sec + '.vec'): os.remove(fn)
    rec['sec'] = round(time.time() - t0, 1)
    open(out, 'a').write(json.dumps(rec) + '\n'); print(json.dumps(rec), flush=True)

if __name__ == '__main__':
    main()

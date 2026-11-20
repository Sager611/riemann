# hdx_rat.py (explore64): rational-cycle ceiling for links almost as long as the window (no integer in (c', c)).
# A path x -> x + p_1 -> ... -> x + p_k, steps +-ln q (q a prime power < c, an atom), endpoint p_k = ln r with c' < r < c, all positions
# distinct and inside the window (span < ln c): no link of diameter <= ln c' holds x and x + ln r, so sum_i sqrt(1 - M) >= 1 along the path
# (triangle inequality for v_x in L^2(nu)); same harmonic bound as hdx_path.py.  Rigorous value at rational alpha (arb).
import sys, math, itertools, numpy as np
from flint import arb, fmpq, ctx
from scipy.optimize import minimize
ctx.prec = 128
sys.path.insert(0, 'Riemann/Deep/ExploreHDX')
from hdx_path import bound_float, ppow

def search(c, cl, kmax=4, top=4):
    lc, lcl = math.log(c), math.log(cl)
    at = [q for q in range(2, int(math.ceil(c))) if q < c and ppow(q)]
    W = {q: math.log(ppow(q)[0])/math.sqrt(q) for q in at}
    cands = []
    for k in range(2, kmax + 1):
        for combo in itertools.combinations_with_replacement([(q, s) for q in at for s in (1, -1)], k):
            L = sum(s*math.log(q) for q, s in combo)
            if not (lcl < L < lc): continue
            if any((q, 1) in combo and (q, -1) in combo for q, _ in combo): continue
            # heuristic value with equal alpha: harmonic mean of steps
            hv = (1/sum(1/(2*W[q]) for q, _ in combo))/(k + 1)
            cands.append((hv, combo))
    cands.sort(reverse=True); out = []
    for hv, combo in cands[:400]:
        best = None
        for perm in set(itertools.permutations(combo)):
            ps = [0.0]
            for q, s in perm: ps.append(ps[-1] + s*math.log(q))
            if max(ps) - min(ps) >= lc - 1e-9: continue
            if len(set(round(p, 9) for p in ps)) < len(ps): continue
            best = perm; break
        if best is None: continue
        ws = [W[q] for q, _ in best]
        f = lambda z: -bound_float(ws, np.exp(z))
        r = minimize(f, np.zeros(len(ws) + 1), method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 6000})
        out.append((-r.fun, best, np.exp(r.x)))
        if len(out) >= 60: break
    out.sort(key=lambda t: -t[0]); return out[:top]

def certify(best, al):
    den = 10**6; a = [fmpq(int(round(x/max(al)*den)), den) for x in al]
    ws = [arb(ppow(q)[0]).log()/arb(q).sqrt() for q, _ in best]
    A = [2*ws[i]*arb(a[i])*arb(a[i + 1]) for i in range(len(ws))]
    return (1/sum(1/x for x in A))/sum(arb(x)*arb(x) for x in a)

if __name__ == '__main__':
    pairs = [(10, 9.5), (10, 9.9), (47, 46.5), (47, 46.9), (47, 46.99), (49, 48.5), (fmpq(38), fmpq(227, 6)), (fmpq(759, 20), fmpq(227, 6))]
    print('# c, c_link, best rational cycles: path (q with sign), bound float >= arb')
    for c, cl in pairs:
        cf, clf = float(c), float(cl)
        res = search(cf, clf)
        txt = ' | '.join(f"{'·'.join((str(q) if s > 0 else '1/' + str(q)) for q, s in b)}: {v:.4f} >= {float(certify(b, al).lower()):.4f}" for v, b, al in res[:2])
        print(f'{cf:7.3f} {clf:8.4f}  {txt if res else "none found (k<=4)"}', flush=True)

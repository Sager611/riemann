# gstate.py (ccmfmm): even ground-state vector g of the truncation QW_lambda^NG (weil_matrix.py, closed forms in Arb), its
# rigorous Rayleigh quotient U = Q(g)/||g||^2 >= lambda_1^even (Galerkin), and approximate low eigenvalues (orientation only:
# e1, e2 of the even block, o1 of the odd block) used to place the certificate levels. g: inverse iteration on mid(T) with the
# seed of ../WeilSimple/cert_simple.py (ground_vec, rayleigh copied unchanged), rounded to exact dyadics.
# Usage: gstate.py c NG bits out.json [NO|0]     (NO = odd-block size for o1; 0 = no eigenvalues)
import sys, os, json, time, random, math
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weil_matrix import WeilMat

def ground_vec(T, n, its=8):   # = cert_simple.py (WS)
    rnd = random.Random(7); v = arb_mat(n, 1)
    for i in range(n): v[i, 0] = arb(rnd.uniform(-1, 1))
    for it in range(its):
        y = T.solve(v, algorithm='approx')
        nrm = sum((y[i, 0]*y[i, 0] for i in range(n)), arb(0)).sqrt()
        for i in range(n): v[i, 0] = (y[i, 0]/nrm).mid()
    return [v[i, 0] for i in range(n)]

def rayleigh(T, g):   # = cert_simple.py (WS): rigorous enclosure of g^T T g and g^T g
    n = len(g); num = arb(0)
    for i in range(n):
        s = arb(0)
        for j in range(n): s += T[i, j]*g[j]
        num += g[i]*s
    return num, sum((x*x for x in g), arb(0))

def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    return fmpq(int(s))

def low_eigs(T, k=2):
    ev = sorted([e.real.mid() for e in acb_mat(T).eig(algorithm="approx")], key=lambda z: z.mid())
    return [e.mid().str(15, radius=False) for e in ev[:k]]

def main():
    cq = parse_c(sys.argv[1]); NG = int(sys.argv[2]); ctx.prec = int(sys.argv[3]); out = sys.argv[4]
    NO = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    t0 = time.time(); WM = WeilMat(cq)
    Tg, _ = WM.block(NG, 'even'); t1 = time.time()
    g = ground_vec(Tg, NG + 1); num, gg = rayleigh(Tg, g); U = num/gg; t2 = time.time()
    res = dict(c=str(cq), NG=NG, prec=ctx.prec, U=U.str(15, radius=True), U_up=U.upper().str(12, radius=False),
               g=[list(map(str, x.mid().man_exp())) for x in g], t_block=round(t1 - t0, 1), t_g=round(t2 - t1, 1))
    if NO:
        res['even'] = low_eigs(Tg, 2); t3 = time.time()
        To, _ = WM.block(NO, 'odd'); res['odd'] = low_eigs(To, 1); res['NO'] = NO
        res['t_eig'] = round(time.time() - t2, 1)
    json.dump(res, open(out, 'w'))
    print(json.dumps({k: v for k, v in res.items() if k != 'g'}), flush=True)

if __name__ == '__main__':
    main()

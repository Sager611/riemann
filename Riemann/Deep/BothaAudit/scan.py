# scan.py (bothaaudit): float-level scan of e1, e2 (even) and o1, o2 (odd) of the truncated Weil form QW_lambda^N, c = lambda^2.
# Matrices: ../CCMFMM/weil_matrix.py's closed forms (copied unchanged as weil_matrix.py), with the per-atom trig values cached (FastWeilMat;
# validated entrywise against WeilMat.block, test mode). Eigenvalues: block inverse iteration (k vectors) + Rayleigh-Ritz, mid-point Arb at
# prec bits (orientation only, no rigour). Usage: scan.py out.jsonl NE NO prec c1 [c2 ...]   |  scan.py test c N prec
import sys, os, json, time, random
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weil_matrix import WeilMat

class FastWeilMat(WeilMat):
    def tables(s, N):
        s.wv = [s.w(j) for j in range(N + 1)]
        s.sn = [[(s.wv[j]*ln).sin() for j in range(N + 1)] for _, ln, _ in s.atoms]
        s.cs = [[(s.wv[j]*ln).cos() for j in range(N + 1)] for _, ln, _ in s.atoms]
    def Cf(s, i, j, k, sector):   # = WeilMat.C(i, j, x_k, sector) from cached sin/cos at atom k
        L = s.L; x = s.atoms[k][1]
        if i == j:
            if i == 0: return 1 - x/L
            w = s.wv[i]
            return (1 - x/L)*s.cs[k][i] + (-1 if sector == 'even' else 1)*s.sn[k][i]/(w*L)
        wi, wj = s.wv[i], s.wv[j]; sg = 1 if (i + j) % 2 == 0 else -1
        if sector == 'even':
            return s.Nj(i)*s.Nj(j)*sg*(wj*s.sn[k][j] - wi*s.sn[k][i])/(wi*wi - wj*wj)
        return (2/L)*sg*(wi*s.sn[k][j] - wj*s.sn[k][i])/(wi*wi - wj*wj)
    def fblock(s, N, sector):
        s.tables(N)
        idx = list(range(0, N+1)) if sector == 'even' else list(range(1, N+1))
        n = len(idx); T = arb_mat(n, n)
        pv = [s.pole_vec(j, sector) for j in idx]; ps = 2 if sector == 'even' else -2
        Sv = {j: s.S(j) for j in idx}
        for a, i in enumerate(idx):
            for b in range(a, n):
                j = idx[b]
                if i == j: ar = s.arch_diag(i, sector)
                else:
                    wi, wj = s.wv[i], s.wv[j]; sg = 1 if (i + j) % 2 == 0 else -1
                    if sector == 'even': ar = -2*s.Nj(i)*s.Nj(j)*sg*(wj*Sv[j] - wi*Sv[i])/(wi*wi - wj*wj)
                    else: ar = -2*(2/s.L)*sg*(wi*Sv[j] - wj*Sv[i])/(wi*wi - wj*wj)
                pr = arb(0)
                for k, (_, ln, wt) in enumerate(s.atoms): pr += wt*s.Cf(i, j, k, sector)
                v = ps*pv[a]*pv[b] + ar - 2*pr
                T[a, b] = v.mid(); T[b, a] = v.mid()
        return T

def lowest(T, k, its=4, seed=7):
    n = T.nrows(); rnd = random.Random(seed); V = arb_mat(n, k)
    for i in range(n):
        for j in range(k): V[i, j] = arb(rnd.uniform(-1, 1))
    hist = []
    for it in range(its):
        Y = T.solve(V, algorithm='approx')
        # Gram-Schmidt (mid-point)
        for j in range(k):
            for l in range(j):
                d = sum((Y[i, j]*Y[i, l] for i in range(n)), arb(0))
                for i in range(n): Y[i, j] = (Y[i, j] - d*Y[i, l]).mid()
            nr = sum((Y[i, j]*Y[i, j] for i in range(n)), arb(0)).sqrt()
            for i in range(n): Y[i, j] = (Y[i, j]/nr).mid()
        TY = T*Y; H = arb_mat(k, k)
        for a in range(k):
            for b in range(k):
                H[a, b] = sum((Y[i, a]*TY[i, b] for i in range(n)), arb(0)).mid()
        ev = sorted([e.real.mid() for e in acb_mat(H).eig(algorithm='approx')], key=lambda z: float(z.mid().log()) if z > 0 else -1e9)
        evs = sorted(ev, key=lambda z: z.mid())
        hist.append([float(e.mid().log()) / 2.302585092994046 for e in evs])  # log10 (assumes > 0)
        V = Y
    return evs, hist

def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    if '.' in s:
        a, b = s.split('.'); return fmpq(int(a + b), 10**len(b))
    return fmpq(int(s))

def main():
    if sys.argv[1] == 'test':
        cq = parse_c(sys.argv[2]); N = int(sys.argv[3]); ctx.prec = int(sys.argv[4])
        W = FastWeilMat(cq)
        for sec in ('even', 'odd'):
            T1, _ = W.block(N, sec); T2 = W.fblock(N, sec); mx = arb(0)
            for i in range(T2.nrows()):
                for j in range(T2.ncols()):
                    dd = abs(T1[i, j].mid() - T2[i, j]);  mx = max(mx, dd) if dd > mx else mx
                    assert T1[i, j].overlaps(T2[i, j]), (sec, i, j)
            print(sec, 'entries agree; max |mid diff| =', mx.mid().str(5), flush=True)
        return
    out = sys.argv[1]; NE = int(sys.argv[2]); NO = int(sys.argv[3]); ctx.prec = int(sys.argv[4])
    for cs in sys.argv[5:]:
        t0 = time.time(); W = FastWeilMat(parse_c(cs))
        Te = W.fblock(NE, 'even'); t1 = time.time()
        ee, he = lowest(Te, 3); t2 = time.time(); del Te
        To = W.fblock(NO, 'odd'); t3 = time.time()
        oo, ho = lowest(To, 2); t4 = time.time(); del To
        f = lambda z: z.mid().str(12, radius=False)
        rec = dict(c=cs, NE=NE, NO=NO, prec=ctx.prec, e1=f(ee[0]), e2=f(ee[1]), e3=f(ee[2]), o1=f(oo[0]), o2=f(oo[1]),
                   r_o1e1=float((oo[0]/ee[0]).mid()), r_e2e1=float((ee[1]/ee[0]).mid()),
                   conv_e=[round(h[0] - he[-1][0], 9) for h in he], conv_o=[round(h[0] - ho[-1][0], 9) for h in ho],
                   sec=[round(t1 - t0, 1), round(t2 - t1, 1), round(t3 - t2, 1), round(t4 - t3, 1)])
        open(out, 'a').write(json.dumps(rec) + '\n'); print(json.dumps(rec), flush=True)

if __name__ == '__main__':
    main()

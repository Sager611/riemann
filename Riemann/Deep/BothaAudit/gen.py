# gen.py (bothaaudit): generators of the truncated Weil blocks (even: n = N + 1, odd: n = N) for deig.c, mid-points of Arb values at prec bits.
# Entry form (algebra in REPORT.md sec. 2, = weil_matrix.py's block term by term):  T_ij = ps p_i p_j + a_i a_j (u_j - u_i)/(x_i - x_j)  (i != j),
# T_ii = d_i;  x_j = w_j^2, R_j = S_j + sum_k wt_k sin(w_j ln n_k);  even: a_j = N_j (-1)^j, u_j = -2 w_j R_j, ps = 2;
# odd: a_j = (2/L)^{1/2} (-1)^j w_j, u_j = -2 R_j/w_j, ps = -2.   Usage: gen.py c N prec outprefix  |  gen.py test c N prec
import sys, os
from flint import arb, arb_mat, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scan import FastWeilMat, parse_c

def gens(W, N, sector):
    W.tables(N); ev = sector == 'even'
    idx = range(0, N + 1) if ev else range(1, N + 1)
    rows = []
    for j in idx:
        w = W.wv[j]; sg = 1 if j % 2 == 0 else -1
        sig = arb(0); pr = arb(0)
        for k, (_, ln, wt) in enumerate(W.atoms):
            sig += wt*W.sn[k][j]; pr += wt*W.Cf(j, j, k, sector)
        R = W.S(j) + sig; p = W.pole_vec(j, sector)
        d = (2 if ev else -2)*p*p + W.arch_diag(j, sector) - 2*pr
        if ev: a = W.Nj(j)*sg; u = -2*w*R
        else: a = (2/W.L).sqrt()*sg*w; u = -2*R/w
        rows.append((d, p, a, u, w*w))
    return rows, (2 if ev else -2)

def hx(v):
    m, e = v.mid().man_exp(); m = int(m); e = int(e)
    return ('-' if m < 0 else '') + format(abs(m), 'x') + ' ' + str(e)

def main():
    if sys.argv[1] == 'test':
        cq = parse_c(sys.argv[2]); N = int(sys.argv[3]); ctx.prec = int(sys.argv[4]); W = FastWeilMat(cq)
        for sec in ('even', 'odd'):
            T, _ = W.block(N, sec); rows, ps = gens(W, N, sec); n = len(rows); mx = arb(0)
            for i in range(n):
                for j in range(n):
                    if i == j: v = rows[i][0]
                    else:
                        di, pi_, ai, ui, xi = rows[i]; dj, pj, aj, uj, xj = rows[j]
                        v = ps*pi_*pj + ai*aj*(uj - ui)/(xi - xj)
                    assert v.overlaps(T[i, j]), (sec, i, j)
                    dd = abs(v.mid() - T[i, j].mid()); mx = dd if dd > mx else mx
            print(sec, n, 'generator entries overlap weil_matrix.py; max |mid diff| =', mx.mid().str(5), flush=True)
        return
    cq = parse_c(sys.argv[1]); N = int(sys.argv[2]); ctx.prec = int(sys.argv[3]); pre = sys.argv[4]
    W = FastWeilMat(cq)
    for sec in ('even', 'odd'):
        rows, ps = gens(W, N, sec)
        with open(pre + '.' + sec, 'w') as f:
            f.write('%d %d\n' % (len(rows), ps))
            for r in rows: f.write(' '.join(hx(v) for v in r) + '\n')

if __name__ == '__main__':
    main()

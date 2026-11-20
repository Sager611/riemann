# hdx_count.py (explore64): rigorous lower bounds on N_c'(theta) = #{k : e_k(c') < theta} for Weil's form on [-(ln c')/2, (ln c')/2],
# both parities, from the arb Galerkin matrix in CCM's Fourier basis (../ExploreFunctional/code/wmat.py) and the inertia of A - theta I
# (symmetric LDL^T in arb, no pivoting; Sylvester).  Galerkin eigenvalues are upper bounds (min-max), so
#   #negative pivots of (A_N - theta I) <= N_c'(theta),  and e_k(c') < theta for every k <= that count.
# zeta: WeilSimple's audited WeilMat (rigorous tail remainder); other kinds: wmat, whose tail truncation (2^-(prec+40), ratio c^-2)
# is covered by padding every entry by 2^-(prec-20) (kept for zeta too: harmless).
import sys, math
sys.path.insert(0, 'Riemann/Deep/ExploreFunctional/code')
from flint import arb, fmpq, ctx
ctx.prec = int(sys.argv[3]) if len(sys.argv) > 3 else 200
import wmat
sys.path.insert(0, 'Riemann/Deep/WeilSimple')
from weil_matrix import WeilMat   # audited (WSGSAudit), rigorous archimedean tail s.trem

def ldl_inertia(T, theta):
    n = T.nrows(); pad = arb(0, arb(2)**(-(ctx.prec - 20)))
    A = [[T[i, j] + pad - (theta if i == j else 0) for j in range(n)] for i in range(n)]
    neg = pos = 0
    for k in range(n):
        d = A[k][k]
        if d.contains(0): return None
        if d < 0: neg += 1
        else: pos += 1
        inv = 1/d
        for i in range(k + 1, n):
            f = A[i][k]*inv
            if f == 0: continue
            for j in range(k + 1, i + 1):
                A[i][j] -= f*A[j][k]
    return neg, pos

def counts(cq, N, thetas, kind='zeta'):
    F = WeilMat(cq) if kind == 'zeta' else wmat.Form(cq, kind); out = {}
    for sec in ['even', 'odd']:
        T, idx = F.block(N, sec)
        out[sec] = [ldl_inertia(T, arb(th)) for th in thetas]
    return out

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    jobs = [(fmpq(47), ['5616/10000', '1023/10000', '1/1000', '1e-6']), (fmpq(227, 6), ['5755/10000', '2017/10000', '1/1000']), (fmpq(71), ['5026/10000']),
            (fmpq(3162277, 10**6), ['16245/10000', '1/1000']), (fmpq(5623413, 10**6), ['10579/10000']), (fmpq(7943282, 10**6), ['5104/10000']),
            (fmpq(6855655, 10**6), ['51538/10000', '1/1000']), (fmpq(1795, 100), ['31687/10000']), (fmpq(3198, 100), ['13715/10000'])]
    if len(sys.argv) > 2 and sys.argv[2] != 'all': jobs = [j for j in jobs if str(float(j[0]))[:5] in sys.argv[2].split(',')]
    print(f'# N = {N} (even: N+1 modes, odd: N modes), prec {ctx.prec}; c_link, theta, (#neg, #pos) even, odd, total lower bound on N_c(theta)')
    for cq, ths in jobs:
        thq = [fmpq(*map(int, t.split('/'))) if '/' in t else fmpq(1, int(round(1/float(t)))) for t in ths]
        res = counts(cq, N, thq)
        for i, t in enumerate(ths):
            ev = res['even'][i]; od = res['odd'][i]
            tot = (ev[0] if ev else '?', od[0] if od else '?')
            print(f'{float(cq):8.4f} theta={t:>12s} even {ev} odd {od}  N_c(theta) >= {sum(x for x in tot if x != "?")}', flush=True)

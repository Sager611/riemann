# scheck.py (bothaaudit): falsification test of Botha's (V1) input s(L) (stage S, Elston's schur_upper; imported trust): s must bound
# sup <C_b h, h>/||h||^2, C_b = sum_{log n < 2b} Lambda(n) n^{-1/2} (tau_{log n} + tau_{-log n}) on L^2(-b, b) (PROOFS P9; = minus the prime part
# of Q, whose form is -2 sum w_n C_f(log n)). Galerkin on our bases (../CCMFMM/weil_matrix.py closed forms C_ij(x), double precision): the top
# eigenvalue of 2 PR_N (PR_ij = sum_k w_k C_ij(log n_k)), even and odd blocks, rises with N to the sup; a value above s refutes s.
# Usage: scheck.py <glue dir> N1,N2,...
import sys, os, json, glob, math
import numpy as np
from fractions import Fraction

def atoms(a):
    out = []
    for n in range(2, int(math.exp(2*a)) + 2):
        if math.log(n) >= 2*a: break
        m = n; p = None
        for q in range(2, n + 1):
            if m % q == 0: p = q; break
        while m % p == 0: m //= p
        if m == 1: out.append((math.log(n), math.log(p)/math.sqrt(n)))
    return out

def top(a, N, sector):
    Lp = 2*a; idx = np.arange(0, N + 1) if sector == 'even' else np.arange(1, N + 1)
    w = 2*np.pi*idx/Lp; Nj = np.where(idx == 0, 1/math.sqrt(Lp), math.sqrt(2/Lp)) if sector == 'even' else np.full(len(idx), math.sqrt(2/Lp))
    sg = np.where((idx[:, None] + idx[None, :]) % 2 == 0, 1.0, -1.0)
    PR = np.zeros((len(idx), len(idx)))
    for x, wt in atoms(a):
        s = np.sin(w*x); c = np.cos(w*x)
        with np.errstate(divide='ignore', invalid='ignore'):
            D = w[:, None]**2 - w[None, :]**2
            if sector == 'even': C = Nj[:, None]*Nj[None, :]*sg*(w[None, :]*s[None, :] - w[:, None]*s[:, None])/D
            else: C = (2/Lp)*sg*(w[:, None]*s[None, :] - w[None, :]*s[:, None])/D
        dg = np.where(w == 0, 1 - x/Lp, (1 - x/Lp)*c + (-1 if sector == 'even' else 1)*np.where(w == 0, 0, s/np.where(w == 0, 1, w*Lp)))
        np.fill_diagonal(C, dg); PR += wt*C
    M = PR + PR.T   # = 2 PR (C_ij(x) symmetric in i, j for these bases)
    return float(np.linalg.eigvalsh(M)[-1])

def main():
    G = sys.argv[1]; Ns = [int(x) for x in sys.argv[2].split(',')]
    for f in sorted(glob.glob(G + '/g_*_odd.json'), key=lambda p: float(os.path.basename(p).split('_')[1])):
        L = os.path.basename(f).split('_')[1]
        if float(L) < 1.70: continue
        d = json.load(open(f)); s = float(Fraction(d['s'])) if '/' in str(d['s']) else float(d['s'])
        b = math.ceil(float(Fraction(L))*2**24)/2**24
        vals = {N: max(top(b, N, 'even'), top(b, N, 'odd')) for N in Ns}
        print(json.dumps(dict(L=L, b=b, s=s, s_float=d.get('s_float'), ritz=vals, margin=s - max(vals.values()))), flush=True)

if __name__ == '__main__':
    main()

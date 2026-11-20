# p3_polefree.py (explore57): Galerkin spectra of Weil's form Q and its pole-free part B = Q - pole (zeta), and of DH's form,
# in CCM's Fourier basis (wmat.py from explore19/39), to compare with the Berezin-Lieb / Weyl phase-space bounds (p2_berezin.py)
# and to check the rank-one (Weinstein-Aronszajn / Cauchy) interlacing of the pole: lam_k(Q^ev) in [mu_k(B^ev), mu_{k+1}(B^ev)],
# lam_k(Q^odd) in [mu_{k-1}(B^odd), mu_k(B^odd)].  Exploration numerics (float eigvalsh of arb midpoints), not a certificate.
import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'ExploreFunctional', 'code'))
from flint import ctx, fmpq
ctx.prec = int(sys.argv[3]) if len(sys.argv) > 3 else 200
import numpy as np
from wmat import Form
N = int(sys.argv[2])
def f2np(T):
    n = T.nrows(); return np.array([[float(T[i, j].mid()) for j in range(n)] for i in range(n)])
for cs in sys.argv[1].split(','):
    cq = fmpq(int(cs)) if '/' not in cs else fmpq(*map(int, cs.split('/')))
    out = []
    for kind in ['zeta', 'dh']:
        F = Form(cq, kind)
        for sec in ['even', 'odd']:
            T, idx = F.block(N, sec); A = f2np(T); ev = np.linalg.eigvalsh(A)
            row = f'{kind} {sec}: Q eig[0:3] = ' + ', '.join(f'{x:.3e}' for x in ev[:3])
            if kind == 'zeta':
                F.pole = False; Tb, _ = F.block(N, sec); F.pole = True
                eb = np.linalg.eigvalsh(f2np(Tb))
                row += ' | B = Q - pole eig[0:3] = ' + ', '.join(f'{x:.3e}' for x in eb[:3])
                row += f' | sum neg(B) = {-eb[eb < 0].sum():.4f}, #neg(B) = {int((eb < 0).sum())}'
            out.append(row)
    print(f'c = {cs}, N = {N}:'); [print('  ' + r) for r in out]; sys.stdout.flush()

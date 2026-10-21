# pyright: reportMissingImports=false
"""Independent consistency check (not part of the proof): the x-space Galerkin matrix of the Weil form Q (weil.py,
explicit formula in x-space, conventions checked against zeta zeros in Riemann/Deep/Weil/cert/check_ef.py) must
dominate the certificate's lower-bound operator M + V = gamma - C_p + poles + band(a - gamma) + U on the span of
the first N0 Legendre modes:  Q_N0 - (M+V)_N0 = compression of (T_tail - U) >= 0.
Usage: python cons16.py <grams.bin> <band.bin> <N0> <u0> <u1>"""
import sys
from flint import arb, ctx
import bc, weil
ctx.prec = 256
gf, bf, N0, u0, u1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
S = bc.Setup('16')
G = bc.load_mats(gf); Bd = bc.load_mats(bf)
L = S.L_(); gamma = Bd['gamma']
Qm = weil.weil_matrix(2*L, N0)            # weil.py uses the full width: support [-L, L]
for p in (0, 1):
    idx = list(range(p, N0, 2)); d = len(idx)
    Qp = weil.sub(Qm, idx)
    A0 = G[f'A0_{p}']; Kb = Bd[f'K_{p}']; Bb = Bd[f'Bb_{p}']
    Av = [Bd['Av'][j] for j in idx]; sg = 2 if p == 0 else -2; uu = u0 if p == 0 else u1
    Eh = [(1 if i == 0 else 0) - Bb[i, 0] for i in range(d)]
    MV = bc.mat_from(lambda i, k: (gamma if i == k else 0) - A0[i, k] + Kb[i, k] + sg*Av[i]*Av[k] + uu*Eh[i]*Eh[k], d)
    D = Qp - MV.mid()
    print(f"parity {p} (N0={N0}): lam_min(Q_N0) = {weil.lam_min(Qp).str(5, radius=False)}  lam_min((M+V)_N0) = {weil.lam_min(MV.mid()).str(5, radius=False)}"
          f"  lam_min(Q_N0 - (M+V)_N0) = {weil.lam_min(D).str(5, radius=False)} (must be >= 0 up to x-space quadrature error)", flush=True)

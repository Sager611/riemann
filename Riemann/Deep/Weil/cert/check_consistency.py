# pyright: reportMissingImports=false
# python-flint lives in a local venv (see README.md), not the system interpreter.

# Consistency: frequency-side envelope matrix M_R (cert.py) must satisfy M_R <= Q (x-space Weil matrix, weil.py),
# with Q - M_R small for large R. Checks normalizations/signs of the certificate independently.
from flint import arb, arb_mat, ctx
import cert, weil
ctx.prec = 192
cert.DIAG = False
cap = []
orig = cert.certified_shift
cert.certified_shift = lambda H: (cap.append(H), None)[1]
L = arb('1.5'); N = 40
for R in ['200', '600']:
    cap.clear()
    cert.certify(L, arb(R), N, 64, [0.0], log=lambda s: None)
    Q = weil.weil_matrix(L, N)
    for p, H in enumerate(cap):
        idx = list(range(p, N, 2))
        G = arb_mat(H.nrows(), H.nrows())
        for i in range(H.nrows()):
            for j in range(H.nrows()): G[i, j] = H[i, j].mid()
            G[i, i] += cert.LAST_TAU.mid()
        H = G
        D = weil.sub(Q, idx) - H
        print('R', R, 'even' if p == 0 else 'odd', 'lam_min(Q - M_R) =', weil.lam_min(D).str(5, radius=False),
              ' lam_min(Q)=', weil.lam_min(weil.sub(Q, idx)).str(5, radius=False), ' lam_min(M_R)=', weil.lam_min(H.mid()).str(5, radius=False))

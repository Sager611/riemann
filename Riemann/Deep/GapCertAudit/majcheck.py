"""gapcertaudit: entrywise test of Lemma T on the computed arrays of the certificate run (kappa = 167/50, r = 2, D = 72, K = 8):
max over M, N < D of |H_n[M][N]| / (4((r-1)!)^4 C(M+2r-1,2r-1) C(N+2r-1,2r-1) beta^(M+N) (2+M+N+4r)^n), n <= 8, and the same
ratio with the sharper (2+M+N)^n (own derivation).  Run from a copy of ExploreShape (imports its certJ.py and ../GapBar/gcert.py).
usage: majcheck.py [D] [K]  (K = 12 also tests Lemma T' of certJ5.py: base 3c/2, factor V_n, (M+N+4r)^n)"""
import sys, os, math
sys.path.insert(0, os.getcwd())
import certJ as CJ
from certJ import GC
from flint import arb, fmpq
from math import comb

D = int(sys.argv[1]) if len(sys.argv) > 1 else 72
r = 2; K = int(sys.argv[2]) if len(sys.argv) > 2 else 8
GC.TH = fmpq(2499, 10000)
c = GC.A_(fmpq(167, 50)) * arb.pi()
beta = int(float((1.5 * c + K + 1).upper())) + 1
GC.set_prec(128 + int(D * math.log2(max(beta / 4, 2))) + 64)
c = GC.A_(fmpq(167, 50)) * arb.pi()
H, imok = CJ.gamma_moments(c, D, beta, r, K)
print('D', D, 'beta', beta, 'imok', imok)
from math import factorial
cu = arb(c.upper())
Tk = [max(arb(factorial(k)), max(arb(1), arb(k) / cu) ** k) for k in range(K + 1)]
Um = [sum((comb(m, k) * Tk[k] * Tk[m - k] for k in range(m + 1)), arb(0)) for m in range(K + 1)]
Vn = [sum((comb(n, j) * Um[n - j] for j in range(n + 1)), arb(0)) for n in range(K + 1)]
for n in range(K + 1):
    wT = 0.0; wh = None
    for M in range(D):
        for N in range(D):
            bT = 4 * Vn[n] * comb(M + 2 * r - 1, 2 * r - 1) * comb(N + 2 * r - 1, 2 * r - 1) * (3 * cu / 2) ** (M + N) * arb(M + N + 4 * r) ** n
            q = float((abs(H[n][M][N]).upper() / bT).upper())
            if q > wT:
                wT = q; wh = (M, N)
    print("n=%d  Lemma T' max ratio %.3e at %s" % (n, wT, wh), flush=True)
for n in range(K + 1):
    worst = [0.0, 0.0]; where = [None, None]; far = 0.0
    for M in range(D):
        for N in range(D):
            base = 4 * comb(M + 2 * r - 1, 2 * r - 1) * comb(N + 2 * r - 1, 2 * r - 1) * arb(beta) ** (M + N)
            v = abs(H[n][M][N]).upper()
            for i, e in enumerate((2 + M + N + 4 * r, 2 + M + N)):
                q = float((v / (base * arb(e) ** n)).upper())
                if q > worst[i]:
                    worst[i] = q; where[i] = (M, N)
                if i == 0 and max(M, N) >= D - 8 and q > far:
                    far = q
    print('n=%d  max ratio (code majorant) %.3e at %s;  (sharper (2+M+N)^n) %.3e at %s;  max over max(M,N) >= D-8: %.3e' % (n, worst[0], where[0], worst[1], where[1], far), flush=True)

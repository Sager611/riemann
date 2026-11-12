"""gapbar: TWISTED Gamma arrays. For complex amplifiers write f = e^{icU/2} g (orientation L; e^{-icU/2} for the other).
Then conj f(U1) f(U2) U_L conj(U_R) = conj g(U1) g(U2) * K~ with the pointwise REAL kernel
  K~ = e^{ic(r2-r1)/2} U_L conj(U_R) = 2cos(c(X1-X2-X3+X4)/2) - 2cos(c(X1-X2+X3-X4)/2 + c(tau1-tau2)),  r1 = X1+X2, r2 = X3+X4,
so sup over complex f of the BM pencil = sup over real g of the real pencil with K~ (a Hermitian form with a real symmetric
kernel has the same top over complex and real vectors).  The twist only shifts the phases in the four terms of
certr.gamma_arrays_r: (X1,X3) += (-c/2,+c/2), (X2,X4) += (-c/2,+c/2); all |phases| <= c, so the tail majorant
|Gam_ij| <= 4(i+j+8)^2((i+2r)(j+2r))^(2r) gam^(i+j), gam >= 2c+2, still holds.  IMCHECK: Im parts contain 0 (Arb)."""
import math
from flint import arb, acb, arb_mat, acb_poly, fmpq, ctx
import gcert as GC
IMCHECK = None

def gamma_arrays_tw(c, D, gam, r):
    I = acb(0, 1); ic = I * c
    h2 = c / 2   # twist e^{ic(X3+X4-X1-X2)/2}: (X1,X3) phases += (-c/2,+c/2), (X2,X4) phases += (-c/2,+c/2)
    terms = [(1, (h2, -h2), (-h2, h2), (acb(0), acb(0)), True), (1, (-h2, h2), (h2, -h2), (acb(0), acb(0)), True),
             (-1, (h2, h2), (-h2, -h2), (ic, -ic), False), (-1, (-h2, -h2), (h2, h2), (-ic, ic), False)]
    wr = [[arb(math.factorial(a + r - 1) * math.factorial(b + r - 1)) / (GC.FAC[a] * GC.FAC[b]) / arb(gam) ** 0 for b in range(D)] for a in range(D)]
    g = [[[acb(0)] * D for _ in range(D)] for _ in range(3)]
    for s, (a1, a3), (a2, a4), (z1, z2), zero in terms:
        A = lambda x: x if isinstance(x, arb) else arb(x)
        us = [GC.uhat(A(a1), A(a3), k, z1, zero, D, gam) for k in range(3)]
        ws = [GC.uhat(A(a2), A(a4), k, z2, zero, D, gam) for k in range(3)]
        us = [[[U[a][b] * wr[a][b] for b in range(D)] for a in range(D)] for U in us]
        ws = [[[W[a][b] * wr[a][b] for b in range(D)] for a in range(D)] for W in ws]
        combos = [[(0, 0, 1)], [(1, 0, 1), (0, 1, 1)], [(2, 0, 1), (1, 1, 2), (0, 2, 1)]]
        for gi, cl in enumerate(combos):
            for (k, l, mult) in cl:
                C = GC.conv2(us[k], ws[l], D)
                for M in range(D):
                    for N in range(D):
                        g[gi][M][N] += s * mult * C[M][N]
    gp = [arb(gam) ** n for n in range(2 * D)]
    g = [[[x[M][N] * gp[M + N] for N in range(D)] for M in range(D)] for x in g]
    def sh(G):
        return [[(G[M - 1][N] * (M + 2 * r - 1) if M > 0 else 0) + (G[M][N - 1] * (N + 2 * r - 1) if N > 0 else 0) for N in range(D)] for M in range(D)]
    a = g[0]
    b = [[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(sh(a), g[1])]
    s1 = sh(g[1]); ss = sh(sh(a))
    cc = [[ss[M][N] + 2 * s1[M][N] + g[2][M][N] for N in range(D)] for M in range(D)]
    sym = lambda G: [[(G[M][N].real + G[N][M].real) / 2 for N in range(D)] for M in range(D)]
    # exact reality check: the twisted kernel is pointwise real, so Im of the Hermitian part must contain 0
    him = lambda G: [[(G[M][N].imag - G[N][M].imag) / 2 for N in range(D)] for M in range(D)]
    global IMCHECK
    IMCHECK = all(x.contains(0) for G in (a, b, cc) for row in him(G) for x in row)
    return {'a': sym(a), 'b': sym(b), 'c': sym(cc)}


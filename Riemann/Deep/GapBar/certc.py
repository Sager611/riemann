"""gapbar: Arb certificate of h1 > 1 (Lambda > kappa) for ONE COMPLEX amplifier, general r (BM, arXiv 1410.3635v3).
Copy of GapCeil certr.py with the Hermitian form: for complex P the main term is sesquilinear (BBLR Thm, BM l.500-517:
sum a_h conj(a_k)), Q = sum_ij Gam_ij conj(I^i phi)^T W_p (I^j phi) with Gam the full Hermitian Taylor array
(GapCeil kept Re Gam only, i.e. real P).  Orientation: if the true form conjugates the other side, conj(P) works instead.
Original header:
Arb certificate of h1 > 1 (Lambda > kappa) for ONE amplifier with general r (BM, arXiv 1410.3635v3, Thm smoothed).
Same reduction as gcert.py (DERIVATION.md 1-2) with the r-weights X^{r^2-1}(X1X2X3X4)^{r-1}:
  phi = I^{2r} f, Q_j = sum_{i,j} Gam^(r)_ij <I^i phi, X^{r^2-1+p} I^j phi>, Gam^(r) = uhat^(r) * what^(r),
  uhat^(r)_ab = (a+r-1)!(b+r-1)! [X1^a X3^b] u, shifts (r1+r2) -> (i+2r-1) G_{i-1,j} + (j+2r-1) G_{i,j-1}.
f(U) = P(1-U/theta) is a polynomial on [0, theta], so phi lies in the Legendre window V and Q_j = phi^T A phi exactly,
up to the operator-norm tail eps |phi|^2 of the truncated series. theta = 2499/10000 < 1/4 (BBLR range).
usage: certr.py logs/r2opt_....json kappa_cert [D=70]"""
import sys, json, math, time, resource
from math import comb, factorial
from flint import arb, acb, arb_mat, acb_poly, fmpq, ctx
import gcert as GC

def gamma_arrays_r(c, D, gam, r):
    I = acb(0, 1); ic = I * c
    terms = [(1, (c, -c), (0, 0), (acb(0), acb(0)), True), (1, (0, 0), (c, -c), (acb(0), acb(0)), True),
             (-1, (c, 0), (0, -c), (ic, -ic), False), (-1, (0, -c), (c, 0), (-ic, ic), False)]
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
    # Hermitian part (the exact Gam is Hermitian by the swap (X1,X2) <-> (X3,X4), rho -> conj rho)
    hre = lambda G: [[(G[M][N].real + G[N][M].real) / 2 for N in range(D)] for M in range(D)]
    him = lambda G: [[(G[M][N].imag - G[N][M].imag) / 2 for N in range(D)] for M in range(D)]
    return {'a': (hre(a), him(a)), 'b': (hre(b), him(b)), 'c': (hre(cc), him(cc))}

def main():
    t0 = time.time()
    js = json.load(open(sys.argv[1])); kap = GC.dec(sys.argv[2])
    r = js['r']; theta = GC.dec(js['theta']); d = js['d']; v = GC.dec(js['v'])
    q15 = lambda x: GC.dec(('%.15f' % float(x)) if abs(float(x)) >= 1 else ('%.16f' % float(x)))
    An = [q15(x) for x in js['Are']]; Bn = [q15(x) for x in js['Aim']]
    GC.TH = theta
    gam = int(float(2 * GC.A_(kap) * arb.pi().mid() + 3)) + 1
    D = int(sys.argv[3]) if len(sys.argv) > 3 else 70
    GC.set_prec(128 + int(D * math.log2(max(gam / 4, 2))) + 32)
    th = GC.th
    c = GC.A_(kap) * arb.pi()
    Gam = gamma_arrays_r(c, D, gam, r)
    N = d + 2 * r + 1
    q0 = r * r - 1
    n2 = N + D + 8
    I, T, G0 = GC.legendre_ops(n2 + q0 + 4)
    G0d = arb_mat([[G0[i] if i == j else arb(0) for j in range(n2 + q0 + 4)] for i in range(n2 + q0 + 4)])
    Wq = {}
    Tp = arb_mat([[arb(1) if i == j else arb(0) for j in range(n2 + q0 + 4)] for i in range(n2 + q0 + 4)])
    for q in range(q0 + 3):
        if q >= q0:
            Wq[q - q0] = GC.sub(G0d * Tp, n2, n2)
        Tp = Tp * T
    Is = GC.sub(I, n2, n2)
    # phi = I^{2r} f in the unnormalised basis; f coefficients An (length d+1)
    f = arb_mat([[GC.A_(An[k]) if k <= d else arb(0)] for k in range(n2)])
    fi = arb_mat([[GC.A_(Bn[k]) if k <= d else arb(0)] for k in range(n2)])
    phi = f; phii = fi
    for _ in range(2 * r):
        phi = Is * phi; phii = Is * phii
    # Q_h,p = sum_ij Gam_ij (I^i phi)^T W_p (I^j phi)
    Ys = [phi]; Zs = [phii]
    for k in range(1, D):
        Ys.append(Is * Ys[-1]); Zs.append(Is * Zs[-1])
    def Q(h, p):
        # Q_L = sum_ij ReH_ij S_ij - ImH_ij T_ij, S_ij = yr_i W yr_j + yi_i W yi_j, T_ij = yr_i W yi_j - yi_i W yr_j
        Gr, Gi = Gam[h]; W = Wq[p]
        WY = [W * Y for Y in Ys]; WZ = [W * Z for Z in Zs]
        tot = arb(0)
        for i in range(D):
            Yt = Ys[i].transpose(); Zt = Zs[i].transpose()
            for j in range(D):
                S = (Yt * WY[j])[0, 0] + (Zt * WZ[j])[0, 0]
                T = (Yt * WZ[j])[0, 0] - (Zt * WY[j])[0, 0]
                tot += Gr[i][j] * S - Gi[i][j] * T
        return tot
    vals = {hp: Q(*hp) for hp in [('a', 0), ('a', 1), ('a', 2), ('b', 0), ('b', 1), ('c', 0)]}
    # tail bound (i or j >= D): |Gam^(r)_ij| <= 4(i+j+8)^2 (i+1)(i+2r)^{2r-2}... use crude 4(i+j+8)^2 ((i+2r)(j+2r))^{2r} gam^{i+j}
    qg = arb(gam) * th
    tail = arb(0)
    for i in range(340):
        for j in range(340):
            if i >= D or j >= D:
                tail += 4 * (i + j + 8) ** 2 * (arb(i + 2 * r) * (j + 2 * r)) ** (2 * r) * qg ** (i + j) / (GC.FAC[i] * GC.FAC[j])
    tail += arb(10) ** -60
    phin2 = sum(((phi[k, 0] ** 2 + phii[k, 0] ** 2) * th / (2 * k + 1) for k in range(n2)), arb(0))
    k2 = GC.A_(kap) ** 2; vv = GC.A_(v)
    Q0 = vals[('a', 0)]
    MS = vals[('a', 1)] + vals[('b', 0)]; MSS = vals[('a', 2)] + 2 * vals[('b', 1)] + vals[('c', 0)]
    Q1 = vv * vv * Q0 - 2 * vv * MS + MSS
    e0 = tail * phin2 * th ** q0
    e1 = tail * phin2 * th ** q0 * (vv * vv + 2 * abs(vv) * (1 + th) + (1 + th) ** 2 + 4)
    Q0b = Q0 + arb(0, e0.upper()); Q1b = Q1 + arb(0, e1.upper())
    h1 = Q0b / (k2 * Q1b)
    H = Q0b - k2 * Q1b
    ok = bool(H > 0 and Q1b > 0 and theta < fmpq(1, 4))
    print('r=%d theta=%s deg P=%d kappa=%s v=%s D=%d prec=%d' % (r, theta, d, kap, v, D, ctx.prec))
    print('  Q0 = %s\n  Q1 = %s\n  h1 = %s\n  H = Q0 - k^2 Q1 = %s\n  tail eps = %s ; float h1 = %.8f' % (Q0b.str(12), Q1b.str(12), h1.str(12), H.str(8), tail.str(3), js['h1_float']))
    print('  CERTIFIED Lambda > %s: %s   (%.0fs, peak %.0f MB)' % (kap, ok, time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
    print('  P(w) = sum_n (-1)^n (A_n + i B_n) psi_n(w); A_n, B_n (exact decimals) =', [str(x) for x in An], [str(x) for x in Bn])

if __name__ == '__main__':
    main()

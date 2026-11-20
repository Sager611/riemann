"""explore53 CERTIFICATE: Arb values of BM's main-term moments <S^m>, m <= K (= 2J), for ONE amplifier, midpoint mollifier,
and the J-th order criterion H = -sum_j a_j kappa^{2j} c_j(v) > 0  (c_j(v) = <(v - S)^{2j}>, CERT.md).
Reduction as GapBar certr.py/gtw.py (read-only import of gcert for Legendre operators and utilities):
  Q[(.)] = sum_{M,N} Gam_MN <I^M phi, X^{r^2-1} (.) I^N phi>,  phi = I^{2r} f,  f(U) = sum_n A_n P_n(2U/theta - 1);
  Gam for the tau-weight (tau1+tau2)^m: sum_{k+l=m} C(m,k) conv2(u_k, w_l) over the four twisted terms (gtw phases);
  S^n-kernels: H_n = sum_j C(n,j) sh^j(g_{n-j}) (sh = multiplication by r1+r2), and <S^m> = sum_i C(m,i) Q(H_{m-i}, X^i).
Truncation M, N < D with the majorant (CERT.md, Lemma T):
  |H_n[M][N]| <= 4 ((r-1)!)^4 C(M+2r-1,2r-1) C(N+2r-1,2r-1) beta^(M+N) (2+M+N+4r)^n,  beta = 3c/2 + K + 1,
  |<I^M phi, X^q I^N phi>| <= theta^q (theta^M/M!)(theta^N/N!) |phi|^2.
usage: certJ.py amp.json   (amp.json: r, theta, d, kappa, v, A (decimal strings), a (rational strings, J+1), D)"""
import sys, os, json, time, resource, math
from math import comb, factorial
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'GapBar'))
from flint import arb, acb, arb_mat, acb_poly, fmpq, ctx
import gcert as GC


def h_series_K(k, z, zero, Q):
    """Taylor coefficients (q < Q) of h(y) = int_0^{1-y} t^k e^{zt} dt (GapBar gcert.h_series, any k)."""
    if zero:
        out = [acb(0)] * Q
        for j in range(k + 2):
            if j < Q:
                out[j] = acb(comb(k + 1, j) * (-1) ** j) / (k + 1)
        return out
    ez = z.exp()
    e = [acb(1)]
    for n in range(1, Q):
        e.append(e[-1] * (-z) / n)
    g = [acb(0)] * Q
    for i in range(k + 1):
        a = acb(comb(k, i) * (-1) ** i)
        for n in range(Q - i):
            g[i + n] += a * e[n]
    g = [x * ez for x in g]
    m = (ez - 1) / z
    for j in range(1, k + 1):
        m = (ez - j * m) / z
    return [m] + [-g[q - 1] / q for q in range(1, Q)]


def uhat_K(al1, al3, k, z, zero, D, gam, cache={}):
    """GapBar gcert.uhat for any k, as matrix products: U = B1 Hank(hh) B3^T, B[a][n] = C(a,n) (i al gam^-1)^(a-n)."""
    from flint import acb_mat
    h = h_series_K(k, z, zero, 2 * D)
    gi = 1 / arb(gam)
    hh = [h[q] * GC.FAC[q] * gi ** q for q in range(2 * D)]
    key = (str(al1), str(al3), D, gam, ctx.prec)
    if key not in cache:
        I = acb(0, 1)
        p1 = [(I * al1 * gi) ** n for n in range(D)]
        p3 = [(I * al3 * gi) ** n for n in range(D)]
        B1 = acb_mat([[comb(a, n) * p1[a - n] if n <= a else acb(0) for n in range(D)] for a in range(D)])
        B3t = acb_mat([[comb(b, n) * p3[b - n] if n <= b else acb(0) for b in range(D)] for n in range(D)])
        cache.clear(); cache[key] = (B1, B3t)
    B1, B3t = cache[key]
    Hk = acb_mat([[hh[n1 + n2] for n2 in range(D)] for n1 in range(D)])
    U = B1 * Hk * B3t
    return [[U[a, b] for b in range(D)] for a in range(D)]


def gamma_moments(c, D, gam, r, K):
    """twisted Gamma arrays g_m (tau-weight (tau1+tau2)^m), m <= K, as in gtw.gamma_arrays_tw; returns real symmetric
    H_n (n <= K) and the imaginary-part check."""
    I = acb(0, 1); ic = I * c; h2 = c / 2
    terms = [(1, (h2, -h2), (-h2, h2), (acb(0), acb(0)), True), (1, (-h2, h2), (h2, -h2), (acb(0), acb(0)), True),
             (-1, (h2, h2), (-h2, -h2), (ic, -ic), False), (-1, (-h2, -h2), (h2, h2), (-ic, ic), False)]
    wr = [[arb(factorial(a + r - 1) * factorial(b + r - 1)) / (GC.FAC[a] * GC.FAC[b]) for b in range(D)] for a in range(D)]
    g = [[[acb(0)] * D for _ in range(D)] for _ in range(K + 1)]
    for s, (a1, a3), (a2, a4), (z1, z2), zero in terms:
        A = lambda x: x if isinstance(x, arb) else arb(x)
        us = [[[U[a][b] * wr[a][b] for b in range(D)] for a in range(D)] for U in (uhat_K(A(a1), A(a3), k, z1, zero, D, gam) for k in range(K + 1))]
        ws = [[[W[a][b] * wr[a][b] for b in range(D)] for a in range(D)] for W in (uhat_K(A(a2), A(a4), k, z2, zero, D, gam) for k in range(K + 1))]
        for m in range(K + 1):
            for k in range(m + 1):
                C = GC.conv2(us[k], ws[m - k], D)
                mult = s * comb(m, k)
                for M in range(D):
                    for N in range(D):
                        g[m][M][N] += mult * C[M][N]
        del us, ws
    gp = [arb(gam) ** n for n in range(2 * D)]
    g = [[[x[M][N] * gp[M + N] for N in range(D)] for M in range(D)] for x in g]
    def sh(G):
        return [[(G[M - 1][N] * (M + 2 * r - 1) if M > 0 else 0) + (G[M][N - 1] * (N + 2 * r - 1) if N > 0 else 0) for N in range(D)] for M in range(D)]
    Hr = []; imok = True
    for n in range(K + 1):
        tot = [[acb(0)] * D for _ in range(D)]
        for j in range(n + 1):
            x = g[n - j]
            for _ in range(j):
                x = sh(x)
            cf = comb(n, j)
            tot = [[tot[M][N] + cf * x[M][N] for N in range(D)] for M in range(D)]
        imok = imok and all(((tot[M][N].imag - tot[N][M].imag) / 2).contains(0) for M in range(D) for N in range(D))
        Hr.append([[(tot[M][N].real + tot[N][M].real) / 2 for N in range(D)] for M in range(D)])
        del tot
    del g
    return Hr, imok


def main():
    t0 = time.time()
    js = json.load(open(sys.argv[1]))
    r = js['r']; theta = GC.dec(js['theta']); d = js['d']; kap = GC.dec(js['kappa']); v = GC.dec(js['v'])
    D = js.get('D', 90)
    a = [fmpq(Fr(x).numerator, Fr(x).denominator) for x in js['a']]
    J = len(a) - 1; K = 2 * J
    An = [GC.dec(x) for x in js['A']]
    GC.TH = theta
    c = GC.A_(kap) * arb.pi()
    beta = int(float((1.5 * c + K + 1).upper())) + 1          # majorant base (Lemma T), integer
    gam = beta                                                   # scaling of the arrays (any positive number works)
    GC.set_prec(128 + int(D * math.log2(max(gam / 4, 2))) + 64)
    th = GC.th
    H, imok = gamma_moments(c, D, gam, r, K)
    t1 = time.time()
    N = d + 2 * r + 1; q0 = r * r - 1
    n2 = N + D + 8; nf = n2 + q0 + K + 4
    I, T, G0 = GC.legendre_ops(nf)
    G0d = arb_mat([[G0[i] if i == j else arb(0) for j in range(nf)] for i in range(nf)])
    Wq = {}
    Tp = arb_mat([[arb(1) if i == j else arb(0) for j in range(nf)] for i in range(nf)])
    for q in range(q0 + K + 1):
        if q >= q0:
            Wq[q - q0] = GC.sub(G0d * Tp, n2, n2)
        Tp = Tp * T
    Is = GC.sub(I, n2, n2)
    f = arb_mat([[GC.A_(An[k]) if k <= d else arb(0)] for k in range(n2)])
    phi = f
    for _ in range(2 * r):
        phi = Is * phi
    cols = [phi]
    for k in range(1, D):
        cols.append(Is * cols[-1])
    Y = arb_mat([[cols[j][i, 0] for j in range(D)] for i in range(n2)])
    Yt = Y.transpose()
    Z = {p: Yt * Wq[p] * Y for p in range(K + 1)}
    def Q(n, p):
        G = H[n]; Zp = Z[p]
        return sum((G[i][j] * Zp[i, j] for i in range(D) for j in range(D)), arb(0))
    s = []
    for m in range(K + 1):
        s.append(sum((comb(m, i) * Q(m - i, i) for i in range(m + 1)), arb(0)))
    t2 = time.time()
    # tail (Lemma T): sum over (M,N) with M >= D or N >= D, up to 600 (beyond: geometric, < 1e-300)
    phin2 = sum((phi[k, 0] ** 2 * th / (2 * k + 1) for k in range(n2)), arb(0))
    bt = arb(beta) * th
    pr = lambda M: arb(comb(M + 2 * r - 1, 2 * r - 1))
    # tail (Lemma T), factorised: (2+M+N+4r)^n <= 2^(n-1)((1+2r+M)^n + (1+2r+N)^n), f(M) = C(M+2r-1,2r-1)(beta th)^M/M!;
    # sum over the tail region {M >= D or N >= D} of f(M)f(N)A_M^n = Fn_t F0 + Fn_D F0_t (no cancellation);
    # beyond M = 700 the ratio f(M+1)/f(M) < 1/40, so that remainder is < 2 f(700) (added as 1e-300).
    MM = 700
    fM = [pr(M) * bt ** M / GC.FAC[M] for M in range(MM)]
    tails = []
    for n in range(K + 1):
        A = [arb(1 + 2 * r + M) ** n for M in range(MM)]
        F0D = sum(fM[:D], arb(0)); F0t = sum(fM[D:], arb(0)) + arb(10) ** -300
        FnD = sum((fM[M] * A[M] for M in range(D)), arb(0)); Fnt = sum((fM[M] * A[M] for M in range(D, MM)), arb(0)) + arb(10) ** -300
        F0 = F0D + F0t
        if n == 0:
            tl = F0t * F0 + F0D * F0t
        else:
            tl = arb(2) ** (n - 1) * 2 * (Fnt * F0 + FnD * F0t)
        tails.append(4 * arb(factorial(r - 1)) ** 4 * tl)
    # error of s_m: sum_i C(m,i) |Q(H_{m-i}, X^i) tail| <= sum_i C(m,i) tails[m-i] theta^(q0+i) |phi|^2
    errs = [sum((comb(m, i) * tails[m - i] * th ** (q0 + i) for i in range(m + 1)), arb(0)) * phin2 for m in range(K + 1)]
    sb = [s[m] + arb(0, errs[m].upper()) for m in range(K + 1)]
    vv = GC.A_(v); k2 = GC.A_(kap) ** 2
    cj = [sum((comb(2 * j, i) * vv ** (2 * j - i) * (-1) ** i * sb[i] for i in range(2 * j + 1)), arb(0)) for j in range(J + 1)]
    aa = [arb(x) for x in a]
    Hval = -sum((aa[j] * k2 ** j * cj[j] for j in range(J + 1)), arb(0))
    ratio = Hval / cj[0]
    ok = bool(Hval > 0 and cj[0] > 0 and imok and theta < fmpq(1, 4))
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20
    out = {'r': r, 'theta': str(theta), 'd': d, 'kappa': str(kap), 'v': str(v), 'D': D, 'prec': ctx.prec, 'beta': beta,
           'J': J, 'a': js['a'], 'imag_check': imok, 's': [x.str(15) for x in sb], 'c': [x.str(15) for x in cj],
           'tails': [x.str(3) for x in tails], 'H': Hval.str(15), 'H_over_c0': ratio.str(15), 'certified': ok,
           'times_s': [round(t1 - t0), round(t2 - t1), round(time.time() - t2)], 'peak_MiB': round(peak)}
    print(json.dumps(out, indent=1))
    fn = js.get('out')
    if fn:
        json.dump(out, open(os.path.join(HERE, fn), 'w'), indent=1)


if __name__ == '__main__':
    main()

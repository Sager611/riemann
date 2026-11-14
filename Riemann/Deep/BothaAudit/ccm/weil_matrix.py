# weil_matrix.py: the Weil form Q (= QW_lambda of Connes-Consani-Moscovici, arXiv 2511.22755, with c = lambda^2) on the
# truncation spaces E_N, as Arb matrices with closed-form entries (derivation: DERIVATION.md sec. 2).
#   I = [-L/2, L/2], L = ln c. Even block: phi_j = N_j cos(w_j x), j = 0..N; odd block: psi_j = (2/L)^{1/2} sin(w_j x), j = 1..N;
#   w_j = 2 pi j/L, N_0 = L^{-1/2}, N_j = (2/L)^{1/2}. These span the even / odd parts of CCM's E_N = span{V_n : |n| <= N}.
#   Q(f) = 2<f,cosh(x/2)>^2 - 2<f,sinh(x/2)>^2 + (1/2 pi) int phi(t) |F(t)|^2 dt - 2 sum_{n<c} Lambda(n) n^{-1/2} C_f(log n),
#   phi(t) = -log pi + Re psi(1/4 + it/2), C_f(x) = int f(y) f(y - x) dy.
from flint import arb, acb, arb_mat, ctx, fmpq

def ppow(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
    return None

class WeilMat:
    def __init__(s, cq):
        s.cq = cq; s.c = arb(cq); s.L = s.c.log(); s.pi = arb.pi()
        s.atoms = []
        n = 2
        while fmpq(n) < cq:
            pe = ppow(n)
            if pe is not None: s.atoms.append((n, arb(n).log(), arb(pe[0]).log()/arb(n).sqrt()))
            n += 1
        # b_k = 2k + 1/2 and c^{-b_k}; truncated when c^{-b_K} < 2^-(prec+40); every tail term below is <= 16 c^{-b_k}
        s.tail = []; k = 0
        while True:
            b = arb(2*k) + arb(1)/2; t = s.c**(-b)
            s.tail.append((b, t)); k += 1
            if t < arb(2)**(-(ctx.prec + 40)) and k > 2: break
        bK = arb(2*k) + arb(1)/2
        s.trem = arb(0, (16*s.c**(-bK)/(1 - s.c**(-2))).upper())
        s._psi = {}
    def w(s, j): return 2*s.pi*j/s.L
    def psi(s, j):   # (Re psi, Im psi, Re psi') at 1/4 + i w_j/2
        if j not in s._psi:
            z = acb(arb(1)/4, s.w(j)/2); p = z.digamma(); p1 = acb(2).zeta(z)
            s._psi[j] = (p.real, p.imag, p1.real)
        return s._psi[j]
    def S(s, j):     # int_0^L e^{x/2} sin(w x)/(e^x - e^{-x}) dx = (1/2) Im psi - sum_k c^{-b_k} w/(b_k^2 + w^2)
        if j == 0: return arb(0)
        w = s.w(j); tot = arb(0)
        for b, t in s.tail: tot += t*w/(b*b + w*w)
        return s.psi(j)[1]/2 - tot - s.trem
    def arch_diag(s, j, sector):
        L = s.L; lp = -s.pi.log()
        if j == 0:
            z = acb(arb(1)/4); tot = arb(0)
            for b, t in s.tail: tot += t/(b*b)
            return lp + z.digamma().real + acb(2).zeta(z).real/(2*L) - 2*tot/L + s.trem
        w = s.w(j); rp, ip, rp1 = s.psi(j); tot = arb(0)
        if sector == 'even':
            for b, t in s.tail: tot += t*b*b/(b*b + w*w)**2
            return lp + rp + (ip/w + rp1/2)/L - 4*tot/L + s.trem
        for b, t in s.tail: tot += t/(b*b + w*w)**2
        return lp + rp - (ip/w - rp1/2)/L + 4*w*w*tot/L + s.trem
    def C(s, i, j, x, sector):   # C_ij(x) = int b_i(y) b_j(y - x) dy, 0 <= x <= L (closed form)
        L = s.L
        if i == j:
            if i == 0: return 1 - x/L
            w = s.w(i)
            return (1 - x/L)*(w*x).cos() + (-1 if sector == 'even' else 1)*(w*x).sin()/(w*L)
        wi, wj = s.w(i), s.w(j); sg = 1 if (i + j) % 2 == 0 else -1
        if sector == 'even':
            return s.Nj(i)*s.Nj(j)*sg*(wj*(wj*x).sin() - wi*(wi*x).sin())/(wi*wi - wj*wj)
        return (2/L)*sg*(wi*(wj*x).sin() - wj*(wi*x).sin())/(wi*wi - wj*wj)
    def Nj(s, j): return (1/s.L).sqrt() if j == 0 else (2/s.L).sqrt()
    def pole_vec(s, j, sector):   # <phi_j, cosh(x/2)> (even) or <psi_j, sinh(x/2)> (odd)
        w = s.w(j); sg = 1 if j % 2 == 0 else -1; sL4 = (s.L/4).sinh()
        if sector == 'even': return s.Nj(j)*sg*sL4/(w*w + arb(1)/4)
        return -2*s.Nj(j)*sg*w*sL4/(w*w + arb(1)/4)
    def block(s, N, sector, parts=False):
        idx = list(range(0, N+1)) if sector == 'even' else list(range(1, N+1))
        n = len(idx); T = arb_mat(n, n)
        pv = [s.pole_vec(j, sector) for j in idx]; ps = 2 if sector == 'even' else -2
        Sv = {j: s.S(j) for j in idx}
        for a, i in enumerate(idx):
            for b in range(a, n):
                j = idx[b]
                if i == j: ar = s.arch_diag(i, sector)
                else:
                    wi, wj = s.w(i), s.w(j); sg = 1 if (i + j) % 2 == 0 else -1
                    if sector == 'even': ar = -2*s.Nj(i)*s.Nj(j)*sg*(wj*Sv[j] - wi*Sv[i])/(wi*wi - wj*wj)
                    else: ar = -2*(2/s.L)*sg*(wi*Sv[j] - wj*Sv[i])/(wi*wi - wj*wj)
                pr = arb(0)
                for _, ln, wt in s.atoms: pr += wt*s.C(i, j, ln, sector)
                v = ps*pv[a]*pv[b] + ar - 2*pr
                T[a, b] = v; T[b, a] = v
        return T, idx

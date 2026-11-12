# wmat.py (explore19): Weil-type form on the window [-L/2, L/2], L = ln c, in CCM's Fourier basis (= Woodward-Lawson basis),
# generalised from Riemann/Deep/WeilSimple/weil_matrix.py (closed forms, DERIVATION.md sec. 2 there) to gamma factor
# (q/pi)^{s/2} Gamma((s + a)/2): archimedean symbol log(q/pi) + Re psi(z0 + it/2), z0 = 1/4 + a/2, kernel exponents b_k = 2k + 2 z0.
#   form: Q(f) = [pole: 2<f,cosh(x/2)>^2 - 2<f,sinh(x/2)>^2 if pole] + (1/2pi) int symbol |F|^2 dt - 2 sum_{n<c} wt(n) C_f(log n).
# Exploration code (midpoint use only; not a certificate).
from flint import arb, acb, arb_mat, ctx, fmpq

def ppow(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
    return None

def zeta_atoms(cq):
    out = []; n = 2
    while fmpq(n) < cq:
        pe = ppow(n)
        if pe is not None: out.append((n, arb(n).log(), arb(pe[0]).log()/arb(n).sqrt()))
        n += 1
    return out

def dh_coeffs(nmax):
    # Davenport-Heilbronn: a(n) periodic mod 5: 1, kappa, -kappa, -1, 0; kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1)
    k = ((10 - 2*arb(5).sqrt()).sqrt() - 2)/(arb(5).sqrt() - 1)
    per = [arb(0), arb(1), k, -k, arb(-1)]
    a = [arb(0)] + [per[n % 5] for n in range(1, nmax+1)]
    # -f'/f = sum b(n) n^-s:  sum_{de=n} a(d) b(e) = a(n) log n, a(1) = 1, b(1) = 0
    b = [arb(0)]*(nmax+1)
    for n in range(2, nmax+1):
        s = a[n]*arb(n).log()
        d = 2
        while d*d <= n:
            if n % d == 0:
                e = n//d
                s -= a[d]*b[e]
                if e != d: s -= a[e]*b[d]
            d += 1
        s -= a[n]*b[1]
        b[n] = s
    return a, b, k

def chi5_atoms(cq):
    # control: L(s, chi) mod 5, chi(2) = i (odd, complex); for real f only Re chi enters: weight log p Re(chi(p)^k)/sqrt n
    chi = {1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0)}
    out = []; n = 2
    while fmpq(n) < cq:
        pe = ppow(n)
        if pe is not None and pe[0] % 5 != 0:
            p, e = pe; z = complex(*chi[p % 5])**e
            re = int(round(z.real))
            if re != 0: out.append((n, arb(n).log(), re*arb(p).log()/arb(n).sqrt()))
        n += 1
    return out

def dh_atoms(cq):
    nmax = int(cq) + 1
    a, b, k = dh_coeffs(nmax)
    out = []; n = 2
    while fmpq(n) < cq:
        out.append((n, arb(n).log(), b[n]/arb(n).sqrt()))
        n += 1
    return out

class Form:
    def __init__(s, cq, kind='zeta'):
        s.cq = cq; s.c = arb(cq); s.L = s.c.log(); s.pi = arb.pi()
        if kind == 'zeta':
            s.q, s.z0, s.pole = arb(1), arb(1)/4, True; s.atoms = zeta_atoms(cq)
        elif kind == 'dh':
            s.q, s.z0, s.pole = arb(5), arb(3)/4, False; s.atoms = dh_atoms(cq)
        elif kind == 'chi5':
            s.q, s.z0, s.pole = arb(5), arb(3)/4, False; s.atoms = chi5_atoms(cq)
        else: raise ValueError(kind)
        s.lp = (s.q/s.pi).log()
        s.tail = []; k = 0
        while True:
            b = arb(2*k) + 2*s.z0; t = s.c**(-b)
            s.tail.append((b, t)); k += 1
            if t < arb(2)**(-(ctx.prec + 40)) and k > 2: break
        s._psi = {}
    def w(s, j): return 2*s.pi*j/s.L
    def psi(s, j):
        if j not in s._psi:
            z = acb(s.z0, s.w(j)/2); p = z.digamma(); p1 = acb(2).zeta(z)
            s._psi[j] = (p.real, p.imag, p1.real)
        return s._psi[j]
    def S(s, j):
        if j == 0: return arb(0)
        w = s.w(j); tot = arb(0)
        for b, t in s.tail: tot += t*w/(b*b + w*w)
        return s.psi(j)[1]/2 - tot
    def Nj(s, j): return (1/s.L).sqrt() if j == 0 else (2/s.L).sqrt()
    def arch_diag(s, j, sector):
        L = s.L; lp = s.lp
        if j == 0:
            z = acb(s.z0); tot = arb(0)
            for b, t in s.tail: tot += t/(b*b)
            return lp + z.digamma().real + acb(2).zeta(z).real/(2*L) - 2*tot/L
        w = s.w(j); rp, ip, rp1 = s.psi(j); tot = arb(0)
        if sector == 'even':
            for b, t in s.tail: tot += t*b*b/(b*b + w*w)**2
            return lp + rp + (ip/w + rp1/2)/L - 4*tot/L
        for b, t in s.tail: tot += t/(b*b + w*w)**2
        return lp + rp - (ip/w - rp1/2)/L + 4*w*w*tot/L
    def C(s, i, j, x, sector):
        L = s.L
        if i == j:
            if i == 0: return 1 - x/L
            w = s.w(i)
            return (1 - x/L)*(w*x).cos() + (-1 if sector == 'even' else 1)*(w*x).sin()/(w*L)
        wi, wj = s.w(i), s.w(j); sg = 1 if (i + j) % 2 == 0 else -1
        if sector == 'even':
            return s.Nj(i)*s.Nj(j)*sg*(wj*(wj*x).sin() - wi*(wi*x).sin())/(wi*wi - wj*wj)
        return (2/L)*sg*(wi*(wj*x).sin() - wj*(wi*x).sin())/(wi*wi - wj*wj)
    def pole_vec(s, j, sector):
        w = s.w(j); sg = 1 if j % 2 == 0 else -1; sL4 = (s.L/4).sinh()
        if sector == 'even': return s.Nj(j)*sg*sL4/(w*w + arb(1)/4)
        return -2*s.Nj(j)*sg*w*sL4/(w*w + arb(1)/4)
    def block(s, N, sector='even'):
        idx = list(range(0, N+1)) if sector == 'even' else list(range(1, N+1))
        n = len(idx); T = arb_mat(n, n)
        pv = [s.pole_vec(j, sector) for j in idx]; ps = 2 if sector == 'even' else -2
        Sv = {j: s.S(j) for j in idx}
        for a_, i in enumerate(idx):
            for b_ in range(a_, n):
                j = idx[b_]
                if i == j: ar = s.arch_diag(i, sector)
                else:
                    wi, wj = s.w(i), s.w(j); sg = 1 if (i + j) % 2 == 0 else -1
                    if sector == 'even': ar = -2*s.Nj(i)*s.Nj(j)*sg*(wj*Sv[j] - wi*Sv[i])/(wi*wi - wj*wj)
                    else: ar = -2*(2/s.L)*sg*(wi*Sv[j] - wj*Sv[i])/(wi*wi - wj*wj)
                pr = arb(0)
                for _, ln, wt in s.atoms: pr += wt*s.C(i, j, ln, sector)
                v = ar - 2*pr
                if s.pole: v += ps*pv[a_]*pv[b_]
                T[a_, b_] = v; T[b_, a_] = v
        return T, idx
    def ft_even(s, j, z):
        # F_j(z) = int phi_j(x) e^{izx} dx, phi_j = N_j cos(w_j x) on [-L/2, L/2]  (z complex acb)
        Lh = s.L/2; w = s.w(j); sg = 1 if j % 2 == 0 else -1
        return s.Nj(j)*sg*(z*Lh).sin()*2*z/(z*z - w*w)
    def ft_odd(s, j, z):
        # psi_j = (2/L)^{1/2} sin(w_j x): int sin(wx) e^{izx} dx = i int sin(wx) sin(zx) dx = i * sg * sin(zL/2) * 2 w/(w^2 - z^2)... (derived below)
        Lh = s.L/2; w = s.w(j); sg = 1 if j % 2 == 0 else -1
        # int_{-L/2}^{L/2} sin(wx) sin(zx) dx = sin((w-z)L/2)/(w-z) - sin((w+z)L/2)/(w+z); with wL/2 = pi j:
        # sin(pi j - zL/2) = -sg sin(zL/2)*(-1)... computed directly below for safety
        a = ((w - z)*Lh).sin()/(w - z) - ((w + z)*Lh).sin()/(w + z)
        return acb(0, 1)*s.Nj(j)*a

def to_mp(T, mp):
    n = T.nrows()
    M = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = mp.mpf(T[i, j].mid().str(radius=False, more=True, condense=0) if False else T[i, j].mid().str(50 + int(ctx.prec*0.302), radius=False))
    return M

"""gapjinf: exact per-gap quadratic forms of order J in sine/boundary coordinates (REPORT.md, 'u-representation').
Piece f on [0,g], f(0)=f(g)=0; w_k = k pi/g; even boundary data a_l = f^(2l)(0), b_l = f^(2l)(g), l = 1..L (L = floor((J-1)/2)).
Odd block coordinates z = (beta_k, k odd <= K; sig_l = a_l + b_l); even block z = (beta_k, k even <= K; del_l = b_l - a_l).
sine coefficient b_k = sum_{l<=L} eps_{k,l} + beta_k (beta_k = 0 for k > K), eps_{k,l} = (-1)^l (2/g) w_k^(-2l-1) tau_l,
tau = sig (k odd), -del (k even);  u_{k,n} = beta_k + sum_{n<l<=L} eps_{k,l}.
  ||f^(j)||^2 = (g/2) sum_k w_k^(2j) u_{k,n(j)}^2 + [j odd, j>=3] del_{(j-1)/2}^2 / g,   n(0)=0, n(j)=floor((j-1)/2);
  T_ml = [f^(2m+1) f^(2l)]_0^g = [m>=1] del_m del_l / g - (g/2) sum_k (-1)^(m+l) w_k^(2(m+l+1)) u_{k,m} eps_{k,l}.
Tails k > K summed exactly with Hurwitz zeta."""
import numpy as np
from scipy.special import zeta

PI = np.pi


def nj(j):
    return 0 if j == 0 else (j - 1) // 2


def tail(s, K, par, g):
    """sum_{k>K, k = par mod 2} (k pi/g)^(-s)"""
    k0 = K + 1 if (K + 1) % 2 == par else K + 2
    return (PI / g) ** (-s) * 2.0 ** (-s) * zeta(s, k0 / 2.0)


def block_mats(g, J, K, par):
    """returns (Ms[j], Ts{(m,l)}) for parity par (1 = odd k with sig, 0 = even k with del)"""
    L = (J - 1) // 2
    ks = np.arange(1, K + 1); ks = ks[ks % 2 == par]
    nb = len(ks); n = nb + L
    w = ks * PI / g
    sgn = 1.0 if par == 1 else -1.0           # tau_l = sgn * coord_l
    C = np.zeros((nb, L + 1))                  # C[k, l] = coefficient of coord_l in eps_{k,l}
    for l in range(1, L + 1):
        C[:, l] = sgn * (-1) ** l * (2 / g) * w ** (-2.0 * l - 1)
    def rvec(i, nn):                           # u_{k_i, nn} as a row vector in z
        r = np.zeros(n); r[i] = 1.0
        for l in range(nn + 1, L + 1):
            r[nb + l - 1] = C[i, l]
        return r
    def tailc(l, s):                            # coefficient of coord_l in eps_{k,l} times w^... for k>K: sgn(-1)^l (2/g) w^(-2l-1)
        return sgn * (-1) ** l * (2 / g)
    Ms = []
    for j in range(J + 1):
        M = np.zeros((n, n)); nn = nj(j)
        R = np.array([rvec(i, nn) for i in range(nb)])
        M += (g / 2) * (R * (w ** (2.0 * j))[:, None]).T @ R
        for l in range(nn + 1, L + 1):        # tail: (g/2) sum_k w^(2j) eps_l eps_l'
            for lp in range(nn + 1, L + 1):
                s = 2 * l + 2 * lp + 2 - 2 * j
                M[nb + l - 1, nb + lp - 1] += (g / 2) * tailc(l, 0) * tailc(lp, 0) * tail(s, K, par, g)
        if j % 2 == 1 and j >= 3 and par == 0:
            mm = (j - 1) // 2
            if mm <= L:
                M[nb + mm - 1, nb + mm - 1] += 1.0 / g
        Ms.append(0.5 * (M + M.T))
    Ts = {}
    for m in range(0, (J - 2) // 2 + 1):
        for l in range(1, L + 1):
            T = np.zeros((n, n))
            R = np.array([rvec(i, m) for i in range(nb)])
            e = np.zeros((nb, n)); e[:, nb + l - 1] = C[:, l]
            T += -(g / 2) * (-1) ** (m + l) * (R * (w ** (2.0 * (m + l + 1)))[:, None]).T @ e
            for lp in range(m + 1, L + 1):     # tail: u_{k,m} = sum_{l'>m} eps_{k,l'}
                s = 2 * lp - 2 * m
                T[nb + lp - 1, nb + l - 1] += -(g / 2) * (-1) ** (m + l) * tailc(lp, 0) * tailc(l, 0) * tail(s, K, par, g)
            if m >= 1 and par == 0:
                T[nb + m - 1, nb + l - 1] += 1.0 / g
            Ts[(m, l)] = 0.5 * (T + T.T)
    return Ms, Ts


def test_poly(J=5, K=400, g=2.7):
    """check against a polynomial piece f(t) = t (g - t) (1 + t^2 + c t^3) computed directly."""
    import numpy.polynomial.polynomial as Pl
    c = 0.3
    p = Pl.polymul([0, g, -1], [1, 0, 1, c])
    L = (J - 1) // 2
    der = [p]
    for _ in range(J + 1):
        der.append(Pl.polyder(der[-1]))
    norms = [Pl.polyval(g, Pl.polyint(Pl.polymul(der[j], der[j]))) for j in range(J + 1)]
    # coordinates
    out = {}
    for par in (1, 0):
        ks = np.arange(1, K + 1); ks = ks[ks % 2 == par]
        w = ks * PI / g
        tt = np.linspace(0, g, 20001)
        fv = Pl.polyval(tt, p)
        bk = np.array([2 / g * np.trapezoid(fv * np.sin(wk * tt), tt) for wk in w])
        al = [Pl.polyval(0, der[2 * l]) for l in range(1, L + 1)]; bl = [Pl.polyval(g, der[2 * l]) for l in range(1, L + 1)]
        tau = np.array([(a + b) if par == 1 else (b - a) for a, b in zip(al, bl)])
        sgn = 1.0 if par == 1 else -1.0
        eps = np.array([[sgn * (-1) ** l * (2 / g) * wk ** (-2.0 * l - 1) * tau[l - 1] for l in range(1, L + 1)] for wk in w])
        beta = bk - eps.sum(axis=1)
        print('parity %d: max |beta_k| (numerical; exact 0 for this degree-5 piece) = %.2e' % (par, np.abs(beta).max()))
        beta = 0 * beta
        out[par] = np.concatenate([beta, tau])
    tot = []
    for j in range(J + 1):
        v = 0
        for par in (1, 0):
            Ms, Ts = block_mats(g, J, K, par)
            v += out[par] @ Ms[j] @ out[par]
        tot.append(v)
    print('norms direct :', ' '.join('%.8g' % x for x in norms))
    print('norms blocks :', ' '.join('%.8g' % x for x in tot))
    for (m, l) in Ts:
        d = Pl.polyval(g, Pl.polymul(der[2 * m + 1], der[2 * l])) - Pl.polyval(0, Pl.polymul(der[2 * m + 1], der[2 * l]))
        v = 0
        for par in (1, 0):
            Ms, Ts2 = block_mats(g, J, K, par)
            v += out[par] @ Ts2[(m, l)] @ out[par]
        print('T_%d%d direct %.8g blocks %.8g' % (m, l, d, v))


if __name__ == '__main__':
    test_poly()

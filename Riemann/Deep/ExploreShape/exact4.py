"""explore53: EXACT per-gap validity for J = 4 even forms (Lemma 2 of REPORT.md), float evaluation.
Q(f) = a0|f|^2 + a1|f'|^2 + a2|f''|^2 + a3|f'''|^2 + a4|f''''|^2 + phi12 [f' f'']_0^g + phi32 [f''' f'']_0^g,
f in H^4(0,g), f(0) = f(g) = 0.  With sine coefficients b_k (w_k = k pi/g) and A0 = f''(0), A1 = f''(g) free:
  Q = (g/2) sum_k [ p(w) b^2 + beta(w) b d_k + gam(w) d_k^2 ] + (a3 + phi32)(A1 - A0)^2/g,  d_k = (2/g)((-1)^k A1 - A0),
  p = sum a_j w^{2j}, beta = -2a3 w^3 - 2a4 w^5 + phi12 w - phi32 w^3, gam = a3 + phi32 + a4 w^2.
Valid on gaps <= pi  <=>  p(w) > 0 for w >= 1, S_odd(g) = sum_{k odd} rho(w_k) >= 0, S_even(g) + (a3+phi32)/2 >= 0,
rho = gam - beta^2/(4p) = O(w^-2), for all g in (0, pi].  (J = 3 is the case a4 = phi32 = 0.)
usage: from exact4 import check4"""
import numpy as np


def parts(w, a, ph12, ph32):
    a = list(a) + [0.0] * (5 - len(a))
    p = sum(a[j] * w ** (2 * j) for j in range(5))
    beta = -2 * a[3] * w ** 3 - 2 * a[4] * w ** 5 + ph12 * w - ph32 * w ** 3
    gam = a[3] + ph32 + a[4] * w ** 2
    return p, beta, gam


def rho_stable(x, a, ph12, ph32):
    """rho = N(x)/(4 p(x)), x = w^2, with N = 4 p gam - x q^2 cancelled analytically (x^5, x^4 terms vanish identically)."""
    c = a[3] + ph32; e = 2 * a[3] + ph32
    n0 = 4 * a[0] * c
    n1 = 4 * (a[0] * a[4] + a[1] * c) - ph12 ** 2
    n2 = 4 * (a[1] * a[4] + a[2] * c) + 2 * ph12 * e
    n3 = 4 * (a[2] * a[4] + a[3] * c) - (e * e - 4 * ph12 * a[4])
    p = a[0] + x * (a[1] + x * (a[2] + x * (a[3] + x * a[4])))
    return (n0 + x * (n1 + x * (n2 + x * n3))) / (4 * p), n3


def check4(a, ph12, ph32, ng=600, KM=20000):
    a = list(a) + [0.0] * (5 - len(a))
    # (i) p > 0 on [1, inf): dense check + leading coefficient
    w = np.concatenate([np.linspace(1, 20, 200001), np.geomspace(20, 1e4, 20000)])
    p, _, _ = parts(w, a, ph12, ph32)
    pmin = p.min()
    k = np.arange(1, KM + 1)
    worst = (np.inf, None, None)
    for g in np.pi * np.concatenate([np.linspace(0.02, 1, ng)]):
        xk = (k * np.pi / g) ** 2
        rho, n3 = rho_stable(xk, a, ph12, ph32)
        # tails k > KM: rho ~ n3 g^2/(4 a4 pi^2 k^2) (a4 > 0); sum over odd / even k > KM ~ 1/(2 KM) each
        tail = (n3 * g * g / (4 * a[4] * np.pi ** 2)) / (2 * KM) if a[4] > 0 else 0.0
        so = rho[0::2].sum() + tail; se = rho[1::2].sum() + tail + 0.5 * (a[3] + ph32)
        m = min(so, se)
        if m < worst[0]:
            worst = (m, g, 'odd' if so < se else 'even')
    return pmin, worst

"""explore53: EXACT per-gap validity for J = 6 even forms (Lemma 4). Same boundary data as Lemma 3 (A = f''(0), f''(g);
B = f''''(0), f''''(g)); new terms a6|f^(6)|^2 (sine coeff of f^(6): -(w e + w^6 b - w^3 d)), phi52[f^(5) f''], phi54[f^(5) f''''],
f^(5)(x) = (B1-B0)/g + sum (e + w^5 b - w^2 d) cos(w x).  Per mode (x = w^2):
  qd = phi12 - (2a3+phi32) x + (phi52 - 2a4) x^2 - 2a5 x^3 - 2a6 x^4,  qe = phi14 - phi34 x + (2a5 + phi54) x^2 + 2a6 x^3,
  gdd = a3 + phi32 + (a4 - phi52) x + a5 x^2 + a6 x^3,  gde = phi34 + phi52 - (2a5 + phi54) x - 2a6 x^2,  gee = a5 + phi54 + a6 x,
  C = [[a3+phi32, (phi34+phi52)/2], [(phi34+phi52)/2, a5+phi54]].  Reduces to Lemma 3 when a6 = phi52 = phi54 = 0."""
from fractions import Fraction as Fr
import numpy as np
from exact5 import pmul, padd, scal, trim, peval


def numer6(a, ph):
    a = [Fr(x) for x in a] + [Fr(0)] * (7 - len(a))
    g = lambda k: Fr(ph.get(k, 0))
    f12, f14, f32, f34, f52, f54 = (g(k) for k in ('12', '14', '32', '34', '52', '54'))
    p = a[:]
    X = [Fr(0), Fr(1)]
    qd = [f12, -(2 * a[3] + f32), f52 - 2 * a[4], -2 * a[5], -2 * a[6]]
    qe = [f14, -f34, 2 * a[5] + f54, 2 * a[6]]
    gdd = [a[3] + f32, a[4] - f52, a[5], a[6]]
    gde = [f34 + f52, -(2 * a[5] + f54), -2 * a[6]]
    gee = [a[5] + f54, a[6]]
    N11 = trim(padd(scal(4, pmul(p, gdd)), scal(-1, pmul(X, pmul(qd, qd)))))
    N12 = trim(padd(scal(2, pmul(p, gde)), scal(-1, pmul(X, pmul(qd, qe)))))
    N22 = trim(padd(scal(4, pmul(p, gee)), scal(-1, pmul(X, pmul(qe, qe)))))
    DET = trim(padd(pmul(N11, N22), scal(-1, pmul(N12, N12))))
    C = [[a[3] + f32, (f34 + f52) / 2], [(f34 + f52) / 2, a[5] + f54]]
    return trim(p), N11, N12, N22, DET, C


def check6(a, ph, ng=400, KM=4000):
    p, N11, N12, N22, DET, C = numer6(a, ph)
    k = np.arange(1, KM + 1)
    worst = (np.inf, None, None)
    for gg in np.pi * np.linspace(0.02, 1, ng):
        x = (k * np.pi / gg) ** 2
        pv = sum(float(c) * x ** i for i, c in enumerate(p))
        r = [sum(float(c) * x ** i for i, c in enumerate(Nn)) / (4 * pv) for Nn in (N11, N12, N22)]
        So = np.array([[r[0][0::2].sum(), r[1][0::2].sum()], [r[1][0::2].sum(), r[2][0::2].sum()]])
        Se = np.array([[r[0][1::2].sum(), r[1][1::2].sum()], [r[1][1::2].sum(), r[2][1::2].sum()]]) + 0.5 * np.array([[float(C[i][j]) for j in range(2)] for i in range(2)])
        eo, ee = np.linalg.eigvalsh(So)[0], np.linalg.eigvalsh(Se)[0]
        if min(eo, ee) < worst[0]:
            worst = (min(eo, ee), gg, 'odd' if eo < ee else 'even')
    return {'worst': worst, 'deg': {n: len(P) - 1 for n, P in (('p', p), ('N11', N11), ('N12', N12), ('N22', N22), ('DET', DET))},
            'lead': {n: float(P[-1]) for n, P in (('N11', N11), ('N22', N22), ('DET', DET))},
            'C_eigs': list(np.linalg.eigvalsh(np.array([[float(C[i][j]) for j in range(2)] for i in range(2)])))}

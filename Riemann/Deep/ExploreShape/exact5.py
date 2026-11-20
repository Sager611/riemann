"""explore53: EXACT per-gap validity for J = 5 even forms (Lemma 3), exact rational numerators, float/Arb evaluation.
Q(f) = sum_{j<=5} a_j|f^(j)|^2 + phi12[f'f''] + phi14[f'f''''] + phi32[f'''f''] + phi34[f'''f''''] on [0,g], f(0)=f(g)=0.
Free boundary values A = (f''(0), f''(g)), B = (f''''(0), f''''(g)); d_k = (2/g)((-1)^k A1 - A0), e_k likewise with B.
Per mode (x = w^2): p b^2 + w qd b d + w qe b e + gdd d^2 + gde d e + gee e^2 with
  p = sum a_j x^j, qd = phi12 - (2a3+phi32)x - 2a4 x^2 - 2a5 x^3, qe = phi14 - phi34 x + 2a5 x^2,
  gdd = a3+phi32 + a4 x + a5 x^2, gde = phi34 - 2a5 x, gee = a5.
Minimising b: 2x2 matrix R = [[N11, N12],[N12, N22]]/(4p), N11 = 4p gdd - x qd^2, N12 = 2p gde - x qd qe, N22 = 4p gee - x qe^2.
min_b Q = (2/g)[ s^T S_odd s + t^T (S_even + C/2) t ],  s = (A0+A1, B0+B1), t = (A1-A0, B1-B0),
S_odd/even = sum_{k odd/even} R(x_k), C = [[a3+phi32, phi34/2],[phi34/2, a5]].
Sufficient for validity on all gaps g <= pi: p > 0 on [1,inf), S_odd(g) >= 0 and S_even(g) + C/2 >= 0 (PSD) for g in (0, pi].
usage: from exact5 import numer5, check5"""
from fractions import Fraction as Fr
import numpy as np


def pmul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def padd(*ps):
    n = max(len(p) for p in ps)
    out = [Fr(0)] * n
    for p in ps:
        for i, x in enumerate(p):
            out[i] += x
    return out


def scal(c, p):
    return [c * x for x in p]


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p


def numer5(a, ph):
    a = [Fr(x) for x in a]; f12, f14, f32, f34 = (Fr(ph[k]) for k in ('12', '14', '32', '34'))
    p = a[:]
    X = [Fr(0), Fr(1)]
    qd = [f12, -(2 * a[3] + f32), -2 * a[4], -2 * a[5]]
    qe = [f14, -f34, 2 * a[5]]
    gdd = [a[3] + f32, a[4], a[5]]
    gde = [f34, -2 * a[5]]
    gee = [a[5]]
    N11 = trim(padd(scal(4, pmul(p, gdd)), scal(-1, pmul(X, pmul(qd, qd)))))
    N12 = trim(padd(scal(2, pmul(p, gde)), scal(-1, pmul(X, pmul(qd, qe)))))
    N22 = trim(padd(scal(4, pmul(p, gee)), scal(-1, pmul(X, pmul(qe, qe)))))
    DET = trim(padd(pmul(N11, N22), scal(-1, pmul(N12, N12))))
    C = [[a[3] + f32, f34 / 2], [f34 / 2, a[5]]]
    return p, N11, N12, N22, DET, C


def peval(p, x):
    out = 0.0
    for c in reversed(p):
        out = out * x + float(c)
    return out


def check5(a, ph, ng=400, KM=4000):
    p, N11, N12, N22, DET, C = numer5(a, ph)
    xs = np.concatenate([np.linspace(1, 50, 200001), np.geomspace(50, 1e8, 20000)])
    pmin = min(peval(p, x) for x in xs[::50])
    roots = {}
    for name, P in (('p', p), ('N11', N11), ('N22', N22), ('DET', DET)):
        rr = np.roots([float(c) for c in reversed(P)])
        real = [z.real for z in rr if abs(z.imag) < 1e-9 * max(1, abs(z))]
        roots[name] = sorted(real)
    lead = {n: float(P[-1]) for n, P in (('N11', N11), ('N22', N22), ('DET', DET))}
    k = np.arange(1, KM + 1)
    worst = (np.inf, None, None)
    for g in np.pi * np.linspace(0.02, 1, ng):
        x = (k * np.pi / g) ** 2
        pv = sum(float(c) * x ** i for i, c in enumerate(p))
        r11 = sum(float(c) * x ** i for i, c in enumerate(N11)) / (4 * pv)
        r12 = sum(float(c) * x ** i for i, c in enumerate(N12)) / (4 * pv)
        r22 = sum(float(c) * x ** i for i, c in enumerate(N22)) / (4 * pv)
        So = np.array([[r11[0::2].sum(), r12[0::2].sum()], [r12[0::2].sum(), r22[0::2].sum()]])
        Se = np.array([[r11[1::2].sum(), r12[1::2].sum()], [r12[1::2].sum(), r22[1::2].sum()]]) + 0.5 * np.array([[float(C[0][0]), float(C[0][1])], [float(C[1][0]), float(C[1][1])]])
        m = min(np.linalg.eigvalsh(So)[0], np.linalg.eigvalsh(Se)[0])
        if m < worst[0]:
            worst = (m, g, 'odd' if np.linalg.eigvalsh(So)[0] < np.linalg.eigvalsh(Se)[0] else 'even')
    Cf = np.array([[float(C[0][0]), float(C[0][1])], [float(C[1][0]), float(C[1][1])]])
    return {'pmin': pmin, 'roots': roots, 'lead': lead, 'worst': worst, 'C_eigs': list(np.linalg.eigvalsh(Cf)),
            'deg': {n: len(P) - 1 for n, P in (('N11', N11), ('N12', N12), ('N22', N22), ('DET', DET))}}

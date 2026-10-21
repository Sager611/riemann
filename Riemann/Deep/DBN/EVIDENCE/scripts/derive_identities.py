"""Symbolic proof of every algebraic step behind the kernel certificates (needs sympy; ~1 min).
Each line prints an expression that must simplify to 0 (or a statement that must be True).
Output: EVIDENCE/04a_identities.txt"""
import signal
import sympy as sp
from common import EVID, Log

signal.alarm(900)
log = Log(EVID / '04a_identities.txt')
log('# Symbolic identities behind the kernel certificates (derive_identities.py, sympy %s)' % sp.__version__)
d, b, y, r, h, z, w, M = sp.symbols('d b y r h z w M', real=True)
D, v, B, R = sp.symbols('D v B R', nonnegative=True)


def K(hh, dd, bb):
    return (hh - bb) / (dd ** 2 + (hh - bb) ** 2) + (hh + bb) / (dd ** 2 + (hh + bb) ** 2)


def J(hh, rr, dd, bb):
    return (K(hh, dd + rr * hh, bb) + K(hh, dd - rr * hh, bb)) / 2


def blocks(Dx, vx, Bx, Rx):
    C = Dx + Rx + 1 + Bx
    A = C ** 2 + 4 * Dx * Rx - 4 * Bx
    return C, A, A ** 2 - 16 * Dx * Rx * C ** 2, A * (C - 2 * Bx) - 8 * Dx * Rx * C, (Dx + vx - Bx) ** 2 + 4 * Dx * Bx


sub = {D: d ** 2, v: y ** 2, B: b ** 2, R: r ** 2}
C, A, D4, G, Qs = blocks(D, v, B, R)


def zero(name, e):
    s = sp.cancel(sp.together(sp.expand(e)))
    log('%-74s -> %s' % (name, s))
    assert s == 0, name


zero('scaling: K_{hy}(hd,hb) - K_y(d,b)/h', K(h * y, h * d, h * b) - K(y, d, b) / h)
zero('scaling: J_{h,r}(hd,hb) - J_{1,r}(d,b)/h', J(h, r, h * d, h * b) - J(1, r, d, b) / h)
zero('evenness: K_y(-d,-b) - K_y(d,b) and J_r(-d,-b) - J_r(d,b)', (K(y, -d, -b) - K(y, d, b)) + (J(1, r, -d, -b) - J(1, r, d, b)) * sp.Symbol('t'))
zero('2K_y(d,b)/y - 4(D+v-B)/Q*   (D=d^2, v=y^2, B=b^2)', 2 * K(y, d, b) / y - (4 * (D + v - B) / Qs).subs(sub))
zero('Q* - (d^2+(y-b)^2)(d^2+(y+b)^2)', Qs.subs(sub) - (d ** 2 + (y - b) ** 2) * (d ** 2 + (y + b) ** 2))
zero('J_{1,r}(d,b) - 2G/D4   (R=r^2)', J(1, r, d, b) - (2 * G / D4).subs(sub))
zero('D4 - prod_{+-}((d+-r)^2+(1-b)^2)((d+-r)^2+(1+b)^2)', D4.subs(sub) - sp.Mul(*[((d + e) ** 2 + (1 - b) ** 2) * ((d + e) ** 2 + (1 + b) ** 2) for e in (r, -r)]))
X, Q, A1, A2, G1, G2, l1, l2, lam = sp.symbols('X Q A1 A2 G1 G2 l1 l2 lam')
zero('one shift: 4X/Q - 2lam G/A1 - 4[X A1 - (Q/2) lam G]/(Q A1)', 4 * X / Q - 2 * lam * G1 / A1 - 4 * (X * A1 - Q / 2 * lam * G1) / (Q * A1))
zero('two shifts: 4X/Q - 2l1G1/A1 - 2l2G2/A2 - 4Pmix/(Q A1 A2)', 4 * X / Q - 2 * l1 * G1 / A1 - 2 * l2 * G2 / A2 - 4 * (X * A1 * A2 - Q / 2 * (l1 * G1 * A2 + l2 * G2 * A1)) / (Q * A1 * A2))
u1 = 1 - z
Dz = z / (1 - z)
Ch = z + u1 * (R + 1 + B)
Ah = Ch ** 2 + 4 * z * u1 * R - 4 * B * u1 ** 2
D4h = Ah ** 2 - 16 * z * u1 * Ch ** 2 * R
Gh = Ah * (Ch - 2 * B * u1) - 8 * z * u1 * Ch * R
Wh = z + u1 * (v - B)
Qh = Wh ** 2 + 4 * z * u1 * B
C2, A2b, D42, G2b, Q2 = blocks(Dz, v, B, R)
for name, e in (('W_h - (1-z)(D+v-B)', Wh - u1 * (Dz + v - B)), ('C_h - (1-z) C', Ch - u1 * C2), ('A_h - (1-z)^2 A', Ah - u1 ** 2 * A2b),
                ('D4_h - (1-z)^4 D4', D4h - u1 ** 4 * D42), ('G_h - (1-z)^3 G', Gh - u1 ** 3 * G2b), ('Q_h - (1-z)^2 Q*', Qh - u1 ** 2 * Q2)):
    zero('homogenised block (Planat check_geometry.poly): ' + name + ', D=z/(1-z)', e)
zero('degree: (1-z)^5 P has W_h D4_h and Q_h G_h of total z-weight 1+4 = 2+3', (Wh * D4h - lam / 2 * Qh * Gh) - u1 ** 5 * ((Dz + v - B) * D42 - lam / 2 * Q2 * G2b))
zero('own pair: J_{1,r}(0,b) - s_r(b^2), s_r(w) = 2(1+R-w)/((1+R-w)^2+4Rw)', J(1, r, 0, b) - (2 * (1 + r ** 2 - b ** 2) / ((1 + r ** 2 - b ** 2) ** 2 + 4 * r ** 2 * b ** 2)))
den = (1 + R - w) ** 2 + 4 * R * w
zero('own pair: M den - 2(1+R-w) - [M w^2 + (M(2R-2)+2) w + M(1+R)^2 - 2(1+R)]', M * den - 2 * (1 + R - w) - (M * w ** 2 + (M * (2 * R - 2) + 2) * w + M * (1 + R) ** 2 - 2 * (1 + R)))
zero('own pair: den - [(1-w)^2 + R^2 + 2R(1+w)]   (so den > 0 for w >= 0, R > 0)', den - ((1 - w) ** 2 + R ** 2 + 2 * R * (1 + w)))
c, t, q0, qe, zz = sp.symbols('c t q0 qe zz', positive=True)
qt = (q0 + 2 / c) * sp.exp(-c * t) - 2 / c
zero("ODE: q(t) = (q0+2/c)e^{-ct} - 2/c solves q' = -2 - c q", sp.diff(qt, t) - (-2 - c * qt))
T = sp.log(1 + c * (q0 - qe) / (2 + c * qe)) / c
zero('terminal time: q(T) - qe with T = ln(1+z)/c, z = c(q0-qe)/(2+c qe)', sp.simplify(qt.subs(t, T) - qe))
wv = zz / (2 + zz)
zero('log bound: (1+w)/(1-w) - (1+z) with w = z/(2+z)  (ln(1+z) = 2 artanh w)', (1 + wv) / (1 - wv) - (1 + zz))
log('positivity facts used: Q* > 0 unless d = 0 and b = y; D4 > 0 because 0 <= b <= y <= rho < 1 (all rho < 1 checked);')
log('(1-z)^n > 0 on [0,1); lcm and binomial scalings and powers of 2 are positive integers. J_r and K_y are even in d and b,')
log('so D = d^2 >= 0 and 0 <= B <= v cover every real d and every other zero a+ib with |b| <= y.')
log('RESULT all identities hold')
log.save()

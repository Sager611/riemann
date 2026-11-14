#!/usr/bin/env python3
"""kcheck.py (explaudit): independent test of the kernel majorants K1 (triangle), K2 (trapezoids w+-), Thm-4 kernel,
on random zeros, with W from closed forms that are first checked against direct quadrature. Not part of the proof."""
import random, mpmath as mp
mp.mp.dps = 30
B = mp.mpf('0.8733819666132076'); b = 1 - B; H = mp.mpf(1) / 2
ell, a = mp.mpf(7), mp.mpf('4.9'); c = ell - a
def Wtri(s):
    m = s - B
    return (mp.exp(-m * ell) / a - (1 / a + 1 / c) * mp.exp(-m * c) + 1 / c) / m ** 2
def g(u): return (u + ell) / a if u <= -c else -u / c
s0 = mp.mpc('0.7', '3.1')
q1 = mp.quad(lambda u: g(u) * mp.exp((s0 - B) * u), [-ell, -c, 0])
print("triangle closed form vs quad:", mp.nstr(abs(Wtri(s0) - q1), 3))
n0, n1, n2 = ell / 2, mp.mpf(2), 2 * (1 / a + 1 / c)
mL = H - B
m0 = mp.quad(lambda u: g(u) * mp.exp(mL * u), [-ell, -c, 0])
J1 = mp.quad(lambda u: (1 / a if u < -c else 1 / c) * mp.exp(mL * u), [-ell, -c, 0])
J2 = mp.exp(-mL * ell) / a + (1 / a + 1 / c) * mp.exp(-mL * c) + 1 / c
m1 = J1 + abs(mL) * m0; m2 = J2 + 2 * abs(mL) * J1 + mL ** 2 * m0
def Wtrap(s, al, d, plus):
    if plus: pts = [(al * (1 - d), 1 / (al * d)), (al, -1 / (al * d)), (mp.mpf(1), -1 / d), (1 + d, 1 / d)]
    else: pts = [(al, 1 / (al * d)), (al * (1 + d), -1 / (al * d)), (1 - d, -1 / d), (mp.mpf(1), 1 / d)]
    return sum(m * t ** (s + 1) for t, m in pts) / (s * (s + 1))
al0, d0 = mp.mpf('0.05'), mp.mpf('0.3')
def wp(t):
    if t <= al0 * (1 - d0) or t >= 1 + d0: return 0
    if t < al0: return (t - al0 * (1 - d0)) / (al0 * d0)
    return 1 if t <= 1 else (1 + d0 - t) / d0
q2 = mp.quad(lambda t: wp(t) * t ** (s0 - 1), [al0 * (1 - d0), al0, 1, 1 + d0])
print("trapezoid closed form vs quad:", mp.nstr(abs(Wtrap(s0, al0, d0, True) - q2), 3))
random.seed(7); worst = {"K1": 0, "K1L": 0, "K2": 0, "K2L": 0, "T4": 0, "T4L": 0}
for it in range(4000):
    gam = mp.mpf(10) ** random.uniform(-3, 3.5)
    left = random.random() < 0.5
    beta = mp.mpf(random.uniform(float(b), 0.5)) if left else mp.mpf(random.uniform(0.5, float(B)))
    # K1
    Y = mp.mpf(10) ** random.uniform(0, 2); X = Y + ell; eY = mp.e * Y
    if left:
        G = min(m0, (m1 + m0 / eY) / gam, (m2 + 2 * m1 / eY + (2 / eY) ** 2 * m0) / gam ** 2); r = mp.exp(beta * X) * abs(Wtri(beta + 1j * gam)) / (mp.exp(H * X) * G); k = "K1L"
    else:
        G = min(n0, (n1 + n0 / eY) / gam, (n2 + 2 * n1 / eY + (2 / eY) ** 2 * n0) / gam ** 2); r = mp.exp(beta * X) * abs(Wtri(beta + 1j * gam)) / (mp.exp(B * X) * G); k = "K1"
    worst[k] = max(worst[k], r)
    # K2
    X = mp.mpf(random.uniform(4.16, 60)); al = mp.exp(-X / 2); d = mp.mpf(10) ** random.uniform(-6, math_log := -0.30103)
    e = H if left else B
    WBe = ((1 + d) ** (e + 1) - 1) / (d * e * (e + 1)); A1 = (1 + d) ** e * (1 + al ** e); A2 = ((1 + d) ** (e + 1) + 1) * (1 + al ** e) / d
    G = min(WBe, A1 / gam, A2 / gam ** 2)
    r = max(abs(Wtrap(beta + 1j * gam, al, d, pl)) for pl in (True, False)) * mp.exp(beta * X) / (mp.exp(e * X) * G)
    k = "K2L" if left else "K2"; worst[k] = max(worst[k], r)
    # Thm 4 kernel: x^{beta-1}|V(rho-1)| <= x^{-b} min(1/b, 1/g, ((1+h)^B+1)/(h g^2)) resp. x^{-1/2} min(2, 1/g, ((1+h)^{1/2}+1)/(h g^2))
    h = mp.mpf(random.uniform(0.3, 1)); X = mp.mpf(random.uniform(8, 200)); sv = beta - 1 + 1j * gam
    V = ((1 + h) ** (sv + 1) - 1) / (h * sv * (sv + 1))
    if left: G = min(mp.mpf(2), 1 / gam, ((1 + h) ** H + 1) / (h * gam ** 2)); r = mp.exp((beta - 1) * X) * abs(V) / (mp.exp(-X / 2) * G); k = "T4L"
    else: G = min(1 / b, 1 / gam, ((1 + h) ** B + 1) / (h * gam ** 2)); r = mp.exp((beta - 1) * X) * abs(V) / (mp.exp(-b * X) * G); k = "T4"
    worst[k] = max(worst[k], r)
print("max ratio actual/majorant (must be <= 1):", {k: mp.nstr(v, 5) for k, v in worst.items()})

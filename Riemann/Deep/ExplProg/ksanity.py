#!/usr/bin/env python3
"""ksanity.py -- float sanity checks of the kernel lemmas and of inequality (1.3) (not part of the proof).
 1. Lemma K1 (triangle kernel): x^beta |W(beta+i g)| <= x^B min(n0, k1/|g|, k2/g^2) on random (beta, g, Y); same for the left kernel.
 2. Lemma K2 (trapezoids w+ and w-): x^beta |W(beta+i g)| <= x^B min(WB, A1/|g|, A2/g^2) (and left) on random samples.
 3. Theorem 4 kernel: x^(beta-1)|V(rho-1)| <= x^-b min(1/b, 1/|g|, ((1+h)^B+1)/(h g^2)) (and left).
 4. Inequality (1.3) for the triangle kernel at x = 1e7: |phi(q) S(x;q,a) - x W(1)| <= Omega_1 + (phi-1) Omega + phi (W(0)c + D_q)
    with the certified bounds, for all q <= 30 and all a (exact sums by a sieve).
Run: uv run --no-project --with python-flint --with mpmath --with numpy python3 ksanity.py"""
import math, cmath, random, sys
import numpy as np
sys.path.insert(0, '.')
import cert
from cert import arb, A, Fr
B = float(cert.B.mid()); b = 1 - B
rnd = random.Random(7)
K = cert.KER
ell, a, c = 7.0, 4.9, 2.1
def W_tri(s):      # W(s) = int g(u) e^{(s-B)u} du (closed form, complex)
    m = s - B
    return (cmath.exp(-m * ell) / a - (1 / a + 1 / c) * cmath.exp(-m * c) + 1 / c) / (m * m)
worst = 0.0
for _ in range(20000):
    Y = rnd.uniform(1.0, 200.0); X = Y + ell
    left = rnd.random() < 0.3
    beta = rnd.uniform(b, 0.5) if left else rnd.uniform(0.5, B)
    g = 10 ** rnd.uniform(-3, 3) * (1 if rnd.random() < 0.5 else -1)
    lhs = math.exp(beta * X) * abs(W_tri(beta + 1j * g))
    (n0, k1, k2), (m0, k1L, k2L) = K.ks(arb(Y))
    if left:
        rhs = math.exp(X / 2) * min(float(m0.mid()), float(k1L.mid()) / abs(g), float(k2L.mid()) / g / g)
    else:
        rhs = math.exp(B * X) * min(float(n0.mid()), float(k1.mid()) / abs(g), float(k2.mid()) / g / g)
    worst = max(worst, lhs / rhs)
print("1. Lemma K1: max lhs/rhs over 20000 samples = %.4f (must be <= 1)" % worst)
assert worst <= 1 + 1e-9

def W_trap(s, u0, u1, v0, v1):
    return ((v1 ** (s + 1) - v0 ** (s + 1)) / (v1 - v0) - (u1 ** (s + 1) - u0 ** (s + 1)) / (u1 - u0)) / (s * (s + 1))
worst = 0.0
for _ in range(20000):
    X = rnd.uniform(math.log(64), 80.0); x = math.exp(X); al = math.exp(-X / 2)
    d = 10 ** rnd.uniform(-4, math.log10(0.5))
    plus = rnd.random() < 0.5
    u0, u1, v0, v1 = (al * (1 - d), al, 1.0, 1 + d) if plus else (al, al * (1 + d), 1 - d, 1.0)
    left = rnd.random() < 0.3
    beta = rnd.uniform(b, 0.5) if left else rnd.uniform(0.5, B)
    g = 10 ** rnd.uniform(-3, 4) * (1 if rnd.random() < 0.5 else -1)
    lhs = math.exp(beta * X) * abs(W_trap(complex(beta, g), u0, u1, v0, v1))
    (WB, A1, A2), (Wh, A1L, A2L) = cert.trap_consts(A(Fr(d).limit_denominator(10 ** 12)), arb(al))
    if left:
        rhs = math.exp(X / 2) * min(float(Wh.mid()), float(A1L.mid()) / abs(g), float(A2L.mid()) / g / g)
    else:
        rhs = math.exp(B * X) * min(float(WB.mid()), float(A1.mid()) / abs(g), float(A2.mid()) / g / g)
    worst = max(worst, lhs / rhs)
print("2. Lemma K2: max lhs/rhs over 20000 samples = %.4f (must be <= 1)" % worst)
assert worst <= 1 + 1e-6

worst = 0.0
for _ in range(20000):
    X = rnd.uniform(math.log(2974), 100.0); h = rnd.uniform(0.05, 1.0)
    left = rnd.random() < 0.3
    beta = rnd.uniform(b, 0.5) if left else rnd.uniform(0.5, B)
    g = 10 ** rnd.uniform(-3, 4) * (1 if rnd.random() < 0.5 else -1)
    s = complex(beta - 1, g)
    V = ((1 + h) ** (s + 1) - 1) / (h * s * (s + 1))
    lhs = math.exp((beta - 1) * X) * abs(V)
    if left:
        rhs = math.exp(-X / 2) * min(2.0, 1 / abs(g), ((1 + h) ** 0.5 + 1) / (h * g * g))
    else:
        rhs = math.exp(-b * X) * min(1 / b, 1 / abs(g), ((1 + h) ** B + 1) / (h * g * g))
    worst = max(worst, lhs / rhs)
print("3. Theorem 4 kernel: max lhs/rhs over 20000 samples = %.4f (must be <= 1)" % worst)
assert worst <= 1 + 1e-9

# 4. inequality (1.3) at x = 1e7
x = 10 ** 7; lo = int(math.ceil(x * math.exp(-ell)))
N = x
sv = np.ones(N + 1, dtype=np.bool_); sv[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
lam = np.zeros(N + 1)
for p in np.nonzero(sv)[0]:
    p = int(p); v = p
    while v <= N:
        lam[v] = math.log(p); v *= p
n = np.arange(lo, N + 1, dtype=np.float64)
u = np.log(n / x)
gval = np.where(u <= -c, (u + ell) / a, -u / c)
gval = np.clip(gval, 0, None)
w = gval * np.exp(-B * u)
wl = lam[lo:N + 1] * w
X = math.log(x); Y = X - ell
W1 = float(K.W1.mid()); W0 = float(K.W0.mid()); winf = float(K.winf.mid())
worst = 0.0
for q in range(1, 31):
    ph = sum(1 for k in range(1, q + 1) if math.gcd(k, q) == 1)
    L = math.log(max(q, 9))
    Z, ZL, z1, z1L = [float(v.mid()) for v in cert.omega_parts(K, arb(Y), arb(L))]
    om1 = x ** B * z1 + x ** 0.5 * z1L
    om = (x ** B * Z + x ** 0.5 * ZL) / 2
    om_q = om if q >= 3 else 0.0
    triv = W0 / (1 - math.exp(-2 * Y))
    omega = len([p for p in range(2, q + 1) if q % p == 0 and all(p % d for d in range(2, int(p ** 0.5) + 1))])
    Dq = winf * (omega * ell + math.log(q)) if q > 1 else 0.0
    bound = om1 + (ph - 1) * om_q + ph * (triv + Dq)
    for aa in range(1, q + 1):
        if math.gcd(aa, q) != 1: continue
        S = wl[(np.arange(lo, N + 1) % q) == (aa % q)].sum()
        dev = abs(ph * S - x * W1)
        worst = max(worst, dev / bound)
print("4. (1.3) at x = 1e7, q <= 30, all a: max |phi S - x W(1)| / bound = %.4f (must be <= 1); x W(1) = %.4g" % (worst, x * W1))
assert worst <= 1

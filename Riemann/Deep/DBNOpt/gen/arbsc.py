# Rigorous (Arb, python-flint) upper bounds for the nine scalar quantities of HB 2610.12803 Prop two-density / Sec kernel-checks.
import math
from fractions import Fraction as F
from flint import arb, ctx, fmpq
ctx.prec = 128
N1 = 65536
logs = [None, arb(0)] + [arb(n).log() for n in range(2, N1 + 1)]
_sv = bytearray([1]) * (N1 + 1); _sv[0] = _sv[1] = 0
for _i in range(2, int(N1 ** 0.5) + 1):
    if _sv[_i]: _sv[_i * _i::_i] = bytearray(len(_sv[_i * _i::_i]))
ppow = []  # (n = l^j, j, l) with n <= N1
for _l in range(2, N1 + 1):
    if _sv[_l]:
        _m, _j = _l, 1
        while _m <= N1:
            ppow.append((_m, _j, _l)); _m *= _l; _j += 1
logN1 = logs[N1]
def A(x): return arb(fmpq(x.numerator, x.denominator))
p = arb(61) / 20; sigp = arb(81) / 40; Nst = arb(69098); Lst = 2 * Nst.log() + 3
def T(j, v):
    s = arb(0)
    for k in range(j + 1):
        s += arb(math.factorial(j) // math.factorial(j - k)) * logN1 ** (j - k) / (v - 1) ** (k + 1)
    return arb(N1) ** (1 - v) * s
def arb_scalars(tm, tp):
    """tm = t^-, tp = t^+ (Fractions). Returns arb enclosures of Z_*, psi_*, psi2_*, delta_0..2, rho_0..2."""
    u = sigp + A(F(557, 100) * tm); u8 = sigp + A(F(278, 100) * tm); e = (p - 1) / 2 - A(F(278, 100) * tm)
    tpa4 = A(tp) / 4
    Z = arb(1); d0 = arb(0); d1 = arb(0); d2 = arb(0)
    for n in range(2, N1 + 1):
        ln = logs[n]; npw = (-u * ln).exp()
        Z += npw
        x = (tpa4 * ln * ln).expm1() * npw
        d0 += x; x = x * ln; d1 += x; d2 += x * ln
    psi = arb(0); psi2 = arb(0)
    for (n, j, l) in ppow:
        ll = logs[l]; npw = (-u * logs[n]).exp()
        psi += ll * npw; psi2 += j * ll * ll * npw
    rho0 = arb(101) / 100 * Nst ** (-p) * (Nst ** (e + 1) / (e + 1) + Nst ** e)
    return dict(Z=Z + T(0, u), psi=psi + T(1, u), psi2=psi2 + T(2, u),
                d0=d0 + T(0, u8), d1=d1 + T(1, u8), d2=d2 + T(2, u8),
                r0=rho0, r1=Lst * rho0, r2=(Lst * Lst + arb(10) ** -9) * rho0)
KEYS = ['Z', 'psi', 'psi2', 'd0', 'd1', 'd2', 'r0', 'r1', 'r2']
DEN = 10 ** 12
def is_strict_ub(U, v): return bool(arb(U) > v * DEN)
def ub_int(v):
    x = v * DEN
    U = int(x.mid().floor().unique_fmpz()) + 1
    while not (arb(U) > x): U += 1
    return U

import sys, random, math
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from lib import *
from c1809b_det import det_bareiss

def fund_unit(p):
    # smallest (t,u) > 0 with t^2 - p u^2 = +-4 ; eps = (t + u sqrt p)/2
    u = 1
    while True:
        for s in (-4, 4):
            t2 = p * u * u + s
            if t2 > 0:
                t = isqrt(t2)
                if t * t == t2:
                    return (t, u)
        u += 1

def h_real(p):
    t, u = fund_unit(p)
    logeps = math.log((t + u * math.sqrt(p)) / 2)
    s = -sum(leg(a, p) * math.log(math.sin(math.pi * a / p)) for a in range(1, (p - 1) // 2 + 1))
    h = s / logeps
    assert abs(h - round(h)) < 1e-6, (p, h)
    return round(h)

def unit_pow(p, e):
    # (t + u sqrt p)/2 ^ e = A + B sqrt p with A, B Fractions
    t, u = fund_unit(p)
    A, B = Fr(1), Fr(0)
    a, b = Fr(t, 2), Fr(u, 2)
    for _ in range(e):
        A, B = A * a + B * b * p, A * b + B * a
    return A, B

def consts(p):
    c = {'l2': leg(2, p)}
    if p % 4 == 1:
        h = h_real(p)
        c['h'] = h
        c['a'], c['b'] = unit_pow(p, h)
        c['a1'], c['b1'] = unit_pow(p, (2 - leg(2, p)) * h)
    else:
        c['hm'] = h_minus_p(p) if p > 3 else 1
    return c

def detM(p, lo, hi, base, x, y, z, w, use_w=True, use_z=True):
    rng = range(lo, hi + 1)
    M = [[x + base(j, k) + leg(j, p) * y + (leg(k, p) * z if use_z else 0) + (leg(j * k, p) * w if use_w else 0) for k in rng] for j in rng]
    return det_bareiss(M)

# Rigorous constants for rhow_exact.c (python-flint arb balls): upper-rounded binary64 weights w_m = Lambda(m)/sqrt(m)
# (smallest double >= upper endpoint of a 256-bit ball) and integer enclosures N_lo < 2^56 e^{j a} < N_hi (j = 1, 2).
from flint import arb, ctx
from fractions import Fraction as F
import math
ctx.prec = 256
a = arb(23029)/10000
lim = (2*a).exp(); assert lim > 100 and lim < 101
pp = []
for m in range(2, 101):
    p = next(d for d in range(2, m+1) if m % d == 0); x = m
    while x % p == 0: x //= p
    if x == 1: pp.append((m, p))
assert len(pp) == 35
def exact(b):          # arb point -> Fraction
    man, ex = b.man_exp(); return F(int(man))*(F(2)**int(ex)) if ex >= 0 else F(int(man), 2**(-int(ex)))
ws = []
for m, p in pp:
    w = arb(p).log()/arb(m).sqrt()
    up = exact(w.upper())
    f = float(up)
    if F(f) < up: f = math.nextafter(f, math.inf)
    assert F(f) >= up
    ws.append(f)
S = 56
cons = []
for j in (1, 2):
    c = (j*a).exp()*arb(2)**S
    lo = exact(c.lower()); hi = exact(c.upper())
    nl = math.floor(lo); assert math.floor(hi) == nl and nl+1 < 2**63
    cons.append((nl, nl+1))
with open('consts.h', 'w') as fh:
    fh.write('static const int NMP=35;\n')
    fh.write('static const unsigned long long MP[]={%s};\n' % ','.join(str(m) for m, _ in pp))
    fh.write('static const double WUP[]={%s};\n' % ','.join(x.hex() for x in ws))
    fh.write('static const int ESHIFT=%d;\n' % S)
    fh.write('static const unsigned long long EA_LO=%dULL, EA_HI=%dULL, E2A_LO=%dULL, E2A_HI=%dULL;\n' % (cons[0][0], cons[0][1], cons[1][0], cons[1][1]))
print('ok', cons, ws[:3])

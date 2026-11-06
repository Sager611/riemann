# Upper-rounded double weights w_m = Lambda(m)/sqrt(m) for prime powers m < e^{2a}, a = 23029/10000.
from mpmath import mp, mpf, log, sqrt, exp
import math, struct
mp.prec = 300
a = mpf(23029)/10000
lim = exp(2*a)
assert 100 < lim < 101
pp = []
for m in range(2, 101):
    p = next(d for d in range(2, m+1) if m % d == 0)
    x = m
    while x % p == 0: x //= p
    if x == 1: pp.append((m, p))
assert len(pp) == 35
out = []
for m, p in pp:
    w = log(p)/sqrt(m)
    f = float(w)
    if mpf(f) < w: f = math.nextafter(f, math.inf)
    assert mpf(f) >= w
    out.append((m, f.hex()))
with open('weights.h', 'w') as fh:
    fh.write('static const int NMP=%d;\n' % len(out))
    fh.write('static const unsigned long long MP[]={%s};\n' % ','.join(str(m) for m,_ in out))
    fh.write('static const double WUP[]={%s};\n' % ','.join(h for _,h in out))
print(len(out), out[:3])


import re, sys, json, time, math
sys.path.insert(0, 'Riemann/Deep/Ingham3')
t0 = time.time()
# (1) KLN Table 1 in the TeX vs cert.py's KLN_ROWS
tex = open('.scratch/pbpaudit/2101.12263/main.tex').read().split('\n')
print('header:', tex[1693][:150])
rows = []
for l in tex[1694:1716]:
    nums = re.findall(r'\$\s*([0-9.]+)\s*\$', l)
    if len(nums) == 8: rows.append(tuple(nums))
src = open('Riemann/Deep/Ingham3/cert.py').read()
m = re.search(r'KLN_ROWS = \[(.*?)\]\nclass', src, re.S)
cert_rows = re.findall(r"\('([^']*)','([^']*)','([^']*)','([^']*)','([^']*)','([^']*)','([^']*)','([^']*)'\)", m.group(1))
print('TeX rows', len(rows), 'cert rows', len(cert_rows))
bad = [ (a, b) for a, b in zip(rows, cert_rows) if tuple(float(x) for x in a) != tuple(float(x) for x in b)]
print('mismatches:', bad)
# (2) max |zeta(1/2+it)| on [0,3] (KLN a2 = 2.851 needs <= a2 - 0.63*6/e = 1.46042 for |t| <= e^-6.. and generally <= 1.461)
from flint import arb, acb, ctx
ctx.prec = 128
need = arb('2.851') - arb('0.63')*6/arb(1).exp()
print('a2 - a1*6/e =', need.str(10))
n = 3000; mx = arb(0)
for i in range(n):
    tb = arb(3*i/n).union(arb(3*(i+1)/n))
    v = abs(acb(arb('0.5'), tb).zeta())
    mx = mx.max(arb(v.upper()))
print('max |zeta(1/2+it)|, t in [0,3], upper bound over 3000 balls:', mx.str(10), ' < need:', bool(mx < need))
print('zeta(1/2) =', acb(0.5).zeta().real.str(12))
# (3) k = 7: a proven prime in (n^7, (n+1)^7) for 1 <= n <= 10^4 (flint fmpz.is_prime: proven primality)
from flint import fmpz
worstgap = (0, 0); cnt = 0
for nn in range(1, 10001):
    a = fmpz(nn)**7; b = fmpz(nn + 1)**7
    p = a + 1
    while not p.is_probable_prime(): p += 1
    assert p < b and p.is_prime(), nn
    g = int(p - a)
    if g > worstgap[0]: worstgap = (g, nn)
    cnt += 1
print('k=7: n = 1..10^4 all have a proven prime in (n^7,(n+1)^7):', cnt, ' largest gap p - n^7 =', worstgap)
print('time %.1fs' % (time.time() - t0))

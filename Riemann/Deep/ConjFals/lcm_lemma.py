# 2511.02106 Remark after Lemma "upper bounds for x_2" (TeX l. 1285-1347): the Lemma x_2 < sqrt(2x) (proved for kappa >= 3/2)
# "seems to hold for all kappa > 1".  By the proof (l. 1289-1292) it is equivalent to H(y) = F(y^2/2, 1) - F(y, 2) > 0 for all y > 2,
# F(x,a) = log[((x^{a+1}-1)^k - (x^a-1)^k)/((x^{a+1}-x)^k - (x^a-x)^k)]/(k log x)  (l. 1050).  Also tests the Remark's sufficient condition
# G(x,k) = [(x^2+2)^k - 2^k][(x+1)^k - 1]^2 - [(x^2+x+1)^k - (x+1)^k]^2 increasing in x.
from mpmath import mp, mpf, log, nstr
mp.dps = 40
def F(x, a, k):
    num = (x**(a+1) - 1)**k - (x**a - 1)**k
    den = (x**(a+1) - x)**k - (x**a - x)**k
    return log(num/den)/(k*log(x))
def H(y, k): return F(y*y/2, 1, k) - F(y, 2, k)
def G(x, k): return ((x*x+2)**k - 2**k)*((x+1)**k - 1)**2 - ((x*x+x+1)**k - (x+1)**k)**2
for k in ['1.001','1.01','1.05','1.1','1.2','1.3','1.4','1.49','1.5','2']:
    k = mpf(k); mn = (mpf(10)**9, None); gdec = None
    ys = [mpf(2) + mpf(10)**(e/mpf(20)) for e in range(-120, 161)]   # y in (2, 2 + 1e8]
    for y in ys:
        h = H(y, k)
        if h < mn[0]: mn = (h, y)
    xs = [mpf(10)**(e/mpf(40)) - 1 + mpf('1e-6') for e in range(0, 400)]
    prev = None
    for x in xs:
        g = G(x, k)
        if prev is not None and g < prev and gdec is None: gdec = x
        prev = g
    print('kappa=%s: min H(y) over y in (2,1e8] = %s at y=%s;  G first decreases at x=%s' % (nstr(k,4), nstr(mn[0],6), nstr(mn[1],6), nstr(gdec,6) if gdec else 'never'))

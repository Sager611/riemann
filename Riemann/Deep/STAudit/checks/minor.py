
import math
def lam(n, th):
    if th <= 1e-12 or abs(th - math.pi) <= 1e-12:
        return 6.0/((n+1)*(n+2)*(2*n+3))
    s = math.sin(th)
    return 4*s**3/((2*n+3)*s - math.sin((2*n+3)*th))
def mu(a, b):  # ST measure of theta in [a,b]: (2/pi) int sin^2 = (1/pi)[t - sin t cos t]
    F = lambda t: (t - math.sin(t)*math.cos(t))/math.pi
    return F(b) - F(a)
worst = 0; worstcase = None; viol = 0; checked = 0
for L in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0]:
    nbound = math.ceil(37.7/L)
    steps = 400
    for i in range(steps+1):
        a = (math.pi - L)*i/steps; b = a + L
        # least n with positive mean
        n = 1
        while mu(a,b) - lam(n,a) - lam(n,b) <= 0:
            n += 1
            if n > 10**6: break
        checked += 1
        r = n*L
        if r > worst: worst = r; worstcase = (L, a, n)
        if n > nbound: viol += 1
print('checked', checked, 'intervals; violations of n <= ceil(37.7/L):', viol)
print('max n*L =', round(worst,3), 'at (L, a, n) =', worstcase, ' (claim: n*L <= 37.7 + L)')

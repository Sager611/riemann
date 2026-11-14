# wyconst.py -- spot check of WY Prop. 3.5 constants (TeX l. 404-437): c, alpha(120409), -0.7702, 7.1476, 10724
from mpmath import mp, mpf, log, sqrt, pi, euler, exp
mp.dps = 30
c = euler + 2 - log(4*pi)
def alpha(x):
    S2 = x*log(x)**4/(64*pi**2)
    return S2*(log(x) + mpf('1.31'))/(2*x**mpf(1.5)) + (c - 2) + (8 - 4*c)/log(x) + 2*log(x)/x**(mpf(1)/6) + log(2*pi)*log(x)/sqrt(x)
print('c =', mp.nstr(c, 8), ' alpha(120409) =', mp.nstr(alpha(mpf(120409)), 8), ' alpha(1.43e47) =', mp.nstr(alpha(mpf('1.43e47')), 6))
print('first-order 2 + c - sqrt2 - 0.99154 sqrt2 =', mp.nstr(2 + c - sqrt(2) - mpf('0.99154')*sqrt(2), 8), ' second-order 2.6616 + 4 + 0.486 =', mpf('2.6616') + 4 + mpf('0.486'))
print('exp(7.1476/0.7702) =', mp.nstr(exp(mpf('7.1476')/mpf('0.7702')), 8))
xs = [mpf(120409) * mpf(1.0001)**k for k in range(0, 2000)]
print('alpha decreasing on [120409, 120409*1.0001^2000]:', all(alpha(xs[i]) > alpha(xs[i+1]) for i in range(len(xs)-1)))

# PYAudit: high-precision checks of the decimal constants in Patel-Yang 2302.13444v1 and Yang 2301.03165
# uv run --no-project --with mpmath --with python-flint python3 consts.py
import mpmath as mp
mp.mp.dps = 40
pi = mp.pi
out = []
def rep(name, val, printed, direction):
    # direction: 'upper' means printed must be >= true value (an upper bound is used), 'lower' the opposite
    ok = (printed >= val) if direction == 'upper' else (printed <= val)
    out.append('%-40s true=%s printed=%s used-as=%s -> %s' % (name, mp.nstr(val, 12), printed, direction, 'OK' if ok else 'ROUNDED WRONG WAY'))
# Lemma 3.4 constants
L = 53632*mp.sqrt(2)*pi**4/729
HC5 = 76545*mp.sqrt(2)/107264
b1 = (1/pi)*L**(mp.mpf(1)/30)*HC5**(mp.mpf(1)/8)
b2 = pi**(-mp.mpf(7)/8)*(1/L)**(mp.mpf(1)/30)
K3 = mp.mpf(1800)/2911*b1*mp.sqrt(pi)
K4 = mp.mpf(28800)/54481*b2*mp.sqrt(pi)
rep('0.4750 = (1800/2911) b1 sqrt(pi) [C1]', K3, mp.mpf('0.4750'), 'upper')
rep('0.2531 = (28800/54481) b2 sqrt(pi) [C2]', K4, mp.mpf('0.2531'), 'upper')
rep('12.496 = (8/3)*4.686 [E1]', mp.mpf(8)/3*mp.mpf('4.686'), mp.mpf('12.496'), 'upper')
rep('4.465 = 0.893*5 [E2]', mp.mpf('0.893')*5, mp.mpf('4.465'), 'upper')
# x-function bounds (l. 700-705): P(x)/(x+2)^{7/2} on (0,1/4]
P = lambda x: (8*x**3 + 36*x**2 + 60*x + 35)/(x + 2)**mp.mpf(3.5)
xs = [mp.mpf(i)/4000 for i in range(0, 1001)]
vals = [P(x) for x in xs]
mono = all(vals[i+1] < vals[i] for i in range(len(vals)-1))
out.append('P(x)/(x+2)^3.5 on [0,1/4]: max %s (35/(8sqrt2)=%s) min %s (6704/2187=%s) decreasing=%s' % (mp.nstr(max(vals),12), mp.nstr(35/(8*mp.sqrt(2)),12), mp.nstr(min(vals),12), mp.nstr(mp.mpf(6704)/2187,12), mono))
dlog = [x*mp.diff(P, x)/P(x) for x in xs[1:]]
out.append('x P1/P on (0,1/4]: min %s max %s (phi5 monotone iff < 9/2)' % (mp.nstr(min(dlog),6), mp.nstr(max(dlog),6)))
# lambda_5 constant (3/(2 sqrt pi))*(6704/2187)*(2pi)^{9/2} == 53632 sqrt2 pi^4/729
out.append('lambda5 const check: %s vs %s' % (mp.nstr(3/(2*mp.sqrt(pi))*mp.mpf(6704)/2187*(2*pi)**mp.mpf(4.5), 20), mp.nstr(L, 20)))
out.append('h5 const check: %s vs %s' % (mp.nstr((35/(8*mp.sqrt(2)))/(mp.mpf(6704)/2187), 20), mp.nstr(HC5, 20)))
# log(2+x) <= log 2 + 0.893 x^{1/3}
f = lambda x: (mp.log(2 + x) - mp.log(2))/mp.cbrt(x)
xm = mp.findroot(lambda x: mp.diff(f, x), 31)
rep('0.893 >= max (log(2+x)-log2)/x^(1/3)', f(xm), mp.mpf('0.893'), 'upper')
out.append('   argmax x = %s' % mp.nstr(xm, 8))
# Kershner constant: sup over intervals |int_u1^u2 e(x^2/2) dx| (f''=1); 1.343 claimed sup for |f''|>=1
F = lambda u: mp.quad(lambda x: mp.expjpi(x*x), [0, u])   # int_0^u e(x^2/2) dx = int exp(i pi x^2)
# extremal for symmetric quadratic: maximize |F(u2)-F(u1)|; by symmetry F odd -> u1=-u2 gives 2|F(u)|; also try asymmetric
best = (0, None)
for i in range(1, 400):
    u = mp.mpf(i)/200
    v = 2*abs(F(u))
    if v > best[0]: best = (v, u)
u0 = best[1]
g = lambda u: -2*abs(F(u))
um = mp.findroot(lambda u: mp.diff(g, u), u0)
kq = 2*abs(F(um))
rep('1.343 >= sup |int e(x^2/2)| symmetric', kq, mp.mpf('1.343'), 'upper')
out.append('   at u = %s (Kershner gamma0/sqrt(2pi): 3.3643/sqrt(2pi) = %s)' % (mp.nstr(um, 10), mp.nstr(mp.mpf('3.3643')/mp.sqrt(2*pi), 10)))
# asymmetric scan
asym = 0
for i in range(-300, 301, 6):
    for j in range(i+6, 601, 6):
        v = abs(F(mp.mpf(j)/200) - F(mp.mpf(i)/200))
        if v > asym: asym = v; aw = (i/200, j/200)
out.append('   asymmetric grid max %s at %s' % (mp.nstr(asym, 10), aw))
# Riemann-Siegel C0 maximum (Hiary/Gabcke): beta = max |cos(pi/2 (z^2+3/4))/cos(pi z)|, z in [0,1]
C0 = lambda z: mp.cos(pi/2*(z*z + mp.mpf(3)/4))/mp.cos(pi*z)
zs = [mp.mpf(i)/2000 for i in range(0, 2001) if i != 1000]
bmax = max(abs(C0(z)) for z in zs)
out.append('RS beta = max|C0| = %s (cos(pi/8)=%s); beta*(2pi)^(1/4) = %s (printed 1.48)' % (mp.nstr(bmax, 12), mp.nstr(mp.cos(pi/8), 12), mp.nstr(bmax*(2*pi)**0.25, 12)))
# Lemma 2.1 constants: (10 log 2 + 1 + gamma + 9/32)/pi + 1 <= 4 ; gamma+2 <= 1/2+3log2
out.append('Lem2.1 constant: %s (<= 4)' % mp.nstr((10*mp.log(2) + 1 + mp.euler + mp.mpf(9)/32)/pi + 1, 10))
out.append('Lem2.1 step gamma+2=%s <= 1/2+3log2=%s' % (mp.nstr(mp.euler + 2, 10), mp.nstr(mp.mpf(1)/2 + 3*mp.log(2), 10)))
out.append('psi(3/2)+gamma=%s <= 2log2=%s' % (mp.nstr(mp.digamma(1.5) + mp.euler, 10), mp.nstr(2*mp.log(2), 10)))
# harmonic bound (Montgomery-Vaughan form) sum_{n<=x} 1/n <= log x + gamma + 9/32 for x >= 2 : scan
worst = -1
for i in range(200, 20001):
    x = mp.mpf(i)/100
    s = mp.harmonic(mp.floor(x))
    d = s - mp.log(x) - mp.euler
    if d > worst: worst = d; xw = x
out.append('max_{2<=x<=200} H(floor x) - log x - gamma = %s at x=%s (9/32 = 0.28125)' % (mp.nstr(worst, 10), xw))
# Lemma 2.3: S3 bound true form: (2/pi)(psi(N+1/2)-psi(1/2)) with N<beta-alpha  vs (2/pi)log(beta-alpha)+2
w = -100
for i in range(101, 100001, 7):
    D = mp.mpf(i)/100   # beta-alpha > 1
    N = mp.ceil(D) - 1
    lhs = 2/pi*(mp.digamma(N + mp.mpf(1)/2) - mp.digamma(mp.mpf(1)/2))
    rhs = 2/pi*mp.log(D) + 2
    if lhs - rhs > w: w = lhs - rhs; Dw = D
out.append('Lem2.3 S3: grid max (b-a = 1.01 + 0.07k <= 1000) of [true bound - (2/pi)log(b-a)-2] = %s at %s (must be <0; exact sup -0.72676 at b-a -> 1+, s3_examples.py)' % (mp.nstr(w, 8), mp.nstr(Dw, 6)))
# Yang Lemma 2.5 edge-case constant: 2^{25/3} (32/(3 sqrt(5 pi)))^2 > 2336
out.append('Yang 2.5: 2^(25/3)*(32/(3sqrt(5pi)))^2 = %s (> 2336 needed)' % mp.nstr(2**(mp.mpf(25)/3)*(32/(3*mp.sqrt(5*pi)))**2, 12))
# Yang Lemma 2.2 lambda_0 and A2
l0 = 1 + 4*(2 - mp.sqrt(4 + pi))/pi
out.append('Yang 2.2 lambda0=%s check (4/sqrtpi+sqrt l0) sqrt l0 = %s' % (mp.nstr(l0, 10), mp.nstr((4/mp.sqrt(pi) + mp.sqrt(l0))*mp.sqrt(l0), 12)))
print('\n'.join(out))

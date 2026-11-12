# Falsification tests for Lemma 5.9 (lem:chi, TeX l. 1252-1273) and Lemma 5.10 (lem:res, l. 1275-1295) of Chourasiya-Simonic.
import mpmath as mp
mp.mp.dps = 25
pi = mp.pi
def chi(z): return (2*pi)**z*mp.rgamma(z)/(2*mp.cos(pi*z/2))
worst = {}
def rec(k, v, info):
    if k not in worst or v > worst[k][0]: worst[k] = (v, info)
ys = [mp.mpf(2) + j*mp.mpf('0.05') for j in range(0, 400)] + [mp.mpf(22) + j for j in range(0, 200)] + [mp.mpf(10)**(mp.mpf(k)/4) for k in range(10, 17)]
xs = [mp.mpf(-1) + j*mp.mpf('0.05') for j in range(0, 41)]
y0 = mp.mpf(2)
for x in xs:
    for y in ys:
        z = mp.mpc(x, y)
        q = abs(mp.gamma(z))/(mp.sqrt(2*pi)*abs(z)**(x - mp.mpf(1)/2)*mp.exp(-pi*y/2))
        rec('(eq:Gamma) upper: q/exp(1/(2y0^2)), y0=2', q/mp.exp(1/(2*y0**2)), (float(x), float(y)))
        rec('(eq:Gamma) lower: exp(-1/(2y0^2))/q, y0=2', mp.exp(-1/(2*y0**2))/q, (float(x), float(y)))
        c = (abs(z)/(2*pi))**(x - mp.mpf(1)/2)*abs(chi(z))
        Up = 1/(1 - mp.exp(-pi*y0))*mp.exp(1/(2*y0**2)); Um = 1/(1 + mp.exp(-pi*y0))*mp.exp(-1/(2*y0**2))
        rec('(eq:chiBound) upper c/U+(2)', c/Up, (float(x), float(y)))
        rec('(eq:chiBound) lower U-(2)/c', Um/c, (float(x), float(y)))
# |chi(z)| <= 0.4|y|^{1/2-x} + 16 for x in [-1/2, 0], all real y
for j in range(0, 11):
    x = mp.mpf(-1)/2 + j*mp.mpf('0.05')
    for y in [mp.mpf(k)/10 for k in range(0, 1200)] + [mp.mpf(120) + 7*k for k in range(0, 300)] + [mp.mpf(10)**k for k in range(4, 9)]:
        z = mp.mpc(x, y); v = abs(chi(z))/(mp.mpf('0.4')*abs(y)**(mp.mpf(1)/2 - x) + 16)
        rec('|chi| <= 0.4|y|^{1/2-x}+16, x in [-1/2,0]', v, (float(x), float(y)))
# |chi(z)| <= 6.31/(1-x), x in [1/2,1), all real y; and |(z-1)chi(z)| <= 6.31 on [1/2,1]x[-2pi,2pi]; |chi| <= 1.1 for |y| >= 2pi
for j in range(0, 51):
    x = mp.mpf(1)/2 + j*mp.mpf('0.01')
    for y in [mp.mpf(k)/50 for k in range(0, 400)] + [mp.mpf(8) + k for k in range(0, 300)] + [mp.mpf(10)**k for k in range(3, 8)]:
        z = mp.mpc(x, y)
        if x == 1 and y == 0: continue
        cz = chi(z)
        if x < 1: rec('|chi| <= 6.31/(1-x), x in [1/2,1)', abs(cz)*(1 - x)/mp.mpf('6.31'), (float(x), float(y)))
        if y <= 2*pi: rec('|(z-1)chi(z)| <= 6.31, |y| <= 2pi', abs((z - 1)*cz)/mp.mpf('6.31'), (float(x), float(y)))
        if y >= 2*pi: rec('|chi| <= 1.1, |y| >= 2pi', abs(cz)/mp.mpf('1.1'), (float(x), float(y)))
    for y in [2*pi]:
        z = mp.mpc(x, y); rec('|(z-1)chi(z)| <= 6.31, |y| <= 2pi', abs((z - 1)*chi(z))/mp.mpf('6.31'), (float(x), float(y)))
# Lemma 5.10: |Res| <= 12 sqrt(T) e^{-pi T/2} log T for t >= T >= 5
for T in [mp.mpf(5), mp.mpf(6), mp.mpf(10), mp.mpf(30), mp.mpf(100)]:
    for dt in [0, mp.mpf('0.1'), 1, 5, 20]:
        t = T + dt
        res = (2*mp.euler + mp.digamma(mp.mpf(1)/2 - 1j*t) + mp.log(T))*T**(mp.mpf(1)/2 - 1j*t)*mp.gamma(mp.mpf(1)/2 - 1j*t)
        rec('Lemma 5.10 |Res|/(12 sqrt(T) e^{-pi T/2} log T)', abs(res)/(12*mp.sqrt(T)*mp.exp(-pi*T/2)*mp.log(T)), (float(T), float(t)))
for k, (v, info) in sorted(worst.items()):
    print('%-48s max = %.6f at %s%s' % (k, float(v), info, '  <-- COUNTEREXAMPLE' if v > 1 else ''))

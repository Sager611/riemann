# T3: EP sensitivity of the approximate functional equation main term at Lehmer's pair 7005
import mpmath as mp
mp.mp.dps = 30
N = 33
def theta(t): return (mp.loggamma(mp.mpf(1)/4 + 0.5j*t) - mp.loggamma(mp.mpf(1)/4 - 0.5j*t))/(2j) - t*mp.log(mp.pi)/2
def Zm(t):
    th = theta(t)
    return 2*mp.fsum(mp.cos(th - t*mp.log(n))/mp.sqrt(n) for n in range(1, N+1))
g1, g2 = mp.mpf('7005.0628661749205814'), mp.mpf('7005.1005646726467216')
tm = mp.findroot(lambda t: mp.siegelz(t, derivative=1), (g1+g2)/2)
Zt, Z2 = mp.siegelz(tm), mp.siegelz(tm, derivative=2)
R = Zt - mp.re(Zm(tm))
print('N=%d; extremum t_m=%s  Z(t_m)=%s  Z\'\'(t_m)=%s' % (N, mp.nstr(tm,12), mp.nstr(Zt,6), mp.nstr(Z2,6)))
print('EP distance |Z(t_m)| = %s vs |Z\'\'|delta^2/8 = %s ; RS remainder R(t_m) = %s' % (mp.nstr(abs(Zt),5), mp.nstr(abs(Z2)*(g2-g1)**2/8,5), mp.nstr(R,5)))
# main term's own extremum and quadratic (EP normal form) prediction
ts = mp.findroot(lambda t: mp.re(mp.diff(Zm, t)), tm)
z0, z2 = mp.re(Zm(ts)), mp.re(mp.diff(Zm, ts, 2))
pred = mp.sqrt(2*z0/z2) if z0*z2 > 0 else None
print('Z_main extremum t*=%s value %s, Z_main\'\'=%s' % (mp.nstr(ts,12), mp.nstr(z0,6), mp.nstr(z2,6)))
root = mp.findroot(Zm, ts + 1j*(pred if pred else mp.mpf('0.01')))
print('complex zero of Z_main: %s ; Puiseux prediction t* +- i*%s' % (mp.nstr(root,12), mp.nstr(pred,6) if pred else 'none'))
# linear (non-EP) sensitivity at a well-separated zero, for contrast
g3 = mp.im(mp.zetazero(6712))
r3 = mp.findroot(lambda t: mp.re(Zm(t)), g3)
Rg = mp.siegelz(g3) - mp.re(Zm(g3))
print('isolated zero 6712: gamma=%s, Z_main zero %s, shift %s, linear prediction -R/Z\' = %s' % (mp.nstr(g3,12), mp.nstr(r3,12), mp.nstr(r3-g3,5), mp.nstr(Rg/mp.siegelz(g3,derivative=1),5)))

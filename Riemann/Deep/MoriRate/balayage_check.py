
import mpmath as mp
mp.mp.dps = 30
def check(y0, x1):
    w0 = (y0/2)/(1 + mp.sqrt(1 - y0**2/4))
    P = lambda th: (1 - w0**2)/abs(mp.expj(th) - w0)**2/(2*mp.pi)
    # x = 2/cos(theta), theta in (-pi, pi); log|x1 - x| integrated against the Poisson kernel (both slit sides)
    f = lambda th: mp.log(abs(x1 - 2/mp.cos(th)))*P(th)
    U = mp.quad(f, [-mp.pi, -mp.pi/2, -mp.acos(2/abs(x1)) if abs(x1) > 2 else -1, 0, mp.acos(2/abs(x1)) if abs(x1) > 2 else 1, mp.pi/2, mp.pi])
    return U, mp.log(abs(x1 - y0))
for y0 in (0.5, 1.3):
    for x1 in (2.5, 3.0, -4.0, 10.0):
        U, L = check(mp.mpf(y0), mp.mpf(x1))
        print(y0, x1, mp.nstr(U, 15), mp.nstr(L, 15), mp.nstr(U - L, 3))
# Green function at 0 vs g2
for y0 in (0.5, 1.3):
    w0 = (y0/2)/(1 + mp.sqrt(1 - y0**2/4))
    print('g(0,y0)', mp.nstr(-mp.log(w0), 15), 'g2', mp.nstr(mp.log((2 + mp.sqrt(4 - y0**2))/y0), 15))

import os
BASE = os.path.abspath(os.environ.get('RAMAUDIT', '.scratch/ramaudit')) + '/'
import mpmath as mp, time, json
mp.mp.dps = 60
def epsR(L):
    L = mp.mpf(L); x = mp.e**L
    li = lambda y: mp.li(y)
    g = (li(x)**2 - (mp.e*x/L)*li(x/mp.e))/((x/L)**2)
    b = 2*li(x)/(x/L) + L/(L-1)
    if g >= 0: return None, g
    e = (-b + mp.sqrt(b*b - 4*g))/2
    return e, g
out = {}
for L in [58, 103, 500, 872, 877, 878, 907, 1000, 2000, 3158, 3158.442, 3159, 3160]:
    e, g = epsR(L)
    out[L] = (mp.nstr(e, 12), mp.nstr(3*mp.mpf(L)**4*e, 8) if e else None, mp.nstr(g, 8))
    print(L, out[L], flush=True)


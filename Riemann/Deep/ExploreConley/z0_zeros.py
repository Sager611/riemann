
import time, mpmath as mp, json
mp.mp.dps = 30
t0=time.time()
zs=[]
for n in range(6700, 6721):
    g = mp.zetazero(n).imag
    zs.append((n, g))
print('zetazero time', round(time.time()-t0,2))
for i,(n,g) in enumerate(zs):
    gap = (zs[i+1][1]-g) if i+1 < len(zs) else None
    sp = 2*mp.pi/mp.log(g/(2*mp.pi))
    print(n, mp.nstr(g,15), '' if gap is None else mp.nstr(gap,8), '' if gap is None else mp.nstr(gap/sp,6))
t0=time.time()
s = mp.mpf(0.5)+1j*mp.mpf(7005.08)
for k in range(5):
    v = mp.zeta(s + 0.01j*k + 0.05)
print('zeta off-line x5 time', round(time.time()-t0,3), mp.nstr(v,8))
t0=time.time()
for k in range(5):
    v = mp.siegelz(mp.mpf(7005.08)+0.01*k)
print('siegelz x5 time', round(time.time()-t0,3))
json.dump([[n, str(g)] for n,g in zs], open('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/zeros_6700_6720.json','w'))

import mpmath as mp
mp.mp.dps = 30
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
a = {1: mp.mpf(1), 2: kap, 3: -kap, 4: mp.mpf(-1), 0: mp.mpf(0)}
def phidh(u):
    u = mp.mpf(u); return 2*mp.e**(1.5*u)*mp.fsum(a[n % 5]*n*mp.e**(-mp.pi*n*n*mp.e**(2*u)/5) for n in range(1, 400))
us = [mp.mpf(k)/100 for k in range(-250, 251)]
vals = [phidh(u) for u in us]
print('kappa', mp.nstr(kap, 12))
print('even check max |Phi(u) - Phi(-u)| / max:', mp.nstr(max(abs(phidh(u) - phidh(-u)) for u in us[::25])/max(vals), 5))
print('min/max Phi_DH on [-2.5, 2.5]:', mp.nstr(min(vals)/max(vals), 6), 'at u =', mp.nstr(us[vals.index(min(vals))], 4))
print('sign changes:', sum(1 for i in range(len(vals) - 1) if vals[i]*vals[i+1] < 0))

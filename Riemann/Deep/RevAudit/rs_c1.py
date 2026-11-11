# Max of the Riemann-Siegel coefficients |C0(p)|, |C1(p)| on [0,1] (mpmath, 30 digits), to pin the scale of Gabcke's 0.127.
# C0(p) = Psi(p) = cos(2 pi (p^2 - p - 1/16)) / cos(2 pi p);  C1(p) = -Psi'''(p) / (96 pi^2)  (Gabcke 1979 / Haselgrove).
import mpmath as mp
mp.mp.dps = 30
def Psi(p):
    return mp.cos(2*mp.pi*(p*p - p - mp.mpf(1)/16)) / mp.cos(2*mp.pi*p)
def C1(p):
    return -mp.diff(Psi, p, 3) / (96*mp.pi**2)
best0 = max((abs(Psi(mp.mpf(k)/2000 + mp.mpf('1e-9'))), k) for k in range(0, 2001) if abs(mp.cos(2*mp.pi*(mp.mpf(k)/2000)))>1e-6)
best1 = max((abs(C1(mp.mpf(k)/400 + mp.mpf('1e-7'))), k) for k in range(0, 401) if abs(mp.cos(2*mp.pi*(mp.mpf(k)/400)))>1e-3)
print('max|C0| ~', mp.nstr(best0[0], 10), 'at p =', best0[1]/2000, '; cos(pi/8) =', mp.nstr(mp.cos(mp.pi/8), 10))
print('max|C1| ~', mp.nstr(best1[0], 10), 'at p =', best1[1]/400)
print('(2pi)^(1/4) cos(pi/8) =', mp.nstr((2*mp.pi)**0.25*mp.cos(mp.pi/8), 10))
print('(2pi)^(3/4) max|C1| =', mp.nstr((2*mp.pi)**0.75*best1[0], 10), ' (bound 0.127 t^-3/4 needs this <= 0.127 asymptotically)')

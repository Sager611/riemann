import numpy as np
from scipy.optimize import curve_fit
Bd = np.array([20,30,40,50,60,70,80,100.])
F = np.array([0.02896023,0.02802508,0.02773068,0.02755454,0.02747949,0.02743165,0.02739549,0.02735516])
n = np.array([118,253,443,700,994,1348,1772,2741.])
for Bpe in [0.0271148, 0.0270648, 0.0271648]:
    E = F - Bpe
    print('B+eps =', Bpe, 'E =', np.array2string(E, precision=6))
    for lo in [2, 3]:
        x = Bd[lo:]; y = F[lo:]
        f3 = lambda b, Einf, c, p: Einf + c * b ** (-p)
        try:
            p3, cv = curve_fit(f3, x, y, p0=[0.0273, 0.05, 1.0], maxfev=20000)
            print(f'  fit Bd>={x[0]:.0f}: Finf={p3[0]:.7f} (+-{np.sqrt(cv[0,0]):.1e}) c={p3[1]:.4f} p={p3[2]:.3f}  -> Finf-(B+eps)={p3[0]-Bpe:+.2e}  resid max {np.abs(f3(x,*p3)-y).max():.1e}')
        except Exception as e:
            print('  fit fail', e)
    # fixed Einf=0 power law on E
    m = np.polyfit(np.log(Bd[3:]), np.log(E[3:]), 1)
    print(f'  pure power E ~ Bd^{m[0]:.3f}, resid {np.abs(np.polyval(m, np.log(Bd[3:])) - np.log(E[3:])).max():.3f}')
# local exponents
F0 = 0.0271148
E = F - F0
for i in range(1, len(Bd)):
    print(f'local p {Bd[i-1]:.0f}->{Bd[i]:.0f}: {-np.log(E[i]/E[i-1])/np.log(Bd[i]/Bd[i-1]):.3f}   in n: {-np.log(E[i]/E[i-1])/np.log(n[i]/n[i-1]):.3f}')
# Richardson-type: assume F = Finf + c/Bd^p with p from consecutive triples
for i in range(2, len(Bd)):
    b0, b1, b2 = Bd[i-2:i+1]; f0, f1, f2 = F[i-2:i+1]
    # solve for p numerically
    from scipy.optimize import brentq
    g = lambda p: (f0 - f1) / (f1 - f2) - (b0 ** -p - b1 ** -p) / (b1 ** -p - b2 ** -p)
    try:
        p = brentq(g, 0.05, 6)
        c = (f0 - f1) / (b0 ** -p - b1 ** -p); Finf = f2 - c * b2 ** -p
        print(f'triple {b0:.0f},{b1:.0f},{b2:.0f}: p={p:.3f} Finf={Finf:.7f} Finf-(B+eps)={Finf-F0:+.2e}')
    except Exception as e:
        print('triple fail', b0, b1, b2, e)

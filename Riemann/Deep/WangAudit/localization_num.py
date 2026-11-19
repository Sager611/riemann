# Numerical check of lem:localization with actual zeros (delta = 0): D(x) = 2 pi F_I(x) - int_I |A(x,t)|^2 dt  vs  x L^3,
# and of the prop:pointwise main term: (2 pi / H) F_I(x) vs L^2/x^2 + log x  (rough at T ~ 600).
import sys, numpy as np
G = np.array([float(l.split()[1]) for l in open(sys.argv[1])])
Gall = np.concatenate([G, -G])
T, th = 600.0, 0.95
H = T**th; L = np.log(T)
g = G[(G > T) & (G <= T+H)]
d = g[:, None] - g[None, :]
w = 4/(4+d**2)
ts = np.arange(T, T+H, 0.01) + 0.005
print(f'T={T} theta={th} H={H:.1f} N(I)={len(g)} Zmax={G[-1]:.1f} grid {len(ts)}')
for x in (1.0, 3.0, 10.0, 40.0):
    F = (np.cos(d*np.log(x))*w).sum()
    acc = 0.0
    for k in range(0, len(ts), 400):
        t = ts[k:k+400]
        A = (2*np.exp(1j*np.log(x)*(Gall[None, :]-t[:, None]))/(1+(t[:, None]-Gall[None, :])**2)).sum(axis=1)
        acc += (np.abs(A)**2).sum()*0.01
    Dl = 2*np.pi*F - acc
    print(f'x={x:5.1f}: 2piF_I={2*np.pi*F:10.2f} int_I|A|^2={acc:10.2f} diff={Dl:9.2f} diff/(x L^3)={Dl/(x*L**3):7.4f} '
          f'| (2pi/H)F_I={2*np.pi*F/H:7.3f} vs L^2/x^2+log x={L*L/x/x+np.log(x):7.3f} (with log(T/2pi): {np.log(T/2/np.pi)**2/x/x+np.log(x):7.3f})')

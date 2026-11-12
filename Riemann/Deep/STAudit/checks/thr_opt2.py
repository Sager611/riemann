
import numpy as np
# exponents (powers of X, O(delta) dropped) for m = 2 route of STLemP Prop F2, with free y (log_X Y), v (log_X V)
def classI(eta, y, v):
    nus1 = np.linspace(v, 1.0, 400)            # blocks V/2 < N <= X: Lemma L (p0=1) + Theorem 1' with theta=1/2
    E1 = 0.5 + np.maximum(1.0, 0.5 + v + 0.5*np.maximum(0, nus1-0.5)) + nus1*(2*eta-1)
    nus2 = np.linspace(1.0, max(y,1.0), 400)   # blocks X < N < 2Z: Corollary L (lift)
    E2 = np.maximum(nus2/2, v + nus2/4) + 2*nus2*eta
    return max(E1.max(), E2.max())
def best(eta):
    b = 9; arg = None
    ys = 15/(8*(1+eta))
    for y in np.linspace(ys-0.03, ys+0.03, 121):
        c2 = 15/8 - y*(0.5-eta)
        if c2 >= b: continue
        for v in np.linspace(0.0, 0.5, 201):
            if v > y: continue
            E = max(c2, classI(eta, y, v))
            if E < b: b, arg = E, (y, v)
    return b, arg
etas = np.array([1/14, 0.1, 0.15, 0.17, 0.1752, 0.18, 0.2, 0.24])
res = [(eta,) + best(eta) for eta in etas]
Areq = max((E-1)/eta for eta,E,_ in res if eta >= 1/14 - 1e-9)
i = int(np.argmax([(E-1)/eta if eta >= 1/14-1e-9 else -9 for eta,E,_ in res]))
print('grid threshold max_eta (E*(eta)-1)/eta over eta in [1/14, 0.26]:', round(Areq,4), 'at eta', round(res[i][0],4), 'y,v =', res[i][2])
print('A2* = (23-3*sqrt5)/8 =', (23-3*5**0.5)/8)
for eta,E,arg in res: print(f' eta={eta:.3f} E*={E:.4f} (y,v)=({arg[0]:.3f},{arg[1]:.3f})  formula F2={15*(1+4*eta)/(16*(1+eta)):.4f}')

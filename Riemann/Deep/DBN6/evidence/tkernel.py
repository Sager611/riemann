# Feasibility scan (non-rigorous): translated derivative kernel for Raymond's probes.
# 2K_y(d,b) >= lam*y*J(d,b) - mu*y*dJ(d,b), h=1, probes at d+-r, y<=rho, 0<=b<=y.
import numpy as np
def K(h, d, b):
    return (h-b)/(d*d+(h-b)**2) + (h+b)/(d*d+(h+b)**2)
def dK(h, d, b):
    s = 0
    for c in (h-b, h+b):
        s = s + (d*d - c*c)/(d*d + c*c)**2
    return s
dg = np.concatenate([np.linspace(0, 3, 601), np.linspace(3, 60, 300), [1e3, 1e4]])
def frontier(r, rho, mus):
    ys = rho*np.linspace(0.02, 1, 50)
    out = []
    LHS=[]; Jv=[]; dJv=[]
    for y in ys:
        for bb in np.linspace(0, 1, 41):
            b = bb*y
            d = dg.copy()
            if bb == 1: d = d[d > 1e-9]
            lhs = 2*K(y, d, b)/y
            J = 0.5*(K(1, d+r, b) + K(1, d-r, b))
            dJ = 0.5*(dK(1, d+r, b) + dK(1, d-r, b))
            LHS.append(lhs); Jv.append(J); dJv.append(dJ)
    LHS = np.concatenate(LHS); Jv = np.concatenate(Jv); dJv = np.concatenate(dJv)
    for mu in mus:
        lam = np.min((LHS + mu*dJv)/Jv)
        out.append((mu, lam))
    return out
mus = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
for r, rho in [(0.575, 0.4), (0.6, 0.4), (0.6, 0.45), (0.6, 0.3), (0.6, 0.2), (0.0, 0.4), (0.3, 0.4), (0.8, 0.4)]:
    fr = frontier(r, rho, mus)
    print('r=%.3f rho=%.2f ' % (r, rho) + ' '.join('mu=%.1f:lam=%.4f' % (m, l) for m, l in fr))

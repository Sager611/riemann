# uc_t4b.py - is (lambda(mu)-lambda(0))/mu ever > 1/(1-rho)? (DBN6 lever kernel table, same formulas as Riemann/Deep/DBN6/evidence/rmodel_lib.py)
import numpy as np
def K(h, d, b): return (h-b)/(d*d+(h-b)**2) + (h+b)/(d*d+(h+b)**2)
def dK(h, d, b): return sum((d*d - c*c)/(d*d + c*c)**2 for c in (h-b, h+b))
dg = np.concatenate([np.linspace(0, 3, 401), np.linspace(3, 60, 200), [1e3, 1e4]])
def lam_table(r, rho, mus):
    LHS=[]; Jv=[]; dJv=[]
    for y in rho*np.linspace(0.02, 1, 40):
        for bb in np.linspace(0, 1, 31):
            b = bb*y; d = dg if bb < 1 else dg[dg > 1e-9]
            LHS.append(2*K(y, d, b)/y); Jv.append(0.5*(K(1, d+r, b)+K(1, d-r, b))); dJv.append(0.5*(dK(1, d+r, b)+dK(1, d-r, b)))
    LHS = np.concatenate(LHS); Jv = np.concatenate(Jv); dJv = np.concatenate(dJv)
    return [np.min((LHS + mu*dJv)/Jv) for mu in mus]
mus = [0, 0.05, 0.1, 0.3, 0.7]
print("r rho lambda(mu=0,.05,.1,.3,.7)  max_mu (lam(mu)-lam(0))/mu  vs Poisson break-even 1/(1-rho)")
for r in [0.0, 0.25, 0.5, 1.0, 1.5]:
    for rho in [0.15, 0.3, 0.45]:
        L = lam_table(r, rho, mus)
        q = max((L[i]-L[0])/mus[i] for i in range(1, len(mus)))
        print(r, rho, " ".join("%.4f" % x for x in L), " %.4f  %.4f %s" % (q, 1/(1-rho), "GAIN POSSIBLE" if q > 1/(1-rho) else "no"))

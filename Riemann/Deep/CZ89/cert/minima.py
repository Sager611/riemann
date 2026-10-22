
import numpy as np, mpmath as mp
import explore as E
mp.mp.dps = 30
b = 1/mp.sqrt(2); a = [mp.mpf(2)/25, mp.mpf(1)/25, -mp.mpf(1)/200]
def s(z): return mp.mpf(1) if z == 0 else mp.sin(z)/z
def K(x):
    x = mp.mpf(x)
    return b/(2*mp.sin(b))*(s(mp.pi*x-b)+s(mp.pi*x+b)) + sum(a[k-1]/2*(s(mp.pi*(x-k))+s(mp.pi*(x+k))) for k in (1,2,3))
Wd = {}
rows = {0:['0.306','0.645','1','1','1','2'],1:['0.347','0.047','0','0','1'],2:['0.347','0.616','0','1'],3:['0.347','0.047','1'],4:['0.347','0.645'],5:['0.306']}
for i,v in rows.items():
    for t,val in enumerate(v): Wd[(i,i+1+t)] = mp.mpf(val)
def S(g):
    y = [mp.mpf(0)]
    for gi in g: y.append(y[-1]+mp.mpf(gi))
    return sum(Wd[(i,j)]*K(y[j]-y[i])**2 for (i,j) in Wd) + sum(mp.mpf(x) for x in g)/200
def grad_hess(g, h=1e-5):
    g = np.array(g, dtype=float)
    f0 = E.S_exact(g)[0]
    G = np.zeros(6); Hm = np.zeros((6,6))
    for k in range(6):
        e = np.zeros(6); e[k] = h
        G[k] = (E.S_exact(g+e)[0]-E.S_exact(g-e)[0])/(2*h)
        for l in range(6):
            f = np.zeros(6); f[l] = h
            Hm[k,l] = (E.S_exact(g+e+f)[0]-E.S_exact(g+e-f)[0]-E.S_exact(g-e+f)[0]+E.S_exact(g-e-f)[0])/(4*h*h)
    return f0, G, Hm
starts = [[1.0837,1.8939,1.0414,1.0249,1.0289,1.0256],[1.0426,1.0469,1.9031,1.0568,1.0297,1.0274],[1.0227,1.0232,1.0184,1.0184,1.0232,1.0227],[1.0447,1.0507,1.912,1.0809,1.9032,1.0855],[1.0847,1.904,1.0464,1.0464,1.904,1.0847],[1.9018,1.0749,1.0349,1.0227,1.0259,1.0214]]
for st in starts:
    g = np.array(st)
    for it in range(30):
        f0, G, Hm = grad_hess(g)
        step = np.linalg.solve(Hm, G)
        g = g - step
        if np.abs(step).max() < 1e-12: break
    f0, G, Hm = grad_hess(g)
    ev = np.linalg.eigvalsh(Hm)
    print("g*", np.array2string(g, precision=10), "sum %.6f" % g.sum())
    print("  S float %.12f  S mp %s  margin %.3e  |grad| %.1e" % (f0, mp.nstr(S(list(g)), 15), f0-0.044, np.abs(G).max()))
    print("  hess eig", np.array2string(ev, precision=4))
gstar = np.array(starts[0])


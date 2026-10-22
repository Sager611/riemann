# Non-rigorous model of Hakimoglu-Brown's derivative chain (arXiv 2610.12803) as a function of the wall X.
# Reproduces their C_i, G_i, lambda_i, mu_i, k_i and final time, then explores higher X and lower final height b.
import numpy as np, mpmath as mp, sys, json
mp.mp.dps = 30
p = 61/20.; B = p*p; sig_p = (1+p)/2
N1 = 65536
n = np.arange(2, N1+1, dtype=np.float64); ln = np.log(n)
# von Mangoldt for psi sums done by mpmath zeta derivatives instead
def zeta_terms(u):
    Z = float(mp.zeta(u)); z1 = float(mp.zeta(u, derivative=1)); z2 = float(mp.zeta(u, derivative=2))
    psi = -z1/Z
    psi2 = z2/Z - (z1/Z)**2   # (zeta'/zeta)' = zeta''/zeta - (zeta'/zeta)^2 >0
    return Z, psi, psi2
def tail_int(j, v, a):
    # int_a^inf (log y)^j y^{-v} dy
    return float(mp.quad(lambda y: mp.log(y)**j * y**(-v), [a, mp.inf]))
def CG(X, tm, tp, eta=(1e-6, 4e-5, 16e-4)):
    L = float(mp.log(X/(4*mp.pi)))
    Nst = int(mp.floor(mp.sqrt(X/(4*mp.pi))))
    u = sig_p + (L/4)*tm; u8 = sig_p + (L/8)*tm
    e = (p-1)/2 - (L/8)*tm
    w = (np.exp(tp*ln*ln/4) - 1) * n**(-u)
    d = [float(np.sum(w*ln**j)) + tail_int(j, u8, N1) for j in range(3)]
    r0 = 1.01*Nst**(-p)*(Nst**(e+1)/(e+1) + Nst**e)
    Ls = 2*np.log(Nst)+3
    r = [r0, Ls*r0, (Ls*Ls+1e-9)*r0]
    Z, psi, psi2 = zeta_terms(u)
    c1 = 0.5+1e-10; c2 = 0.25+1e-10
    e0 = d[0]+r[0]+eta[0]; e1 = c1*d[1]+r[1]+eta[1]; e2 = c2*d[2]+r[2]+eta[2]+1e-9
    E0 = Z*e0; E1 = Z*(e1+c1*psi*e0); E2 = Z*(e2+2*c1*psi*e1+c2*(psi2+psi*psi)*e0)
    C = L/4 - 1e-6 - c1*psi - E1/(1-E0)
    G = c2*psi2 + E2/(1-E0) + (E1/(1-E0))**2 + 2e-10
    return C, G, dict(L4=L/4, Z=Z, psi=psi, psi2=psi2, d=d, r=r, eps0=E0)
# kernel LP: maximize lam*a - mu*bcoef subject to Phi(D,r)>=0 on grid
zz = np.concatenate([np.linspace(0,0.999,1500), 1-np.logspace(-3,-9,200)])
Dg = zz/(1-zz)
def lp_kernel(A2, a, bcoef, rgrid=60):
    rs = np.linspace(0, A2, rgrid)
    D, R = np.meshgrid(Dg, rs, indexing='ij'); D = D.ravel(); R = R.ravel()
    w = D + B - R; Dl = w*w + 4*D*R
    f = (D + A2 + 3*R)
    cl = f*2*w*Dl                       # coefficient of lam
    cm = -f*((2*w+4*B)*Dl - 8*B*w*w)    # coefficient of mu (N = 2 lam w Dl - mu((2w+4B)Dl - 8 B w^2))
    rhs = 2*Dl*Dl
    # constraint: cl*lam + cm*mu <= rhs (cl>0). For each mu, lam_max = min (rhs - cm*mu)/cl
    best = (-1e9, 0, 0)
    for mu in np.linspace(0, 2.5, 1251):
        lam = np.min((rhs - cm*mu)/cl)
        val = lam*a - mu*bcoef
        if val > best[0]: best = (val, lam, mu)
    return best
def k_of(C, G, q, lam, mu):
    return 4*(lam*C/p - mu*G - 2*lam/(B-q) - 2*mu*(B+q)/(B-q)**2)
def chain(X, b, dt=1/200., q0=9/16, tmax=0.2, verbose=False, fixed=None):
    t = 0.0; q = q0; rows = []
    i = 0
    while q > b + 1e-15:
        if fixed and i < len(fixed):
            C, G, lam, mu = fixed[i]
        else:
            C, G, info = CG(X, t, t+dt)
            a = C/p - 2/(B-q); bc = G + 2*(B+q)/(B-q)**2
            val, lam, mu = lp_kernel(q, a, bc)
            lam *= (1-1e-9); mu = mu
        k = k_of(C, G, q, lam, mu)
        qn = (q + 2/k)*np.exp(-k*dt) - 2/k
        if qn <= b:
            d = np.log((2+k*q)/(2+k*b))/k
            rows.append((i, t, q, C, G, lam, mu, k, d)); t += d; q = b; break
        rows.append((i, t, q, C, G, lam, mu, k, dt))
        t += dt; q = qn; i += 1
        if t > tmax: return None, rows
    if verbose:
        for r in rows: print('%2d t=%.4f q=%.6f C=%.5f G=%.5f lam=%.4f mu=%.4f k=%.4f d=%.5f' % r)
    return t + b/2, rows
if __name__ == '__main__':
    X = 60000083951.5
    for t in [0, 0.06, 0.115]:
        C, G, info = CG(X, t, t+1/200)
        print('X=6e10 t=%.3f C=%.5f G=%.5f eps0=%.4g' % (t, C, G, info['eps0']))
    val = chain(X, 1/25, verbose=True)[0]
    print('model Lambda bound at HB X, b=1/25:', val)

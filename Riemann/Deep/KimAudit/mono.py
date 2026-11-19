#!/usr/bin/env python3
"""kimaudit T5: own implementation of the DERIVATION Sec. 8-9 majorants (from the text, not from certify.py);
dense grid u in [12, 300] (+ larger u): monotone directions of Sec. 9 and FINAL > 0. Not a proof (floats at 40 digits).
Run: bash .scratch/memcap.sh 290 600 uv run --no-project --with mpmath python3 mono.py [a]"""
import sys
import mpmath as mp
mp.mp.dps = 40
mpf, pi, e, L, sq = mp.mpf, mp.pi, mp.exp, mp.log, mp.sqrt
a = mpf(sys.argv[1]) if len(sys.argv) > 1 else mpf('0.37338196661320752')
th, thU, dR, De, K, ee = mpf(9)/10, mpf(49)/50, mpf(1), mpf(17)/25, mpf(4), mpf(1)/1000
beta = 3/(2*pi); rho = sq(thU/th); Ch = th*(2-th)/(4*(1-th)**2)
eps0 = mp.nsum(lambda k: k**2*(2*pi*k**2-3)/(2*pi-3)*e(-pi*(k**2-1)), [2, mp.inf])
def maj(u):
    n = u*(pi*e(2*u) - mpf(9)/4)
    M1 = 4 - 2*beta*u*e(-2*u)/(1-beta) - u*(u+1)*a**2/(2*n)                    # 8.1 (Q1)*u
    kl = 4*(n/u)*e(-2*(De+2*ee)) - a**2/(2*n); eta = sq(K/kl)                   # 8.3
    Phi = 2*a*(1+Ch)*eta*(sq(th) + a*(u+1)/n)                                   # 8.3
    E = a**2*rho**4*(th*(u+1)**2/n + a**2*(u+1)**4/n**3)/(8*(1-thU)**2)         # 8.4
    y = (rho-1)*u + rho; F = (e(2*y)-1-2*y)/u                                    # 8.5
    ep = 1/n + 2*eta; uL = u - De
    dh1 = De**2/(uL*u) + (2/u)*(1+9*u/(4*n))*(De-(1-e(-2*De))/2)
    dh2 = ep**2/u**2 + (1/u)*(1+9*u/(4*n))*(e(2*ep)-1-2*ep)
    G = thU/4 + a*rho**2*(u+1)*(sq(th)+a*(u+1)/n)/(2*n) + dh1 + dh2 + th/4 + th**2/(32*(1-th))
    c6 = F - G
    Hcp = (n/u)*(2*(e(2*rho)-1) - rho*th/2) - a*rho*(sq(th)+a)
    lr = L(1+eps0) - L(1-beta*e(-2*uL)) + L(5*n/K)/2 - c6*n
    fin = 2*mp.cos(Phi) - (1+eps0)*e(E)/(2*K) - e(lr)
    hcc = 4*n/u - (th*n/(2*(u+1)**2) + a*(sq(th)/(u+1) + a/n))
    return dict(n=n, M1=M1, eta=eta, M3a=ee-eta, M3b=dR-1/n-2*eta, Phi=Phi, E=E, F=F, G=G, c6=c6, Hcp=Hcp, lr=lr, fin=fin, hcc=hcc, c6m=c6-1/(2*n))
up_ = ['M1', 'M3a', 'M3b', 'F', 'c6', 'Hcp', 'fin', 'hcc']; dn_ = ['eta', 'Phi', 'E', 'G', 'lr']
us = [mpf(12) + mpf(j)/20 for j in range(0, 5761)] + [mpf(x) for x in (400, 600, 1000, 3000, 10000)]
prev = None; viol = {}; mins = {}
for u in us:
    m = maj(u)
    for k in ('M1', 'M3a', 'M3b', 'c6', 'Hcp', 'fin', 'hcc', 'c6m'):
        mins[k] = min(mins.get(k, mpf(10)**99), m[k] if k != 'Hcp' else m[k]-1)
    if prev is not None:
        for k in up_:
            if m[k] < prev[k]: viol[k] = viol.get(k, 0) + 1
        for k in dn_:
            if m[k] > prev[k]: viol[k] = viol.get(k, 0) + 1
    prev = m
m12 = maj(mpf(12))
print('a =', a, ' eps0 =', mp.nstr(eps0, 12))
print('u=12:', ', '.join('%s=%s' % (k, mp.nstr(m12[k], 10)) for k in ('n', 'M1', 'eta', 'Phi', 'E', 'F', 'G', 'c6', 'Hcp', 'lr', 'fin', 'hcc')))
print('Q2 slack =', mp.nstr(1-e(-2*De) - th/4 - th**2/(16*(1-th)), 12), ' rho =', mp.nstr(rho, 12), ' 1/(2(rho-1)) =', mp.nstr(1/(2*(rho-1)), 8))
print('grid points:', len(us), ' monotonicity violations:', viol if viol else 'none')
print('minima over grid (must be > 0; Hcp shown as Hcp-1):', ', '.join('%s=%s' % (k, mp.nstr(v, 6)) for k, v in mins.items()))

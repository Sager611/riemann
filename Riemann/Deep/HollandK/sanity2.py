#!/usr/bin/env python3
# Sanity check 2: (a) Jensen hyperbolicity by Sturm for small d,n; (b) the parameter pipeline at n ~ 1e45 computed
# from the saddle main term at high precision (error terms of the main term are invisible at 4th differences).
import mpmath as mp, sympy as sp, sys
pi = mp.pi
mode = sys.argv[1] if len(sys.argv) > 1 else 'b'
if mode == 'a':
    from flint import fmpz_poly, arb as farb
    mp.mp.dps = 50
    def ureal(n):   # real saddle of 2n log u + 9u/2 - pi e^{2u}: 2n/u + 9/2 = 2 pi e^{2u}  (bisection)
        if n == 0: return mp.mpf('0.3')
        a, b_ = mp.mpf('1e-6'), mp.mpf(50)
        for _ in range(200):
            m = (a + b_)/2
            if 2*n/m + mp.mpf(9)/2 - 2*pi*mp.exp(2*m) > 0: a = m
            else: b_ = m
        return a
    def phi(v, z):
        s = 0
        for k in range(1, 4):
            s += (2*pi**2*k**4*mp.exp(9*v/2) - 3*pi*k**2*mp.exp(5*v/2))*mp.exp(-pi*k**2*mp.exp(2*v))
        return s*mp.exp(2*z*mp.log(v)) if v > 0 else mp.mpf(0)
    def Mt(n):
        U = ureal(n); w = 1/mp.sqrt(8*n + 20)
        pts = sorted(set([mp.mpf(0)] + [max(U + k*w, mp.mpf('1e-9')) for k in (-40, -12, -6, -3, -1.5, 0, 1.5, 3, 6, 12, 40)] + [U + 2, mp.mpf(8)]))
        return mp.quad(lambda v: phi(v, n), pts)
    NS = (0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 300)
    need = sorted(set(n + j for n in NS for j in range(31)))
    gam = {}
    for i, m in enumerate(need):
        gam[m] = mp.factorial(m)/mp.factorial(2*m)*Mt(m)   # = gamma_true(m) / 4^{m+1}; hyperbolicity unaffected
        if i % 40 == 0: print("quad", i, "/", len(need), flush=True)
    bad = 0; cnt = 0
    for n in NS:
        for d in range(2, 31):
            sc = gam[n+d]
            co = [int(mp.nint(mp.binomial(d, j)*gam[n+j]/sc*mp.mpf(10)**45)) for j in range(d+1)]
            P = fmpz_poly(co)
            rts = P.complex_roots()
            nreal = sum(mult for (r, mult) in rts if r.imag.is_zero())   # Arb returns certified real roots with exactly zero imaginary part
            cnt += 1
            if nreal != d: bad += 1; print("NOT hyperbolic?", n, d, nreal, flush=True)
    print("(a) Jensen J^{d,n}, d=2..30, n in", NS, ":", cnt, "polynomials,", bad, "non-hyperbolic (certified real-root isolation, Arb)")
else:
    mp.mp.dps = 600
    def u0(z):
        u = mp.log(z)/2
        for _ in range(400): u = mp.log((z/u + mp.mpf(9)/4)/pi)/2
        return u
    def logmain(z):
        u = u0(z)
        f = 2*z*mp.log(u) + 9*u/2 - pi*mp.exp(2*u)
        fpp = -2*z/u**2 - 4*pi*mp.exp(2*u)
        return mp.log(2*pi**2) + f + mp.log(2*pi/(-fpp))/2
    for Lt in (90, 100, 150, 300):
        a0 = mp.mpf(Lt)/2
        n = mp.floor(a0*(pi*mp.exp(2*a0) - mp.mpf(9)/4)) + 1
        u = u0(n); L = 2*u; eps = 1/L
        lg = [logmain(n+j) - (n+j)*mp.log(4) - mp.loggamma(n+j+mp.mpf(1)/2) for j in range(5)]   # log gamma up to const
        R = [mp.exp(lg[j]-lg[0]) for j in range(5)]
        q = [R[k+1]**2/(R[k]*R[k+2]) for k in range(3)]
        b = n + mp.mpf(1)/2
        delta0 = (lg[2]-2*lg[1]+lg[0]) + 0  # this is -(log q0) + log(1+1/b) etc.; recompute from h:
        h = [logmain(n+j) for j in range(5)]
        d0_ = h[2]-2*h[1]+h[0]; d1_ = h[3]-2*h[2]+h[1]
        X0 = n*L*d0_; Y = n**2*L*(d1_-d0_)
        q1L = 2 - 1/q[0]; Dn = mp.log(q1L/q[1])
        T = (q[0]-1)/(q[1]-1); AJ = 2*q[0]/(T-q[0]); BJ = 1/(q[0]*(1+1/AJ)-1)
        aJ = AJ/(n*L); betaJ = BJ/n; bJ = (BJ-b)*L/n
        RJ4 = R[1]**4*0  # placeholder
        # Jacobi model ratio at j=4: S^4 (A)_4/(A^4 (B)_4), S = B R_1
        S = BJ*R[1]
        RJ = lambda j: S**j*mp.rf(AJ, j)/(AJ**j*mp.rf(BJ, j))
        Qn = -mp.log(R[4]/RJ(4))
        tQ = n**3*L*Qn
        D = BJ/2
        f = lambda k, U: mp.log((U+k)/(U+k+1))
        def eqs(a, bb, c):
            A = a*AJ; B = BJ + bb*n/L; C = D + c*n/L
            return [f(k, A)-f(k, B)+f(k, C)-f(k, D) - mp.log(q[k]) for k in range(3)]
        sol = mp.findroot(lambda a, bb, c: eqs(a, bb, c), (mp.mpf(6)/5, mp.mpf(1)/3, mp.mpf(1)/24), tol=mp.mpf(10)**(-500))
        print("L=%d eps=%s: X0=%s Y=%s n2LDn=%s aJ=%s betaJ=%s bJ=%s tQ=%s sol=(%s,%s,%s)" % (Lt, mp.nstr(eps, 5),
              mp.nstr(X0, 8), mp.nstr(Y, 8), mp.nstr(n**2*L*Dn, 8), mp.nstr(aJ, 8), mp.nstr(betaJ, 10), mp.nstr(bJ, 8), mp.nstr(tQ, 8),
              mp.nstr(sol[0], 8), mp.nstr(sol[1], 8), mp.nstr(sol[2], 8)))

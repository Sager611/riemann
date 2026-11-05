
# ef_check.py: explicit formula for L(s,chi), real primitive chi, checked against zeros and against the lattice identity.
import mpmath as mp, sys, time
mp.mp.dps = 30
CH = {'m4': ([0,1,0,-1], 4, 1), 'm3': ([0,1,-1], 3, 1), 'p5': ([0,1,-1,-1,1], 5, 0), 'p8': ([0,1,0,-1,0,-1,0,1], 8, 0)}
def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else 0
    return 0
def Z(t, tab, q, a):
    s = mp.mpf(1)/2 + 1j*t
    th = t/2*mp.log(q/mp.pi) + mp.im(mp.loggamma((mp.mpf(1)/2 + a + 1j*t)/2))
    v = mp.exp(1j*th)*mp.dirichlet(s, tab)
    return v
def zeros(tab, q, a, T, h=0.05):
    zs = []; t0 = mp.mpf('0.01'); z0 = mp.re(Z(t0, tab, q, a)); t = t0
    while t < T:
        t1 = t + h; z1 = mp.re(Z(t1, tab, q, a))
        if mp.sign(z1) != mp.sign(z0):
            r = mp.findroot(lambda u: mp.re(Z(u, tab, q, a)), (t, t1), solver='anderson')
            zs.append(r)
        t, z0 = t1, z1
    return zs
def run(key, A, m, T):
    tab, q, a = CH[key]
    t0 = time.time()
    zs = zeros(tab, q, a, T)
    imax = max(abs(mp.im(Z(g, tab, q, a))) for g in zs[:5])
    nu = mp.mpf(m) + mp.mpf(1)/2; fm = mp.factorial(m)
    Fe = lambda t: A*mp.sqrt(mp.pi)*fm*2**nu*(A*t)**(-nu)*mp.besselj(nu, A*t) if t != 0 else A*mp.sqrt(mp.pi)*fm*mp.gamma(nu+1)**-1*mp.mpf(2)**0*  mp.mpf(1)*mp.beta(mp.mpf(1)/2, m+1)/ (A*mp.sqrt(mp.pi)*fm) * A*mp.sqrt(mp.pi)*fm / A * A
    def Fe(t):
        if abs(t) < mp.mpf('1e-20'): return A*mp.beta(mp.mpf(1)/2, m+1)
        return A*mp.sqrt(mp.pi)*fm*2**nu*(A*t)**(-nu)*mp.besselj(nu, A*t)
    def Fo2(t):   # |F|^2 for f = (x/A)(1-(x/A)^2)^m
        if abs(t) < mp.mpf('1e-20'): return mp.mpf(0)
        return (A*mp.sqrt(mp.pi)*fm*2**nu*(A*t)**(-nu)*mp.besselj(nu+1, A*t))**2
    fe = lambda x: (1-(x/A)**2)**m if abs(x) <= A else mp.mpf(0)
    fo = lambda x: (x/A)*(1-(x/A)**2)**m if abs(x) <= A else mp.mpf(0)
    phi = lambda t: mp.log(q/mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + mp.mpf(a)/2 + 1j*t/2))
    out = {}
    for sec, H2, f in (('even', lambda t: Fe(t)**2, fe), ('odd', Fo2, fo)):
        zs_sum = 2*sum(H2(g) for g in zs)
        arch = mp.quad(lambda t: phi(t)*H2(t), mp.linspace(0, 4*T, 81))/mp.pi
        c = mp.e**(2*A); pr = mp.mpf(0)
        for n in range(2, int(mp.floor(c))+1):
            if n >= c: break
            w = vm(n); x = tab[n % q]
            if w == 0 or x == 0: continue
            s = mp.log(n)
            g = mp.quad(lambda y: f(y)*f(y-s), [s-A, (s-A+A)/2, A])
            pr += x*w/mp.sqrt(n)*2*g
        Qform = arch - pr
        # lattice: L = 2A, r_k = pi k / L, mu_k
        L = 2*A; cc = mp.e**L
        def mu(k):
            r = mp.pi*k/L
            tail = sum(cc**(-(2*n + mp.mpf(1)/2 + a))*(2*n + mp.mpf(1)/2 + a)/((2*n + mp.mpf(1)/2 + a)**2 + r*r) for n in range(40))
            ps = 0
            for n in range(2, int(mp.floor(cc))+1):
                if n >= cc: break
                w = vm(n); x = tab[n % q]
                if w and x: ps += x*w/mp.sqrt(n)*mp.cos(r*mp.log(n))
            return phi(r) + 2*(-1)**k*tail - 2*ps
        K = 3000
        lat = (mu(0)*H2(0) + 2*mp.fsum(mu(k)*H2(mp.pi*k/L) for k in range(1, K)))/(2*L)
        out[sec] = (zs_sum, Qform, lat)
        print(f"{key} A={A} m={m} T={T} zeros={len(zs)} first={mp.nstr(zs[0],12)} sector={sec}: zero-sum={mp.nstr(zs_sum,18)} formula={mp.nstr(Qform,18)} lattice={mp.nstr(lat,18)} |zs-form|={mp.nstr(abs(zs_sum-Qform),3)} |lat-form|={mp.nstr(abs(lat-Qform),3)} imZmax={mp.nstr(imax,2)}", flush=True)
    print(f"  time {time.time()-t0:.0f}s", flush=True)
if __name__ == '__main__':
    run(sys.argv[1], mp.mpf(sys.argv[2]), int(sys.argv[3]), mp.mpf(sys.argv[4]))

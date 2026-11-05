# ef_check_zeta.py: the Weil explicit formula for zeta (DERIVATION.md sec. 1-2), checked against a direct sum over the zeros
# (mpmath.zetazero) and against the lattice identity (sec. 3).  Test functions f = (1 - (x/A)^2)^m (even) and
# (x/A)(1 - (x/A)^2)^m (odd) on [-A, A] (C^{m-1}, so the explicit formula applies).  Usage: ef_check_zeta.py A m T
import mpmath as mp, sys, time
mp.mp.dps = 30

def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else 0
    return 0

def run(A, m, T):
    t0 = time.time()
    zs = []; n = 1
    while True:
        g = mp.im(mp.zetazero(n))
        if g > T: break
        zs.append(g); n += 1
    nu = mp.mpf(m) + mp.mpf(1)/2; fm = mp.factorial(m)
    def Fe(t):     # int (1-(x/A)^2)^m cos(tx) dx
        if abs(t) < mp.mpf('1e-20'): return A*mp.beta(mp.mpf(1)/2, m+1)
        return A*mp.sqrt(mp.pi)*fm*2**nu*(A*t)**(-nu)*mp.besselj(nu, A*t)
    def Fo2(t):    # |int (x/A)(1-(x/A)^2)^m e^{itx} dx|^2
        if abs(t) < mp.mpf('1e-20'): return mp.mpf(0)
        return (A*mp.sqrt(mp.pi)*fm*2**nu*(A*t)**(-nu)*mp.besselj(nu+1, A*t))**2
    fe = lambda x: (1-(x/A)**2)**m if abs(x) <= A else mp.mpf(0)
    fo = lambda x: (x/A)*(1-(x/A)**2)**m if abs(x) <= A else mp.mpf(0)
    phi = lambda t: -mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2))
    c = mp.e**(2*A); L = 2*A
    for sec, H2, f in (('even', lambda t: Fe(t)**2, fe), ('odd', Fo2, fo)):
        # sum over zeros rho = 1/2 + i gamma, both signs of gamma; |F(-gamma)| = |F(gamma)| for real f
        zs_sum = 2*mp.fsum(H2(g) for g in zs)
        # expected omitted part: 2 int_T^inf (1/2 pi) log(t/(2 pi)) |F|^2 dt  (Riemann-von Mangoldt density)
        tail = mp.quad(lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)*H2(t), [T, 2*T, 8*T, mp.inf])*2
        # pole terms h(i/2) + h(-i/2) = 2 Re[F(i/2) conj F(-i/2)]: +2 (int f cosh(x/2))^2 (even f), -2 (int f sinh(x/2))^2 (odd f)
        if sec == 'even': pole = 2*mp.quad(lambda x: f(x)*mp.cosh(x/2), [-A, 0, A])**2
        else: pole = -2*mp.quad(lambda x: f(x)*mp.sinh(x/2), [-A, 0, A])**2
        arch = mp.quad(lambda t: phi(t)*H2(t), mp.linspace(0, 4*T, 81) + [mp.inf])/mp.pi
        pr = mp.mpf(0)
        for n in range(2, int(mp.floor(c)) + 1):
            if n >= c: break
            w = vm(n)
            if w == 0: continue
            s = mp.log(n)
            g = mp.quad(lambda y: f(y)*f(y - s), [s - A, s/2, A])
            pr += w/mp.sqrt(n)*2*g
        Qform = pole + arch - pr
        def mu(k):
            r = mp.pi*k/L
            tl = mp.fsum(c**(-(2*j + mp.mpf(1)/2))*(2*j + mp.mpf(1)/2)/((2*j + mp.mpf(1)/2)**2 + r*r) for j in range(40))
            ps = 0
            for nn in range(2, int(mp.floor(c)) + 1):
                if nn >= c: break
                w = vm(nn)
                if w: ps += w/mp.sqrt(nn)*mp.cos(r*mp.log(nn))
            return phi(r) + 2*(-1)**k*tl - 2*ps
        K = 3000
        lat = pole + (mu(0)*H2(0) + 2*mp.fsum(mu(k)*H2(mp.pi*k/L) for k in range(1, K)))/(2*L)
        print(f"zeta A={A} m={m} T={T} zeros={len(zs)} first={mp.nstr(zs[0], 12)} sector={sec}: zero-sum={mp.nstr(zs_sum, 20)} "
              f"formula={mp.nstr(Qform, 20)} (pole {mp.nstr(pole, 8)}, arch {mp.nstr(arch, 8)}, primes {mp.nstr(pr, 8)}) "
              f"|zs-form|={mp.nstr(abs(zs_sum - Qform), 3)} expected tail={mp.nstr(tail, 3)} |zs+tail-form|={mp.nstr(abs(zs_sum + tail - Qform), 3)} "
              f"lattice={mp.nstr(lat, 20)} |lat-form|={mp.nstr(abs(lat - Qform), 3)}", flush=True)
    print(f"  time {time.time()-t0:.0f}s", flush=True)

if __name__ == '__main__':
    run(mp.mpf(sys.argv[1]), int(sys.argv[2]), mp.mpf(sys.argv[3]))

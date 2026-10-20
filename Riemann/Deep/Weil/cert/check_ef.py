
import mpmath as mp
mp.mp.dps = 30
def vm(n):
    # von Mangoldt
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else mp.mpf(0)
    return mp.mpf(0)
for s in [mp.mpf('0.25'), mp.mpf('0.15')]:
    g = lambda x: mp.e**(-x**2/(2*s**2))
    h = lambda r: s*mp.sqrt(2*mp.pi)*mp.e**(-s**2*r**2/2)
    pole = 2*s*mp.sqrt(2*mp.pi)*mp.e**(s**2/8)
    archF = mp.quad(lambda r: h(r)*mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2)), [-mp.inf,0,mp.inf])/(2*mp.pi) - g(0)*mp.log(mp.pi)
    archX = -g(0)*(mp.euler+mp.log(mp.pi)) + mp.quad(lambda x: (2*g(0)*mp.e**(-2*x) - mp.e**(-x/2)*(g(x)+g(-x)))/(1-mp.e**(-2*x)), [0,1,mp.inf])
    N = int(mp.e**(s*mp.sqrt(2*80)))+2
    primes = 2*sum(vm(n)/mp.sqrt(n)*g(mp.log(n)) for n in range(2, N))
    rhs = pole + archF - primes
    zs = 0; k = 1
    while True:
        gm = mp.im(mp.zetazero(k)); t = 2*h(gm); zs += t; k += 1
        if t < mp.mpf(10)**-25: break
    print('s',s,'archF-archX',mp.nstr(archF-archX,5),'rhs',mp.nstr(rhs,20),'zeros',mp.nstr(zs,20),'diff',mp.nstr(rhs-zs,5),'nz',k-1)

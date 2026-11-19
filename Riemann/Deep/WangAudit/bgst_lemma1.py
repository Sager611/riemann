# Check BGST Lemma 1 (Wang (eq:explicit)) with actual zeros: the EXACT identity behind it,
# sum_rho 2x^{i(g-t)}/(1+(t-g)^2) = x R(3/2+it) - x^{-1} R(-1/2+it),
# R(s) = sum_{n<=x} Lambda(n) n^{-s} + zeta'/zeta(s) - x^{1-s}/(1-s) - sum_k x^{-2k-s}/(2k+s)   (Landau; x not a prime power),
# and the size of E_x(t) = A + D_x - log(t+2)/x  (claimed << 1/x + x^{1/2}/t^2 + x^{-5/2}/t).
import sys, mpmath as mp
mp.mp.dps = 20
G = [mp.mpf(l.split()[1]) for l in open(sys.argv[1])]
Zmax = G[-1]
def vonmangoldt_upto(N):
    lam = [0.0]*(N+1)
    sieve = list(range(N+1))
    for p in range(2, N+1):
        if sieve[p] == p:
            for q in range(p*p, N+1, p):
                if sieve[q] == q: sieve[q] = p
    for n in range(2, N+1):
        p = sieve[n]; m = n
        while m % p == 0: m //= p
        if m == 1: lam[n] = mp.log(p)
    return lam
NL = 200000
lam = vonmangoldt_upto(NL)
def A(x, t):  # zeros +-gamma, delta = 0 (all computed zeros are on the line)
    s = 0
    for g in G:
        for gg in (g, -g):
            s += 2*mp.power(x, 1j*(gg-t))/(1+(t-gg)**2)
    return s
def tailbound(t):  # crude bound for the omitted zeros |gamma| > Zmax: 2 sum 1/(gamma-t)^2 ~ (1/pi) log(Z/2pi)/(Z-t)
    return 2*(mp.log(Zmax/(2*mp.pi))/(2*mp.pi))*(1/(Zmax-t) + 1/(Zmax+t))
def D(x, t):  # sum_{n>=2} Lambda(n) n^{-1/2-it} min(n/x, x/n); tail n > NL by PNT integral
    s = 0
    for n in range(2, NL+1):
        if lam[n]: s += lam[n]*mp.power(n, -0.5-1j*t)*min(n/x, x/n)
    s += x*mp.power(NL, -0.5-1j*t)/(0.5+1j*t)   # int_NL^inf x u^{-3/2-it} du
    return s
def R(s, x):
    out = 0
    for n in range(2, int(x)+1):
        if lam[n]: out += lam[n]*mp.power(n, -s)
    out += mp.zeta(s, derivative=1)/mp.zeta(s) - mp.power(x, 1-s)/(1-s)
    out -= mp.nsum(lambda k: mp.power(x, -2*k-s)/(2*k+s), [1, mp.inf])
    return out
print('zeros', len(G), 'Zmax', mp.nstr(Zmax,8))
for x in (mp.mpf('1.5'), mp.mpf('6.5'), mp.mpf('30.5')):
    for t in (mp.mpf(300), mp.mpf(800)):
        a = A(x, t); exact = x*R(mp.mpf(1.5)+1j*t, x) - R(mp.mpf(-0.5)+1j*t, x)/x
        dd = D(x, t)
        e = a + dd - mp.log(t+2)/x            # E_x from the zero side
        ep = exact + dd - mp.log(t+2)/x       # E_x from the prime side (Landau)
        print(f'x={mp.nstr(x,4)} t={int(t)}: |A - (xR(3/2+it) - R(-1/2+it)/x)| = {mp.nstr(abs(a-exact),3)} (omitted-zero bound {mp.nstr(tailbound(t),3)}); '
              f'x|E_x| zeros-side = {mp.nstr(x*abs(e),4)}, primes-side = {mp.nstr(x*abs(ep),4)}; |A| = {mp.nstr(abs(a),4)}')

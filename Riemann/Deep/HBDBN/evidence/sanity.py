
# Sanity: exact t=0 values S(x,p,0) = (1/2) Re (xi'/xi)(w), dS/dp = (1/4) Re (xi'/xi)'(w), w=(1+p)/2 - i x/2
import mpmath as mp
mp.mp.dps = 30
N = 10**6
sieve = bytearray([1])*(N+1); sieve[0]=sieve[1]=0
for i in range(2, int(N**0.5)+1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
pp = []
for l in range(2, N+1):
    if sieve[l]:
        n = l
        while n <= N:
            pp.append((n, mp.log(l))); n *= l
print("prime powers:", len(pp))
X = mp.mpf('60000083951.5'); p = mp.mpf(61)/20
for dx in (0, mp.mpf('0.37'), mp.mpf('12.9'), mp.mpf('1000.25')):
    x = X + dx
    w = (1+p)/2 - 1j*x/2
    rho = mp.mpc(0); rhop = mp.mpc(0)
    for (n, ll) in pp:
        ln = mp.log(n); e = mp.exp(-w*ln)
        rho += ll*e; rhop += ll*ln*e
    zl = -rho          # zeta'/zeta(w)
    zlp = rhop         # (zeta'/zeta)'(w)
    xil = 1/w + 1/(w-1) - mp.log(mp.pi)/2 + mp.digamma(w/2)/2 + zl
    xilp = -1/w**2 - 1/(w-1)**2 + mp.psi(1, w/2)/4 + zlp
    S = mp.re(xil)/2; dS = mp.re(xilp)/4
    print("x=X+%s: S(x,61/20,0)=%s  dS/dp=%s  L/4=%s" % (mp.nstr(dx,8), mp.nstr(S,10), mp.nstr(dS,10), mp.nstr(mp.log(x/(4*mp.pi))/4,10)))
print("paper row 0: C_0 = 5.2902, G_0 = 0.2340; tail of truncation < 1e-6")


# Bare diagonal Q(e_k) of the collar ground vector (paper's "target diagonal" B_k = 1/2 + log(1/w_k)).
import mpmath as mp, sys
mp.mp.dps = 20
def run(k):
    a = mp.log(k+1)/2; ak = mp.log(k)/2
    nrm = mp.sqrt(2*(k**-0.5 - (k+1)**-0.5))
    def e(y):
        ay = abs(y)
        return (mp.sign(y)*mp.e**(-ay/2)/nrm) if (ak < ay < a) else 0
    def g(t):
        # int e(y) e(y-t) dy ; both supports: y in +-(ak,a), y-t in +-(ak,a)
        tot = 0
        for (lo1, hi1) in ((ak, a), (-a, -ak)):
            for (lo2, hi2) in ((ak, a), (-a, -ak)):
                lo = max(lo1, lo2 + t); hi = min(hi1, hi2 + t)
                if hi > lo: tot += mp.quad(lambda y: e(y)*e(y-t), [lo, hi])
        return tot
    g0 = g(0)
    kinks = sorted(set([mp.mpf(0), a-ak, 2*ak, a+ak, 2*a]))
    arch = g0*(-(mp.euler + mp.log(mp.pi) + mp.log(1-mp.e**(-4*a)))) + mp.quad(lambda t: (2*mp.e**(-2*t)*g0 - 2*g(t)*mp.e**(-t/2))/(1-mp.e**(-2*t)), kinks)
    A = mp.quad(lambda y: e(y)*mp.e**(y/2), [-a, -ak]) + mp.quad(lambda y: e(y)*mp.e**(y/2), [ak, a])
    pole = -2*A**2
    pr = 0
    for n in range(2, k+1):
        p = next(q for q in range(2, n+1) if n % q == 0); m = n
        while m % p == 0: m //= p
        if m == 1: pr += mp.log(p)/mp.sqrt(n)*2*g(mp.log(n))
    Q = arch + pole - pr
    Bk = mp.mpf(1)/2 + mp.log(1/mp.log(1+mp.mpf(1)/k))
    print(f"k={k}: ||e||^2={mp.nstr(g0,10)} Q(e_k)={mp.nstr(Q,10)} (arch {mp.nstr(arch,8)}, pole {mp.nstr(pole,8)}, primes {mp.nstr(-pr,8)}); paper B_k={mp.nstr(Bk,10)}, 1/2+log k={mp.nstr(0.5+mp.log(k),8)}", flush=True)
for k in sys.argv[1:]: run(int(k))

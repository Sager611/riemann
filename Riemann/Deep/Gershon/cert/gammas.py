
import mpmath as mp, sys, time
mp.mp.dps = 60
pi = mp.pi
def Phi(u):
    s = mp.mpf(0)
    for n in range(1, 8):
        s += (2*pi**2*n**4*mp.e**(9*u) - 3*pi*n**2*mp.e**(5*u))*mp.e**(-pi*n**2*mp.e**(4*u))
    return s
K = int(sys.argv[1]) if len(sys.argv)>1 else 160
t0=time.time()
lg = []
for k in range(K+1):
    f = lambda u: Phi(u)*u**(2*k)
    # peak location roughly
    val = mp.quad(f, [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3])
    lg.append(mp.log(val) - mp.loggamma(2*k+1))
print("time", time.time()-t0, file=sys.stderr)
with open("gammas_log.txt","w") as fo:
    for k,v in enumerate(lg): fo.write(f"{k} {mp.nstr(v,50)}\n")

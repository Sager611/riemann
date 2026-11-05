# /// script
# dependencies = ["mpmath", "numpy"]
# ///
# nblp2 checks (single-threaded):
# (1) exact N=10 L^p norms via periodicity (period q*lcm(1..10), digamma weights) for chi_{-3}, chi_{-4}, chi0 mod 3, chi0 mod 4;
# (2) A4(b) residue identity J = 2pi[R(2s')/(2s') - delta c_q/((2s'-1)L(1+eps))] on Re s = s' = 1/p - eps/2, incl. a complex chi mod 5;
#     and B(eps) = (1/2pi) int |1-R|^2/|s|^2 (truncated) vs K_{p,chi} eps;
# (3) illustration ||1 + H_eps||_{1.1} for chi_{-3} (H_eps = -sum_{j<=m} chi(j) phi_eps(j) on I_m), M = 2e6 (tail not bounded).
import mpmath as mp, numpy as np, time, sys
from fractions import Fraction as Fr
t0 = time.time()
PARTS = sys.argv[1] if len(sys.argv) > 1 else '123'
def mobius_upto(n):
    mu=[1]*(n+1); isp=[True]*(n+1)
    for i in range(2,n+1):
        if isp[i]:
            for j in range(2*i,n+1,i): isp[j]=False
            for j in range(i,n+1,i): mu[j]*=-1
            for j in range(i*i,n+1,i*i): mu[j]=0
    mu[0]=0; return mu
mu = mobius_upto(200)
CH = {'chi-3': (3,[0,1,-1],0), 'chi-4': (4,[0,1,0,-1],0), 'chi0mod3': (3,[0,1,1],1), 'chi0mod4': (4,[0,1,0,1],1)}
# ---------- (1) exact N = 10
mp.mp.dps = 25
N = 10; L10 = 2520
for name,(q,vals,dl) in (CH.items() if '1' in PARTS else []):
    chi = lambda n: vals[n % q]
    P = q*L10
    a = {k: Fr(mu[k]*(N-k), N) for k in range(1,N+1)}
    PN1 = sum(a[k]*chi(k)/Fr(k) for k in range(1,N+1))
    v = []; cum = Fr(0); S = 0
    for n in range(1, 2*P+1):
        A = sum(a[k] for k in range(1,N+1) if n % k == 0)
        cum += chi(n)*A; S += chi(n)
        v.append(1 + (dl*S*PN1 if dl else 0) - cum)
    per = all(v[i] == v[i+P] for i in range(P))
    out = []
    for p in [mp.mpf('1.1'), mp.mpf(2)]:
        tot = mp.mpf(0)
        for r in range(1, P+1):
            w = (mp.digamma(mp.mpf(r+1)/P) - mp.digamma(mp.mpf(r)/P))/P
            x = v[r-1]; tot += abs(mp.mpf(x.numerator)/x.denominator)**p * w
        out.append(mp.nstr(tot**(1/p), 12))
    print("(1) %s N=10 period=%d periodic=%s  exact norm_1.1=%s  norm_2=%s  max|v|=%s  [%.0fs]" % (name, P, per, out[0], out[1], mp.nstr(max(abs(float(x)) for x in v),6), time.time()-t0)); sys.stdout.flush()
# ---------- (2) residue identity (Simpson on a grid, step h, |t| <= T; tail ~ 2/T since R averages to 1)
mp.mp.dps = 12
CH2 = dict(CH); CH2['chi5cplx'] = (5,[0,1,1j,-1j,-1],0)
fG = lambda t, q: (8 + 7*mp.log(q) + mp.log(2+abs(t)))*mp.sqrt(2+abs(t))/(mp.mpf(1)/4+t*t)
T = 100; h = mp.mpf('0.05'); n = int(T/h)
def simpson(vals):  # vals at t_j = j*h, j=0..n (n even)
    return h/3*(vals[0] + vals[-1] + 4*sum(vals[1:-1:2]) + 2*sum(vals[2:-1:2]))
for name in (['chi-3','chi0mod3','chi5cplx'] if '2' in PARTS else []):
    q, vals, dl = CH2[name]
    p, eps = mp.mpf(8)/7, mp.mpf('0.1')
    sp = 1/p - eps/2
    Lf = lambda s: mp.dirichlet(s, vals)
    signs = [1] if all(isinstance(x,int) for x in vals) else [1,-1]
    J = 0; B = 0
    for sg in signs:
        Rv = [Lf(mp.mpc(sp, sg*j*h))/Lf(mp.mpc(sp+eps, sg*j*h)) for j in range(n+1)]
        w = [1/(sp**2 + (j*h)**2) for j in range(n+1)]
        Jh = simpson([Rv[j]*w[j] for j in range(n+1)]); Bh = simpson([abs(1-Rv[j])**2*w[j] for j in range(n+1)])
        if len(signs) == 1: J = 2*mp.re(Jh); B = 2*Bh
        else: J += Jh; B += Bh
    J += 2/mp.mpf(T); B = B/(2*mp.pi)
    cq = mp.mpf(len([x for x in vals if x != 0]))/q
    Jf = 2*mp.pi*(Lf(2*sp)/Lf(2*sp+eps)/(2*sp) - (dl*cq/((2*sp-1)*Lf(1+eps)) if dl else 0))
    Jwrong = 2*mp.pi*(-Lf(2*sp)/Lf(2*sp+eps)/(2*sp) - (dl*cq/((2*sp-1)*Lf(1+eps)) if dl else 0))
    CGq = mp.e**4 * q**mp.mpf(3.5) * 2*mp.quad(lambda t: fG(t,q), [0,1,10,100,1e4,1e6,mp.inf])
    u0 = mp.mpf(1)/2 + 1/p
    K = (4*p/(2-p) if dl else 2*abs(mp.zeta(u0, derivative=1))*mp.zeta(u0)) + CGq/(2*mp.pi)
    print("(2) %s p=8/7 eps=0.1 sigma'=%s: J_num=%s  J_formula=%s  |diff|=%.1e  (flipped-sign formula would give %s) | B_num(|t|<=%d)=%s  B/eps=%s  K_bound=%s  [%.0fs]" % (name, mp.nstr(sp,6), mp.nstr(J,8), mp.nstr(Jf,8), float(abs(J-Jf)), mp.nstr(Jwrong,5), T, mp.nstr(B,5), mp.nstr(B/eps,4), mp.nstr(K,5), time.time()-t0)); sys.stdout.flush()
# ---------- (3) illustration, chi_{-3}
M = 2*10**6
sieve = np.ones(M+1, dtype=bool); sieve[:2] = False
for i in range(2, int(M**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]
m = np.arange(1, M+1, dtype=np.float64); w = 1.0/(m*(m+1.0))
r3 = np.arange(1, M+1) % 3; chiv = np.where(r3 == 0, 0.0, np.where(r3 == 1, 1.0, -1.0))
for eps in ([0.4, 0.2, 0.1, 0.05] if '3' in PARTS else []):
    phi = np.ones(M+1)
    for pr in primes: phi[pr::pr] *= (1.0 - float(pr)**(-eps))
    g = 1.0 - np.cumsum(chiv*phi[1:])
    print("(3) chi-3 eps=%.2f ||1+H_eps||_1.1 ~ %.5f  (M=%d, tail est %.1e)  [%.0fs]" % (eps, np.sum(np.abs(g)**1.1*w)**(1/1.1), M, np.mean(np.abs(g[M//2:])**1.1)/M, time.time()-t0)); sys.stdout.flush()
print("done %.0fs" % (time.time()-t0))

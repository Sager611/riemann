# e1phase ls.py — direct test of the load-bearing input (*) of ceilaudit10's route (a) at its binding block (theta=0.6: P=Q=y, A=1, N=N0):
#   S(C) = sum_{p~P} sum_{q~Q,(p,q)=1} | sum_{n~N} C(n) w(n/N) e(n qbar/p) exp(i X0 sqrt((n/N) PQ/(pq))) |   (Bessel phase KEPT)
# exact Kloosterman fractions qbar = q^{-1} mod p; sizes tied as in e1rate at theta=0.6, nu=(2/7)(3theta-1): N = y^{n0/nu}, X0* = y^{(1-theta)/nu}.
# Outputs per X0: l2 large-sieve constant Delta = s_max(F)^2, l1 sums for C = d(n) (divisor), random unimodular, adversarial (Boyd inf->1 ascent).
import sys, math, time
import numpy as np
th = 0.6; nu = 2*(3*th-1)/7
y = int(sys.argv[1]); P = Q = y
N = int(round(y**((2*nu+1-2*th)/nu))); X0s = y**((1-th)/nu)
xs = [float(v) for v in sys.argv[2].split(',')]          # X0 values as multiples of X0*
rng = np.random.default_rng(1)
pr, qr, rr = [], [], []
for p in range(P, 2*P):
    for q in range(Q, 2*Q):
        if math.gcd(p, q) == 1:
            pr.append(p); qr.append(q); rr.append(pow(q, -1, p))
pr = np.array(pr, float); qr = np.array(qr, float); rr = np.array(rr, float)
R = len(pr); n = np.arange(N, 2*N); u = n/N; w = np.sin(np.pi*(u-1))**2
alpha = rr/pr; kap = np.sqrt(P*Q/(pr*qr))
d = np.array([sum(1 + (k*k != m) for k in range(1, int(math.isqrt(m))+1) if m % k == 0) for m in n], float); d /= d.max()
T0 = R*w.sum()
print('y=%d P=Q=%d rows=%d N=%d X0*=%.1f  trivial l1 T0=%.4g  theory LS N+4P^2=%d  DI(A=1) rel gain %.3f' % (
    y, P, R, N, X0s, T0, N+4*P*P, ((P*Q)**0.5 + (N+1)**0.25*(P*Q*(N+1)*(P+1) + N*Q*Q)**0.25)/(P*Q*N)**0.5), flush=True)
def build(X0, al):
    F = np.exp(2j*np.pi*np.outer(al, n))
    if X0 > 0: F *= np.exp(1j*X0*np.outer(kap, np.sqrt(u)))
    F *= w
    return F
def boyd(F, C0, its=40):
    C = C0/np.maximum(np.abs(C0), 1e-300)
    best = 0.0
    for _ in range(its):
        v = F @ C; s = np.abs(v).sum(); best = max(best, s)
        eta = v/np.maximum(np.abs(v), 1e-300)
        g = F.conj().T @ eta; C = g/np.maximum(np.abs(g), 1e-300)
    return best
def measure(X0, al, tag):
    t = time.time(); F = build(X0, al)
    G = F.conj().T @ F; lam = np.linalg.eigvalsh(G)[-1]; del G
    Sd = np.abs(F @ d).sum()
    Sr = np.mean([np.abs(F @ np.exp(2j*np.pi*rng.random(N))).sum() for _ in range(3)])
    starts = [d.astype(complex), np.exp(-1j*X0*np.sqrt(u))*d] + [np.exp(2j*np.pi*rng.random(N)) for _ in range(2)] + [F[i].conj() for i in rng.integers(0, R, 2)]
    Sa = max(boyd(F, c0) for c0 in starts)   # local ascent; starts incl. dechirped d(n) and single-row resonators
    print('%s X0=%9.1f (%.3g X0*)  Delta=%.4g (Delta/(N+4P^2)=%.3f, /N=%.2f)  l1/T0: d(n)=%.4f rand=%.4f adv=%.4f  [%.1fs]' % (
        tag, X0, X0/X0s, lam, lam/(N+4*P*P), lam/N, Sd/T0, Sr/T0, Sa/T0, time.time()-t), flush=True)
    return lam, Sd, Sr, Sa
for m in xs:
    measure(m*X0s, alpha, 'kloost')
if len(sys.argv) > 3:
    al2 = rng.random(R)
    for m in xs: measure(m*X0s, al2, 'randfr')

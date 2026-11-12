# explore44 c8: Prop TL at scale. 2-D W1 between mu_T (zeta(sigma+it), DH(sigma+it), t uniform in [1000, 2000]) and the Bohr lift, by exact
# assignment on n-point samples, with same-law noise floors; then the TL bound min_R e^R W1 + A(R) vs the actual zero sum.
import numpy as np, math, sys, time
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
sys.path.insert(0, '/Users/adriansager/playground/riemann/.scratch/explore44')
from dhlib import hurwitz, f as fdh
t0 = time.time(); rng = np.random.default_rng(4444)
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
P = primes_upto(30000); Pf = P.astype(float); Pbig = primes_upto(10**7).astype(float); Pbig = Pbig[Pbig > 30000]
chi = {1: 1, 2: 1j, 3: -1j, 4: -1, 0: 0}; cp = np.array([chi[int(p) % 5] for p in P])
kk = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1); c1 = (1 - 1j*kk)/2
n = 2500; L2 = {0.6: 2.77591, 0.7: 1.77818}
def lift(sig, m):
    x = Pf**(-sig); s2 = float(np.sum(Pbig**(-2*sig))); Z = []; F = []
    for i in range(0, m, 100):
        z = np.exp(2j*np.pi*rng.random((100, len(P))))
        g = -np.log(1 - x*z).sum(1) + rng.normal(0, np.sqrt(s2/2), 100) + 1j*rng.normal(0, np.sqrt(s2/2), 100)
        Z.append(np.exp(g))
        G1 = -np.log(1 - cp*x*z).sum(1) + rng.normal(0, np.sqrt(s2/2), 100) + 1j*rng.normal(0, np.sqrt(s2/2), 100)
        G2 = -np.log(1 - np.conj(cp)*x*z).sum(1) + rng.normal(0, np.sqrt(s2/2), 100) + 1j*rng.normal(0, np.sqrt(s2/2), 100)
        F.append(c1*np.exp(G1) + np.conj(c1)*np.exp(G2))
    return np.concatenate(Z), np.concatenate(F)
def W1(a, b):
    C = cdist(np.c_[a.real, a.imag], np.c_[b.real, b.imag]); r, c = linear_sum_assignment(C); return float(C[r, c].mean())
for sig in [0.6, 0.7]:
    t = rng.uniform(1000, 2000, 2*n); s = sig + 1j*t
    zt = np.concatenate([hurwitz(s[i:i+500], 1.0, 1060) for i in range(0, 2*n, 500)])
    ft = np.concatenate([fdh(s[i:i+500], M=1060) for i in range(0, 2*n, 500)])
    zl, fl = lift(sig, 2*n)
    wz = W1(zt[:n], zl[:n]); wz0 = W1(zl[:n], zl[n:]); wzt = W1(zt[:n], zt[n:])
    wf = W1(ft[:n], fl[:n]); wf0 = W1(fl[:n], fl[n:]); wft = W1(ft[:n], ft[n:])
    print('sigma', sig, 'W1 zeta: T-vs-lift %.4f  floors lift-lift %.4f T-T %.4f | DH: T-vs-lift %.4f floors %.4f %.4f' % (wz, wz0, wzt, wf, wf0, wft), round(time.time() - t0, 1), flush=True)
    for w, name in [(wz, 'raw'), (max(wz0, wzt), 'floor')]:
        best = min((math.exp(R)*w + math.sqrt(math.pi*L2[sig]/2)*math.erfc(R/math.sqrt(2*L2[sig])), R) for R in np.linspace(0, 6, 601))
        print('   TL bound for zeta (W1 %s %.4f): min_R e^R W1 + A_Herbst(R) = %.4f at R = %.2f ; actual zero sum (RH) = 0' % (name, w, best[0], best[1]))
print('done', round(time.time() - t0, 1))

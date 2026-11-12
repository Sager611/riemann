# explore44 c9: the quantity Prop TL needs, measured: 1-D W1 between the law of log|F(sigma+it)| (t uniform on [T,2T]) and its Bohr lift,
# for T = 250 ... 4000 (zeta and DH, sigma = 0.6), against the same-size lift-vs-lift noise floor.
import numpy as np, sys, time
sys.path.insert(0, '/Users/adriansager/playground/riemann/.scratch/explore44')
from dhlib import hurwitz, f as fdh
t0 = time.time(); rng = np.random.default_rng(99)
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
P = primes_upto(20000); Pf = P.astype(float); Pbig = primes_upto(10**7).astype(float); Pbig = Pbig[Pbig > 20000]
chi = {1: 1, 2: 1j, 3: -1j, 4: -1, 0: 0}; cp = np.array([chi[int(p) % 5] for p in P])
kk = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1); c1 = (1 - 1j*kk)/2
sig = 0.6; x = Pf**(-sig); s2 = float(np.sum(Pbig**(-2*sig)))
NL = 60000; LZ = []; LF = []
for i in range(0, NL, 100):
    z = np.exp(2j*np.pi*rng.random((100, len(P))))
    LZ.append(-np.log(np.abs(1 - x*z)).sum(1) + rng.normal(0, np.sqrt(s2/2), 100))
    G1 = -np.log(1 - cp*x*z).sum(1) + rng.normal(0, np.sqrt(s2/2), 100) + 1j*rng.normal(0, np.sqrt(s2/2), 100)
    G2 = -np.log(1 - np.conj(cp)*x*z).sum(1) + rng.normal(0, np.sqrt(s2/2), 100) + 1j*rng.normal(0, np.sqrt(s2/2), 100)
    LF.append(np.log(np.abs(c1*np.exp(G1) + np.conj(c1)*np.exp(G2))))
LZ = np.concatenate(LZ); LF = np.concatenate(LF)
print('lift samples', NL, 'means zeta %.4f DH %.4f' % (LZ.mean(), LF.mean()), round(time.time() - t0, 1), flush=True)
q = np.linspace(0.0005, 0.9995, 4000)
def w1(a, b): return float(np.mean(np.abs(np.quantile(a, q) - np.quantile(b, q))))
n = 20000
print('floor (lift %d vs lift %d): zeta %.4f DH %.4f' % (n, NL - n, w1(LZ[:n], LZ[n:]), w1(LF[:n], LF[n:])))
for T in [250.0, 500.0, 1000.0, 2000.0, 4000.0]:
    t = rng.uniform(T, 2*T, n); s = sig + 1j*t; M = int(T) + 60
    zt = np.concatenate([hurwitz(s[i:i+256], 1.0, M) for i in range(0, n, 256)])
    ft = np.concatenate([fdh(s[i:i+256], M=M) for i in range(0, n, 256)])
    a, b = np.log(np.abs(zt)), np.log(np.abs(ft))
    print('T=%6g  W1 zeta %.4f (mean %.4f)  DH %.4f (mean %.4f)' % (T, w1(a, LZ), a.mean(), w1(b, LF), b.mean()), round(time.time() - t0, 1), flush=True)
print('done', round(time.time() - t0, 1))

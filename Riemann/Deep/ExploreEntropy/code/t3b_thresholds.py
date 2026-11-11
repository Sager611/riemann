# t3b_thresholds.py (explore32): thresholds of the Beckner-type bounds and of the exact window minima; entropic (Hirschman-Beckner) bounds.
import sys, numpy as np, mpmath as mp
from scipy.optimize import brentq
from scipy.special import digamma
from scipy.integrate import quad
from flint import fmpq, ctx
ctx.prec = 128
sys.argv = ['x']
import importlib.util
spec = importlib.util.spec_from_file_location('t3', 't3_windows.py')
src = open('t3_windows.py').read().split('def K(u)')[0]   # reuse parts() only
exec(src)
D = digamma(0.25) - np.log(np.pi)
K = lambda u: sum(np.exp(-(2*m + 0.5)*u)/(2*m + 0.5) for m in range(int(60/(2*u)) + 2))
aK = brentq(lambda a: D + 2*K(a), 1e-3, 1)
aKo = brentq(lambda a: D + 2*K(a) - 2*(np.sinh(a) - a), 1e-3, 1)
print('killing bound zero: even a=%.6f (c=%.5f); odd a=%.6f (c=%.5f); e^{-gamma}/(2 pi) = %.6f' % (aK, np.exp(2*aK), aKo, np.exp(2*aKo), np.exp(-np.euler_gamma)/(2*np.pi)))
N = 80
def lm(c, sec, which):
    cq = fmpq(int(round(c*10**6)), 10**6)
    A, P, R = parts(cq, N, sec)
    M = {'arch': A, 'arch_pole': A + P, 'Q': A + P + R}[which]
    return float(np.linalg.eigvalsh(M)[0])
for sec, which, lo, hi in [('even', 'arch', 1.15, 1.35), ('odd', 'arch', 2.0, 2.2), ('even', 'arch_pole', 2.2, 2.5), ('odd', 'arch_pole', 2.05, 2.2)]:
    cz = brentq(lambda c: lm(c, sec, which), lo, hi, xtol=1e-6)
    print('%s %s: zero at c=%.6f a=%.6f (N=%d)' % (sec, which, cz, np.log(cz)/2, N), flush=True)
for c in [2.0]:
    print('c=2 margins (N=%d): even arch+pole %.6g, odd arch+pole %.6g; n=2 weight ln2/sqrt2 = %.6f' % (N, lm(c, 'even', 'arch_pole'), lm(c, 'odd', 'arch_pole'), np.log(2)/np.sqrt(2)))
# Hirschman-Beckner (continuous): f supported in [-a, a] => H(|f~|^2) >= ln(e/2) - ln(2a); Arch = int phi(2 pi xi) |f~|^2 dxi;
# Gibbs: Arch/||f||^2 >= sup_{beta>1} (h0 - ln Z(beta))/beta, Z(beta) = int exp(-beta phi(2 pi xi)) dxi
phi = lambda t: -np.log(np.pi) + float(mp.re(mp.digamma(0.25 + 0.5j*t)))
def lnZ(beta):
    g = lambda xi: np.exp(-beta*phi(2*np.pi*xi))
    v = 2*(quad(g, 0, 1, limit=200)[0] + quad(g, 1, 50, limit=400)[0] + quad(g, 50, np.inf, limit=400)[0])
    return np.log(v)
betas = [1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2, 2.5, 3, 4, 6]
LZ = {b: lnZ(b) for b in betas}
def HB(a):
    h0 = np.log(np.e/2) - np.log(2*a)
    return max((h0 - LZ[b])/b for b in betas)
for a in [0.001, 0.005, 0.01, 0.02, 0.05, 0.0894, 0.1733, 0.3466]:
    print('Hirschman-Beckner Gibbs bound a=%.4f: Arch/||f||^2 >= %.4f   (killing bound %.4f)' % (a, HB(a), D + 2*K(a)))
aH = brentq(HB, 1e-4, 0.3)
print('Hirschman-Beckner bound zero: a=%.6f (c=%.6f)' % (aH, np.exp(2*aH)))
# Lattice (Maassen-Uffink/Hirschman on the circle of length 2L, f supported on half): H(q) >= ln 2;
# sum_k mu(k) q_k >= sup_beta (ln 2 - ln sum_k e^{-beta mu(k)})/beta  (Chekmachev's lattice symbol mu, c <= 2: no primes)
def mu(k, c):
    L = np.log(c); r = np.pi*k/L
    t = sum(c**(-(2*n + 0.5))*(2*n + 0.5)/((2*n + 0.5)**2 + r*r) for n in range(400))
    return phi(r) + 2*(-1)**k*t
for c in [1.2, 1.5, 2.0]:
    ks = np.arange(-20000, 20001); m = np.array([mu(int(k), c) for k in range(0, 20001)])
    mm = np.concatenate([m[:0:-1], m])
    best = max((np.log(2) - (np.log(np.sum(np.exp(-b*(mm - mm.min())))) - b*mm.min()))/b for b in [1.05, 1.1, 1.25, 1.5, 2, 3, 5, 8])
    print('lattice entropic bound c=%.2f: sum mu q >= %.4f (mu(0) = %.4f, min mu = %.4f)' % (c, best, m[0], mm.min()), flush=True)

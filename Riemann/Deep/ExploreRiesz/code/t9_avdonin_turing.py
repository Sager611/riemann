# T9: exact Avdonin-Turing identity for unfolded zeros lam_n = n - S_n (S_n = S(gamma_n) midpoint):
#   sum_{n=a}^{b-1} S_n = int_{gamma_a}^{gamma_b} S dlam + (S_b^2 - S_a^2)/2 - (S_b - S_a)/2,   dlam = theta'(t)/pi dt,
# so Avdonin block means = -(1/N)[L*Delta S_1 + boundary]; check numerically and compare block-sum sd with sqrt(2 Var S_1) L.
import numpy as np
from common import zeta_unfolded, theta, PI
g, lam = zeta_unfolded()
n = np.arange(1, len(g) + 1); S = n - lam                     # S_n
# exact int of S dlam over each gap: on (g_n, g_{n+1}) S = S_n + 1/2 - (lam(t) - lam_n)  => int = G(S_n + 1/2) - G^2/2
G = np.diff(lam); I = G*(S[:-1] + 0.5) - G**2/2
cum = np.concatenate([[0], np.cumsum(I)])                      # int_{g_1}^{g_n} S dlam
# also int S dt via antiderivative of theta (independent check of the dlam integral's size)
def Th(t): return (t**2/4)*np.log(t/(2*PI)) - 3*t**2/8 - PI*t/8 + np.log(t)/48
It = (n[:-1] - 1)*np.diff(g) + 0.5*np.diff(g) - (Th(g[1:]) - Th(g[:-1]))/PI  # S(t) = N(t) - 1 - theta/pi, N = n on (g_n, g_{n+1}); midpoint +1/2 shift
cumt = np.concatenate([[0], np.cumsum(It)])
rng = np.random.default_rng(3); errs = []
for _ in range(2000):
    a = rng.integers(1000, len(g) - 600); N = rng.integers(1, 512); b = a + N
    lhs = S[a:b].sum(); rhs = (cum[b] - cum[a]) + (S[b]**2 - S[a]**2)/2 - (S[b] - S[a])/2
    errs.append(abs(lhs - rhs))
print('identity max error over 2000 random blocks: %.3e' % max(errs))
# S_1 at zeros (int S dt from g_1) and block-sum decomposition
for N in [8, 64, 512]:
    idx = np.arange(1000, len(g) - N, N)
    bs = np.array([S[i:i+N].sum() for i in idx]); bs -= bs.mean()
    s1 = np.array([cumt[i+N] - cumt[i] for i in idx]); L = (np.log(g[idx]/(2*PI))/(2*PI))
    print(f'N={N:4d}: sd(block sum) {bs.std():.4f}  sd(L*dS1) {np.std(L*s1):.4f}  sd(boundary) {np.std([(S[i+N]**2-S[i]**2)/2-(S[i+N]-S[i])/2 for i in idx]):.4f}')
v = (1/(2*PI**2))*sum(1/(m**4*p**m*np.log(p)**2) for p in [q for q in range(2, 100000) if all(q % d for d in range(2, int(q**0.5) + 1))] for m in range(1, 8))
print('Var S_1 (random model, p < 1e5) = %.4f, sd = %.4f; sqrt(2 Var S1)*L(1e4) = %.4f' % (v, v**0.5, (2*v)**0.5*np.log(1e4/(2*PI))/(2*PI)))
s1all = cumt - np.mean(cumt[1000:])
print('sd of S_1(gamma_n) - mean over n > 1000: %.4f' % np.std(cumt[1000:] - np.polyval(np.polyfit(g[1000:], cumt[1000:], 1), g[1000:])))

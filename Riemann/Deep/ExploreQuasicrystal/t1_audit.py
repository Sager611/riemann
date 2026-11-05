# t1_audit.py -- explore5 T-A: which Fourier-quasicrystal hypotheses fail for the zeta zero measure, and by how much
# run: OMP_NUM_THREADS=1 uv run --no-project --with numpy python3 t1_audit.py
import numpy as np, math, sys
XMAX = 17.0                      # spectral side: n <= e^17 ~ 2.4e7
N = int(math.exp(XMAX)) + 1
s = np.ones(N + 1, dtype=bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
primes = np.nonzero(s)[0]
print("primes up to e^%.0f: %d" % (XMAX, len(primes)))
# prime powers with Lambda(n)/sqrt(n)
pp = []; lw = []
for p in primes:
    q = int(p); L = math.log(p)
    while q <= N:
        pp.append(q); lw.append(L / math.sqrt(q)); q *= int(p)
pp = np.array(pp, dtype=np.float64); lw = np.array(lw)
o = np.argsort(pp); pp = pp[o]; lw = lw[o]; logn = np.log(pp)
# (a) density: zeta zeros (1/2pi) log(T/2pi) vs prime-torus Kurasov-Sarnak density theta(P)/2pi (one unit of degree per prime)
print("\n(a) density match: theta(P) = sum_{p<=P} log p must equal log(T/2pi)")
th = np.cumsum(np.log(primes[:200].astype(float)))
for T in [1e3, 1e6, 1e12, 1e24, 1e100]:
    need = math.log(T / (2 * math.pi)); k = int(np.searchsorted(th, need)); P = int(primes[k])
    print("  T=%.0e  log(T/2pi)=%.2f  smallest P with theta(P)>=: P=%d (pi(P)=%d variables, theta=%.2f)" % (T, need, P, k + 1, th[k]))
# (b) spectral masses in unit windows: |mu^|((X,X+1]) = sum Lambda(n)/sqrt n vs pole term 2(e^{(X+1)/2}-e^{X/2})
print("\n(b) atomic Fourier mass of the zero measure in (X, X+1]  (FQ needs this translation-bounded)")
print("   X   mass        pole 2(e^((X+1)/2)-e^(X/2))  ratio   log(mass)/X")
for X in range(0, 16):
    m = lw[(logn > X) & (logn <= X + 1)].sum(); pole = 2 * (math.exp((X + 1) / 2) - math.exp(X / 2))
    print("  %2d  %10.3f  %10.3f  %6.4f  %s" % (X, m, pole, m / pole, ("%.3f" % (math.log(m) / X)) if X > 0 else "-"))
# (c) residual after subtracting the pole: R(X) = sum_{log n<=X} Lambda(n)/sqrt n - 2 e^{X/2}; tempered (polynomial in X) <=> RH
cum = np.cumsum(lw); R = cum - 2 * np.exp(logn / 2)
Rm = cum - lw - 2 * np.exp(logn / 2)      # just before each atom
absR = np.maximum(np.abs(R), np.abs(Rm))
print("\n(c) residual R(X) = sum_{log n<=X} Lambda(n)n^-1/2 - 2e^{X/2}  (RH <=> |R(X)| = O(X^2); unconditional only O(e^{(B-1/2+eps)X}), B = 0.87338...)")
for X in [5, 8, 10, 12, 14, 16, 17]:
    sel = logn <= X
    j = np.argmax(absR[sel]); 
    print("  max_{log n<=%2d} |R| = %8.4f at log n=%.3f ;  X^2/8 = %7.2f ; e^{0.3734 X} = %10.1f" % (X, absR[sel][j], logn[sel][j], X * X / 8, math.exp(0.37338196661 * X)))
# (d) spectral counting: rank-r torus spectra grow polynomially, zeta's spectrum {log n} exponentially
print("\n(d) spectral counting #{frequencies <= X}: zeta (all n) e^X vs P-smooth (rank pi(P) torus)")
def smooth_count(P, X):
    ps = [int(p) for p in primes if p <= P]; L = [math.log(p) for p in ps]
    cnt = 0; stack = [(0, 0.0)]
    while stack:
        i, acc = stack.pop()
        if i == len(ps): cnt += 1; continue
        a = acc
        while a <= X + 1e-12:
            stack.append((i + 1, a)); a += L[i]
    return cnt
for P, X in [(31, 10), (31, 20), (31, 30), (7, 20), (7, 30)]:
    c = smooth_count(P, X); r = len([p for p in primes if p <= P]); prodL = math.prod(math.log(p) for p in primes if p <= P)
    print("  P=%2d (rank %2d) X=%2d: #P-smooth n<=e^X = %9d ; leading X^r/(r! prod log p) = %11.1f ; all n: e^X = %.3e" % (P, r, X, c, X ** r / math.factorial(r) / prodL, math.exp(X)))

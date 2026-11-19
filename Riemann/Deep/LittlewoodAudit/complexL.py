# Falsification (ii): |L(1,chi)| for ALL nonprincipal characters of every prime conductor 3 <= q <= Q (all primitive),
# complex and real, via L(1,chi) = -(1/q) sum_a chi(a) psi(a/q) and one FFT per q, against Theorem 1(a):
# |L(1,chi)| >= 0.11/max(3, loglog q) (= 0.11/3 here).  Run: uv run --no-project --with numpy --with scipy python3 complexL.py Q
import sys, math, time
import numpy as np
from scipy.special import digamma
Q = int(sys.argv[1]); t0 = time.time()
s = np.ones(Q + 1, dtype=bool); s[:2] = False
for p in range(2, int(Q**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = [int(p) for p in np.nonzero(s)[0] if p >= 3]
def primroot(q):
    n = q - 1; fs = []; m = n; f = 2
    while f*f <= m:
        if m % f == 0:
            fs.append(f)
            while m % f == 0: m //= f
        f += 1
    if m > 1: fs.append(m)
    for g in range(2, q):
        if all(pow(g, n//f, q) != 1 for f in fs): return g
nchar = 0; minc = (9, None); minr = (9, None); maxc = (0, None)
chk = {}
for q in primes:
    g = primroot(q); a = np.empty(q - 1, dtype=np.int64); a[0] = 1
    for j in range(1, q - 1): a[j] = a[j-1]*g % q
    Lv = -(q - 1)*np.fft.ifft(digamma(a/q))/q          # Lv[k] = L(1, chi_k), chi_k(g^j) = e(jk/(q-1))
    A = np.abs(Lv); kr = (q - 1)//2
    if A[kr] < minr[0]: minr = (A[kr], (q, kr))
    Ac = A.copy(); Ac[0] = np.inf; Ac[kr] = np.inf; k = int(np.argmin(Ac))
    if Ac[k] < minc[0]: minc = (Ac[k], (q, k))
    k = int(np.argmax(A[1:])) + 1
    if A[k] > maxc[0]: maxc = (A[k], (q, k))
    nchar += q - 2
    if q in (5, 23, 163): chk[q] = Lv[kr]
print(f'Q = {Q}: {len(primes)} prime conductors, {nchar} nonprincipal (primitive) characters, t = {time.time()-t0:.1f}s')
print(f'sanity: L(1,(./5)) = {chk[5].real:.12f} (2 log((1+sqrt5)/2)/sqrt5 = {2*math.log((1+5**0.5)/2)/5**0.5:.12f}); '
      f'L(1,chi_-23) = {chk[23].real:.12f} (3 pi/sqrt23 = {3*math.pi/23**0.5:.12f}); L(1,chi_-163) = {chk[163].real:.12f} (pi/sqrt163 = {math.pi/163**0.5:.12f})')
b = 0.11/3
print(f'min |L(1,chi)| complex chi: {minc[0]:.6f} at (q, k) = {minc[1]};  real chi: {minr[0]:.6f} at {minr[1]};  bound 0.11/3 = {b:.6f}: {"OK" if min(minc[0], minr[0]) >= b else "VIOLATED"}')
print(f'max |L(1,chi)|: {maxc[0]:.4f} at {maxc[1]} (side result (e) needs loglog q >= 3: not testable here)')

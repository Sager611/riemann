# t3c_scale.py (explore32): the parts of Q at scale (double precision, CCM basis N = 80): archimedean and prime parts are O(1)-O(10),
# Q's minimum is ~e^{-4 pi c}: every bound that treats the parts separately loses.
import numpy as np, json
from flint import fmpq, ctx
src = open('t3_windows.py').read().split('def K(u)')[0]
exec(src)
D = -5.372183419225666
def ppow(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
K = lambda u: sum(np.exp(-(2*m + 0.5)*u)/(2*m + 0.5) for m in range(int(60/(2*u)) + 2))
for c in [5, 10, 20, 47]:
    a = np.log(c)/2
    W = sum(np.log(ppow(n)[0])/np.sqrt(n) for n in range(2, c) if ppow(n))
    row = {'c': c, 'D': D, 'killing': round(D + 2*K(a), 4), 'D_c = D - 2 sum w_n': round(D - 2*W, 4), 'sum w_n': round(W, 4)}
    for sec in ['even', 'odd']:
        A, P, R = parts(fmpq(c), 80, sec)
        row[sec] = {'arch': round(float(np.linalg.eigvalsh(A)[0]), 4), 'arch_pole': round(float(np.linalg.eigvalsh(A + P)[0]), 4), 'primes_min': round(float(np.linalg.eigvalsh(R)[0]), 4), 'primes_max': round(float(np.linalg.eigvalsh(R)[-1]), 4)}
    print(json.dumps(row), flush=True)

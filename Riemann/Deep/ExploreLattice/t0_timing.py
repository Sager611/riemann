# t0_timing.py (explore48): python-flint API check + timing of Dirichlet L at height t for conductors 20..5460.
import time, math
from flint import acb, arb, dirichlet_char, dirichlet_group, ctx
def kron(d, a):
    # Kronecker symbol (d/a), a > 0
    from math import gcd
    if gcd(d, a) != 1: return 0
    res = 1; aa = a
    while aa % 2 == 0:
        aa //= 2
        if d % 8 in (3, 5, -3, -5): res = -res
    # Jacobi (d/aa) for odd aa
    x, n = d % aa, aa; j = 1
    while x:
        while x % 2 == 0:
            x //= 2
            if n % 8 in (3, 5): j = -j
        x, n = n, x
        if x % 4 == 3 and n % 4 == 3: j = -j
        x %= n
    return res * (j if n == 1 else 0)
def real_char(d):
    q = abs(d)
    test = [a for a in range(2, 4*q) if math.gcd(a, q) == 1][:40]
    for n in range(1, q):
        if math.gcd(n, q) != 1: continue
        if (n*n) % q != 1: continue  # real characters: order <= 2 in Conrey labels
        chi = dirichlet_char(q, n)
        ok = True
        for a in test:
            v = chi(a)
            if abs(float(v.real.mid()) - kron(d, a)) > 1e-9: ok = False; break
        if ok: return chi
    raise ValueError(d)
if __name__ != '__main__': pass
else:
  print('ctx', ctx)
for d in (() if __name__ != '__main__' else (-20, 5, -4, -427, -1435, -5460, -3)):
    t0 = time.time(); chi = real_char(d); t1 = time.time()
    for t in (100.0, 300.0):
        s = acb(0.5, t); t2 = time.time(); v = s.dirichlet_l(chi); t3 = time.time()
        print(d, 'cond', chi.conductor(), 't', t, 'L', v.str(radius=False)[:40], 'time %.4f s' % (t3-t2), 'charsearch %.2f s' % (t1-t0))

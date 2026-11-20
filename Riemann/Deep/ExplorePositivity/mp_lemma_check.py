# mp_lemma_check.py — sanity: recurrence = Romik's 2F1 definition; MP expansion reproduces Xi;
# Romik's Polya-type lemma (real-rooted p => sum k!/(3/2)_k a_k f_k real-rooted) on random p.
import random, mpmath as mp
from mp_core import *
mp.mp.dps = 40
F = f_coeffs(12)
err = 0
for n in range(0, 9):
    for x in [mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('-2.2')]:
        a = sum(c * x**k for k, c in enumerate(F[n]))
        b = f_hyp(n, x)
        err = max(err, abs(a - b))
print('recurrence vs 2F1 max err', mp.nstr(err, 3))
# lemma test
random.seed(1)
worst = 0
for trial in range(200):
    d = random.randint(2, 10)
    roots = [random.uniform(-3, 3) for _ in range(d)]
    p = [mp.mpf(1)]
    for r in roots:
        q = [mp.mpf(0)] * (len(p) + 1)
        for k, v in enumerate(p):
            q[k+1] += v; q[k] -= r * v
        p = q
    T = [mp.mpf(0)] * (d + 1)
    for k, a in enumerate(p):
        w = mp.factorial(k) / mp.rf(mp.mpf(3)/2, k) * a
        for j, c in enumerate(F[k]):
            T[j] += w * c
    rr = mp.polyroots(T[::-1], maxsteps=200, extraprec=200)
    worst = max(worst, max(abs(mp.im(z)) for z in rr))
print('lemma: 200 random real-rooted p (deg 2..10): max |Im root| of T[p] =', mp.nstr(worst, 3))
# negative control: p with non-real roots
p = [mp.mpf(1), 0, mp.mpf(1)]  # 1 + z^2
T = [mp.mpf(0)] * 3
for k, a in enumerate(p):
    w = mp.factorial(k) / mp.rf(mp.mpf(3)/2, k) * a
    for j, c in enumerate(F[k]):
        T[j] += w * c
print('control 1+z^2 -> T coeffs', [mp.nstr(v, 6) for v in T], 'roots', [mp.nstr(z, 6) for z in mp.polyroots(T[::-1])])

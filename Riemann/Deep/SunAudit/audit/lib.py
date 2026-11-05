# sunaudit common exact helpers (stdlib only)
from math import gcd, isqrt

def primes_upto(n):
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

def leg(a, p):
    """Legendre symbol (a/p), p odd prime, via Euler's criterion."""
    a %= p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return 1 if r == 1 else -1

def jacobi(a, n):
    """Jacobi symbol (a/n), n odd positive."""
    assert n > 0 and n % 2 == 1
    a %= n
    res = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0

def h_forms(D):
    """Class number h(D) for negative discriminant D (counts primitive reduced forms ax^2+bxy+cy^2,
    |b|<=a<=c, b>=0 if |b|==a or a==c)."""
    assert D < 0 and D % 4 in (0, 1)
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if b < 0 and a == c:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h

def h_dirichlet_minus_p(p):
    """h(-p) for prime p = 3 mod 4, p > 3, via Dirichlet: h = -(1/p) sum a (a/p)."""
    assert p % 4 == 3 and p > 3
    s = sum(a * leg(a, p) for a in range(1, p))
    assert s % p == 0
    return -s // p

def h_minus_p(p):
    """class number of Q(sqrt(-p)), p prime = 3 mod 4 (disc -p)."""
    h1 = h_forms(-p)
    if p > 3:
        h2 = h_dirichlet_minus_p(p)
        assert h1 == h2, (p, h1, h2)
    return h1

def inv_count(seq):
    """number of pairs i<j with seq[i] > seq[j] (strict)."""
    n = len(seq)
    c = 0
    for i in range(n):
        si = seq[i]
        for j in range(i + 1, n):
            if si > seq[j]:
                c += 1
    return c

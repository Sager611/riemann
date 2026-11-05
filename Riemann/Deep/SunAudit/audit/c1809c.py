import sys
sys.path.insert(0, '.')
from lib import *
def inv_strict(seq):
    # merge-sort inversion count (strict >)
    def rec(a):
        if len(a) <= 1: return a, 0
        m = len(a)//2
        l, x = rec(a[:m]); r, y = rec(a[m:])
        out = []; i = j = 0; c = x + y
        while i < len(l) and j < len(r):
            if l[i] <= r[j]: out.append(l[i]); i += 1
            else: out.append(r[j]); j += 1; c += len(l) - i
        out += l[i:]; out += r[j:]
        return out, c
    return rec(list(seq))[1]
bad61 = []; bad612 = []
for p in primes_upto(3000):
    if p < 3: continue
    n = (p - 1)//2
    if p % 4 == 1:
        s = inv_strict([k*k % p for k in range(1, n+1)])
        t = 0
        for j in range(1, n+1):
            jj = j*j
            for k in range(j+1, n+1):
                if 2*((k*k - jj) % p) > p: t += 1
        rhs = sum(1 for k in range(1, p) if 4*k < p and leg(k, p) == 1)
        if (s + t - rhs) % 2: bad61.append(p)
    c8 = inv_strict([pow(j, 8, p) for j in range(1, n+1)])
    if p % 8 == 1: r8 = sum(1 for k in range(1, p) if 4*k < p and leg(k, p) == 1)
    elif p % 8 == 3: r8 = 0
    elif p % 8 == 5: r8 = (p - 5)//8
    else: r8 = (h_minus_p(p) + 1)//2
    if (c8 - r8) % 2: bad612.append(p)
print('1809.07766 Conj 6.1 failures p<3000:', bad61, '| (6.12) failures p<3000:', bad612)

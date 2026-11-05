import sys, cmath, math
sys.path.insert(0, '.')
from lib import *

def s_p(p):
    n = (p - 1) // 2
    return inv_count([k * k % p for k in range(1, n + 1)])

def t_p(p):
    n = (p - 1) // 2
    c = 0
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            if 2 * ((k * k - j * j) % p) > p:
                c += 1
    return c

def cnt_qr_below_quarter(p, val):
    # #{1<=k<p/4 : (k/p)=val}
    return sum(1 for k in range(1, p) if 4 * k < p and leg(k, p) == val)

def R(k, n):
    r = k % n
    return min(r, n - r)

out = []
P = [p for p in primes_upto(400) if p > 2]
# Conj 6.1
bad61 = []
for p in P:
    if p % 4 == 1:
        lhs = (s_p(p) + t_p(p)) % 2
        rhs = cnt_qr_below_quarter(p, 1) % 2
        if lhs != rhs:
            bad61.append((p, s_p(p), t_p(p), cnt_qr_below_quarter(p, 1)))
print('6.1 failures p<400:', bad61)
print('6.1 detail p=5,13,17:', [(p, s_p(p), t_p(p), cnt_qr_below_quarter(p,1)) for p in (5,13,17,29)])
# sanity Theorem 1.4: (-1)^s = (-1)^t = 1 if p=3 mod 8, (-1)^((h+1)/2) if p=7 mod 8
for p in P:
    if p % 8 == 3:
        assert s_p(p) % 2 == 0 and t_p(p) % 2 == 0, p
    if p % 8 == 7:
        e = ((h_minus_p(p) + 1) // 2) % 2
        assert s_p(p) % 2 == e and t_p(p) % 2 == e, p
print('Thm1.4 sanity ok')

# Conj 6.2 for all a
bad62a = []; bad62b = []
for p in [q for q in P if q < 200]:
    n = (p - 1) // 2
    for a in range(1, p):
        seq = [R(a * i * i, p) for i in range(1, n + 1)]
        if inv_count(seq) % 2 != ((p + 1) // 8) % 2:
            bad62a.append((p, a))
        c = sum(1 for i in range(n) for j in range(i + 1, n) if 2 * (seq[i] + seq[j]) > p)
        if p % 4 == 1:
            rhs = (-1) ** cnt_qr_below_quarter(p, -1) * leg(a, p) ** ((1 - leg(2, p)) // 2)
        else:
            rhs = 1
        if (-1) ** c != rhs:
            bad62b.append((p, a))
print('6.2 (6.2) failures p<200 all a:', bad62a[:20], len(bad62a))
print('6.2 (6.3) failures p<200 all a:', bad62b[:20], len(bad62b))

# Conj 6.3
bad64 = []; bad65 = []
for p in [q for q in P if q > 3]:
    n = (p - 1) // 2
    seq = [(j * (j + 1) // 2) % p for j in range(1, n + 1)]
    if p % 4 == 3:
        lhs = (-1) ** inv_count(seq)
        rhs = (-1) ** ((h_minus_p(p) + 1) // 2 + sum(1 for k in range(1, (p + 1) // 8 + 1) if leg(k, p) == 1))
        if lhs != rhs: bad64.append(p)
    c = sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] + seq[j] > p)
    if p % 8 == 1: rhs = (-1) ** ((p - 1) // 8)
    elif p % 8 == 5: rhs = (-1) ** cnt_qr_below_quarter(p, -1)
    else: rhs = (-1) ** ((h_minus_p(p) + 1) // 2 + sum(1 for k in range(1, (p + 1) // 8 + 1) if leg(k, p) == -1))
    if (-1) ** c != rhs: bad65.append(p)
print('6.3 (6.4) failures:', bad64)
print('6.3 (6.5) failures:', bad65)

# Conj 6.4
bad66 = []; bad67 = []
for p in P:
    n = (p - 1) // 2
    seq = [(j * (j + 1)) % p for j in range(1, n + 1)]
    if p % 4 == 3:
        if (-1) ** inv_count(seq) != (-1) ** ((p + 1) // 8): bad66.append(p)
    c = sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] + seq[j] > p)
    if p % 4 == 1: rhs = (-1) ** ((p - 1) // 8)
    elif p % 8 == 3 and p > 3: rhs = (-1) ** ((h_minus_p(p) + 1) // 2)
    elif p % 8 == 7: rhs = 1
    else: rhs = None
    if rhs is not None and (-1) ** c != rhs: bad67.append(p)
print('6.4 (6.6) failures:', bad66)
print('6.4 (6.7) failures:', bad67)

# Conj 6.5(i)
bad68 = []; bad69 = []
for p in P:
    if p % 6 == 5:
        n = (p - 1) // 2
        d = sum(1 for k in range(1, n + 1) if 2 * (pow(k, 3, p)) > p) - (p + 1) // 6
        if d < 0 or d % 2: bad68.append((p, d))
        if inv_count([pow(j, 3, p) for j in range(1, p)]) % 2 != ((p + 1) // 6) % 2: bad69.append(p)
print('6.5 (6.8) failures:', bad68)
print('6.5 (6.9) failures:', bad69)

# Conj 6.6
bad611 = []; bad612 = []
for p in P:
    n = (p - 1) // 2
    c4 = inv_count([pow(j, 4, p) for j in range(1, n + 1)])
    rhs = (p + 1) // 8 + (((h_minus_p(p) + 1) // 2) if p % 8 == 7 else 0)
    if c4 % 2 != rhs % 2: bad611.append((p, c4, rhs))
    c8 = inv_count([pow(j, 8, p) for j in range(1, n + 1)])
    if p % 8 == 1: r8 = cnt_qr_below_quarter(p, 1)
    elif p % 8 == 3: r8 = 0
    elif p % 8 == 5: r8 = (p - 5) // 8
    else: r8 = (h_minus_p(p) + 1) // 2
    if c8 % 2 != r8 % 2: bad612.append((p, c8, r8))
print('6.6 (6.11) failures:', bad611)
print('6.6 (6.12) failures:', bad612)

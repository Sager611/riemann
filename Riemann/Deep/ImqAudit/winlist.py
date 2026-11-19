# winlist.py LO HI: every fundamental d in [LO, HI] (d = |disc|), one per line
import sys
LO, HI = int(sys.argv[1]), int(sys.argv[2])
n = HI - LO + 1
bad = bytearray(n)  # bad[i]: some odd p^2 | LO+i
import math
lim = math.isqrt(HI) + 1
sieve = bytearray(lim + 1)
for p in range(3, lim + 1, 2):
    if sieve[p]: continue
    for j in range(p * p, lim + 1, 2 * p): sieve[j] = 1
    q = p * p
    if q > HI: break
    st = (-LO) % q
    bad[st::q] = b'\x01' * len(range(st, n, q))
out = []
for i in range(n):
    d = LO + i; r = d % 16
    if r in (3, 7, 11, 15) and not bad[i]: out.append(d)
    elif r in (4, 8) and not bad[i]:   # d = 4m, odd p^2 | d iff odd p^2 | m; 4 | m excluded by r
        out.append(d)
sys.stdout.write('\n'.join(map(str, out)) + '\n')

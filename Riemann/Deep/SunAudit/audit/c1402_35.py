import sys
sys.path.insert(0, '.')
from lib import *
N = 200000
P = primes_upto(4 * 10**6)
Pset = set(P)
pr = [None] + P  # 1-indexed
assert len(P) > N
reps = {}
missing = []
for s in range(3, 20001):
    found = None
    for k in range(1, (s + 1) // 2):
        m = s - k
        if k >= m: break
        q = 2 * pr[m] - pr[k]
        if q in Pset:
            found = (k, m, pr[k], pr[m], q); break
    if found is None:
        missing.append(s)
    elif s < 12:
        reps[s] = found
print('s in 3..20000 without representation:', missing)
print('small reps:', reps)

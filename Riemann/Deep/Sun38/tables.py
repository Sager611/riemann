#!/usr/bin/env python3
"""tables.py (stdlib): small case tables behind the proofs in REPORT.md (output: tables.txt)."""
import itertools

out = []
# TeX Conj3.5(i) (printed 3.4(i)): A = {0,1,2,3,4} in Z/7, all 24 cyclic orders (a1 = 0)
out.append('Z/7, A = {0,1,2,3,4}: every cyclic order repeats a sum a_i + 2a_{i+1} (mod 7)')
for rest in itertools.permutations([1, 2, 3, 4]):
    a = (0,) + rest
    s = [(a[i] + 2 * a[(i + 1) % 5]) % 7 for i in range(5)]
    rep = next((i, j) for i in range(5) for j in range(i + 1, 5) if s[i] == s[j])
    out.append('  order %s  sums %s  s%d = s%d = %d' % (''.join(map(str, a)), ''.join(map(str, s)), rep[0] + 1, rep[1] + 1, s[rep[0]]))
# TeX Conj3.3 (printed 3.2): Z/6, n = 4, failing (A, a1)
out.append('Z/6, n = 4: all (A, a1) with no permutation b (b1 = a1) having distinct differences b_i - b_{i+1}')
for A in itertools.combinations(range(6), 4):
    for a1 in A:
        rest = [x for x in A if x != a1]
        ok = any(len({(b[i] - b[i + 1]) % 6 for i in range(3)}) == 3 for b in ((a1,) + p for p in itertools.permutations(rest)))
        if not ok: out.append('  A = %s, a1 = %d (A minus a1 = %s, a coset of {0,2,4})' % (set(A), a1, set(rest)))
open('tables.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))

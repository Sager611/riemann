"""gapfamily: exact criterion for GL(2): H = -sum_j a_j Lam^{2j}/(2j+1) with the exact rational form certified by cert5.py
(logs/cert5_<tag>.json). H > 0 and the per-gap certificate give Lambda_f > Lam (Criterion J, uniform measure, 1410.7765 Thm mixed)."""
import sys, json
from fractions import Fraction as Fr
tag = sys.argv[1]; lam = Fr(sys.argv[2])
c = json.load(open('logs/cert5_%s.json' % tag))
a = [Fr(x) for x in c['a']]
assert a[1] == 1, a[1]
H = -sum(a[j] * lam ** (2 * j) / (2 * j + 1) for j in range(len(a)))
print('tag %s: per-gap ok = %s; Lam = %s; H = %s = %.10g (exact rational) -> %s' % (tag, c['ok'], lam, H if H.denominator < 10**6 else '<rational>', float(H), 'CRITERION HOLDS' if (H > 0 and c['ok']) else 'FAILS'))
json.dump({'tag': tag, 'lam': str(lam), 'H': str(H), 'Hfloat': float(H), 'ok': bool(H > 0 and c['ok'])}, open('logs/gl2crit_%s.json' % tag, 'w'), indent=1)

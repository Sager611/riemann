import json
from fractions import Fraction as Fr
C = json.load(open('../IntCheb/CERTIFICATE_LOWER_new.json'))
lib = [[int(a) for a in q] for q in C['library']]; cs = C['multipliers']
def ev(q, x):
    v = Fr(0)
    for a in q: v = v * x + a
    return v
act = [(q, c) for q, c in zip(lib, cs) if c > 0]
for name, q0 in (('x', [1, 0]), ('x-1', [1, -1]), ('2x-1', [2, -1])):
    root = Fr(-q0[1], q0[0]); a = q0[0]
    vals = []
    for q, c in act:
        if q == q0: continue
        R = a ** (len(q) - 1) * ev(q, root)      # Res(q0, q) = a^{deg q} q(root)
        vals.append(abs(R))
    print(name, 'max |Res(.,Q_j)| over active j != i:', max(vals), ' all equal 1:', all(v == 1 for v in vals))

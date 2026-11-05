#!/usr/bin/env python3
# extract sofa_certificate.json into plain files for the independent C checker (sidecheck)
import json, base64, sys, os
from fractions import Fraction as F
src = sys.argv[1]; out = sys.argv[2]; os.makedirs(out, exist_ok=True)
d = json.load(open(src))
T = F(d['T']); W = F(d['W'])
def quadrant(j, c, s): return [(c, s), (-s, c), (-c, -s), (s, -c)][j]
with open(os.path.join(out, 'cases.txt'), 'w') as f:
    f.write('%d %.17g %.17g\n' % (len(d['cases']), float(T), float(W)))
    for idx, case in enumerate(d['cases']):
        corr = [(int(j), F(c), F(s)) for j, c, s in case['corridors']]
        raw = base64.b64decode(case['tree_b64'])
        nbits = int.from_bytes(raw[:8], 'little')
        open(os.path.join(out, 'tree%d.bin' % idx), 'wb').write(raw[8:])
        f.write('%s|%d|%d\n' % (case['name'].replace(' ', '_'), len(corr), nbits))
        for (j, c, s), dom in zip(corr, case['domain']):
            cc, ss = quadrant(j, c, s)
            dv = [F(v) for v in dom]
            f.write('%d %.17g %.17g %s\n' % (j, float(cc), float(ss), ' '.join('%.17g' % float(v) for v in dv)))
            assert all(F(float(v)) == v for v in dv), 'domain not exactly representable'
print('ok')

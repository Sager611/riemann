#!/usr/bin/env python3
"""Search code (not needed for verification): builds the branch-and-bound trees with the C program sofabb
and writes sofa_certificate.json.  Usage: python3 make_certificate.py [T_num T_den]"""
import json, base64, math, os, subprocess, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
def quadrant(j, c, s):
    return [(c, s), (-s, c), (-c, -s), (s, -c)][j]
def dyadic_up(x, den):  return F(math.ceil(x * den), den)
def dyadic_down(x, den): return F(math.floor(x * den), den)
def main():
    T = F(int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else F(19, 10)
    phi = (F(21, 29), F(20, 29)); phi2 = (F(20, 29), F(21, 29))
    cases = []
    for (j1, j2) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
        cases.append(('pair %d%d' % (j1, j2), [(j1,) + phi, (j2,) + phi]))
    for j in (1, 3):
        cases.append(('same quadrant %d' % j, [(j,) + phi, (j,) + phi2]))
    W = dyadic_up(min((1 + c) / s + (1 + s) / c for c, s in (phi, phi2)), 32)
    Kv = [(F(0), F(0)), (W, F(0)), (W, F(1)), (F(0), F(1))]
    out = []
    for name, corrs in cases:
        rot = [quadrant(j, c, s) for j, c, s in corrs]
        dom = []
        for c, s in rot:
            ps = [c * x + s * y for x, y in Kv]; qs = [-s * x + c * y for x, y in Kv]
            dom.append((dyadic_down(min(ps) - 1, 64), dyadic_up(max(ps), 64), dyadic_down(min(qs) - 1, 64), dyadic_up(max(qs), 64)))
        cfg = os.path.join(HERE, 'tmp_cfg.txt'); bits = os.path.join(HERE, 'tmp_tree.bits')
        with open(cfg, 'w') as f:
            f.write('%d %r %r %r %d\n' % (len(rot), float(W), float(T), 1e-9, 60))
            for c, s in rot:
                f.write('%d %d %d %d\n' % (c.numerator, c.denominator, s.numerator, s.denominator))
            for d in dom:
                f.write('%r %r %r %r\n' % tuple(float(v) for v in d))
        r = subprocess.run([os.path.join(HERE, 'sofabb'), cfg, bits], capture_output=True, text=True)
        print(name, r.stdout.strip())
        if r.returncode != 0: raise SystemExit('search failed for ' + name)
        raw = open(bits, 'rb').read()
        nb = int.from_bytes(raw[:8], 'little'); raw = raw[:8 + (nb + 7) // 8]
        out.append({'name': name, 'corridors': [[j, str(c), str(s)] for j, c, s in corrs],
                    'domain': [[str(v) for v in d] for d in dom], 'tree_b64': base64.b64encode(raw).decode()})
        os.remove(cfg); os.remove(bits)
    cert = {'constant': 'C_41b (ambidextrous moving sofa constant), upper bound',
            'T': str(T), 'W': str(W), 'cases': out}
    json.dump(cert, open(os.path.join(HERE, 'sofa_certificate.json'), 'w'))
    print('wrote sofa_certificate.json')
if __name__ == '__main__':
    main()

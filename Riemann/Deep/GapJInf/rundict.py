"""gapjinf: run the dictionary bound h_dict(J) for one amplifier's moments.
usage: rundict.py forms.npz amp.json dictname Jmax [nth] [Xi] [scale s: dilate nu by s]"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ampmom import load
import dictsdp as DS

PI = np.pi


def dictionary(name):
    pats = []
    oms = [1.0, 1.01, 1.03, 1.06, 1.1, 1.15, 1.2, 1.3, 1.4, 1.55, 1.75, 2.0]
    if name in ('sep', 'mix', 'mix2'):
        pats += [[PI / w] for w in oms]
    if name in ('mix', 'mix2'):
        ts = [0.15, 0.3, 0.45, 0.6, 0.75, 0.9]
        pats += [[PI, PI * t] for t in ts]
        pats += [[PI, PI, PI * t] for t in ts]
        pats += [[PI, PI * t, PI * t] for t in ts]
        pats += [[PI, PI, PI, PI * t] for t in ts]
    if name == 'mix2':
        rng = np.random.default_rng(1)
        for m in (2, 3, 4, 5):
            for _ in range(12):
                g = PI * rng.uniform(0.2, 1.0, m); g[0] = PI
                pats.append(list(g))
    return pats


if __name__ == '__main__':
    npz, amj, name, Jmax = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    nth = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    Xi = float(sys.argv[6]) if len(sys.argv) > 6 else 3.6
    s = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    sv, am, kap, v = load(npz, amj)
    m = sv[0::2] * s ** (2 * np.arange(len(sv[0::2])))
    pats = dictionary(name)
    t0 = time.time()
    blocks = []
    for g in pats:
        blocks += DS.pattern_blocks(g, nth, Xi)
    print('kappa=%.4f v=%.6f dilation s=%.4f (kappa_eff %.4f) dict=%s: %d patterns, %d blocks, max size %d (%.1f s)' % (
        kap, v, s, kap * s, name, len(pats), len(blocks), max(N.shape[1] for _, N in blocks), time.time() - t0), flush=True)
    for J in range(1, Jmax + 1):
        t0 = time.time()
        val, a, st = DS.h_dict(m, blocks, J, ret=True)
        print('J=%d: h_dict = %.6f  [%s, %.1f s]  a = %s' % (J, val, st, time.time() - t0, ' '.join('%.4g' % x for x in a / a[1])), flush=True)

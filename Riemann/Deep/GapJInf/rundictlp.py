"""usage: rundictlp.py forms.npz amp.json dict Jmax [nth] [Xi] [s]"""
import sys, os, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ampmom import load
import dictlp as DL
PI = np.pi


def dictionary(name, seed=1):
    pats = [[PI / w] for w in [1.0, 1.01, 1.03, 1.06, 1.1, 1.15, 1.2, 1.3, 1.4, 1.55, 1.75, 2.0]]
    if name == 'sep':
        return pats
    ts = [0.15, 0.3, 0.45, 0.6, 0.75, 0.9]
    pats += [[PI, PI * t] for t in ts] + [[PI, PI, PI * t] for t in ts] + [[PI, PI * t, PI * t] for t in ts]
    pats += [[PI, PI, PI, PI * t] for t in ts]
    if name == 'mix':
        return pats
    rng = np.random.default_rng(seed)
    nrand = {'mix2': 120, 'mix3': 600}[name]
    for i in range(nrand):
        m = int(rng.integers(2, 8))
        g = PI * rng.uniform(0.1, 1.0, m)
        if rng.random() < 0.7:
            g[rng.integers(0, m)] = PI
        pats.append(list(g))
    return pats


if __name__ == '__main__':
    npz, amj, name, Jmax = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    nth = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    Xi = float(sys.argv[6]) if len(sys.argv) > 6 else 4.0
    s = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    Jmin = int(os.environ.get('JMIN', '1'))
    sv, am, kap, v = load(npz, amj)
    m = sv[0::2] * s ** (2 * np.arange(len(sv[0::2])))
    pats = dictionary(name)
    t0 = time.time()
    fib = []
    for g in pats:
        fib += DL.pattern_fibres(g, nth, Xi)
    print('kappa=%.4f v=%.6f s=%.4f (kappa_eff %.4f) dict=%s: %d patterns, %d fibres, max dim %d (%.1f s)' % (
        kap, v, s, kap * s, name, len(pats), len(fib), max(N.shape[1] for _, N in fib), time.time() - t0), flush=True)
    for J in range(Jmin, Jmax + 1):
        t0 = time.time()
        L = DL.DictLP(m, J, fib)
        val, a, it = L.run(verbose=False)
        print('J=%d: h_dict = %.6f  [%d rounds, %d cuts, %.1f s]  a/a1 = %s' % (J, val, it, len(L.cuts), time.time() - t0, ' '.join('%.4g' % x for x in a / a[1])), flush=True)

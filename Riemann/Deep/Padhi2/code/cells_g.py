"""padhi2: Padhi's c5 uses Bellotti's |G(1/2+it)| <= d5 W/log W (proof invalid). The certified bound is 9.75707 W/log^2 W
(BellottiFix cert.py), so the printed c5 is valid iff d5*wL >= 9.75707 in every cell (W = T^w, log W = wL). Same for any
cell parameters actually used (the certified cells take p from code/lemmaD_tableGH0.01_*.pkl)."""
import pickle, glob, os, math, sys
sys.path.insert(0, os.path.dirname(__file__))
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
REPO = os.environ.get('PADHI_REPO', os.path.join(ROOT, '.scratch', 'padhi2', 'repo', 'ramanujans-inequality-main'))
tab = {}
for f in sorted(glob.glob(os.path.join(REPO, 'code', 'lemmaD_tableGH0.01_*.pkl'))): tab.update(pickle.load(open(f, 'rb')))
used = [(k, v) for k, v in tab.items() if v is not None and v[0] <= 1e8]
wl = [(float(v[1][2])*float(k[1]), k) for k, v in used]
d5wl = []
for x, k in wl:
    d5 = 1.0061*(1 + (1.333 + 4*math.exp(-x/2))/x); d5wl.append((d5*x, k))
print('cells used:', len(used), ' wL range: %.4f (%s) .. %.4f' % (min(wl)[0], min(wl)[1], max(wl)[0]))
print('min d5*wL = %.4f at %s (need >= 9.75707): %s' % (min(d5wl)[0], min(d5wl)[1], min(d5wl)[0] >= 9.75707))

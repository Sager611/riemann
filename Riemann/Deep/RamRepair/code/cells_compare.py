# cells_compare.py -- compare our certified cells (out/cells_final.json) with Padhi's published cells (certify/arb/cells_arb_*.json)
import json, glob, os, math
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'repo', 'ramanujans-inequality-main', 'certify', 'arb')
old = {}
for f in glob.glob(os.path.join(R, 'cells_arb_*.json')):
    for k, v in json.load(open(f)).items():
        if v is not None: old[k] = float(v['logN_up'])
mine = json.load(open(os.path.join(HERE, 'out', 'cells_final.json')))
tot = len(mine); none = [k for k, v in mine.items() if v is None or v.get('logN_up') is None]
d = []
for k, v in mine.items():
    if v is None or v.get('logN_up') is None or k not in old: continue
    d.append((float(v['logN_up']) - old[k], k, float(v['logN_up']), old[k]))
d.sort()
print('cells', tot, 'dropped (inadmissible/c1 fails)', len(none), none[:10])
print('compared', len(d), 'mine - published logN: min %.4f  max %.4f  mean %.4f' % (d[0][0], d[-1][0], sum(x[0] for x in d)/len(d)))
print('largest increases', [(k, round(x, 4)) for x, k, a, b in d[-5:]])
print('largest decreases', [(k, round(x, 4)) for x, k, a, b in d[:5]])
missing = [k for k in old if k not in mine]
print('published cells not in ours', len(missing), missing[:10])
d4s = [float(v['d4'].split()[0].strip('[')) for v in mine.values() if v is not None]
print('d4 range %.4f .. %.4f' % (min(d4s), max(d4s)))
dens = [float(v['den'].split()[0].strip('[')) for v in mine.values() if v is not None]
print('c1^2 - d4 c5r range %.4f .. %.4f' % (min(dens), max(dens)))


# merge5.py: shards -> Riemann/Deep/DBN/EVIDENCE2/floors_dbn5.json (+ certs_dbn5.json.gz with per-subrange (c, Ctop, K, floor))
import json, gzip, glob, sys, os
import lib5 as L
from fractions import Fraction as F
E2 = os.path.join(L._repo(), 'Riemann/Deep/DBN/EVIDENCE2') + '/'
recs = {}
for fn in sorted(glob.glob('out/shard_*.jsonl')):
    for line in open(fn):
        r = json.loads(line); recs[(r['cell'], r['h'])] = r
SRC = os.path.join(os.environ.get('DBN4V_WORK', os.path.join(L._repo(), '.scratch', 'dbn4v')), 'arxiv/anc/certificate/inputs/refined') + '/'
plans = [[F(r['height']) for r in json.load(open(SRC + n))['plan']['rows']] for n in ('SOURCE_OLD_ENHANCED.json', 'SOURCE_RETUNED.json')]
need = sum(len(set([plans[0][i], plans[1][i]])) for i in range(608))
missing = [(i, str(h)) for i in range(608) for h in set([plans[0][i], plans[1][i]]) if (i, str(h)) not in recs]
if missing: sys.exit('missing combos: %s' % missing[:5])
if len(recs) != need: sys.exit('incomplete: %d of %d' % (len(recs), need))
fl = []; ce = []
for (cell, h) in sorted(recs, key=lambda k: (k[0], F(k[1]))):
    r = recs[(cell, h)]
    assert len(r['shifts']) == 10
    fl.append({'cell': cell, 'h': h, 'floor': {s: v['floor'] for s, v in r['shifts'].items()}})
    ce.append({'cell': cell, 'h': h, 'shifts': {s: {'Ep_upper': v['Ep'], 'subs': v['subs']} for s, v in r['shifts'].items()}})
json.dump({'schema': 'DBN5_PAIRED_FLOORS_V1', 'meaning': 'certified lower bound of the translated-probe average for t in cell i, x >= X, height h, shift r (premises: see REPORT6)', 'floors': fl}, open(E2 + 'floors_dbn5.json', 'w'), separators=(',', ':'), sort_keys=True)
with gzip.open(E2 + 'certs_dbn5.json.gz', 'wt') as f: json.dump(ce, f, separators=(',', ':'))
print('merged', len(fl), 'records')

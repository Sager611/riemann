# summary of gram_win -dyn on 200 isolated bad points of Jerby's Zenodo table (data as in jerby_data.py); run from JerbyConf/
import re, collections, csv
rows = [l for l in open('logs/dyn_gbg_sample.log') if l.startswith('DYNSUM')]
conf = sum(1 for l in rows if 'confined_along_diagonal=yes' in l)
coll = [l for l in rows if re.search(r'real_range=\[0,', l)]
print('sample', len(rows), 'confined along diagonal', conf, '; with some r where the real sign changes vanish (pair non-real):', len(coll))
print(collections.Counter(re.search(r'count_range=\[(-?\d+),(-?\d+)\]', l).groups() for l in rows))
bad = [int(re.search(r'n=(\d+)', l).group(1)) for l in rows if 'confined_along_diagonal=no' in l]
collided = set(int(re.search(r'n=(\d+)', l).group(1)) for l in coll)
mu = {}
with open('/tmp/jc_zen/Zonedo/250000-GramPoints_Dataset/GramPoints_CORE_1_250000.csv') as f:
    r = csv.reader(f); next(r)
    want = set(int(re.search(r'n=(\d+)', l).group(1)) for l in rows)
    for row in r:
        n = int(row[0])
        if n in want: mu[n] = float(row[4])
print('not confined along diagonal:', bad[:30])
ns = sorted(mu)
for lo, hi in [(0, 6), (6, 10), (10, 20), (20, 1e9)]:
    sel = [n for n in ns if lo <= mu[n] < hi]
    print(f'mu in [{lo},{hi}): {len(sel)} points, collide along the diagonal: {sum(1 for n in sel if n in collided)}, not confined: {sum(1 for n in sel if n in bad)}')

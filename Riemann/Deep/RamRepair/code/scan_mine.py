# scan_mine.py -- certified rescan of Proposition master (Padhi Sec. 3, Prop. 3.5 + Lemma 3.6) on every interval [L1, L2] of the
# published cover, with Padhi's own ball-arithmetic evaluator q_arb (certify/arb/scan_arb.py) and three replacements:
#   (a) N4: our density cells out/cells_final.json (far term on both sides, kap2 factor, rigorous G bound, own d4; no Thm 5.2);
#   (b) zero-free region: our certified table out/zfr_table.json (per height, from Mossinghoff-Trudgian-Yang 1/5.558691), no BTY floor;
#   (c) N3 (Bellotti's tables and uniform bounds) removed.  N1 (Hasanalizade-Shen-Wong) and N2 (Kadiri-Lumley-Ng Table 1) kept.
# usage: scan_mine.py PART NPARTS [kernel-file]; output out/scan_mine_PART.txt, lines 'L1 L2 q_upper c=.. f=.. {parts}' (resumable).
import sys, os, json, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'repo', 'ramanujans-inequality-main')
sys.path.insert(0, os.path.join(R, 'certify', 'arb')); sys.path.insert(0, os.path.join(R, 'certify'))
os.chdir(os.path.join(R, 'certify', 'arb'))
import scan_arb as SA
from flint import arb
# (a) cells
CF = json.load(open(os.path.join(HERE, 'out', 'cells_final.json')))
SA.CELL.clear()
for k, v in CF.items():
    if v is not None and v.get('logN_up') is not None: SA.CELL[k] = arb(v['logN_up'])
# (b) zero-free region: lam(l) = table value at the smallest grid point >= l; 1/5.558691 (rounded down) above the table
ZT = json.load(open(os.path.join(HERE, 'out', 'zfr_table.json')))
GRID = sorted((float(k), float(v)) for k, v in ZT['table'].items()); ABOVE = float(ZT['above'])
def lam_mine(lt):
    for L, v in GRID:
        if L >= lt - 1e-12: return arb(v)
    return arb(ABOVE)
def edge_mine(lt): return 1 - lam_mine(lt)/arb(lt)
SA.lam_at = lam_mine; SA.edge = edge_mine
# (c) no N3
SA.ROWS = []; SA.ROWLOG = []; SA.LOG_C1 = arb(10)**9; SA.LOG_C2 = arb(10)**9
part, nparts = int(sys.argv[1]), int(sys.argv[2])
kfile = sys.argv[3] if len(sys.argv) > 3 else None
ivs = []
if kfile:
    for line in open(kfile):
        p = line.split()
        if len(p) >= 4: ivs.append((float(p[0]), float(p[1]), int(p[2]), float(p[3])))
else:
    for fn in ['scan_arb_0.txt', 'scan_arb_1.txt']:
        for line in open(fn):
            p = line.split()
            ivs.append((float(p[0]), float(p[1]), int(p[3].split('=')[1]), float(p[4].split('=')[1])))
ivs.sort()
out = os.path.join(HERE, 'out', 'scan_mine_%d.txt' % part if not kfile else 'scan_mine_k%d.txt' % part)
done = set()
if os.path.exists(out):
    for line in open(out):
        p = line.split()
        if len(p) >= 3: done.add((float(p[0]), float(p[1])))
t0 = time.time()
with open(out, 'a') as fo:
    for i, (a, b, c, f) in enumerate(ivs):
        if i % nparts != part or (a, b) in done: continue
        q, parts = SA.q_arb(a, b, 3.0e12, c, f)
        fo.write('%r %r %r  c=%d f=%.2f %s\n' % (a, b, SA.fup_safe(q), c, f, json.dumps(parts))); fo.flush()
        print('%g %g q=%.6f c=%d f=%.2f %.0fs' % (a, b, SA.fup(q), c, f, time.time() - t0), flush=True)
print('done %.0fs' % (time.time() - t0), flush=True)


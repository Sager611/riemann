"""Step 5: exact certificates for the 345 two-shift mixtures of the 0.148369 history.
Per mixture (parameters read from the history's mixed_kernels table, which the auditor uses):
  * polynomial: leafcert construction A == construction B == the producer's polynomial (dbn4/certify_mix.py
    build + homogenize); degree <= 9 in each variable; identity exact at 3 rational points;
  * coefficient signs lam1, lam2 >= 0; r1 != r2 both among the 10 declared shifts; 0 < rho < 1;
  * leafcert prover records leaf boxes; its leaf/node/depth counts are compared with the producer's record;
  * certificate check: exact tiling of [0,1]^3 and all directly recomputed Bernstein coefficients >= 0;
  * the table is identical to dbn4/mixed_kernels.json.
Negative controls: lam1 raised by 1e-2 on 5 mixtures -> their boxes fail.
Outputs: EVIDENCE/05_mixtures.tsv, EVIDENCE/05_mixtures_summary.txt, EVIDENCE/certs/mix_leaves.txt"""
import collections, random, time
from fractions import Fraction as F
from common import *
import leafcert as L

guard(5400)
log = Log(EVID / '05_mixtures_summary.txt')
log('# Exact certificates for the 345 two-shift mixtures (step s5_mix_certs.py)')
H = load_json(HIST_MIX)
MK = H['mixed_kernels']
tab = load_json(DBN4 / 'mixed_kernels.json')
log('history', HIST_MIX.name, 'sha256', sha(HIST_MIX), '; mixtures', len(MK), '; equal to dbn4/mixed_kernels.json:', tab['mixed_kernels'] == MK)
assert tab['mixed_kernels'] == MK and len(MK) == 345
sys.path.insert(0, str(DBN4))
import certify_mix as CM  # pyright: ignore[reportMissingImports]
SHIFTS = {'9/20', '19/40', '49/100', '1/2', '51/100', '21/40', '107/200', '113/200', '23/40', '3/5'}
use = collections.Counter(rec['mix'] for rec in H['records'] if rec['type'] == 'mix')
rnd = random.Random(20261022)
rows, certs = [], []
t0 = time.time()
tot = 0
same_counts = 0
for m in MK:
    assert m['r1'] in SHIFTS and m['r2'] in SHIFTS and m['r1'] != m['r2']
    r1, r2, rho, l1, l2 = (F(m[x]) for x in ('r1', 'r2', 'rho', 'lam1', 'lam2'))
    assert l1 >= 0 and l2 >= 0 and 0 < rho < 1
    sh = [(r1, l1), (r2, l2)]
    a = L.poly_A(sh, rho)
    c3 = {k: F(v) for k, v in CM.homogenize(CM.build(m['r1'], m['r2'], m['rho'], m['lam1'], m['lam2']), 9).items()}
    same = a == L.poly_B(sh, rho) == c3
    deg = L.degrees(a)
    pts = []
    for case in range(3):
        y = rho * F(rnd.randint(1, 1000), 1000)
        b = F(0) if case == 0 else y * F(rnd.randint(0, 999), 1000)
        d = F(0) if case == 1 else F(rnd.randint(-4000, 4000), 1000)
        pts.append(L.point_identity(a, sh, rho, d, y, b))
    aint, den = L.to_int(a)
    c = L.dense(aint, 9)
    boxes, nodes, md = L.prove_boxes(L.bern_on_box(c, 9, ((0, 0), (0, 0), (0, 0))), 9)
    assert boxes is not None, m['id']
    ok, why, mn = L.check_certificate(a, 9, boxes)
    pc = m['cert']
    sc = (len(boxes), nodes, md) == (pc['leaves'], pc['nodes'], pc['maxdepth']) and pc['geometry']
    same_counts += sc
    assert same and max(deg) <= 9 and all(pts) and ok, (m['id'], same, deg, pts, why)
    tot += len(boxes)
    rows.append('%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s' % (
        m['id'], m['cell'], m['provider'], m['r1'], m['r2'], m['rho'], m['lam1'], m['lam2'], use.get(m['id'], 0),
        len(boxes), nodes, md, 'same' if sc else 'differs', 'tiles+all>=0', 'A=B=producer,3/3pts'))
    certs.append('mix\t%s|%s|%s|%s|%s\t%s' % (r1, r2, rho, l1, l2, L.encode_boxes(boxes)))
    if m['id'] % 50 == 0:
        print('mixture', m['id'], '%.0fs' % (time.time() - t0), flush=True)
(EVID / 'certs').mkdir(exist_ok=True)
(EVID / 'certs' / 'mix_leaves.txt').write_text(
    '# leaf boxes per mixture: kind<TAB>r1|r2|rho|lam1|lam2<TAB>boxes; box = kz.dz,kp.dp,ku.du in (z,p,u)\n'
    + '\n'.join(certs) + '\n')
(EVID / '05_mixtures.tsv').write_text(
    'id\tcell\tprovider\tr1\tr2\trho\tlam1\tlam2\tsteps_in_0.148369\tleaves\tnodes\tmax_depth\tvs_producer_record\tcoverage+positivity\tpolynomial\n'
    + '\n'.join(rows) + '\n')
log('all 345 mixtures: polynomial A == B == producer polynomial; degree <= 9; identity exact at 3 points each (1035 points)')
log('all 345 mixtures: lam1, lam2 >= 0; r1 != r2 among the declared shifts; 0 < rho < 1 (max rho %s)' % max((F(m['rho']) for m in MK)))
log('all 345 mixtures: boxes tile [0,1]^3 exactly; all directly recomputed Bernstein coefficients >= 0')
log('leaf boxes total', tot, '; leaf/node/depth counts equal to the producer\'s recorded certificate for', same_counts, 'of 345')
log('mixtures used by the history:', len(use), '; mixture steps', sum(use.values()))
log('\n## negative controls')
for m in MK[:5]:
    sh = [(F(m['r1']), F(m['lam1']) + F(1, 100)), (F(m['r2']), F(m['lam2']))]
    a = L.poly_A(sh, F(m['rho']))
    ok, why, mn = L.check_certificate(a, 9, L.decode_boxes(certs[m['id']].split('\t')[2]))
    log('  mixture %d lam1+1e-2: its boxes %s (%s)' % (m['id'], 'ACCEPT' if ok else 'reject', why[:45]))
    assert not ok
log('certificate file EVIDENCE/certs/mix_leaves.txt sha256', sha(EVID / 'certs' / 'mix_leaves.txt'))
log('RESULT all 345 mixtures certified; %.0f s' % (time.time() - t0))
log.save()

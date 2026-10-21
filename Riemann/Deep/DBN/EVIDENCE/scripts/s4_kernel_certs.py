"""Step 4: exact certificates for all 843 single kernels in the evidence sources (57 Planat + 786 appended).
Per kernel:
  * polynomial: leafcert construction A == Planat's unchanged check_geometry.poly == leafcert construction B;
    degree <= 5 in each variable; the identity 2K_y/y - lam J = 4P/(Q* D4) holds exactly at 3 rational points;
  * Planat's unchanged prove() (depth cap 30, node cap 15000): pass/fail;
    for each failure, his prove() with ONLY the depth cap replaced by the kernel's leaf depth must pass;
  * leafcert prover (same split rule) records the leaf boxes; when Planat passes, leaf and node counts agree;
  * certificate check: boxes tile [0,1]^3 exactly and every Bernstein coefficient, recomputed on each box
    directly from the power basis, is >= 0;
  * own pair: cap valid on [0, rho^2] (Planat's formula and leafcert.cap_ok); 0 < rho < 1 (so D4 > 0); 0 < lam <= 2.
Negative controls: lam raised by 1e-4 on the 5 most used kernels -> their boxes fail, Planat's prove fails;
the forged kernel of step 3 (lambda 21/10) has no certificate.
Outputs: EVIDENCE/04_single_kernels.tsv, EVIDENCE/04_single_kernels_summary.txt, EVIDENCE/certs/single_leaves.txt"""
import collections, inspect, random, time
from fractions import Fraction as F
from common import *
import leafcert as L

guard(3600)
log = Log(EVID / '04_single_kernels_summary.txt')
log('# Exact certificates for the 843 single kernels (step s4_kernel_certs.py)')
cg, hc = planat_module('check_geometry')
log('Planat code/check_geometry.py sha256', hc, '(equals his SHA256SUMS.json entry)')
S = [load_json(WORK / 'evidence' / n) for n in NAMES]
for pi in (0, 1):
    assert sha(WORK / 'evidence' / NAMES[pi]) == EVIDENCE_SHA256[pi]
assert S[0]['kernels'] == S[1]['kernels'] and len(S[0]['kernels']) == 843
log('kernels read from the hashed evidence sources (identical in both providers): 843')
H = load_json(HIST_SINGLE)
used = collections.Counter(rec['kernel'] for rec in H['records'])
src = inspect.getsource(cg.prove)
assert src.count('depth>=30') == 1


def planat_prove_with_cap(cap):
    ns = dict(cg.__dict__)
    exec(src.replace('depth>=30', 'depth>=%d' % cap), ns)
    return ns['prove']


rnd = random.Random(20261021)
rows = []
certs = []
t0 = time.time()
fails = []
tot_leaves = 0
minmargin = None
for idx, k in enumerate(S[0]['kernels']):
    r, rho, lam, M = (F(k[x]) for x in ('r', 'rho', 'lambda', 'self_cap'))
    assert 0 < rho < 1 and 0 < lam <= 2
    a = L.poly_A([(r, lam)], rho)
    same = a == cg.poly(r, rho, lam) == L.poly_B([(r, lam)], rho)
    deg = L.degrees(a)
    pts = []
    for case in range(3):
        y = rho * F(rnd.randint(1, 1000), 1000)
        b = F(0) if case == 0 else y * F(rnd.randint(0, 999), 1000)
        d = F(0) if case == 1 else F(rnd.randint(-4000, 4000), 1000)
        pts.append(L.point_identity(a, [(r, lam)], rho, d, y, b))
    pp, pleaves, pnodes, pdepth = cg.prove(cg.bernstein(cg.poly(r, rho, lam)))
    aint, den = L.to_int(a)
    c = L.dense(aint, 5)
    boxes, nodes, md = L.prove_boxes(L.bern_on_box(c, 5, ((0, 0), (0, 0), (0, 0))), 5)
    assert boxes is not None, idx
    agree = (len(boxes), nodes, md) == (pleaves, pnodes, pdepth) if pp else None
    capped = ''
    if not pp:
        p2 = planat_prove_with_cap(md)(cg.bernstein(cg.poly(r, rho, lam)))
        capped = 'cap%d:%s' % (md, 'PASS' if p2[0] and p2[1] == len(boxes) else 'FAIL')
        fails.append((idx, k['r'], k['rho'], k['lambda'], md, nodes, len(boxes), p2[0], used.get(idx, 0)))
    ok, why, mn = L.check_certificate(a, 5, boxes)
    R = r * r; qa = M; qb = M * (2 * R - 2) + 2; qc = M * (1 + R) ** 2 - 2 * (1 + R)
    cand = [F(0), rho * rho]; vx = -qb / (2 * qa)
    if 0 < vx < rho * rho: cand.append(vx)
    capP = min(qa * w * w + qb * w + qc for w in cand) >= 0
    capL = L.cap_ok(M, r, rho * rho)
    good = same and max(deg) <= 5 and all(pts) and ok and capP and capL and (pp or capped.endswith('PASS')) and agree is not False
    assert good, (idx, same, deg, pts, ok, why, capP, capL, pp, capped, agree)
    tot_leaves += len(boxes)
    rows.append('%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s' % (
        idx, k['r'], k['rho'], k['lambda'], k['self_cap'], used.get(idx, 0), len(boxes), nodes, md,
        'PASS' if pp else 'FAIL(cap30)', capped or '-', 'tiles+all>=0', 'A=P=B,3/3pts', 'ok'))
    certs.append('single\t%s|%s|%s\t%s' % (r, rho, lam, L.encode_boxes(boxes)))
    if idx % 100 == 0:
        print('kernel', idx, '%.0fs' % (time.time() - t0), flush=True)
(EVID / 'certs').mkdir(exist_ok=True)
(EVID / 'certs' / 'single_leaves.txt').write_text(
    '# leaf boxes per kernel: kind<TAB>r|rho|lambda<TAB>boxes; box = kz.dz,kp.dp,ku.du = prod [k/2^d,(k+1)/2^d] in (z,p,u)\n'
    + '\n'.join(certs) + '\n')
(EVID / '04_single_kernels.tsv').write_text(
    'index\tr\trho\tlambda\tself_cap\tsteps_in_0.148409\tleaves\tnodes\tmax_depth\tplanat_check_geometry_unchanged\tplanat_prove_depth_cap_raised\tcoverage+positivity\tpolynomial\town_pair_cap\n'
    + '\n'.join(rows) + '\n')
log('all 843 kernels: polynomial A == Planat poly == B; degree <= 5; identity exact at 3 rational points each (2529 points)')
log('all 843 kernels: certificate boxes tile [0,1]^3 (exact volume 1, disjoint interiors) and all Bernstein')
log('coefficients recomputed directly on each box are >= 0; own-pair cap valid on [0, rho^2]; 0 < rho < 1; 0 < lam <= 2')
log('total leaf boxes', tot_leaves, '(Planat 57:', sum(int(x.split('\t')[6]) for x in rows[:57]), '; appended 786:', sum(int(x.split('\t')[6]) for x in rows[57:]), ')')
log('appended kernels used by the 0.148409 history:', sum(1 for i in used if i >= 57), 'distinct, all certified; original kernels used:', sum(1 for i in used if i < 57))
log('Planat\'s unchanged prove() (depth cap 30): passes %d of 843; agrees with leafcert leaf and node counts on all of them' % (843 - len(fails)))
log('kernels needing a deeper subdivision than his cap 30 (his prove with only "depth>=30" -> "depth>=max_depth"):')
log('  %5s %7s %10s %12s %9s %6s %6s %12s %s' % ('index', 'r', 'rho', 'lambda', 'max_depth', 'nodes', 'leaves', 'passes then', 'steps used'))
for f in fails:
    log('  %5d %7s %10s %12s %9d %6d %6d %12s %d' % (f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8]))
log('  (a cap is a resource limit: raising it cannot admit a wrong kernel, since every leaf still needs all coefficients >= 0)')
log('\n## negative controls')
top = [i for i, _ in used.most_common(5)]
for idx in top:
    k = S[0]['kernels'][idx]
    r, rho, lam = F(k['r']), F(k['rho']), F(k['lambda']) + F(1, 10 ** 4)
    a = L.poly_A([(r, lam)], rho)
    boxes = L.decode_boxes(certs[idx].split('\t')[2])
    ok, why, mn = L.check_certificate(a, 5, boxes)
    pp = cg.prove(cg.bernstein(cg.poly(r, rho, lam)))[0]
    log('  kernel %d lambda+1e-4: its boxes %s (%s); Planat prove %s' % (idx, 'ACCEPT' if ok else 'reject', why[:40], 'PASS' if pp else 'fails'))
    assert not ok and not pp
k = S[0]['kernels'][67]
a = L.poly_A([(F(k['r']), F(21, 10))], F(k['rho']))
ok = L.check_certificate(a, 5, L.decode_boxes(certs[67].split('\t')[2]))[0]
pp = cg.prove(cg.bernstein(cg.poly(F(k['r']), F(k['rho']), F(21, 10))))[0]
log('  forged kernel of step 3 (kernel 67 with lambda 21/10): boxes', 'ACCEPT' if ok else 'reject', '; Planat prove', 'PASS' if pp else 'fails')
assert not ok and not pp
log('certificate file EVIDENCE/certs/single_leaves.txt sha256', sha(EVID / 'certs' / 'single_leaves.txt'))
log('RESULT all 843 single kernels certified; %.0f s' % (time.time() - t0))
log.save()

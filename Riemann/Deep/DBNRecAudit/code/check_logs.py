#!/usr/bin/env python3
# dbnrecaudit check_logs.py: consistency of every logged shell of rows 1-8 (DBNRec rows 1-4, DBNRec2 rows 4-8) with its cert:
# unresolved-cell count (bb DONE line) = lines of the dump; thr/K of bb = prep.json; CERT present, no ABORT; shells contiguous from 0;
# cert tags = logged shells; DBNRec2 row-4 shells 00-05 = DBNRec copies; rebuild chain F's source from Raymond's file + cert strings
# (mkrec.py logic re-implemented) and compare its sha256 with the pinned 6321ae69...; same for chain B (04f16da4...).
import json, re, glob, hashlib, copy, os
from fractions import Fraction as Fr
R = 'Riemann/Deep/'; out = []
def shells(d, row): return sorted(glob.glob(R + d + '/logs/r%d_*.bb.out' % row))
bad = []
for d, rows in (('DBNRec', (1, 2, 3, 4)), ('DBNRec2', (4, 5, 6, 7, 8))):
    for row in rows:
        cert = json.load(open(R + d + '/logs/cert_r%d.json' % row))
        tags = [s['tag'] for s in cert['shells']]
        files = [os.path.basename(f)[:-7] for f in shells(d, row)]
        if tags != files: bad.append((d, row, 'tags', tags, files))
        prev = Fr(0)
        for s in cert['shells']:
            tag = s['tag']; pj = json.loads(open(R + d + '/logs/' + tag + '.prep.json').readline())
            bb = open(R + d + '/logs/' + tag + '.bb.out').read()
            m = re.search(r'DONE K (\d+) .*? unres (\d+) \|', bb); hd = re.search(r'thr ([0-9.]+) Gneed', bb)
            nlines = sum(1 for _ in open(R + d + '/logs/' + tag + '.unres.txt'))
            ok = (m is not None and int(m.group(1)) == pj['K'] and int(m.group(2)) == nlines == s['cells'] and 'ABORT' not in bb and 'CERT sup' in bb
                  and abs(float(hd.group(1)) - pj['thr']) < 5e-10 and Fr(pj['u0']) == prev and pj['row'] == row)
            prev = Fr(pj['u1'])
            T = Fr(pj['T']); sq = Fr(pj['s'])
            out.append(dict(d=d, row=row, tag=tag, u=[pj['u0'], pj['u1']], K=pj['K'], unres=int(m.group(2)) if m else None, lines=nlines, thr=pj['thr'],
                            thr_minus_su0_minus_T=float(Fr(pj['thr']) - sq*Fr(pj['u0']) - T), ok=ok))
            if not ok: bad.append((d, row, tag))
        if Fr(cert['far']['uF']) != prev: bad.append((d, row, 'uF', cert['far']['uF'], str(prev)))
# DBNRec2 row 4 shells 00-05 are tranche-1 copies
for j in range(6):
    for ext in ('.prep.json', '.bb.out', '.unres.txt'):
        a = open(R + 'DBNRec/logs/r4_%02d' % j + ext, 'rb').read(); b = open(R + 'DBNRec2/logs/r4_%02d' % j + ext, 'rb').read()
        if a != b: bad.append(('r4 copy differs', j, ext))
# rebuild chain sources (mkrec.py 'source' logic) from Raymond's source and the cert strings
src_raw = open('.scratch/dbnrecaudit/dl/Proof_and_Replay__Lambda_Peer_Review__inputs__workspace__lambda_round15__flow_main_v2__MAIN_CLOCK_SOURCE.json', 'rb').read()
assert hashlib.sha256(src_raw).hexdigest().startswith('b59905d7')
chains = {}
for name, certs in (('B', [('DBNRec', r) for r in (1, 2, 3, 4)]), ('F', [('DBNRec', 1), ('DBNRec', 2), ('DBNRec', 3)] + [('DBNRec2', r) for r in (4, 5, 6, 7, 8)])):
    src = json.loads(src_raw); new = copy.deepcopy(src); rows = []
    for d, row in certs:
        c = json.load(open(R + d + '/logs/cert_r%d.json' % row)); r = new['records'][row]
        A0, B0, A1, B1 = Fr(r['first_hi']), Fr(r['last_hi']), Fr(c['A_new']), Fr(c['B_new'])
        assert A1 <= A0 and B1 <= B0 and A1 >= B1 > 0
        r['first_hi'], r['last_hi'] = c['A_new'], c['B_new']; rows.append((row, float(A0 - A1), float(B0 - B1)))
    h = hashlib.sha256(json.dumps(new).encode()).hexdigest()
    chains[name] = dict(sha256=h, rows=rows)
mk = {'B': json.load(open(R + 'DBNRec/logs/mkrec_B.json'))['sha256'], 'F': json.loads(open(R + 'DBNRec2/logs/mkrec_F.json').read())['sha256']}
ver = {}
for name, f in (('B', R + 'DBNRec/logs/EXACT_ARC_VERIFICATION_B.json'), ('F', R + 'DBNRec2/logs/EXACT_ARC_VERIFICATION_F.json')):
    v = json.load(open(f)); ver[name] = {k: v.get(k) for k in ('PASS', 'arcs_checked', 'exact_checks', 'exact_full_cleanup_endpoint', 'headline', 'all_rows_covered')}
    ver[name]['source_sha256'] = v['input_sha256']['source']; ver[name]['rebuilt_equals_verified'] = v['input_sha256']['source'] == chains[name]['sha256']
print(json.dumps(dict(nshells=len(out), all_ok=all(o['ok'] for o in out) and not bad, bad=bad, chains=chains, mkrec=mk, verifier=ver), indent=0)[:6000])
json.dump(out, open('Riemann/Deep/DBNRecAudit/logs/check_logs_shells.json', 'w'), indent=0)


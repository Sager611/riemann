#!/usr/bin/env python3
# dbnrec mkrec.py (I/O and exact rational bookkeeping only).
#  mkrec.py source cert_r<i>.json ...  : MAIN_CLOCK_SOURCE.json -> ws/dbnrec/MAIN_CLOCK_SOURCE_REC.json with first_hi/last_hi of the certified rows
#        replaced by the certified (A', B') (only if A' <= first_hi and B' <= last_hi); asserts that nothing else changed; writes make_flow_rec.py
#        (= HP141 make_flow_b2.py with the SOURCE path/hash swapped).
#  mkrec.py verifier                   : writes verify_arcs_rec.py (= HP141 verify_arcs_141.py with the source/flow paths and hashes and the pinned
#        terminal constants taken from the new flow summary; no mathematical check changed) into .scratch/dbnrec3/verify/.
import sys, json, hashlib, copy
from pathlib import Path
from flint import fmpq
WS = Path('.scratch/dbnrec3/ws'); C = Path('Riemann/Deep/DBNRec3/code')
OLD = 'lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json'; NEW = 'dbnrec/MAIN_CLOCK_SOURCE_REC.json'; FLOW = 'dbnrec/FLOW_REC.json.gz'
def Q(v): return fmpq(*map(int, v.split('/'))) if '/' in v else fmpq(int(v))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def swap(text, old, new):
    assert text.count(old) == 1, old[:60]; return text.replace(old, new)
if sys.argv[1] == 'source':
    raw = (WS/OLD).read_bytes(); assert hashlib.sha256(raw).hexdigest() == 'b59905d7b994d3424ba14f9347be75a2429ae4d4b0de2c2404971ffda1d4d536'
    src = json.loads(raw); new = copy.deepcopy(src); log = []
    for cf in sys.argv[2:]:
        c = json.load(open(cf)); i = c['row']; r = new['records'][i]
        assert int(r['row']) == i and r['mode'] == 'retained'
        A0, B0, A1, B1 = Q(r['first_hi']), Q(r['last_hi']), Q(c['A_new']), Q(c['B_new'])
        assert Q(c['A_raymond']) == A0 and Q(c['B_raymond']) == B0 and A1 >= B1 > 0
        if A1 <= A0 and B1 <= B0:
            r['first_hi'], r['last_hi'] = str(A1), str(B1); log.append(dict(row=i, dA=str(A0 - A1), dB=str(B0 - B1), dA_float=float(int((A0-A1).p))/int((A0-A1).q), dB_float=float(int((B0-B1).p))/int((B0-B1).q)))
        else: log.append(dict(row=i, skipped=True))
    for a, b in zip(src['records'], new['records']):   # nothing else changed
        for k in set(a) | set(b):
            if k not in ('first_hi', 'last_hi'): assert a.get(k) == b.get(k), (a['row'], k)
    for k in set(src) | set(new):
        if k != 'records': assert src.get(k) == new.get(k), k
    (WS/'dbnrec').mkdir(exist_ok=True); (WS/NEW).write_text(json.dumps(new))
    h = sha(WS/NEW)
    t = (C/'make_flow_b2.py').read_text()
    t = swap(t, "SOURCE = ('lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json', 'b59905d7b994d3424ba14f9347be75a2429ae4d4b0de2c2404971ffda1d4d536')",
             "SOURCE = ('%s', '%s')  # dbnrec: certified head fields (mkrec.py)" % (NEW, h))
    (C/'make_flow_rec.py').write_text(t)
    print(json.dumps(dict(source=str(WS/NEW), sha256=h, rows=log)))
else:
    summ = json.loads((WS/'dbnrec/FLOW_REC_SUMMARY.json').read_text())
    hsrc = sha(WS/NEW); hflow = sha(WS/FLOW)
    assert summ['source_sha256'] == hsrc
    t = (C/'verify_arcs_141.py').read_text()
    t = swap(t, "'source':('lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json','b59905d7b994d3424ba14f9347be75a2429ae4d4b0de2c2404971ffda1d4d536'),", "'source':('%s','%s'),  # dbnrec" % (NEW, hsrc))
    t = swap(t, "'flow':('dbnb2/FLOW_THETA141.json.gz','ea2ab380b81b16251e6096039bc0807c2d80491167399c8a4e875744f54bf32f'),  # hp141: theta141 flow", "'flow':('%s','%s'),  # dbnrec flow" % (FLOW, hflow))
    t = swap(t, "require(len(arcs)==flow['steps']==220162,", "require(len(arcs)==flow['steps']==%d," % summ['steps'])
    tt, tq = Q(summ['terminal_t']), Q(summ['terminal_q']); e = Q(summ['endpoint']); hd = Q(summ['headline'])
    t = swap(t, "Q(110081,1024000)", "Q(%d,%d)" % (int(tt.p), int(tt.q)))
    t = swap(t, "Q(7905864869759814112092869,79228162514264337593543950336)", "Q(%d,%d)" % (int(tq.p), int(tq.q)))
    t = swap(t, "Q(2130264443492687883657300572721,19807040628566084398385987584000)", "Q(%d,%d)" % (int(e.p), int(e.q)))
    t = swap(t, "Q(21510173917,200000000000)", "Q(%d,%d)" % (int(hd.p), int(hd.q)))
    t = swap(t, "'headline':'21510173917/200000000000',", "'headline':'%s'," % str(hd))
    mc = Q(summ['minimum_kernel_cap_margin'])
    t = swap(t, "Q(33041217100869987898812442879269,30948500982134506872478105600000000)", "Q(%d,%d)" % (int(mc.p), int(mc.q)))   # producer's min cap - q
    t = swap(t, "require([r['arcs']for r in rowstats]==[5120]*20+[2560]*46+[2],'all actual timecells completelyowned')",
             "require([r['arcs']for r in rowstats][:last_row]==[int(((rows[k]['end']-rows[k]['start'])/step).floor())for k in range(last_row)]and 0<rowstats[last_row]['arcs']<=int(((rows[last_row]['end']-rows[last_row]['start'])/step).floor())and all(r['arcs']==0 for r in rowstats[last_row+1:]),'all actual timecells completelyowned')  # dbnrec: same check for any stopping row")
    vd = Path('.scratch/dbnrec3/verify'); vd.mkdir(exist_ok=True); (vd/'verify_arcs_rec.py').write_text(t)
    (C/'verify_arcs_rec.py').write_text(t)
    print(json.dumps(dict(flow_sha256=hflow, source_sha256=hsrc, steps=summ['steps'], endpoint=str(e), headline=str(hd))))

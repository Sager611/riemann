#!/usr/bin/env python3
"""Make dbnb2 copies of Raymond's arc verifier and compact-witness script.

Only the hard-coded constants of the old 7/8 run are replaced (flow path/hash, initial
square, arc count, terminal pair, endpoint, headline, per-row arc counts, minimum cap
clearance). Mode A also switches off the vestigial row-55 tube check (see REPORT.md);
mode B keeps it. Every mathematical check is otherwise byte-identical.
"""
import gzip, hashlib, json, sys
from pathlib import Path
from fractions import Fraction as F

mode, flowpath, outdir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
SRC = Path('pr/Proof_and_Replay/Lambda_Peer_Review/verification/arcs')
v = (SRC/'verify_arcs_fresh.py').read_text()
c = (SRC/'compact_arc_witnesses.py').read_text()
raw = flowpath.read_bytes()
digest = hashlib.sha256(raw).hexdigest()
flow = json.loads(gzip.decompress(raw))
n = flow['steps']
counts = [0]*67
for r in flow['records']:
    counts[r['row']] += 1
tq = F(flow['terminal_q']); tt = F(flow['terminal_t']); ep = F(flow['endpoint']); hd = F(flow['headline'])
mc = F(flow['minimum_kernel_cap_margin'])
assert counts[:20] == [5120]*20 and counts[20:66] == [2560]*46, counts
last = counts[66]

def rep(text, old, new, count=1):
    assert text.count(old) == count, (old, text.count(old))
    return text.replace(old, new)

v = rep(v, "'flow':('lambda_round16/ordinary_7_8_flow/full_flow/FLOW.json.gz','a664e829a3681d96e474bcf875d80cba7b16fe8f65e050c69f01b5ea07de4737'),",
           "'flow':('dbnb2/%s','%s'),  # dbnb2: new flow" % (flowpath.name, digest))
v = rep(v, "require(plan['initial_sigma']==flow['initial_sigma']=='7/8','ordinary7/8 initialization')",
           "require(plan['initial_sigma']=='7/8' and flow['initial_sigma']=='5458637291332547/6250000000000000','dbnb2: unchanged 7/8-era plan rows; flow initialised at Lean-checked b0=0.87338196661320752')")
v = rep(v, "require(len(arcs)==flow['steps']==221010,'every expected actual arc supplied')",
           "require(len(arcs)==flow['steps']==%d,'every expected actual arc supplied')" % n)
v = rep(v, "initial=Q(9,16)+Q(1,10**12)", "initial=(2*Q(87338196661320752,10**17)-1)**2+Q(1,10**12)  # dbnb2: (2b0-1)^2+1e-12")
v = rep(v, "require(t==Q(flow['terminal_t'])==Q(22101,204800),", "require(t==Q(flow['terminal_t'])==Q(%d,%d)," % (tt.numerator, tt.denominator))
v = rep(v, "require(q==Q(flow['terminal_q'])==Q(7902377963830536747158643,2**96),", "require(q==Q(flow['terminal_q'])==Q(%d,%d)," % (tq.numerator, tq.denominator))
v = rep(v, "Q(427693072077942473661578098491,3961408125713216879677197516800)", "Q(%d,%d)" % (ep.numerator, ep.denominator))
v = rep(v, "require(endpoint<Q(flow['headline'])==Q(107964910079,10**12)<endpoint+eta,", "require(endpoint<Q(flow['headline'])==Q(%d,%d)<endpoint+eta," % (hd.numerator, hd.denominator))
v = rep(v, "[5120]*20+[2560]*46+[850]", "[5120]*20+[2560]*46+[%d]" % last)
v = rep(v, "require(minimums['cap_clearance'][0]==Q(21059999,10**12),", "require(minimums['cap_clearance'][0]==Q(%d,%d)," % (mc.numerator, mc.denominator))
v = rep(v, "'headline':'107964910079/1000000000000',", "'headline':'%s'," % flow['headline'])
if mode == 'A':
    v = rep(v, "   if i==55:\n", "   if CHECK_RETAINED_ROW55_TUBE and i==55:  # dbnb2 mode A: vestigial, see REPORT.md\n")
    v = rep(v, "COUNT=0\nCURRENT=None\n", "COUNT=0\nCURRENT=None\nCHECK_RETAINED_ROW55_TUBE=False  # dbnb2 mode A\n")
c = rep(c, " assert count==221010\n", " assert count==%d\n" % n)
outdir.mkdir(parents=True, exist_ok=True)
(outdir/'verify_arcs_b2.py').write_text(v)
(outdir/'compact_arc_witnesses_b2.py').write_text(c)
print(json.dumps({'mode': mode, 'flow': str(flowpath), 'sha256': digest, 'steps': n, 'last_row_arcs': last}))

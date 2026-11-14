#!/usr/bin/env python3
"""hp141: copy of DBNB2/evidence/patch_verifier.py (mode A) for the theta141 flow, plus one I/O-only change in the author's verifier:
the flow is STREAMED (StreamFlow: header parsed from the prefix, records decoded one at a time; a counting pass first, so len(arcs) and
arcs[0] behave as before) instead of json.loads of the whole 220k-record file (708,952,064 B max RSS in DBNB2/evidence/A_verify_run.log).
Every mathematical check of verify_arcs_fresh.py is otherwise byte-identical; the replaced constants are the run's own (from the flow).
usage: patch_verifier141.py SRC_ARCS_DIR FLOW.json.gz OUTDIR
"""
import gzip, hashlib, json, sys
from pathlib import Path
from fractions import Fraction as F
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import theta141 as T

STREAM = r'''
class _Arcs:
 def __init__(self,path):
  self.path=path;self.n=0;self.first=None
  for a in self:
   if self.n==0:self.first=a
   self.n+=1
 def __len__(self):return self.n
 def __getitem__(self,k):
  if k!=0:raise IndexError(k)
  return self.first
 def __iter__(self):
  dec=json.JSONDecoder()
  with gzip.open(self.path,'rt',encoding='utf-8')as f:
   buf=''
   while ',"records":[' not in buf:
    s=f.read(1<<16)
    if not s:raise ValueError('no records')
    buf+=s
   pos=buf.index(',"records":[')+12;eof=False
   while True:
    if not eof and len(buf)-pos<(1<<15):
     s=f.read(1<<16);eof=not s;buf=buf[pos:]+s;pos=0
    if buf[pos]==',':pos+=1
    if buf[pos]==']':return
    obj,pos=dec.raw_decode(buf,pos)
    yield obj
class StreamFlow:   # hp141: streamed flow (I/O only)
 def __init__(self,path):
  with gzip.open(path,'rt',encoding='utf-8')as f:
   buf=''
   while ',"records":[' not in buf:buf+=f.read(1<<16)
  self.header=json.loads(buf[:buf.index(',"records":[')]+'}')
  self.arcs=_Arcs(path)
 def __getitem__(self,k):return self.arcs if k=='records' else self.header[k]
'''

def rep(text, old, new, count=1):
    assert text.count(old) == count, (old, text.count(old))
    return text.replace(old, new)

def main():
    src, flowpath, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    v = (src/'verify_arcs_fresh.py').read_text()
    digest = hashlib.sha256(flowpath.read_bytes()).hexdigest()
    ns = {}; exec('import gzip, json\n' + STREAM, ns)
    flow = ns['StreamFlow'](flowpath)
    n = len(flow['records']); counts = [0]*67
    for r in flow['records']: counts[r['row']] += 1
    tq, tt, ep, hd = F(flow['terminal_q']), F(flow['terminal_t']), F(flow['endpoint']), F(flow['headline'])
    mc = F(flow['minimum_kernel_cap_margin'])
    assert n == flow['steps'] and counts[:20] == [5120]*20 and counts[20:66] == [2560]*46, counts
    assert F(flow['initial_sigma']) == T.TH
    last = counts[66]
    v = rep(v, "'flow':('lambda_round16/ordinary_7_8_flow/full_flow/FLOW.json.gz','a664e829a3681d96e474bcf875d80cba7b16fe8f65e050c69f01b5ea07de4737'),",
               "'flow':('dbnb2/%s','%s'),  # hp141: theta141 flow" % (flowpath.name, digest))
    v = rep(v, "require(plan['initial_sigma']==flow['initial_sigma']=='7/8','ordinary7/8 initialization')",
               "require(plan['initial_sigma']=='7/8' and flow['initial_sigma']=='%s','hp141: unchanged 7/8-era plan rows; flow initialised at PR 141 theta141')" % flow['initial_sigma'])
    v = rep(v, "require(len(arcs)==flow['steps']==221010,'every expected actual arc supplied')",
               "require(len(arcs)==flow['steps']==%d,'every expected actual arc supplied')" % n)
    v = rep(v, "initial=Q(9,16)+Q(1,10**12)", "initial=(2*Q(%d,%d)-1)**2+Q(1,10**12)  # hp141: (2 theta141 - 1)^2 + 1e-12" % (T.TH_NUM, T.TH_DEN))
    v = rep(v, "require(t==Q(flow['terminal_t'])==Q(22101,204800),", "require(t==Q(flow['terminal_t'])==Q(%d,%d)," % (tt.numerator, tt.denominator))
    v = rep(v, "require(q==Q(flow['terminal_q'])==Q(7902377963830536747158643,2**96),", "require(q==Q(flow['terminal_q'])==Q(%d,%d)," % (tq.numerator, tq.denominator))
    v = rep(v, "Q(427693072077942473661578098491,3961408125713216879677197516800)", "Q(%d,%d)" % (ep.numerator, ep.denominator))
    v = rep(v, "require(endpoint<Q(flow['headline'])==Q(107964910079,10**12)<endpoint+eta,", "require(endpoint<Q(flow['headline'])==Q(%d,%d)<endpoint+eta," % (hd.numerator, hd.denominator))
    v = rep(v, "[5120]*20+[2560]*46+[850]", "[5120]*20+[2560]*46+[%d]" % last)
    v = rep(v, "require(minimums['cap_clearance'][0]==Q(21059999,10**12),", "require(minimums['cap_clearance'][0]==Q(%d,%d)," % (mc.numerator, mc.denominator))
    v = rep(v, "'headline':'107964910079/1000000000000',", "'headline':'%s'," % flow['headline'])
    v = rep(v, "   if i==55:\n", "   if CHECK_RETAINED_ROW55_TUBE and i==55:  # dbnb2 mode A: vestigial, see DBNB2/REPORT.md\n")
    v = rep(v, "COUNT=0\nCURRENT=None\n", "COUNT=0\nCURRENT=None\nCHECK_RETAINED_ROW55_TUBE=False  # dbnb2 mode A\n" + STREAM)
    v = rep(v, "  result[kind]=json.loads(gzip.decompress(data)if kind=='flow'else data)",
               "  result[kind]=StreamFlow(path)if kind=='flow'else json.loads(data)  # hp141: streamed flow (I/O only)")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir/'verify_arcs_141.py').write_text(v)
    print(json.dumps({'mode': 'A', 'flow': str(flowpath), 'sha256': digest, 'steps': n, 'last_row_arcs': last}))

if __name__ == '__main__':
    main()

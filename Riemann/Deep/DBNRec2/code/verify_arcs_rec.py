#!/usr/bin/env python3
"""Fresh exact arc verification, independent of the historical producers.

Only the actual plan, source numbers and complete arc list are inputs. No
saved acceptance receipt, producer or old verifier is loaded. Exact FLINT
rationals implement two different algebraic checks:

* Joint own derivatives use Newton sums of the reciprocal conjugate roots.
* exp(-x) uses an even diagonal [6/6] Pade upper bound, with an integral
  remainder, rather than historical alternating Taylor/reciprocal-Taylor
  implementations.

This proves the finite comparison inequalities. The all-x source theorem,
kernel/own monotonicity, initial ordinary-zeta theorem and analytic global
contact argument are independent mathematical interfaces, not inferred
from the supplied row metadata.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
import time

sys.dont_write_bytecode=True
sys.path.insert(0,'/tmp/lambda_publication_flint_runtime')
from flint import fmpq as Q

OUT=Path(__file__).resolve().parent
BOUND={
 'plan':('lambda_round14/material/m8_late/M8_LATE_PLAN.json','ec346b7c1136eda38021eb6323033b3b48ed21e6f7458eeeac7cc9579db2497a'),
 'source':('dbnrec/MAIN_CLOCK_SOURCE_REC.json','6321ae690ef0760caa98808bb5fc05af5e08a0113b96046046a538e5e4109aab'),  # dbnrec
 'flow':('dbnrec/FLOW_REC.json.gz','bca7c24873395ffee60c21258d3e5b58216477e46681d93c1e07f9b5a9fceb4f'),  # dbnrec flow
}
COUNT=0
CURRENT=None
CHECK_RETAINED_ROW55_TUBE=False  # dbnb2 mode A

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

def require(condition,explanation):
 global COUNT
 COUNT+=1
 if not condition:
  raise ArithmeticError(json.dumps({'arc':CURRENT,'failure':explanation}))

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb')as f:
  for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
 return h.hexdigest()

def read_inputs(workspace):
 result={}
 for kind,(relative,digest)in BOUND.items():
  path=workspace/relative
  data=path.read_bytes()
  require(hashlib.sha256(data).hexdigest()==digest,'input bytes: '+kind)
  result[kind]=StreamFlow(path)if kind=='flow'else json.loads(data)  # hp141: streamed flow (I/O only)
 return result

def own_via_newton(q,row):
 """Power sums satisfy S_n=(u+v)S_(n-1)-uv*S_(n-2).

 For u=1/(p-sqrt(q)), v=1/(p+sqrt(q)), u+v=2p/(p²-q),
 uv=1/(p²-q), and (d/dp)^m[2p/(p²-q)]/m!=(-1)^m S_(m+1).
 Thus neither sqrt(q), binomial expansion nor Taylor-division is used.
 """
 result=Q(0)
 for h,v,w in row['rational_nodes']:
  hh=h*h
  vv=v*v
  denominator=(hh-q)*(hh-q)+vv*(2*hh+2*q+vv)
  require(hh>q and denominator>0,'whole ordinary pair has no pole')
  result+=2*w*h*(hh+vv-q)/denominator
 if row['overlay'] is None:return result
 alpha,centres=row['overlay']
 result*=alpha
 for p,coefficients in centres:
  pole=p*p-q
  require(pole>0,'joint pole outside contact')
  root_sum=2*p/pole
  root_product=1/pole
  previous=Q(2)
  current=root_sum
  scale=1/p
  for m,b in enumerate(coefficients):
   result+=2*b*scale*current*(-1 if m%2 else 1)
   previous,current=current,root_sum*current-root_product*previous
   scale*=p
 return result

def pade_negative_upper(x):
 """P6(-x)/P6(x)>=exp(-x) for x>=0.

 P6(z)=sum ((12-k)!6!)/((12)!k!(6-k)!) z^k.
 With the harmless factor665280 its ascending coefficients are below.
 The identity
 exp(-x)P6(x)-P6(-x)
 = -x^13/12! * integral_0^1 exp(-tx)t^6(1-t)^6 dt
 proves the sign; P6(x)>0. It follows by integration by parts.
 """
 require(x>=0,'nonnegative Pade argument')
 positive=Q(1)
 negative=Q(1)
 # Horner from the x^6 coefficient; all are exact integer coefficients.
 for coefficient in (42,840,10080,75600,332640,665280):
  positive=positive*x+coefficient
  negative=-negative*x+coefficient
 require(positive>0,'positive Pade denominator')
 return negative/positive

def prepare_rows(plan,source):
 require(len(plan['rows'])==len(source['records'])==67,'all67 rows supplied')
 rows=[]
 public=[]
 previous_end=Q(0)
 for i,(p,s)in enumerate(zip(plan['rows'],source['records'])):
  require(p['row']==s['row']==i,'ordered source/plan row')
  require(p['mode']==s['mode'],'same source/kernel functional units')
  start,end=map(Q,p['time'])
  cap0,cap=map(Q,p['q_interval'])
  require(start==previous_end<end,'closed thermal cells join')
  require(cap0==0<cap,'complete contact interval starts atzero')
  require([Q(s['t_lo']),Q(s['t_hi'])]==[start,end],'source and kernel have same times')
  require(list(map(Q,s['q_interval']))==[cap0,cap],'source and kernel have same cap')
  main,debit,first,last=(Q(s[k])for k in ('main_lo','native_debit','first_hi','last_hi'))
  require(first>=last>0,'source lower chord nondecreasing')
  intercept=main-debit-first
  slope=(first-last)/(end-start)
  nodes=[]
  for node in p['nodes']:
   h,v,w=(Q(node[k])for k in ('height','shift','weight'))
   require(h>0 and h*h>cap and abs(v)<=Q(39,20),'all ordinary poles and translations valid on fullcap')
   nodes.append((h,v,w))
  overlay=None
  if p.get('joint_jet_overlay')is not None:
   o=p['joint_jet_overlay']
   alpha=Q(o['alpha'])
   require(alpha>=0,'nonnegative scalar-parent coefficient')
   centres=[]
   for c in o['centres']:
    p0=Q(c['p']);bs=list(map(Q,c['b']))
    require(p0>0 and p0*p0>cap and bs,'joint centres outside fullcap')
    centres.append((p0,bs))
   overlay=alpha,centres
  rows.append({'start':start,'end':end,'cap':cap,'intercept':intercept,
               'slope':slope,'rational_nodes':nodes,'overlay':overlay})
  public.append({'row':i,'mode':p['mode'],'time':[str(start),str(end)],
    'q_interval':['0',str(cap)],'main_lo':str(main),'native_debit':str(debit),
    'first_hi':str(first),'last_hi':str(last),
    'source_lower_intercept_at_row_start':str(intercept),'source_lower_slope':str(slope),
    'nodes':p['nodes'],'joint_jet_overlay':p.get('joint_jet_overlay')})
  previous_end=end
 require(previous_end==Q(87,800),'complete thermal-plan horizon')
 (OUT/'ACTUAL_67_ROW_ARC_INPUTS.json').write_text(json.dumps(public,indent=2)+'\n')
 return rows

def main():
 global CURRENT
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--workspace',type=Path,default=Path('/workspace/lambda_publication_20261028/evidence/workspace'))
 args=parser.parse_args()
 if not __debug__:raise RuntimeError('Run with assertions enabled.')
 started=time.monotonic()
 data=read_inputs(args.workspace.resolve())
 plan,source,flow=(data[k]for k in ('plan','source','flow'))
 require(plan['initial_sigma']=='7/8' and flow['initial_sigma']=='9908119238306612133369374975000622988619121777878229799765263366946606045502455625911569003364955967645745531862430068266267717/11345012023616779600791582323373494582287250283647959965785871011074347975117040113347744503437234376916434871293982689424818700','hp141: unchanged 7/8-era plan rows; flow initialised at PR 141 theta141')
 require(flow['schema']=='MULTIPLICATIVE_PACKED_FLOW_V1','actual arc schema')
 rows=prepare_rows(plan,source)
 arcs=flow['records']
 require(len(arcs)==flow['steps']==219999,'every expected actual arc supplied')
 require(flow['step_denominator']==2048000,'actual step granularity')
 initial=(2*Q(9908119238306612133369374975000622988619121777878229799765263366946606045502455625911569003364955967645745531862430068266267717,11345012023616779600791582323373494582287250283647959965785871011074347975117040113347744503437234376916434871293982689424818700)-1)**2+Q(1,10**12)  # hp141: (2 theta141 - 1)^2 + 1e-12
 require(Q(flow['initial_q'])==initial,'strict initial square9/16+1e-12')
 require(flow['plan_sha256']==BOUND['plan'][1]and flow['source_sha256']==BOUND['source'][1],'exact plan/source identity')
 step=Q(1,2048000)
 eta=Q(1,10**12)
 t=Q(0);q=initial;last_row=-1
 rowstats=[{'row':i,'arcs':0,'minima':{},'maxima':{}}for i in range(67)]
 minimums={}
 maximums={}
 first_witness=None
 trace=OUT/'ALL_221010_ARC_INEQUALITY_WITNESSES.jsonl.gz'
 transcript=hashlib.sha256()
 with gzip.open(trace,'wt',compresslevel=1,encoding='utf-8')as stream:
  for index,a in enumerate(arcs):
   CURRENT=index
   i=a['row']
   require(isinstance(i,int)and 0<=i<67 and i in (last_row,last_row+1),'valid contiguous source-row sequence')
   row=rows[i]
   incoming,outgoing,r,lower,upper,dt=(Q(a[k])for k in ('qin','qout','rate','source_lower','own_upper','dt'))
   require(Q(a['t'])==t and incoming==q,'exact incoming state and time join')
   require(row['start']<=t<t+dt<=row['end']and dt==step,'whole arc in one exact thermalcell')
   require(0<q<=row['cap'],'whole decreasing arc in contactcap')
   source_value=row['intercept']+row['slope']*(t-row['start'])
   source_slack=source_value-lower
   require(source_slack>=0,'stored source bound rounded downward')
   own=own_via_newton(q,row)
   own_slack=upper-own
   require(own_slack>=0,'stored full signed ownbound rounded upward')
   strict_margin=2*(lower-upper)-r
   require(r>0 and strict_margin>=eta,'strict differential comparison rate')
   positive=q-(2+r*q)*dt
   require(positive==Q(a['linear_positive_lower'])and positive>0,'entire closed ideal arc positive')
   exp_upper=pade_negative_upper(r*dt)
   outgoing_upper=(q+2/r)*exp_upper-2/r
   exp_slack=outgoing-outgoing_upper
   require(exp_slack>=0 and 0<outgoing<q,'Pade upper enclosure and decreasing outwardreset')
   denominator=int(outgoing.denominator)
   require(denominator&(denominator-1)==0 and denominator<=2**96,'outgoing96bit dyadic')
   if CHECK_RETAINED_ROW55_TUBE and i==55:  # dbnb2 mode A: vestigial, see DBNB2/REPORT.md
    require(positive>=Q(289,10000)and outgoing>=Q(289,10000)
            and q<=Q(361,10000)and outgoing<=Q(361,10000),'whole incoming arc and outgoing jump fit retainedtube')
   measures={'source_downward_slack':source_slack,'own_upward_slack':own_slack,
             'rate_margin_after_required_eta':strict_margin-eta,
             'pade_outgoing_slack':exp_slack,'whole_arc_positive_floor':positive,
             'cap_clearance':row['cap']-q,'rate':r}
   witness={'arc':index,'row':i,'t':str(t),'dt':str(dt),'qin':str(q),'qout':str(outgoing),
            **{k:str(v)for k,v in measures.items()}}
   rendered=json.dumps(witness,separators=(',',':'))+'\n'
   stream.write(rendered);transcript.update(rendered.encode())
   if first_witness is None:first_witness=witness
   rs=rowstats[i];rs['arcs']+=1
   if rs['arcs']==1:rs['first_arc']=index;rs['entry_t']=str(t);rs['entry_q']=str(q)
   rs['last_arc']=index;rs['exit_t']=str(t+dt);rs['exit_q']=str(outgoing)
   for k,v in measures.items():
    if k not in minimums or v<minimums[k][0]:minimums[k]=v,index,i
    if k not in maximums or v>maximums[k][0]:maximums[k]=v,index,i
    if k not in rs['minima']or v<Q(rs['minima'][k]['value']):rs['minima'][k]={'value':str(v),'arc':index}
    if k not in rs['maxima']or v>Q(rs['maxima'][k]['value']):rs['maxima'][k]={'value':str(v),'arc':index}
   t+=dt;q=outgoing;last_row=i
   if (index+1)%20000==0:print('CHECKED',index+1,'row',i,'elapsed',round(time.monotonic()-started,2),flush=True)
 CURRENT=None
 require(t==Q(flow['terminal_t'])==Q(219999,2048000),'actual full history reaches statedterminaltime')
 require(q==Q(flow['terminal_q'])==Q(7855654131117049873834299,79228162514264337593543950336),'actual full history reaches statedterminalsquare')
 require(0<q<=Q(1,10000),'positive complete terminalstrip')
 tau=Q(flow['startup_tau'])
 require(tau==Q(1,10**14)and 0<tau<step,'positive startup before firstarc ends')
 startup=eta-Q(arcs[0]['rate'])*initial*tau
 require(startup==Q(flow['startup_gap_lower'])and startup>=Q(1,2*10**12),'strict classical-startup separation')
 endpoint=t+q/2
 require(endpoint==Q(flow['endpoint'])==Q(2128681727881580061683701397967,19807040628566084398385987584000),'entire q/2 cleanup exact identity')
 require(endpoint<Q(flow['headline'])==Q(107470962867,1000000000000)<endpoint+eta,'upward decimal ceiling')
 require(t<Q(19,125)and endpoint<Q(41,250),'complete source and wall horizons')
 require([r['arcs']for r in rowstats][:last_row]==[int(((rows[k]['end']-rows[k]['start'])/step).floor())for k in range(last_row)]and 0<rowstats[last_row]['arcs']<=int(((rows[last_row]['end']-rows[last_row]['start'])/step).floor())and all(r['arcs']==0 for r in rowstats[last_row+1:]),'all actual timecells completelyowned')  # dbnrec: same check for any stopping row
 require(minimums['cap_clearance'][0]==Q(38274296043905262789039278035519,30948500982134506872478105600000000),'minimum actual capclearance')
 report={'PASS':True,'fresh_verification':True,'saved_PASS_receipts_loaded':False,
   'producer_or_historical_checker_imported':False,'arcs_checked':len(arcs),'exact_checks':COUNT,
   'all_rows_covered':67,'input_sha256':{k:v[1]for k,v in BOUND.items()},
   'checker_sha256':sha(__file__),'arithmetic':'FLINT exact rationals; no floating evaluations',
   'independent_algorithms':{'joint_own':'Newton power sums for reciprocal conjugate roots',
     'exponential_upper':'Even diagonal [6/6] Pade with sign-definite integral remainder'},
   'startup_tau':str(tau),'startup_gap_lower':str(startup),'terminal_t':str(t),
   'terminal_q':str(q),'exact_full_cleanup_endpoint':str(endpoint),
   'headline':'107470962867/1000000000000',
   'minimums':{k:{'value':str(v[0]),'arc':v[1],'row':v[2]}for k,v in minimums.items()},
   'maximums':{k:{'value':str(v[0]),'arc':v[1],'row':v[2]}for k,v in maximums.items()},
   'first_arc_witness':first_witness,'row_statistics':rowstats,
   'all_arc_witnesses':{'path':trace.name,'bytes':trace.stat().st_size,'sha256':sha(trace),
                        'uncompressed_canonical_transcript_sha256':transcript.hexdigest()},
   'scope':'Every finite actual arc inequality, closed time/contactdomain, exact initialization, strictstartup and fullcleanup checked. Globalphysicalsource proof, continuouskernel/ownmonotonicity and ordinaryinitialtheorem remain separate mathematical gates.',
   'seconds':time.monotonic()-started}
 (OUT/'EXACT_ARC_VERIFICATION_FRESH_V1.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in ('PASS','arcs_checked','exact_checks','terminal_t','terminal_q','exact_full_cleanup_endpoint','all_arc_witnesses','seconds')},indent=2))

if __name__=='__main__':main()

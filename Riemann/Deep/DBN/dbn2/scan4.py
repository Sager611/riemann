
import subprocess, math, sys
def run(N,t,y2,w,ps,ls=None):
    args=['nice','-n','10','./floor',str(N),repr(t),repr(y2),w]+[str(p) for p in ps]+([] if ls is None else ['@'+repr(ls)])
    out=subprocess.run(args,capture_output=True,text=True).stdout
    d=dict(kv.split('=') for kv in out.split() if '=' in kv)
    return float(d['floor'])
N0=690988
P11=[2,3,5,7,11]
t=0.1692; y2=0.005
for ls in [0.8,0.9,1.0,1.1,1.2,1.4]:
    print('lamscale',ls, run(N0,t,y2,'S',P11,ls), flush=True)
def frontier(y2,ps,thr=1e-6,lo=0.15,hi=0.175,it=12):
    for _ in range(it):
        mid=(lo+hi)/2
        if run(N0,mid,y2,'S',ps)>=thr: hi=mid
        else: lo=mid
    return hi
P23=[2,3,5,7,11,13,17,19,23]
for y2 in [0.003,0.006]:
    tt=frontier(y2,P23)
    print('P23',y2,tt,tt+y2/2,flush=True)

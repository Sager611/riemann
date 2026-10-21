
import subprocess, math, sys, json
def fl(N,t,y2,w,ps):
    out=subprocess.run(['nice','-n','10','./floor',str(N),repr(t),repr(y2),w]+[str(p) for p in ps],capture_output=True,text=True).stdout
    d=dict(kv.split('=') for kv in out.split() if '=' in kv)
    return float(d['floor'])
N=690988
ps=[2,3,5,7,11]
res=[]
for y2 in [0.025,0.03,0.0327249,0.0350708,0.04,0.05]:
    lo,hi=0.12,0.17
    if fl(N,hi,y2,'S',ps)<1e-6: res.append((y2,None)); continue
    for _ in range(14):
        mid=(lo+hi)/2
        if fl(N,mid,y2,'S',ps)>=1e-6: hi=mid
        else: lo=mid
    res.append((y2,hi,hi+y2/2))
    print(y2,hi,hi+y2/2,flush=True)


import subprocess, math, sys
def fl(N,t,y2,w,ps):
    out=subprocess.run(['nice','-n','10','./floor',str(N),repr(t),repr(y2),w]+[str(p) for p in ps],capture_output=True,text=True).stdout
    d=dict(kv.split('=') for kv in out.split() if '=' in kv)
    return float(d['floor'])
N=690988
def frontier(y2,ps,w='S',thr=1e-6,lo=0.12,hi=0.175,it=13):
    if fl(N,hi,y2,w,ps)<thr: return None
    for _ in range(it):
        mid=(lo+hi)/2
        if fl(N,mid,y2,w,ps)>=thr: hi=mid
        else: lo=mid
    return hi
for ps in ([2,3,5,7,11],[2,3,5,7,11,13,17]):
    for y2 in [0.002,0.005,0.01,0.015,0.02]:
        t=frontier(y2,ps)
        print(len(ps),y2, t, None if t is None else t+y2/2, flush=True)

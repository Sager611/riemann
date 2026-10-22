import subprocess, sys, json, time
sys.setrecursionlimit(10000)
from ref import RM, Halt, m,j,x,f,v
TM=sys.argv[1] if len(sys.argv)>1 else '../repro/machines/RH_120_paper.tm'
ROLES=sys.argv[2] if len(sys.argv)>2 else 'roles120.txt'
def q_of(N):
    q=0
    while 4**(q+1)<=3*N: q+=1
    return q
def thr(N):
    q=q_of(N); return (q+1)*2**q
def beta(b,B):
    c=0
    while B>0: B=(B-1)//b; c+=1
    return c
def formula(point,regs):
    M,J,X,F,V=regs
    if point=='width': return F>=thr(J+1)
    if point=='digits':
        N=J+V+1; assert J==0 and F==N-1
        return max(N-1-beta(N,X),0)>=thr(N)
    if point=='seed':
        N=J+1; return max(N-1-beta(N,X),0)>=thr(N)
    if point=='compare': return F>X
def reference(point,regs):
    rm=RM(regs)
    try:
        if point=='width': rm.width()
        elif point=='digits': rm.digits(); rm.width()
        elif point=='seed': rm.seed(); rm.digits(); rm.width()
        elif point=='compare':
            while rm.dec(x): rm.dec(f)
            if rm.dec(f): raise Halt()
        elif point=='stage': rm.stage()
        halted=False
    except Halt: halted=True
    return halted, rm.ops, '%016x'%rm.h, tuple(rm.r)
def tm(point,regs,maxsteps=10**11,timeout=900):
    import os
    cmd=['nice','./tmsim',TM,ROLES,'--load',os.environ.get('SNAP','snap_')+'%s.bin'%point,*map(str,regs),'--maxsteps',str(maxsteps)]
    out=subprocess.run(cmd,capture_output=True,text=True,timeout=timeout).stdout.strip().splitlines()[-1]
    kv=dict(p.split('=',1) for p in out.split()[1:] if '=' in p)
    kind=out.split()[0]
    if kind=='HALT': return True, int(kv['stage_ops']), kv['hash'], None, kv
    if kind=='RESTART': return False, int(kv['stage_ops']), kv['hash'], tuple(map(int,kv['regs'].split(','))), kv
    return None, None, None, None, kv
def run(cases):
    bad=0
    for point,regs in cases:
        t0=time.time()
        fh=formula(point,regs); rh,rops,rhash,rregs=reference(point,regs)
        th,tops,thash,tregs,kv=tm(point,regs)
        ok = (th==rh==fh) and tops==rops and thash==rhash and (th or tregs==rregs)
        bad+= not ok
        print('OK ' if ok else 'BAD', point, regs, 'formula_halt=%s ref_halt=%s tm_halt=%s ops=%s/%s hash=%s regs=%s steps=%s from_state=%s %.1fs'%(fh,rh,th,rops,tops,'eq' if thash==rhash else (rhash,thash),tregs,kv.get('steps'),kv.get('from_state'),time.time()-t0), flush=True)
    return bad
if __name__=='__main__':
    which=sys.argv[3] if len(sys.argv)>3 else 'small'
    cases=[]
    if which=='small':
        for N in [2,3,4,5,6,7,8,20,21,22,23,84,85,86,87]:
            T=thr(N)
            for F in sorted({0,max(T-1,0),T,T+1,2*T+3}):
                cases.append(('width',(0,N-1,0,F,0)))
        for Q,F in [(0,0),(0,1),(1,1),(1,2),(3,3),(3,4),(7,9),(12,12),(12,13)]:
            cases.append(('compare',(5,0,Q,F,0)))
        for N,X in [(5,0),(5,1),(5,5),(5,6),(13,0),(13,1),(14,14),(14,15),(15,240),(15,241),(34,34),(34,35),(35,1260),(35,1261),(21,3),(86,1000),(342,50),(2,3),(3,35),(4,143),(5,3599)]:
            cases.append(('digits',(0,0,X,N-1,N-1)))
        for N,X in [(5,0),(14,14),(14,15),(15,241)]:
            cases.append(('seed',(0,N-1,X,0,0)))
    elif which=='mini':
        for N in [5,6,21,22,85,86,341,342]:
            T=thr(N)
            for F in [T-1,T]: cases.append(('width',(0,N-1,0,F,0)))
        for N,X in [(5,0),(5,1),(14,14),(14,15),(15,240),(15,241),(35,1260),(35,1261),(86,1000)]:
            cases.append(('digits',(0,0,X,N-1,N-1)))
        for Q,F in [(3,3),(3,4)]: cases.append(('compare',(5,0,Q,F,0)))
    elif which=='big':
        for N in [341,342,1365,1366,2656,2657]:
            T=thr(N)
            for F in [T-1,T]: cases.append(('width',(0,N-1,0,F,0)))
        for N,X in [(16,4368),(16,4369)]:
            cases.append(('digits',(0,0,X,N-1,N-1)))
    print('BAD COUNT',run(cases))

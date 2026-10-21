
"""Generate an extended fixed-kernel catalogue for the 10 declared shifts, select the kernels used by a
float replay, and certify each selected kernel exactly (certify_kernel.py). Output: catalogue.json"""
import numpy as np, math, json, sys, time
from fractions import Fraction as F
sys.path.insert(0,'.')
from replay import load, kernels, replay
from certify_kernel import certify, s_cap_ok
SHIFTS=['9/20','19/40','49/100','1/2','51/100','21/40','107/200','113/200','23/40','3/5']
def ratio(D,v,B,R):
    C=D+R+1+B; A=C*C+4*D*R-4*B
    D4=A*A-16*D*R*C*C; G=A*(C-2*B)-8*D*R*C
    Qs=(D+v-B)**2+4*D*B
    return 2*(D+v-B)*D4/(Qs*G)
zz=np.linspace(1e-7,1-1e-9,200001); DD=zz/(1-zz)
def lam_max(r,rho):
    R=r*r; v=rho*rho
    vals=ratio(DD,v,v*(1-1e-13),R); k=int(np.nanargmin(vals)); m=vals[k]
    # local refine
    lo=DD[max(k-2,0)]; hi=DD[min(k+2,len(DD)-1)]
    xs=np.linspace(lo,hi,2001); m=min(m,np.nanmin(ratio(xs,v,v*(1-1e-13),R)))
    return float(m)
def s_r(w,R): return 2*(1+R-w)/((1+R-w)**2+4*R*w)
def grid():
    g=[F(k,400) for k in range(1,80)]+[F(k,1000) for k in range(200,400)]+[F(k,10000) for k in range(4000,4801)]
    return g
def make_kernel(rs,rho,shrink=F(3,10**6)):
    r=F(rs); lm=lam_max(float(r),float(rho))
    if lm>=2-1e-9: lam=F(2)
    else:
        lam=F(math.floor((F(lm)-shrink)*10**7),10**7)
    R=float(r*r); V=float(rho*rho)
    w=np.linspace(0,V,20001); M=float(np.max(s_r(w,R)))
    Mr=F(math.ceil(M*10**9)+1,10**9)
    assert s_cap_ok(Mr,r,rho)
    return dict(r=rs,rho=str(rho),lam=lam,M=Mr,lam_max=lm)
if __name__=='__main__':
    t0=time.time()
    G=grid(); table={}
    for rs in SHIFTS:
        for rho in G:
            k=make_kernel(rs,rho)
            if F(k['lam'])>0: table[(rs,rho)]=k
    print('table',len(table),'%.1fs'%(time.time()-t0),flush=True)
    # float replay with original 57 + table; record used
    S=load(); K57=kernels(S)
    H=[[float(F(S[p]['plan']['rows'][i]['height'])) for i in range(608)] for p in (0,1)]
    modes={rs:'r_'+rs.replace('/','_') for rs in SHIFTS}
    ELL=[[{rs:float(F(S[p]['rows'][i][modes[rs]])) for rs in SHIFTS} for i in range(608)] for p in (0,1)]
    cand_all=[(rs,float(F(rs)),float(rho),float(k['lam']),float(k['M']),('T',rs,rho)) for (rs,rho),k in table.items()]
    cand_all+=[(None,k['r'],k['rho'],k['lam'],k['M'],('O',ki)) for ki,k in enumerate(K57)]
    for c in cand_all[-57:]:
        pass
    k57modes=[k['mode'] for k in K57]
    used=set()
    def br(i,j,q):
        best=(0,None)
        for p in (0,1):
            h=H[p][i]
            for (rs,r,rho,lam,M,key) in cand_all:
                if q>h*h*rho*rho or abs(h*r)>=2: continue
                ell=ELL[p][i][rs] if rs else ELL[p][i][ [s for s in SHIFTS if modes[s]==k57modes[key[1]]][0] ]
                c=2*lam/h*(ell-M/h)
                if c>best[0]: best=(c,key)
        used.add(best[1]); return best[0]
    U=replay(S,K57,rate_fn=br)
    print('float replay U=%.9f'%U,'used',len(used),'%.1fs'%(time.time()-t0),flush=True)
    sel=sorted([k for k in used if k[0]=='T'],key=lambda k:(SHIFTS.index(k[1]),k[2]))
    out=[]; nfail=0
    for (_,rs,rho) in sel:
        k=table[(rs,rho)]
        lam=k['lam']
        for attempt in range(6):
            res=certify(rs,rho,lam,k['M'])
            if res['geometry'] and res['cap']: break
            nfail+=1; lam=lam-F(1,10**5)*(2**attempt)
        else:
            raise SystemExit('cannot certify '+str((rs,rho)))
        res['lam_max_float']=k['lam_max']
        out.append(res)
    print('certified',len(out),'retries',nfail,'leaves',sum(o['leaves'] for o in out),'%.1fs'%(time.time()-t0),flush=True)
    json.dump(dict(schema='DBN4_EXTENDED_KERNELS_V1',count=len(out),kernels=out),open('catalogue.json','w'),indent=1)

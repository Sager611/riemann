
import numpy as np, math, sys, json, time
from fractions import Fraction as F
from scipy.optimize import linprog
sys.path.insert(0,'.')
from replay import load, kernels, replay
S=load(); K57=kernels(S)
SH=['9/20','19/40','49/100','1/2','51/100','21/40','107/200','113/200','23/40','3/5']
shifts=[float(F(x)) for x in SH]
modes=['r_'+x.replace('/','_') for x in SH]
H=[[float(F(S[p]['plan']['rows'][i]['height'])) for i in range(608)] for p in (0,1)]
ELL=[[[float(F(S[p]['rows'][i][m])) for m in modes] for i in range(608)] for p in (0,1)]
def Kp(d,b,y): return (y-b)/(d*d+(y-b)**2)+(y+b)/(d*d+(y+b)**2)
def J(d,b,h,r): return 0.5*(Kp(d+r*h,b,h)+Kp(d-r*h,b,h))
def s_r(w,r):
    R=r*r; return 2*(1+R-w)/((1+R-w)**2+4*R*w)
# exact-rho single probe lam_max
def ratio(D,v,B,R):
    C=D+R+1+B; A=C*C+4*D*R-4*B
    D4=A*A-16*D*R*C*C; G=A*(C-2*B)-8*D*R*C
    Qs=(D+v-B)**2+4*D*B
    return 2*(D+v-B)*D4/(Qs*G)
zz=np.linspace(1e-7,1-1e-9,40001); DD=zz/(1-zz)
def lam1(r,rho):
    v=rho*rho; return min(2.0,float(np.nanmin(ratio(DD,v,v*(1-1e-13),r*r))))
def single(i,q):
    best=0
    for p in (0,1):
        h=H[p][i]; rho=math.sqrt(q)/h
        if rho>=0.48: continue
        for si,r in enumerate(shifts):
            M=max(s_r(w,r) for w in np.linspace(0,rho*rho,60))
            c=2*lam1(r,rho)/h*(ELL[p][i][si]-M/h)
            best=max(best,c)
    return best
dgrid=np.unique(np.concatenate([np.linspace(0,4,500),np.geomspace(4,1e4,150)]))
def lpmix(i,q):
    y=math.sqrt(q); probes=[(p,si) for p in (0,1) for si in range(10) if y/H[p][i]<0.5]
    rows=[];rhs=[]
    for yy in (y,0.97*y,0.9*y,0.75*y,0.5*y,0.25*y):
        for b in np.linspace(0,yy*(1-1e-9),14):
            Kv=2*Kp(dgrid,b,yy)
            A=np.array([(yy/H[p][i])*J(dgrid,b,H[p][i],shifts[si]) for (p,si) in probes]).T
            m=~((dgrid<1e-12)&(abs(b-yy)<1e-6*yy))
            rows.append(A[m]); rhs.append(Kv[m])
    A=np.vstack(rows); bv=np.concatenate(rhs)
    obj=[]
    for (p,si) in probes:
        h=H[p][i]; r=shifts[si]
        M=max(s_r(w,r) for w in np.linspace(0,q/h/h,60))
        obj.append(-(2/h)*(ELL[p][i][si]-M/h))
    res=linprog(obj,A_ub=A,b_ub=bv,bounds=[(0,None)]*len(probes),method='highs')
    return -res.fun,[(probes[k],round(res.x[k],4)) for k in range(len(probes)) if res.x[k]>1e-5]
# get q per cell from the phase-1 exact history
import gzip
Hh=json.load(gzip.open('histories/classical_dbn4.json.gz','rt'))
qcell={}
for rec in Hh['records']:
    if rec['sub']==0: qcell[rec['cell']]=float(F(rec['q_in']))
out={}
t=time.time()
for i in list(range(0,600,25))+[590,594]:
    if i not in qcell: continue
    q=qcell[i]; s=single(i,q); m,w=lpmix(i,q)
    out[i]=(q,s,m)
    print('cell',i,'y=%.4f'%math.sqrt(q),'single=%.5f'%s,'mix=%.5f'%m,'rel gain %.2e'%((m-s)/s),w,'%.0fs'%(time.time()-t),flush=True)
json.dump(out,open('mixgain.json','w'))

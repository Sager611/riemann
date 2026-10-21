
import json, math, sys
from fractions import Fraction as F
base='/Users/adriansager/playground/riemann/.scratch/dbn3/src/anc/certificate/'
def load(kind='refined'):
    return [json.load(open(base+'inputs/'+kind+'/'+n)) for n in ['SOURCE_OLD_ENHANCED.json','SOURCE_RETUNED.json']]
def kernels(S):
    out=[]
    for k in S[0]['kernels']:
        out.append(dict(r=float(F(k['r'])),rho=float(F(k['rho'])),lam=float(F(k['lambda'])),M=float(F(k['self_cap'])),mode='r_'+k['r'].replace('/','_')))
    return out
def cand_cell(S,K,i,floor_adj=None):
    out=[]
    for pi in (0,1):
        h=float(F(S[pi]['plan']['rows'][i]['height']))
        for ki,k in enumerate(K):
            ell=float(F(S[pi]['rows'][i][k['mode']]))
            if floor_adj: ell+=floor_adj(pi,i,k)
            rate=2*k['lam']/h*(ell-k['M']/h)
            if rate>0 and abs(h*k['r'])<2:
                out.append((rate,pi,ki,h*h*k['rho']**2,h,ell))
    out.sort(reverse=True)
    return out
def replay(S,K,floor_adj=None,rate_fn=None,parts=32,q0=1.000001**2,qe=1e-8,log=None):
    dt=1.0/(4000*parts); q=q0
    for i in range(608):
        cand=None if rate_fn else cand_cell(S,K,i,floor_adj)
        for j in range(parts):
            if rate_fn:
                c=rate_fn(i,j,q)
            else:
                e=next((c for c in cand if q<=c[3]),None)
                c=e[0]
            nxt=(q+2/c)*math.exp(-c*dt)-2/c
            now=i/4000+j*dt
            if log is not None: log.append((i,j,q,c))
            if nxt<=qe:
                lag=math.log(1+c*(q-qe)/(2+c*qe))/c
                return now+lag+qe/2
            q=nxt
    return None
if __name__=='__main__':
    S=load(); K=kernels(S)
    lg=[]
    U=replay(S,K,log=lg)
    print('U float replay', U)
    # per-cell summary
    import collections
    cells=collections.OrderedDict()
    for (i,j,q,c) in lg:
        cells.setdefault(i,[]).append((q,c))
    for i in [0,50,100,150,200,250,300,350,400,450,500,550,594]:
        v=cells.get(i)
        if v: print(i, 'q_in %.5f'%v[0][0], 'rate %.4f'%v[0][1], 'y %.4f'%math.sqrt(v[0][0]))

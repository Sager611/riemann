
# Reproduce DPC large-gap quadratic forms (TeX 905-966): I, J_K, max J/I in nested polynomial spaces.
import sys, time, itertools, math
import numpy as np
from flint import arb, arb_mat, ctx
t0=time.time()
PREC=int(sys.argv[1]) if len(sys.argv)>1 else 320
RMAX=int(sys.argv[2]) if len(sys.argv)>2 else 12
DMAX=int(sys.argv[3]) if len(sys.argv)>3 else 5
PHIS=[s for s in (sys.argv[4].split(',') if len(sys.argv)>4 else ['363846723/100000000'])]
KT=40
ctx.prec=PREC
def q(s):
    a,b=(s.split('/')+['1'])[:2]; return arb(int(a))/int(b)
ell=arb(39327)/20000; th=ell*ell
E=[arb(5)/4, arb(3)/2, arb(2), arb(3)]
def edot(n): return sum((E[i]*n[i] for i in range(4)), arb(0))
def vecs(maxdeg):
    out=[]
    for d in range(maxdeg+1):
        for n in itertools.product(range(d+1),repeat=4):
            if sum(n)==d: out.append(n)
    return out
fact=[math.factorial(i) for i in range(120)]
def mfact(n): 
    p=1
    for x in n: p*=fact[x]
    return p
# --- moments of P via EGF exp(sum_m th Gamma(e.m) z^m/m!) ---
def moments(Evec, maxdeg):
    def ed(n): return sum((Evec[i]*n[i] for i in range(4)), arb(0))
    allv=vecs(maxdeg)
    G={n: th*ed(n).gamma()/mfact(n) for n in allv if sum(n)>0}
    Fh=[{(0,0,0,0):arb(1)}]
    Gh=[{n:G[n] for n in G if sum(n)==d} for d in range(maxdeg+1)]
    for d in range(1,maxdeg+1):
        acc={}
        for k in range(1,d+1):
            for m,gv in Gh[k].items():
                for n,fv in Fh[d-k].items():
                    key=tuple(m[i]+n[i] for i in range(4))
                    acc[key]=acc.get(key,arb(0))+k*gv*fv
        Fh.append({key:val/d for key,val in acc.items()})
    M={}
    for d in range(maxdeg+1):
        for n,fv in Fh[d].items():
            M[n]=th.gamma()/(th+ed(n)).gamma()*fv*mfact(n)
    return M
# sanity: e=(1,1,1,1) must give all moments 1
Mt=moments([arb(1)]*4, 6)
print("sanity e=1: max|M-1| =", max(abs(float((v-1).mid())) for v in Mt.values()))
M=moments(E, 2*DMAX)
print("moments done", len(M), "t=%.1f"%(time.time()-t0))
print(" E t1..t4:", [float(M[n].mid()) for n in [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]])
idt=vecs(DMAX); NT=len(idt); NV=RMAX+1
At=arb_mat(NT,NT)
for a,j in enumerate(idt):
    for b,k in enumerate(idt):
        At[a,b]=M[tuple(j[i]+k[i] for i in range(4))]
Av=arb_mat(NV,NV)
for r in range(NV):
    for s in range(NV): Av[r,s]=1/(th+r+s)
def chol(A):
    n=A.nrows(); L=arb_mat(n,n); piv=[]
    for i in range(n):
        s=A[i,i]
        for k in range(i): s-=L[i,k]*L[i,k]
        piv.append(s)
        if not (s>0): raise SystemExit("Cholesky failed at %d: %s"%(i,s.str(5)))
        L[i,i]=s.sqrt()
        for j in range(i+1,n):
            s=A[j,i]
            for k in range(i): s-=L[j,k]*L[i,k]
            L[j,i]=s/L[i,i]
    return L,piv
Lv,pv=chol(Av); Lt,pt=chol(At)
print("min pivot ratio v: %.3e  t: %.3e"%(min(float(p.mid()) for p in pv)/float(pv[0].mid()), min(float(p.mid()) for p in pt)/float(pt[0].mid())), "t=%.1f"%(time.time()-t0))
Lvi=Lv.inv(); Lti=Lt.inv(); LtT=Lt.transpose(); LtiT=Lti.transpose()
pairs=[]
for b,k in enumerate(idt):
    for a,l in enumerate(idt):
        if all(l[i]<=k[i] for i in range(4)):
            bc=1
            for i in range(4): bc*=math.comb(k[i],l[i])
            pairs.append((a,b,l,k,bc))
print("sparse pairs", len(pairs))
res={}
for ph in PHIS:
    phi=q(ph); c=arb.pi()*phi/2
    # per pair: 13x13 matrix W[r,s]
    Wp=[]
    cpow=[c**(2*kk+1)/fact[2*kk+1] for kk in range(KT+1)]
    for (a,b,l,k,bc) in pairs:
        al=edot(l); be=edot(k)
        g2=[(be-al+2*kk+1).gamma() for kk in range(KT+1)]
        W=arb_mat(NV,NV)
        for r in range(NV):
            g1=(th+r+al).gamma()
            g3=[(th+r+be+2*kk+1).gamma() for kk in range(KT+1)]
            coef=[(-1)**kk*cpow[kk]*g1*g2[kk]/g3[kk] for kk in range(KT+1)]
            for s in range(NV):
                W[r,s]=bc*sum((coef[kk]/(th+r+s+2*kk+1) for kk in range(KT+1)), arb(0))
        Wp.append(Lvi*W*Lvi.transpose())
    print(" W done phi=%s t=%.1f"%(ph,time.time()-t0))
    Jt=np.zeros((NV*NT,NV*NT)); maxrad=0.0
    for rp in range(NV):
        for sp in range(NV):
            B=arb_mat(NT,NT)
            for (a,b,l,k,bc),W in zip(pairs,Wp):
                B[a,b]=W[rp,sp]
            X=LtT*B*LtiT
            for a in range(NT):
                for b in range(NT):
                    x=X[a,b]; Jt[rp*NT+a,sp*NT+b]=float(x.mid()); maxrad=max(maxrad,float(x.rad()))
    print(" Jtilde done t=%.1f maxrad=%.2e"%(time.time()-t0,maxrad))
    Js=(Jt+Jt.T)/2
    np.save('/Users/adriansager/playground/riemann/.scratch/dpcaudit/Jsym_%s.npy'%ph.replace('/','_'), Js)
    fl=float((4*ell/arb.pi()).mid()); phf=float(phi.mid())
    out=[]
    for R in sorted(set([0,2,4,6,8,10,RMAX])):
        if R>RMAX: continue
        for D in range(DMAX+1):
            nt=sum(1 for j in idt if sum(j)<=D)
            ix=[r*NT+a for r in range(R+1) for a in range(nt)]
            lam=np.linalg.eigvalsh(Js[np.ix_(ix,ix)])[-1]
            out.append((R,D,len(ix),lam,1-phf+fl*lam))
            print("  R=%2d D=%d dim=%4d  max J/I=%.15f  1-phi+(4l/pi)J/I=%+.6e"%(R,D,len(ix),lam,1-phf+fl*lam))
    res[ph]=out
print("total t=%.1f"%(time.time()-t0))

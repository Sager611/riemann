
# test_lambda_promise.py -- Prop div:divide promise on genus-1 curves Y^3=F(X) over K=F_p (p = 1 mod 3):
# G = c*j (pushforward of D = P - inf_0, normalized at the ramification point P_1) is a norm
#   <=>  [P - inf_0] in lambda E(F_p).   '=>' is the paper's claim; the norm solver must succeed on it.
import sys, random, itertools, time
sys.path.insert(0, '.')
from genus1 import GF, Curve
from ffpoly import *
from nsolve import NormSolver, Frac

def primary_gen(q,p):
    n=p-1; s=0; t=n
    while t%q==0: t//=q; s+=1
    for g in range(2,p):
        w=pow(g,t,p)
        if pow(w,q**(s-1),p)!=1: return w

def run(p, nF, seed):
    random.seed(seed)
    k=GF(p,[0,1])          # k = F_p, elements are residues
    zeta=k.exp[(p-1)//3]
    omega=primary_gen(3,p)
    res={'curves':0,'points':0,'pred_div':0,'norm_ok_on_div':0,'norm_fail_on_div':0,'norm_ok_on_nondiv':0,'norm_fail_on_nondiv':0}
    for _ in range(nF):
        roots=random.sample(range(p),3)
        cur=Curve(k,roots,zeta)
        cur.zeta=zeta
        cur.sigma=lambda P,cur=cur: cur.norm((P[0],k.mul(zeta,P[1]),P[2]))
        cur.inf=[cur.norm((1,k.pow(zeta,j),0)) for j in range(3)]; cur.O=cur.inf[0]
        pts=cur.points()
        image=set(cur.lam(P) for P in pts)
        F=[1]
        for r in roots: F=pmul(F,[(-r)%p,1],p)
        x1=roots[0]
        res['curves']+=1
        for P in pts:
            # D = P - inf_0 ; pushforward j: affine P -> X - x(P); inf_j -> 1
            if P[2]==0:
                A=[1]; b=0
            else:
                xP=P[0]
                A=[(-xP)%p,1]
                b=1 if xP==x1 else 0
            # c = ((j/F^b)(x1))^{-1}
            jv=peval(A,x1,p)
            if b==1:
                Fb=pdiv_exact(F,[(-x1)%p,1],p)  # F/(X-x1)
                val=pow(peval(Fb,x1,p),p-2,p)  # (X-x1)/F at x1 = 1/F'(x1)
                c=pow(val,p-2,p)
            else:
                c=pow(jv,p-2,p)
            G=pscal(A,c,p)
            pred = P in image
            res['points']+=1; res['pred_div']+=int(pred)
            ok=False
            try:
                ns=NormSolver(p,3,zeta,omega,F,G,[1])
                h=ns.run()
                ok = (ns.norm_of(h)==Frac(ns.A,ns.B,p))
            except Exception as ex:
                ok=False
            key=('norm_ok_on_div' if ok else 'norm_fail_on_div') if pred else ('norm_ok_on_nondiv' if ok else 'norm_fail_on_nondiv')
            res[key]+=1
            if pred and not ok: print('PROMISE VIOLATION? p',p,'roots',roots,'P',P,flush=True)
            if ok and not pred: print('norm but not divisible?? p',p,roots,P,flush=True)
    return res

if __name__=='__main__':
    for p in [int(a) for a in sys.argv[1].split(',')]:
        t0=time.time()
        print('p',p,run(p,int(sys.argv[2]),p),'%.0fs'%(time.time()-t0),flush=True)


#!/usr/bin/env python3
"""alptoy item 1: function-field toy of Alpoge remark fn 3.
Sum_{h,x1..xk in F_q} t(x1)...t(xk) e_q(h(x1x2+...+x_{k-1}x_k)) = qN = M^k + R.
N = sum over the path quadric Q=0 of prod t(x_i) (exact DP, integers); M = sum t; R = sum_{h!=0} S(h)
(transfer matrices, complex floats); Parseval bound |S(h)| <= |t|_2^2 q^{(k-1)/2} |t|_inf^{k-2}."""
import numpy as np, math, sys, json
def is_prime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
def legendre(q): 
    sq={(x*x)%q for x in range(1,q)}; return np.array([0]+[1 if x in sq else -1 for x in range(1,q)],dtype=np.int64)
def polymod(a,b,q):  # a mod b over F_q, lists low->high
    a=a[:]; 
    while len(a)>=len(b) and any(a):
        if a[-1]%q==0: a.pop(); continue
        c=a[-1]*pow(b[-1],q-2,q)%q; s=len(a)-len(b)
        for i,bi in enumerate(b): a[s+i]=(a[s+i]-c*bi)%q
        a.pop()
    while a and a[-1]%q==0: a.pop()
    return a
def squarefree(f,q):
    df=[(i*f[i])%q for i in range(1,len(f))]
    a,b=f[:],df[:]
    while b and any(b): a,b=b,polymod(a,b,q)
    return len([c for c in a if c%q])>0 and len(a)==1
CURVES={'E1: y^2=x^3+x+1':[1,1,0,1],'E2: y^2=x^3+2x+3':[3,2,0,1],'C2: y^2=x^5+x+3 (genus 2)':[3,1,0,0,0,1]}
def ev(f,x,q): return sum(c*pow(x,i,q) for i,c in enumerate(f))%q
def N_exact(t,q,k):
    A=np.zeros((q,q),dtype=object)  # A[x_last, s]
    for x in range(q): A[x,0]=int(t[x])
    for _ in range(k-1):
        B=np.zeros((q,q),dtype=object)
        for x in range(q):
            for s in range(q):
                v=A[x,s]
                if v==0: continue
                for y in range(q):
                    if t[y]: B[y,(s+x*y)%q]+=v*int(t[y])
        A=B
    return int(sum(A[x,0] for x in range(q)))
def S_all(t,q,k):
    xs=np.arange(q); out=[]
    for h in range(q):
        F=np.exp(2j*np.pi*((h*np.outer(xs,xs))%q)/q)
        v=t.astype(complex)
        for _ in range(k-1): v=t*(F@v)
        out.append(v.sum())
    return np.array(out)
rows=[]
for name,f in CURVES.items():
    for q in [p for p in range(3,32) if is_prime(p)]:
        if not squarefree(f,q) or f[-1]%q==0: continue
        ch=legendre(q); fx=np.array([ev(f,x,q) for x in range(q)]); c=ch[fx]
        a=int(c.sum())  # a = sum chi(f(x)) = -(sum of inverse roots); affine count = q + a
        for tname,t in [('1+chi(f)',1+c),('chi(f)',c)]:
            M=int(t.sum()); n2=int((t*t).sum()); ninf=int(abs(t).max())
            for k in range(2,7):
                N=N_exact(t,q,k); S=S_all(t,q,k)
                lhs=S.sum(); R=q*N-M**k; Rd=S[1:].sum()
                B=n2*q**((k-1)/2)*ninf**(k-2)
                if B==0:
                    print('B=0',name,q,tname,k,n2,ninf,t.tolist()); continue
                rows.append(dict(curve=name,q=q,t=tname,k=k,a=a,M=M,N=N,qN=q*N,Mk=M**k,R=R,
                    err_lhs=abs(lhs-q*N),err_R=abs(Rd-R),Rimag=abs(Rd.imag),
                    maxS_over_B=float(np.abs(S[1:]).max()/B),R_over_qB=abs(R)/((q-1)*B),
                    R_over_sqrt=abs(R)/q**((k+3)/2)))
json.dump(rows,open(sys.argv[1],'w'),default=int)
mx=lambda key: max(r[key] for r in rows)
print('cases',len(rows),'max|LHS-qN|',mx('err_lhs'),'max|R_direct-(qN-M^k)|',mx('err_R'),'max|Im R|',mx('Rimag'))
print('max |S(h)|/ParsevalBound',mx('maxS_over_B'),' max |R|/((q-1)B)',mx('R_over_qB'))
for tn in ['1+chi(f)','chi(f)']:
    for k in range(2,7):
        sel=[r for r in rows if r['t']==tn and r['k']==k]
        sg=[ (1 if r['R']>0 else (-1 if r['R']<0 else 0))*(-1)**(k-1) for r in sel]
        sgN=[ (1 if r['qN']>0 else (-1 if r['qN']<0 else 0))*(-1)**(k-1) for r in sel]
        print(f"t={tn:9s} k={k}: n={len(sel)} max|R|/q^((k+3)/2)={max(r['R_over_sqrt'] for r in sel):.3f} "
              f"sign((-1)^(k-1)R): +{sg.count(1)} -{sg.count(-1)} 0:{sg.count(0)} | sign((-1)^(k-1)qN): +{sgN.count(1)} -{sgN.count(-1)} 0:{sgN.count(0)}")


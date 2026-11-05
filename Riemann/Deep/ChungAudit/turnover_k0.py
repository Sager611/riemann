# turnover_k0.py: Lemma 4.0 turnover index r*(z_m) at k=0 vs claim r*(z_m)=m (m in {d-2,d-1}; used for m=d-3 in Thm 5.3)
import sys, json, mpmath as mp
mp.mp.dps=40
pi=mp.pi
def Phi_and_d(u):
    s=mp.mpf(0); sd=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,30):
        a=pi*n*n*e2u
        t=a*(2*a-3)*mp.e**(-a)
        # d/du of 2e^{u/2} a(2a-3)e^{-a} with a'=2a: factor 2e^{u/2}e^{-a}[ (1/2)a(2a-3) + 2a(2a-3) + 4a^2 - 2a*a(2a-3) ]
        td=( mp.mpf(1)/2*a*(2*a-3) + 2*a*(2*a-3) + 4*a*a - 2*a*a*(2*a-3) )*mp.e**(-a)
        s+=t; sd+=td
        if abs(t)<mp.mpf(10)**(-60) and n>2: break
    return 2*mp.e**(u/2)*s, 2*mp.e**(u/2)*sd
def kappa0(u):
    P,Pd=Phi_and_d(u); return -Pd/(u*P)
def ystar(z):
    f=lambda y: 1/y**2+mp.mpf(1)/2-(2/z**2)*kappa0(mp.sqrt(2)*y/z)
    # f decreasing in y; bracket
    lo=mp.mpf('0.01'); hi=mp.mpf(1)
    while f(hi)>0: hi*=2
    return mp.findroot(f,(lo,hi),solver='anderson')
def hermite_pos_zeros(n):
    # probabilists He_n zeros: eigenvalues of Jacobi matrix with off-diag sqrt(i)
    A=mp.matrix(n,n)
    for i in range(1,n):
        A[i,i-1]=mp.sqrt(i); A[i-1,i]=mp.sqrt(i)
    ev=mp.eigsy(A)[0]
    return sorted([e for e in ev if e>mp.mpf(10)**-20])
if __name__=='__main__':
    R=json.load(open(sys.argv[1]))
    viol=[]
    for d in list(range(4,31)):
        key=f"{d},0"
        if key not in R: continue
        z=[mp.mpf(s) for s in R[key]['z']]
        yz=hermite_pos_zeros(2*d-1)  # d-1 positive zeros
        rs=[]
        for m in range(1,d):
            ys=ystar(z[m-1])
            r=1+sum(1 for yy in yz if yy<ys)  # lobe index
            rs.append(r)
        edge_ok = (rs[d-2]==d-1) and (rs[d-3]==d-2)
        m3 = rs[d-4] if d>=4 else None
        print('d',d,'r*(z_m) for m=1..d-1:',rs,' edge claim r*(z_{d-1})=d-1,r*(z_{d-2})=d-2:','PASS' if edge_ok else 'FAIL',' r*(z_{d-3})=',m3,'(Thm5.3 needs d-3=%d)'%(d-3))

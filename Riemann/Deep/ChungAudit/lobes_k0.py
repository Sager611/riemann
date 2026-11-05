# lobes_k0.py: Lemma 4.0(1)-(2) and Thm 4.1(2) checks at k=0 using Phi_0 (exact theta series)
# L_j(z)=int_{xi_{j-1}}^{xi_j} u Phi(u)|He_{2d-1}(z u/sqrt2)| du, xi_j=sqrt2 y_j/z
# Claims tested: (P1) L_1<...<L_{r*} and L_j decreasing for j>r* (peak exactly at r*=lobe of y*_turn)
#                (P2) M_r=(-1)^{d-1-r} W(xi_r) > 0 for all r (m in {d-2,d-1})
import sys, json, mpmath as mp
sys.path.insert(0,'.')
from turnover_k0 import Phi_and_d, ystar, hermite_pos_zeros
mp.mp.dps=30
def He(n,y):
    a,b=mp.mpf(1),y
    if n==0: return a
    for k in range(1,n): a,b=b,y*b-k*a
    return b
def gl(n):
    return mp.gauss_legendre if False else None
X,Wt=None,None
def glnodes(n):
    # Gauss-Legendre nodes/weights on [-1,1] via mpmath
    xs=[]; ws=[]
    for i in range(1,n+1):
        x=mp.cos(mp.pi*(i-mp.mpf(1)/4)/(n+mp.mpf(1)/2))
        for _ in range(100):
            p0,p1=mp.mpf(1),x
            for k in range(2,n+1): p0,p1=p1,((2*k-1)*x*p1-(k-1)*p0)/k
            dp=n*(x*p1-p0)/(x*x-1)
            dx=p1/dp; x-=dx
            if abs(dx)<mp.mpf(10)**(-mp.mp.dps+2): break
        xs.append(x); ws.append(2/((1-x*x)*dp*dp))
    return xs,ws
X,Wt=glnodes(40)
def integ(f,a,b):
    h=(b-a)/2; c=(a+b)/2
    return h*sum(w*f(c+h*x) for x,w in zip(X,Wt))
def lobes(d,z):
    yz=hermite_pos_zeros(2*d-1)
    xi=[mp.mpf(0)]+[mp.sqrt(2)*y/z for y in yz]
    f=lambda u: u*Phi_and_d(u)[0]*abs(He(2*d-1,z*u/mp.sqrt(2)))
    L=[]
    for j in range(1,d):
        L.append(integ(f,xi[j-1],xi[j]))
    # last lobe [xi_{d-1}, inf): split until negligible
    a=xi[d-1]; s=mp.mpf(0); step=mp.mpf('0.25')
    while True:
        part=integ(f,a,a+step); s+=part; a+=step
        if part < mp.mpf(10)**(-25)*s or a>6: break
    L.append(s)
    return L,yz
if __name__=='__main__':
    R=json.load(open(sys.argv[1]))
    D=list(map(int,sys.argv[2].split(',')))
    for d in D:
        z=[mp.mpf(s) for s in R[f"{d},0"]['z']]
        for m in range(1,d):
            L,yz=lobes(d,z[m-1])
            ys=ystar(z[m-1]); rstar=1+sum(1 for y in yz if y<ys)
            peak=1+max(range(d),key=lambda j:L[j])
            inc=all(L[j]<L[j+1] for j in range(rstar-1))
            dec=all(L[j]>L[j+1] for j in range(rstar-1,d-1))
            alt=sum(((-1)**(d-j))*L[j-1] for j in range(1,d+1))
            M=[sum(((-1)**(j-r-1))*L[j-1] for j in range(r+1,d+1)) for r in range(1,d)]
            Mpos=all(x>0 for x in M)
            print(('PASS ' if (inc and dec) else 'FAIL ')+'d=%d m=%d r*=%d argmax L=%d  P1 inc-to-r*:%s dec-after-r*:%s  M_r>0 all:%s  alt/maxL=%s'%(d,m,rstar,peak,inc,dec,Mpos,mp.nstr(alt/max(L),3)))

# omega_k0.py: Thm 4.1 Part 2 (p.20) check: outermost root m=d-1, k=0: claim 0<Omega_1<...<Omega_{d-1},
# Omega_r = int_{tau_{r-1}}^{tau_r} |W(u)| kappa_0'(u) du, W(u)=int_u^inf t Phi(t) He_{2d-1}(z t/sqrt2) dt, tau_r zeros of W
# also independent check of L_7 < L_8 at (9,0), z_8 with mp.quad
import sys, json, mpmath as mp
sys.path.insert(0,'.')
from turnover_k0 import Phi_and_d, kappa0, hermite_pos_zeros
mp.mp.dps=25
def He(n,y):
    a,b=mp.mpf(1),y
    if n==0: return a
    for k in range(1,n): a,b=b,y*b-k*a
    return b
R=json.load(open(sys.argv[1])); d=int(sys.argv[2]); m=int(sys.argv[3])
z=mp.mpf(R[f"{d},0"]['z'][m-1]); n=2*d-1
yz=hermite_pos_zeros(n); xi=[mp.mpf(0)]+[mp.sqrt(2)*y/z for y in yz]
g=lambda t: t*Phi_and_d(t)[0]*He(n,z*t/mp.sqrt(2))
# W at nodes (backward cumulative)
Wn=[None]*d; Wn_last=mp.quad(g,[xi[d-1],xi[d-1]+1,xi[d-1]+2,xi[d-1]+4])
Wn[d-1]=Wn_last
for r in range(d-2,-1,-1):
    Wn[r]=Wn[r+1]+mp.quad(g,[xi[r],xi[r+1]])
def W(u):
    # locate interval
    r=max(i for i in range(d) if xi[i]<=u)
    if r==d-1: return mp.quad(g,[u,u+1,u+2,u+4])
    return Wn[r+1]+mp.quad(g,[u,xi[r+1]])
print('W(0)=f(z_m) ~', mp.nstr(Wn[0],5), ' scale', mp.nstr(max(abs(x) for x in Wn),5))
print('signs of W at nodes xi_1..xi_{d-1}:', [int(mp.sign(Wn[r])) for r in range(1,d)])
# zeros tau_r in (xi_r, xi_{r+1}) for r=1..d-2
tau=[mp.mpf(0)]
for r in range(1,d-1):
    a,b=xi[r],xi[r+1]
    fa,fb=Wn[r],Wn[r+1]
    assert fa*fb<0, ('no sign change',r)
    t=mp.findroot(W,(a,b),solver='illinois')
    tau.append(t)
dk=lambda u: mp.diff(kappa0,u)
Om=[]
for r in range(1,d):
    a=tau[r-1]; b=tau[r] if r<=d-2 else None
    pts=[a]+[x for x in xi if x>a and (b is None or x<b)]
    if b is not None: pts.append(b)
    else: pts+= [pts[-1]+mp.mpf('0.5'), pts[-1]+1, pts[-1]+2]
    val=mp.mpf(0)
    for i in range(len(pts)-1):
        val+=mp.quad(lambda u: abs(W(u))*dk(u),[pts[i],pts[i+1]],maxdegree=6)
    Om.append(val)
print('Omega_r (r=1..d-1):',[mp.nstr(x,8) for x in Om])
print('Thm 4.1(2) Omega strictly increasing:', 'PASS' if all(Om[i]<Om[i+1] for i in range(len(Om)-1)) else 'FAIL')
A=sum(((-1)**(d-1-r))*Om[r-1] for r in range(1,d))
print('A = sum (-1)^{d-1-r} Omega_r =', mp.nstr(A,8))
# independent mp.quad lobes
f=lambda u: u*Phi_and_d(u)[0]*abs(He(n,z*u/mp.sqrt(2)))
L=[mp.quad(f,[xi[j-1],xi[j]]) for j in range(1,d)]
print('L_j (mp.quad) j=1..d-1:',[mp.nstr(x,10) for x in L])

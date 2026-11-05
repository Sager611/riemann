# lobes_k0_small.py: as lobes_k0.py but computes roots itself for given d (k=0), prints L_j table for chosen cases
import sys, json, mpmath as mp
sys.path.insert(0,'.')
from turnover_k0 import Phi_and_d, ystar, hermite_pos_zeros
from lobes_k0 import lobes
mp.mp.dps=30
g=[mp.mpf(s) for s in json.load(open(sys.argv[1]))]
D=list(map(int,sys.argv[2].split(',')))
show=sys.argv[3] if len(sys.argv)>3 else ''
for d in D:
    k=0; N=d-1
    coeffs=[mp.binomial(N,j)*g[k+1+j] for j in range(N+1)]
    rts=mp.polyroots(coeffs[::-1],maxsteps=400,extraprec=200)
    xs=sorted([mp.re(x) for x in rts],reverse=True)
    z=[mp.sqrt(-2*x) for x in xs]
    fails=[]
    for m in range(1,d):
        L,yz=lobes(d,z[m-1])
        ys=ystar(z[m-1]); rstar=1+sum(1 for y in yz if y<ys)
        peak=1+max(range(d),key=lambda j:L[j])
        inc=all(L[j]<L[j+1] for j in range(rstar-1))
        dec=all(L[j]>L[j+1] for j in range(rstar-1,d-1))
        M=[sum(((-1)**(j-r-1))*L[j-1] for j in range(r+1,d+1)) for r in range(1,d)]
        if not(inc and dec) or (m>=d-2 and rstar!=m): fails.append((m,rstar,peak,inc,dec,all(x>0 for x in M)))
        if show==f"{d},{m}":
            print('  z_m=',mp.nstr(z[m-1],15),' y*=',mp.nstr(ys,12),' Hermite zeros y_j:',[mp.nstr(y,10) for y in yz])
            print('  L_j:',[mp.nstr(x,10) for x in L])
            print('  M_r:',[mp.nstr(x,6) for x in M])
    print(('PASS ' if not fails else 'FAIL ')+'d=%d k=0: failures (m, r*, argmax L, inc-to-r*, dec-after-r*, M_r>0):'%d, fails)

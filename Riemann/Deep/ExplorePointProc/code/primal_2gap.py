
# explore40: continuum primal LP (float, sampled) for Lamzouri's small-gap criterion, plain vs 2-gap realizability-augmented.
# f(a) = sum_n c_n h_{2n}(sqrt(2pi) a/s)  (h_n orthonormal Hermite functions), fhat(u) = s sum_n c_n (-1)^n h_{2n}(sqrt(2pi) s u).
# Plain (Thm Main2): f>=0, fhat(0)=1, fhat<=0 on [lam,inf), maximise J_theta(f) = f(0)+2 int_0^theta f(a) a da ; works iff J>1.
# Augmented (hard-core set Z_0, gaps g_i>=lam, consecutive-gap realizability, U<=3 lam): fhat<=0 on [U,inf) only, and
#   G(g1)+G(g1+g2)+A(g1)-A(g2)+beta(1-g1) <= c  for g1,g2>=lam  (G = fhat on [lam,U), 0 beyond; A const beyond U; beta>=0),
#   then sum_{j!=l in Z_0} fhat <= 2c N_0 + o(N); works iff J > 1 + 2c.
import numpy as np, sys, time, threading, os, resource, json
def _wd(cap=285e6):
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > cap:
            print("MEMCAP exceeded, abort", flush=True); os._exit(3)
        time.sleep(0.1)
threading.Thread(target=_wd, daemon=True).start()
from scipy.optimize import linprog
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6

def hermfun(nmax, y):
    y=np.asarray(y,float); H=np.zeros((nmax+1,)+y.shape)
    H[0]=np.pi**-0.25*np.exp(-y*y/2)
    if nmax>=1: H[1]=np.sqrt(2)*y*H[0]
    for n in range(1,nmax):
        H[n+1]=np.sqrt(2/(n+1))*y*H[n]-np.sqrt(n/(n+1))*H[n-1]
    return H
S2P=np.sqrt(2*np.pi)
def basis_f(D, s, a):   # (len a, D+1)
    H=hermfun(2*D, S2P*np.asarray(a)/s); return H[0::2].T
def basis_fh(D, s, u):
    H=hermfun(2*D, S2P*s*np.asarray(u)); sg=(-1.0)**np.arange(D+1); return (s*H[0::2]*sg[:,None]).T

def solve(lam, theta, D, s, U=None, aug=True, dg=float(os.environ.get('DG','0.01')), verbose=False):
    nc=D+1
    al=np.arange(0,8.0,0.004); Bf=basis_f(D,s,al)
    x,wq=np.polynomial.legendre.leggauss(240); aq=(x+1)*theta/2; wq=wq*theta/2
    Jvec=basis_f(D,s,[0.0])[0]+2*(basis_f(D,s,aq)*(wq*aq)[:,None]).sum(0)
    fh0=basis_fh(D,s,[0.0])[0]
    if not aug: U=lam
    uu=np.arange(U,7.0,0.002); Bu=basis_fh(D,s,uu)
    rows=[]; rhs=[]
    # variables: c (nc), then if aug: A at g-grid nodes (ng), beta, cc
    if aug:
        g=np.arange(lam, U+1e-12, dg); g=np.unique(np.r_[g,U]); ng=len(g)
        nv=nc+ng+2
    else:
        ng=0; nv=nc
    import scipy.sparse as sp
    def pad(M):
        return sp.hstack([sp.csr_matrix(M), sp.csr_matrix((M.shape[0], nv-nc))]).tocsr()
    A_ub=[pad(-Bf)]; b_ub=[np.zeros(len(al))]           # f>=0
    A_ub.append(pad(Bu)); b_ub.append(np.zeros(len(uu)))  # fhat<=0 on [U,7)
    if aug:
        Bg=basis_fh(D,s,g)           # G at nodes (g<U); at g=U index: G=0 handled
        Gmask=(g<U-1e-12).astype(float)
        iA=nc; ib=nc+ng; icc=nc+ng+1
        dense=np.zeros((ng*ng,nc)); ri=[]; ci=[]; vv=[]
        for i1 in range(ng):
            g1=g[i1]
            s12=g1+g
            m=s12<U-1e-12
            B12=np.zeros((ng,nc))
            if m.any(): B12[m]=basis_fh(D,s,s12[m])
            dense[i1*ng:(i1+1)*ng]=Gmask[i1]*Bg[i1][None,:]+B12
            for i2 in range(ng):
                rr=i1*ng+i2
                if i2!=i1:
                    ri+= [rr,rr]; ci+=[i1, i2]; vv+=[1.0,-1.0]
                ri+=[rr,rr]; ci+=[ng, ng+1]; vv+=[1-g1, -1.0]
        sparse_part=sp.csr_matrix((vv,(ri,ci)),shape=(ng*ng,ng+2))
        A_ub.append(sp.hstack([sp.csr_matrix(dense), sparse_part]).tocsr()); b_ub.append(np.zeros(ng*ng))
        del dense
    A_ub=sp.vstack(A_ub).tocsr(); b_ub=np.concatenate(b_ub)
    A_eq=np.zeros((1,nv)); A_eq[0,:nc]=fh0; b_eq=[1.0]
    cobj=np.zeros(nv); cobj[:nc]=-Jvec
    bounds=[(None,None)]*nc
    if aug:
        cobj[icc]=2.0
        bounds+= [(None,None)]*ng+[(0,None),(None,None)]
        A_eq=np.hstack([A_eq]); 
        # fix gauge A(U)=0
        row=np.zeros((1,nv)); row[0,iA+ng-1]=1.0; A_eq=np.vstack([A_eq,row]); b_eq=[1.0,0.0]
    t0=time.time()
    r=linprog(cobj,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=bounds,method='highs')
    if r.status!=0: return None, r.message
    cvec=r.x[:nc]; J=Jvec@cvec
    out={'lam':lam,'D':D,'s':s,'U':U,'aug':aug,'J':J,'obj':-r.fun,'time':time.time()-t0,'rows':A_ub.shape[0]}
    if aug: out.update({'c':r.x[icc],'beta':r.x[ib],'g':g.tolist(),'A':r.x[iA:iA+ng].tolist()})
    out['coef']=cvec.tolist()
    return out, None

if __name__=='__main__':
    theta=eval(sys.argv[1]); D=int(sys.argv[2]); s=float(sys.argv[3]); mode=sys.argv[4]
    lams=[float(x) for x in sys.argv[5].split(',')]; Uf=float(sys.argv[6]) if len(sys.argv)>6 else 2.99
    for lam in lams:
        aug=(mode=='aug')
        o,err=solve(lam,theta,D,s,U=Uf*lam,aug=aug)
        if o is None: print('lam',lam,mode,'FAIL',err,flush=True); continue
        print(f"{mode} theta={theta:.5f} D={D} s={s} lam={lam} U={o['U']:.4f} J={o['J']:.6f} obj(J-2c or J)={o['obj']:.6f}  margin={o['obj']-1:+.6f}"+(f" c={o['c']:.5f} beta={o['beta']:.4f}" if aug else '')+f" rows={o['rows']} ({o['time']:.1f}s) rss={rss():.0f}MB",flush=True)
        json.dump(o,open(f"logs/primal_{mode}_D{D}_lam{lam}.json","w"))

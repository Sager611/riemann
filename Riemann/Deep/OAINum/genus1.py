
# genus1.py -- Test of Section 5 (geometry) of the factorization paper on C: Y^3 = F(X), deg F = 3
# (genus 1, q = N = 3, r = 1, m = 3, k = cubic extension of K).
# Checks: (a) Frobenius lemma: every lambda^m-division of [Ne] and every lift chain d_{i,t} is k-rational;
# (b) label formula (Lemma geom:ramification) via explicit functions (Miller) vs. brute-force T_1 decomposition;
# (c) Prop geom:separation: for EVERY admissible d_m some reached test has non-constant labels.
import sys, random, itertools

def factor_int(n):
    fs=[]; d=2
    while d*d<=n:
        while n%d==0: fs.append(d); n//=d
        d+=1
    if n>1: fs.append(n)
    return sorted(set(fs))

class GF:
    def __init__(self, p, modpoly):
        self.p=p; self.d=len(modpoly)-1; self.mod=modpoly; self.Q=p**self.d
        d=self.d; Q=self.Q
        self.pw=[p**i for i in range(d)]
        # find generator
        fs=factor_int(Q-1)
        for cand in range(2,Q):
            v=self.dec(cand)
            if all(self.vpow(v,(Q-1)//r)!=[1]+[0]*(d-1) for r in fs):
                g=v; break
        self.exp=[0]*(Q-1); self.log={}
        cur=[1]+[0]*(d-1)
        for i in range(Q-1):
            e=self.enc(cur); self.exp[i]=e; self.log[e]=i
            cur=self.vmul(cur,g)
        assert len(self.log)==Q-1
    def enc(self,v): return sum(c*w for c,w in zip(v,self.pw))
    def dec(self,a):
        v=[]
        for _ in range(self.d): v.append(a%self.p); a//=self.p
        return v
    def vmul(self,a,b):
        p=self.p; d=self.d
        r=[0]*(2*d-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b):
                    r[i+j]+=x*y
        for k in range(2*d-2,d-1,-1):
            c=r[k]%p
            if c:
                for j in range(d): r[k-d+j]-=c*self.mod[j]
            r[k]=0
        return [x%p for x in r[:d]]
    def vpow(self,a,e):
        r=[1]+[0]*(self.d-1)
        while e:
            if e&1: r=self.vmul(r,a)
            e>>=1
            if e: a=self.vmul(a,a)
        return r
    def mul(self,a,b):
        if a==0 or b==0: return 0
        return self.exp[(self.log[a]+self.log[b])%(self.Q-1)]
    def inv(self,a): return self.exp[(-self.log[a])%(self.Q-1)]
    def add(self,a,b):
        p=self.p; r=0; w=1
        while a or b:
            r+=((a%p+b%p)%p)*w; a//=p; b//=p; w*=p
        return r
    def neg(self,a):
        p=self.p; r=0; w=1
        while a:
            r+=((-(a%p))%p)*w; a//=p; w*=p
        return r
    def sub(self,a,b): return self.add(a,self.neg(b))
    def const(self,c): return c%self.p
    def pow(self,a,e):
        if a==0: return 0 if e else 1
        return self.exp[(self.log[a]*e)%(self.Q-1)]

class Curve:
    # Y^3 = X^3 + f2 X^2 Z + f1 X Z^2 + f0 Z^3 over field k (GF); F = (X-x1)(X-x2)(X-x3)
    def __init__(self,k,roots,zeta):
        self.k=k; p=k.p
        x1,x2,x3=roots
        self.roots=roots
        self.f2=(-(x1+x2+x3))%p; self.f1=(x1*x2+x1*x3+x2*x3)%p; self.f0=(-(x1*x2*x3))%p
        self.zeta=zeta
        self.inf=[self.norm((1,k.pow(zeta,j),0)) for j in range(3)]
        self.O=self.inf[0]
    def C(self,P):
        k=self.k; X,Y,Z=P
        X2=k.mul(X,X); Z2=k.mul(Z,Z)
        t=k.mul(k.mul(Y,Y),Y)
        t=k.sub(t,k.mul(X2,X))
        t=k.sub(t,k.mul(k.const(self.f2),k.mul(X2,Z)))
        t=k.sub(t,k.mul(k.const(self.f1),k.mul(X,Z2)))
        t=k.sub(t,k.mul(k.const(self.f0),k.mul(Z2,Z)))
        return t
    def norm(self,P):
        k=self.k; X,Y,Z=P
        if Z: i=k.inv(Z); return (k.mul(X,i),k.mul(Y,i),1)
        if X: i=k.inv(X); return (1,k.mul(Y,i),0)
        assert Y; return (0,1,0)
    def lin(self,P,a,Q,b):
        k=self.k
        return tuple(k.add(k.mul(a,P[i]),k.mul(b,Q[i])) for i in range(3))
    def grad(self,P):
        k=self.k; X,Y,Z=P
        f2,f1,f0=k.const(self.f2),k.const(self.f1),k.const(self.f0)
        # C = Y^3 - X^3 - f2 X^2 Z - f1 X Z^2 - f0 Z^3
        three=k.const(3); two=k.const(2)
        dX=k.neg(k.add(k.add(k.mul(three,k.mul(X,X)),k.mul(two,k.mul(f2,k.mul(X,Z)))),k.mul(f1,k.mul(Z,Z))))
        dY=k.mul(three,k.mul(Y,Y))
        dZ=k.neg(k.add(k.add(k.mul(f2,k.mul(X,X)),k.mul(two,k.mul(f1,k.mul(X,Z)))),k.mul(three,k.mul(f0,k.mul(Z,Z)))))
        return (dX,dY,dZ)
    def third(self,P,Q):
        k=self.k
        inv2=k.inv(k.const(2)); one=1; m1=k.neg(1)
        if P!=Q:
            cp=self.C(self.lin(P,one,Q,one)); cm=self.C(self.lin(P,one,Q,m1))
            c1=k.mul(k.add(cp,cm),inv2); c2=k.mul(k.sub(cp,cm),inv2)
            return self.norm(self.lin(P,c1,Q,k.neg(c2)))
        g=self.grad(P)
        # pick T on tangent line, not proportional to P
        cands=[(0,g[2],k.neg(g[1])),(g[2],0,k.neg(g[0])),(g[1],k.neg(g[0]),0)]
        T=None
        for c in cands:
            if c==(0,0,0): continue
            # check not proportional to P
            if self.norm(c)!=P: T=c; break
        assert T is not None
        cp=self.C(self.lin(P,one,T,one)); cm=self.C(self.lin(P,one,T,m1))
        c1=k.mul(k.add(cp,cm),inv2); c0=k.mul(k.sub(cp,cm),inv2)
        return self.norm(self.lin(P,c0,T,k.neg(c1)))
    def add(self,P,Q): return self.third(self.O,self.third(P,Q))
    def negp(self,P): return self.third(self.third(self.O,self.O),P)
    def sub(self,P,Q): return self.add(P,self.negp(Q))
    def mulint(self,n,P):
        if n<0: return self.mulint(-n,self.negp(P))
        R=self.O
        for _ in range(n): R=self.add(R,P)
        return R
    def sigma(self,P):
        k=self.k; return self.norm((P[0],k.mul(k.const(self.zeta),P[1]),P[2]))
    def sigma_star(self,P): return self.sub(self.sigma(P),self.inf[1])
    def lam(self,P): return self.sub(P,self.sigma_star(P))
    def points(self):
        k=self.k
        cubes={}
        for y in range(k.Q): cubes.setdefault(k.mul(k.mul(y,y),y),[]).append(y)
        pts=list(self.inf)
        for x in range(k.Q):
            v=k.add(k.mul(k.mul(x,x),x),k.add(k.mul(k.const(self.f2),k.mul(x,x)),k.add(k.mul(k.const(self.f1),x),k.const(self.f0))))
            for y in cubes.get(v,[]): pts.append((x,y,1))
        return pts
    def divclass(self,W):
        # W = (a0,a1,a2) degree-zero infinity divisor -> group element sum a_j inf_j
        R=self.O
        for j,a in enumerate(W): R=self.add(R,self.mulint(a,self.inf[j]))
        return R

def lam_lattice_inv(W):
    # solve (1-sigma)W' = W in Lambda: b_j - b_{j-1} = a_j ; sum b = 0
    a0,a1,a2=W
    pre=[0,a1,a1+a2]   # b_j - b_0
    s=sum(pre)
    if s%3: return None
    b0=-s//3
    return (b0+pre[0],b0+pre[1],b0+pre[2])

def lam_lattice(W):
    a=W; return (a[0]-a[2],a[1]-a[0],a[2]-a[1])

# ---- Miller-type construction of h with given principal divisor; values at P_i via leading coefficients in Y
class LinForm:
    def __init__(self,c): self.c=c   # (alpha,beta,gamma): alpha X + beta Y + gamma Z

def line_through(cur,P,Q):
    k=cur.k
    if P==Q:
        g=cur.grad(P); return g
    # cross product
    return (k.sub(k.mul(P[1],Q[2]),k.mul(P[2],Q[1])), k.sub(k.mul(P[2],Q[0]),k.mul(P[0],Q[2])), k.sub(k.mul(P[0],Q[1]),k.mul(P[1],Q[0])))

def lead_at_ram(cur,L,xi):
    # leading coeff and valuation of L(X,Y,1) at P=(xi,0) in uniformizer Y; X = xi + Y^3/F'(xi) + ...
    k=cur.k
    al,be,ga=L
    c0=k.add(k.mul(al,xi),ga)
    if c0: return (0,c0)
    if be: return (1,be)
    x1,x2,x3=[k.const(r) for r in cur.roots]
    Fp=1
    for r in (x1,x2,x3):
        if r!=xi: Fp=k.mul(Fp,k.sub(xi,r))
    return (3,k.mul(al,k.inv(Fp)))

def h_values(cur,divpts,xs):
    # divpts: list of (point,mult) summing to principal divisor. Build h = prod of g_{T,Q}^{+-1}
    # g_{P,Q} = L_{P,Q}/L_{O,P*Q}: div = P + Q - (P+Q) - O.  Track (valuation, leading coeff) at each P_i.
    k=cur.k
    acc={xi:(0,1) for xi in xs}
    def mulL(L,sgn):
        for xi in xs:
            v,c=lead_at_ram(cur,L,xi); v0,c0=acc[xi]
            acc[xi]=(v0+sgn*v, k.mul(c0, c if sgn>0 else k.inv(c)))
    # express D = sum n_P P (deg 0) as sum n_P (P - O) ; handle negative via (-P) - O = -(P - O) + div(g_{P,-P})
    pts=[]
    for P,n in divpts:
        if P==cur.O: continue
        if n>0: pts+= [P]*n
        elif n<0:
            nP=cur.negp(P)
            for _ in range(-n):
                # -(P-O) = (nP - O) - div(g_{P,nP})  => divisor accumulates; h gets factor g^{-1}
                pts.append(nP)
                L1=line_through(cur,P,nP); R=cur.third(P,nP); L2=line_through(cur,cur.O,R)
                mulL(L1,-1); mulL(L2,+1)
    # now D = sum (Q-O) - div(prod g^{-1} parts) ... sequential addition
    T=cur.O
    for Qp in pts:
        if T==cur.O:
            T=Qp; continue
        L1=line_through(cur,T,Qp); R=cur.third(T,Qp); L2=line_through(cur,cur.O,R)
        # div(L1/L2) = T + Qp - (T+Qp) - O ; sum (Q-O) = (T_last - O) + div(prod g)
        mulL(L1,+1); mulL(L2,-1)
        T=cur.add(T,Qp)
    assert T==cur.O, 'divisor not principal'
    # We accumulated factors so that h = prod(...) has div(h) = D: check valuations zero later
    return acc

def run(p, triples, Kdeg=1, seed=0, verbose=False):
    random.seed(seed)
    # field k = F_p^{3*Kdeg}; zeta in K
    d=3*Kdeg
    # modulus: random irreducible of degree d over F_p (test by gcd-free brute force using GF construction attempt)
    while True:
        modp=[random.randrange(p) for _ in range(d)]+[1]
        try:
            k=GF(p,modp)
            break
        except Exception:
            continue
    # zeta: primitive cube root of unity in k (lies in K)
    zeta=k.exp[(k.Q-1)//3]
    # verify zeta in K = F_{p^Kdeg}: zeta^(p^Kdeg) == zeta
    assert k.pow(zeta,p**Kdeg)==zeta
    stats={'F':0,'rational27':0,'chain_rational':0,'sep_ok':0,'sep_fail':0,'label_formula_ok':0,'label_formula_bad':0,'first_t_hist':{}}
    for roots in triples:
        cur=Curve(k,[r%p for r in roots],zeta)
        # need zeta as an element of k; Curve.sigma uses k.const(zeta) -> fix: store raw
        cur.zeta=zeta
        cur.sigma=lambda P,cur=cur: cur.norm((P[0],k.mul(zeta,P[1]),P[2]))
        cur.inf=[cur.norm((1,k.pow(zeta,j),0)) for j in range(3)]; cur.O=cur.inf[0]
        pts=cur.points()
        lamt={P:cur.lam(P) for P in pts}
        Ne=cur.mulint(3,cur.inf[1])     # [3e] with e = inf_1 - inf_0
        sols=[P for P in pts if lamt[lamt[lamt[P]]]==Ne]
        stats['F']+=1
        if len(sols)==27: stats['rational27']+=1
        elif verbose: print('only',len(sols),'rational lambda^3-divisions',roots)
        # lift chains of P_i: all d with lambda^2 d = P_i rational (9 each)
        Pis=[cur.norm((k.const(r),0,1)) for r in roots]
        ok_chain=all(sum(1 for P in pts if lamt[lamt[P]]==Pi)==9 for Pi in Pis)
        stats['chain_rational']+=int(ok_chain)
        W={1:(-3,3,0)}            # W_1 = 3e = 3 inf_1 - 3 inf_0
        W[2]=lam_lattice_inv(W[1]); W[3]=lam_lattice_inv(W[2])
        assert W[3] is not None and lam_lattice_inv(W[3]) is None
        # epsilon_i classes
        for d3 in sols:
            dch={3:d3,2:lamt[d3],1:lamt[lamt[d3]]}
            first=None
            for t in (1,2,3):
                if lamt[dch[t]]!=cur.divclass(W[t]):
                    break  # invariant not reached (should not happen before first non-constant test)
                # brute-force T_1 labels
                labs=[]
                for a in itertools.product(range(3),repeat=3):
                    Wd=tuple(W[t][j]-(sum(a) if j==1 else (-sum(a) if j==0 else 0)) for j in range(3))
                    Wp=lam_lattice_inv(Wd)
                    if Wp is None: continue
                    S=dch[t]
                    for ai,Pi in zip(a,Pis): S=cur.sub(S,cur.mulint(ai,Pi))
                    if S==cur.divclass(Wp): labs.append(a)
                assert len(labs)==3, ('T1 decomposition count',len(labs))
                const = any(a[0]==a[1]==a[2] for a in labs)
                # function-based labels: D = d_t - O ; div h = (1-sigma)D - W_t
                dt=dch[t]
                divpts=[]
                if dt!=cur.O:
                    divpts+=[(dt,1),(cur.O,-1),(cur.sigma(dt),-1),(cur.inf[1],1)]
                for j in range(3): divpts.append((cur.inf[j],-W[t][j]))
                xs=[k.const(r) for r in roots]
                acc=h_values(cur,divpts,xs)
                # h is unit at P_i
                vals=[acc[xi] for xi in xs]
                if all(v==0 for v,_ in vals):
                    hv=[c for _,c in vals]
                    g=hv[0]
                    logs=[]
                    for c in hv:
                        r=k.mul(c,k.inv(g)); lg=None
                        for e in range(3):
                            if k.pow(zeta,e)==r: lg=e
                        logs.append(lg)
                    if None in logs: stats['label_formula_bad']+=1
                    else:
                        vP=[1 if dt==Pi else 0 for Pi in Pis]
                        lf=tuple((vP[i]-logs[i])%3 for i in range(3))
                        if any(all((lf[i]-a[i])%3==(lf[0]-a[0])%3 for i in range(3)) for a in labs): stats['label_formula_ok']+=1
                        else: stats['label_formula_bad']+=1
                else:
                    stats['label_formula_bad']+=1
                if not const:
                    first=t; break
            if first is None: stats['sep_fail']+=1
            else:
                stats['sep_ok']+=1; stats['first_t_hist'][first]=stats['first_t_hist'].get(first,0)+1
    return stats

if __name__=='__main__':
    p=int(sys.argv[1]); mode=sys.argv[2]; Kdeg=int(sys.argv[3]) if len(sys.argv)>3 else 1
    rng=random.Random(p)
    allt=[t for t in itertools.combinations(range(p),3)]
    if mode=='all': tri=allt
    else: tri=rng.sample(allt,min(int(mode),len(allt)))
    print('p',p,'Kdeg',Kdeg,'#F',len(tri),run(p,tri,Kdeg),flush=True)


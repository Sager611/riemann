# P1-P4: Theorem 12.5.1.1 (Masullo p.519) with a WIDE strip (all zeros) = Connes' form of Weil's explicit formula.
# Checks the constants 1/2 (real), 1/(2 pi) (complex), log q_v/mu_v(O_v^x) (finite) and the normalisation of M(phi)
# for the trivial character of Q, Q(i), Q(sqrt-3), Q(x^3-x-1). phi(g) = f(log|g|), f smooth bump on (A,B), A>0,
# so phi(C_F^1)=0. Measure on C_F: vol(C_F^1)=1, dt/t on R>0 (Connes' normalisation).
import cypari2, sys
pari = cypari2.Pari(); pari.allocatemem(200*10**6); pari.set_real_precision(38)
A, B, H = 0.3, 2.3, 260   # support of f; zeros up to height H
gp = pari
gp(f"A={A}; B={B}; c=(A+B)/2; w=(B-A)/2;")
gp("f(u)=my(x=(u-c)/w, y=1-x^2); if(y<1/200,0,exp(-1/y));")
gp("TAB=intnumgaussinit(700); Fm(z)=intnumgauss(u=A,B,f(u)*exp(z*u),TAB);")  # M(phi)(|.|^z), Gauss-Legendre 700 nodes
gp("Re_t=intnumgauss(u=A,B,f(u)*exp(u/2)*(1/(exp(u)-1)+1/(exp(u)+1)),TAB)/2;")  # real place
gp("Cx_t=intnumgauss(u=A,B,f(u)*exp(u/2)/(exp(u)-1),TAB);")                     # complex place
def run(name, pol, facs):
    # zeros of zeta_F from its factors (lfunzeros on the product misses close pairs: dbg2, 278 vs 114+170 for Q(i))
    gp(f"nf=nfinit({pol}); r=nf.sign; Z=concat([lfunzeros(Lf,{H}) | Lf<-[{facs}]]); Zd=lfunzeros(lfuncreate({pol}),{H});")
    gp("fin=0; forprime(p=2,exp(B)+1, P=idealprimedec(nf,p); for(j=1,#P, q=p^P[j].f; m=1; while(m*log(q)<B, fin+=log(q)*q^(-m/2)*f(m*log(q)); m++)));")
    gp("lhs=r[1]*Re_t + r[2]*Cx_t + fin;")
    gp("zs=sum(k=1,#Z, 2*real(Fm(I*Z[k])));")      # zeros come in pairs +-gamma (trivial char, real field data)
    gp("rhs=Fm(1/2)+Fm(-1/2)-zs;")
    out = gp("[#Z, lhs, rhs, lhs-rhs]")
    print(f"{name:14s} zeros<{H}: {out[0]} (product lfunzeros: {gp('#Zd')})  LHS={float(out[1]):.12f}  RHS={float(out[2]):.12f}  diff={float(out[3]):.2e}")
    sys.stdout.flush()
for name, pol, facs in [("Q","x","lfuncreate(1)"), ("Q(i)","x^2+1","lfuncreate(1),lfuncreate(-4)"), ("Q(sqrt-3)","x^2+x+1","lfuncreate(1),lfuncreate(-3)"),
                        ("Q(x^3-x-1)","x^3-x-1","lfuncreate(1),lfundiv(lfuncreate(x^3-x-1),lfuncreate(1))"), ("Q(sqrt5)","x^2-x-1","lfuncreate(1),lfuncreate(5)")]:
    run(name, pol, facs)
# sanity: wrong constant at the real place (1 instead of 1/2) must break the identity for Q
gp("nf=nfinit(x); Z=lfunzeros(lfuncreate(x),%d);" % H)
gp("fin=0; forprime(p=2,exp(B)+1, m=1; while(m*log(p)<B, fin+=log(p)*p^(-m/2)*f(m*log(p)); m++));")
print("Q with real-place factor 1 instead of 1/2: diff =", float(gp("2*Re_t+fin-(Fm(1/2)+Fm(-1/2)-sum(k=1,#Z,2*real(Fm(I*Z[k]))))")))

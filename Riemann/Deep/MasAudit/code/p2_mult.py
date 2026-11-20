# P6-P10: multiplicity / central-point / polar-multiplicity probes (PARI via cypari2).
import cypari2, sys
gp = cypari2.Pari(); gp.allocatemem(250*10**6); gp.set_real_precision(30)
def say(*a): print(*a); sys.stdout.flush()
# P6: Artin-forced multiplicities. S3 field K6 = splitting field of x^3-x-1: zeta_K6 = zeta*L(chi_-23)*L(rho)^2.
gp("K6=polredabs(nfsplitting(x^3-x-1)); Lr=lfundiv(lfuncreate(x^3-x-1),lfuncreate(1)); zr=lfunzeros(Lr,12);")
say("P6 S3: K6 =", gp("K6"), " first zeros of L(rho)=zeta_F3/zeta:", gp("zr"))
gp("L6=lfuncreate(K6); r1=1/2+I*zr[1];")
say("P6 S3: |zeta_K6(r1)|, |zeta_K6'(r1)|, |zeta_K6''(r1)| =", gp("[abs(lfun(L6,r1)), abs(lfun(L6,r1,1)), abs(lfun(L6,r1,2))]"), "(double zero expected)")
# A4: quartic x^4+8x+12 (Galois group A4): zeta_F4 = zeta*L(rho3); in the degree-12 closure zeta_K12 = zeta_k3 * L(rho3)^3.
say("P6 A4: galois group of x^4+8x+12:", gp("polgalois(x^4+8*x+12)"))
gp("L3=lfundiv(lfuncreate(x^4+8*x+12),lfuncreate(1)); z3=lfunzeros(L3,10);")
say("P6 A4: zeros of L(rho3) with 0<t<10 (each a zero of multiplicity >=3 of zeta_K12, [K12:Q]=12):", gp("z3"))
# P7: F=Q consequence of Thm 13.3.1.1 + Weil: every non-central zero of zeta and of Dirichlet L simple. Data check.
for D in [1,-4,-3,5,-7,8]:
    v = gp(f"L=lfuncreate({D}); Z=lfunzeros(L,100); [#Z, vecmin(vector(#Z,k,abs(lfun(L,1/2+I*Z[k],1))))]")
    say(f"P7 L(s,({D}/.)): zeros 0<t<100: {v[0]}, min |L'(rho)| = {float(v[1]):.3e}")
# P8: CM rank-3 central zero: E_1254: y^2 = x^3 - 1254^2 x (CM by Z[i]).
gp("E=ellinit([0,0,0,-1254^2,0]);")
say("P8 E_1254: root number", gp("ellrootno(E)"), " conductor", gp("ellglobalred(E)[1]"))
r = gp("R=ellrank(E); R")
say("P8 ellrank [lower, upper]:", gp("[R[1],R[2]]"), " points:", gp("R[4]"))
say("P8 regulator of the found points (det of height matrix):", gp("matdet(ellheightmatrix(E,R[4]))"))
say("P8 ellanalyticrank:", gp("ellanalyticrank(E)"))
# P9: central-point definition (p.528-529): unique point allowed as a pole of multiplicity >1 of partial Euler products.
# Q(i), chi=1, S={inf, (2+i), (2-i)}: prod_{v|5} L_v(s) = (1-5^-s)^-2.
gp("s1=2*Pi*I/log(5);"); gp("g(s)=(1-5^(-s))^2;")  # separate calls: a GP function body swallows the rest of the line
say("P9 Q(i), p=5: (1-5^-s)^2 and its derivative at s=2 pi i/log 5 (=", gp("s1"), "):", gp("[abs(g(s1)), abs(derivnum(t=s1,g(t))), abs(derivnum(t=s1,g(t),2))]"), "-> double pole of prod_{v|5}L_v at a NON-central point")
# P10: m_p for the CM Hecke character psi of Q(i) attached to y^2=x^3-x (conductor 32): split p=a^2+b^2 has a_p = 2 Re(pi);
# unitary values psi_v(pi_v) = pi/|pi|, conj; the two poles of L_v L_vbar on iR coincide iff (pi/|pi|)^2 = 1, i.e. b=0. Count.
gp("E32=ellinit([0,0,0,-1,0]);")
say("P10 CM psi (E32): #split p<10^5 with coincident poles (a_p^2 = 4p):", gp("c=0; n=0; forprime(p=3,10^5, if(p%4==1, n++; a=ellap(E32,p); if(a^2==4*p, c++))); [n,c]"), "-> m_p(psi)=1 for all p<10^5, m_F(psi)=1 < [F:Q]=2")
say("P10 chi=1 over Q(i): m_5 = 2 = [F:Q] (P9). Galois-type characters reach [F:Q] at completely split primes.")

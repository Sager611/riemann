
# Part D: remaining scalar certificates of HB 2610.12803 (exact rationals; arb only for reference logs)
from fractions import Fraction as F
from math import isqrt, comb, factorial
from flint import arb, ctx, fmpq
ctx.prec = 200
def E(x, n):
    s = F(0); term = F(1)
    for j in range(n+1):
        if j > 0: term = term * x / j
        s += term
    return s
def Ebar(x, n):
    assert x < n + 2
    return E(x, n) + x**(n+1) / factorial(n+1) / (1 - x/(n+2))
X = F(600000839515, 10)
out = []
def chk(name, cond, val=None): out.append((name, bool(cond), val))
chk("Ebar160(111/5) < 7*6e10/88", Ebar(F(111,5),160) < F(7*60000000000,88))
k = F(4675,988)
chk("k=21/4-8/(16-9/16)", k == F(21,4) - 8/(16-F(9,16)))
chk("E20(k*4/25) > 73691/34608", E(k*F(4,25),20) > F(73691,34608))
z = 1 + sum(F(10**12, n*n*isqrt(n*10**24)) for n in range(2,65)) + F(2, 3*64*8)
chk("zeta(5/2) cert < 269/200", z < F(269,200), float(z))
def A12(r):
    return 2*sum(r**(2*j+1)/(2*j+1) for j in range(13)) + 2*r**27/(27*(1-r*r))
def ell(n):
    m = n.bit_length()-1; r = F(n - 2**m, n + 2**m)
    return m*A12(F(1,3)) + A12(r)
J = sum(F(10**12,1)*ell(n)/(n*n*isqrt(n*10**24)) for n in range(2,65)) + F(1,512)*(F(2,3)*ell(64)+F(4,9))
chk("J cert < 39/100", J < F(39,100), float(J))
chk("log n <= ell_n (n=2..64, arb)", all(arb(n).log() < arb(fmpq(ell(n).numerator, ell(n).denominator)) for n in range(2,65)))
chk("E60(12) > 60000", E(F(12),60) > 60000)
Sval = F(111,20) - F(1,10**6) - ((F(1,2)+F(1,10**10))*F(39,100) + F(4,10**6))/(2 - F(269,200) - F(2,10**6))
chk("S(x,4,t) bound > 5.25228 > 21/4", Sval > F(525228,100000), float(Sval))
chk("F'/F bound < 0.297717", ((F(1,2)+F(1,10**10))*F(39,100) + F(4,10**6))/(2 - F(269,200) - F(2,10**6)) < F(297717,10**6))
# old bound
tstar = (1/arb(fmpq(4675,988))) * (arb(73691)/34608).log()
chk("t* < 4/25 and Lambda<=0.17972754...", tstar < arb(4)/25, (tstar + arb(1)/50).str(12))
# step-table normalizer
chk("Ebar160(557/25) < 7X/88", Ebar(F(557,25),160) < 7*X/88)
chk("Ebar160(13929/625) < 113X/1420", Ebar(F(13929,625),160) < 113*X/1420)
chk("69098^2 < 113X/1420", 69098**2 < 113*X/1420)
Larb = (arb(fmpq(X.numerator, X.denominator))/(4*arb.pi())).log()
chk("direct: L/4 at X (arb)", Larb/4 > arb(fmpq(55716,10000)), (Larb/4).str(12))
chk("direct: sqrt(X/4pi) >= 69098", (arb(fmpq(X.numerator,X.denominator))/(4*arb.pi())).sqrt() > 69098, ((arb(fmpq(X.numerator,X.denominator))/(4*arb.pi())).sqrt()).str(12))
J2 = sum(arb(n).log()/(n*n) for n in range(2,17)) + (arb(16).log()+1)/16
chk("J(2) bound < 1", J2 < 1, J2.str(10))
chk("log 65536 < 12", arb(65536).log() < 12)
v = arb(81)/40; N1 = arb(65536)
tail = N1**(1-v)*(N1.log()/(v-1) + 1/(v-1)**2)
chk("int_{N1}^inf log y y^{-81/40} < 1e-3", tail < arb(10)**-3, tail.str(10))
# Gaussian checks at A = 3e10-3, t=1/5 (rational exps bounded by Ebar100)
A = F(30000000000-3); t = F(1,5); D = 1 - t/A
g1 = (1/D)*Ebar(17/A + 16*t/(A*A*D), 100) - 1
g2 = F(2,10)*(1/D)*Ebar(17/A + t*(F(11,10)+4/A)**2/D, 100)
g3 = F(29,1000)*(1/D)*Ebar(17/A + t*(F(7,20)+4/A)**2/D, 100)
chk("gauss1 < 1e-8", g1 < F(1,10**8), float(g1)); chk("gauss2 < 3/10", g2 < F(3,10), float(g2)); chk("gauss3 < 4/100", g3 < F(4,100), float(g3))
chk("log3<11/10, log2/2<7/20", arb(3).log() < arb(11)/10 and arb(2).log()/2 < arb(7)/20)
# eps_sigma < 1/25 at a>=69098 (a-1 >= 69097)
chk("eps_sigma worst < 1/25", (F(3,10)*9**4 + F(4,100)*8)/69097 + F(1,10**8) < F(1,25))
P0 = arb(2).sqrt()/(4*arb.pi()**2); k0 = 1/((3-2*arb(2).log())*arb.pi()).sqrt()
chk("P0<0.036, kappa0<0.445", P0 < arb(36)/1000 and k0 < arb(445)/1000)
a = arb(60000)
chk("first-100 coeff sum < 0.029", P0*(arb.pi().sqrt()*k0 + k0**2/a)/(1-50*k0**2/a**2) < arb(29)/1000)
# d(x) < 1e-10 at x=6e10
x0 = arb(6*10**10); L0 = (x0/(4*arb.pi())).log()
chk("d(x)<1e-10", (L0**2/400 + arb(626)/1000)/(x0 - arb(666)/100) < arb(10)**-10, ((L0**2/400 + arb(626)/1000)/(x0 - arb(666)/100)).str(5))
# Remark thresholds: F4(sigma)=1/5 bracket and T4 = 1/5 bracket
def T4(sig):
    a_ = (2*arb(fmpq(sig.numerator,sig.denominator))-1)**2; b_=arb(1)/25
    k4 = arb(21)/4 - 8/(16-a_)
    return (1/k4)*((2+k4*a_)/(2+k4*b_)).log()
chk("F4 bracket (0.90664,0.90665)", T4(F(90664,100000))+arb(1)/50 < arb(1)/5 and T4(F(90665,100000))+arb(1)/50 > arb(1)/5)
chk("T4 bracket (0.93835,0.93836)", T4(F(93835,100000)) < arb(1)/5 and T4(F(93836,100000)) > arb(1)/5)
def T5(sig):
    a_ = (2*arb(fmpq(sig.numerator,sig.denominator))-1)**2
    k5 = 4 - 8/(25-a_)
    return (1/k5)*((2+k5*a_)/(2+k5/25)).log()
chk("F5 bracket (0.88328,0.88330)", T5(F(88328,100000))+arb(1)/50 < arb(1)/5 and T5(F(88330,100000))+arb(1)/50 > arb(1)/5)
chk("T5 bracket (0.91043,0.91045)", T5(F(91043,100000)) < arb(1)/5 and T5(F(91045,100000)) > arb(1)/5)
# height-5 bound S(x,5,t) > 111/20 - 1e-6 - 7/13 > 5
chk("S(x,5)>5", F(111,20) - F(1,10**6) - F(7,13) > 5)
# 0.19387... height 5 value: k5(9/16)
a_ = arb(9)/16; k5 = 4 - 8/(25-a_)
chk("height5 Lambda value", True, ((1/k5)*((2+k5*a_)/(2+k5/25)).log() + arb(1)/50).str(12))

# triangle and Euler step tables
pp = F(301,100); BB = pp*pp
tri = [(0,562500,42786,F(1,100),516903),(1,516903,44307,F(1,100),472419),(2,472419,45557,F(1,100),429333),(3,429333,46600,F(1,100),387828),(4,387828,47483,F(1,100),348016),(5,348016,48237,F(1,100),309958),(6,309958,48888,F(1,100),273677),(7,273677,49455,F(1,100),239169),(8,239169,49952,F(1,100),206410),(9,206410,50391,F(1,100),175361),(10,175361,50780,F(1,100),145974),(11,145974,51127,F(1,100),118194),(12,118194,51439,F(1,100),91960),(13,91960,51720,F(1,100),67208),(14,67208,51974,F(1,100),43874),(15,43874,52204,F(43,25000),40000)]
eul = [(0,562500,52708,F(1,100),510004),(1,510004,52851,F(1,100),460480),(2,460480,53006,F(1,100),413770),(3,413770,53166,F(1,100),369731),(4,369731,53323,F(1,100),328228),(5,328228,53475,F(1,100),289133),(6,289133,53619,F(1,100),252324),(7,252324,53754,F(1,100),217683),(8,217683,53881,F(1,100),185096),(9,185096,53999,F(1,100),154454),(10,154454,54108,F(1,100),125652),(11,125652,54209,F(1,100),98588),(12,98588,54303,F(1,100),73166),(13,73166,54390,F(1,100),49294),(14,49294,54470,F(51,12500),40000)]
triU = [l.split() for l in open('tri.txt').read().strip().split('\n')]
eulU = [l.split() for l in open('eul.txt').read().strip().split('\n')]
nu = F(1,2)+F(1,10**10); e0 = F(2,10**5); e1 = F(66,10**5)
def gate_rows(rows, Cstar):
    mC = None; mG = None; ok = True
    for idx, r in enumerate(rows):
        i, qi, Ci, di, qn = r; qi = F(qi,10**6); qn = F(qn,10**6); Ci = F(Ci,10**4)
        kk = 4*Ci/pp - 8/(BB-qi)
        ok = ok and kk > 0 and BB >= 5*qi and qn >= F(1,25) and qn < qi
        g = E(kk*di,20)*(2+kk*qn) - (2+kk*qi)
        c = Cstar[idx] - Ci
        mC = c if mC is None or c < mC else mC; mG = g if mG is None or g < mG else mG
    return ok, mC, mG
Ctri = []
for r in triU:
    UZ = F(int(r[1]),10**10); UJ = F(int(r[2]),10**10)
    Ctri.append(F(557,100) - F(1,10**6) - (nu*UJ + e1)/(2 - UZ - e0))
ok, mC, mG = gate_rows(tri, Ctri)
chk("triangle: chain/k>0/p^2>=5q", ok); chk("triangle min(C*-C) > 8.5e-7", mC > F(85,10**8), float(mC)); chk("triangle min gate > 3.48e-7", mG > F(348,10**9), float(mG))
tsum = sum(r[3] for r in tri); chk("triangle final time 3793/25000", tsum == F(3793,25000), float(tsum + F(1,50)))
CE = []
for r in eulU:
    UZ = F(int(r[1]),10**10); Up = F(int(r[2]),10**10); U0 = F(int(r[3]),10**10); U1 = F(int(r[4]),10**10)
    assert 1 - (U0+e0)*UZ > 0
    CE.append(F(557,100) - F(1,10**6) - (nu*(Up + U1*UZ) + e1*UZ)/(1 - (U0+e0)*UZ))
ok, mC, mG = gate_rows(eul, CE)
chk("Euler: chain/k>0/p^2>=5q", ok); chk("Euler min(C*-C) > 1.48e-6", mC > F(148,10**8), float(mC)); chk("Euler min gate > 1.25e-7", mG > F(125,10**9), float(mG))
esum = sum(r[3] for r in eul); chk("Euler final time 1801/12500", esum == F(1801,12500), float(esum + F(1,50)))
# Euler tails used in Prop euler-step
v2 = arb(401)/200
chk("tail n>60000 n^-2.005 < 1.7e-5", arb(60000)**(1-v2)/(v2-1) < arb(17)/10**6)
chk("tail n>60000 log n n^-2.005 < 2.1e-4", arb(60000)**(1-v2)*(arb(60000).log()/(v2-1)+1/(v2-1)**2) < arb(21)/10**5)
chk("tail n>65536 n^-2.005 < 1.6e-5", N1**(1-v2)/(v2-1) < arb(16)/10**6)
chk("tail n>65536 log < 1.9e-4", N1**(1-v2)*(N1.log()/(v2-1)+1/(v2-1)**2) < arb(19)/10**5)
for o in out: print(("OK  " if o[1] else "FAIL"), o[0], "" if o[2] is None else o[2])
print("ALL OK:", all(o[1] for o in out))


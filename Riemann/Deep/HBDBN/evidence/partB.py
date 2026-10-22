
# Part B/C: Arb verification that printed table entries are rigorous upper bounds (HB 2610.12803)
import sys, time
from flint import arb, ctx, fmpq
ctx.prec = 128
t0 = time.time()
N1 = 65536
def A(x): return arb(fmpq(x.numerator, x.denominator)) if hasattr(x,'numerator') else arb(x)
from fractions import Fraction as F
logs = [None, arb(0)] + [arb(n).log() for n in range(2, N1+1)]
# primes
sieve = bytearray([1])*(N1+1); sieve[0]=sieve[1]=0
for i in range(2, int(N1**0.5)+1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
primes = [i for i in range(N1+1) if sieve[i]]
ppow = []  # (n, j, log l)
for l in primes:
    n = l; j = 1
    while n <= N1:
        ppow.append((n, j, l)); n *= l; j += 1
logN1 = logs[N1]
def T(j, v):
    s = arb(0); fact = 1
    import math
    for k in range(j+1):
        s += arb(math.factorial(j)//math.factorial(j-k)) * logN1**(j-k) / (v-1)**(k+1)
    return arb(N1)**(1-v) * s
def sums(u, tp, need):
    # returns dict of finite sums n<=N1
    Z = arb(1); J = arb(0); d0 = arb(0); d1 = arb(0); d2 = arb(0)
    for n in range(2, N1+1):
        ln = logs[n]
        npw = (-u*ln).exp()
        Z += npw
        if 'J' in need: J += ln*npw
        if tp is not None:
            hm1 = (tp*ln*ln/4).expm1()
            x = hm1*npw
            d0 += x
            x = x*ln; d1 += x
            if 'd2' in need: d2 += x*ln
    psi = arb(0); psi2 = arb(0)
    if 'psi' in need:
        for (n, j, l) in ppow:
            ll = logs[l]; npw = (-u*logs[n]).exp()
            psi += ll*npw; psi2 += j*ll*ll*npw
    return dict(Z=Z, J=J, d0=d0, d1=d1, d2=d2, psi=psi, psi2=psi2)
def ub_ok(U_int, den, val):
    # True iff U_int/den > val rigorously
    return (arb(U_int) - val*arb(den)) > 0
def margin(U_int, den, val):
    return (arb(U_int)/arb(den) - val)

bad = []
# ---------- kernel table (main result) ----------
zeta_tab = [l.split() for l in open('kz.txt').read().strip().split('\n')]
pert_tab = [l.split() for l in open('kp.txt').read().strip().split('\n')]
p = arb(61)/20; sigp = arb(81)/40
Nst = arb(69098); Lst = 2*Nst.log() + 3
minrel = {}
for i in range(24):
    tm = F(i,200); d = F(1,200) if i < 23 else F(463,100000); tp = tm + d
    u = sigp + A(F(557,100)*tm); u8 = sigp + A(F(278,100)*tm)
    e = (p-1)/2 - A(F(278,100)*tm)
    s = sums(u, A(tp), {'d2','psi'})
    Zs = s['Z'] + T(0,u); psis = s['psi'] + T(1,u); psi2s = s['psi2'] + T(2,u)
    d0 = s['d0'] + T(0,u8); d1 = s['d1'] + T(1,u8); d2 = s['d2'] + T(2,u8)
    rho0 = arb(101)/100 * Nst**(-p) * (Nst**(e+1)/(e+1) + Nst**e)
    rho1 = Lst*rho0; rho2 = (Lst*Lst + arb(10)**-9)*rho0
    vals = dict(Z=Zs, psi=psis, psi2=psi2s, d0=d0, d1=d1, d2=d2, r0=rho0, r1=rho1, r2=rho2)
    U = dict(Z=int(zeta_tab[i][1]), psi=int(zeta_tab[i][2]), psi2=int(zeta_tab[i][3]),
             d0=int(pert_tab[i][1]), d1=int(pert_tab[i][2]), d2=int(pert_tab[i][3]),
             r0=int(pert_tab[i][4]), r1=int(pert_tab[i][5]), r2=int(pert_tab[i][6]))
    for k in vals:
        ok = ub_ok(U[k], 10**12, vals[k])
        m = margin(U[k], 10**12, vals[k])
        rel = float((m/vals[k]).mid())
        if k not in minrel or rel < minrel[k][0]: minrel[k] = (rel, i, float(m.mid()))
        if not ok: bad.append(('kernel', i, k, U[k], vals[k].str(20)))
    if i in (0, 5, 23):
        print("kernel row", i, {k: vals[k].str(15) for k in vals}, flush=True)
print("kernel table: min relative margin (U-val)/val per column:", minrel, flush=True)
print("kernel bad entries:", bad, "  elapsed %.1fs"%(time.time()-t0), flush=True)

# ---------- triangle table (p=301/100) ----------
tri = [l.split() for l in open('tri.txt').read().strip().split('\n')]
sigp2 = arb(401)/200
minrel = {}
for i in range(16):
    u = sigp2 + A(F(278,100)*F(i,100))
    s = sums(u, None, {'J'})
    Zs = s['Z'] + T(0,u); Js = s['J'] + T(1,u)
    for k, v, Ui in (('Z', Zs, int(tri[i][1])), ('J', Js, int(tri[i][2]))):
        if not ub_ok(Ui, 10**10, v): bad.append(('tri', i, k, Ui, v.str(20)))
        rel = float((margin(Ui,10**10,v)/v).mid())
        if k not in minrel or rel < minrel[k][0]: minrel[k] = (rel, i)
print("triangle table min rel margins:", minrel, " elapsed %.1fs"%(time.time()-t0), flush=True)

# ---------- Euler table ----------
eul = [l.split() for l in open('eul.txt').read().strip().split('\n')]
teu = [F(i,100) for i in range(15)] + [F(7,50)+F(51,12500)]
minrel = {}
for i in range(15):
    u = sigp2 + A(F(557,100)*teu[i]); tp = A(teu[i+1])
    s = sums(u, tp, {'psi'})
    Zs = s['Z'] + T(0,u); psis = s['psi'] + T(1,u)
    d0 = s['d0'] + arb(4)/10**5; d1 = s['d1'] + arb(5)/10**4
    for k, v, Ui in (('Z',Zs,int(eul[i][1])),('psi',psis,int(eul[i][2])),('d0',d0,int(eul[i][3])),('d1',d1,int(eul[i][4]))):
        if not ub_ok(Ui, 10**10, v): bad.append(('eul', i, k, Ui, v.str(20)))
        rel = float((margin(Ui,10**10,v)/v).mid())
        if k not in minrel or rel < minrel[k][0]: minrel[k] = (rel, i)
print("Euler table min rel margins:", minrel, " elapsed %.1fs"%(time.time()-t0), flush=True)
print("ALL BAD:", bad)


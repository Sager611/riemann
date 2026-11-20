# elem_arb.py -- own Arb certificate (stc2) of sup_{T>=1} max(F_u, F_l) (mode N) and of max(F_u^S, F_l^S) (mode S),
# Casana-Giner (Zenodo 23866074) (4.2), (4.3), (4.8), (4.16), Prop. 4.9, Lemmas 4.7, 4.8.  Usage: elem_arb.py NAME A1 N|S [T2]
# [1, T2] (T2 = 1e8): adaptive bisection of a geometric partition (ratio 1.01); on [a,b] the derivative F'([a,b]) is enclosed
# by ball arithmetic (natural extension of the exact derivative, polygammas on the complex ball); if F' < 0 (> 0) on [a,b] the
# bound is F(a) (F(b)), else the centred form F(m) + |F'| (b-a)/2; intervals are bisected until the bound is within 1e-9
# of the best certified point value. T >= T2: Lemma 4.7 (the paper's K of (4.15), evaluated at T2).
import json, sys, time
from fractions import Fraction as Fr
from math import factorial
from flint import arb, acb, ctx
ctx.prec = 96
KJ = json.load(open('Riemann/Deep/STC2/kernels.json'))
name = sys.argv[1]; A1 = arb(Fr(sys.argv[2]).numerator)/Fr(sys.argv[2]).denominator; MODE = sys.argv[3]
T2 = float(sys.argv[4]) if len(sys.argv) > 4 else 1e8
k = KJ[name]; dq = Fr(k['d']).limit_denominator(10**6); d = arb(dq.numerator)/dq.denominator
B = [(X, int(j)) for X, j in k['B']]; C = [arb(c) for c in k['c']]; J = max(j for _, j in B)
PI = arb.pi(); LOGPI = PI.log(); s1r = arb(1)/2 + d; s2r = arb(1)/2 + 2*d; cR1 = C[0]; C1 = cR1/(2*PI)
assert B[0] == ('R', 1) and A1 >= C1
def Dj(j, s, ps):
    if j == 1: return 1/s + 1/(s-1) - LOGPI/2 + ps[0]/2
    return s**(-j) + (s-1)**(-j) + arb((-1)**(j-1))/factorial(j-1)/2**j*ps[j-1]
def dDj(j, s, ps):
    if j == 1: return -s**-2 - (s-1)**-2 + ps[1]/4
    return -j*s**(-j-1) - j*(s-1)**(-j-1) + arb((-1)**(j-1))/factorial(j-1)/2**(j+1)*ps[j]
def combine(vals):
    R = arb(0); I = arb(0)
    for v, (X, j), c in zip(vals, B, C):
        if X == 'R': R += c*v.real
        else: I -= c*v.imag
    return R, I
def F(T):   # (F_u, F_l) at T (arb; a point)
    s1 = acb(s1r, T); s2 = acb(s2r, T)
    Ael = lambda s, sr: (T/sr).atan() + (T/(sr-1)).atan() - T/2*LOGPI + (s/2).lgamma().imag
    base = (2*Ael(s1, s1r) - Ael(s2, s2r))/PI
    if MODE == 'N': base -= T/(2*PI)*(T/(2*PI*arb(1).exp())).log()
    else: base -= (acb(arb(1)/4, T/2).lgamma().imag - T/2*LOGPI)/PI + 1
    ps = [(s1/2).polygamma(kk) for kk in range(J)]
    R, I = combine([Dj(j, s1, ps) for (_, j) in B])
    return base + (R + I)/PI - A1*T.log(), -(base - (R - I)/PI) - A1*T.log()
def dF(T):  # (F_u', F_l') on a ball T
    s1 = acb(s1r, T); s2 = acb(s2r, T)
    dA = lambda s: (1/s).real + (1/(s-1)).real - LOGPI/2 + (s/2).digamma().real/2
    db = (2*dA(s1) - dA(s2))/PI
    if MODE == 'N': db -= (T/(2*PI)).log()/(2*PI)
    else: db -= (acb(arb(1)/4, T/2).digamma().real/2 - LOGPI/2)/PI
    ps = [(s1/2).polygamma(kk) for kk in range(J+1)]
    R, I = combine([acb(0, 1)*dDj(j, s1, ps) for (_, j) in B])
    return db + (R + I)/PI - A1/T, -(db - (R - I)/PI) - A1/T
t0 = time.time(); cache = {}
def Fp(x):
    if x not in cache: cache[x] = F(arb(x))
    return cache[x]
edges = [1.0]
while edges[-1] < T2: edges.append(min(edges[-1]*1.01, T2))
best = [None, None]; worst = [None, None]; where = [None, None]; n = 0
for x in edges:
    v = Fp(x)
    for i in (0, 1):
        lb = float(v[i].lower())
        if best[i] is None or lb > best[i]: best[i] = lb; where[i] = x
stack = [(edges[i], edges[i+1], 0) for i in range(len(edges)-1)]
while stack:
    a, b, dep = stack.pop(); n += 1
    Tb = arb((a + b)/2, (b - a)/2*(1 + 1e-15)); g = dF(Tb); m = (a + b)/2
    ubs = []
    for i in (0, 1):
        if g[i] < 0: ub = float(Fp(a)[i].upper())
        elif g[i] > 0: ub = float(Fp(b)[i].upper())
        else:
            v = Fp(m)[i]; ub = float((v + abs(g[i])*arb((b - a)/2)).upper())
            if float(v.lower()) > best[i]: best[i] = float(v.lower()); where[i] = m
        ubs.append(ub)
    if all(ubs[i] <= best[i] + 1e-9 or ubs[i] <= max(best) - 1e-6 for i in (0, 1)) or dep >= 40:
        for i in (0, 1):
            if worst[i] is None or ubs[i] > worst[i]: worst[i] = ubs[i]
        continue
    stack += [(a, m, dep+1), (m, b, dep+1)]
# T >= T2: Lemma 4.7 with the paper's K (4.15)
T = arb(T2); s1 = s1r; s2 = s2r
K = (2*s2-1)/PI + (s1*s1/2 + (s2-1)*s2/2 + 1)/PI + (2*s2*s2/4 + (s1-1)*s1)/PI
K += abs(cR1)/PI*((2*s1-1) + s1*s1/4 + arb(8)/3)/T
for (X, j), cc in zip(B, C):
    if j == 1: continue
    kk = j - 1
    ps = arb(factorial(kk-1))*(2/T)**kk + arb(factorial(kk))*(2/T)**(kk+1)/2 + arb(factorial(kk+1))*(4/(3*T))*(4/T)**(kk+1)
    K += abs(cc)/PI*(2/T**j + ps/(factorial(j-1)*2**j))*T
if MODE == 'N': Fi = (arb(7)/8 - C1*(2*PI).log(), -arb(7)/8 - C1*(2*PI).log())
else: Fi = (-C1*(2*PI).log(), -C1*(2*PI).log()); K = K + (arb(1)/16 + arb(1)/8 + arb(1)/3)/PI
tail = [float((Fi[i] + K/T).upper()) for i in (0, 1)]
sup = max(worst[0], worst[1], tail[0], tail[1])
print('%s mode %s A1=%s: [1,%.0e] %d intervals; F_u <= %.10f (point max %.10f at T=%.6g), F_l <= %.10f (point max %.10f at T=%.6g); T>=%.0e (Lemma 4.7, K=%.4f): F_u <= %.6f, F_l <= %.6f; SUP <= %.10f  %.0f s' % (
    name, MODE, sys.argv[2], T2, n, worst[0], best[0], where[0], worst[1], best[1], where[1], T2, float(K.mid()), tail[0], tail[1], sup, time.time()-t0))
json.dump(dict(name=name, mode=MODE, A1=sys.argv[2], sup_upper=sup, Fu_upper=worst[0], Fl_upper=worst[1], Fu_point=best[0], Fl_point=best[1],
               T_Fu=where[0], T_Fl=where[1], tail=tail, K=float(K.mid()), intervals=n), open('Riemann/Deep/STC2/elem_%s_%s.json' % (name, MODE), 'w'))


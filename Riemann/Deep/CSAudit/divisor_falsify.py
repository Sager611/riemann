# Falsification tests for the explicit divisor-sum lemmas of Chourasiya-Simonic §5.2 (TeX l. 955-1243):
# Lemma 5.2 (lem:d), 5.3 (lem:dExp), 5.4 (lem:Sumd2, Cully-Hugill-Trudgian), 5.5 (lem:d2), 5.6 (lem:d2Exp1), 5.7 (lem:d2Exp2), 5.8 (lem:d2Exp3).
# For each lemma: max over a grid of (x0, sigma) and all integers x in [x0, N] (worst x for step sums; both one-sided limits for tails)
# of LHS/RHS. A value > 1 is a counterexample. numpy float64; tails for sigma > 1 via exact zeta(s)^2, zeta(s)^4/zeta(2s) (mpmath).
import numpy as np, mpmath as mp, sys, math
mp.mp.dps = 30
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
g = float(mp.euler); pi2 = float(mp.pi**2)
d = np.zeros(N + 1, dtype=np.int32)
for k in range(1, N + 1): d[k::k] += 1
n = np.arange(N + 1, dtype=np.float64); n[0] = 1.0
d = d.astype(np.float64); d2 = d*d; d[0] = 0; d2[0] = 0
L = np.log(n)
def D1(x0, s0):
    l = np.log(x0); r = np.sqrt(x0)
    return 1 + (2*g - 1 - s0/(1 - s0))/l + (1 - s0)*(0.5 + abs(s0) - s0)/r + (1 - 2*g*s0)/(r*l)
def D2(x0, s0):
    l = np.log(x0); return s0 - 1 + (1 + 2*g*(s0 - 1))/l + (4*s0 - 1)/((2*s0 - 1)*np.sqrt(x0)*l)
def D3(x0):
    l = np.log(x0); e = np.e; return 1/e + (1 + 2*g)/(e*l) + 5/(2*e*np.sqrt(x0)*l)
def D41(x0, s):
    l = np.log(x0)
    return (1/((1 - s)*pi2) + (-3*s/((1 - s)**2*pi2) + 0.744342/(1 - s))/l
            + (6*s/((1 - s)**3*pi2) - 1.488684*s/(1 - s)**2 + 0.823266/(1 - s))/l**2
            + (-6*s/((1 - s)**4*pi2) + 1.488684*s/(1 - s)**3 - 0.823266*s/(1 - s)**2 + 0.460324/(1 - s))/l**3 + 17.52/(x0**0.25*l**2))
def D42(x0, s):
    l = np.log(x0)
    return (1/pi2 + 0.744342/l + 0.823266/l**2 + 0.460324/l**3 + 9.73*s*(1 - s)/(x0**min(0.25, 1 - s)*l) + 9.73*(1 - s)/(x0**0.25*l**2))
def D4(x0, s): return D41(x0, s) if s < 0 else D42(x0, s)/(1 - s)
def D5(x0, s):
    l = np.log(x0); a = s - 1
    return (a**3/pi2 + a**2*(3*s/pi2 + 0.744342*a)/l + a*(6*s/pi2 + 1.488684*s*a + 0.823266*a**2)/l**2
            + (6*s/pi2 + 1.488684*s*a + 0.823266*s*a**2 + 0.460324*a**3)/l**3 + (1 + s/(s - 0.75)*(1 + 1/((s - 0.75)*l)))*9.73*a**4/(x0**0.25*l**2))
def D6(x0):
    l = np.log(x0); return 0.349436 + 1.155975/l + 1.28359/l**2 + 154/l**3 + 9.73/(x0**0.25*l**2)
def D7(x0):
    l = np.log(x0); return 3/(2*np.e*np.pi)**2 + 0.12441/l + 0.364/l**2 + 0.3/l**3 + 4.61/(x0**0.25*l**2)
X0S = [2, 2.5, 3, 4, 5, 7, 10, 20, 50, 100, 1000, 1e4, 1e5]
def worst(ratio, xs):
    if len(ratio) == 0: return -1.0, None
    i = int(np.argmax(ratio)); return float(ratio[i]), float(xs[i])
res = {}
def rec(name, val, info):
    if name not in res or val > res[name][0]: res[name] = (val, info)
xi = np.arange(2, N + 1)                   # integer x (worst case for sums over n <= x)
xf = xi.astype(np.float64)
for s in [-3, -2.5, -2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.1, 0.25, 0.4, 0.5]:
    ps = np.cumsum(d*n**(-s))              # sum_{n<=x} d(n)/n^s
    for s0 in [s0 for s0 in [-3, -2, -1, -0.5, -0.25, 0] if s0 <= min(s, 0)]:
        for x0 in X0S:
            m = xi >= x0
            r = ps[xi[m]]/(D1(x0, s0)/(1 - s)*xf[m]**(1 - s)*np.log(xf[m]))
            v, xw = worst(r, xf[m]); rec('L5.2a sum_{n<=x} d/n^s', v, (x0, s0, s, xw))
    if -3 <= s < 1:
        p2 = np.cumsum(d2*n**(-s))
        for x0 in X0S:
            m = xi >= x0
            r = p2[xi[m]]/(D4(x0, s)*xf[m]**(1 - s)*np.log(xf[m])**3)
            v, xw = worst(r, xf[m]); rec('L5.5a sum_{n<=x} d^2/n^s', v, (x0, s, xw))
for s in [1.001, 1.01, 1.05, 1.1, 1.25, 1.5, 2, 3, 4]:
    a1t = d*n**(-s); a2t = d2*n**(-s)
    tN1 = float(mp.zeta(s)**2 - mp.mpf(math.fsum(a1t[1:]))); tN2 = float(mp.zeta(s)**4/mp.zeta(2*s) - mp.mpf(math.fsum(a2t[1:])))
    R1 = np.append(np.cumsum(a1t[::-1])[::-1], 0.0); R2 = np.append(np.cumsum(a2t[::-1])[::-1], 0.0)   # R[m] = sum_{m<=n<=N}, R[N+1] = 0
    for x0 in X0S:
        idx0 = xi[xi >= x0]; xs = idx0.astype(np.float64)
        for side, idx in (('x->m-', idx0 - 1), ('x=m', idx0)):
            if side == 'x->m-':          # x -> m^- must stay >= x0: drop m <= x0
                keep = xs > x0; idx, xs_s = idx[keep], xs[keep]
            else: xs_s = xs  # tail sum_{n>x}: x -> m^- includes n = m
            t1 = R1[idx + 1] + tN1; t2 = R2[idx + 1] + tN2      # sum_{n > idx}
            for s0 in [s0 for s0 in [1.001, 1.01, 1.05, 1.1, 1.25, 1.5, 2, 3, 4] if s0 <= s]:
                r = t1/(D2(x0, s0)/(s0 - 1)**2*xs_s**(1 - s)*np.log(xs_s))
                v, xw = worst(r, xs_s); rec('L5.2b sum_{n>x} d/n^s', v, (x0, s0, s, side, xw))
            r = t2/(D5(x0, s)/(s - 1)**4*xs_s**(1 - s)*np.log(xs_s)**3)
            v, xw = worst(r, xs_s); rec('L5.5b sum_{n>x} d^2/n^s', v, (x0, s, side, xw))
# Lemma 5.4: |D0(x)| <= 9.73 x^{3/4} log x, x >= 2 (check at integers from both sides)
S2 = np.cumsum(d2); xs = xi.astype(np.float64); lx = np.log(xs)
main = xs*lx**3/pi2 + 0.7443415*xs*lx**2 + 0.8232655*xs*lx + 0.4603235*xs
for side, val in (('x=m', S2[xi]), ('x->m-', S2[xi - 1])):
    kk = xs > 2 if side == 'x->m-' else xs >= 2
    r = np.abs(val[kk] - main[kk])/(9.73*xs[kk]**0.75*lx[kk]); v, xw = worst(r, xs[kk]); rec('L5.4 |D0(x)|/(9.73x^{3/4}log x)', v, (side, xw))
# Lemma 5.6: |sum_{n<=x} d^2/n - log^4 x/(4 pi^2)| <= D6(x0) log^3 x
p1 = np.cumsum(d2/n)
for x0 in X0S:
    xs_ = xi[xi >= x0].astype(np.float64)
    for side, val in (('x=m', p1[xi[xi >= x0]]), ('x->m-', p1[xi[xi >= x0] - 1])):
        kk = (xs_ > x0) if side == 'x->m-' else (xs_ >= x0)
        r = np.abs(val[kk] - np.log(xs_[kk])**4/(4*pi2))/(D6(x0)*np.log(xs_[kk])**3); v, xw = worst(r, xs_[kk]); rec('L5.6 d^2/n sum', v, (x0, side, xw))
# Lemma 5.3: first: sum_{n<=x} d/n^s (1-e^{-n/x}) <= D1(x0, s-1)/(2-s) x^{1-s} log x, s in [0,1]; second: sum_{n>x} d n^{-s} e^{-n/x} <= D3 x^{1-s} log x
# Lemma 5.7: sum_{n>x} d^2 n^{-s} e^{-2n/x} <= D7(x0) x^{1-s} log^3 x, x >= x0 >= e^4, s >= -1
# Lemma 5.8: sum_{n<=x} d^2 n^{-s} (e^{-n/x}-1)^2 <= D4(x0, s-2) x^{1-s} log^3 x, s in [-1, 2]
xgrid = np.unique(np.concatenate([np.arange(2, 400), np.geomspace(400, N/45, 300)]))
for x in xgrid:
    k = np.arange(1, int(min(N, 45*x)) + 1).astype(np.float64); ki = k.astype(np.int64); w = np.exp(-k/x); lo = k <= x; hi = ~lo
    for s in [0, 0.25, 0.5, 0.75, 1]:
        a = (d[ki]*k**(-s)*(1 - w))[lo].sum(); b = (d[ki]*k**(-s)*w)[hi].sum()
        for x0 in [x0 for x0 in X0S if x0 <= x]:
            rec('L5.3a', a/(D1(x0, s - 1)/(2 - s)*x**(1 - s)*np.log(x)), (x0, s, float(x)))
            rec('L5.3b', b/(D3(x0)*x**(1 - s)*np.log(x)), (x0, s, float(x)))
    for s in [-1, -0.5, 0, 0.5, 1, 1.5, 2]:
        c = (d2[ki]*k**(-s)*(1 - w)**2)[lo].sum()
        for x0 in [x0 for x0 in X0S if x0 <= x]:
            rec('L5.8', c/(D4(x0, s - 2)*x**(1 - s)*np.log(x)**3), (x0, s, float(x)))
    if x >= np.exp(4):
        for s in [-1, -0.5, 0, 0.5, 1, 2, 3]:
            e2 = (d2[ki]*k**(-s)*w**2)[hi].sum()
            for x0 in [x0 for x0 in [np.exp(4), 60, 100, 1000, 1e4] if x0 <= x]:
                rec('L5.7', e2/(D7(x0)*x**(1 - s)*np.log(x)**3), (x0, s, float(x)))
print('N =', N)
for k_, (v, info) in sorted(res.items()):
    print('%-36s max LHS/RHS = %.6f  at %s%s' % (k_, v, info, '   <-- COUNTEREXAMPLE' if v > 1 else ''))

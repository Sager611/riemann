# Recompute E of Bi-Hou Thm 4.8 / Lemma 4.5 at their parameters, check hypotheses, window, margins. mpmath, 60 digits.
from mpmath import mp, mpf, log, exp, sqrt, tan, pi, ceil, floor, nstr
mp.dps = 60
def gamma6(x):  # upper incomplete Gamma(6, x), exact closed form
    return exp(-x) * (x**5 + 5*x**4 + 20*x**3 + 60*x**2 + 120*x + 120)
def logB(u, k, T0, T1, delta, M):
    c0 = tan(1); t6log = 3*u/10 - mpf(3)/5
    lB1 = log(4) - u/4 + log(4*exp(150) + T0**(1 + c0*log(4) + 11) * exp(14))
    lB2 = -k*T0**2 + c0*log(T0)*log(1/delta) + 11*log(T0) + 16
    lB3 = log(u) + M + 2 + log(gamma6(k*T0**2)) - log(2*k**6)
    lB4 = -k*T1**2 + 16*log(T1) + c0*log(T1)*log(u) + 18
    lB5 = log(2) + 18 - k*T1**2/2 - log(k*T1)
    # k*t6^2 = exp(2*t6log + log k) is astronomically large; keep log B6 <= -(k t6^2)/2 symbolically as -exp(L6) with L6 below
    L6 = 2*t6log + log(k) - log(2)   # valid since 1 + 3u/2 - log(2k) - t6log <= k t6^2 / 2 here
    assert 1 + 3*u/2 - log(2*k) - t6log < exp(min(L6, mpf(50)))
    lB6 = -exp(min(L6, mpf(10**6))) if L6 < 10**6 else -mpf(10)**30  # placeholder far below -1e6; true value -exp(L6)
    if L6 > 50 and not getattr(logB, 'said', False): print('  (log B6 <= -exp(%s); B6 reported as 0)' % nstr(L6, 8)); logB.said = True
    return [lB1, lB2, lB3, lB4, lB5, lB6]
def ex(x):  # exp with clamp: terms below e^-1e6 are reported as 0 (they are < 10^-434294)
    return exp(x) if x > -10**6 else mpf(0)
def report(tag, v, k, T0, T1, delta, M, K2low):
    u = v - k; c0 = tan(1); A = c0*log(u) + 16
    conds = {'u>16': u > 16, '0<k<1/4': 0 < k < mpf(1)/4, '4<=T0<T1': 4 <= T0 < T1, 'T1+c0<=3e12': T1 + c0 <= mpf('3e12'),
             '(5logT1+3)/u<1/2': (5*log(T1) + 3)/u < mpf(1)/2, 'T1>=sqrt(A/k)': T1 >= sqrt(A/k), '(k/2)T1^2>=A logT1': k*T1**2/2 >= A*log(T1),
             'uT1>64': u*T1 > 64, '0<delta<=1': 0 < delta <= 1}
    lB = logB(u, k, T0, T1, delta, M)
    tail_log = log(2) - u/2 + 2*sqrt(k*u)
    E = sum(ex(x) for x in lB)/pi + ex(tail_log)
    print('==', tag); print('  conditions', all(conds.values()), {a: bool(b) for a, b in conds.items() if not b})
    print('  A =', nstr(A, 8), ' sqrt(A/k) =', nstr(sqrt(A/k), 8), ' (k/2)T1^2 - A logT1 =', nstr(k*T1**2/2 - A*log(T1), 8), ' kT0^2 =', nstr(k*T0**2, 8))
    print('  log B1..B6 =', [nstr(x, 6) for x in lB]); print('  log tail =', nstr(tail_log, 6))
    print('  E =', nstr(E, 6), '  K2low - 1 - E =', nstr(K2low - 1 - E, 8))
    w = 2*sqrt(k*u); print('  u - 2sqrt(ku) =', nstr(u - w, 25), '  u + 2sqrt(ku) =', nstr(u + w, 25), '  ceil =', int(ceil(u + w)))
    # margins: largest M (resp. smallest kT0^2 shift) keeping E < K2low - 1
    room = log((K2low - 1) * pi) - lB[2]
    print('  slack in M (B3 alone, other terms negligible):', nstr(room, 6))
    return E
v = mpf('3484971332733132396.5383341489440461'); k = mpf('1.5e-10')
K2low = mpf('1.000305250377074143011122')
report('Bi-Hou Sec.7', v, k, mpf(1300000), mpf('4.3e6'), mpf('0.18286'), mpf('22.7671'), K2low)
report('Bi-Hou Rem.7.2 (StR M=23.05, T1=1e7)', v, k, mpf(1300000), mpf('1e7'), mpf('0.18286'), mpf('23.05'), K2low)
# Table 2 spot check: smallest integer T0 with E <= 1e-6 at k=1.5e-10 (T1=1e7, M=23.05, delta=0.1)
for kk, T0tab in (('1.5e-10', 1272608), ('3e-9', 273596), ('1e-10', 1566567)):
    kk = mpf(kk); u = v - kk
    def Ef(T0):
        lB = logB(u, kk, mpf(T0), mpf('1e7'), mpf('0.1'), mpf('23.05')); return sum(ex(x) for x in lB)/pi + ex(log(2) - u/2 + 2*sqrt(kk*u))
    print('Table2 k=%s T0=%d: E(T0)=%s E(T0-1)=%s' % (nstr(kk, 3), T0tab, nstr(Ef(T0tab), 6), nstr(Ef(T0tab - 1), 6)))
# Gallego's parameters, Bi-Hou's B-form (k=5.7e-12, T0=7999035.172, T1=22494898, M=24.25): what E does Bi-Hou's Lemma 4.5 give?
yg = mpf('3484971332733132396.538348536063730810969705971170')
report('Gallego params in Bi-Hou form (delta=0.1 nominal)', yg, mpf('5.7e-12'), mpf('7999035.172'), mpf('22494898'), mpf('0.1'), mpf('24.25'), mpf('1.0020694979178448'))

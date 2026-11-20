# Independent evaluation of Kadiri-Lumley-Ng Lemma 4.14 (arXiv 2101.12263) from the TeX definitions, mpmath 256 bits.
# Not using cert.py's formulas: I(A,n) by direct quadrature, every constant typed from the TeX (line numbers in comments).
# Part 1: H0 = 3.0610046e10 -> compare with Table 1 (l. 1692-1717). Part 2: H0 = 3e12 -> admissibility, A, B; compare with cert.py.
import sys, os, json
from mpmath import mp, mpf, sqrt, log, exp, pi, zeta, quad, inf, euler, findroot, binomial
mp.prec = 256
b1, b2, b3, b4 = mpf('0.62'), mpf('1.048'), mpf('0.605'), mpf('0.529')       # l. 429-434
a1, a2 = mpf('0.63'), mpf('2.851')                                           # l. 416-418
m0 = sqrt(1 + mpf(2)/3*sqrt(mpf(6)/5))                                       # l. 559-560
g = +euler; beta = 2
ROWS = [('0.60','0.5','1.251','0.288','0.3140','0.341','2.177','5.663'),('0.65','0.6','1.249','0.256','0.3070','0.340','2.963','5.249'),
('0.70','0.8','1.247','0.222','0.3040','0.339','3.983','4.824'),('0.75','1.0','1.245','0.189','0.3030','0.338','5.277','4.403'),
('0.80','1.0','1.245','0.160','0.3030','0.337','6.918','3.997'),('0.85','1.0','1.245','0.133','0.3030','0.336','8.975','3.588'),
('0.86','1.0','1.245','0.127','0.3030','0.335','9.441','3.514'),('0.87','1.0','1.245','0.122','0.3030','0.335','9.926','3.430'),
('0.88','1.0','1.245','0.116','0.3030','0.335','10.431','3.346'),('0.89','1.0','1.245','0.111','0.3030','0.335','10.955','3.262'),
('0.90','1.0','1.245','0.105','0.3030','0.334','11.499','3.186'),('0.91','1.0','1.245','0.100','0.3030','0.334','12.063','3.102'),
('0.92','1.0','1.245','0.095','0.3030','0.334','12.646','3.017'),('0.93','1.0','1.245','0.089','0.3030','0.333','13.250','2.941'),
('0.94','1.0','1.245','0.084','0.3030','0.333','13.872','2.856'),('0.95','1.0','1.245','0.079','0.3030','0.333','14.513','2.772'),
('0.96','1.0','1.245','0.074','0.3030','0.332','15.173','2.694'),('0.97','1.0','1.245','0.069','0.3030','0.332','15.850','2.609'),
('0.98','1.0','1.245','0.064','0.3030','0.331','16.544','2.532'),('0.99','1.0','1.245','0.060','0.3030','0.331','17.253','2.446')]
_Ic = {}
def I(A, n, al):            # l. 815: int_0^oo x^A e^{-2 al x^2} (log x)^n dx, by quadrature
    key = (A, n, al)
    if key not in _Ic:
        f = lambda x: x**A*exp(-2*al*x**2)*log(x)**n
        _Ic[key] = quad(f, [0, mpf(1)/4, 1, 2, 4, inf])
    return _Ic[key]
def run(H0, eta, rows, mode):
    L0 = log(H0); ll0 = log(L0); H = H0 - 1
    out = []
    for r in rows:
        s0, k, mu, al, de, d, At, Bt = [mpf(v) for v in r]
        X0 = k*H0; LX0 = log(X0)
        C1k = 6/pi**2 + b2/LX0                                              # l. 644
        C2k = pi*m0*b1/LX0 + 6*m0/(pi*X0) + pi*m0*b2/(X0*LX0)               # l. 647-649
        a3 = -6*a1/exp(1) + a2; C3 = a3**2*C1k*LX0                          # l. 683-686
        C4 = C1k*a1**2*(1 + 1/sqrt(C3))**2                                   # l. 689-691
        rr = C2k/C1k*k
        third, sixth = mpf(1)/3, mpf(1)/6
        coef = [I(beta+third,0,al) + rr*I(beta-2*third,0,al),               # J (l. 829-840): coefficient of 1, 1/LT, 1/LT^2, 1/(T^(1/6)LT), 1/(T^(1/6)LT^2), 1/(T^(1/3)LT^2)
                2*I(beta+third,1,al) + 2*rr*I(beta-2*third,1,al),
                I(beta+third,2,al) + rr*I(beta-2*third,2,al),
                2*a2*(I(beta+sixth,0,al) + rr*I(beta-5*sixth,0,al))/a1,
                2*a2*(I(beta+sixth,1,al) + rr*I(beta-5*sixth,1,al))/a1,
                a2**2*(I(beta,0,al) + rr*I(beta-1,0,al))/a1**2]
        def J(T):
            LT = log(T); T6 = T**sixth
            return coef[0] + coef[1]/LT + coef[2]/LT**2 + coef[3]/(T6*LT) + coef[4]/(T6*LT**2) + coef[5]/(T6**2*LT**2)
        om1 = lambda s, T: exp(al*(s/T)**2)                                  # l. 592
        U = lambda T: 4*al*beta*C4*om1(mpf(1)/2, T)**2*J(T)                  # l. 843-845
        C5 = pi*m0*b4/(2*de)*(1 + 2*de/LX0)**2*exp(2*de*g/LX0)               # l. 734-735
        C6 = b4/(5*de*exp(de))*(1 + de/LX0)**2 + b3*exp(-2*de)/LX0**2        # l. 737-738
        K = lambda T: (C5 + C6*pi*m0/(k*T))*I(1,0,al) + C6/k*I(2,0,al)       # l. 948-953
        V = lambda T: 8*al*om1(1 + de/log(k*T), T)**2*K(T)                   # l. 944-946
        Vinf = 8*al*(C5*I(1,0,al) + C6/k*I(2,0,al))
        M = max(L0/(LX0 + 2*de), mpf(1))                                     # l. 1014-1015
        b12 = 1/(2*(1 - 1/H)**2)                                             # l. 1654-1655
        def logC1(s):                                                        # l. 1576-1580
            e1 = mpf(8)/3*de*(2*s - 1)*M + 4*de*(2*s - 1)*ll0/(LX0 + 2*de)
            eU = 2*(1 - s) + 2*d/L0 + 2*de*(2*s - 1)/(LX0 + 2*de)
            e3 = 2*d*(2*ll0 - log(LX0))/L0 + mpf(8)*d/3 + 2*al
            return log(b12) + e1 + eU*log(U(H0)) + (2*s - 1)*log(V(H0)) + e3
        z = zeta(1 + eta)
        b5 = zeta(1 + eta)**4/zeta(2 + 2*eta)**2 + 2*zeta(1 + eta)**2/zeta(2 + 2*eta)   # l. 1291
        def b6(X, et):                                                       # l. 1304-1305
            LX = log(X); return (1 + et)*LX/(et*X**et)*(1 + 1/(et*LX) + g/LX + 7*et/(12*(1 + et)*X*LX))
        b7 = (1 + 2/(mpf('3.006')*z*X0**(1 + eta)))*(mpf('3.006')*z)**2       # l. 1345-1346
        b8 = sqrt((2 + eta)**2/H**2 + ((1 + 2*eta)/H + 1)**2)                # l. 1354-1355
        l2 = log(2)
        b9 = pi*log(b7)/l2 + pi*log(b5)/l2 - 2*pi*log(1 - b6(mpf(10)**9, eta)**2)/l2 + pi + 2*(1 + 2*eta)/l2*log(b8/(2*pi))   # l. 1431-1436
        C7 = (2*(1 + 2*eta) + 2*pi*(1 + eta))/l2 + b9/L0                     # l. 1405-1406
        bb6 = b6(X0, mu - 1); b10 = -log(1 - bb6**2)/bb6**2                 # l. 1497-1498
        b11 = lambda X, tau: 1 + 3/((tau - 1)*log(X)) + 6/((tau - 1)**2*log(X)**2) + 6/((tau - 1)**3*log(X)**3)   # l. 1526-1527
        C8 = b10*LX0**2/X0**(2*mu - 2)*(4*mu*b11(X0, 2*mu)/(k*(2*mu - 1)) + 2*pi*m0*(2*mu - 1)*b11(X0, 2*mu - 1)/(mu - 1))   # l. 1474-1476
        Bfun = lambda s: (C7*(mu - s + d/L0) + C8)/(2*pi*d)                  # l. 1581-1582
        Afun = lambda s: exp(logC1(s))/(2*pi*d)
        eta0 = findroot(lambda e: b6(mpf(10)**9, e) - 1, mpf('0.236'))
        adm = dict(U_H0=U(H0), U_gt1=U(H0) > 1, Jcoef_pos=all(c > 0 for c in coef), Vinf=Vinf, Vinf_gt1=Vinf > 1,
                   k_ok=(k >= mpf(10)**9/H0 and k <= 1), eta_ok=(eta0 <= eta < mpf(1)/2), mu_ok=(1 + eta0 <= mu <= 1 + eta),
                   delta_ok=(0 < de < L0*(ll0 - 1)/2), H_ok=(1002 <= H <= H0), b9_pos=b9 > 0, smin=mpf(1)/2 + d/L0)
        SIG = ['0.83', '0.86', '0.8733819666132076', r[0]]     # evaluated here: the closures bind loop variables late
        out.append(dict(row=r[0], A=Afun(s0), B=Bfun(s0), At=At, Bt=Bt, adm=adm, Av={s: Afun(mpf(s)) for s in SIG}, Bv={s: Bfun(mpf(s)) for s in SIG}, b9=b9, eta0=eta0))
    return out
H0a = mpf('3.0610046e10'); H0b = mpf('3e12'); ETA = mpf('0.2561882')
res = {}
P1 = run(H0a, ETA, ROWS, 'old')
worstA = max(abs(p['A']/p['At'] - 1) for p in P1); worstB = max(abs(p['B']/p['Bt'] - 1) for p in P1)
print('eta0 (root of b6(1e9, eta) = 1):', mp.nstr(P1[0]['eta0'], 12))
print('PART 1, H0 = 3.0610046e10: max |A/A_table - 1| = %.2e, max |B/B_table - 1| = %.2e over 20 rows; all admissible: %s' % (worstA, worstB, all(all(v for kk, v in p['adm'].items() if isinstance(v, bool)) for p in P1)))
for p in P1[:1] + P1[5:8] + P1[18:]:
    print('  row %s: A %s (table %s)  B %s (table %s)  U(H0) %s' % (p['row'], mp.nstr(p['A'], 7), p['At'], mp.nstr(p['B'], 7), p['Bt'], mp.nstr(p['adm']['U_H0'], 6)))
P2 = run(H0b, ETA, ROWS, 'new')
print('PART 2, H0 = 3e12 (H = H0 - 1, same Table 1 parameters):')
for p in P2:
    a = p['adm']; fails = [kk for kk, v in a.items() if isinstance(v, bool) and not v]
    print('  row %s: A(s0) %s  B(s0) %s  U(H0) %s  Vinf %s  b9 %s  smin %s  %s' % (p['row'], mp.nstr(p['A'], 8), mp.nstr(p['B'], 8), mp.nstr(a['U_H0'], 6), mp.nstr(a['Vinf'], 6), mp.nstr(p['b9'], 6), mp.nstr(a['smin'], 6), 'FAILS ' + ','.join(fails) if fails else 'admissible'))
# U(T) decreasing for row 0.60 between the two heights
p60 = P2[0]
print('row 0.60: U(3.0610046e10) = %s, U(3e12) = %s (same k, alpha; U depends on H0 only via C1(k), C2(k), C3(k))' % (mp.nstr(P1[0]['adm']['U_H0'], 8), mp.nstr(p60['adm']['U_H0'], 8)))
# compare with cert.py (copy) at H0 = 3e12: A(sigma), B(sigma) at several sigma
sys.argv = ['x']; sys.dont_write_bytecode = True
import importlib.util
spec = importlib.util.spec_from_file_location('c4', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Ingham4', 'cert.py'))   # read-only import, no bytecode written
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
PC = C.build_kln_H0('3e12')
names = [p.name for p in PC]
print('cert.py build_kln_H0(3e12) keeps rows:', len(PC), '; dropped:', [r[0] for r in ROWS if r[0] not in names])
worst = mpf(0)
for p in P2:
    if p['row'] not in names: continue
    q = PC[names.index(p['row'])]
    for s in ['0.83', '0.86', '0.8733819666132076', p['row']]:
        sv = mpf(s)
        Ac = exp(mpf(q.logA(C.arb(s)).mid().str(40, radius=False))); Bc = mpf(q.Bc(C.arb(s)).mid().str(40, radius=False))
        worst = max(worst, abs(Ac/p['Av'][s] - 1), abs(Bc/p['Bv'][s] - 1))
print('cert.py vs independent at H0 = 3e12, 19 rows x sigma in {0.83, 0.86, B, sigma0}: max relative difference in A or B = %.2e' % worst)
for nm in ['0.85', '0.87', '0.98']:
    p = P2[[r[0] for r in ROWS].index(nm)]; q = PC[names.index(nm)]
    print('  row %s at sigma0: A indep %s cert %s; B indep %s cert %s' % (nm, mp.nstr(p['A'], 10), q.logA(C.arb(nm)).exp().str(10, radius=False), mp.nstr(p['B'], 10), q.Bc(C.arb(nm)).str(10, radius=False)))

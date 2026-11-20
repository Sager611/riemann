# KLN Lemma 4.14 rows from the TeX (same formulas and line references as kln_indep.py), returning affine coefficients.
from mpmath import mp, mpf, sqrt, log, exp, pi, zeta, quad, inf, euler, findroot
def kln_rows(H0, eta):
    b1, b2, b3, b4 = mpf('0.62'), mpf('1.048'), mpf('0.605'), mpf('0.529'); a1, a2 = mpf('0.63'), mpf('2.851')
    m0 = sqrt(1 + mpf(2)/3*sqrt(mpf(6)/5)); g = +euler; beta = 2
    ROWS = [('0.60','0.5','1.251','0.288','0.3140','0.341'),('0.65','0.6','1.249','0.256','0.3070','0.340'),('0.70','0.8','1.247','0.222','0.3040','0.339'),
        ('0.75','1.0','1.245','0.189','0.3030','0.338'),('0.80','1.0','1.245','0.160','0.3030','0.337'),('0.85','1.0','1.245','0.133','0.3030','0.336'),
        ('0.86','1.0','1.245','0.127','0.3030','0.335'),('0.87','1.0','1.245','0.122','0.3030','0.335'),('0.88','1.0','1.245','0.116','0.3030','0.335'),
        ('0.89','1.0','1.245','0.111','0.3030','0.335'),('0.90','1.0','1.245','0.105','0.3030','0.334'),('0.91','1.0','1.245','0.100','0.3030','0.334'),
        ('0.92','1.0','1.245','0.095','0.3030','0.334'),('0.93','1.0','1.245','0.089','0.3030','0.333'),('0.94','1.0','1.245','0.084','0.3030','0.333'),
        ('0.95','1.0','1.245','0.079','0.3030','0.333'),('0.96','1.0','1.245','0.074','0.3030','0.332'),('0.97','1.0','1.245','0.069','0.3030','0.332'),
        ('0.98','1.0','1.245','0.064','0.3030','0.331'),('0.99','1.0','1.245','0.060','0.3030','0.331')]
    def I(A, n, al): return quad(lambda x: x**A*exp(-2*al*x**2)*log(x)**n, [0, mpf(1)/4, 1, 2, 4, inf])
    L0 = log(H0); ll0 = log(L0); H = H0 - 1; out = []
    for r in ROWS:
        k, mu, al, de, d = [mpf(v) for v in r[1:]]
        X0 = k*H0; LX0 = log(X0)
        C1k = 6/pi**2 + b2/LX0; C2k = pi*m0*b1/LX0 + 6*m0/(pi*X0) + pi*m0*b2/(X0*LX0)
        a3 = -6*a1/exp(1) + a2; C3 = a3**2*C1k*LX0; C4 = C1k*a1**2*(1 + 1/sqrt(C3))**2; rr = C2k/C1k*k
        t3, t6 = mpf(1)/3, mpf(1)/6
        coef = [I(beta+t3,0,al) + rr*I(beta-2*t3,0,al), 2*I(beta+t3,1,al) + 2*rr*I(beta-2*t3,1,al), I(beta+t3,2,al) + rr*I(beta-2*t3,2,al),
                2*a2*(I(beta+t6,0,al) + rr*I(beta-5*t6,0,al))/a1, 2*a2*(I(beta+t6,1,al) + rr*I(beta-5*t6,1,al))/a1, a2**2*(I(beta,0,al) + rr*I(beta-1,0,al))/a1**2]
        LT = L0; T6 = H0**t6
        J = coef[0] + coef[1]/LT + coef[2]/LT**2 + coef[3]/(T6*LT) + coef[4]/(T6*LT**2) + coef[5]/(T6**2*LT**2)
        UH0 = 4*al*beta*C4*exp(al*(mpf(1)/2/H0)**2)**2*J
        C5 = pi*m0*b4/(2*de)*(1 + 2*de/LX0)**2*exp(2*de*g/LX0); C6 = b4/(5*de*exp(de))*(1 + de/LX0)**2 + b3*exp(-2*de)/LX0**2
        VH0 = 8*al*exp(al*((1 + de/log(k*H0))/H0)**2)**2*((C5 + C6*pi*m0/(k*H0))*I(1,0,al) + C6/k*I(2,0,al))
        Vinf = 8*al*(C5*I(1,0,al) + C6/k*I(2,0,al))
        M = max(L0/(LX0 + 2*de), mpf(1)); b12 = 1/(2*(1 - 1/H)**2)
        def logC1(s):
            e1 = mpf(8)/3*de*(2*s - 1)*M + 4*de*(2*s - 1)*ll0/(LX0 + 2*de)
            eU = 2*(1 - s) + 2*d/L0 + 2*de*(2*s - 1)/(LX0 + 2*de)
            e3 = 2*d*(2*ll0 - log(LX0))/L0 + mpf(8)*d/3 + 2*al
            return log(b12) + e1 + eU*log(UH0) + (2*s - 1)*log(VH0) + e3
        z = zeta(1 + eta); b5 = z**4/zeta(2 + 2*eta)**2 + 2*z**2/zeta(2 + 2*eta)
        def b6(X, et):
            LX = log(X); return (1 + et)*LX/(et*X**et)*(1 + 1/(et*LX) + g/LX + 7*et/(12*(1 + et)*X*LX))
        b7 = (1 + 2/(mpf('3.006')*z*X0**(1 + eta)))*(mpf('3.006')*z)**2; b8 = sqrt((2 + eta)**2/H**2 + ((1 + 2*eta)/H + 1)**2); l2 = log(2)
        b9 = pi*log(b7)/l2 + pi*log(b5)/l2 - 2*pi*log(1 - b6(mpf(10)**9, eta)**2)/l2 + pi + 2*(1 + 2*eta)/l2*log(b8/(2*pi))
        C7 = (2*(1 + 2*eta) + 2*pi*(1 + eta))/l2 + b9/L0
        bb6 = b6(X0, mu - 1); b10 = -log(1 - bb6**2)/bb6**2
        b11 = lambda X, tau: 1 + 3/((tau - 1)*log(X)) + 6/((tau - 1)**2*log(X)**2) + 6/((tau - 1)**3*log(X)**3)
        C8 = b10*LX0**2/X0**(2*mu - 2)*(4*mu*b11(X0, 2*mu)/(k*(2*mu - 1)) + 2*pi*m0*(2*mu - 1)*b11(X0, 2*mu - 1)/(mu - 1))
        l0_, l1_ = logC1(mpf(0)), logC1(mpf(1))
        ok = UH0 > 1 and all(c > 0 for c in coef) and Vinf > 1 and 0 < de < L0*(ll0 - 1)/2 and mpf(10)**9/H0 <= k <= 1 and mu <= 1 + eta
        out.append(dict(name=r[0], k=k, c0=l0_ - log(2*pi*d), c1=l1_ - l0_, q0=(C7*(mu + d/L0) + C8)/(2*pi*d), q1=C7/(2*pi*d), smin=mpf(1)/2 + d/L0, ok=ok))
    return out

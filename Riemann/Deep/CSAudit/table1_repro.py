# Independent recomputation of Table 1 of Chourasiya-Simonic (arXiv 2507.15184 v2) from Thm 2.1, Lemmas 4.1-4.5, Cor 2.4
# (TeX l. 228-247, 256-259, 627-760, 346-356). Written from the TeX, not from Ingham4/cs_repro.py. mpmath, 30 digits.
import sys
import mpmath as mp
mp.mp.dps = 30
pi, EG, ln2 = mp.pi, mp.euler, mp.log(2)
b1, b2, b3 = mp.mpf('0.62'), mp.mpf('1.048'), mp.mpf('0.605')
m0 = mp.sqrt(1 + mp.mpf(2)/3*mp.sqrt(mp.mpf(6)/5))
HRH = mp.mpf('3e12'); T0 = H0 = HRH; h = mp.mpf(1)
M1_10, M2_10, MM2_10 = mp.mpf('34.842'), mp.mpf('1.3897e16'), mp.mpf('9.5079e15')
TAB = [  # s1, s2, B1, B2, B3, d, kappa, delta, log10 A0  (TeX l. 113-128)
 ('0.50000','0.53125','5.069','9.461','167.8','0.28','0.34','0.012','11.25'),
 ('0.53125','0.56250','6.046','9.168','162.6','0.28','0.32','0.025','11.26'),
 ('0.56250','0.59375','7.077','8.875','157.4','0.28','0.30','0.038','11.27'),
 ('0.59375','0.62500','8.185','8.582','152.2','0.28','0.29','0.052','11.27'),
 ('0.62500','0.65625','9.387','8.290','147.0','0.28','0.27','0.066','11.28'),
 ('0.65625','0.68750','10.70','7.997','141.9','0.28','0.25','0.081','11.29'),
 ('0.68750','0.71875','12.11','7.704','136.7','0.28','0.23','0.097','11.30'),
 ('0.71875','0.75000','13.66','7.411','131.5','0.28','0.21','0.114','11.31'),
 ('0.75000','0.78125','15.35','7.118','126.3','0.28','0.20','0.132','11.32'),
 ('0.78125','0.81250','17.21','6.826','121.1','0.28','0.18','0.150','11.33'),
 ('0.81250','0.84375','19.31','6.533','115.9','0.28','0.16','0.170','11.35'),
 ('0.84375','0.87500','21.77','6.240','110.7','0.28','0.14','0.192','11.36'),
 ('0.87500','0.90625','24.83','5.745','101.9','0.29','0.13','0.215','11.37'),
 ('0.90625','0.93750','29.01','5.463','96.87','0.29','0.11','0.241','11.40'),
 ('0.93750','0.96875','35.35','5.180','91.86','0.29','0.09','0.270','11.42'),
 ('0.96875','1.00000','46.06','4.897','86.84','0.29','0.07','0.303','11.46')]
def C1(X0): return 6/pi**2 + b2/mp.log(X0)
def C2(X0): return pi*m0*b1/mp.log(X0) + 6*m0/(pi*X0) + pi*m0*b2/(X0*mp.log(X0))
def C3(X0, de, b4): return pi*m0*b4/(2*de)*mp.exp(2*de*EG/mp.log(X0))
def C4(X0, de, b4): return b4/(5*de*mp.exp(de)) + b3*mp.exp(-2*de)/mp.log(X0)
def b6sq(u): return 9*(1 + (2 + EG + mp.mpf(7)/(36*u))/mp.log(u))**2*mp.log(u)**2/u
def C5(u, hh):
    q = 9*(1 + (2 + EG + mp.mpf(7)/(36*u))/mp.log(u))**2
    return -mp.log(1 - b6sq(u))*(3/hh + 8*pi*m0)*(1 + mp.mpf('3.31')/mp.log(u))/q*mp.log(u)
def C6(u, Hz):
    return (220 + pi/ln2*mp.log(1 + u**mp.mpf(-1.5)) + 4/ln2*mp.log(mp.sqrt((mp.mpf(5)/(2*Hz))**2 + (1 + 2/Hz)**2)/(2*pi))
            - 2*pi/ln2*mp.log(1 - b6sq(u)))
def a1(A0): return A0*mp.log(A0/2)**4/(2*pi**2) + M1_10*A0*mp.log(A0/2)**3 + M2_10
def a2(A0):
    L = mp.log(A0); return (1 - ln2/L)**4/(2*pi**2) + M1_10/L*(1 - ln2/L)**3 + MM2_10/(A0*L**4)
def K1(A0, ka):
    X0 = h*T0; c1, c2 = C1(X0), C2(X0); q = 3/(4*ka)
    f = lambda u: mp.exp(-u)*u**(mp.mpf(1)/6)*(1 + mp.log(1 + q*u)/(2*mp.log(T0)))**(mp.mpf(4)/3)*(c1*mp.sqrt(q*u) + h*c2)**(mp.mpf(2)/3)
    I = mp.quad(f, [0, 1, 10, 50, mp.inf])
    g = lambda v: 6*v**5*f(v**6)          # u = v^6: smooth at 0, cross-check
    I2 = mp.quad(g, [0, 1, mp.mpf(10)**(mp.mpf(1)/6), 2, mp.inf])
    assert abs(I - I2) < mp.mpf('1e-20')*I, (I, I2)
    return mp.cbrt(2*a2(A0)*mp.sqrt(q))*I + mp.cbrt(2*a1(A0))*(h*c2 + A0*c1/T0)**(mp.mpf(2)/3)/(mp.cbrt(T0)*mp.log(T0)**(mp.mpf(4)/3))
def K2(A0, ka): return mp.cbrt(2)*mp.sqrt(3/(4*ka))*mp.gamma(mp.mpf(1.5)) + mp.cbrt(2)*4*ka/3*(A0/T0)**3
def K3(ka, de, b4):
    X0 = h*T0; return C3(X0, de, b4) + C4(X0, de, b4)/h*(mp.gamma(mp.mpf(1.5))/mp.sqrt(2*ka) + pi*m0/T0)
def L0(ka, de): return mp.exp(ka*(1 + (1/(4*T0)*(1 + 2*de/mp.log(h*T0)))**2))/(1 - 1/H0)
def L1(k1, k3, de):
    e = 3*de/mp.log(h*T0); lt = mp.log(T0); llt = mp.log(lt)
    return max(1, k1**e)*max(1, k3**(-e))*mp.exp(3*de*max(1, lt/mp.log(h*T0)) + 4*de*max(llt/mp.log(h*T0), llt/lt))
def L2(k2, k3, de):
    e = 3*de/mp.log(h*T0); return max(1, k2**e)*max(1, k3**(-e))*mp.exp(3*de*max(1, mp.log(T0)/mp.log(h*T0)))
def row_consts(r, b4):
    s1, s2, d, ka, de = [mp.mpf(x) for x in (r[0], r[1], r[5], r[6], r[7])]
    A0 = mp.mpf(10)**mp.mpf(r[8])
    k1, k2, k3, l0 = K1(A0, ka), K2(A0, ka), K3(ka, de, b4), L0(ka, de)
    l1, l2 = L1(k1, k3, de), L2(k2, k3, de)
    lt = mp.log(T0); llt = mp.log(lt); e = lambda s: 3*(1 - s)/(2 - s); f = lambda s: (2*s - 1)/(2 - s)
    common = max(1, l0**(-2*d/lt))*max(1, k3**(-3*d/lt))*max(l0**(2/(2 - s1)), l0**(2/(2 - s2)))*max(k3**f(s1), k3**f(s2))
    A1 = l1*mp.exp(d*(3 + 4*llt/lt))*max(1, k1**(3*d/lt))*common*max(k1**e(s1), k1**e(s2))
    A2 = l2*mp.exp(3*d)*max(1, k2**(3*d/lt))*common*max(k2**e(s1), k2**e(s2))
    w = (3 - 2*s1)/(2*d) + 1/lt
    A3, A4 = w/(pi*ln2), 3*w/(2*ln2)
    A5 = C6(h*T0, H0)/(2*pi)*w + C5(h*T0, h)/(2*pi*d)
    B1 = (A1 + A2*mp.log(HRH)**(-6*(1 - s2)/(2 - s2)))/(2*pi*d)
    return dict(B1=B1, B2=A3 + A4, B3=A3*mp.log(HRH) + A5, K1=k1, K2=k2, K3=k3, L0=l0, L1=l1, L2=l2, A1=A1, A2=A2,
                C5=C5(h*T0, h), C6=C6(h*T0, H0), eps=2*de/mp.log(h*T0), ratio21=A2*mp.log(HRH)**(-6*(1 - s2)/(2 - s2))/A1)
if __name__ == '__main__':
    print('m0 = %s  C1 = %s  C2 = %s' % (mp.nstr(m0, 10), mp.nstr(C1(T0), 10), mp.nstr(C2(T0), 10)))
    print('row            B1tab   B1(0.54)  ratio   B2tab  B2      B3tab   B3      | eps=2d/logX  B1(0.5501) B1(b4pp)  B1(0.5936)  K1     K3(0.54)  A2/A1-share')
    for r in TAB:
        c = row_consts(r, mp.mpf('0.54'))
        eps = c['eps']
        cA = row_consts(r, mp.mpf('0.5501'))
        bpp = mp.mpf('0.5501')*(1 + (1 - EG)*eps)
        cB = row_consts(r, bpp)
        cC = row_consts(r, mp.mpf('0.5936'))
        print('%s-%s  %6s  %8.4f  %6.4f  %5s  %6.4f  %6s  %7.3f | %.5f      %8.4f  %8.4f  %8.4f  %7.3f  %6.3f  %.4f' % (
            r[0], r[1], r[2], c['B1'], c['B1']/mp.mpf(r[2]), r[3], c['B2'], r[4], c['B3'], eps, cA['B1'], cB['B1'], cC['B1'], c['K1'], c['K3'], c['ratio21']))
    sys.stdout.flush()

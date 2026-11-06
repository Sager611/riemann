#!/usr/bin/env python3
# cs_repro.py -- recompute the constants B1, B2, B3 of Chourasiya-Simonic (arXiv 2507.15184 v2) Corollary 1.1 / Table 1 from the
# formulas of their Theorem 2.1, Lemmas 4.1-4.5 and Corollary 2.4 (TeX l. 228-247, 249-275, 627-860), with their b4 = 0.54 and with the
# b4 supported by Ramare-Zuniga Alterman (arXiv 2603.25961 v3, Thm 1.2). A consistency check of their table (floats, 40 digits), not a proof.
# usage: uv run --no-project --with mpmath python3 Riemann/Deep/Ingham4/cs_repro.py
import mpmath as mp
mp.mp.dps = 40
pi, EG = mp.pi, mp.euler
b1, b2, b3 = mp.mpf('0.62'), mp.mpf('1.048'), mp.mpf('0.605')
m0 = mp.sqrt(1 + mp.mpf(2)/3*mp.sqrt(mp.mpf(6)/5))
HRH = mp.mpf('3e12'); T0 = H0 = HH = HRH; h = mp.mpf(1)
LT0 = mp.log(T0)
ROWS = [('0.50000','0.53125','5.069','9.461','167.8','0.28','0.34','0.012','11.25'), ('0.53125','0.56250','6.046','9.168','162.6','0.28','0.32','0.025','11.26'),
('0.56250','0.59375','7.077','8.875','157.4','0.28','0.30','0.038','11.27'), ('0.59375','0.62500','8.185','8.582','152.2','0.28','0.29','0.052','11.27'),
('0.62500','0.65625','9.387','8.290','147.0','0.28','0.27','0.066','11.28'), ('0.65625','0.68750','10.70','7.997','141.9','0.28','0.25','0.081','11.29'),
('0.68750','0.71875','12.11','7.704','136.7','0.28','0.23','0.097','11.30'), ('0.71875','0.75000','13.66','7.411','131.5','0.28','0.21','0.114','11.31'),
('0.75000','0.78125','15.35','7.118','126.3','0.28','0.20','0.132','11.32'), ('0.78125','0.81250','17.21','6.826','121.1','0.28','0.18','0.150','11.33'),
('0.81250','0.84375','19.31','6.533','115.9','0.28','0.16','0.170','11.35'), ('0.84375','0.87500','21.77','6.240','110.7','0.28','0.14','0.192','11.36'),
('0.87500','0.90625','24.83','5.745','101.9','0.29','0.13','0.215','11.37')]
def C1(X0): return 6/pi**2 + b2/mp.log(X0)
def C2(X0): return pi*m0*b1/mp.log(X0) + 6*m0/(pi*X0) + pi*m0*b2/(X0*mp.log(X0))
def C3(X0, de, b4): return pi*m0*b4/(2*de)*mp.exp(2*de*EG/mp.log(X0))
def C4(X0, de, b4): return b4/(5*de*mp.exp(de)) + b3*mp.exp(-2*de)/mp.log(X0)
def a1(A0): return A0*mp.log(A0/2)**4/(2*pi**2) + mp.mpf('34.842')*A0*mp.log(A0/2)**3 + mp.mpf('1.3897e16')
def a2(A0):
    l = mp.log(A0); return (1 - mp.log(2)/l)**4/(2*pi**2) + mp.mpf('34.842')/l*(1 - mp.log(2)/l)**3 + mp.mpf('9.5079e15')/(A0*l**4)
def K1(A0, ka):
    X0 = h*T0
    f = lambda u: mp.exp(-u)*u**(mp.mpf(1)/6)*(1 + mp.log(1 + 3*u/(4*ka))/(2*LT0))**(mp.mpf(4)/3)*(C1(X0)*mp.sqrt(3*u/(4*ka)) + h*C2(X0))**(mp.mpf(2)/3)
    I = mp.quad(f, [0, 1, 10, 60, mp.inf])
    return (2*a2(A0)*mp.sqrt(3/(4*ka)))**(mp.mpf(1)/3)*I + (2*a1(A0))**(mp.mpf(1)/3)*(h*C2(X0) + A0*C1(X0)/T0)**(mp.mpf(2)/3)/(T0**(mp.mpf(1)/3)*LT0**(mp.mpf(4)/3))
def K2(A0, ka): return mp.cbrt(2)*mp.sqrt(3/(4*ka))*mp.gamma(mp.mpf(3)/2) + mp.cbrt(2)*4*ka/3*(A0/T0)**3
def K3(ka, de, b4):
    X0 = h*T0; return C3(X0, de, b4) + C4(X0, de, b4)/h*(mp.gamma(mp.mpf(3)/2)/mp.sqrt(2*ka) + pi*m0/T0)
def L0(ka, de): return mp.exp(ka*(1 + (1/(4*T0)*(1 + 2*de/mp.log(h*T0)))**2))/(1 - 1/H0)
def L1(k1, k3, de):
    lh = mp.log(h*T0); return max(1, k1**(3*de/lh))*max(1, k3**(-3*de/lh))*mp.exp(3*de*max(1, LT0/lh) + 4*de*max(mp.log(LT0)/lh, mp.log(LT0)/LT0))
def L2(k2, k3, de):
    lh = mp.log(h*T0); return max(1, k2**(3*de/lh))*max(1, k3**(-3*de/lh))*mp.exp(3*de*max(1, LT0/lh))
def C5(u, hh):
    q = (1 + (2 + EG + 7/(36*u))/mp.log(u))**2
    return -mp.log(1 - 9*q*mp.log(u)**2/u)*((3/hh + 8*pi*m0)*(1 + mp.mpf('3.31')/mp.log(u)))/(9*q)*mp.log(u)
def C6(u, H0_):
    q = (1 + (2 + EG + 7/(36*u))/mp.log(u))**2
    return (220 + pi/mp.log(2)*mp.log(1 + u**(-mp.mpf(3)/2)) + 4/mp.log(2)*mp.log(mp.sqrt((5/(2*H0_))**2 + (1 + 2/H0_)**2)/(2*pi))
            - 2*pi/mp.log(2)*mp.log(1 - 9*q*mp.log(u)**2/u))
def consts(row, b4):
    s1, s2, B1t, B2t, B3t, d, ka, de, lA0 = [mp.mpf(v) for v in row]
    A0 = mp.mpf(10)**lA0
    k1, k2, k3 = K1(A0, ka), K2(A0, ka), K3(ka, de, b4)
    l0, l1, l2 = L0(ka, de), L1(k1, k3, de), L2(k2, k3, de)
    e = lambda s: (2*s - 1)/(2 - s); g = lambda s: 3*(1 - s)/(2 - s)
    common = max(1, l0**(-2*d/LT0))*max(1, k3**(-3*d/LT0))*max(l0**(2/(2 - s1)), l0**(2/(2 - s2)))*max(k3**e(s1), k3**e(s2))
    A1 = l1*mp.exp(d*(3 + 4*mp.log(LT0)/LT0))*max(1, k1**(3*d/LT0))*common*max(k1**g(s1), k1**g(s2))
    A2 = l2*mp.exp(3*d)*max(1, k2**(3*d/LT0))*common*max(k2**g(s1), k2**g(s2))
    Bb1 = (A1 + A2*mp.log(HRH)**(-6*(1 - s2)/(2 - s2)))/(2*pi*d)
    w = (3 - 2*s1)/(2*d) + 1/LT0
    A3 = w/(pi*mp.log(2)); A4 = 3*w/(2*mp.log(2)); A5 = C6(h*T0, H0)/(2*pi)*w + C5(h*T0, h)/(2*pi*d)
    return dict(B1=Bb1, B2=A3 + A4, B3=A3*mp.log(HRH) + A5, K1=k1, K2=k2, K3=k3, L0=l0, B1t=B1t, B2t=B2t, B3t=B3t, s2=s2, de=de)
if __name__ == '__main__':
    print('row            B1(table) B1(b4=0.54) ratio   B2(table) B2(ours)  B3(table) B3(ours)  K1      K3     | b4 supported  B1(recomputed)  B1(scaled, cert.py)')
    for row in ROWS:
        c = consts(row, mp.mpf('0.54'))
        eps = 2*c['de']/LT0
        b4s = mp.mpf('0.5501')*(1 + (1 - EG)*eps)
        cs = consts(row, b4s)
        scaled = c['B1t']*(b4s/mp.mpf('0.54'))**((2*c['s2'] - 1)/(2 - c['s2']))
        print('%s-%s  %7.3f  %9.4f  %6.4f  %7.3f  %8.4f  %7.2f  %8.3f  %6.3f  %6.3f | %.5f   %9.4f   %9.4f' % (row[0], row[1], c['B1t'], c['B1'], c['B1']/c['B1t'], c['B2t'], c['B2'], c['B3t'], c['B3'], c['K1'], c['K3'], b4s, cs['B1'], scaled))

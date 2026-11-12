# Recompute F1(T0'', c, sigma, a, b) of Chourasiya-Simonic Thm 2.2 (TeX l. 296-312) at the printed Table 3 parameters (l. 396-406)
# from the J-formulas (cJ0, cJ1, cJ21, cJ22 = G1..G6 via g1..g8, cJ3..cJ6; l. 1300-1650) and compare with Table 2's M1(T0) (l. 186-194).
# Then redo the recursion of l. 418-427 for M2(10^n, 2 10^n) and frak M2(10^n). mpmath 25 digits.
import mpmath as mp
mp.mp.dps = 25
pi, g, ln = mp.pi, mp.euler, mp.log
m0 = mp.sqrt(1 + mp.mpf(2)/3*mp.sqrt(mp.mpf(6)/5))
def eta(T0): return 1 + pi*m0/T0
def Um(y): return 1/(1 + mp.exp(-pi*y))*mp.exp(-1/(2*y**2))
def Up(y): return 1/(1 - mp.exp(-pi*y))*mp.exp(1/(2*y**2))
def D1(x0, s0):
    l = ln(x0); r = mp.sqrt(x0)
    return 1 + (2*g - 1 - s0/(1 - s0))/l + (1 - s0)*(mp.mpf(1)/2 + abs(s0) - s0)/r + (1 - 2*g*s0)/(r*l)
def D2(x0, s0):
    l = ln(x0); return s0 - 1 + (1 + 2*g*(s0 - 1))/l + (4*s0 - 1)/((2*s0 - 1)*mp.sqrt(x0)*l)
def D3(x0):
    l = ln(x0); e = mp.e; return 1/e + (1 + 2*g)/(e*l) + 5/(2*e*mp.sqrt(x0)*l)
def D4(x0, s):
    l = ln(x0)
    assert -3 <= s < 1, s
    if s < 0:
        return (1/((1 - s)*pi**2) + (-3*s/((1 - s)**2*pi**2) + mp.mpf('0.744342')/(1 - s))/l
                + (6*s/((1 - s)**3*pi**2) - mp.mpf('1.488684')*s/(1 - s)**2 + mp.mpf('0.823266')/(1 - s))/l**2
                + (-6*s/((1 - s)**4*pi**2) + mp.mpf('1.488684')*s/(1 - s)**3 - mp.mpf('0.823266')*s/(1 - s)**2 + mp.mpf('0.460324')/(1 - s))/l**3
                + mp.mpf('17.52')/(x0**mp.mpf('0.25')*l**2))
    d42 = (1/pi**2 + mp.mpf('0.744342')/l + mp.mpf('0.823266')/l**2 + mp.mpf('0.460324')/l**3
           + mp.mpf('9.73')*s*(1 - s)/(x0**min(mp.mpf('0.25'), 1 - s)*l) + mp.mpf('9.73')*(1 - s)/(x0**mp.mpf('0.25')*l**2))
    return d42/(1 - s)
def D5(x0, s):
    l = ln(x0); a = s - 1; assert s > 1
    return (a**3/pi**2 + a**2*(3*s/pi**2 + mp.mpf('0.744342')*a)/l + a*(6*s/pi**2 + mp.mpf('1.488684')*s*a + mp.mpf('0.823266')*a**2)/l**2
            + (6*s/pi**2 + mp.mpf('1.488684')*s*a + mp.mpf('0.823266')*s*a**2 + mp.mpf('0.460324')*a**3)/l**3
            + (1 + s/(s - mp.mpf('0.75'))*(1 + 1/((s - mp.mpf('0.75'))*l)))*mp.mpf('9.73')*a**4/(x0**mp.mpf('0.25')*l**2))
def D6(x0):
    l = ln(x0); return mp.mpf('0.349436') + mp.mpf('1.155975')/l + mp.mpf('1.28359')/l**2 + 154/l**3 + mp.mpf('9.73')/(x0**mp.mpf('0.25')*l**2)
def D7(x0):
    l = ln(x0); return 3/(2*mp.e*pi)**2 + mp.mpf('0.12441')/l + mp.mpf('0.364')/l**2 + mp.mpf('0.3')/l**3 + mp.mpf('4.61')/(x0**mp.mpf('0.25')*l**2)
_IG = {}
def IG(c):   # int_R |Gamma(c + iu)| du
    if c not in _IG: _IG[c] = 2*mp.quad(lambda u: abs(mp.gamma(mp.mpc(c, u))), [0, 1, 5, 20, 60])
    return _IG[c]
def F1(T0, c1, c2, s1, s2, s3, s4, s5, a1, a2, b1, b2):
    T0 = mp.mpf(T0); c1, c2, s1, s2, s3, s4, s5, a1, a2, b1, b2 = [mp.mpf(v) for v in (c1, c2, s1, s2, s3, s4, s5, a1, a2, b1, b2)]
    e, U, L = eta(T0), Um(T0), ln(T0); A0, B0 = a1*L**a2, b1*L**b2; ac = abs(c1)
    # conditions (2)-(6) of Thm 2.2
    ok = (-1 < c1 < -mp.mpf(1)/2 and 0 < c2 < mp.mpf(1)/2 and 0 <= s1 < mp.mpf(1)/2 and mp.mpf(1)/2 < s2 <= 1 and 0 <= s3 < mp.mpf(1)/2
          and mp.mpf(1)/2 < s4 < ac and mp.mpf(1)/2 - c2 < s5 < mp.mpf(1)/2 and a1 > 0 and a2 > 1 and b1 > 0 and b2 > 1
          and T0 - A0 >= 2 and a1*L**(a2 - 1) >= 1 and T0 - B0 >= 2 and b1*L**(b2 - 1) >= 2 and s4 + c2 < 1)
    UpA, UpB = Up(T0 - A0), Up(T0 - B0); IG1, IG2 = IG(c1), IG(c2)
    J1 = 2*e*D6(T0) + 4*pi*m0*D4(T0, 0) + m0*L/(2*pi*T0)
    J21 = (2*(2*pi)**(1 - 2*s1)/U**2*(e*D4(T0, 2*s1) + 2*pi*m0*D4(T0, 2*s1 - 1))
           + 2/L*(D1(T0, 0)/U)**2*mp.quad(lambda u: ((2*pi)**(1 - 2*u) + pi**(1 - 2*u))/(1 - u)**2, [s1, mp.mpf(1)/2]))
    q7 = mp.sqrt((e + 2*pi*m0)*D7(T0))
    G1 = (1/U**2*(2 + 1/T0)**(s2 - mp.mpf(1)/2)*mp.sqrt(e*D4(T0, 2 - 2*s2) + 2*pi*m0*D4(T0, 1 - 2*s2))
          *(q7 + mp.sqrt(e*D4(T0, 2*s2 - 2) + 2*pi*m0*D4(T0, 2*s2 - 3)))
          + (2*pi)**(1 - 2*s3)/U**2*mp.sqrt(e*D4(T0, 2*s3) + 2*pi*m0*D4(T0, 2*s3 - 1))
          *(q7 + mp.sqrt(e*D4(T0, 2*s3 - 2) + 2*pi*m0*D4(T0, 2*s3 - 3))))
    h = lambda u: D1(T0, u - 1)/(2 - u) + D3(T0)
    G2 = D1(T0, 0)/U**2*(mp.quad(lambda u: ((1 + 1/T0)**(u - mp.mpf(1)/2) + (1 + 1/(2*T0))**(u - mp.mpf(1)/2))*h(u)/u, [mp.mpf(1)/2, s2])
                          + mp.quad(lambda u: ((2*pi)**(1 - 2*u) + pi**(1 - 2*u))*h(u)/(1 - u), [s3, mp.mpf(1)/2]))
    def g1(s4_):
        return (1/U**2*(1/pi + 1/(2*pi*T0))**(2*s4_ - 1)*(1/pi + (mp.mpf(1)/2 + A0)/(2*pi*T0))**(2*(1 - 2*s4_) + 4*ac)*UpA**4*IG1**2
                *(e*D5(T0, 2 - 2*s4_ + 2*ac)/((2*pi)**2*(1 - 2*s4_ + 2*ac)**4) + m0*D5(T0, 1 - 2*s4_ + 2*ac)/(2*pi*(2*ac - 2*s4_)**4)))
    def g2(s4_):
        return (1/U**2*(1/pi + 1/(2*pi*T0))**(2*s4_ - 1)*(1/(2*pi) - B0/(2*pi*T0))**(2*(1 - 2*(s4_ + c2)))*UpB**4*IG2**2
                *(e*D4(T0, 2*(1 - s4_ - c2))/(2*pi)**2 + m0*D4(T0, 1 - 2*(s4_ + c2))/(2*pi)))
    g3 = (1/((2*pi)**(1 + 2*s5)*U**2)*(1/pi + (1 + A0)/(2*pi*T0))**(2*(1 - 2*s5) + 4*ac)*UpA**4*IG1**2
          *(e*D5(T0, 2 - 2*s5 + 2*ac)/(1 - 2*s5 + 2*ac)**4 + 2*pi*m0*D5(T0, 1 - 2*s5 + 2*ac)/(2*ac - 2*s5)**4))
    g4 = (1/((2*pi)**(1 + 2*s5)*U**2)*(1/(2*pi) - B0/(2*pi*T0))**(2*(1 - 2*(s5 + c2)))*UpB**4*IG2**2
          *(e*D4(T0, 2*(1 - s5 - c2)) + 2*pi*m0*D4(T0, 1 - 2*(s5 + c2))))
    g5 = mp.quad(lambda s: ((1 + 1/T0)**(s - mp.mpf(1)/2)*(1/(2*pi) + (mp.mpf(1)/2 + A0)/(2*pi*T0))**(1 - 2*(s - ac))
                            + (1 + 1/(2*T0))**(s - mp.mpf(1)/2)*(1/pi + (mp.mpf(1)/2 + A0)/(2*pi*T0))**(1 - 2*(s - ac)))*D2(T0, 1 - s + ac)/(s*(s - ac)**2), [mp.mpf(1)/2, s4])
    g6 = mp.quad(lambda s: ((1 + 1/T0)**(s - mp.mpf(1)/2)*(1/(2*pi) - B0/(2*pi*T0))**(1 - 2*(s + c2))
                            + (1 + 1/(2*T0))**(s - mp.mpf(1)/2)*(1/pi - B0/(2*pi*T0))**(1 - 2*(s + c2)))/(s*(s + c2)), [mp.mpf(1)/2, s4])
    g7 = mp.quad(lambda s: ((2*pi)**(1 - 2*s)*(1/(2*pi) + (1 + A0)/(2*pi*T0))**(1 - 2*(s - ac))
                            + pi**(1 - 2*s)*(1/pi + (1 + A0)/(2*pi*T0))**(1 - 2*(s - ac)))*D2(T0, 1 - s + ac)/((1 - s)*(s - ac)**2), [s5, mp.mpf(1)/2])
    g8 = mp.quad(lambda s: ((2*pi)**(1 - 2*s)*(1/(2*pi) - B0/(2*pi*T0))**(1 - 2*(s + c2))
                            + pi**(1 - 2*s)*(1/pi - B0/(2*pi*T0))**(1 - 2*(s + c2)))/((1 - s)*(s + c2)), [s5, mp.mpf(1)/2])
    G3 = (2*pi)**s4/U*mp.sqrt(e/(2*pi)*D4(T0, 2*(1 - s4)) + m0*D4(T0, 1 - 2*s4))*(mp.sqrt(g1(s4)) + mp.sqrt(g2(s4)))
    G4 = (2*pi)**(mp.mpf(1)/2 - s5)/U*mp.sqrt(e*D4(T0, 2*s5) + 2*pi*m0*D4(T0, 2*s5 - 1))*(mp.sqrt(g3) + mp.sqrt(g4))
    G5 = D1(T0, 0)/(2*pi)*(UpA/U)**2*g5*IG1 + D1(T0, 0)*D1(T0, -mp.mpf(1)/2)/(2*pi)*(UpB/U)**2*g6*IG2
    G6 = D1(T0, 0)/(2*pi)*(UpA/U)**2*g7*IG1 + 1/(2*pi)*(D1(T0, 0)*UpB/U)**2*g8*IG2
    J22 = G1 + G3 + G4 + (G2 + G5 + G6)/L
    J3 = (e + 2*pi*m0)*D7(T0); J4 = e*D4(T0, -1) + 2*pi*m0*D4(T0, -2); J5 = g1(mp.mpf(1)/2); J6 = g2(mp.mpf(1)/2)
    J01 = (2**(2*ac)*mp.sqrt(2/pi)*(mp.mpf('0.4') + 16/(2*T0 + A0)**ac)**2*mp.exp(1/(2*A0**2))*D2(T0, mp.mpf(1)/2 + ac)/(ac - mp.mpf(1)/2)**2
           *mp.quad(lambda u: (1 + u/(2*T0))**(2*ac)*mp.exp(-u/2)/u**(ac + mp.mpf(1)/2), [4, 20, 100, mp.inf]))
    J02 = 2/mp.e**4*mp.sqrt(2/pi)*(mp.mpf('6.31')/(mp.mpf(1)/2 - c2))**2*mp.exp(1/(2*B0**2))*B0**(c2 - mp.mpf(1)/2)*D1(T0, 0)/(mp.mpf(1)/2 + c2)
    J0 = mp.mpf(7)/2*((1 + mp.mpf(3)/14*J01)*J01*T0**(1 - a1*L**(a2 - 1)) + (1 + mp.mpf(3)/14*J02)*J02*T0**(2 - b1*L**(b2 - 1)) + mp.mpf('12.1')*T0*mp.exp(-pi*T0/2))
    Js = [J3, J4, J5, J6]
    mixed = 2*mp.fsum(mp.sqrt(Js[i]*Js[j]) for i in range(4) for j in range(i + 1, 4))
    F = J0 + J1 + J21 + 2*J22 + mp.fsum(Js) + mixed
    return F, ok, dict(J0=J0, J1=J1, J21=J21, J22=J22, J3=J3, J4=J4, J5=J5, J6=J6, mixed=mixed, G1=G1, G2=G2, G3=G3, G4=G4, G5=G5, G6=G6)
ROWS = [  # T0, c, sigma, a, b (TeX l. 396-404), printed M1 (l. 186-194)
 (1e4,  -0.879, 0.1809,  0.3545, 0.7898, 0.4008, 0.5456, 0.3933, 0.339, 1.7553, 1.7215, 1.303, '59.300'),
 (1e5,  -0.938, 0.194,   0.368,  0.718,  0.435,  0.629,  0.3654, 1.99, 1.44, 1.33, 1.45, '48.801'),
 (1e6,  -0.9295, 0.1721, 0.4077, 0.6957, 0.3413, 0.5672, 0.4288, 1.6769, 1.5632, 1.6891, 1.6748, '43.230'),
 (1e7,  -0.9273, 0.2061, 0.4104, 0.6263, 0.3536, 0.5743, 0.41, 1.8473, 1.8288, 1.649, 1.9094, '39.720'),
 (1e8,  -0.9266, 0.1637, 0.3561, 0.7103, 0.3381, 0.5645, 0.4285, 0.6164, 1.9123, 1.0018, 1.9929, '37.550'),
 (1e9,  -0.834, 0.22985, 0.3964, 0.9479, 0.3927, 0.5846, 0.4235, 1.3609, 1.1367, 0.7112, 1.8455, '36.211'),
 (1e10, -0.8303, 0.1676, 0.3985, 0.668,  0.4088, 0.551,  0.386,  1.4356, 1.28, 1.5563, 1.6332, '34.842'),
 (1e15, -0.813, 0.221,   0.3464, 0.8268, 0.4263, 0.5698, 0.3511, 1.0851, 1.0245, 0.8137, 1.8276, '31.170'),
 (1e20, -0.875, 0.1936,  0.3785, 0.8615, 0.45,   0.619,  0.393,  0.1917, 1.525, 1.1358, 1.717, '29.364')]
if __name__ == '__main__':
    M1 = {}
    for r in ROWS:
        T0 = mp.mpf(10)**int(round(mp.log10(r[0])))
        F, ok, parts = F1(T0, *[str(v) for v in r[1:12]])
        M1[int(round(mp.log10(r[0])))] = F
        print('T0 = 1e%-2d  F1 = %.5f  printed M1 = %s  ratio %.6f  conditions(2)-(6) & s4+c2<1: %s' % (round(mp.log10(r[0])), float(F), r[12], float(F/mp.mpf(r[12])), ok))
        print('      ' + '  '.join('%s=%.4g' % (k, float(v)) for k, v in parts.items()))

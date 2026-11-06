# Re-derive the coefficients of arXiv 2602.05614 (CorputSum1 for region 1 and 10, CorputSum2 at t0=8.97e17)
# with the paper's own formulas (TeX l. 590-945, 1083-1230). mpmath, 50 digits. Sanity check only.
from mpmath import mp, mpf, sqrt, log, pi, asinh, ceil, floor, zeta
mp.dps = 50
def H(N):  # sum_{n=1}^{N} n^{-1/2}, exact via Hurwitz zeta
    return zeta(mpf(1)/2) - zeta(mpf(1)/2, N+1)
def consts(t0, t1, c, phi, r0=4):
    t0, t1, c, phi = mpf(t0), mpf(t1), mpf(c), mpf(phi)
    R0 = int(ceil((sqrt(t0/(2*pi)) - 1)/(t0**phi + 1) - 1))
    R1 = int(floor(t1**(mpf(1)/2-phi)/sqrt(2*pi)))
    mu0 = (1 + mpf(1)/r0)*(2 + mpf(1)/r0 + 1/(t0**phi - 1))/2
    e = t0**(phi - mpf(1)/3)
    M1 = 1/c + mu0*(32/(15*sqrt(pi))*sqrt(c + 1/(5*pi**(mpf(1)/3)*e)) + c/(15*pi**(mpf(1)/3)*e) + 1/(75*pi**(mpf(2)/3)*e**2))
    M2 = (1/t0**(phi-mpf(1)/6))/(sqrt(2)*pi**(mpf(1)/6))*(1+mpf(1)/R0) + mu0*c*(16*sqrt(2)/(15*pi**(mpf(2)/3))/t0**(phi-mpf(1)/6)*(1+mpf(1)/R0)*sqrt(c + 1/(5*pi**(mpf(1)/3)*e)) + c/(3*t0**phi) + 1/(15*pi**(mpf(1)/3)*t0**(2*phi-mpf(1)/3)))
    M3 = (1/t0**phi)/c*(6/pi-1) + 16/(3*pi**(mpf(5)/6)*sqrt(c))/t0**(phi-mpf(1)/3)*(1+mpf(1)/R0)**2*(1+1/t0**phi) + (1-6/pi)*5*pi**(mpf(1)/3)/t1**(mpf(1)/3)
    M4 = (1/t0**(2*phi-mpf(1)/6))/(sqrt(2)*pi**(mpf(1)/6))*(6/pi-1)*(1+mpf(1)/R0) + (1/t0**(2*phi-mpf(1)/2))*8*sqrt(2*c)/(3*pi)*(1+mpf(1)/R0)**3*(1+1/t0**phi) + (1-6/pi)/t1**(phi+mpf(1)/6)*5*pi**(mpf(1)/6)*c/sqrt(2)*(1+mpf(1)/R1)
    M5 = log(1+1/(2*mpf(R0))) + 2*log(1+sqrt(1+1/(R0+mpf(1)/2))) + 1/sqrt(mpf(r0*(r0+1))) - log(sqrt(2*pi)) - 2*asinh(sqrt(r0+mpf(1)/2))
    M6 = 4*sqrt(r0*(1+1/t0**phi))
    N0 = int(ceil(t0**phi*r0))
    M7 = 2*H(N0-1) - 4*sqrt(mpf(N0-1))
    return dict(R0=R0, R1=R1, mu0=mu0, M1=M1, M2=M2, M3=M3, M4=M4, M5=M5, M6=M6, M7=M7)
def corput1(t0, t1, c, phi):
    d = consts(t0, t1, c, phi); t0 = mpf(t0)
    S = sqrt(4/pi**(mpf(1)/3)*(d['M1']+d['M2']+d['M3']+d['M4']))
    a1 = (mpf(1)/2-mpf(phi))*S; a2 = d['M5']*S; a3 = d['M6']; a4 = d['M7'] + mpf('1.463')*t0**(-mpf(1)/4) + mpf('0.127')*t0**(-mpf(3)/4)
    return a1, a2, a3, a4
def corput2(t0, c):
    phi = mpf(1)/3; d = consts(t0, t0, c, phi); t0 = mpf(t0); c = mpf(c)
    M8 = sqrt(5)*pi**(mpf(1)/6)*sqrt((32/(3*sqrt(pi*c)) + (6/pi-1)/c/(25*pi**(mpf(2)/3)))*(1 + 1/t0**(mpf(1)/3) + 5*c/(4*sqrt(2)*pi**(mpf(1)/6)*t0**(mpf(1)/6))))
    b1 = sqrt(d['M1']+d['M2'])/(3*pi**(mpf(1)/6))
    b2 = d['M6'] + 2*d['M5']/pi**(mpf(1)/6)*sqrt(d['M1']+d['M2']) + M8/sqrt(2*pi)
    b3 = d['M7'] - 3*M8 + mpf('1.463')*t0**(-mpf(1)/4) + mpf('0.127')*t0**(-mpf(3)/4)
    return b1, b2, b3, M8
def ratio(coef, t, phi):
    t = mpf(t); a1,a2,a3,a4 = coef
    return (a1*t**(mpf(1)/6)*log(t) + a2*t**(mpf(1)/6) + a3*t**(mpf(phi)/2) + a4)/(t**(mpf(1)/6)*log(t))
regions = [(3.80e7,3.85e7,1.82,0.33559,(0.622447,-8.21438,8.01143,-2.87533)),
           (3.85e7,3.96e7,1.82,0.33576,(0.621314,-8.20793,8.01135,-2.8755)),
           (3.96e7,4.12e7,1.82,0.33599,(0.619747,-8.19871,8.01119,-2.8758)),
           (4.12e7,4.43e7,1.82,0.33643,(0.61681,-8.1818,8.01096,-2.87626)),
           (4.43e7,4.96e7,1.82,0.33711,(0.612254,-8.15526,8.01057,-2.87706)),
           (4.96e7,5.82e7,1.81,0.33816,(0.605265,-8.11448,8.00999,-2.87828)),
           (5.82e7,7.06e7,1.81,0.33962,(0.595642,-8.05816,8.00922,-2.87995)),
           (7.06e7,8.24e7,1.80,0.34098,(0.58663,-8.00412,8.00842,-2.8818)),
           (8.24e7,3.92e8,1.83,0.33536,(0.617933,-8.20176,8.00885,-2.88182)),
           (3.92e8,8.97e17,1.83,0.33711,(0.597618,-8.09885,8.00507,-2.89251))]
print("region  ours(a1,a2,a3,a4)  paper  max ratio at t0,t1")
for (t0,t1,c,phi,paper) in regions:
    co = corput1(t0,t1,c,phi)
    print("%.3g-%.3g" % (t0,t1), " ".join("%.6g" % float(x) for x in co), "| paper", paper, "| ratio t0 %.5f t1 %.5f" % (float(ratio(co,t0,phi)), float(ratio(co,t1,phi))))
b1,b2,b3,M8 = corput2(8.97e17, 1.34)
print("CorputSum2 t0=8.97e17 c=1.34: b1=%.6g b2=%.6g b3=%.6g (paper 0.469028 3.99399 -21.4623) M8=%.6g" % (b1,b2,b3,M8))
t = mpf('8.97e17'); print("  ratio at t0: %.6f" % float(b1 + b2/log(t) + b3/(t**(mpf(1)/6)*log(t))))
b1,b2,b3,M8 = corput2(7.00e11, 1.34)
print("CorputSum2 t0=7e11 c=1.34: b1=%.6g b2=%.6g b3=%.6g (paper 0.470795 4.04972 -21.5437)" % (b1,b2,b3))
t = mpf('7e11'); print("  ratio at t0: %.6f" % float(b1 + b2/log(t) + b3/(t**(mpf(1)/6)*log(t))))
# RS trivial bound ratio at 3.3e7 and 3.8e7
for t,cst in ((3.3e7,2.0895),(3.8e7,2.8805)):
    t = mpf(t); print("RS trivial ratio at %.3g: %.6f" % (t, float((4*t**(mpf(1)/4)/(2*pi)**(mpf(1)/4) - cst)/(t**(mpf(1)/6)*log(t)))))

# Independent evaluation of Bellotti-Wong's leading constant C1 = Cbar1/(2 pi log(r/(c-1/2))) (BW arXiv 2412.15470,
# Sec. 3.4.1, transcribed from the TeX, not from Ismail's code) at the nine rows of Ismail's Table 1 (n = 4).
import mpmath as mp
mp.mp.dps = 30
ROWS = [('HPY',1.1948780,1.667319,0.194877,0.14074976),('HPY',1.1688080,1.581471,0.168807,0.13479015),('HPY',1.1450330,1.501535,0.145032,0.12934985),
        ('HPY',1.1150630,1.395716,0.115062,0.12238338),('HPY',1.0994560,1.339248,0.099455,0.11879775),('HPY',1.0832160,1.279677,0.083215,0.11514373),
        ('SW',1.0660270,1.215723,0.066026,0.11125810),('SW',1.0677880,1.135577,0.067787,0.10810835),('SW',1.0428930,1.085787,0.042892,0.10523393)]
for reg, c, r, eta, C1paper in ROWS:
    c, r, eta = mp.mpf(c), mp.mpf(r), mp.mpf(eta); n = 4
    k2 = mp.mpf(1)/6 if reg == 'HPY' else mp.mpf(27)/164
    th = lambda y: mp.acos((y - c)/r)
    sg = lambda w: c + r*mp.cos(w)
    S = lambda k: 1 - mp.mpf(k)/(2**k - 2)
    q = lambda f, a, b: mp.quad(f, [a, b])
    N4 = 2**(n+4)
    C = q(lambda w: (N4-1)*(1-sg(w)) + (N4-2)*(sg(w)-S(n+4)), th(1), th(S(n+4)))/((N4-2)*(1-S(n+4)))
    for h in range(n):
        A, B = 2**(h+4)-2, 2**(h+5)-2
        C += q(lambda w: B*(A+1)*(S(5+h)-sg(w)) + A*(B+1)*(sg(w)-S(4+h)), th(S(5+h)), th(S(4+h)))/(A*B*(S(5+h)-S(4+h)))
    C += q(lambda w: 14*(k2+1)*(S(4)-sg(w)) + 15*(sg(w)-mp.mpf(1)/2), th(S(4)), th(mp.mpf(1)/2))/(14*(S(4)-mp.mpf(1)/2))
    C += q(lambda w: (1-2*sg(w)+4*k2*sg(w))/2, th(mp.mpf(1)/2), th(0))
    C += q(lambda w: (1-2*sg(w))/2, th(-eta), mp.pi)
    C += th(1) - th(mp.mpf(1)/2)
    C += q(lambda w: -sg(w)*(1+2*eta)/(2*eta) + (sg(w)+eta)/(2*eta), th(0), th(-eta))
    C1 = C/(2*mp.pi*mp.log(r/(c-mp.mpf(1)/2)))
    print(reg, float(c), float(r), float(eta), 'C1 indep = %.8f  paper raw %.8f  diff %.1e' % (C1, C1paper, C1-C1paper))

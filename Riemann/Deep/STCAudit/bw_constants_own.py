# Own transcription of Bellotti-Wong's constants C2 (and C1) and C3, C~3 (arXiv 2412.15470v2 TeX, Sec. 3.4, l. 838-912, 943)
# at Ismail's nine rows; compare with Ismail's raw values; then the corrected D3 with zeta(1+eta) in place of (1+eta).
import mpmath as mp
mp.mp.dps = 30
f = mp.mpf
T0 = f(30610046000); J1, J2 = 64, 39
HPY = (f('0.618'), f(1)/6, f(1)); SW = (f('66.7'), f(27)/164, f(0))
ROWS = [(HPY,1.1948780,1.667319,0.194877,0.15095051,1.37729567),(HPY,1.1688080,1.581471,0.168807,0.15960505,1.50581607),
        (HPY,1.1450330,1.501535,0.145032,0.16931805,1.64349485),(HPY,1.1150630,1.395716,0.115062,0.18602913,1.85901425),
        (HPY,1.0994560,1.339248,0.099455,0.19717111,1.99600992),(HPY,1.0832160,1.279677,0.083215,0.21117818,2.16357331),
        (SW,1.0660270,1.215723,0.066026,0.15456953,2.73093591),(SW,1.0677880,1.135577,0.067787,0.20755749,2.79552145),
        (SW,1.0428930,1.085787,0.042892,0.21440763,3.17543128)]
Q0 = f(10)**9; c1 = f('0.5') + f('0.6633')/mp.log(Q0); c2 = f(1); n = 4
Q2 = f('1.5'); Q10 = Q11 = Q0; Q0n = Q0; Qy = f(1)
def consts(k, c, r, eta, fix_zeta=False, loglog=True):
    k1, k2, k3 = k; c, r, eta = f(c), f(r), f(eta)
    th = lambda y: f(0) if y >= c + r else (mp.pi if y <= c - r else mp.acos((y - c)/r))
    sg = lambda w: c + r*mp.cos(w)
    q = lambda g, a, b: mp.quad(g, [a, b])
    S = lambda kk: 1 - f(kk)/(2**kk - 2)
    L = mp.log(r/(c - f(1)/2)); s1 = c + (c - f(1)/2)**2/r; dl = 2*c - s1 - f(1)/2
    N4 = 2**(n+4); h2 = f(1)/2
    C2b = q(lambda w: (1 - sg(w)) + c2*(sg(w) - S(n+4)), th(1), th(S(n+4)))/(1 - S(n+4))
    C2b += sum(th(S(4+h)) - th(S(5+h)) for h in range(n))
    C2b += q(lambda w: k3*(S(4) - sg(w)) + (sg(w) - h2), th(S(4)), th(h2))/(S(4) - h2)
    C2b += (c2/eta)*q(lambda w: 1 + eta - sg(w), th(1 + eta), th(1))
    C2b += q(lambda w: c2*(1 - 2*sg(w)) + 2*k3*sg(w), th(h2), th(0))
    C2b += q(lambda w: (sg(w) + eta)/eta*c2, th(0), th(-eta))
    C2 = C2b/(2*mp.pi*L)
    lc = mp.log(c1*f('1.00212')**c2); lY = mp.log(f('1.546')*f('1.00212')); lk = mp.log(k1*f('1.00212')**k3)
    D3 = lc/eta*q(lambda w: 1 + eta - sg(w), th(1+eta), th(1)) + mp.log(mp.zeta(1+eta))/eta*q(lambda w: sg(w) - 1, th(1+eta), th(1))
    D3 += lY/(1 - S(n+4))*q(lambda w: 1 - sg(w), th(1), th(S(n+4))) + lc/(1 - S(n+4))*q(lambda w: sg(w) - S(n+4), th(1), th(S(n+4)))
    D3 += sum(lY*(th(S(4+h)) - th(S(5+h))) for h in range(n))
    D3 += lk/(S(4) - h2)*q(lambda w: S(4) - sg(w), th(S(4)), th(h2)) + lY/(S(4) - h2)*q(lambda w: sg(w) - h2, th(S(4)), th(h2))
    D3 += mp.log(c1*f('1.00212')**c2/mp.sqrt(2*mp.pi))*q(lambda w: 1 - 2*sg(w), th(h2), th(0)) + 2*lk*q(lambda w: sg(w), th(h2), th(0))
    num = (mp.zeta(1 + eta)/f('1.00212')**c2) if fix_zeta else (1 + eta)
    D3 += q(lambda w: -(sg(w)/eta)*mp.log(num/(c1*(2*mp.pi)**eta)) + mp.log(c1*f('1.00212')**c2/mp.sqrt(2*mp.pi)), th(0), th(-eta))
    D3 += -mp.log(2*mp.pi)*q(lambda w: (1 - 2*sg(w))/2, th(-eta), mp.pi)
    lz = (mp.log(mp.zeta(1 + eta)) + mp.log(mp.zeta(c)))/2
    D3 += lz*(th(1 + eta) - mp.pi/2) + mp.pi/(4*J1)*mp.log(mp.zeta(c)) + lz*(th(1 - c) - th(-eta)) + (mp.pi - th(1 - c))/(2*J2)*mp.log(mp.zeta(c))
    k1s = mp.pi/(4*J1)*(mp.log(mp.zeta(c + r)) + 2*mp.fsum(mp.log(mp.zeta(c + r*mp.cos(mp.pi*j/(2*J1)))) for j in range(1, J1)))
    t1c = th(1 - c)
    k2s = (mp.pi - t1c)/(2*J2)*(mp.log(mp.zeta(1 - c + r)) + 2*mp.fsum(mp.log(mp.zeta(1 - c - r*mp.cos(mp.pi*j/J2 + (1 - f(j)/J2)*t1c))) for j in range(1, J2)))
    Ls = lambda w, j: (j + c + r*mp.cos(w))**2/T0 + (r*mp.sin(w))**2/T0 + 2*r*mp.sin(w)
    M1 = q(lambda w: Ls(w, -1), 0, th(1 + eta)) + q(lambda w: Ls(w, Q0), th(1 + eta), th(1))
    M1 += q(lambda w: Ls(w, -1) + (1 - 2*sg(w) + 4*k2*sg(w))/2*Ls(w, Q11), th(h2), th(0))
    M1 += q(lambda w: ((N4 - 1)*(1 - sg(w)) + (N4 - 2)*(sg(w) - S(n+4)))*Ls(w, Q0n), th(1), th(S(n+4)))/((N4 - 2)*(1 - S(n+4)))
    for h in range(n):
        A, B = 2**(h+4) - 2, 2**(h+5) - 2
        M1 += q(lambda w: (B*(A + 1)*(S(5+h) - sg(w)) + A*(B + 1)*(sg(w) - S(4+h)))*Ls(w, Qy), th(S(5+h)), th(S(4+h)))/(A*B*(S(5+h) - S(4+h)))
    M1 += q(lambda w: (14*(k2 + 1)*(S(4) - sg(w)) + 15*(sg(w) - h2))*Ls(w, Q2), th(S(4)), th(h2))/(14*(S(4) - h2))
    M1 += q(lambda w: Ls(w, -1) + (h2 - sg(w))*Ls(w, Q10), th(0), th(-eta)) + q(lambda w: Ls(w, -1) + (1 - 2*sg(w))/2*Ls(w, 1), th(-eta), mp.pi)
    M2 = q(lambda w: c2/eta*(1 + eta - sg(w))*Ls(w, Q0), th(1 + eta), th(1)) + q(lambda w: ((1 - sg(w)) + c2*(sg(w) - S(n+4)))*Ls(w, Q0n), th(1), th(S(n+4)))/(1 - S(n+4))
    M2 += sum(q(lambda w: Ls(w, Qy), th(S(5+h)), th(S(4+h))) for h in range(n))
    M2 += q(lambda w: (k3*(S(4) - sg(w)) + (sg(w) - h2))*Ls(w, Q2), th(S(4)), th(h2))/(S(4) - h2)
    M2 += q(lambda w: (c2*(1 - 2*sg(w)) + 2*k3*sg(w))*Ls(w, Q11), th(h2), th(0)) + q(lambda w: (sg(w) + eta)/eta*c2*Ls(w, Q10), th(0), th(-eta))
    k3s = max(0, M1)/(2*T0) + max(0, M2)/(2*T0*(mp.log(mp.log(T0)) if loglog else mp.log(T0)))
    C3 = f(7)/8 + f(1)/4 + 1/(50*T0) + mp.log(mp.zeta(s1))/mp.pi + mp.log(mp.zeta(c)/mp.zeta(2*c))/(2*L) + ((640*dl - 112)/(1536*(3*T0 - 1)) + f(2)**-10)/2 + (D3 + k1s + k2s + k3s)/(2*mp.pi*L)
    C3T = C3 - f(7)/8 - 1/(50*T0) + (mp.atan((s1 - 1)/T0) + mp.atan(1/(2*T0)))/mp.pi
    arc = q(lambda w: -sg(w)/eta, th(0), th(-eta))
    return C2, C3T, arc, L
print('row  C2(own)      C2(IS raw)   d        C~3(own,TeX loglogT0)  C~3(own,logT0)  C~3(IS raw)  d(logT0)   C~3 with zeta(1+eta)  increase')
for i, (k, c, r, eta, C2p, C3p) in enumerate(ROWS, 1):
    C2, C3a, arc, L = consts(k, c, r, eta, loglog=True)
    _, C3b, _, _ = consts(k, c, r, eta, loglog=False)
    _, C3z, _, _ = consts(k, c, r, eta, fix_zeta=True, loglog=False)
    print('%d  %.8f  %.8f  %+.1e  %.8f  %.8f  %.8f  %+.1e  %.8f  %+.6f' % (i, C2, C2p, C2 - C2p, C3a, C3b, C3p, C3b - C3p, C3z, C3z - C3b))

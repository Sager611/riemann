# Size of the HSW slip (log(1+eta) for log zeta(1+eta) on -eta <= sigma <= 0) for the published parameter rows of
# Hasanalizade-Shen-Wong (Table 2, no 1.00212 factor), Bellotti-Wong (Table 2, with 1.00212^c2, c2 = 1) and Ismail rows.
import mpmath as mp
mp.mp.dps = 30
def delta(c, r, eta, fudge):
    c, r, eta = mp.mpf(c), mp.mpf(r), mp.mpf(eta)
    th = lambda y: mp.acos((y - c)/r)
    arc = mp.quad(lambda w: -(c + r*mp.cos(w))/eta, [th(0), th(-eta)])
    L = mp.log(r/(c - mp.mpf(1)/2))
    gap = mp.log(mp.zeta(1 + eta)) - mp.log(1 + eta) - (mp.log(mp.mpf('1.00212')) if fudge else 0)
    return gap*arc/(2*mp.pi*L)
print('HSW Table 2 (C3 printed):')
for c, r, eta, C3 in [('1.000011314','1.064340602','4.2826451e-6',9.367419),('1.025253504','1.182375395','0.009944751381',4.030486),('1.035766557','1.229059659','0.014325507360',3.746756)]:
    dd = delta(c, r, eta, False); print('  c=%s r=%s eta=%s  dC3 = %.7f  C3 %.6f -> %.6f' % (c, r, eta, dd, C3, C3 + dd))
print('BW Table 2 (C3, C~3 printed):')
for c, r, eta, C3, C3t in [('1.000225','1.000605','0.000158',8.08344,7.20844),('1.070007','1.182997','0.069901',3.71067,2.83567),('1.043400','1.250450','0.040000',3.77417,2.89916),('1.000060','1.499556','1.542440e-5',6.25796,5.38296),('1.499159','1.998357','0.499050',1.96334,1.08834)]:
    dd = delta(c, r, eta, True); print('  c=%s r=%s eta=%s  dC3 = %.7f  C3 %.5f -> %.5f  C~3 %.5f -> %.5f' % (c, r, eta, dd, C3, C3 + dd, C3t, C3t + dd))

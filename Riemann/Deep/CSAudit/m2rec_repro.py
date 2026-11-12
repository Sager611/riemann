# (a) C1(3000) = 20.7225 (TeX l. 371, Cor 2.3) from F2 and J_n at T0' = 3000; (b) the recursion of l. 418-427 for M2(10^n, 2 10^n), M2(10^n), frak M2(10^n).
import mpmath as mp, f1_repro as f
mp.mp.dps = 25
pi, ln = mp.pi, mp.log
# (a) at T0' = 3000, c = (-0.923364, 0.161176), a = (0.531241, 1.52906), b = (0.754804, 1.68921); sigma's are irrelevant for F2 (dummy valid values)
T0 = mp.mpf(3000)
F, ok, P = f.F1(T0, '-0.923364', '0.161176', '0.4', '0.7', '0.4', '0.55', '0.45', '0.531241', '1.52906', '0.754804', '1.68921')
Js = [P['J3'], P['J4'], P['J5'], P['J6']]
F2 = P['J0'] + 2*P['J1'] + mp.fsum(Js) + P['mixed']
C1 = 2*mp.fsum(mp.sqrt(J/pi**2 + 2*P['J1']*J/ln(T0/2)) for J in Js) + F2/mp.sqrt(ln(T0/2))
print('C1(3000) = %.6f (printed 20.7225); conditions ok: %s; J1=%.4f J3..J6=%s J0=%.3g' % (float(C1), ok, float(P['J1']), [round(float(j), 4) for j in Js], float(P['J0'])))
# (b) recursion with the printed frak M1 (Table 2) and with our recomputed F1 (both give upper bounds since ours <= printed)
M1p = {4: '59.300', 5: '48.801', 6: '43.230', 7: '39.720', 8: '37.550', 9: '36.211', 10: '34.842', 15: '31.170', 20: '29.364'}
M1p = {k: mp.mpf(v) for k, v in M1p.items()}
pr_MM2 = {4: '9.9857e8', 5: '1.9276e10', 6: '2.8034e11', 7: '3.9692e12', 8: '5.4791e13', 9: '7.2965e14', 10: '9.5079e15', 15: '3.0154e21', 20: '6.5442e26'}
pr_M2 = {4: '3.8095e8', 5: '8.7748e9', 6: '1.3688e11', 7: '1.8478e12', 8: '2.5002e13', 9: '3.3317e14', 10: '4.3883e15', 15: '1.4174e21', 20: '3.1233e26'}
pr_fM2 = {4: '1.3796e9', 5: '2.8051e10', 6: '4.1722e11', 7: '5.8170e12', 8: '7.9793e13', 9: '1.0629e15', 10: '1.3897e16', 15: '4.4328e21', 20: '9.6675e26'}
def C2(x): return 2/pi**2*x*ln(x)**4 + mp.mpf('20.7225')*(2*x*ln(x)**mp.mpf(3.5)) + mp.mpf('1.25e6')
MM2, M2, fM2 = {}, {}, {}
MM2[4] = C2(mp.mpf(10)**4); M2[4] = C2(mp.mpf(10)**4/2) + mp.mpf('7.032e5'); fM2[4] = M2[4] + MM2[4]
A = {n: list(range(1, n - 3)) for n in range(5, 11)}; A[15] = list(range(5, 12)); A[20] = [5] + list(range(10, 17))
for n in [5, 6, 7, 8, 9, 10, 15, 20]:
    X = mp.mpf(10)**n
    c1 = [X/pi**2*ln(X)**4 + 2*X*M1p[n - i]*ln(X)**3 + MM2[n - i] for i in A[n]]
    c2 = [X/(2*pi**2)*ln(X/2)**4 + M1p[n - i]*X*ln(X/2)**3 + fM2[n - i] for i in A[n]]
    MM2[n] = min([C2(X)] + c1); M2[n] = min([C2(X/2) + mp.mpf('7.032e5')] + c2); fM2[n] = M2[n] + MM2[n]
for n in [4, 5, 6, 7, 8, 9, 10, 15, 20]:
    print('n=%-2d  M2(10^n,2 10^n) <= %.5e (printed %s)  M2(10^n) <= %.5e (printed %s)  frak M2 >= %.5e (printed %s)  %s' % (
        n, float(MM2[n]), pr_MM2[n], float(M2[n]), pr_M2[n], float(fM2[n]), pr_fM2[n],
        'OK' if MM2[n] <= mp.mpf(pr_MM2[n]) and M2[n] <= mp.mpf(pr_M2[n]) and fM2[n] <= mp.mpf(pr_fM2[n]) else 'MISMATCH'))
# Cor 2.4 constants at A0 = 10^11.36: a2 and its pieces
A0 = mp.mpf(10)**mp.mpf('11.36'); L = ln(A0)
print('a2(10^11.36) pieces: main %.5f, M1 term %.5f, M2 term %.5f' % (float((1 - ln(2)/L)**4/(2*pi**2)), float(mp.mpf('34.842')/L*(1 - ln(2)/L)**3), float(mp.mpf('9.5079e15')/(A0*L**4))))

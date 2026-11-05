# Run: uv run --no-project --with mpmath python3 xcheck.py   (output: xcheck_output.txt)

# Independent re-implementation of Proposition 10 (PROOF.md) in mpmath, to cross-check cert.py values.
import mpmath as mp
mp.mp.dps = 50
B = mp.mpf('0.87338196661320752'); b = 1 - B
g = mp.euler; E = mp.mpf('-1.33258227573')
def F(ell, L, h, L1, h1):
    L, h, L1, h1, ell = map(mp.mpf, (L, h, L1, h1, ell))
    x = mp.e**L; x1 = mp.e**L1
    g0 = lambda be: (mp.e**(-(1+h)*L*(1-be))/(1+h) + mp.e**(-L*(1-be)))/(h*L**2)
    g1 = lambda be: (mp.e**(-(1+h1)*L1*(1-be)) + mp.e**(-L1*(1-be)))/(h1*L1)
    k0 = g0(B)*B/b + g0(mp.mpf(1)/2); k1 = g1(B)*B/b + g1(mp.mpf(1)/2)
    kp = k0/(1-k1)
    tau = lambda y: y**(-mp.mpf(1)/2) + ((mp.log(y)/mp.log(2))*(y**(mp.mpf(1)/3)+3)/2 + 8)/y
    T0 = mp.e**(-L)/(L*(1-mp.e**(-2*L))); T1 = mp.e**(-L1)/(1-mp.e**(-2*L1))
    Eh = ((1+h)*mp.log(1+h) - h)/h
    epsD = 1/(10*L**2) + mp.mpf(4)/(15*L**3)
    W = L1*(1+h1/2) + E + mp.mpf('0.2')/L1 + mp.mpf('0.2')/L1**2 + L1*tau(x1) + (1+h1)*L1/(2*(x1-1)) + (ell - mp.log(mp.pi))/2 - g/2 + T1
    C0 = mp.log(mp.pi**2/6) - g
    return C0 - mp.log(L) - Eh - (epsD + 1/(x-1) + tau(x)/2 + T0) - kp*W, k1
# cert values: lower ends of the arb balls printed in cert_output.txt, section (D)
for (Lam, P, cert) in [(3, (12.854, 1.016, 8.0, 1.832), -3.287947688),
                       (5, (19.746, 0.753, 13.297, 1.399), -3.584157302),
                       (40, (287.412, 0.065, 286.0, 0.080), -5.798049278)]:
    v, k1 = F(mp.e**Lam, *P)
    print(Lam, P, 'F =', mp.nstr(v, 10), ' |L(1,chi)| >=', mp.nstr(mp.e**v, 8), ' cert:', cert, ' diff:', mp.nstr(v - cert, 3), ' kappa1 =', mp.nstr(k1, 6))
# tightest cell of (B): [3.13, 3.135], P = (13.133, 1.003, 8.0, 1.879): need F(e^3.135) > log(0.115/3.13)
v, _ = F(mp.e**mp.mpf('3.135'), 13.133, 1.003, 8.0, 1.879)
print('B tight cell margin:', mp.nstr(v - mp.log(mp.mpf('0.115')/mp.mpf('3.13')), 8))

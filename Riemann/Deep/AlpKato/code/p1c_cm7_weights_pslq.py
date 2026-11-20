# P1c: as p1b, e = 1..12, pslq at pi^(2e-j)
# P1b (alpkato): CM forms f_e = theta(psi^{2e}) of weight 2e+1, level 7, nebentypus (-7/.), CM by Q(sqrt-7).
# a_n = (1/2) sum_{N alpha = n} alpha^{2e}; critical values L(f_e, j), j = e+1..2e; algebraic parts vs Omega(49a1)^{2e}.
import mpmath as mp, sympy as sp
mp.mp.dps = 70
M = 900
s7 = mp.sqrt(7)
Om = mp.mpf('1.933311705616811546733077')  # real period of 49a1 (from p1_cm7.log, 25 digits)
# recompute Omega to full precision via AGM as in p1 (one real root)
x = sp.symbols('x'); rts = sp.Poly(4*x**3 - 3*x**2 - 8*x - 4, x).nroots(n=60)
e1 = [r for r in rts if abs(sp.im(r)) < 1e-40][0]; e2 = [r for r in rts if sp.im(r) > 1e-40][0]
e1 = mp.mpf(str(sp.re(e1))); e2 = mp.mpc(str(sp.re(e2)), str(sp.im(e2)))
beta = abs(mp.sqrt((e1-e2)*(e1-mp.conj(e2)))); alpha = 3*e1 - mp.mpf(3)/4
Om = 2*mp.pi/mp.agm(2*mp.sqrt(beta), mp.sqrt(2*beta + alpha))
R = int(2*M**0.5)+3
pts = []
for xx in range(-R, R+1):
    for yy in range(-R, R+1):
        N = xx*xx + xx*yy + 2*yy*yy
        if 0 < N <= M: pts.append((N, sp.Rational(2*xx+yy, 2), sp.Rational(yy, 2)))  # alpha = u + v sqrt(-7)
def coeffs(e):
    a = [0]*(M+1)
    for (N, u, v) in pts:
        # Re((u + v i sqrt7)^(2e)) exactly
        re = sum(sp.binomial(2*e, 2*k) * u**(2*e-2*k) * (-7*v*v)**k for k in range(e+1))
        a[N] += re
    h = [sp.Rational(c)/2 for c in a]
    return [int(c) if c.q == 1 else None for c in h]
out = []
for e in range(1, 13):
    a = coeffs(e); k = 2*e+1
    assert all(c is not None for c in a[1:]), "non-integral"
    # epsilon from theta(1/t) = eps t^k theta(t)
    th = lambda t: mp.fsum(a[n]*mp.e**(-2*mp.pi*n*t/s7) for n in range(1, M+1))
    t = mp.mpf('1.21'); eps = th(1/t)/(t**k*th(t))
    epsr = int(mp.nint(eps))
    def L(s):
        x0 = 2*mp.pi/s7
        lam = mp.fsum(a[n]*(mp.gammainc(s, x0*n)*(s7/(2*mp.pi*n))**s + epsr*mp.gammainc(k-s, x0*n)*(s7/(2*mp.pi*n))**(k-s)) for n in range(1, M+1))
        return lam/((s7/(2*mp.pi))**s*mp.gamma(s))
    print("e=%d k=%d eps=%s (%s) a_2..a_11=%s" % (e, k, epsr, mp.nstr(eps, 12), a[2:12]))
    for j in range(e+1, 2*e+1):
        v = L(j)
        best = None
        pa = 2*e - j
        for b in (0, 1):
            val = v * mp.pi**pa * s7**b / Om**(2*e)
            rel = mp.pslq([val, 1], maxcoeff=10**30, maxsteps=10**6)
            if rel and rel[0] != 0:
                fr = sp.Rational(-rel[1], rel[0])
                if abs(mp.mpf(fr.p)/fr.q - val) < mp.mpf(10)**-50:
                    best = (pa, b, fr); break
        if best:
            pa, b, fr = best
            print("   L(f,%d) * pi^%d * sqrt7^%d / Omega^%d = %s   numerator factors %s  denominator factors %s" % (j, pa, b, 2*e, fr, sp.factorint(abs(sp.numer(fr))), sp.factorint(sp.denom(fr))))
        else:
            print("   L(f,%d) = %s  (no rational identification at pi^%d)" % (j, mp.nstr(v, 20), pa))

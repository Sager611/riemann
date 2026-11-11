# chaudit2: fallbacks and margins of Cor. mertensimple (2511.14736v1) if parts of the square-free chain fail
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi; T = mp.mpf(10)**10 + 1; C0 = mp.mpf('11.3505136')   # 2 x Platt's one-sided sum at 1e10+1
P = print
# (1) c = 1 from Lemma oldtrivial (|R(x)| <= sqrt x for all x >= 0): Cor. metamert for x >= max(e^3 T^2, 4T), Cor. mertensimplon (garance2) below
X = mp.e**3*T**2
C_low = (pi/2 - 3/pi)*mp.mpf(10)**-10*mp.sqrt(X) + mp.mpf('11.350514')      # garance2 written as 3x/(pi 1e10) + C sqrt x, worst at x = X
C_high = C0 + mp.mpf('2.9')*1 + 2/mp.sqrt(X)
P("c = 1 fallback: hand-over e^3 T^2 =", mp.nstr(X, 6), "; C below =", mp.nstr(C_low, 8), "; C above =", mp.nstr(C_high, 8), "=> |M(x)| <= 3x/(pi 1e10) +", mp.nstr(max(C_low, C_high), 6), "sqrt x for all x >= 1")
# (2) MOST-free: c from Prop. qmediano (nopgik's inputs, valid 4.05e6 < x <= 1e32) and gopnik (x >= 1e32)
def Bq(x):   # nopgik's bound / sqrt(x), exact coefficients
    c1, c2 = mp.mpf(16)/25, mp.mpf(114)/25; km = mp.mpf('0.570591'); kap = mp.mpf('6.738093'); eps = pi/(2*mp.mpf(10)**7)
    c3 = (1/(2*pi**2/6) - c1/4)/2
    return (2*eps*mp.sqrt(x) + 5*c1**mp.mpf(0.6)/3*(km/2)**mp.mpf(0.4)*x**mp.mpf(0.4) + 3*c2**(mp.mpf(1)/3)*c3**(mp.mpf(2)/3)*x**(mp.mpf(1)/3)
            + 4*kap/3*x**mp.mpf(0.25) + (c2/c3)**(mp.mpf(1)/6)/2*x**(mp.mpf(1)/6))/mp.sqrt(x)
def need(c):  # sup over omega >= (cT)^2 of the available bound / sqrt(omega)
    w0 = (c*T)**2
    if w0 > mp.mpf(10)**32: return mp.mpf('0.00161')
    s = max(Bq(w0*mp.mpf(10)**(j/50.0)) for j in range(0, int(50*(32 - mp.log10(w0))) + 1))
    return max(s, mp.mpf('0.00161'))
lo, hi = mp.mpf('0.0134'), mp.mpf('0.05')
for _ in range(50):
    mid = (lo + hi)/2
    if need(mid) <= mid: hi = mid
    else: lo = mid
cM = hi
CM = C0 + mp.mpf('2.9')*cM
P("MOST-free: smallest self-consistent c =", mp.nstr(cM, 6), "; sqrt-constant =", mp.nstr(CM, 8), "; hand-over e^3(cT)^2 =", mp.nstr(mp.e**3*(cM*T)**2, 6))
xg = ((CM - mp.mpf('9.7588'))/(pi/(2*(mp.mpf(10)**9 - 1)) - 3/(pi*mp.mpf(10)**10)))**2
P("   garadiol branch (T = 1e9) stays below 3x/(pi 1e10) + C sqrt x up to x =", mp.nstr(xg, 6), "(must exceed the hand-over)")
# (3) where MOST is used: omega in [(0.0134 T)^2, 1e18]
P("MOST range actually used by bataclar's role in mertosimp: omega in [", mp.nstr((mp.mpf('0.0134')*T)**2, 6), ", 1e18 ]; there 1.12543 omega^(-1/4) <=", mp.nstr(mp.mpf('1.12543')*((mp.mpf('0.0134')*T)**2)**mp.mpf(-0.25), 6))
P("qmediano alone gives |R|/sqrt(x) <= 0.0136 from x =", end=" ")
lo, hi = mp.mpf(10)**16, mp.mpf(10)**18
for _ in range(80):
    mid = mp.sqrt(lo*hi)
    if Bq(mid) <= mp.mpf('0.013616'): hi = mid
    else: lo = mid
P(mp.nstr(hi, 6))
# (4) T^2-slack in metamert: (3/pi)/T + A/T^2 + 39/(64T^3) < (3/pi)/(T-1) holds for A up to
for TT in [mp.mpf(50), mp.mpf(10)**10 + 1]:
    Amax = ((3/pi)/(TT - 1) - (3/pi)/TT - mp.mpf(39)/(64*TT**3))*TT**2
    P("metamert T =", mp.nstr(TT, 6), ": allowed coefficient of 1/T^2:", mp.nstr(Amax, 8), "vs used (pi/4)(10/9) =", mp.nstr(pi/4*mp.mpf(10)/9, 8), " slack factor", mp.nstr(Amax/(pi/4*mp.mpf(10)/9), 6))
# (5) MertensCH: CH enters only for x <= X0 ~ 1e240; its sqrt-constant C matters only at x = 1: need 3/(pi 1e10) + C <= smallest record constant 1.54e4
P("MertensCH tolerance: any sqrt-constant C <= 1.54e4 leaves c_CH = max(CH(1), CH(X0)/X0^theta) unchanged; Platt's one-sided sum may be as large as", mp.nstr((mp.mpf('1.54e4') - 2.9*1)/2, 6), "(claimed 5.675256)")
# (6) mertensimplon linear term instead of 3/(pi 1e10): factor
P("if the square-free gain were lost entirely: linear coefficient pi/(2 1e10) instead of 3/(pi 1e10): factor", mp.nstr((pi/2)/(3/pi), 8), "= pi^2/6")

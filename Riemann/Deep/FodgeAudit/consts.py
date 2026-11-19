# FodgeAudit consts.py: Lemma C.3 constants and Prop. C.1 constants recomputed from the stated ingredients in Arb (python-flint), 256 bits.
from flint import arb, acb, ctx
ctx.prec = 256
pi = arb.pi(); z3 = arb(3).zeta(); half = arb(1)/2
def up(x): return x.upper()
def show(name, x, claim=None):
    s = f'{name}: [{x.lower().str(15, radius=False)}, {x.upper().str(15, radius=False)}]'
    if claim is not None:
        ok = x < arb(claim); s += f'  <= {claim}: {ok}  (margin {(arb(claim)-x).upper().str(5, radius=False)})'
        assert ok, name
    print(s)
# --- Lemma C.3 ---
N = arb(10**4); Np = N - half; rL = (arb(1)/75)*N/Np         # L/N' <= (1/75) N/(N-1/2), decreasing in N
hx = ((1 + arb(1)/75)*N/Np)**3; show('h = ((1+1/75) 1e4/(1e4-1/2))^3', hx, '1.0407')
e3 = 3/Np; show('eps3 = 3/N\'', e3, '3.0002e-4'); e4 = 12/Np**2
def c3(h, e3, e4, rL):
    X = e3**2*(h/6 + h**1.5/2)**2 + h.sqrt()*(e3**2*(7*h*h/36 + h**2.5/2) + e4*(h/12 + h**1.5/3))  # Xi
    Np_ = 3/e3; amp = h**2.5*X*(rL*Np_)**2/(4*pi)             # h^(5/2) Xi L^2/(4 pi) with L <= rL N', N' = 3/eps3
    K = 8/pi.sqrt() + (h.sqrt() - 1)/2
    C0 = 1 + e3/8 + amp
    E1 = 7*z3/pi**2*h**1.5 + 21*z3/(2*pi**2)*h
    return K, C0, E1, amp
for lab, h in (('h = 1.0407', arb('1.0407')), ('h = exact bound', hx)):
    K, C0, E1, amp = c3(h, e3, e4, rL)
    print(f'-- {lab}'); show(' 8/sqrt(pi) + (sqrt h - 1)/2', K, '4.524'); show(' amplitude h^(5/2) Xi L^2/(4pi)', amp, '0.00027')
    show(' 1 + eps3/8 + amp', C0, '1.0004'); show(' 7 z3 h^1.5/pi^2 + 21 z3 h/(2 pi^2)', E1, '2.24')
# monotonicity in N of the amplitude: Xi L^2 uses (L/N')^2 and eps4 L^2 = 12 (L/N')^2, eps3^2 L^2 = 9 (L/N')^2: all decrease in N. Check at N = 2e4:
N2 = arb(2*10**4); Np2 = N2 - half
print('   amp at N=2e4 (L = N/75):', c3(((1+arb(1)/75)*N2/Np2)**3, 3/Np2, 12/Np2**2, (arb(1)/75)*N2/Np2)[3].upper().str(8, radius=False), '(decreases in N)')
# mu' identity: mu' = ((N+L)/(N'+L))^(3/2) phi(L/N'), phi(v) = 2(1+v)/(1+sqrt(1+v))
import random; random.seed(1); worst = arb(0)
for _ in range(2000):
    n_ = random.randint(10**4, 10**9); Nn = arb(n_); L = arb(random.randint(2, n_//75))
    Npn = Nn - half
    mu1 = 2*(Nn+L)**1.5/L*(1/Npn.sqrt() - 1/(Npn+L).sqrt()); v = L/Npn
    mu2 = ((Nn+L)/(Npn+L))**1.5*2*(1+v)/(1+(1+v).sqrt())
    d = abs(mu1 - mu2); worst = d if d.upper() > worst.upper() else worst
print('mu\' identity: max |difference| over 2000 random (N, L) <=', worst.upper().str(3, radius=False))

# --- Proposition C.1 constants, independent recomputation (c, r, nu, T1) ---
def propC1(c, r, nu, T1):
    c = arb(c); nu = arb(nu); T1 = arb(T1); r = arb(r)
    t13 = T1**(arb(1)/3); N0 = nu*t13; u = 1/r + 1/N0
    assert N0 > 10**4 and u < arb(1)/75
    N0p = N0 - half; v = u*(N0p + half)/N0p
    mup = (1 + 1/(2*N0p))**1.5*2*(1+v)/(1+(1+v).sqrt())
    K, C0, E1 = arb('4.524'), arb('1.0004'), arb('2.24')
    A = 8*mup/15; ainf = 1/c + A*c.sqrt(); sa = ainf.sqrt()
    abar = lambda y: 1/c + A*(c + 1/y).sqrt()
    ys = pi**(arb(1)/3)*nu; y0 = ys*(1 + 1/r)
    # I_inf = int_{ys}^inf (sqrt(abar(y)) - sqrt(ainf)) dy/y; y = e^s; rigorous acb integration to Y = e^40, tail A/(4 sqrt(c) sa Y)
    Ac, cc, sac = acb(A), acb(c), acb(sa)
    def f(s, analytic):
        return ((1/cc + Ac*(cc + (-s).exp()).sqrt(analytic=analytic)).sqrt(analytic=analytic) - sac)
    s0 = acb(ys.log()); s1 = acb(40)
    Ii = acb.integral(f, s0, s1).real
    Iinf = Ii + A/(4*c.sqrt()*sa*arb(40).exp())
    kap = (2+u)/(2*(1+u).sqrt()); kpp = u/(1+u).log()
    bbar = 8*K/(3*c.sqrt()) + C0/y0 + E1*(c + 1/y0)/(3*y0**3)
    e1 = bbar*pi**(arb(1)/3)*(1+u)*r + c*bbar*pi**(arb(2)/3)*(1+u)**2*r*r/t13
    e0 = c*abar(y0)*pi**(arb(1)/3)*(1+u)*r
    p16 = pi**(arb(1)/6); Cp = ((1+u)/((2*pi).sqrt()*nu)).log()
    b1 = kap*sa/(3*p16); b2 = 4*nu.sqrt() + 2/p16*kap*(sa*Cp + Iinf)
    S = p16*e1*(1+u)**2/(sa*(2*pi).sqrt()) + e0*kpp*(T1.log()/6 + Cp)/(p16*sa*T1**(arb(1)/6))
    assert T1.log()/6 + Cp > 1 and 1.48/T1**0.25 + 0.127/T1**0.75 < 2
    return dict(u=u, mup=mup, kap=kap, ainf=ainf, Iinf=Iinf, bbar=bbar, e1=e1, e0=e0, Cp=Cp, b1=b1, b2=b2, S=S)
for (c_, r_, nu_, T1_, claims) in (('2.44', 80, '0.95', '2.9e12', ('0.3080532', '2.405897', '506.6')),
                                   ('2.44', 80, '1.15', arb('27.4').exp(), ('0.3080571', '2.428866', '485.31'))):
    d = propC1(c_, r_, nu_, T1_)
    print(f'-- Prop C.1 constants c={c_} r={r_} nu={nu_} T1={arb(T1_).str(8, radius=False)}')
    for k in ('u', 'mup', 'kap', 'ainf', 'Iinf', 'bbar', 'e1', 'e0', 'Cp'): show(' '+k, d[k])
    show(' b1', d['b1'], claims[0]); show(' b2', d['b2'], claims[1]); show(' S', d['S'], claims[2])

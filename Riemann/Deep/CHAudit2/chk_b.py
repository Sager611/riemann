# chaudit2: square-free part of 2511.14736v1 (Props. gould l.4356, qmediano l.4538; Lemmas nopgik l.4589, gopnik l.4607;
# Cor. bataclar l.4634; Cor. metamert l.4744, mertosimp l.4791) -- recomputation and falsification.
import math, random, sys
import numpy as np
import mpmath as mp
mp.mp.dps = 40
Z2 = mp.pi**2/6
out = []
def P(*a):
    s = " ".join(str(t) for t in a); print(s); sys.stdout.flush()

# ---------------- Part 1: constants of nopgik / gopnik / bataclar / mertosimp ----------------
def c3f(c1): return (1/(2*Z2) - mp.mpf(c1)/4)/2
def gould_a(eps, kap, c1, c2):
    c3 = c3f(c1)
    return dict(x35=4*eps/(16*kap/c1)**mp.mpf(0.4), x25=5*c1**mp.mpf(0.6)/3*(kap/2)**mp.mpf(0.4),
                x13=3*c2**(mp.mpf(1)/3)*c3**(mp.mpf(2)/3), x14=5*kap/7, x16=(c2/c3)**(mp.mpf(1)/6)/2,
                cond1=64*10**5*(kap/c1)**4, cond2=(c2/c3)**5/(16*kap/c1)**6, cond_omit=mp.mpf(2)**11*73**5*(kap/c1)**4)
def gould_b(eps, kap, c1, c2):
    c3 = c3f(c1)
    return dict(x12=2*mp.sqrt(eps*c1), x13=3*c2**(mp.mpf(1)/3)*c3**(mp.mpf(2)/3), x14=kap/3*(c1/eps)**mp.mpf(0.75),
                x16=(c2/c3)**(mp.mpf(1)/6)/2, cond=(c2/c3)**2*(c1/(16*eps))**3, epsmax=c1/584)
def qmed(eps, kap, km, v0, v1, c1, c2):
    c3 = c3f(c1)
    return dict(x12=2*eps, x25=5*c1**mp.mpf(0.6)/3*(km/2)**mp.mpf(0.4), x13=3*c2**(mp.mpf(1)/3)*c3**(mp.mpf(2)/3),
                x14=4*kap/3, x16=(c2/c3)**(mp.mpf(1)/6)/2,
                conds=[64*10**5*(km/c1)**4, (c2/c3)**5/(16*km/c1)**6, (2*v0**2)**mp.mpf(1.25)/(16*km/c1), 2*v0**2], xmax=v1**2)
f = lambda v: mp.nstr(v, 10)
P("== nopgik (qmediano: kappa_-=0.570591, v0=33, v1=1e16, eps=pi/2e7, kappa=6.738093, c1=16/25, c2=114/25) ==")
q = qmed(mp.pi/(2*10**7), mp.mpf('6.738093'), mp.mpf('0.570591'), mp.mpf(33), mp.mpf(10)**16, mp.mpf(16)/25, mp.mpf(114)/25)
P("coefficients x^1/2, x^2/5, x^1/3, x^1/4, x^1/6:", f(q['x12']), f(q['x25']), f(q['x13']), f(q['x14']), f(q['x16']))
P("paper: pi/1e7 =", f(mp.pi/10**7), "0.772103 0.860845 8.98413 0.999")
P("lower conditions (must be < 1e18):", [mp.nstr(c, 6) for c in q['conds']], " x <= v1^2 =", mp.nstr(q['xmax'], 4))
x0, x1 = mp.mpf(10)**18, mp.mpf(10)**32
def ratio25(x):
    return (q['x12']*x**mp.mpf(0.1) + q['x25'] + q['x13']*x**(mp.mpf(1)/3-mp.mpf(2)/5) + q['x14']*x**(mp.mpf(1)/4-mp.mpf(2)/5) + q['x16']*x**(mp.mpf(1)/6-mp.mpf(2)/5))
bound_paper = q['x12']*x1**mp.mpf(0.1) + q['x25'] + q['x13']*x0**(-mp.mpf(1)/15) + q['x14']*x0**(-mp.mpf(3)/20) + q['x16']*x0**(-mp.mpf(7)/30)
P("paper's endpoint bound for |R|/x^(2/5) =", f(bound_paper), "(<= 0.84491 claimed); true max over [1e18,1e32] =", f(max(ratio25(mp.mpf(10)**(18+14*j/2000)) for j in range(2001))))
P("0.845 x^(2/5) <= 0.0134 x^(1/2) iff x >= (0.845/0.0134)^10 =", mp.nstr((mp.mpf('0.845')/mp.mpf('0.0134'))**10, 8))
cmin_nop = max(ratio25(mp.mpf(10)**(18+14*j/2000))*mp.mpf(10)**((18+14*j/2000)*(-0.1)) for j in range(2001))
P("sup over [1e18,1e32] of the nopgik bound / sqrt(x) (exact coefficients) =", f(cmin_nop), "(this is the c that nopgik really gives)")

P("== gopnik (gould: eps=pi/2e10, kappa=11.350514, c1=768/1225, c2=9458/1225) ==")
eps, kap, c1, c2 = mp.pi/(2*10**10), mp.mpf('11.350514'), mp.mpf(768)/1225, mp.mpf(9458)/1225
ga = gould_a(eps, kap, c1, c2); gb = gould_b(eps, kap, c1, c2)
P("(a) coeffs x^3/5, x^2/5, x^1/3, x^1/6:", f(ga['x35']), f(ga['x25']), f(ga['x13']), f(ga['x16']), " paper 6.51e-11 2.5222 1.042 1.086")
P("(a) conditions:", mp.nstr(ga['cond1'], 6), mp.nstr(ga['cond2'], 6), " omit-5k/7 condition:", mp.nstr(ga['cond_omit'], 6), "(all < 1e32)")
P("(a) at 1e32: 2.5222 + 1.042 x^(-1/15) + 1.086 x^(-7/30) =", f(ga['x25'] + ga['x13']*mp.mpf(10)**(-mp.mpf(32)/15) + ga['x16']*mp.mpf(10)**(-mp.mpf(32)*7/30)), "(<= 2.53 claimed)")
P("(a) |R|/sqrt x on [1e32,1e50] <=", f(ga['x35']*mp.mpf(10)**5 + mp.mpf('2.53')*mp.mpf(10)**(-mp.mpf('3.2'))), "(<= 0.00161 claimed); exact-coefficient sup:",
  f(max((ga['x35']*mp.mpf(10)**((32+18*j/1000)*0.1) + ga['x25']*mp.mpf(10)**(-(32+18*j/1000)*0.1) + ga['x13']*mp.mpf(10)**(-(32+18*j/1000)/6) + ga['x16']*mp.mpf(10)**(-(32+18*j/1000)/3)) for j in range(1001))))
P("(b) coeffs x^1/2, x^1/3, x^1/4, x^1/6:", f(gb['x12']), f(gb['x13']), f(gb['x14']), f(gb['x16']), " paper 1.99e-5 1.042 6.01e7 1.1")
P("(b) eps < c1/584:", eps < gb['epsmax'], " condition x >=", mp.nstr(gb['cond'], 6), "(< 1e50)")
xb = mp.mpf(10)**50
P("(b) at 1e50: |R|/sqrt x <=", f(mp.mpf('1.99e-5') + mp.mpf('1.042')*xb**(-mp.mpf(1)/6) + mp.mpf('6.01e7')*xb**(-mp.mpf(1)/4) + mp.mpf('1.1')*xb**(-mp.mpf(1)/3)), "(<= 3.9e-5 claimed)")
P("== bataclar: jbond 1.12543 x^(1/4) <= 0.0134 sqrt(x) iff x >=", mp.nstr((mp.mpf('1.12543')/mp.mpf('0.0134'))**4, 8), "(<= 1e8 needed)")
P("== mertosimp ==")
T = mp.mpf(10)**10 + 1
P("2.9*0.0134 =", f(mp.mpf('2.9')*mp.mpf('0.0134')), "; 11.350514 + 0.0389 =", f(mp.mpf('11.350514') + mp.mpf('0.0389')))
P("absorb 1/|zeta(0)| = 2: 11.389414 sqrt(x) + 2 <= 11.39 sqrt(x) iff x >=", mp.nstr((2/(mp.mpf('11.39') - mp.mpf('11.389414')))**2, 6))
P("e^3 (cT)^2 with c = 0.0134:", mp.nstr(mp.e**3*(mp.mpf('0.0134')*T)**2, 8), "(< 3.61e17 needed); largest c with e^3(cT)^2 <= 3.61e17:", f(mp.sqrt(mp.mpf('3.61e17')/mp.e**3)/T))
P("largest c keeping 11.39: (11.39 - 11.350514)/2.9 =", f((mp.mpf('11.39') - mp.mpf('11.350514'))/mp.mpf('2.9')))
P("garadiol branch below 11.39 bound up to x =", mp.nstr(((mp.mpf('11.39') - mp.mpf('9.7588'))/(mp.pi/(2*(mp.mpf(10)**9 - 1)) - 3/(mp.pi*mp.mpf(10)**10)))**2, 8))
P("metamert: 3/(pi T) + (pi/4)(10/9)/T^2 + 39/(64T^3) < 3/(pi(T-1)) at T=50:", 3/(mp.pi*50) + mp.pi/4*mp.mpf(10)/9/2500 + mp.mpf(39)/(64*125000) < 3/(mp.pi*49),
  "; trivial-zero constant (1/(2-1) + 2*pi/(12 pi))(2pi)^2/zeta(3) =", f((1 + mp.mpf(1)/6)*(2*mp.pi)**2/mp.zeta(3)), "< 39")

# ---------------- Part 2: proof-as-computation of gould (a), (b) and qmediano ----------------
z2 = math.pi**2/6
def bounds_from_proof(x, eps, kap, c1, c2, case, km=None, v0=None):
    c3 = 0.5*(1/(2*z2) - c1/4)
    alpha, beta = c2, 2*c3*math.sqrt(x)
    ups0 = (beta/(2*alpha))**(2/3)
    Kp = math.ceil(ups0 - 0.5)
    if case in ('a', 'q'):
        kk = km if case == 'q' else kap
        a = 8*kk/3*x**0.25; b = c1/4*math.sqrt(x)
        y0 = (2*b/(3*a))**0.8
        K = math.ceil(y0 - 0.5)
    else:
        a = 4*eps; b = c1/4
        y0 = math.sqrt(b/a)
        K = math.ceil(y0**2 - 0.5)
    ok = (K >= 3) and (Kp >= K)
    ks = np.arange(1, K+1, dtype=float)
    U = K + 0.5
    if case == 'q':
        # main part with |M(v)| <= km sqrt v on [v0, v1] for u >= 1, eps v + kap sqrt v for u < 1
        main = np.sum(km*(x/ks)**0.25) + km*x**0.25*(4/3)*(U**0.75 - 1) + 2*eps*math.sqrt(x) + 4/3*kap*x**0.25
        ok = ok and (math.sqrt(x/U) >= v0)
    else:
        main = np.sum(eps*np.sqrt(x/ks) + kap*(x/ks)**0.25) + 2*math.sqrt(U)*eps*math.sqrt(x) + 4/3*U**0.75*kap*x**0.25
    err = c1*math.sqrt(x)/4*(1/math.sqrt(K+0.5) - 1/math.sqrt(Kp+0.5)) + c2*(Kp-K) + 0.5*(math.sqrt(x/(Kp+0.5))/z2 + (x/(Kp+0.5))**0.25)
    return main + err, ok, K, Kp
def closed(x, eps, kap, c1, c2, case, km=None):
    c3 = 0.5*(1/(2*z2) - c1/4)
    if case == 'a':
        r = 4*eps/(16*kap/c1)**0.4*x**0.6 + 5*c1**0.6/3*(kap/2)**0.4*x**0.4 + 3*c2**(1/3)*c3**(2/3)*x**(1/3) + 0.5*(c2/c3)**(1/6)*x**(1/6)
        if not x > 2**11*73**5*(kap/c1)**4: r += 5*kap/7*x**0.25
        return r
    if case == 'b':
        return 2*math.sqrt(eps*c1)*x**0.5 + 3*c2**(1/3)*c3**(2/3)*x**(1/3) + kap/3*(c1/eps)**0.75*x**0.25 + 0.5*(c2/c3)**(1/6)*x**(1/6)
    return 2*eps*x**0.5 + 5*c1**0.6/3*(km/2)**0.4*x**0.4 + 3*c2**(1/3)*c3**(2/3)*x**(1/3) + 4*kap/3*x**0.25 + 0.5*(c2/c3)**(1/6)*x**(1/6)
pairs = [(1, 1), (3/4, 3/2), (2/3, 8/3), (16/25, 114/25), (768/1225, 9458/1225), (18432/29645, 361192/29645), (442368/715715, 14328304/715715)]
random.seed(7)
worst = {'a': (0,), 'b': (0,), 'q': (0,)}; fails = {'a': 0, 'b': 0, 'q': 0}; n_case = {'a': 0, 'b': 0, 'q': 0}
for it in range(30000):
    case = random.choice('abq')
    if random.random() < 0.5:
        c1, c2 = random.choice(pairs)
    else:
        c1 = random.uniform(0.01, 1.0); c2 = 10**random.uniform(-3, 1.5)
    c3 = 0.5*(1/(2*z2) - c1/4)
    kap = 10**random.uniform(-3, 2)
    if case == 'a':
        eps = random.choice([0.0, 10**random.uniform(-12, -1)])
        lo = max(64e5*(kap/c1)**4, (c2/c3)**5/(16*kap/c1)**6)
        x = lo*(1 + 1e-9) if random.random() < 0.2 else lo*10**random.uniform(0, 12)
        if x > 1e40: continue
        B, ok, K, Kp = bounds_from_proof(x, eps, kap, c1, c2, 'a')
    elif case == 'b':
        eps = c1/584*random.choice([0.999999, random.uniform(0, 1)])
        if eps <= 0: continue
        lo = (c2/c3)**2*(c1/(16*eps))**3
        x = lo*(1 + 1e-9) if random.random() < 0.2 else lo*10**random.uniform(0, 12)
        if x > 1e40 or x < 1: continue
        B, ok, K, Kp = bounds_from_proof(x, eps, kap, c1, c2, 'b')
    else:
        km = 10**random.uniform(-2, 1); v0 = 10**random.uniform(0, 3); eps = random.choice([0.0, 10**random.uniform(-12, -1)])
        lo = max(64e5*(km/c1)**4, (c2/c3)**5/(16*km/c1)**6, (2*v0**2)**1.25/(16*km/c1), 2*v0**2)
        x = lo*(1 + 1e-9) if random.random() < 0.2 else lo*10**random.uniform(0, 12)
        if x > 1e40: continue
        B, ok, K, Kp = bounds_from_proof(x, eps, kap, c1, c2, 'q', km=km, v0=v0)
    C = closed(x, eps, kap, c1, c2, case, km=km if case == 'q' else None)
    n_case[case] += 1
    if not ok: fails[case] += 1; continue
    r = B/C
    if r > worst[case][0]: worst[case] = (r, x, eps, kap, c1, c2, K, Kp)
P("== Part 2: proof-as-computation (sum form of Prop. andalas + Lemma oldtrivial with the proof's K, K'; exact sums, no Lemma ramando) vs closed forms ==")
for cs_ in 'abq':
    P("case", cs_, ": cases", n_case[cs_], " structural failures (K<3/37, K'<K or v0 check):", fails[cs_], " max (proof bound)/(closed form) =", worst[cs_][0], "at", worst[cs_][1:])

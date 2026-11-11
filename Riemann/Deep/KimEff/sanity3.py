
# sanity3: numerical checks of the Bessel facts (G1)-(G6) of DERIVATION Sec. 5 for moderate b (not part of the proof)
import mpmath as mp, random
mp.mp.dps = 30
random.seed(1)
bad = 0
for nn in [5, 20, 60]:
    b = nn + mp.mpf(1)/2; nu = b - 1
    j1 = mp.besseljzero(nu, 1)
    c1 = j1**2/4
    print('n=%d nu=%s j_{nu,1}=%s > nu: %s ; c1/(nu^2/4)=%s' % (nn, nu, mp.nstr(j1,10), j1 > nu, mp.nstr(c1/(nu**2/4),6)))
    # Rayleigh sums from zeros (first 400 zeros + tail estimate) vs 1/b and 1/(b^2(b+1))
    ck = [mp.besseljzero(nu, k)**2/4 for k in range(1, 401)]
    S1 = sum(1/c for c in ck); S2 = sum(1/c**2 for c in ck)
    print('   S1(400 zeros)=%s vs 1/b=%s ; S2=%s vs 1/(b^2(b+1))=%s' % (mp.nstr(S1,8), mp.nstr(1/b,8), mp.nstr(S2,10), mp.nstr(1/(b**2*(b+1)),10)))
    G = lambda x: mp.hyp0f1(b, x)
    S1e = 1/b; S2e = 1/(b**2*(b+1))
    for trial in range(40):
        tau = mp.mpf(random.uniform(-0.5, 0.97)); X = tau*nu**2/4  # X <= 0.97 nu^2/4 <= 0.97 c1
        Y = mp.mpf(random.uniform(-1, 1))*mp.sqrt(abs(X)+1)*3
        x = mp.mpc(-X, Y)
        t_true = X/c1; tp = max(t_true, 0)
        lhs = abs(G(x))
        ok1 = lhs >= G(-X)*(1-mp.mpf(10)**-20)
        ok2 = mp.log(lhs) <= mp.log(G(-X)) + Y**2*S2e/(2*(1-tp)**2) + mp.mpf(10)**-20
        lG = mp.log(G(-X))
        if X >= 0:
            ok3 = (-X*S1e - X**2*S2e/(2*(1-t_true)) <= lG + 1e-20) and (lG <= -X*S1e - X**2*S2e/2 + 1e-20)
            g = mp.diff(lambda y: mp.log(G(y)), -X)
            ok4 = S1e - 1e-20 <= g <= S1e + X*S2e/(1-t_true) + 1e-20
        else:
            ok3 = lG >= 0
            g = mp.diff(lambda y: mp.log(G(y)), -X)
            ok4 = 0 < g <= S1e + 1e-20
        # phase derivative bound: d/dt arg G(t^2 w/4) at t=1 with w = 4x  ->  (|v|/2) * sum c_k/((c_k-X)^2+Y^2) <= (|v|/2)(S1 + max(X,0) S2 (2-tp)/(1-tp)^2)
        w = 4*x; v = mp.im(w)
        dphase = mp.im(mp.diff(lambda t: mp.log(G(t**2*w/4)), 1))
        ok5 = abs(dphase) <= (abs(v)/2)*(S1e + max(X,0)*S2e*(2-tp)/(1-tp)**2) + 1e-15
        if not (ok1 and ok2 and ok3 and ok4 and ok5):
            bad += 1; print('   FAIL', nn, mp.nstr(tau,4), mp.nstr(Y,4), ok1, ok2, ok3, ok4, ok5)
print('sanity3 failures:', bad)

"""liuaudit2 T9: Prop 5.2 (l.2812-2866) local Mellin factor and multiplier growth, r_F = 2 (Q(zeta12)).
(a) int_0^oo x K_{1/3}(c sqrt x) x^t dx/x = 2^{2t+1} c^{-2t-2} Gamma(5/6+t) Gamma(7/6+t)   (Re t > -5/6)
(b) |[Gamma(5/6+t)Gamma(7/6+t)/(Gamma(5/6-t)Gamma(7/6-t))]^{r_F}| ~ C |tau|^{4 r_F sigma}
(c) the multiplier is holomorphic on Re t > -5/6: numerator poles at t = -5/6 - k, -7/6 - k only.
Run: uv run --no-project --with mpmath python3 mellin_q12.py
"""
import mpmath as mp
mp.mp.dps = 30
fails = 0
for c in (mp.mpf(1), mp.mpf('2.5'), mp.mpf('0.7')):
    for t in (mp.mpf('-0.6'), mp.mpf('-0.33'), mp.mpf(0), mp.mpf('0.5'), mp.mpc('0.25', '3.0'), mp.mpc('-0.5', '-1.5')):
        lhs = mp.quad(lambda x: x*mp.besselk(mp.mpf(1)/3, c*mp.sqrt(x))*x**t/x, [0, mp.mpf('1e-12'), mp.mpf('1e-8'), mp.mpf('1e-4'), mp.mpf('1e-2'), 1, 10, 100, mp.inf], maxdegree=10)
        rhs = 2**(2*t+1)*c**(-2*t-2)*mp.gamma(mp.mpf(5)/6+t)*mp.gamma(mp.mpf(7)/6+t)
        rel = abs(lhs-rhs)/abs(rhs)
        ok = rel < mp.mpf('1e-9')  # quadrature limit near the endpoint singularity x^(t-1/6)
        fails += (not ok)
        print("(a) c=%s t=%s rel.err=%.2e %s" % (mp.nstr(c, 3), mp.nstr(t, 4), float(rel), "ok" if ok else "FAIL"))
rF = 2
def ratio(t): return (mp.gamma(mp.mpf(5)/6+t)*mp.gamma(mp.mpf(7)/6+t)/(mp.gamma(mp.mpf(5)/6-t)*mp.gamma(mp.mpf(7)/6-t)))**rF
for sigma in (mp.mpf('-0.3'), mp.mpf('0.25'), mp.mpf(1), mp.mpf(3)):
    vals = [abs(ratio(mp.mpc(sigma, tau)))/mp.mpf(tau)**(4*rF*sigma) for tau in (10**3, 10**4, 10**5)]
    drift = abs(mp.log(vals[2]/vals[0]))
    ok = drift < mp.mpf('1e-2')
    fails += (not ok)
    print("(b) sigma=%s |m|/tau^(4 r_F sigma) at tau=1e3,1e4,1e5: %s %s" % (mp.nstr(sigma, 3), [mp.nstr(v, 8) for v in vals], "ok" if ok else "FAIL"))
# (c) no pole of the ratio in -5/6 < Re t: check finiteness on a grid near the strip boundary
mx = 0
for s_ in [mp.mpf(-5)/6 + mp.mpf('1e-3'), mp.mpf('-0.5'), mp.mpf('-1')/3, mp.mpf(0)]:
    for tau in [0, mp.mpf('0.1'), 1, 5]:
        v = abs(ratio(mp.mpc(s_, tau)))
        mx = max(mx, v)
print("(c) max |ratio| on grid in -5/6 < Re t <= 0: %s (finite)" % mp.nstr(mx, 6))
print("RESULT:", "PASS" if fails == 0 else "FAIL (%d)" % fails)

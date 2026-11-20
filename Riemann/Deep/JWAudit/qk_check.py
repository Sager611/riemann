# JWAudit independent algebra checks.
# (1) Q_k (k = 2..6): derived directly in u (Psi' = 2 log u0, d/dz = Z'(u)^{-1} d/du), compared with the (t,s)
#     construction of certify6.py l. 106-116 and with numerical derivatives of Psi(x) = f(u0(x), x) at x = 1e8.
# (2) Savva's T_k recurrence (13.6) with P1's factor (m - d + 1), exact rationals; (m - d) must fail.
# (3) Jacobi diagonal identity (SD l. 87): (U/2)(1 - b_k) - V = [V S (d-1) + 2k(S+k+1)(U-2V)]/((2k+S)(2k+S+2)).
# Run: uv run --no-project --with sympy --with mpmath python3 qk_check.py
import sympy as sp, mpmath as mp
u = sp.symbols('u', positive=True)
E = sp.exp(2 * u)
Zu = u * (sp.pi * E - sp.Rational(9, 4)); Zpu = sp.diff(Zu, u)
P = {1: 2 * sp.log(u)}
for k in range(2, 7): P[k] = sp.diff(P[k - 1], u) / Zpu
Qm = {k: P[k] * Zu ** (k - 1) * 2 * u for k in range(2, 7)}
t, s = sp.symbols('t s', positive=True)
Zp = sp.Rational(9, 4) * (t + 2 - t * s) / (t * s); Zs = sp.Rational(9, 4) * (1 - s) / (t * s)
def Dz(F): return sp.cancel((sp.diff(F, t) * (-t ** 2) + sp.diff(F, s) * (-2 * s)) / Zp)
Ps = {2: sp.cancel(2 * t / Zp)}
for k in range(2, 6): Ps[k + 1] = Dz(Ps[k])
Qs = {k: sp.cancel(Ps[k] * Zs ** (k - 1) * 2 / t) for k in range(2, 7)}
mp.mp.dps = 50
ok = True
for k in range(2, 7):
    fm = sp.lambdify(u, Qm[k], 'mpmath'); fs = sp.lambdify((t, s), Qs[k], 'mpmath')
    worst = max(abs(fm(mp.mpf(uv)) / fs(1 / mp.mpf(uv), 9 / (4 * mp.pi * mp.exp(2 * mp.mpf(uv)))) - 1) for uv in ('0.7', '1.9', '3.3', '7.1', '45'))
    ser = sp.series(Qs[k].subs(s, 0), t, 0, 3).removeO()
    print("Q_%d: u-derivation vs (t,s) code: max rel diff %s; Q_%d(t,0) = %s + O(t^3)" % (k, mp.nstr(worst, 3), k, sp.expand(ser)))
    ok &= worst < mp.mpf(10) ** -40
# numerical derivatives of Psi(x) = f(u0(x), x) at x = 1e8 (independent of all algebra)
mp.mp.dps = 120
x0 = mp.mpf(10) ** 8
def u0of(x):
    v = mp.log(x) / 2
    for _ in range(300):
        dv = (v * (mp.pi * mp.exp(2 * v) - mp.mpf(9) / 4) - x) / (mp.pi * mp.exp(2 * v) * (1 + 2 * v) - mp.mpf(9) / 4); v -= dv
        if abs(dv) < mp.mpf(10) ** -110: break
    return v
Psi = lambda x: 2 * x * mp.log(u0of(x)) + mp.mpf(9) * u0of(x) / 2 - mp.pi * mp.exp(2 * u0of(x))
uu = u0of(x0)
print("Psi'(x0) numeric %s vs 2 log u0 %s" % (mp.nstr(mp.diff(Psi, x0, 1), 20), mp.nstr(2 * mp.log(uu), 20)))
for k in range(2, 7):
    num = mp.diff(Psi, x0, k)
    fs = sp.lambdify((t, s), Qs[k], 'mpmath')
    pred = fs(1 / uu, 9 / (4 * mp.pi * mp.exp(2 * uu))) / (x0 ** (k - 1) * 2 * uu)
    print("Psi^(%d)(1e8): numeric/formula - 1 = %s" % (k, mp.nstr(num / pred - 1, 3)))
    ok &= abs(num / pred - 1) < mp.mpf(10) ** -30
# (2) recurrence, exact
d = 6; A, B, C, D = sp.Rational(37, 2), sp.Rational(23, 3), sp.Rational(29, 4), sp.Rational(13, 5); y0 = sp.Rational(3, 7)
Y = sp.symbols('Y')
kap = D * Y / (A * C)
p = sum(sp.rf(-d, j) * sp.rf(A, j) * sp.rf(C, j) / (sp.rf(B, j) * sp.rf(D, j) * sp.factorial(j)) * kap ** j for j in range(d + 1))
T = [sp.nsimplify((Y ** k * sp.diff(p, Y, k) / p).subs(Y, y0)) for k in range(d + 2)]
kv = kap.subs(Y, y0)
for shift, lab in ((1, "(m-d+1)"), (0, "(m-d) [App l. 688 misprint]")):
    res = []
    for m in range(d - 1):
        P3 = 1 - kv; P2 = B + D + 2 * m + 1 - kv * (A + C + 3 * m - d + 3)
        P1 = (B + m) * (D + m) - kv * ((A + m) * (C + m) + (m - d + shift) * (A + C + 2 * m + 1)); P0 = kv * (d - m) * (A + m) * (C + m)
        res.append(sp.simplify(P3 * T[m + 3] + P2 * T[m + 2] + P1 * T[m + 1] + P0 * T[m]))
    print("recurrence with %s: residuals m=0..%d: %s" % (lab, d - 2, res))
    if shift == 1: ok &= all(r == 0 for r in res)
# (3) Jacobi diagonal identity
U, V, dd, k = sp.symbols('U V d k')
al = V - 1; be = U - V - dd; S = al + be
bk = (be ** 2 - al ** 2) / ((2 * k + S) * (2 * k + S + 2))
lhs = U * (1 - bk) / 2 - V; rhs = (V * S * (dd - 1) + 2 * k * (S + k + 1) * (U - 2 * V)) / ((2 * k + S) * (2 * k + S + 2))
z3 = sp.simplify(lhs - rhs)
print("Jacobi diagonal identity: lhs - rhs = %s" % z3)
ok &= z3 == 0
print("ALL CHECKS:", "PASS" if ok else "FAIL")

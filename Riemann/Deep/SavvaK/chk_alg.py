# Audit checks of Savva's algebraic claims (exact/sympy + mpmath). Not part of the proof.
import sympy as sp, random, mpmath as mp
random.seed(1)
# (2) finite-free factorization p_F = Q^{(1)}_{A,B} [x] Q^{(D)}_{C,D} (ascending, constant term one): coefficientwise
ok2 = True
for trial in range(5):
    d = random.randint(3, 9)
    A, B, C, D = [sp.Rational(random.randint(50, 400), random.randint(2, 7)) for _ in range(4)]
    for j in range(d + 1):
        bn = sp.binomial(d, j)
        q1 = sp.rf(-d, j) * sp.rf(A, j) / (sp.rf(B, j) * sp.factorial(j)) / A ** j
        q2 = sp.rf(-d, j) * sp.rf(C, j) / (sp.rf(D, j) * sp.factorial(j)) * (D / C) ** j
        pf = sp.rf(-d, j) * sp.rf(A, j) * sp.rf(C, j) / (sp.rf(B, j) * sp.rf(D, j) * sp.factorial(j)) * (D / (A * C)) ** j
        # multiplicative finite free convolution in the normalized-coefficient sense: c_j/binom(d,j) multiply (with sign (-1)^j)
        if sp.simplify(pf / bn * (-1) ** j - (q1 / bn * (-1) ** j) * (q2 / bn * (-1) ** j)) != 0: ok2 = False
print("(2) p_F = Q_{A,B}^{(1)} box Q_{C,D}^{(D)} coefficientwise:", "OK" if ok2 else "FAIL")

# (3) limiting system (11.1)-(11.4): zero, Jacobian, det, inverse row sums
al, t, w, de = sp.symbols('alpha t w delta')
F = sp.Matrix([1 / al + w / t ** 2 + de - 2, w / t ** 3 + de - 1, 3 * w / t ** 4 + 3 * de - 2, 4 * w / t ** 5 + 4 * de - 2])
ys = {al: 3, t: 2, w: sp.Rational(16, 3), de: sp.Rational(1, 3)}
J = F.jacobian([al, t, w, de]).subs(ys)
Ji = J.inv()
print("(3) F(y*) =", list(F.subs(ys)), " DF =", J.tolist(), " det =", J.det(), " inverse row sums =", [sum(abs(x) for x in Ji.row(i)) for i in range(4)])
# uniqueness in the positive orthant: 3F1-F2 and 4/3 F2 - F3 identities
print("    (11.2) with 0-based components: 3F_1 - F_2 - (3w(t-1)/t^4 - 1) =", sp.simplify(3 * F[1] - F[2] - (3 * w * (t - 1) / t ** 4 - 1)),
      ";  (4/3)F_2 - F_3 - (4w(t-1)/t^5 - 2/3) =", sp.simplify(sp.Rational(4, 3) * F[2] - F[3] - (4 * w * (t - 1) / t ** 5 - sp.Rational(2, 3))))
print("    with 1-based components 3F_1 - F_2 would be", sp.simplify(3 * F[0] - F[1]), "(not the printed identity)")
# (4) leading E6 at the limit: 48 - 120 (delta + w/t^6)
print("(4) leading n^5 L E6 at y* =", 48 - 120 * (sp.Rational(1, 3) + sp.Rational(16, 3) / 64))
# (5) Jacobi localization (Savva app:jacobi-eight): roots of 2F1(-d,U;V;Y/U) in [V-8 sqrt(Vd), V+8 sqrt(Vd)] for U>=V+d, V>=32d
mp.mp.dps = 60
worst = 0; worstsharp = 0
for (d, V, U) in [(5, 160, 165), (5, 160, 320), (10, 320, 640), (10, 320, 330), (20, 640, 1280), (20, 2000, 4000), (8, 256, 10 ** 6), (30, 960, 1920), (12, 5000, 10000)]:
    cs = [mp.mpf(1)]
    for k in range(d):
        cs.append(cs[-1] * (k - d) * (k + U) / ((k + 1) * (k + V) * U))
    roots = mp.polyroots(cs[::-1], maxsteps=400, extraprec=400)
    dev = max(abs(mp.re(r) - V) for r in roots) / mp.sqrt(V * d)
    sharp = 2 * mp.sqrt(mp.mpf(U - V) / U)
    worst = max(worst, dev); worstsharp = max(worstsharp, dev / sharp)
    print("    d=%d V=%d U=%d: max|root-V|/sqrt(Vd) = %.4f  (2 sqrt((U-V)/U) = %.4f)  max|Im| = %.1e" % (d, V, U, dev, sharp, max(abs(mp.im(r)) for r in roots)))
print("(5) Jacobi localization: worst ratio %.4f <= 8: %s; worst / (2 sqrt((U-V)/U)) = %.4f" % (worst, worst <= 8, worstsharp))

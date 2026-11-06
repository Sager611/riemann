# Which P_{1,m} satisfies the four-term recurrence?  (exact rational arithmetic)
import sympy as sp, random
random.seed(2)
def pF_coeffs(d, A, B, C, D):
    c = [sp.Integer(1)]
    for k in range(d):
        c.append(c[-1] * D * (k - d) * (k + A) * (k + C) / (A * C * (k + 1) * (k + B) * (k + D)))
    return c
y = sp.symbols('y')
tot = {'expanded(Savva app)': 0, 'decomposed(Savva 13.6)': 0, 'corrected': 0}
cnt = 0
for trial in range(8):
    d = random.choice([4, 6, 7, 9])
    A, B, C, D = [sp.Rational(random.randint(200, 900), random.randint(3, 9)) for _ in range(4)]
    p = sum(cf * y ** j for j, cf in enumerate(pF_coeffs(d, A, B, C, D)))
    der = [p]
    for k in range(d + 3): der.append(sp.diff(der[-1], y))
    yv = sp.Rational(random.randint(1, 60), 7)
    pv = p.subs(y, yv)
    T = [yv ** k * der[k].subs(y, yv) / pv for k in range(d + 3)]
    for m in range(0, d - 1):
        Dy = D * yv / (A * C)
        P3 = 1 - Dy
        P2 = B + D + 2 * m + 1 - Dy * (A + C + 3 * m - d + 3)
        P0 = -Dy * (m - d) * (A + m) * (C + m)
        P1e = (B + m) * (D + m) - Dy * ((A + m) * (C + m) + (m - d) * (A + C + 2 * m + 1))
        a = 1 - yv / A; ep = (C - D) / C
        bm = lambda mm: B + mm - yv + (d - 1 - 2 * mm) * yv / A
        cm = lambda mm: (d - mm) * (1 + mm / A) * yv
        beta = A * (2 * m + 1 - d) - d * (2 * m + 1) + 3 * m * m + 3 * m + 1
        P1d = cm(m + 1) + (D + m) * bm(m) + ep * (yv / A) * beta
        P1c = (B + m) * (D + m) - Dy * ((A + m) * (C + m) + (m - d + 1) * (A + C + 2 * m + 1))
        cnt += 1
        for nm, P1 in (('expanded(Savva app)', P1e), ('decomposed(Savva 13.6)', P1d), ('corrected', P1c)):
            if sp.simplify(P3 * T[m + 3] + P2 * T[m + 2] + P1 * T[m + 1] + P0 * T[m]) == 0: tot[nm] += 1
print("recurrence satisfied (out of %d (d,m,y) cases):" % cnt, tot)
# symbolic difference of the forms
A, B, C, D, m, d, Y = sp.symbols('A B C D m d y')
Dy = D * Y / (A * C)
P1e = (B + m) * (D + m) - Dy * ((A + m) * (C + m) + (m - d) * (A + C + 2 * m + 1))
P1c = (B + m) * (D + m) - Dy * ((A + m) * (C + m) + (m - d + 1) * (A + C + 2 * m + 1))
a = 1 - Y / A; ep = (C - D) / C
bm = lambda mm: B + mm - Y + (d - 1 - 2 * mm) * Y / A
cm = lambda mm: (d - mm) * (1 + mm / A) * Y
beta = A * (2 * m + 1 - d) - d * (2 * m + 1) + 3 * m * m + 3 * m + 1
P1d = cm(m + 1) + (D + m) * bm(m) + ep * (Y / A) * beta
print("P1c - P1e =", sp.factor(sp.simplify(P1c - P1e)))
print("P1c - P1d =", sp.factor(sp.simplify(P1c - P1d)))
P2 = B + D + 2 * m + 1 - Dy * (A + C + 3 * m - d + 3)
P2d = bm(m + 1) + (D + m) * a + ep * (Y / A) * (A - d + 3 + 3 * m)
P0 = -Dy * (m - d) * (A + m) * (C + m)
Gam = m * (m - d) * (m + A)
P0d = (D + m) * cm(m) + ep * (Y / A) * Gam
print("P2 - P2d =", sp.factor(sp.simplify(P2 - P2d)), "   P0 - P0d =", sp.factor(sp.simplify(P0 - P0d)))

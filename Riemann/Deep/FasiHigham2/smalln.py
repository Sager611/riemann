
# small n: c_n(m,m) as rational functions of m; prove c_n > 0 for all real m >= M0 by checking the
# numerator and denominator have only positive coefficients after substituting m = M0 + x.
import sympy as sp, sys
M0 = int(sys.argv[1]) if len(sys.argv) > 1 else 3
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 12
m, x = sp.symbols('m x')
def q(j):
    # q_j = (-1)^j C(m,j)(2m-j)!/(2m)! = (-1)^j / j! * prod_{i<j} (m-i)/(2m-i)
    e = sp.Integer((-1)**j) / sp.factorial(j)
    for i in range(j): e *= (m - i) / (2*m - i)
    return sp.together(e)
qs = [None] + [q(j) for j in range(1, NMAX + 1)]
c = [sp.Integer(1)]
ok = True
for n in range(1, NMAX + 1):
    cn = sp.together(-sum(qs[j] * c[n - j] for j in range(1, n + 1)))
    cn = sp.factor(cn)
    c.append(cn)
    num, den = sp.fraction(sp.together(cn))
    pn = sp.Poly(sp.expand(num.subs(m, M0 + x)), x); pd = sp.Poly(sp.expand(den.subs(m, M0 + x)), x)
    cn_num = pn.all_coeffs(); cd = pd.all_coeffs()
    allpos = (all(cf >= 0 for cf in cn_num) and all(cf >= 0 for cf in cd)) or (all(cf <= 0 for cf in cn_num) and all(cf <= 0 for cf in cd))
    nonzero = pn.degree() >= 0 and any(cf != 0 for cf in cn_num)
    print('n=%d  c_n =' % n, cn if n <= 4 else '(deg num %d, den %d)' % (sp.Poly(num, m).degree(), sp.Poly(den, m).degree()), ' positive for m>=%d:' % M0, allpos and nonzero, flush=True)
    ok = ok and allpos and nonzero
print('ALL small-n positive:', ok)

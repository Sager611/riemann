# cert.py: rigorous certificate for Lamzouri 2611.14486 Thm Main2 (exact rationals + Arb balls, python-flint)
# usage: uv run --no-project --with python-flint python3 cert.py cert.json
# f(x) = P(sqrt(2 pi) x/s) e^{-pi x^2/s^2}, P(u) = sum_k B_k H_{2k}(u) = p(u^2);  fhat(xi) = s Q(sqrt(2 pi) s xi) e^{-pi s^2 xi^2},
# Q(v) = sum_k (-1)^k B_k H_{2k}(v) = q(v^2)   [H_n(sqrt(2pi)x) e^{-pi x^2} has Fourier transform (-i)^n H_n(sqrt(2pi)xi) e^{-pi xi^2}]
import sys, json
from fractions import Fraction
from flint import fmpz_poly, fmpq_poly, fmpq, fmpz, arb, ctx
ctx.prec = 4000
C = json.load(open(sys.argv[1]))
d = C['d']; sF = Fraction(C['s']); lamF = Fraction(C['lam']); thF = Fraction(C['theta'])
q_ = lambda F: fmpq(F.numerator, F.denominator)
s, lam, th = q_(sF), q_(lamF), q_(thF)
B = [q_(Fraction(b)) for b in C['B']]
ok = True
def CHECK(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name); ok = ok and bool(cond)
CHECK('theta = 4/7 - 10^-5 and 1/2 < theta < 4/7', thF == Fraction(4, 7) - Fraction(1, 10**5) and Fraction(1, 2) < thF < Fraction(4, 7))
CHECK('0 < lambda <= 1', 0 < lamF <= 1)
# physicists Hermite polynomials (exact integers): H_{n+1} = 2u H_n - 2n H_{n-1}
H = [fmpz_poly([1]), fmpz_poly([0, 2])]
for n in range(1, 2*d):
    H.append(fmpz_poly([0, 2])*H[n] - 2*n*H[n-1])
P = fmpq_poly([0]); Q = fmpq_poly([0])
for k in range(d+1):
    P += B[k]*fmpq_poly(H[2*k]); Q += (B[k] if k % 2 == 0 else -B[k])*fmpq_poly(H[2*k])
# even polynomials -> p(w), q(w) with P(u) = p(u^2), Q(v) = q(v^2)
Pc, Qc = P.coeffs(), Q.coeffs()
CHECK('P, Q even', all(Pc[i] == 0 for i in range(1, len(Pc), 2)) and all(Qc[i] == 0 for i in range(1, len(Qc), 2)))
p = fmpq_poly(Pc[0::2]); q = fmpq_poly(Qc[0::2])
print('deg p =', p.degree(), ' deg q =', q.degree())
def to_int(r):
    den = fmpz(1)
    for c in r.coeffs(): den = den.lcm(c.q) if hasattr(den, 'lcm') else den*c.q//fmpz(__import__('math').gcd(int(den), int(c.q)))
    return fmpz_poly([int(c*den) for c in r.coeffs()])
def no_real_root_at_or_above(r, w0, label):
    """all real roots of r (rigorous Arb isolation, flint complex_roots) lie strictly below w0; non-real roots have Im != 0"""
    R = to_int(r).complex_roots()
    tot = sum(m for _, m in R)
    real = [(c, m) for c, m in R if c.imag == 0]
    bad = [c for c, m in real if not (c.real < arb(w0))]
    nonreal_ok = all((c.imag > 0) or (c.imag < 0) for c, m in R if not (c.imag == 0))
    print(f'  {label}: degree {r.degree()}, roots counted {tot}, real roots {len(real)}: {sorted(float(c.real.mid()) for c, _ in real)}')
    return tot == r.degree() and nonreal_ok and not bad
def sturm_count(r, a):
    """number of distinct real roots of r in (a, oo) by Sturm's theorem, exact (primitive integer remainders)"""
    def prim(f):
        g = to_int(f); c = abs(g.content()) if g.degree() >= 0 else fmpz(1)
        return fmpq_poly([fmpq(int(x), int(c)) for x in g.coeffs()])
    seq = [prim(r), prim(r.derivative())]
    while True:
        rem = seq[-2] % seq[-1]
        if rem.degree() < 0: break
        seq.append(-prim(rem))   # positive rescaling of -rem keeps the Sturm sign pattern
    def V(vals):
        vals = [x for x in vals if x != 0]
        return sum(1 for i in range(len(vals)-1) if (vals[i] > 0) != (vals[i+1] > 0))
    va = V([f(a) for f in seq]); vinf = V([f.coeffs()[-1] for f in seq])
    return va - vinf
# (1) f >= 0 on R  <=>  p(w) > 0 for w >= 0
CHECK('p(0) > 0', p(fmpq(0)) > 0)
CHECK('f > 0 on R: p has no real root w >= 0 (Arb root isolation)', no_real_root_at_or_above(p, fmpq(-1, 10**40), 'p') )
n1 = sturm_count(p, fmpq(0)); CHECK(f'f > 0 on R: Sturm count of roots of p in (0, oo) = {n1}', n1 == 0)
# (2) fhat(xi) <= 0 for |xi| >= lam  <=>  q(w) < 0 for w >= 2 pi s^2 lam^2 ; use rational w0 <= 2 pi s^2 lam^2
wl = 2*arb.pi()*arb(s)**2*arb(lam)**2
w0 = fmpq(int((wl.lower()*10**30).floor().unique_fmpz()), 10**30)
CHECK('w0 <= 2 pi s^2 lambda^2 (rational lower bound)', arb(w0) < wl)
CHECK('q(w0) < 0', q(w0) < 0)
CHECK('fhat < 0 on |xi| >= lambda: q has no real root >= w0 (Arb root isolation)', no_real_root_at_or_above(q, w0, 'q'))
n2 = sturm_count(q, w0); CHECK(f'fhat < 0 on |xi| >= lambda: Sturm count of roots of q in (w0, oo) = {n2}', n2 == 0)
# (3) normalisation and functional: fhat(0) = s q(0);  J = f(0) + 2 int_0^theta f(a) a da
#   = p(0) + (s^2/pi) sum_n p_n 2^n n! (1 - e^{-X} sum_{j<=n} X^j/j!),  X = pi theta^2 / s^2
fh0 = s*q(fmpq(0)); CHECK('fhat(0) = s q(0) > 0', fh0 > 0)
X = arb.pi()*arb(th)**2/arb(s)**2
eX = (-X).exp(); term = arb(1); en = arb(0); fact = arb(1); Jint = arb(0)
pc = p.coeffs()
for n in range(len(pc)):
    if n > 0: term = term*X/n; fact = fact*n
    en = en + term
    Jint = Jint + arb(pc[n])*arb(2)**n*fact*(1 - eX*en)
J = arb(p(fmpq(0))) + arb(s)**2/arb.pi()*Jint
ratio = J/arb(fh0)
print('f(0)/fhat(0) =', (arb(p(fmpq(0)))/arb(fh0)).str(15))
print('J(f)/fhat(0) =', ratio.str(20))
CHECK('normalised functional f(0) + 2 int_0^theta f(a) a da > 1', ratio > 1)
print('lambda =', C['lam'], ' s =', C['s'], ' d =', d)
print('ALL CHECKS PASS' if ok else 'SOME CHECK FAILED')

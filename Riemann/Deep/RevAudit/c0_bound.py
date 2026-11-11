# Rigorous (Arb) proof that |C0(z)| <= 0.924062 (= 1.463/(2 pi)^{1/4}, the constant behind 1.463 t^{-1/4}) for |z| <= 1, where
# C0(z) = cos(pi(z^2/2 + 3/8))/cos(pi z) (Gabcke). Gabcke's Satz 2.3.4(b) states only |C0| < 0.93 (which would give 1.4725 > 1.463).
# C0 is even; on [0,1] bisect adaptively (exact dyadic endpoints); on [1/2 - 1/64, 1/2 + 1/64] (removable singularity) use Cauchy's
# mean value theorem: C0(z) = N'(xi)/D'(xi) with xi between z and 1/2.
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 80
PI = arb.pi(); BOUND = arb('1.463') / (2 * PI).root(4)
def A(q): return arb(q.numerator) / q.denominator          # exact for dyadic q
def N(z): return (PI * (z * z / 2 + arb(3) / 8)).cos()
def D(z): return (PI * z).cos()
def Np(z): return -PI * z * (PI * (z * z / 2 + arb(3) / 8)).sin()
def Dp(z): return -PI * (PI * z).sin()
L0, L1 = Fr(1, 2) - Fr(1, 64), Fr(1, 2) + Fr(1, 64)
def upper_abs(lo, hi):
    z = A(lo).union(A(hi))
    if lo >= L0 and hi <= L1: return abs(Np(z) / Dp(z)).upper()
    if hi > L0 and lo < L1: return None
    d = D(z)
    if d.contains(0): return None
    return abs(N(z) / d).upper()
stack = [(Fr(0), L0), (L0, L1), (L1, Fr(1))]; n = 0; worst = arb(0); worst_at = None
while stack:
    lo, hi = stack.pop(); u = upper_abs(lo, hi)
    if u is not None and u < BOUND:
        n += 1
        if u > worst: worst, worst_at = u, (float(lo), float(hi))
        continue
    if hi - lo < Fr(1, 2 ** 40): raise SystemExit("FAIL near z=%s" % float(lo))
    mid = (lo + hi) / 2; stack += [(lo, mid), (mid, hi)]
print("|C0(z)| < %s on [-1,1]: %d intervals certified; largest certified upper bound %s on [%.9f, %.9f]; cos(pi/8) = %s"
      % (BOUND.str(10), n, worst.str(10), worst_at[0], worst_at[1], (PI / 8).cos().str(10)))

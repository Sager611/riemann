# gl2audit: python-flint arb semantics that replay.py relies on ('>' is the certain predicate; lower(); fmpq -> arb encloses)
from flint import arb, fmpq, ctx
ctx.prec = 320
a = arb(1, 0.5); print("arb(1 +/- 0.5) > 0.9 :", a > arb(0.9), "| > 0.4 :", a > arb(0.4), "| < 1.6 :", a < arb(1.6), "| < 1.2 :", a < arb(1.2))
b = arb(fmpq(1, 3)); print("arb(1/3) contains 1/3:", b.contains(fmpq(1, 3)), "rad", b.rad())
l = (arb(2, 1e-10)).lower(); print("lower() of [2 +/- 1e-10]:", l, type(l).__name__, "exact:", l.rad() == 0)
print("bool of uncertain compare:", bool(arb(0, 1) > arb(0)), bool(arb(0, 1) < arb(0)))
e = arb(fmpq(-111, 128)).exp(); print("exp ball rad < 1e-90:", e.rad() < arb(10)**-90)

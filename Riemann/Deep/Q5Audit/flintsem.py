#!/usr/bin/env python3
"""q5audit: python-flint semantics used by q5cert.py (enclosures and strict comparisons)."""
import flint
from flint import arb, fmpq, fmpq_poly, ctx
print("python-flint", flint.__version__)
ctx.prec = 53
x = arb(fmpq(1, 3)); print("arb(1/3) at 53 bits:", x.str(20, radius=True), " rad>0:", x.rad() > 0, " contains 1/3:", x.contains(fmpq(1, 3)) if hasattr(x,'contains') else 'n/a')
e = arb(fmpq(15491, 10000)).exp(); print("exp(1.5491):", e.str(20), " rad:", e.rad())
print("strict < on overlapping balls:", bool(arb(1) < arb(1)), bool(arb("1 +/- 0.1") < arb("1.05")), " certain:", bool(arb(1) < arb(2)))
X = fmpq_poly([0, 1]); p = fmpq_poly([1, 2, 3]); print("fmpq_poly call with poly = composition:", p(1 - 2 * X), "| expected 12x^2-16x+6 ->", fmpq_poly([6, -16, 12]))
ctx.prec = 4000
s35 = arb(fmpq(3, 5)).sqrt(); print("4000-bit sqrt(3/5) radius:", s35.rad().str(3), "; 2sqrt(3/5) =", (2 * s35).str(15))
print("1.5491 < 2 sqrt(3/5) certain:", bool(arb(fmpq(15491, 10000)) < 2 * s35))

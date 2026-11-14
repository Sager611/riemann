#!/usr/bin/env python3
"""certcells.py (explaudit): run cert.py's own pnt_cell over its own cell loop and write every cell value to a file.
Run in a scratch dir: cp Riemann/Deep/ExplProg/cert.py cert_copy.py; uv run --no-project --with python-flint --with mpmath python3 certcells.py cells.txt (17 s, 15 MB).
Then: indep.py pnt cells.txt [n_random] (INDEP_PREC=64 for a fast run)."""
import sys
sys.path.insert(0, ".")
import cert_copy as C
from fractions import Fraction as Fr
out = open(sys.argv[1], "w")
X0 = C.A(Fr(416, 100)); XT = C.A(Fr(60)); n = 0
while X0 < XT:
    X1 = X0 + C.amax(C.A(Fr(1, 20)), X0 / 50)
    for kind in ("q1", "q2"):
        R = C.pnt_cell(kind, None, None, X0, X1)
        out.write("%s %s %s %r %r %.10g\n" % (kind, "-", "-", X0.mid().__float__(), X1.mid().__float__(), float(R.upper()))); n += 1
    L0 = C.arb(3).log()
    while L0 < X1:
        L1 = L0 + C.amax(C.A(Fr(1, 10)), L0 / 20)
        R = C.pnt_cell("gen", L0, L1, X0, X1)
        out.write("gen %r %r %r %r %.10g\n" % (L0.mid().__float__(), L1.mid().__float__(), X0.mid().__float__(), X1.mid().__float__(), float(R.upper()))); n += 1
        L0 = L1
    X0 = X1
out.close(); print("cells", n)

# Record: covering factor of non-separable families of homothetic triangles

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** (A) An explicit non-separable family of 6003 homothetic triangles with covering factor 500000000000000/470233160167623 = 1.0633022984. (B) Theorem: λ_simplex ≥ 1.0633083237, where λ_simplex is the supremum of the covering factors of non-separable families of homothetic copies of a triangle.

**Before.** 250000000/235141779 = 1.0631883498 with 303 copies (Lian–Xue–Yuan, arXiv 2609.06558, Sept 2026). Their §7 continuum model gave the value 1.06330832375 as a candidate, "not a proved lower bound", and left open whether finite families reach it.

**Method.** The paper's cyclic chain recurrence, optimized with 2000 bridges per chain at 60 digits, then every coordinate rounded so that all 12003 contacts hold strictly (A). Equal-step chains at the continuum parameters exist for every large number of copies; their covering factors converge, and an interval-arithmetic computation bounds the limit (B).

**Certificate.** verify.py checks A exactly in 0.2 s (standard library only); limit.py (mpmath interval arithmetic, 5 s) certifies B.

**Review.** Both results correct; significant but incremental. An earlier review flagged one rounding direction in limit.py; it now rounds outward and the bound is unchanged.

**Novelty.** No source proves a covering factor above 1.0631883; the value 1.0633083 was known only as an unproved candidate ([novelty](../Riemann/Deep/TriCover/NOVELTY.md)).

**Files.** [report](../Riemann/Deep/TriCover/REPORT.md), [note](../Riemann/Deep/TriCover/NOTE.md), [configuration](../Riemann/Deep/TriCover/config.txt), [verify.py](../Riemann/Deep/TriCover/verify.py), [limit.py](../Riemann/Deep/TriCover/limit.py)

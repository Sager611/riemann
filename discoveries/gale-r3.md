# Record: a stronger counterexample to Gale's conjecture in R³

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** An explicit 170-point set V ⊂ R³ with β₄(V)/D(V) = √(7932058815501876/10000000211897483) = 0.8906210556…: every partition of V into 4 parts has a part of diameter ≥ 0.89062·D(V).

**Before.** Gale conjectured that the ball, with √((3+√3)/6) = 0.8880738, is the worst case. The previous record was 0.8892811624 (the 101-point set of arXiv 2610.28828 and its code repository, Oct 2026).

**Method.** Interior refinement and pruning around the paper's construction. Not converged: more points would likely give more.

**Certificate.** Exact integer arithmetic: two exhaustive 4-colouring searches plus a DRAT proof checked by drat-trim; verify.py runs in 1 s.

**Review.** Certificate correct. An independent recomputation from the points alone agrees.

**Novelty.** Nothing higher found ([novelty](../Riemann/Deep/Gale/NOVELTY.md)).

**Files.** [report](../Riemann/Deep/Gale/REPORT.md), [points](../Riemann/Deep/Gale/points.txt), [verify.py](../Riemann/Deep/Gale/verify.py)

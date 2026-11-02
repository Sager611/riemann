# Record: a larger lower bound for the integer Chebyshev constant of [0, 1]

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** t_Z([0, 1]) ≥ 0.4222354990, where t_Z([0, 1]) is the integer Chebyshev constant (integer transfinite diameter) of [0, 1]. The bracket is now 0.4222354990 ≤ t_Z([0, 1]) ≤ 0.4226846976.

**Before.** 0.4222286000 ≤ t_Z([0, 1]) ≤ 0.4226846976 (arXiv 2610.00386, Oct 2026; earlier 0.4213 ≤ t_Z ≤ 0.422685).

**Method.** The paper's own lower-bound theorem and 206-polynomial library, with a new discrete measure (1832 cells: each of its support intervals re-graded with twice as many Chebyshev cells) and new multipliers from a re-solved linear program. A refinement of the existing method, not a new method.

**Certificate.** verify.py (about 40 s; C with FLINT/Arb, exact integers): inf A ≥ −0.862192066 on [0, 1/2] (20,328 intervals, by the exact symmetry x ↦ 1 − x), the G_i bounds by exact resultants, the mass condition B < M; both published certificates are re-verified as well.

**Review.** Correct; significant as a certified numerical record.

**Novelty.** Nothing at or above 0.42223549 found (arXiv 2610.00386 is at v3 with the old bracket), checked 2 Nov 2026.

**Files.** [report](../Riemann/Deep/IntCheb/REPORT.md), [note](../Riemann/Deep/IntCheb/NOTE.md), [certificate](../Riemann/Deep/IntCheb/CERTIFICATE_LOWER_new.json), [verify.py](../Riemann/Deep/IntCheb/verify.py)

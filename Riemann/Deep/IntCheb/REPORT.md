RECORD: t_Z([0,1]) >= 0.4222354990 (was 0.4222286000, arXiv 2610.00386, Oct 2026)

What t is: the integer Chebyshev constant (integer transfinite diameter) of [0,1].
Bracket now 0.4222354990 <= t_Z([0,1]) <= 0.4226846976 (width 4.49e-4, was 4.56e-4).

1. Reproduction of both published certificates (rigorous, FLINT/arb balls, exact integers):
 - Lower: inf_[0,1] A >= V = -0.86220840540944 (27 522 intervals, 31 s; true margin 1.2e-11),
   min G_i = -0.86220828217 (arb, exact resultants), total mass M-1 = -2e-17, B < M,
   library irreducible/primitive/distinct  =>  t_Z >= 0.42222860001915. Claim holds.
 - Upper: max_[0,1/4] f <= -1.72225754796595 < T (5 752 intervals, 12 s)  =>  t_Z <= 0.422684697505.
   Claim holds. (Their own double-precision scripts also pass: 125 s and 569 s.)
2. New lower bound: same theorem, same 206-polynomial library, new measure and multipliers.
 - LP re-solved by our own dense interior-point method with symmetric cutting planes;
   1832 cells (each of the paper's 23 support intervals on [0,1/2] re-graded with twice as many
   Chebyshev cells, gaps graded), 17 active multipliers.
 - Certified: inf A >= -0.862192066 (20 328 intervals on [0,1/2] by exact symmetry, arb),
   min G_i >= -0.86219116325, M = 1 + 6.8e-17, B - M = -1.08e-9 < 0
   =>  t_Z >= exp(-0.862192066/M) = 0.42223549904149.
 - Gain +6.9e-6. On the paper's own 947 cells the LP optimum is already -0.8622082833
   (theirs -0.8622084054), so most of the gain is the finer cells.
3. Upper bound: reproduced only; no new factor search (time).
Check: python3 verify.py  (about 55 s; needs cc, FLINT 3, GMP, MPFR, mpmath, sympy).
  Add --upper/--lower <published JSON> to re-check the paper's certificates too.
Novelty: arXiv search (Oct 2026) finds no bound newer than 2610.00386 v3.
  Incremental (finer discretization of the same method, which the paper predicts) -> reviewer to judge.
Next: 3x-4x cells or a Chebyshev-weighted basis per support interval (both the measure and the dual
  have inverse-square-root edges); grow the library via dual reduced costs + LLL; LLL factor search for the upper side.
Files: NOTE.md, verify.py, lowcert.c, gcert.c, upcert.c, CERTIFICATE_LOWER_new.json.

# Theorem: the diagonal Fasi–Higham positivity holds for every odd m

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** For every odd m ≥ 1 and every n ≥ 0, the n-th Taylor coefficient of 1/q_{m,m}(z) at 0 is positive, where q_{m,m} is the denominator of the diagonal [m/m] Padé approximant to e^z. This is the whole diagonal case k = m of Fasi and Higham's conjecture (SIAM J. Matrix Anal. Appl. 40 (2019), Conj. 3.2).

**Before.** Known: k = 0 for all odd m and eventual positivity for large n (Scardino, arXiv 2610.39051, 28 Oct 2026, which also shows the conjecture fails off the diagonal, e.g. m = 5, 7, 9 with k ≫ m, and calls the diagonal open); each odd m ≤ 361 certified one by one (Goodchild, Zenodo 23919117, 3 Nov 2026, which lists the general diagonal statement as open); m = 3 for every k (our earlier result, [fasi-higham-m3](fasi-higham-m3.md)).

**Method.** Computer-assisted proof. With ν = m + ½, write 1/q = (e^z/p)·Σ_l X^l, X = πI_ν(z/2)/K_ν(z/2) (Padé remainder with positive coefficients, Bessel-function identities). The coefficient c_n becomes I_L + Σ_l κ^l J_l, contour integrals with their own saddle points. For m ≥ 203 (ν ≥ 203.5), Liouville–Green bounds with explicit constants and saddle-point lower bounds give positivity in the three ranges n ≤ 12, 8 ≤ n < 60ν and n ≥ 60ν; every finite inequality is checked in Arb ball arithmetic on intervals, with monotonicity arguments covering whole parameter ranges. Odd 3 ≤ m ≤ 201 are certified individually (root balls, partial fractions, Arb power series); m = 1 gives 2^{−n}.

**Certificate.** [verify.py](../Riemann/Deep/FasiHigham2/verify.py) re-runs every step (about 25 minutes on one core; python-flint). Run by the author agent and again by the reviewer: all steps pass. Orchestrator spot check: m = 205, all c_n > 0 for n ≤ 2500 (mpmath at 1500 and 2500 digits agree); reviewer: exact rational c_n > 0 for m = 203 (n ≤ 6000) and m = 205, 251, 301, 351, 401 (n ≤ 5000), ball arithmetic to n = 40000.

**Review.** Independent adversarial review (agent fhrev, 5 Nov 2026): VERDICT CORRECT; every lemma re-derived, coverage of all n and all odd m checked, no gap; findings cosmetic ([review](../Riemann/Deep/FasiHigham2/REVIEW.md)). Significant: settles the diagonal case of an open conjecture whose off-diagonal form is false.

**Novelty.** No earlier proof found (5 Nov 2026): Goodchild's Zenodo record and GitHub repository (diagonal "still open"), Scardino 2610.39051 v1 ("remains an open problem"), arXiv math.CA, math.NA, math.CV listings for Oct–Nov 2026, citations of Fasi–Higham 2019.

**Files.** [report](../Riemann/Deep/FasiHigham2/REPORT.md), [proof](../Riemann/Deep/FasiHigham2/PROOF.md), [review](../Riemann/Deep/FasiHigham2/REVIEW.md), [verify.py](../Riemann/Deep/FasiHigham2/verify.py)


# Record: a larger lower bound for the integer Chebyshev constant of [0, 1]

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** t_Z([0, 1]) ≥ 0.4222354990, where t_Z([0, 1]) is the integer Chebyshev constant (integer transfinite diameter) of [0, 1]. The bracket is now 0.4222354990 ≤ t_Z([0, 1]) ≤ 0.4226846976.

**Before.** 0.4222286000 ≤ t_Z([0, 1]) ≤ 0.4226846976 (arXiv 2610.00386, Oct 2026; earlier 0.4213 ≤ t_Z ≤ 0.422685).

**Method.** The paper's own lower-bound theorem and 206-polynomial library, with a new discrete measure (1832 cells: each of its support intervals re-graded with twice as many Chebyshev cells) and new multipliers from a re-solved linear program. A refinement of the existing method, not a new method.

**Certificate.** verify.py (about 51 s for the new certificate, 105 s with both published ones; C with FLINT/Arb, exact integers): inf A ≥ −0.862192066 on [0, 1/2] (20,328 bisection nodes, by the exact symmetry x ↦ 1 − x), the G_i bounds by exact resultants, the mass condition B < M; both published certificates are re-verified as well.

**Review.** Correct; significant as a certified numerical record.

**Novelty.** Nothing at or above 0.42223549 found (arXiv 2610.00386 is at v3 with the old bracket), checked 2 Nov 2026.

**Audit.** AUDIT PASSED (intchebaudit, 19 Nov 2026, [IntChebAudit](../Riemann/Deep/IntChebAudit/REPORT.md)): every item holds as written; Theorem thm:cert of arXiv 2610.00386 (v3, still the newest) re-proved by hand (its middle term matters only for B = 1), the certificate rerun, and independent code finds the true min A = V + 1.7·10⁻⁹ and all 206 G_i within 10⁻¹⁶ of the resultant bounds; the certified value is exp(V/M) = 0.42223549904…, so the dependents' '0.4222355' (RECORDS, TREE) was rounded up and now reads 0.4222354990; orchestrator rerun of verify.py identical (PROVED ≥ 0.42223549904149, 20,328 nodes, peak 53 MB).

**RH relevance:** N10 (barrier record): t_Z([0, 1]) > 1/e already shows that integer Chebyshev polynomials cannot prove the prime number theorem (Pritsker's 0.4213 suffices); the record moves the constant, not the barrier.

**Files.** [report](../Riemann/Deep/IntCheb/REPORT.md), [note](../Riemann/Deep/IntCheb/NOTE.md), [certificate](../Riemann/Deep/IntCheb/CERTIFICATE_LOWER_new.json), [verify.py](../Riemann/Deep/IntCheb/verify.py)

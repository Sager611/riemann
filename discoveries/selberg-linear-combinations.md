# Explicit Selberg-method proportion for linear combinations of N ≥ 4 Dirichlet L-functions, of the order c/N (nodes N12, N8)

**Claim.** Let q ≥ 3, let χ₁, …, χ_N be distinct even primitive characters mod q (odd characters by Dousselin's remark that they are "handled similarly"), c_j ∈ ℝ∖{0}, and F(s) = Σ_j c_j ε_j q^{s/2−1/4} L(s, χ_j) (a function with a functional equation and no Euler product, of Davenport–Heilbronn type). Then the proportion κ_F of zeros of F on the critical line (liminf over T, zeros of odd order counted once) satisfies, unconditionally,

| N | κ_F ≥ | before (Dousselin) |
|---|---|---|
| 4 | 6.8158·10⁻⁴ | 3.68·10⁻⁹ |
| 5 | 5.7578·10⁻⁴ | 2.94·10⁻⁹ |
| 6, 7, 8, 9 | 5.0081, 4.4455, 4.0060, 3.6520·10⁻⁴ | — |
| 10 | 3.3601·10⁻⁴ | 1.46·10⁻⁹ |
| 100 | 5.1386·10⁻⁵ | 1.45·10⁻¹⁰ |
| 1000 | 6.8655·10⁻⁶ | 1.43·10⁻¹¹ |

and κ_F ≥ 3.425·10⁻³/N for every N ≥ 11 (9.696·10⁻³/N for N ≥ 10⁶). This is the order c/N in the number of summands: Selberg proved c/N² and announced c/(N log N), which Dousselin made explicit only for log N > 862.

**Before.** Dousselin (arXiv 2311.10285 = Forum Math. 37 (2025) 1871–1912, Theorem "application", TeX l. 619–640; thesis 2025): 7.38, 4.91, 3.68, 2.94, 1.46·10⁻⁹ for N = 2, 3, 4, 5, 10 and 2.161·10⁻⁶/(N log N)·(1 − … − 862/log N) for huge N. For N = 2, 3 our earlier discovery ([page](linear-combinations-three.md): 34.71%, 2.061%) is far better, so the record is N ≥ 4, where that method gives nothing (it needs a non-concentration lemma). No other explicit bound for N ≥ 4 was found.

**Method.** Dousselin's chain with his analytic inputs as proved in his paper (dominance decomposition, moment lemmas, Bauer's twisted second moment, his Lemma 12 and eqintero, Gabriel convexity) and these changes, each re-derived in [DERIVATION.md](../Riemann/Deep/SelbergN4/DERIVATION.md): (1) counting sign changes by meas(G)/H (main term 1/A instead of 1/(2A)); (2) his N = 1 Cauchy–Schwarz argument inside each dominance set, giving N(√C₁ + √C₂)² instead of 4N(C₁ + C₂); (3) a Minkowski M-estimate with the first moment, C₂ = 4(C₃ − 1)/(θϰ)²; (4) exact main terms from his own Lemma 12 (C₅′ = 2/√(πϰ(1 − θ)), 2.66 times below his contour bound); (5) a weight lemma Φ(t) = sup_{u+u′≥t} f(u)g(u′) from the divisor condition; (6)–(7) exact cancellation of local factors for small v and integration by parts twice for large v, which makes C₁ linear in A and turns c/(N log N) into c/N; (8) signed local δ-sums at every prime p ≤ 10⁶ (equal to −1/p or 0), so their Euler factors become exactly 1. Ablation at N = 4: 4.54·10⁻⁸ (steps 1–2), 4.15·10⁻⁵, 9.73·10⁻⁵, 6.82·10⁻⁴ (all steps).

**Certificate.** [selbergn4.py](../Riemann/Deep/SelbergN4/selbergn4.py): `uv run --no-project --with python-flint --with mpmath --with numpy python3 selbergn4.py` (300 s, 117 MB, one thread): Arb at 128 bits; Part 1 reproduces Dousselin's printed constants to all digits (3.682008·10⁻⁹ at N = 4); Part 2 certifies the new values with special functions as Arb balls, sups over boxes and upward rounding; asserts κ_F(4), κ_F(10) > 10⁵ × Dousselin ([cert.log](../Riemann/Deep/SelbergN4/logs/cert.log), [result.json](../Riemann/Deep/SelbergN4/logs/result.json)).

**Checks.** Rerun by the orchestrator on a copy in /tmp: identical log (only the run time differs) and identical result.json. The agent: normalisations of the exact main terms checked numerically (sum/main term 1.03–1.04 at q = 1, ξ = 10⁶; [mtcheck.log](../Riemann/Deep/SelbergN4/logs/mtcheck.log)); every constant of DERIVATION §11 found in the log (65 of 65); two slips in Dousselin's paper found (a conservative (log ϰ − ϰ) for (log ϰ − 1), and a dropped factor 2 at l. 1487 that changes his table only at N = 1); neither affects this result.

**Review.** Seven verdict-only reviews by the agent on pasted claims ([review.log](../Riemann/Deep/SelbergN4/logs/review.log)); fixes: a missing (q/φ(q)) in the text, an additive error made explicit, two exact constants restated as the upper bounds used. Orchestrator reviewer (report and certificate scheme inline, DERIVATION readable): CORRECT, SIGNIFICANT: "the c/N dependence improves on c/(N log N), not merely its numerical constant"; §8 (bounded d, integrable d′, Λ ≫ log T) justifies the linear-in-A estimate; §10's fixed-prime cancellation before absolute values is valid and the ordered limits control the tails; §9's negative terms are legitimate inside the positive weights; numerical enclosures accepted as certified.

**Novelty.** Agent ([NOVELTY.md](../Riemann/Deep/SelbergN4/NOVELTY.md)): arXiv listings math.NT/CV/CA October–November 2026, the arXiv API (Dousselin; "linear combination" with "critical line"; "Selberg's method"), HAL (the paper and the thesis), Zenodo, one GXL query: no explicit bound for N ≥ 3 besides Dousselin's, and no c/N anywhere; 2609.15329 (an optimised Selberg method for ζ only) is the closest modern work. Orchestrator: lst.py (math.NT October–November 2026: the only related title, Lamzouri's 2610.26841, treats N = 2 with 31%).

**Caveat.** It rests on Dousselin's analytic inputs as proved in his paper, which were not re-proved here; the new steps were refereed line by line by selbaudit (12 Nov 2026: AUDIT PASSED, every step (1)–(8) and §9 re-proved, the certificate code read against DERIVATION and rerun identically, falsification tests negative; §10 needed explicit uniformity statements, written in its FULLPROOF (viii) and confirmed by a reviewer in the stated iterated limits; [SelbAudit](../Riemann/Deep/SelbAudit/REPORT.md); orchestrator rerun of its δ-sum check identical), not by a specialist; §8 (integration by parts) and §10 (small primes) carry the largest gains. κ_F is a liminf over T for fixed q, χ_j, c_j, with no uniformity in q. The certified computation is for even characters (the constants do not depend on the parity). The factor N is structural: removing it needs a joint non-concentration lemma.

**RH relevance:** N12, N8: critical-line proportions for functions with a functional equation and no Euler product (the Davenport–Heilbronn family, where RH fails), and the first explicit bound of the right order c/N in the number of summands.

**Files.** [Riemann/Deep/SelbergN4/](../Riemann/Deep/SelbergN4/): REPORT.md, DERIVATION.md, NOVELTY.md, selbergn4.py (the certificate; `--search` runs the parameter search), checks/mtcheck.py, logs/.

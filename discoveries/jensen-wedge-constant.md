# Explicit constant in Holland's Jensen hyperbolicity wedge (node N6)

**Claim.** Let J^{d,n} be the degree-d Jensen polynomial of Riemann's Xi shifted by n. J^{d,n} has d distinct real (negative) zeros whenever n^3 log^2(n+2) >= 7000 d^5, for every d >= 1 (degrees d <= D_all = 9.0·10^24 are covered by the known all-n range from Chasse's sector lemma and RH to height 3·10^12; our d0 = 1.6141·10^25 also works). Without that range, K = 2.43·10^114 works for every d >= 1 (self-contained form of Holland's Theorem 1.1). With the all-n range this is the first explicit hyperbolic region unbounded in d.

**Before.** Holland (arXiv 2608.08682, 9 Aug 2026, TeX l. 81-89) proves the wedge n^3 log^2(n+2) >= K d^5 with an absolute K that is never computed (l. 1814-1823, 2042-2043); Griffin-Ono-Rolen-Zagier, Griffin-Ono-Rolen-Thorner-Tripp-Wagner (c e^{d/2}, compactness) and Kim-Lee are ineffective; the only explicit region was the all-n range, bounded in d. Savva (Zenodo 22293642, unrefereed) claims a two-thirds wedge with constants left symbolic; Planat (HAL hal-05675741) fits C0 ≈ 0.020 numerically.

**Method.** Holland's chain is kept (Laguerre -> Jacobi -> finite-free comparison model matching R_0..R_4, Gershgorin localization, fifth-order residual multiplier, perturbed-Jacobi derivative ratios, stability proposition). Its analytic input, the complex GORTTW saddle asymptotic, is stated without proof or constant in GORTTW (1910.01227, TeX l. 514-518); we replace it by an explicit complex saddle bound on |z - n| <= n/10 (relative error <= 5.33·10^-20 for n >= N* ≈ 1.725·10^41). Every "sufficiently large n" becomes n >= N*, every O(.) an Arb ball enclosure over 1000 cells. Sharpenings: Brouwer instead of a 1/30 contraction, Gershgorin constant 2 instead of 8, a signed fifth-derivative bound (12.769 instead of about 84), Cauchy radius 1.4r, a + b + theta < 1 instead of each <= 1/4. Kres = (Cfin/thetaE)^2 (1 + lamL)^2 <= 6416.59; K = max(Kres, KN(9e24)) <= 7000.

**Certificate.** certify.py (python-flint arb, 4 s): "K = max(Kres, KN(9e24)) = 6416.58.. <= 7000: PASS" (L* = 90) and KSELF = 2.429·10^114 (L* = 80). Rerun by the orchestrator on a copy: both logs byte-identical.

**Review.** Orchestrator's verdict-only reviewer: first GAP (DERIVATION.md printed c3 <= 7.63 and beta <= 0.0533, rounded the wrong way; the certificate itself always used Arb ball endpoints); after the outward rounding of every quoted constant and recertification: CORRECT (C3 s1 < kappa0/4, deltaA and etaG bounds, a + b + theta = 0.97084 < 1, Kres < 6416.68 even from rounded displays), SIGNIFICANT (first explicit wedge constant valid for unbounded degree).

**Novelty.** lst.py 'Jensen|hyperbolic|Turan|Laguerre' (math.NT, CA, CV; Oct-Nov 2026): no explicit Jensen region (2610.38870 is a GRH criterion); Zenodo "Jensen polynomials" (60 hits): Savva 22293642 symbolic, Chen 22951677 and Chung 23987495 RH claims, Schieber 22936342 conditional, Dumas 22897068 transfers published bounds; HAL: Planat fitted, not rigorous; our files: R374 candidate only. Recheck 6 Nov 2026: unchanged.

**Caveat.** Holland's preprint is unrefereed; its chain is used as written, with the analytic input replaced by our own bound. Item 1 needs the all-n range d <= 9.0·10^24 (Chasse's lemma with Platt-Trudgian's RH height). Not formalized in Lean.

**RH relevance:** N6 (RH iff J^{d,n} is hyperbolic for all d, n; this is the first explicit hyperbolic region unbounded in d).

**Files.** [Riemann/Deep/HollandK/](../Riemann/Deep/HollandK/) (REPORT.md, DERIVATION.md, certify.py, certify.log, certify_self.log, check_identities.py, check_jacobi.py, sanity*.py, sanity.log).

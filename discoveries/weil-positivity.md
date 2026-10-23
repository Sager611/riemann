# Record: certified Weil positivity on half-width 1.55

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** Q(f) ≥ 10^{−112}‖f‖² for every complex C² test function f supported in [−1.55, 1.55], where Q is Weil's explicit-formula quadratic form (window length 3.1, x = e^{3.1}; both parities).

**Before.** 1.5 (ours, earlier); 1.495 (Elston, Zenodo 23174403, audited rigorous).

**Certificate.** Arb ball arithmetic; the full chain is in [PROOF.md](../Riemann/Deep/WeilRecord/PROOF.md). Parameters: prime powers 2..19, ‖C_p‖ ≤ 4.591281212, γ = 109/20, Ω = 1640, N = 3400 Legendre modes, 448 bits. After all quadrature and truncation errors, the smallest eigenvalue is ≥ 1.2e−111 (even part) and ≥ 4.0e−107 (odd part).

**Status.** Not refereed; not in Lean; not progress toward RH.

**Review.** Rigorous (source and log audit of cert4.py and run4/).

**Novelty (2026-10-22).** Zenodo, GitHub (8 repositories) and arXiv monthly listings show no larger full-window certificate. Finite-matrix positivity results (e.g. a 47×47 truncated form at µ = 11) are not window certificates.

**Files.** [report](../Riemann/Deep/WeilRecord/REPORT3.md), [1.5 report](../Riemann/Deep/WeilRecord/REPORT2.md), [Elston audit](../Riemann/Deep/WeilRecord/ELSTON.md)

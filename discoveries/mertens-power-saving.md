# First explicit power-saving bound for the Mertens function (node N10)

**Claim.** For every x ≥ 1, |M(x)| ≤ 1.84·10^6 · x^0.9382, where M(x) = Σ_{n≤x} μ(n). Harnack-only version (one-line certificate): |M(x)| ≤ 7.91·10^5 · x^0.9409. Certified family down to |M(x)| ≤ 1.01·10^382 · x^0.8992; |M(x)| ≤ x^0.95 for x ≥ 10^395. The only input about zeros is the proved half-plane: ζ(s) ≠ 0 for Re s > 7/8.

**Before.** No explicit bound C·x^θ with θ < 1. A power saving needs a zero-free half-plane (M(x) ≪ x^θ makes 1/ζ analytic on Re s > θ), available only since OpenAI's 7/8 theorem (Oct 2026). Explicit records: Ramaré |M(x)| ≤ 0.013 x/log x (x ≥ 1,078,853), Cohen–Dress–El Marraki x/4345, the sub-power savings of arXiv 2208.06141; non-explicit M(x) = O(x^{7/8+ε}) (Lean, [PNTError](../Riemann/Deep/PNTError/)). The ψ analogue is explicit (Fiori–Jaskari, arXiv 2609.23222, with 7/8); it needs ζ′/ζ, not 1/ζ. The new bound beats Ramaré's from x ≈ 10^174 on. Missed at first (flagged by harvest5, 6 Nov 2026): Chirre–Helfgott, arXiv 2511.14736 (Nov 2025, TeX l. 132, 393), |M(x)| ≤ 3x/(π·10^10) + 11.39√x for every x ≥ 1 (from Platt's residues to height 10^10 + 1); it is smaller than ours below x ≈ e^{607} ≈ 10^{263.5}, ours beats it from there on; it is not a power saving, so the first-power-saving claim stands.

**Method.** u = log|(s − 1)ζ(s)/s| is harmonic on Re s > 7/8; Rademacher's explicit convexity bound gives a harmonic majorant h, so v = h − u ≥ 0. Harnack's inequality, or finite signed Poisson-kernel certificates found by linear programming (justified by the Herglotz representation), bound v on Re s = σ₁ by its values on lines Re s > 1, giving |1/ζ(σ₁ + it)| ≤ C|t|^a. A smoothed Perron formula moved to Re s = σ₁ (no residue: 1/ζ(1) = 0), unsmoothing over a window of length h ~ x^θ and a monotonicity argument give the bound for all x ≥ x_a; |M(x)| ≤ x covers the rest. θ = (σ₁ + a)/(1 + a).

**Certificate.** `batch.py all.json` (Arb ball arithmetic via python-flint, 120 bits, about 12 min) checks the certificate inequalities by adaptive bisection, the constants, the integral of |1/ζ| on 0 ≤ t ≤ 1000 and C at x_a for nine configurations. The orchestrator reran the rows H093, ML093g5 and ML092g8 on a copy (95 s): identical output.

**Review.** Reviewer (PROOF.md: majorant on the whole half-plane, directions of the Harnack and certificate inequalities incl. the linear term and negative weights, contour shift, unsmoothing, monotonicity): CORRECT; significant as the first explicit all-x power saving (the qualitative saving follows from the half-plane). The agent's own review pass had asked for three missing hypotheses (0 < a < 1; σ₁ ≥ 0.89), now stated.

**Novelty.** The agent scanned the arXiv math.NT listings for October (854 entries) and November 2026, keyword searches, Zenodo and our files. The orchestrator repeated (5 Nov 2026): listing titles (no Mertens bound; the Möbius titles concern short intervals, function fields, GCD forms, Chowla), an arXiv keyword search (only 2208.06141), Zenodo (no Mertens deposit since October). Nearest: Zenodo "Kernel-checked explicit Chebyshev bounds from a zero-free half-plane" (19 Oct 2026) and Fiori–Jaskari, both for ψ only.

**Caveat.** Uses 7/8; the Lean-checked 0.87338196661320752 would lower θ by less than 0.001 (not done).

**RH relevance:** N10 (explicit prime-error bounds): RH ⟺ M(x) = O(x^{1/2+ε}); this is the first explicit exponent below 1 for M, from the proved half-plane alone.

**Files.** [Riemann/Deep/Mertens78/](../Riemann/Deep/Mertens78/) (REPORT.md, PROOF.md, cert.py, run.py, batch.py, chk.py, all.json, c_*.json, results.json, verify.log).

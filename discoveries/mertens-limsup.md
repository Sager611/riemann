# Mertens oscillations: lim sup M(x)/√x > 2.07390 and lim inf M(x)/√x < −2.07641, given Weisz's zero data (node N10)

**Claim.** lim sup_{x→∞} M(x)/√x > 2.0739001417 and lim inf_{x→∞} M(x)/√x < −2.0764178466, where M is the Mertens function. Condition: the first 14,400 zeros of ζ and their coefficients 1/(ρζ′(ρ)) lie in the Arb balls published by B. Weisz (Zenodo 23651459); his own bounds rest on the same data. A conditional improvement over Weisz's reported bounds.

**Before.** 2.0617536515 / −2.0647155435 (Weisz, Zenodo 23651459, 25 Oct 2026, unrefereed, AI-produced, Jurkat–Peyerimhoff kernel with support γ₁₄₄₀₀ = 13550.64 at two lattice-found points y± with |y| ≈ 10^9626); refereed: 1.826054 / −1.837625 (Hurst 2018). No later claim on arXiv, Zenodo (concept 23298194 has v1 and v2 only) or HAL as of 6 Nov 2026.

**Method.** Ingham's theorem in the form of Odlyzko–te Riele: for a nonnegative kernel K whose Fourier transform k has k(0) = 1 and support in [−T, T], every value h_K(y) = 2Σ_{0<γ<T} k(γ) Re(e^{iγy}/(ρζ′(ρ))) lies between lim inf and lim sup of M(x)/√x, unconditionally (if RH fails or a zero is multiple, both limits are infinite). At Weisz's own points:
- The kernel is the best positive-definite one for the given phases: k = (g∗g̃)/‖g‖² for a step function g with 1,350 cells; the value is a Rayleigh quotient, so the best g is the top eigenvector of a Toeplitz matrix.
- The support grows from 13550.64 to 18292.5. The 5,911 extra zeros (14401–20311) have unknown phases at y±, so each enters with its worst case −2|k(γ)||c|; they cost only 0.0026–0.0028, against a gain of +0.0121 (lim sup) and +0.0117 (lim inf).

**Certificate.** certify.py (python-flint Arb, 282 s, 0.9 GB): phases γy mod 2π at 32,819 bits from the input balls, the integer kernels, the sums; zeros 14401–20311 isolated by acb.zeta_zeros with N(18292.5) = 20311 by Turing's method, and rigorous bounds for their |c|. Result: h_K(y₊) ≥ 2.07390014175980 ± 3.2·10⁻²⁴, h_K(y₋) ≤ −2.07641784661092 ± 4.5·10⁻²⁴. The same run reproduces Weisz's certified values to 10⁻⁴⁰. The orchestrator reran it on a copy (6 Nov 2026): "ALL CLAIMS VERIFIED".

**Data checks.** SHA-256 of Weisz's files; 80 zeros at 256 bits (index, γ, c); 50 zeros at 33,010 bits (the 40 heaviest by kernel weight, zeros 1–72, and 10 random ones up to γ = 13,192): one interval Newton step proves each stored ball holds the zero, and the c balls contain 1/(ρζ′(ρ)); 50/50. A full recheck of the other 14,350 zeros would cost about 6 core-days and was not done.

**Review.** Reviewer, pass 1: kernel admissibility, the worst-case tail and unconditionality valid; GAP: the imported balls are not fully re-certified. Pass 2, after the claim was restated as conditional with the trust chain documented: CORRECT, SIGNIFICANT ("improves both extremal bounds through kernel optimization"), to be listed as a conditional improvement.

**Novelty.** The agent: arXiv listings and API (Mertens, Ingham kernels, Jurkat–Peyerimhoff), Zenodo, HAL, Weisz's GitHub repository. The orchestrator repeated the listing and Zenodo scan (6 Nov 2026): nothing beyond Weisz. Kernels used before: Fejér (Ingham), Jurkat–Peyerimhoff (Odlyzko–te Riele, Best–Trudgian, Hurst, Weisz), Gaussian (counterexample size). Optimizing the kernel for fixed phases and charging the worst case beyond the high-precision range were not used for these constants.

**Caveat.** Conditional on Weisz's data (unrefereed, AI-produced); spot checks do not discharge the condition. New lattice points (aligning more than ~1,050 zeros) are out of reach on one core.

**RH relevance:** N10 (explicit oscillations of the Mertens function): quantitative Ω± bounds for M(x)/√x; the Mertens conjecture |M(x)| < √x fails by these bounds, and larger values constrain how the zeros' phases can align.

**Files.** [Riemann/Deep/MLimsup/](../Riemann/Deep/MLimsup/) (REPORT.md, certify.py, kernel.py, ext_zeros.py, spotcheck.py, kernel_plus.txt, kernel_minus.txt, certify.log, spotcheck.log, spotcheck2.log).

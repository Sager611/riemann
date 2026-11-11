BARRIER: Weil's zero sum is a quadrature whose node kernels all integrate to 2π‖f‖², so off-line zeros act only through the irregularity of their ordinates, and an explicit block of off-line pairs above 3·10¹² that respects every unconditional zero bound used here (B_2 half-plane, Trudgian's S and S₁ bounds, 2606.31441's count) makes it negative on every window with c ≥ 30 (c ≥ 47 at any verification height): Weil window records are not corollaries of zero verification (node N1).

# explore30 — quasi-Monte Carlo and kernel quadrature → N1

## Answer
1. Transfer (exact). With g = f * f^*, every zero ρ = 1/2 + δ + iγ is a node with kernel K(t, δ) = ∫ g(y) cosh(δy) e^{ity} dy:
   **W(f) = Σ_ρ K(γ_ρ, δ_ρ)**. Every kernel has the same integral 2π‖f‖² and is bandlimited to [−2a, 2a]. So by Poisson summation an off-line lattice finer than π/a gives exactly the on-line value, and off-line zeros enter only through ordinate irregularity (a Koksma–Hlawka split with discrepancy S(t)).
2. Barrier (information-based). Zero-side information determines W only up to an adversary. Two-bump test functions have negative phases (fraction arccos(sech 2δa)/π of each period π/a). Filling only those phases with off-line pairs and leaving gaps elsewhere is consistent with every unconditional zero bound used here, and it gives W < 0:
   - at T = 3·10¹², for a ≥ 1.7 (c ≥ 30);
   - at every height, for a ≥ 1.925 (c ≥ 47).

   So no zero-side certificate (zeros as quadrature nodes) can reach the window records c = 47, 71 or 100.
3. Converse. The certified c = 47 window does exclude such blocks at every height. It therefore carries Euler-product information beyond zero verification. This corrects TREE's "window records are finite zero verification": that holds for single pairs (explore19), not for patterns. The excluded patterns are exotic (hundreds of consecutive off-line zeros), so this is no RH progress.
4. N5 (Báez-Duarte as a worst-case error) is a re-encoding and was not pursued (DERIVATION §9): the target measure is the Möbius measure on the nodes, the fooling functions are Burnol's vectors, the prediction form is BSY, and every finite N is DH-blind.

## Transferred statement
Let a ≥ 1.7 and T = 3·10¹². There is a multiset Z of points in 0 < Re s < 1 and a test function f ∈ C_c^∞(−a, a) with Σ_{ρ∈Z} ĝ(γ_ρ) < 0, where Z has all of these properties:
- it equals ζ's zeros below T;
- it is symmetric under ρ ↦ 1 − ρ̄ and ρ ↦ ρ̄;
- it has no point with Re s > B_2;
- it satisfies |S| ≤ Trudgian 2014, the count bound of 2606.31441, and |S₁| ≤ Trudgian 2011 at every t, each with 1 unit kept as margin;
- it lies far inside the explicit zero-density bounds.

For a ≥ 1.925 the same holds with any T′ ≥ T in place of T. Hence for these a no proof of Weil positivity on [−a, a] can use the zeros alone, and such a proof must use the prime side.

Construction: f = φ((x − a′)/w) + e^{iθ}φ((x + a′)/w), φ = (1 − s²)⁸, w = 0.03a–0.06a, modulated to t_0 > T. Every zero in [t_0 − 40/w, t_0 + 40/w] becomes an off-line pair (δ = 0.37) placed in the negative intervals of K. Counts per period are D·π/a′, with second-order rounding that keeps S and ∫S near 0 (DERIVATION §5–6).

## Tests (Riemann/Deep/ExploreQMC/logs; numbers are W/(D·2π‖f‖²), where on-line zeros of density D give ≈ 1)
| test | result |
|---|---|
| exactness ∫K(t,δ)dt = 2π‖f‖², δ = 0, 0.1, 0.3734 | relative error 4.5·10⁻¹³ |
| off-line lattice (δ = 0.3734), pair step 0.467 / 1.538 (< π/a = 1.632) | 1.000000000000 / 1.000000000000 |
| same, step 2.0 / 3.64 (> π/a) | 1 + 2.8·10⁻⁷ / 2.77 (aliasing) |
| adversary at T = 3·10¹², all constraints, a = 1.7 / 1.8 / 1.925 / 2.131 / 2.303 | −0.53 / −0.66 / −0.65 / −0.79 / −0.91 (max|S| ≤ 4.67 of 5.66, max|∫S| ≤ 2.61 of 2.76) |
| same, a = 1.4–1.6 | not admissible (∫S 2.94–3.53 > 2.76); without the S₁ constraint admissible from a = 1.2 |
| t = 10¹⁰⁰, a = 1.85 / 1.925 / 2.131 / 2.303 | −0.69 / −0.75 / −0.95 / −1.14 (|S| 22.0 of 26.5, ∫S 8.9 of 14.7 at a = 1.925) |
| t → ∞, a = 1.925 (w = 0.03a) | needed \|S\|/D = 0.559 vs 0.594 allowed; \|∫S\|/log t = 0.037 vs 0.059; narrow-bump limit a ≥ 1.78 |
| smaller δ at T = 3·10¹² (a = 1.925 / 2.131 / 2.303) | δ = 0.2: −0.21 / −0.24 / −0.27; δ = 0.1: — / −0.055 / −0.068; δ = 0.05: — / −0.013 / −0.016 |

## Adversarial test (real zeros and DH)
| configuration (a = 1.925, localized two-bump f) | min W |
|---|---|
| ζ's zeros near t = 15000 (1286 zeros) | +1.071 |
| one count-preserving merged pair (δ = 0.3), worst of 72 phases | +0.898 (invisible) |
| block of K off-line zeros (δ = 0.3), K = 27 / 51 / 99 / 150 / 250 / 743 | +0.58 / +0.21 / −0.28 / −0.46 / −0.49 / −0.47 (count deviation ≤ 3) |
| Davenport–Heilbronn zeros to 420, 9 real off-line pairs (δ ≤ 0.37), a = 1.5 / 1.71 / 1.925 / 2.303 | +0.49 / +0.44 / +0.38 / +0.26 (W/2π‖f‖², min over phases, centres, widths) |

The window separates ζ's zeros from a patterned off-line block of ≈ 100 zeros at any height. Localized test functions do not see an isolated moved pair or DH's isolated pairs. Those are seen only by the low-frequency near-null modes: explore19's reach m* ≈ 4.4πc, and DH fails at c ≈ 30.7. The two mechanisms are complementary.

## Literature and novelty (quick)
- arXiv keyword searches (Weil positivity with small support, zero verification, quadrature with zero nodes) found nothing on this mechanism.
- TREE and NOGO carry only explore19's single-pair reach law.
- Inputs used:
  - Platt–Trudgian 2021 (RH to 3·10¹²);
  - Richey's Lean B_2;
  - Trudgian 2014 (S), Trudgian 2011 (S₁);
  - Hasanalizade–Shen–Wong 2022 and arXiv 2606.31441 Cor. improvedzeta (counts).
- No reviewer verdict yet.

## Caveats
- The barrier is relative to the explicit constants used. Asymptotically this adversary needs leading coefficients ≥ 0.089 (count) and ≥ 0.037 (S₁) at a = 1.925, so better bounds would remove it for large T, not at 3·10¹², where the constant terms dominate.
- The adversary family (two bumps, full gaps) is not optimized, so the true thresholds are lower. Only W < 0 is claimed.
- Floating point only, not certified. |W| is 0.2–1.1 of the main term. The Gauss–Legendre transform is exact to about 10⁻⁹⁰ for |uw| ≤ 120, and the outside contribution is ≤ 5·10⁻¹⁵.
- The 1-unit margins for the true S and S₁ at the block ends are a modelling choice: the ends are chosen where S_true ≈ 0.
- My first qmc1 run was unchunked and reached 5.5 GB. The orchestrator killed it. All kept runs are chunked, with a 250 MB self-abort: peak RSS qmc1 78 MB (494 s), qmc2 56 MB, qmc3 54, qmc4 54, qmc5 64, qmc6 45, qmc7 47 (/usr/bin/time -l, top samples in logs/*.top).

## Files
Riemann/Deep/ExploreQMC/:
- REPORT.md, DERIVATION.md (lemmas, construction, admissibility, thresholds, N5 note).
- qmc1_adversary.py: library, exactness, lattices, unbalanced scan.
- qmc2_real.py: δ-dependence; real zeros near 15000.
- qmc3_gap.py: needed gaps vs 4πA.
- qmc4_balanced.py, qmc5_thresh.py: balanced S/S₁ thresholds.
- qmc6_dh.py: Davenport–Heilbronn zeros.
- qmc7_delta.py: smaller δ.
- logs/.

Run: uv run --no-project --with numpy python3 <script>, run from this directory (it reads ../ExploreHelio/results).

## Suggested NOGO line
- explore30 (Orca, exploration: quasi-Monte Carlo and kernel quadrature — Koksma–Hlawka, Poisson/trapezoid exactness, RKHS worst-case error, information-based adversaries → N1, N5) BARRIER: Weil's zero sum is a quadrature W(f) = Σ_ρ K(γ_ρ, δ_ρ), K = (g·cosh δy)^, every node kernel integrating to 2π‖f‖² (off-line lattices finer than π/a are exactly invisible), so off-line zeros act only through ordinate irregularity; an off-line block above 3·10¹² in the negative phases of a two-bump f obeys B_2, Trudgian's S and S₁ and 2606.31441's count and gives W/(D·2π‖f‖²) = −0.53 … −0.91 for every c ≥ 30 (c ≥ 47 at any verification height), so the window records are not zero verification (they exclude such blocks at all heights), while an isolated moved pair and DH's isolated pairs stay invisible to localized f (W ≥ 0.26; near-null modes only, explore19); N5's Báez-Duarte worst-case-error form is a re-encoding (Möbius target on the nodes, Burnol = fooling functions, finite N DH-blind) [ExploreQMC](Riemann/Deep/ExploreQMC/REPORT.md)

RH relevance: N1 (barrier for zero-side certificates; corrects "window records are finite zero verification" for patterned off-line blocks); N5 examined, re-encoding.

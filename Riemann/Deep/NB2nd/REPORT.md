BARRIER: under RH, simple zeros, BCF's Σ|ζ′(ρ)|⁻² ≪ T^{3/2−δ} and (H′), the optimum's mean second-order coefficient is C₂ = (J − Σ_∞) + D₀ − κ̄. Here J − Σ_∞ ∈ [0.034, 0.0397] is explicit, so no N-independent small-k correction cancels J, and D₀ = −0.0131. The uncontrolled term is κ̄, the mean gain of N-dependent corrections near k = N (G_N⁻¹ applied to the zero-sum vector e_N). C₂ = 0 ⟺ κ̄ = J − Σ_∞ + D₀ ∈ [0.021, 0.027]. This is undecided; the value measured at N ≤ 640 is κ̄ ≈ 0.029 ± 0.004, which favours C₂ ≈ −0.005 ± 0.006.

# nb2nd: the second-order term of the optimal Báez-Duarte distance (node N5, rate)

## Answer (Pyramid)
1. **J is not cancelled by the head alone.** Take any fixed (N-independent) Dirichlet polynomial B, put A = V_N + B/log N, and write L = log N. Then L²I(A) − CL = J − 2ℓ(B) + ‖ζB‖² + D(N) + o(1), with the explicit functional ℓ_n = ⟨ζ′/ζ, ζn^{−s}⟩ from eq. (2) of derivation.md. The best such B of length K leaves J − Σ_K, where Σ_K = ℓᵀG_K⁻¹ℓ. Computed: J − Σ_K ≥ 0.0397 for every K ≤ 640. Extrapolated: J − Σ_∞ ∈ [0.034, 0.038]. The best fixed construction therefore has second-order mean +0.021…+0.027, the K = 60 one measures 0.0295, and the optimum measures ≈ −0.002 (N ≤ 640).
2. **The remaining ≈ 0.03 is N-dependent and sits near k = N.** The parts identified so far:
   - a tilt/jump P(1) = j/L, with first-order gain −jC (exact residue computation, confirmed by local Fourier analysis) and cost j²X(N), where X(N) = ‖1 − ζM_N‖² ranges over 0.03–0.30;
   - end coefficients, squarefree and non-squarefree.
   The tilt alone recovers 43–60 % of this gain, and tilt plus the support [N e^{−1/2}, N] recovers 75–88 %.
   The mean of this gain, κ̄, is a nonlinear (inverse-Gram) functional of almost periodic zero sums, and no theory controls it. This is the barrier.
   Measured: κ(N) = 0.0275, 0.0311, 0.0273 (N ∈ [100,250], [250,450], [450,640]; s.e. 0.003–0.006 from correlated fluctuations), with no significant trend (a fit κ_∞ − a/L gives 0.032).
   This gives a conditional estimate (assuming no trend in κ beyond N = 640): C₂ = [0.034, 0.0397] − 0.0131 − [0.025, 0.033] ∈ [−0.012, +0.002]. The upper end J − Σ_∞ ≤ J − Σ_640 = 0.0397 holds because Σ_K is non-decreasing in K (Σ_640 computed in double precision, error < 1e−8). C₂ = 0 needs κ̄ ∈ [0.021, 0.027], within about 1σ of the measured 0.026–0.031.
   Both finite-N trends lower C₂(N) as N grows: J − Σ_N falls from 0.041 toward J − Σ_∞, and the ramp-corrected end gain κ(N) − (R(N) − D₀) rises (0.0255, 0.0281, 0.0284). So C₂ < 0 is likely, but not shown; the value +0.0025 at N = 2560 (from 2611.03924's V2) shows that the oscillation persists.
3. **Oscillation.** The optimum's C₂(N) has sd 0.003 (N ∈ [100,640]); nbrate certified sd 0.0009 on [800,1280]. The ramp's R(N) has sd 0.014. The optimum's oscillation fits the zero-pair frequencies γ₂−γ₁ = 6.887 (amplitude 0.0026) and γ₃−γ₁ = 10.876 (amplitude 0.0019), with R² 0.75 together with C₂ + C₃/L. The ramp's scatter is instead the single-zero remainder 0.29 N^{−1/2} cos(γ₁ log N + φ) (R² 0.95), and a fixed head correction removes it (sd 0.0005).
4. **The data cannot fix the sign of C₂.** At N ≤ 640 the finite-N couplings have the size of C₂ itself: κ(N) has sd 0.011–0.017, and the head–E coupling is 0.019 at N = 640.
   - Plain mean over N ∈ [200, 5120]: −0.003 ± 0.004.
   - Fit C₂ + C₃/L: C₂ = −0.021, C₃ = 0.12, R² 0.35, unstable.
   Both are consistent with nbrate's estimate C₂ ∈ [−0.01, +0.005] and with C₂ = 0.

## Derivation (details, tags per step: derivation.md)
| step | statement | status |
|---|---|---|
| 1 | L²d_N² − CL = J + R(N) − S_N, S_N = (Lr)ᵀG_N⁻¹(Lr), r = G_N c^ramp − b_N, r_n = ⟨R_N, ζn^{−s}⟩ | rigorous (finite N) |
| 2 | ℓ_n = −[log 2π (1 + log n − γ) + α₁]/n − Σ_{k≥1} ζ(2k+1)/(2k(2k+1)n^{2k+1}), α₁ = π²/12 − γ² − 2γ₁; L r_n → ℓ_n | rigorous under RH (residues at 0, −2k); limit checked to 0.003 |
| 2′ | Σ μ(n) log n · ℓ_n = J (formally: ζ′/ζ = ζ·Σμ(n) log n n^{−s}); conditionally convergent, ‖ζ′/ζ‖ = ∞ | formal identity, verified term by term |
| 3 | fixed B: L²I(V_N + B/L) − CL = J(B) + D(N) + O(N^{−ε}); B invisible to the zero residues (ζ(1−s) kills E's poles) | proof sketch (as nbrate: needs H′) |
| 4 | jump/tilt P(1) = j/L: (1 − j/L)²J + (1 − j/L)CL·(cross) + D + j²X − … ⇒ gain −jC, cost j²X(N) | heuristic; residues exact, local Fourier check |
| 5 | C₂(N) = (J − Σ_N) + R(N) − κ(N), κ(N) = −2eᵀG_N⁻¹ℓ + eᵀG_N⁻¹e, e = ℓ − Lr | rigorous identity; mean of κ uncontrolled |

## Numerical tests (double precision Gram, N ≤ 640, cond ≤ 1.6e6; d_N² at N = 640 agrees with nbrate's certified value: C₂(640) = −0.0055)
| quantity | value |
|---|---|
| Σ_K (K = 1, 5, 20, 80, 320, 640) | 2.4419, 3.4753, 3.5161, 3.5201, 3.5224, 3.5230 (J = 3.5627) |
| fits of Σ_K | Σ_∞ − 0.036/log K: Σ_∞ = 3.5287; Σ_∞ − A K^{−0.43}: 3.5249; log log K growth fits worse |
| fixed-B construction K = 5, 20, 60 (N ∈ [450,640]) | 0.0744, 0.0337, 0.0298 (sd ≤ 0.0015) vs J − Σ_K + D₀ ⇒ D₀ = −0.0131 ± 0.0003 (nbrate: −0.0142 ± 0.003) |
| head optimizer | β*_1 = −0.152 (P_opt(1) = 0.979 at N = 1280, nbrate 0.981); β*/μ from 0.63 (n ≤ 7) down to 0.19 at n ≈ K; Σβ*_n/n = −1.018 (A(1) ≈ 0, hence V1 ≈ V2) |
| optimum C₂(N) = L²d² − CL | +0.0050 (N ∈ [100,200]), −0.0015 ([200,400]), −0.0015 ([400,640]); nbrate −0.0062 ([800,1280]); 2611.03924 V2: +0.0028 (2560), −0.0091 (5120) |
| N-dependent gain κ(N) | +0.028 ± 0.017, +0.031 ± 0.011, +0.028 ± 0.015 (N ∈ [100,250], [250,450], [450,640]) |
| tilt check (N ∈ [400,640]) | L²·2I₁ = −1.146 ± 0.218 (theory −C − 2J/L = −1.15); L²I₂ = 0.161 ± 0.040 (theory X + J/L² = 0.161) |
| endpoint profile of the optimum (N = 640) | L(P_opt − ramp) = 1.48 (k ≥ 0.98N) → 1.0 (k ≈ N/10): tilt j ≈ 1.2–1.4 + head; P(1) ≈ 0.2 |
Test against 2611.03924's V2 values: V2 − V1 is 2.6e−4 in these units (Σβ*_n/n ≈ −1 forces A(1) ≈ 0). Their values at 2560 and 5120 straddle the
N ≤ 1280 mean, so they are consistent with an oscillating C₂(N) of amplitude ≈ 0.006, and they do not discriminate C₂ = 0 from C₂ = −0.02.

## Literature and novelty
- Read in TeX:
  - BCF 1211.5191: first order only.
  - Burnol math/0103058: liminf lower bound through orthogonal vectors (first order only).
    - His MA(q) toy model has an exact Gram formula. For Q = 1 − z our computation gives the ramp exactly 1/N and the optimum 1/(N + 2), so the optimum beats the ramp by 2 at second order. The reason is that the jump sits outside the ramp's band (degree N + 1).
    - For ζ the jump sits at the band edge, so the first-order gain is −jC with a half-delta overlap. The toy model is therefore not a model of the ζ endpoint.
  - 2611.03924: numerics only; its least-squares fit of BBLS §8 is quoted.
- Scan:
  - `lst.py 'Nyman|Beurling|Baez|Báez'` over math.NT/CA/CV/SP/PR/FA/CO, Oct and Nov 2026: 5 title hits, all Beurling algebras or densities.
  - Zenodo '"Nyman-Beurling" OR "Baez-Duarte"', 67 records, those since Sep 2026 checked. 24124207 (11 Nov) has a Pick-quantity "floor", not the d_N second-order term. 24054105 and 24088197 have small-N certificates. 22647730 (Veltri) and 22650745 were already handled by nbrate.
  - Our files (TREE, COROLLARIES, NOGO, Deep reports): only nbrate.
- New as far as found:
  - (i) the explicit head functional ℓ_n, and the fixed-correction formula (3) with J(B);
  - (ii) the computation that small-k corrections leave J − Σ_∞ ≈ 0.036 > 0, i.e. the optimum's near-cancellation is an N-dependent end effect;
  - (iii) the jump mechanism (gain −jC, cost j²‖1 − ζM_N‖²);
  - (iv) D₀ = −0.0131 ± 0.0003 from Gram data (spread over K and N; the finite-N bias is not bounded), 10× sharper than nbrate's;
  - (v) the barrier: C₂ = 0 is equivalent to an identity between an arithmetic constant and an endpoint zero-sum mean.
  All are conditional or numerical, and none is a theorem about C₂^opt.
- Suggested filing (orchestrator decides):
  - NOGO: "N5 rate: the mean second-order coefficient of d_N² is (J − Σ_∞) + D₀ − κ̄; the head leaves 0.036, and κ̄ (N-dependent end corrections through G_N⁻¹) is uncontrolled; finite N ≤ 5120 cannot fix its sign."
  - COROLLARIES: formula (3) with ℓ_n (routine refinement of BCF/nbrate).
  - RECORDS: D₀ = −0.0131 ± 0.0003.

RH relevance: N5

## Caveats
- Reviewer (nb2nd, algebra only): items (a) ℓ_n residues, (b) no residues in ⟨E, ζB⟩, (c) the tilt residue with the external (1 − j/L) factor (linear gain −jC), (d) S_N = Σ_N + κ(N): all CORRECT; one wording fix applied (the pole at s = 0 is double, not triple). Contour-limit justifications were not assessed.
- Formula (3) and the jump analysis inherit nbrate's (H′): uniformity of the O(N^{−ε}) leftover integrals. The pair-term domination for profiles with a 1/z tail (the jump) is not checked, since Σ_ρ|ζ′(ρ)|⁻¹|ρ − s|⁻¹ is borderline.
- Σ_K is computed only for K ≤ 640. Its convergence is an extrapolation: two fits agree, and a slow divergence fits worse but is not excluded.
  - A local Hardy-space argument says the poles of ζ′/ζ cannot be tracked profitably: tracking them with the one-sided family B ≈ 1/(ζ′(ρ)(s − ρ + η)²) costs +π/(2η) per zero, in units of 1/(2π|ρ|²).
  - The slow rate A/log K (fitted A = 0.036) has the order of A_th = Σ_ρ|ζ″/(2ζ′)(ρ)|²/|ρ|², which is 0.041 after 100 zeros (t8) and ≈ 0.075 with the tail. Local approximation of the regular part of ζ′/ζ at each zero by ζB would cost |a_ρ|²/log K, so a factor ≈ 2 is unexplained.
- κ̄ is measured only at N ≤ 640. Its fluctuations (sd 0.011–0.017) are correlated over about 20 consecutive N, so the quoted standard errors are rough.
- All numbers are double precision. For the quadratic forms used, errors are below 1e−8, far below the effects reported.

## Files
Riemann/Deep/NB2nd/:
- REPORT.md, derivation.md.
- Scripts (uv run --no-project --with numpy --with scipy --with mpmath). Run order: t0 → t1, t2, t3, t5 → t4, t6, t7, t8; intermediates go to .scratch/nb2nd/ and are regenerated in seconds. Peak RSS from /usr/bin/time -l is given per script, all ≤ 70 MB:
  - t0_gram.py: Vasyunin Gram G_640 → .scratch/nb2nd/G640.npz, 49 MB.
  - t1_gains.py: ℓ_n vs L r_n, Σ_K, gains by support, 56 MB.
  - t2_profile.py: head optimizer vs finite N, endpoint profile, all-N C₂(N), 56 MB.
  - t3_tilt.py: tilt/jump family, X(N), 50 MB.
  - t4_tail.py: tail-gain decomposition, 63 MB.
  - t5_head.py: Σ_K fits and head-optimizer shape, 50 MB.
  - t6_freq.py: coupling κ(N) and frequency fits, 69 MB.
  - t7_fixedB.py: fixed-B constructions vs formula (3), 41 MB.
  - t8_kappa.py: κ(N) for every N and A_th from 100 zeros, 70 MB.
- logs/: outputs.
- State: .scratch/nb2nd/STATE.md.

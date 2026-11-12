IDEA: At finite λ = 1/N Burnol's zero vectors give an explicit Christoffel-type bound B(N) ≤ D(1/N)² ≤ d_N² that matches the continuous Nyman–Beurling distance D within 4.5 %, so ζ's second-order term splits into an explicit zero sum (limit −0.1035) plus a flat integer-dilation penalty (log²N·(d_N² − B) = 0.064 ± 0.003, hence C₂ ≈ −0.04 < 0), while Davenport–Heilbronn's penalty is first order (0.020/log N) (node N5).

# explore38: Christoffel functions → N5 (Orca exploration, 12 Nov 2026, 07:05–08:30 UTC)

Derivation: [DERIVATION.md](DERIVATION.md). Code: [code/](code/). Logs: [logs/](logs/). L = log N, C = 2 + γ − log 4π = 0.0461914179.
d_N = Báez-Duarte distance (dilations 1/k, k ≤ N, L²(0,∞), unconstrained V1); D(λ) = BBLS's continuous distance (all dilations θ ∈ [λ,1]); d_N ≥ D(1/N).

## 1. Bridge (measured before choosing)
| node | Christoffel reading | verdict |
|---|---|---|
| N5 | Burnol's lower bound (math/0103058) is Grenander–Rosenblatt's Christoffel–Darboux method: his toy Y^N_α are CD kernels at the roots of an MA(q) symbol, and his ζ vectors Y^λ_ρ ⊥ C_λ are their analogue. He takes λ → 0 only. The transferable tool is the finite-n error term, i.e. the exact finite-λ kernel. | chosen |
| N1 | λ_min(window) = min over band-limited F of Σ_γ\|F(γ)\|²/‖F‖²: a lower sampling constant, bounded above by a Christoffel ratio in a zero gap. Its potential theory is Mori's g₂ (MoriRate, LostLaw, explore35). | already covered |
| N3 | z_ρ = 1 − 1/ρ maps the zeros to \|z\| = 1 (RH) with weights \|1 − z_ρ\|² = 1/\|ρ\|² (total mass C). Li's λ_n are moments of this measure; its Toeplitz positivity fails at 2 × 2 (R395). | no new quantity |

## 2. Transfer and statement
Theorem (unconditional; proof DERIVATION §1–4; reviewer CORRECT on claims 1–5, see §7). Let N ≥ 2 and let R be a finite, conjugation-closed set of simple zeros on the line. Then
  d_N² ≥ D(1/N)² ≥ B_R(N) := c_Rᴴ H_R⁻¹ c_R,  H_ij = (Y_j, Y_i), c_i = (χ₁, Y_i),
  (Y_ρ, Y_ρ') = (1 − e^{i(Φ_N(γ) − Φ_N(γ'))})/(i(γ' − γ)) − E_ρρ'(N),  (Y_ρ, Y_ρ) = Φ_N'(γ) − E_ρρ(N),  c_ρ = (ρ − 1)/ρ² − e_ρ(N),
  Φ_N(t) = t log N + 2θ(t) + 6 arctan 2t,  Φ_N'(γ) = log N + 2θ'(γ) + 3/|ρ|² = log N + log(γ/2π) + O(γ⁻²),
with E, e = O(N^{−5/2}) the explicit pole terms of the multiplier V(s) = (s/(1−s))³/χ(s) at s = −2, −4, … (DERIVATION §2–3).
- The kernel has the form of a model-space reproducing kernel with phase Φ_N sampled at the zeros (de Branges-type by analogy only: V has poles).
  A zero at height γ acts as a Christoffel function with effective bandwidth log N + log(γ/2π): the zero density adds to the bandwidth.
- Expansion (fixed R, L → ∞): B_R = Σ_R|ρ|⁻²/L − [S₁(R) + P_R(L)]/L² + O(L⁻³), S₁ = Σ_ρ(2θ'(γ) + 3/|ρ|²)/|ρ|², P_R(L) = Re Σ_{ρ≠ρ'} c̄_ρ c_ρ' E_ρρ'.
  P has the mean P₀ = Re Σ c̄_ρ c_ρ' i/(γ' − γ) and is almost periodic in L with frequencies γ_j − γ_k. All zeros: S₁ = 0.100516, P₀ = +0.0030 ⇒ c₂^B(∞) = −0.1035.
  The error constants depend on R; uniformity over growing R is not proved, so the limit is computed, not proved.
- Bug found on the way: Burnol's φ_{w,0} needs the boundary term sin(2πt)/(πt). His Thm pretheo assumes k ≥ 1, so this is a k = 0 extension, not an erratum. Without the term the Gram matrix is indefinite (eigenvalue −4.65 at N = 2).

## 3. Tests (ζ)
Smallest case R = {ρ₁, ρ̄₁}: B·L = 0.00778, 0.00844, 0.00897, 0.00965 at N = 10, 100, 1280, 10¹⁰ → 2/|ρ₁|² = 0.009998 (logs/s3_K.log).
At scale (K = 400 zero pairs exact + diagonal tail over 22,491 zeros + density; K = 200 → 600 changes B·L by ≤ 5·10⁻⁵; H positive definite, λ_min ≥ 1.04 at K = 50 for every tested N):
- B ≤ d_N² for every N ≤ 1280 (d_N from the Vasyunin Gram, equal to nbrate's certified values to 12 digits). B/d² = 0.49 (N = 10), 0.72 (100–250), 0.79 (1280).
- B·L = 0.0268 (10), 0.0326 (100), 0.0360 (1280), 0.0376 (10⁴), 0.0400 (10⁶), 0.0423 (10¹⁰), 0.0448 (10³⁰). These are explicit unconditional lower bounds for d_N² log N, by monotonicity of d_N also between grid points. Before: liminf ≥ C − 6·10⁻¹² only.
- Continuous problem. The Farey relaxation (all θ = a/b ∈ [1/N, 1] with b ≤ Bd, plus 1/k) is an upper bound for D(1/N)². Double precision, matching explore13's 40-digit 0.0290 and 0.0280.

| N, Bd | B·L | Farey·L ≥ D²L | d²L | B/Farey | share of d² − B that is continuous |
|---|---|---|---|---|---|
| 10, 30 | 0.02677 | 0.02803 | 0.05474 | 0.955 | ≤ 4.5 % |
| 20, 40 | 0.02918 | 0.03009 | 0.04944 | 0.970 | ≤ 4.5 % |
| 50, 60 | 0.03070 | 0.03151 | 0.04643 | 0.974 | ≤ 5.1 % |

  So D(λ)² ≈ B(λ): Burnol's guess that only the discrete spectrum survives holds at finite λ to within 2.6–4.5 %. The rest of d_N² − B is the price of integer dilations.
- Deficit δ(N) := L²(d_N² − B) = (integer-dilation penalty) + (continuous part ≤ 0.003):

| N | 5–20 | 20–50 | 50–100 | 100–250 | 250–450 | 450–640 | 640–1280 |
|---|---|---|---|---|---|---|---|
| δ mean | 0.0593 | 0.0618 | 0.0655 | 0.0664 | 0.0662 | 0.0673 | 0.0663 |
| sd | | | | 0.0016 | 0.0007 | 0.0008 | 0.0011 |

  Fit on [50,1280]: δ = 0.0644 + 0.0003 L (residual sd 0.0012). Flat.
- The N-dependence of the optimum's C₂(N) = L²d² − CL is that of B. Detrended correlation of L²B and L²d²: 0.993 on [100,1280] (0.986 on [100,400], 0.995 on [400,1280]); sd 0.0018 vs 0.0029. Regression C₂ = 0.1263 + 1.586(L²B − CL) − 0.118/L, residual sd 0.00046 (fitted on N ≤ 1280).
- Out of sample (2611.03924's V2; V2 − V1 ≈ +0.0003):
  - N = 2560: predicted +0.0033 (regression), −0.0017 (flat δ); actual +0.0028.
  - N = 5120: predicted −0.0086 (regression), −0.0100 (flat δ); actual −0.0091.
- Consequence: C₂ = c₂^B + δ_∞ with c₂^B(∞) = −0.1035 explicit and δ ≥ 0, so C₂ < 0 ⟺ δ_∞ < 0.1035. With δ flat at 0.064–0.067: C₂(∞) ≈ −0.04, approached slowly. The flat-δ model gives L²d² − CL = −0.013 (10⁴), −0.019 (10⁵), −0.024 (10¹⁰), −0.029 (10³⁰). Mean d_N² log N then lies below C. Prediction, not proof: δ_∞ is measured, not derived.

## 4. Adversarial
Same code (generic bound(); on ζ it reproduces s3 to 10 digits).
- Davenport–Heilbronn: 2406 on-line zeros and 109 off-line pairs below 2500 (ExplorePotential). Each off-line zero ρ₀ gets the interior vector V⁻¹Q_λVψ_{1−ρ̄₀} ⊥ C_λ (Beurling's floor vector, truncated). Poles of V_f at −1, −3, …
- L(s, χ₋₄): an Euler product with GRH true numerically; 767 zeros below 900 (smooth count 766.3).

| | ζ | L(s, χ₋₄) | DH |
|---|---|---|---|
| B/Farey (N = 10, 20, 40/50) | 0.955, 0.970, 0.974 | 0.980, 0.985, 0.986 | 0.982, 0.982, 0.982 |
| δ = L²(d² − B), N ≈ 10 → 600 | 0.059 → 0.067, flat from N ≈ 50 | 0.050 → 0.085, δ/L falls 0.021 → 0.013 | 0.043 → 0.126, δ/L = 0.020 at every N |
| fit | 0.0644 + 0.0003 L | 0.1134 − 0.189/L or 0.040 + 0.007 L | 0.0040 + 0.0190 L (N to 1000) |
| detrended corr(B, d²) | 0.993 | 0.962 | 0.987 |

- For DH, (d_N² − B)·log N = 0.0196–0.0204 in every window from N = 5 to 1000: a first-order gap. The continuous problem is matched as well as for ζ (B/Farey = 0.982), so ≥ 89 % of the gap (N = 10, 20, 40) is DH's integer-dilation penalty.
- d_N² log N alone cannot separate them: DH also climbs toward its constant (0.1669 at N = 1000), and its Beurling floor 2.5·10⁻⁴ shows only at N ≈ 10³³⁴ (explore13). The decomposition separates ζ and χ₋₄ (penalty second order or saturating) from DH (first order) at N ≤ 10³.
- Moved ζ pair: B treats it through the interior Christoffel quantity. Its change is invisible while (2β − 1) log N ≪ 1, and an off-line ζ zero above 3·10¹² has 1/|ρ|² < 10⁻²⁵. B itself carries no RH information beyond Beurling's floor.
- Heuristic, not tested: integer dilations need the Dirichlet coefficients of 1/f. These are bounded for ζ and χ₋₄ (μ, μχ) and grow for DH (abscissa of 1/f is max Re ρ > 1/2). If so, the penalty's order is Littlewood-type information. The RH and Euler-product effects are not separated: DH lacks both.

## 5. What it gives on N5 (TREE N5 line: "second-order coefficient ≈ 0 ± 0.01, sign open"; NB2nd: C₂ = (J − Σ_∞) + D₀ − κ̄, κ̄ uncontrolled)
- New bound: an explicit unconditional finite-N lower bound d_N² ≥ B(N) at every N (before: only liminf d_N² log N ≥ C); values above.
- Rate: C₂ = c₂^B + δ. The first part is an explicit zero sum, −0.1035 in the limit, and carries all the N-dependence (trend and oscillation). The second is sign-definite (δ ≥ 0) and flat (0.064 ± 0.003); at N = 10, 20, 50 at least 94.9 % of it is integer-dilation penalty.
  This replaces NB2nd's κ̄ (sd 0.011–0.017, sign of C₂ undecided) with a one-sided question: is δ_∞ < 0.1035? Data say yes by a margin of 0.037 ≈ 30 sd. Predicted sign: C₂ < 0.
- Rigorous range: C₂ ∈ [−0.1035 − osc, +0.027]. The lower end is unconditional once uniformity in R is proved; the upper end is under RH and NB2nd's fixed-B construction.
- Barrier side: Christoffel/zero-kernel bounds are RH-blind beyond Beurling's floor. They hold verbatim for DH, and a moved pair enters only through an interior Christoffel function. The RH-sensitive part of d_N is the integer-dilation penalty, which no kernel at the zeros controls.

## 6. Novelty
- Read in TeX: Burnol math/0103058 (λ → 0 only; Gram → identity; no finite-λ kernel, no second order).
- Checked BCF 1211.5191 (first order), NBRate, NB2nd, explore13 (Farey relaxation values only, below C), our TREE, NOGO and COROLLARIES.
- lst.py 'Nyman|Beurling|Báez|Burnol|Christoffel|de Branges|reproducing kernel…' over math.NT/CA/CV/SP/FA, Oct and Nov 2026: 9 title hits, none on NB rates. Zenodo (69 records): 24124207, 24054105, 20277868 (Gram splitting at M ≤ 60), 17957671 — none has a finite-λ Burnol kernel, a continuous-vs-integer split or second-order structure.
- arXiv search (3 queries) found the same.
- New as far as found:
  1. the exact finite-λ Gram (phase Φ_N, Christoffel bandwidth log N + log(γ/2π));
  2. the unconditional finite-N bound and the c₂^B = −(S₁ + P₀) = −0.1035 limit;
  3. D ≈ B to 2–4.5 %;
  4. the flat integer-dilation penalty and the C₂ < 0 prediction;
  5. the ζ/χ₋₄ vs DH contrast.
  (1)–(2) are routine refinements of Burnol; (3)–(5) are numerical.
- Side note for the orchestrator: Zenodo 17957671 (abstract) claims ‖(ζA_N − 1)/s‖_{H²} ≪ (log N)^C/√N for A_N supported on a dyadic prime block of size N. If its support length is polynomial in N, this contradicts BBLS/Burnol (distance² ≥ (C − o(1))/log(length), unconditional). Not read beyond the abstract: candidate for REFUTATIONS.

## 7. Caveats
- Reviewer (verdict only, claims inline; transcript in this run's subagent session):
  - (1)–(5) CORRECT with qualifications: k = 0 extension; removable limits at w = −2j; constants not uniform in R; H = conj of the (Y_ρ, Y_ρ') display, now fixed in DERIVATION §3.
  - (6) INCORRECT in the stated generality. The residue-only expansion fails for general X (counterexample (1 + 2^{−s})², double poles), the target needs V_f(0) = 0, and the norm interpolation was stated reciprocally. All fixed: §6 is restricted to gamma-factor X (DH, χ₋₄), where 1/X has simple poles −(2k+1) and V_f(0) = 0.
  - "De Branges-type" is an analogy.
- All B values are double precision, well conditioned (λ_min(H) ≈ L). The tail beyond K zero pairs is diagonal; its interaction error is ≤ 10⁻⁴ in δ (K-convergence). Farey values are double precision with cond up to 5·10⁶ (agree with 40-digit explore13 at N = 10). Nothing is interval-certified.
- δ's flatness is measured on N ≤ 1280 (L ≤ 7.2). A slow drift (|b| ≲ 0.02 in a + b/L) is not excluded, and nothing proves δ_∞ < 0.1035.
- DH's on-line zeros come from sign changes (ExplorePotential; count deficit matches its 218 off-line zeros). Zeros above 2500 are treated as on-line density (effect < 10⁻⁴). χ₋₄'s zeros come from sign changes of an Euler–Maclaurin Z (matches mpmath to 12 digits).
- Memory: every kept script peaked ≤ 120 MB (/usr/bin/time -l), one process at a time. Copying nb2nd/nbrate's Vasyunin outer-product loop peaked at 390 MB at N = 1280; chunking fixed it (78 MB). This is likely the cause of nbrate's cap violations.

## Suggested NOGO line
- explore38 (Orca, exploration: Christoffel functions — Grenander–Rosenblatt/Burnol CD kernels at finite λ, Máté–Nevai–Totik error terms → N5; N1 is Mori's sampling problem, N3 has no new quantity) IDEA, not a proof: Burnol's zero vectors at λ = 1/N have the exact kernel (1 − e^{i(Φ_N(γ)−Φ_N(γ'))})/(i(γ'−γ)), Φ_N = t log N + 2θ + 6 arctan 2t (diagonal log N + log(γ/2π)). This gives an unconditional d_N² ≥ B(N) (B log N = 0.0360 at 1280, 0.0423 at 10¹⁰) that matches the continuous NB distance within 2.6–4.5 % (Farey brackets, N ≤ 50). ζ's C₂(N) is −(S₁ + P₀) = −0.1035 in the limit, carrying all N-dependence (detrended corr 0.993; predicts 2611.03924's N = 2560, 5120 within 0.0008), plus a flat integer-dilation penalty 0.064 ± 0.003, so C₂ ≈ −0.04 if it stays flat (unproved). The kernel part is RH-blind (DH, moved pairs); the penalty is second order for ζ and L(s, χ₋₄) but first order for Davenport–Heilbronn (0.020/log N, N ≤ 1000) [ExploreChristoffel](Riemann/Deep/ExploreChristoffel/REPORT.md)

RH relevance: N5 (rate of the Nyman–Beurling/Báez-Duarte criterion: second-order term, its sign, and where arithmetic enters).

## Files
- code/: scripts, run with `uv run --no-project --with numpy --with scipy --with mpmath python3 <script>`, OMP etc. = 1. Intermediates go to .scratch/explore38/ and regenerate in minutes. Peak RSS is in parentheses.
  - s0_cache.py: zeros and θ cache (43 MB, 80 s).
  - s1_dN.py 1280: Vasyunin Gram and every d_N² (78 MB, 7 s).
  - s2_burnol.py: first version and the PD check.
  - s3_scan.py K Nspec TAG: B_K(N) with tail (≤ 102 MB).
  - s4_compare.py: B vs d_N, deficit, correlations.
  - s5_predict.py: regressions and predictions to 10³⁰.
  - s6_asym.py: S₁, P₀ (80 MB).
  - s7_dh.py 1000 300: generic bound(), ζ check, DH (94 MB, 13 min).
  - s8_dhcmp.py: windows and fits.
  - s9_chi4.py 600 900 300: L(s, χ₋₄) (120 MB, 3 min).
  - s10_farey.py which N Bd: Farey brackets for D (≤ 117 MB).
- logs/: all outputs (s10_farey.log carries peak RSS per run).

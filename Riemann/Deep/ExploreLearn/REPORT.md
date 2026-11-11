NO IDEA: a learning-curve formula is a function of the feature pool and the size n only, but Báez-Duarte's integers and random subsets of the same size and density differ 5–15× (d² ≈ 0.1, stalled); the formula fails at N = 2 (2.2×), and its continuum limit, Burnol's D(1/N), misses an N-independent second-order head constant ≈ 0.06 (node N5).

# explore31: learning theory (kernel-regression learning curves, omniscient risk, spectral bias) → N5

Orca exploration, 11 Nov 2026. Derivations: [derivations.md](derivations.md). Logs: [logs/](logs/).

## Answer (pyramid)
1. **No new bound.** The transferred claim is false at every N tested, N = 2 to 2880. The claim: a learning curve built from the pool's Gram spectrum and the target's coefficients predicts d_N, its constant C and its second-order term.
   - The formula is implemented correctly: Gaussian sketches of the same pool, its own setting, match it to 0–4% for N ≥ 10 (s8).
   - It fails at first order in both directions: it over-predicts the natural design 1.25–4.5×, and under-predicts random subsets of the same pool up to 9.5×.
   - The argument does not depend on any fit: such a formula depends only on (pool, n), so it must give one value to every n-subset. The natural subset and random subsets differ 4.8–15.5× (N = 10–640), so no such formula can come within a factor √ratio = 2.2–3.9 of both.
2. **Where it breaks: alignment (an explanation, not proved).**
   - Every natural feature ρ(1/(kx)) jumps on the lattice {1/m}. Möbius-type coefficients cancel the jumps at 1/m for 2 ≤ m ≤ N, so what remains sits at x < 1/N (D5).
   - A dilation θ ∉ {1/k} jumps on {θ/n}, which no other design point shares. At the Báez-Duarte density (θ ∝ θ⁻²) about one point falls in each [1/(k+1), 1/k] for every N, so a random design never aligns its head.
   - Designs with the leverage density of {1/k} therefore behave differently:
     - random designs, θ = 1 pinned: d² ≈ 0.09–0.11 from N = 10 to 640, no decrease;
     - the continuum θ ∈ [1/N, 1]: d² = D(1/N)², Burnol/BBLS's quantity, with second-order coefficient −0.06 against ≈ 0 for d_N.
3. **The second-order gap is a head constant.** c_head = C₂(natural) − C₂(continuum) ≈ 0.06–0.067 is N-independent (N = 10–160; its θ ∈ [1/5, 1] part, 0.055–0.061, stays flat to N = 2880). About 90% of it sits in θ ∈ [1/5, 1], mostly [1/3, 1/2]. nb2nd's fixed-head functional, evaluated at real dilations, reproduces 70% (N = 40) to 80% (N = 2880) of it, and that share rises slowly. In the continuum the head overshoots J: Σ ≥ 3.597 > J = 3.5627. So nb2nd's near-cancellation J − Σ_∞ ≈ 0.036 is specific to the integers.
4. **No RH information.** Davenport–Heilbronn shows the same picture: c_head 0.03–0.06, and the formula's misfit has the same size. The formula's only zero-sensitive content is already known: Burnol's per-zero cost (the leading term C) and Beurling's off-line floor (L → ∞). For DH that floor shows only near N ≈ 10^334 (explore13).
5. **N1 is rejected on paper.** Eigenvalue decay of windowed Weil forms is set by the gamma factor and the conductor (explore19, Mori). Any lower bound is a positivity certificate, i.e. finite zero verification.
6. **Side correction (known, not learning theory).** liminf d_N² log N ≥ C holds exactly and unconditionally (D6).
   - Burnol states this himself: "If the Riemann Hypothesis fails this result is true but trivial" (math/0103058, TeX l. 196–204).
   - His B_λ is span{ρ(θ/t) : λ ≤ θ ≤ 1} in L²(0,∞), so d_N ≥ D(1/N) transfers the bound to d_N.
   - TREE and COROLLARIES state the weaker C − 6·10⁻¹².
   - reviewer: CORRECT, ROUTINE; their applicability question is answered by the TeX check above.

## Transfer
- **Tool.** The deterministic-equivalent ("omniscient risk") learning curve for least squares with n features drawn from a pool. This is the same self-consistent equation as Bordelon–Canatar–Pehlevan, Simon et al. (eigenlearning) and the random-feature equivalents of Bach (2023) and Atanasov et al. (2024):
  - E_LC(n) = ‖χ‖² − βᵀ(Γ + ε)⁻¹β, with ε fixed by tr Γ(Γ + ε)⁻¹ = n;
  - Γ is the pool Gram and β the target coefficients (D2).
- **Node N5.** In the Mellin picture, d_N² is a least-squares fit of 1 by ζ(s)k^{−s} in L²(μ). The natural dilations have point density e^u and feature weight e^{−u} in u = log k, so they have **uniform leverage on [0, log N]** (D1). Their population covariance is that of the continuous family θ ∈ [1/N, 1]: a truncated Wiener–Hopf operator with symbol |ζ(½+it)|²/(2π|s|²).
- **Statement tested (LC).** For designs at the Báez-Duarte density (θ ∝ θ⁻² on [1/N, 1]):
  - (a) d² ≈ E_LC(N) (the omniscient-risk prediction);
  - (b) E_LC(N) → δ(log N)² + O(log N/N), where δ(L)² = dist(χ, span{ρ_θ : θ ∈ [e^{−L}, 1]})² = D(e^{−L})², the continuous Nyman–Beurling distance of BBLS and Burnol;
  - (c) hence d_N² = δ(log N)²(1 + o(1/log N)): same C and same second-order coefficient C₂ = lim L(d²L − C) as the continuum.
- **Not an analogy.** Every quantity is a functional of ζ on the critical line, computed exactly:
  - Vasyunin's Gram for rational dilations (D3), q-designs {q/k} → continuum;
  - nb2nd's ℓ(x) extended to real x (D4).

## Tests
**Smallest case, N = 2** (design {1, 1/2}):

| quantity | value |
|---|---|
| d_2² | 0.289965 |
| E_LC, with ρ_1 pinned | 0.6354 (2.2× too high) |
| random {1, θ}, mean ± sd | 0.574 ± 0.144 |
| continuum δ(log 2)² | 0.23015 (q → ∞, converged to 10⁻⁵) |

(b) already fails here: d_2² − δ² = 0.060.

**At scale.** Columns: natural d_N²; continuum δ² (q-design extrapolated a + b/q + c/q²); E_LC with ρ_1 pinned; random designs (θ = 1 + N − 1 uniform pool points, 200 reps unless noted); the ratio natural/E_LC. Logs: s1, s2b, s6.

| N | natural | δ² (q → ∞) | E_LC | random, mean (sd) | natural/E_LC |
|---|---|---|---|---|---|
| 3 | 0.093837 | 0.05850 | 0.3715 | 0.307 (0.117) | 0.25 |
| 5 | 0.036102 | 0.01751 | 0.1640 | 0.139 (0.039) | 0.22 |
| 10 | 0.023772 | 0.01184 | 0.0857 | 0.115 (0.033) | 0.28 |
| 20 | 0.016504 | 0.00980 | 0.0579 | 0.108 (0.039) | 0.29 |
| 40 | 0.012715 | 0.00825 | 0.0391 | 0.109 (0.041) | 0.33 |
| 80 | 0.010983 | 0.00748 | 0.0257 | 0.109 (0.040) | 0.43 |
| 160 | 0.009346 | 0.00673 | 0.0160 | 0.090 (0.034), Q = 64, 20 reps | 0.58 |
| 320 | 0.007954 | – | 0.0099 | 0.094 (0.046), coarse pool Q = 5 | 0.80 |
| 640 | 0.007017 | – | – | 0.109 (0.040), Q = 16, 6 reps | – |

- Standard errors of the random means: ≤ 0.003 for 200 reps, 0.008 at N = 160, 0.016 at N = 640.
- E_LC tracks random designs only up to N ≈ 5, then predicts convergence that random designs do not show.
- Without pinning θ = 1, random designs sit at 0.22–0.24 for every N ≥ 10, because missing the edge θ ≈ 1 costs ≈ 0.26 (s2).

**Control (s8, MMAX 1200, 50 reps, θ = 1 pinned).** Gaussian sketches ψ_i = Σ_j S_ij φ_j of the same pool:

| N | E_LC | sketch mean (s.e.) | natural | random subset |
|---|---|---|---|---|
| 10 | 0.08568 | 0.0892 (0.0049) | 0.02377 | 0.0991 |
| 20 | 0.05779 | 0.0596 (0.0016) | 0.01650 | 0.1093 |
| 40 | 0.03857 | 0.0377 (0.0005) | 0.01272 | 0.1028 |
| 80 | 0.02512 | 0.02509 (0.0003) | 0.01098 | 0.1078 |
| 160 | 0.015108 | 0.015112 (0.00011) | 0.00935 | 0.1041 |

The formula is right for generic features. The natural design beats these generic N-dimensional subspaces 1.6–3.8×, and random subsets lose 1.1–6.9×.

**Second order.** C₂(N) = L(d²L − C); c_head = C₂^nat − C₂^cont (s1, s3b).

| N | 2 | 5 | 10 | 20 | 40 | 80 | 160 | 640 | 2880 |
|---|---|---|---|---|---|---|---|---|---|
| C₂ natural | +0.107 | +0.019 | +0.020 | +0.010 | +0.003 | +0.008 | +0.006 | −0.0055 | −0.0050 |
| C₂ continuum | +0.079 | −0.029 | −0.044 | −0.050 | −0.058 | −0.059 | −0.061 | – | – |
| c_head | 0.029 | 0.048 | 0.063 | 0.060 | 0.061 | 0.067 | 0.067 | – | – |
| C₂, head θ ∈ [1/5, 1] continuum, Q = 64 | | | | | −0.051 | | −0.053 | −0.064 | −0.066 |

**Localisation (s3, N = 160, Q = 16).** Gain L²(d_nat² − d²) from adding all θ in:

| range | gain |
|---|---|
| [1/2, 1] | 0.013 |
| [1/3, 1] | 0.045 |
| [1/5, 1] | 0.054 |
| [1/10, 1] | 0.057 |
| [1/20, 1] | 0.058 |

N = 10 and N = 40 give the same profile.

**Fixed-head check (s7 against s3q/s3b), head θ ∈ [1/5, 1].** Measured gain against nb2nd's prediction Σ_{nat∪H} − Σ_nat, which uses ℓ(x) for real x:

| | Q = 8 | 16 | 32 | 64 |
|---|---|---|---|---|
| prediction | 0.0615 | 0.0698 | 0.0740 | 0.0759 |
| measured, N = 640 | 0.0469 | 0.0535 | 0.0569 | 0.0584 |

- The ratio measured/predicted is 0.70 (N = 40), 0.78 (160), 0.77 (640), 0.79 (1440), 0.80 (2880), the same for every Q. This is consistent with an O(1/L) finite-N coupling (nb2nd's head–E coupling).
- Σ_{nat∪H} = 3.5974 at Q = 64 exceeds J = 3.5627, while Σ_nat = 3.5215 (K = 160).

**Spectral bias (s5, G_1440).**
- Eigenvalues: λ_j j² ∈ [1.5, 10] (capacity exponent 2); cond 7.8·10⁶.
- 91% of the explained target variance lies in the top 5 eigenmodes.
- The top-N eigenmodes fit worse than the first N dilations for N ≤ 100 (ratio 1.3–2.5) and only slightly better for N ≥ 200 (0.9).
- Source/capacity rate theorems describe a fixed RKHS; here the span grows with N and the target lies outside every finite span. The rate 1/log N is a logarithmic source condition, and its constant is Burnol's boundary-zero cost, not a spectral quantity.

## Adversarial: Davenport–Heilbronn
Setup: exact digamma Gram (explore13), q-designs and random designs as for ζ, C_f = 0.19527 (s4).

| N | natural d² | continuum d² | c_head^DH | random+ρ_1 | E_LC | random/natural, DH | random/natural, ζ |
|---|---|---|---|---|---|---|---|
| 5 | 0.065911 | ≈ 0.0539 | 0.031 | 0.0762 | 0.0774 | 1.16 | 3.8 |
| 10 | 0.052861 | ≈ 0.0450 | 0.042 | 0.0714 | 0.0672 | 1.35 | 4.8 |
| 20 | 0.046131 | ≈ 0.0398 | 0.057 | 0.0701 | 0.0576 | 1.52 | 6.5 |

- Same failure: E_LC misses DH's natural design by 1.2–1.3×.
- The smaller alignment penalty comes from DH's coefficients: its jumps at 1/n have sizes a_n ∈ {1, 0.284, −0.284, −1, 0}, so the first misaligned jumps are small (D5). This is a coefficient effect, not a zero effect.
- The off-line zeros enter only through the Beurling floor, 2.5·10⁻⁴, which appears near N ≈ 10^334.
- A ζ pair moved off the line would change the formula and d_N alike, by its floor and a lost 1/(|ρ|² log N) term. Both are invisible at accessible N, so the test is moot once (LC) is false.

## N1 (not pursued)
Analytic-kernel eigen-decay theorems (Widom, Little–Reade) need an analytic kernel with decaying spectrum. On a window, Weil's form is an archimedean symbol growing like log|t|, minus a finite shift sum, plus a rank-2 polar part.
- Its near-null depth e^{−4πμ} (Mori) is fixed by the gamma factor and the conductor, so DH has the same mechanism (explore19).
- A lower bound on λ_min is a positivity certificate, i.e. finite zero verification (TREE N1 barriers).

## Review
reviewer (verdict only, evidence inline):
- **A, the side correction: CORRECT, ROUTINE.** Its one gap, whether Burnol's bound applies to this d_N, is closed by the TeX: d_N ≥ D(1/N).
- **B, the negative result as first phrased: NOT SUPPORTED.** Its points and how they were handled:
  - the alignment mechanism is not proved — it is now labelled an explanation;
  - finite Möbius sums leave the jumps at 1/m, m > N — wording fixed in D5;
  - the formula sees target projections, not only the density — wording fixed;
  - uncertainty estimates were missing — standard errors added.
- The remaining claim is the subset-independence argument (Answer 1), which needs only the measured ensemble gap. The reviewer called that ensemble mismatch plausible.

## Literature
- arXiv search (approach arxiv, four queries): no learning-curve, kernel-regression or deterministic-equivalent treatment of Nyman–Beurling/Báez-Duarte.
- Darses–Hillion (1805.06733, 2006.02953) randomize the θ_k, but inside expectation-type criteria, not as i.i.d. designs.
- BBLS 2000 and Burnol bound the continuum problem D(λ) from below only. Our q-design values are D(λ)² log(1/λ) = 0.1595, 0.0643, 0.0282, 0.0273, 0.0293, 0.0304, 0.0327, 0.0337 at 1/λ = 2, 3, 5, 10, 20, 40, 80, 160. They lie far below C and approach it slowly from below; Burnol's toy model expects order C/log(1/λ).
- explore13's Farey relaxation is the same continuum problem at N ≤ 20. Its 0.0280 at N = 10, B = 30 was unconverged; the q → ∞ value is δ²L = 0.0273.
- New as far as found, all numerical: the continuum second-order coefficient C₂^cont ≈ −0.06 of D(λ), the head constant c_head ≈ 0.06 with its θ-profile, the random-design stall, and the real-dilation head functional ℓ(x) (D4).

## Suggested NOGO line
- explore31 (Orca, exploration: learning theory — kernel-regression learning curves, omniscient risk/deterministic equivalents, spectral bias → N5, N1) NO IDEA: a learning-curve formula depends only on the feature pool and n (correct for generic features: Gaussian sketches match it to 0–4%), but Báez-Duarte's 1/k (uniform leverage, like all dilations in [1/N, 1]) beat it 1.25–4.5× (2.2× off already at N = 2) while random subsets of the same pool stall at d² ≈ 0.1 to N = 640 (5–15× worse than 1/k) — an aligned design (jumps on {1/m}); the continuum limit, Burnol's D(1/N), has C₂ ≈ −0.06 against ≈ 0 for d_N (N-independent head constant 0.06 from θ ∈ [1/5, 1], 70–80% of nb2nd's fixed-head functional at real dilations; continuum head overshoots J: Σ ≥ 3.597 > 3.563); DH alike; N1 eigen-decay is a gamma-factor effect [ExploreLearn](Riemann/Deep/ExploreLearn/REPORT.md)

## Suggested TREE/COROLLARIES fix (known, routine; not learning theory)
- N5: liminf d_N² log N ≥ C = 2 + γ − log 4π exactly, unconditionally. This replaces C − 6·10⁻¹² and needs no verification height:
  - Burnol (math/0103058, TeX l. 196–204) notes the bound is trivially +∞ if RH fails;
  - under RH, Σ m_ρ²/|ρ|² ≥ Σ m_ρ/|ρ|² = C;
  - d_N ≥ D(1/N) carries it from his B_λ to the integers;
  - equality forces RH and simple zeros.

RH relevance: N5

## Caveats
- All numbers are double precision. d_N matches nbrate's certified values (d_640² log 640 = 0.0453405057). The Gram conditioning is ≤ 10⁸, so quadratic-form errors are ≤ 10⁻⁸.
- Continuum values are q → ∞ extrapolations: a + b/q fits, consistent to ~3·10⁻⁵ in d² at N ≤ 20 and ~10⁻⁴ at N = 160. For N ≥ 320 only head refinements were run.
- The random-design stall is measured to N = 640 only (6 reps there). Under RH a Poisson design accumulating at 0 may still converge, extremely slowly; not tested.
- The fixed-head ratio 0.70 → 0.80 is extrapolated to 1 only by assumption (O(1/L)); not shown.
- Memory is measured as physical footprint (top's mem). RSS can be higher because it counts the pages of the memory-mapped Gram file: s2b at MMAX 1600 had footprint 246 MB and RSS 329 MB.
- **Memory-cap breach.** Two early runs of s3b_big.py at N = 2880 (Q = 16, 37 s; Q = 64, 9.4 min) peaked at 311–337 MB physical footprint, above the 300 MB cap. They were fixed (in-place LAPACK, one copy) and rerun at 173 MB (Q = 16) and ≤ 200 MB (Q = 64), with identical results.

## Files (Riemann/Deep/ExploreLearn/)
- REPORT.md, derivations.md, logs/ (outputs of every run below, including s8.log).
- Scripts:
  - run with uv run --no-project --with numpy --with scipy (s7 adds --with mpmath);
  - single-threaded; intermediates go to .scratch/explore31/;
  - peak memory is physical footprint from memcap.py / maximum RSS from /usr/bin/time -l.

| script | purpose | time | peak |
|---|---|---|---|
| s0_gram.py | Vasyunin Gram G_M (copied from nb2nd t0_gram.py); 2880 → G2880.npy | 25 s | ≤ 230 MB RSS |
| s1_qdesign.py | q-designs → continuum δ² | 6 s | 167 MB footprint |
| s2_lc.py, s2b_pinned.py | E_LC vs natural vs random | 15 s | ≤ 250 MB footprint (RSS 329 MB) at MMAX 1600; use 1200 for margin |
| s3_hybrid.py | head refinements; qext mode = Q-extrapolation | 30–50 s | ≤ 243 MB |
| s3b_big.py | head refinement at N ≤ 2880 | 35 s (Q = 16), 9 min (Q = 64) | ≤ 200 MB after the fix |
| s4_dh.py | DH adversarial | 16 s | 68 MB |
| s5_spectral.py | Gram spectrum, spectral truncation | 3 s | 126 MB |
| s6_random.py | random designs with a fine pool | 40 s (N = 160), 3 min (N = 640) | ≤ 116 MB |
| s7_headfun.py | ℓ(x) at real x, fixed-head prediction | 10 s | 66 MB |
| s8_sketch.py | control: Gaussian sketches against E_LC | 10 s | 100 MB footprint, 157 MB RSS |
| memcap.py | footprint check and abort | – | – |

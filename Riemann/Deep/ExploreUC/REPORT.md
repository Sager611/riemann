NO IDEA: quantitative unique continuation controls how small a solution can be on sets and how often it can vanish; for ξ and the de Bruijn–Newman flow H_t that is Jensen's formula in z and Angenent's zero-number theorem in (x, t), both already in use, and neither sees whether a zero is on the line or whether a zero is born at t > 0. Tests at Lehmer's pair (γ ≈ 7005): moving a zero pair off the line changes doubling indices by 0.01–1.6 %; Poon's frequency is ≥ 1 at creation points, but an ordinary zero already reaches 1.13 at s = 0.14; FE-symmetric Carleman windows see an off-line pair at height Y only if the window is narrower than about Y (0.78–1.00 Y); and the Harnack/Poisson bound never beats the break-even point of Raymond's derivative lever.

# Exploration explore2: quantitative unique continuation → RH tree (2026-11-05)

Field: Carleman estimates, Almgren/Poon frequency and the doubling index, Logunov–Malinnikova, the Landis conjecture, and quantitative and backward uniqueness for parabolic equations.
Sources: arXiv site search plus abs pages (no API), two Evidence-API claims (the third hit the daily quota), one scout (three more timed out), our NOGO/TREE/DBN6/HBDBN.
Four transfers were tested numerically and five were argued. None survives.

## Survey (2020–2026)
- **Landis, elliptic.**
  - LMNN 2007.07034 (one of R298's two citations): if Δu + Vu = 0 in ℝ² with real |V| ≤ 1 and u decays faster than exp(−C|x|√log|x|), then u ≡ 0.
  - **Frank–Ivanisvili 2608.00802 (Aug 2026): Landis is false in d ≥ 3.** There is a real bounded V with a solution decaying like exp(−c|x|^{4/3}), Meshkov's rate.
  - Further: 2512.09839 (plane, decaying real V); 2608.22698 (unique continuation at infinity for ∂̄u = Vu); 2609.21048 (reverse-Hölder V); 2506.19130 (a frequency-function approach to quantitative unique continuation); 2412.00788 (survey).
- **Nodal sets, elliptic.**
  - Logunov 1605.02587 and 1605.02589: polynomial upper bound; Nadirashvili's conjecture and Yau's lower bound. Both rest on doubling-index combinatorics, whose hard case is d ≥ 3.
  - 2609.11444: sharp volume bounds for the nodal and positivity sets of harmonic functions.
  - Logunov–Papazov 2012.07169 (R298's other citation): an elliptic Carleman–Domar / Levinson log-log theorem proved via three balls. It concerns growth near a line, not zeros of heat solutions.
- **Parabolic.**
  - Poon's frequency (1996); Colding–Minicozzi 2002.11015; Baldauf–Kim 2201.05505; ancient caloric functions on graphs 2412.15145; Sun–Wang 2512.10139 (frequency for the Ornstein–Uhlenbeck flow).
  - Huang–Jiang 2406.05877: H^{n−1}(Z_t) < ∞ in general and the nodal dimension is monotone. Measure monotonicity fails for n ≥ 2; for n = 1 it is Angenent's zero number.
  - Tao 1908.04958: quantitative Carleman / backward uniqueness for Navier–Stokes. 1809.00309: the zero number under general boundary conditions.
- **Bridge and novelty.**
  - arXiv "unique continuation" + zeta: only 2112.05521 (Oukil, a complex ODE that characterises zeros; not unique continuation in our sense).
  - "Carleman" + "zeta function" zeros: no hits. Evidence API: no supporting passage for unique continuation applied to dBN zeros or to zero creation under heat flow.
  - The nearest dBN work uses the off-zero log-derivative. Planat 2609.37164 gives d/dτ y² ≤ −2 − 4y²𝒢 with 𝒢 ≥ −Im(H'/H)(x + iη), the spatial half of a frequency function (Hakimoglu-Brown's mechanism).
  - **New, 4 Nov 2026: Planat–Raymond 2611.06245, unconditional Λ < 0.156821** by certified collision exclusion (unrefereed). It is not a record for us: 7/8 gives 0.116398 (Raymond) and 0.13963.

## Why the field stops short
1. **One complex variable.** The sharp tools are Jensen, Hadamard's three circles and the two-constants theorem. Logunov–Malinnikova machinery exists to recover these in d ≥ 3, where no conformal structure exists. The doubling index of ξ on a disc is ≈ (r/4) log(x/4π), the Jensen scale.
2. **(x, t) for H_t.** One-dimensional parabolic unique continuation is exact: Angenent's zero number, i.e. de Bruijn's theorem that real zeros stay real as t increases and new ones appear only at double zeros.
3. **Wrong kind of conclusion.** Unique continuation bounds vanishing order from above and gives lower bounds somewhere on balls. RH needs pointwise exclusion of simple off-line zeros; Λ ≤ 0 needs exclusion of order-2 points at t > 0. Both lie inside every such bound.
4. **Scale.** These quantities are computable only where we have control (Euler product for σ > 1, or t ≥ the known Λ bound). At those scales they average over many zeros, and a displaced pair is invisible (tests T1–T4).

## Transfers and tests
All tests use one dataset. Cauchy jets of H₀ = ξ((1+iz)/2)/8 come from 384 values on the circle |z − 14010.163| = 7 at 80 digits. The jets reproduce H₀ to 10⁻⁴². The flow is H_t = Σ_k (−t)^k ∂^{2k}H₀/k!, and it matches H_t(x) = E H₀(x + i√(2t) Y) to 5.6·10⁻¹⁶ at t = 0.1. Recovered Lehmer pair: γ = 7005.062866, 7005.100565.
Planted models replace two real zeros a, b by m ± iY, i.e. F = H₀·(z−w)(z−w̄)/((z−a)(z−b)). LP = Lehmer pair, Y = 0.05, 0.15, 0.40. NP = a normal pair (gap 2.70), Y = 0.30, 0.75 (the 7/8 maximum).

### T1 (N7, N2): doubling-index zero counts of H_t in discs or boxes
- **Transfer.** Bound the zeros in a box by the doubling index DI(r) = ln(M(2r)/M(r)); a Logunov–Malinnikova refinement would count sub-boxes.
- **Test** (base point 1.7, between zeros).
  - True H: DI(2) = 3.3569, Jensen bound 7.77, 3 zeros, 3 of them real.
  - LP0.15 (pair inside the disc): DI 3.3566, Jensen 7.75, 3 zeros, 1 real.
  - NP0.75: DI 3.3033 (−1.6 %).
  - The real-zero count changes by 2 while the doubling index moves by 0.01 %.
- **Verdict.** Blind. Doubling-index counts in boxes are Littlewood's lemma (Jensen on rectangles), which zero-density arguments (N7) already use. The sub-box combinatorics is trivial in one complex variable. Dies.

### T2 (N2): Poon frequency / quantitative backward uniqueness for zero creation
- **Transfer.** v(x, τ) = H_{−τ}(x) is an eternal solution of v_τ = v_xx.
  - At a base point (x₀, t₀), N(s) = s ∫H'_{t₀+s}² G_s / ∫H_{t₀+s}² G_s, with G_s the Gaussian of variance 2s, is nondecreasing in s, and d log ∫H²G/d log s = 2N (parabolic doubling).
  - N(0+) = ½ at simple zeros and 1 at collision or creation points (local model (x−x*)² − 2(t−t*)).
  - So N(s) ≥ 1 at every creation point for every s. Idea: certify Λ ≤ 0 by showing N < 1 at a scale where H_t is known.
- **Test.** N(s) at s = 0.1, 0.1164, 0.14, 0.2; all seven profiles are monotone, with N(0+) = ½ or 1 as predicted.

| base point | N(0.1) | N(0.1164) | N(0.14) | N(0.2) |
|---|---|---|---|---|
| H, ordinary zero, gaps 3.3/2.8 | 0.545 | 0.556 | 0.572 | 0.623 |
| H, ordinary zero in a cluster (gaps 0.99, 0.68) | 0.932 | **1.016** | **1.130** | 1.373 |
| H, Lehmer collision (h = 0.00043, t_c = −7.1120·10⁻⁴; −g²/8 = −7.1059·10⁻⁴) | 1.044 | 1.060 | 1.086 | 1.170 |
| LP0.15 creation (t* = 0.01110) | 1.044 | 1.059 | 1.085 | 1.169 |
| NP0.75 creation (t* = 0.2076; de Bruijn Y²/2 = 0.28) | 1.064 | 1.085 | 1.120 | 1.235 |
| model cos(kx), k = ¼ log(x/4π) = 1.754 | 0.562 | 0.583 | 0.618 | 0.730 |

- **Reading the table.**
  1. At s ≈ the known bound on Λ, an ordinary zero in a cluster (1.13) outranks two creation points (1.085, 1.120), even at γ = 7000.
  2. A past collision (t < 0) and a future creation (t > 0) have the same profile, so the frequency carries no sign of t.
  3. In the model, N_ordinary(s) = sk² coth(2k²s) > 1 once s > 0.9575/k². That holds at s = 0.14 above γ ≈ 2·10⁵ and at s = 0.1164 above γ ≈ 6·10⁵, while RH is verified to 3·10¹².
  4. At γ = 3·10¹² a certificate needs s < 0.021, i.e. H inside the uncertified time band at x-resolution ≈ 0.2. That is zero tracking, no cheaper than verification.
- **Verdict.** Dies. This is the quantitative form of R71's "backward stability ≡ RH". Log-convexity (Agmon–Nirenberg) bounds norms only.

### T3 (N7, N6): Carleman weights adapted to the functional equation
- **Transfer.** |H(x+iy)| = |H(x−iy)|, so I_w(y) = ∫|H(x+iy)|² w(x) dx is even in y and I_w''(0) = 2∫L₁w, with L₁ = H'² − HH''.
  - With w ≡ 1 it holds unconditionally (Plancherel: I = (π/2)∫e^{2tu²}Φ² cosh 2yu).
  - Pointwise it is the Laguerre inequality, part of the Jensen criterion (N6).
  - Question: is there a usable Carleman window in between?
- **Test** (Gaussian windows at the planted centre).
  - L₁ < 0 at every planted centre and > 0 for the true H.
  - The window average turns positive once the window's standard deviation exceeds 0.0500, 0.149, 0.379 (LP, Y = 0.05, 0.15, 0.40) and 0.290, 0.587 (NP, Y = 0.30, 0.75), found by bisection. That is 0.999, 0.993, 0.947, 0.966 and 0.782 Y.
  - Local model: f = (x−m)² + Y² gives L₁ = 2u² − 2Y², whose Gaussian average 2(v − Y²) changes sign at sd = Y exactly. Neighbouring zeros only lower the threshold.
- **Verdict.** An integrated, Carleman-type inequality sees an off-line pair only at windows narrower than its distance to the line, so it is the pointwise Jensen/Laguerre criterion in disguise. Dies.

### T4 (N2, record lever): Harnack/Poisson bound as Raymond's missing certificate G
- **Transfer.** DBN6 found that a certified uniform bound G ≥ mean ∂_hS at the probes (h = 0.69–1.88, t ≤ 0.117) could give Λ ≈ 0.115 < 0.116398.
  - Above all zeros, S = −Im(H'/H)(x+ih) = Σ (h−y_k)/|x+ih−z_k|² is positive harmonic. Harnack/Poisson gives a certified, uniform ∂_hS ≤ S/(h − y_max), with no arithmetic needed.
- **Test.**
  1. In the lever c = (2λ(μ)/h)(ℓ − M/h) − 2μ(G − g/h²), a gain needs (λ(μ) − λ(0))/μ > h(G − g/h²)/(ℓ − M/h).
  2. With the Poisson G the right side is ≥ h/(h − y)·S_max/S_min ≥ 1, and > 1 once y > 0.
  3. DBN6's kernel-table formulas on the grid r ∈ {0, 0.25, 0.5, 1, 1.5}, ρ ∈ {0.15, 0.3, 0.45}, μ ∈ {0.05, …, 0.7} give max (λ(μ) − λ(0))/μ = 1.0000. That is the far field, the same extremal configuration as Poisson (zeros horizontally far). So the gain condition (> 1 needed, ≤ 1 available) can never hold.
  4. At γ ≈ 7005 the true h∂_hS/S reaches 0.78, 0.62, 0.40, 0.29 at h = 0.69, 1, 1.5, 1.88 (above the 3.28 gap) against the Poisson bound 1. The usable margin is arithmetic.
- **Verdict.** Dies: at best break-even, never better. Any certificate for G must be sub-Poisson, i.e. use the Dirichlet series (DBN6's G_ideal).

### Argued only
- **Landis and Meshkov (Frank–Ivanisvili 2026).** No ζ-object solves a Schrödinger equation with bounded potential (−Z''/Z has poles at the zeros). The decay of Ξ is fixed by Γ, e^{−π|t|/4}. Dies at once.
- **Logunov–Papazov / Levinson log-log.** This needs a function holomorphic across the line whose bounds degenerate toward it from both sides. ξ is entire with no such degeneration, while Ξ'/Ξ and 1/Ξ have poles on the line. Dies.
- **Turán–Nazarov lemma (propagation of smallness for exponential sums).** This is Turán's power-sum method, used for zero density near σ = 1 since the 1950s. At most a constant gain near σ = 1, which is not RH-relevant. Not pursued.
- **Huang–Jiang and Sun–Wang (Ornstein–Uhlenbeck frequency).** In one dimension this is Angenent's theorem. The Ornstein–Uhlenbeck flow is explore1's Hermite attractor (wrong direction). Dies.
- **ℂ² structure.** H(z, t) is entire in (z, t). Creation points are real branch points of {H = 0} → t. Bezout/Bernstein-type counts bound how many there are, never zero. Dies.

## Ranking (RH relevance × chance × gain)
| rank | idea | node | relevance | chance | gain | result |
|---|---|---|---|---|---|---|
| 1 | T2: frequency < 1 certifies "no creation" | N2 | high (Λ ≤ 0 ⟺ no order-2 point at t > 0) | 0: an ordinary zero reaches 1.13 at s = 0.14 (γ = 7000); ordinary zeros exceed 1 everywhere above γ ≈ 6·10⁵ | RH | dies |
| 2 | T4: Harnack/Poisson G for the derivative lever | N2 | medium (record ≈ 0.115) | 0: ratio ≤ 1 ≤ h/(h−y), never beats break-even | ≤ 0.0017 | dies |
| 3 | T1: doubling-index zero counts | N7, N2 | low | 0: 0.01–1.6 % change | — | blind, = Littlewood |
| 4 | T3: FE-symmetric Carleman windows | N7, N6 | low | 0: needs window < Y | — | = Laguerre/Jensen |
| 5 | Landis, Levinson, Turán–Nazarov, ℂ² branch points | — | low | 0 | — | argued dead |

## Side notes for node owners
- **N2.** Planat–Raymond 2611.06245 (4 Nov 2026), unconditional Λ < 0.156821, is above the 7/8-based bounds, so it needs no tree change. Their "collision exclusion" is the T2 problem done by certified numerics, not by unique continuation.
- **N2 record line.** TREE.md still names 0.13963 as the record. DBN6 (2026-10-21) found Raymond's unrefereed 0.116398 under the same 7/8 hypothesis. The N2 owner should reconcile the two.
- **N2 tooling.** Cauchy jets on one circle (384 ξ-values, 45 s, 41 MB) give the whole local (x, t)-field: zeros, collisions, landings and frequencies. The isolated-pair landing estimate Y²/2 is beaten by the neighbours by 0.15–26 % at γ = 7000.

## Suggested NOGO line (orchestrator)
- explore2, quantitative unique continuation (Carleman, Almgren/Poon frequency, doubling index, Landis, backward uniqueness): blind to on/off-line and to zero creation. Moving a pair off the line changes doubling indices ≤ 1.6 %. Poon's frequency is ≥ 1 at creation points, but ordinary zeros reach 1.13 at s = 0.14 (γ = 7000) and exceed 1 above γ ≈ 6·10⁵. FE-symmetric Carleman windows need width < Y. The Harnack/Poisson bound for ∂_hS never beats the break-even of Raymond's derivative lever. [report](Riemann/Deep/ExploreUC/REPORT.md)

## Files
- `uc_main.py` (45 s, 41 MB): circle values, jets, validation and real zeros. It writes `jets.pkl` (deleted; rerun to regenerate).
- `uc_tests.py` (55 s, 43 MB): Lehmer collision, planted pairs, T1–T3 (T3 thresholds by bisection).
- `uc_t4.py` (5 s, 39 MB): Poisson ratios at γ ≈ 7005.
- `uc_t4b.py` (4 s, 73 MB): the DBN6 lever kernel table against the Poisson break-even.
- `outputs.txt`: all outputs, with times and peak memory.
- Runs were single-threaded, at nice 19, under a hard alarm ≤ 1500 s, one process at a time.

RH relevance: N2 (de Bruijn–Newman); secondarily N7 (zero density) and N6 (Laguerre/Jensen).

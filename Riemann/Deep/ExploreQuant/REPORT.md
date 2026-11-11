BARRIER: Lloyd–Max/Zador quantizer conditions see zeros only through their gaps (an off-line pair is one merged gap: 2.65–3.24 mean gaps at all 9 Davenport–Heilbronn pairs below height 420, below ζ's own unconditional Λ > 3.2648) or through log-concavity of |Ξ| between zeros (Fleischer's hypothesis = the first Laguerre inequality, N6), and quantizer uniqueness is weaker than that (node N8).

# ExploreQuant — vector quantization and rate–distortion → N8 (explore29, 11 Nov 2026)

## Transfer
- Tool: Lloyd–Max necessary conditions for MSE quantizers (Voronoi cells, centroid condition); Fleischer's uniqueness for log-concave densities (Fleischer 1964, Trushkin 1982, Kieffer 1983); Zador / Fejes Tóth 1-D optimum D_N ≈ N⁻²(∫p^{1/3})³/12 with point density λ ∝ p^{1/3}.
- Node: N8 (gaps between consecutive critical zeros; records λ ≥ 0.7851, Λ > 3.2648, μ ≤ 1/2). Bridge: zeros = quantizer points. Its only RH-sensitive input, log-concavity of |Ξ| between zeros, is N6's first Laguerre inequality.

## Statements (proofs: [DERIVATION.md](DERIVATION.md))
y_1 < … < y_N a Lloyd–Max stationary MSE quantizer of a continuous density p (cells V_n = [m_{n−1}, m_n], m_n = (y_n + y_{n+1})/2, y_n = centroid of p on V_n); gaps g_n = y_{n+1} − y_n.
- L1: p nondecreasing on an interior cell V_n ⇒ g_n ≤ g_{n−1} (equality iff p a.e. constant there); nonincreasing ⇒ g_n ≥ g_{n−1}.
- L2: a strict gap peak g_{j−1} < g_j > g_{j+1} forces a dip of p inside V_j ∪ V_{j+1} (x1 < a < x4 there with p(a) < min(p(x1), p(x4))); dips of distinct peaks are disjoint, so #peaks ≤ #dips(p); unimodal (every log-concave) p has no dip, hence no peak.
- L3: osc_{V_n} log p ≥ 2|log(g_n/g_{n−1})| on every interior cell (bang-bang extremal; = 4 artanh o_n, o_n = |g_n − g_{n−1}|/(g_n + g_{n−1})); so TV(log p) ≥ 2 Σ_n |log(g_n/g_{n−1})|: log p must vary twice as much as the log-gaps.
- Z: for unfolded zeros x_n = θ(γ_n)/π + 1, Voronoi MSE / Fejes Tóth optimum = κ₂ = ⟨δ³⟩ + O(1/N); in t the same ratio for Zador's p ∝ log³(t/2π), the only density whose optimal point density is the zero density.
- F (ζ): RH ⇒ |Ξ| log-concave on every (γ_n, γ_{n+1}) (Laguerre: Ξ'² − ΞΞ'' ≥ 0) ⇒ for all α > 0, k ≥ 1 and sub-intervals J the k-point MSE quantizer of |Ξ|^α on J is unique (Fleischer). Contrapositive test: a non-unique quantizer between consecutive real zeros certifies an off-line zero.

## Tests
### Smallest case
γ_2..γ_5 = 21.0220, 25.0109, 30.4249, 32.9351: gaps 3.9888 < 5.4140 > 2.5102, a strict peak, so by L2 no unimodal (no log-concave) density has these four zeros as consecutive stationary points. Three points never fail; γ_1..γ_4 (gaps 6.8873 > 3.9888 < 5.4140) do not. F holds on every true gap tested (max (log|Ξ|)''·ḡ² from −4.52 to −9.93; quantizers unique for k = 2, 3, 4, α = 1, 2).

### At scale (Odlyzko's zeros; quant_stats.py, baselines.py)
| set | n | ⟨δ²⟩ | κ₂ = ⟨δ³⟩ (dB) | ⟨δ⁴⟩ | corr(δ_n, δ_{n+1}) | strict-peak fraction | L3: mean 2\|log(δ_n/δ_{n−1})\| |
|---|---|---|---|---|---|---|---|
| zeros 1–10³ | 10³ | 1.1444 | 1.4587 (1.64) | 2.029 | −0.404 | 0.357 | |
| 10³–10⁴ | 9·10³ | 1.1552 | 1.4960 (1.75) | 2.127 | −0.366 | 0.355 | |
| 10⁴–10⁵ | 9·10⁴ | 1.1615 | 1.5155 (1.81) | 2.175 | −0.356 | 0.355 | 1.168 |
| 10¹²+ | 10⁴ | 1.1767 | 1.5651 (1.95) | 2.295 | −0.332 | 0.354 | 1.218 |
| 10²¹+ | 10⁴ | 1.1795 | 1.5744 (1.97) | 2.317 | −0.317 | 0.351 | |
| 10²²+ | 10⁴ | 1.1758 | 1.5642 (1.94) | 2.301 | −0.305 | 0.351 | 1.204 |
| sine kernel (exact) | | 1.179994 | 1.577939 (1.98) | 2.333092 | | | |
| CUE, N = 200 | 2·10⁵ | 1.1806 | 1.5798 | 2.337 | −0.307 | 0.3515 | 1.221 |
| Poisson | 2·10⁵ | 1.993 | 5.930 | 23.46 | 0.002 | 0.333 | median 2.20 |

Sine kernel by Nyström (48 Gauss nodes): ⟨s^k⟩ = k(k−1)∫ s^{k−2}E(0;s) ds, ⟨log s⟩ = −0.103463. Reading: the zeros are a 1.94–1.97 dB-inefficient quantizer of their own density at 10¹²–10²², approaching the sine kernel's 1.98 dB from below; a strict gap peak sits at 35 % of positions, so by L2–L3 a density with Lloyd–Max-stationary zeros has a dip every 2.8 zeros and TV(log p) ≥ 1.2 per mean gap, ≥ 0.19 T log T on [T, 2T], against 3 log 2/log T for Zador's log³ density. At 10²² every statistic agrees with CUE to 1 %.

### Adversarial
(a) GUE/CUE samples: the table — the gap statements do not tell ζ from random matrices.

(b) A ζ pair moved off the line (adv2.py, adv3.py): Ξ_δ = Ξ·((t − g0)² + δ²)/((t − γ_a)(t − γ_b)), g0 = (γ_a + γ_b)/2, δ in mean gaps ḡ; L = Lehmer pair #6709/#6710 (7005.0629, 7005.1006; merged on-line gap 3.011 ḡ), G = generic pair #6607/#6608 (6913.3036, 6914.1875; 3.497 ḡ). Columns: Laguerre (max (log|Ξ_δ|)''ḡ² on the merged interval; > 0 = fails), modes of |Ξ_δ|, sub-intervals (of 72) with a non-unique 2-point quantizer for α = 1, 2, distinct Lloyd limits for k = 3, 4 (α = 1; α = 2; 150 random starts), Gaussian water-filling rate R at θ = 10⁻² (S = |Ξ_δ|² normalised by the true mean).

| δ/ḡ | Laguerre L / G | modes | k = 2 non-unique L / G | k = 3, 4 limits L / G | R(θ = 10⁻²) L / G |
|---|---|---|---|---|---|
| 0.02 | +4988 / +4968 | 2 | 0, 0 / 0, 1 | 1,2; 2,2 / 2,2; 2,3 | 1.7987 / 1.9570 |
| 0.1 | +198.0 / +198.3 | 2 | 0, 0 / 0, 1 | 1,2; 2,2 / 2,2; 2,3 | 1.8378 / 1.9930 |
| 0.2 | +48.0 / +48.3 | 2 | 0, 0 / 0, 1 | 1,2; 1,2 / 2,1; 2,2 | 1.9941 / 2.1342 |
| 0.35 | +14.4 / +14.6 | 2 | 0, 0 / 0, 0 | 1,1; 1,1 / 1,1; 1,2 | 2.2686 / 2.3762 |
| 0.5 | +6.05 / +6.31 | 2 | 0, 0 / 0, 0 | 1,1; 1,1 / 1,1; 1,1 | 2.5232 / 2.6039 |
| 0.75 | +1.62 / +1.87 | 2 | 0, 0 / 0, 0 | 1,1; 1,1 / 1,1; 1,1 | 2.9071 / 2.9530 |
| 1.0 | +0.083 / +0.311 | 1 | 0, 0 / 0, 0 | 1,1; 1,1 / 1,1; 1,1 | 3.2464 / 3.2673 |
| true Ξ (per gap) | −9.86, −4.52 / −7.64, −9.93, −6.61 | 1 | unique | 1,1; 1,1 | 1.7954 / 1.9017 |

(c) Davenport–Heilbronn (adv.py, DH section; zeros from ../ExploreHelio): merged on-line gaps at its 9 off-line zeros below 420: 3.050, 3.118, 3.237, 2.809, 2.782, 3.013, 2.717, 2.650, 2.876 mean gaps (density log(5t/2π)/2π); DH's ordinary gaps (t > 20): max 1.928, 99th percentile 1.730. Sine-kernel tail P(s ≥ 2.65) = 5.2·10⁻⁴, P(s ≥ 3.24) = 8.5·10⁻⁶, P(s ≥ 3.2648) = 6.8·10⁻⁶; ζ near 10²²: 7·10⁻⁴ and 0 of 10⁴; ζ has infinitely many gaps > 3.2648 unconditionally (N8). First pair 0.8085 + 85.6993i: Laguerre fails on (83.1088, 87.6475) by +44.6 (ordinary interval: −8.8); 2 modes; k = 2 unique for α = 1 (whole + 52 sub-intervals), α = 2 non-unique on 4 of 52; Lloyd limits k = 3, 4: 1, 1 (α = 1), 2, 1 (α = 2).

Verdict. L1–L3 and Z read only gaps: blind to (a); they see (b) and (c) as a deleted point pair, i.e. as counting (Turing's method), and the merged-gap size certifies nothing because ζ has larger genuine gaps. F (k ≤ 4, α = 1, 2) sees (b) only for δ ≤ 0.2 ḡ (Lehmer) and δ ≤ 0.35 ḡ (generic), and (c) only at (α, k) = (2, 3) or on truncated sub-intervals, while its own hypothesis (log-concavity) fails in every case δ ≤ ḡ: F is strictly weaker than the first Laguerre inequality, which is already known on |x| ≤ T − 1/2 for all orders (arXiv 2611.02944). No new bound on any N8 or N6 quantity.

## Where it breaks
1. Continuity. An off-line pair is born in a collision: (t − g0)² + δ² with δ² passing from < 0 (two real zeros) through 0 to > 0. Quantizer distortions, Zador ratios, Lloyd residuals and Gaussian rate–distortion functions are continuous in δ² except through the number of real points; the only discontinuous ingredients are the zero count (Turing, Backlund) and the sign of (log|Ξ|)'' at g0 (−2/|δ²| − O(ḡ⁻²) for δ² < 0, a pole at δ² = 0, +2/δ² − O(ḡ⁻²) for δ² > 0), i.e. Laguerre.
2. No density. By L2–L3 the zeros are Lloyd–Max stationary for no explicit density, and asymptotically D/D_opt → 1.5779, not 1.
3. Direction. The mechanisms that make zeros quantizer-optimal (κ → 1) are forward flows: repeated differentiation (Farmer–Rhoades 2005) and the de Bruijn–Newman heat flow at t > 0, under which zeros relax to arithmetic progressions; Rodgers–Tao's proof of Λ ≥ 0 is this mechanism (Λ < 0 would make ζ's zeros the relaxed zeros of H_Λ, locally an arithmetic progression, against Montgomery's pair correlation, valid since Λ < 0 implies RH). Quantizer optimality therefore bounds Λ only from below (explore7's filter).

## Other tools screened
- Rate–distortion (reverse water-filling, Shannon lower bound) → N5, N1: R(D) and D(R) depend only on the distribution function of the spectral density S, exist for every S ≥ 0 and move smoothly with δ (table: R from 1.7954 to 1.7987 when the Lehmer pair moves 0.02 ḡ off the line). N5: Nyman–Beurling is the log-scale process with S = |ζ(1/2+it)|²/|s|²; RH ⟺ ζ(s)/s outer in H²(Re s > 1/2), a phase statement, and the Kolmogorov–Szegő/Jensen defect of S is the Balazard–Saias–Yor sum Σ_{Re ρ > 1/2} log|ρ/(1−ρ)| (1999): re-encoding. N1: a Gaussian source with Weil's covariance exists iff Weil positivity, so rate–distortion presupposes it.
- Lattice quantizers / Gersho → N12: the optimal 2-D quantizer A₂ has Epstein zeta 6ζ(s)L(s, χ₋₃), the non-optimal Z² has 4ζ(s)L(s, χ₋₄); forms of class number > 1 have off-line zeros (Davenport–Heilbronn 1936), elongated forms real zeros in (1/2, 1) (Bateman–Grosswald 1964): RH for Epstein zeta follows the class number, not quantizer optimality.

## Review and novelty
Reviewer (math only, inline brief): L1, L2, F and the smallest case OK; L3 needed the closed form ε(o) = 4 artanh o = 2|log(g_n/g_{n−1})| and the p > 0 caveat (fixed). Every table number re-checked against the logs by rg (rounding only). L1–L3 are elementary (L1 is folklore for monotone densities); the DH merged gaps and the Lloyd-versus-Laguerre comparison are tests, not discoveries. Nothing for DISCOVERIES.

## Suggested NOGO line
- explore29 (Orca, exploration: vector quantization and rate–distortion — Lloyd–Max/Fleischer, Zador/Fejes Tóth, reverse water-filling, lattice quantizers → N8, N6, N5, N12) BARRIER: quantizer conditions see zeros only through gaps, which match the sine kernel (⟨δ³⟩ = 1.564 at 10²² vs 1.5779; strict gap peaks at 35 %, so a density with Lloyd–Max-stationary zeros needs a dip every 2.8 zeros), and an off-line pair is one merged gap (2.65–3.24 mean gaps at all 9 Davenport–Heilbronn pairs below 420, below ζ's unconditional Λ > 3.2648); the RH-sensitive part, log-concavity of |Ξ| between zeros, is Fleischer's hypothesis = the first Laguerre inequality, and quantizer uniqueness is weaker (a ζ pair moved 0.35–0.5 mean gaps off the line keeps unique quantizers for k ≤ 4 while Laguerre fails); Gaussian rate–distortion is phase-blind, its RH trace is Balazard–Saias–Yor [ExploreQuant](Riemann/Deep/ExploreQuant/REPORT.md)

RH relevance: N8 (gap-based certificates), N6 (Fleischer = Laguerre), N5 (rate–distortion = Balazard–Saias–Yor).

## Files
- quant_stats.py: smallest case, κ_r, gap peaks, Lloyd tilts by height (logs/quant_stats.log; 6 s, 87 MB).
- baselines.py: exact sine-kernel moments, CUE N = 200 × 1000, Poisson, L3 bound (logs/baselines.log; 141 s, 52 MB).
- adv.py: DH merged gaps, DH Laguerre and 2-point uniqueness (logs/adv.log, DH section; its ζ section is void: |Γ(1/4 + it/2)| underflows at t ≈ 7000, redone in adv2.py; 546 s, 49 MB).
- adv2.py: ζ off-line-pair model, Laguerre, 2-point uniqueness on 72 sub-intervals, water-filling (logs/adv2.log; 336 s, 50 MB).
- adv3.py: Lloyd limits for k = 2, 3, 4 (logs/adv3.log; 65 s with the cache, 55 MB).
- DERIVATION.md: proofs of L1–L3, Z, F and the rate–distortion screen.
- Rerun: put Odlyzko's zeros1.gz, zeros3, zeros4, zeros5 (https://www-users.cse.umn.edu/~odlyzko/zeta_tables/) in DIR, then
  `uv run --no-project --with numpy --with scipy python3 quant_stats.py DIR` (same for baselines.py),
  `uv run --no-project --with mpmath --with numpy python3 -u adv2.py DIR CACHE`, `... adv3.py DIR CACHE ../ExploreHelio/results/dh_zeros.json`,
  `... adv.py DIR ../ExploreHelio/results/dh_zeros.json`; single-threaded (OMP/VECLIB = 1).

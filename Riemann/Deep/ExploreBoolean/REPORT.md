BARRIER: Boolean analysis cannot see M(x), the top Walsh coefficient of the squarefree knapsack down-set: level-k, hypercontractivity, Peres, KKL and Friedgut bound it by at best |D_x| ≈ 0.61x (topology by 0.20x), and a Beurling system with a planted zero 0.8 + 10i matches every measured Boolean statistic of the primes within 3.1 % (Hamming distance 5.4 %) while |M| ≍ x^{0.8} (node N10).

# ExploreBoolean — analysis of Boolean functions → N10 (explore33, 12 Nov 2026)

Brief: [.scratch/orca/explore33-brief.md](../../../.scratch/orca/explore33-brief.md). Derivations: [derivations.md](derivations.md) (D1–D8). Notation: n = π(x), cube {0,1}^{P_x}, D_x = {S : Σ_{p∈S} log p ≤ log x} ≅ squarefree m ≤ x, F(T) = 2^n \hat 1_{D_x}(T), Q(x) = |D_x|.

## Answer
1. **No new bound.** Best Boolean-type bound: |M(x)| ≤ Σβ_k(Δ_x) = #{odd squarefree m ∈ (x/2, x]} = (2/π²)x + O(x^θ) (Björner's Betti count; it is the identity M(x) = Σ_{x/2<m≤x, 2∤m} μ(m)). The four named tools give ≥ Q(x) at the top level. N10 record stands: |M(x)| ≤ 3x/(π·10^{10}) + 11.39√x (Chirre–Helfgott) and our 9.77·10^4 x^{0.9382}.
2. **Collapse.** F(T) = Σ_{m≤x} μ²(m) g_T(m) (g_T(p) = −1 on T, +1 off T): the spectrum is the set of sign-twisted squarefree sums, Parseval RMS √Q(x), so RH ⟺ the all-minus point is not an x^ε-outlier. The top coefficient is permutation-invariant, so it is Σ_k (−1)^k π_k^*(x) = P_x(−1) (explore1's prime-factor gas), i.e. Selberg–Delange at z = −1 = Perron's formula with Π_{p≤x}(1 − p^{−s}). Pakianathan–Winfree (1104.4324) already reduce random-weight knapsack complexes to that Perron integral.
3. **Barrier (rigorous, D5).** Planted Beurling primes with density (1 − 2u^{β−1}cos(γ log u))/log u: Π(1 − q^{−s}) = (s−1)e^{A}/((s−ρ)(s−ρ̄)) on σ > 1/2, so M_q ≠ O(x^{β−ε}), while every bounded-level Walsh coefficient has on σ > 1/2 only the pole s = 1, as for the primes (reviewer: D5(a) OK; D5(b) wording fixed). Fit at β = 0.8, γ = 10: M_q/x^{0.8} = 0.0201 cos(10 log x − 1.456), R² = 0.988, best grid (β, γ) = (0.79, 10.00).
4. **Adversarial.** The bounds are functions of (n, α) and LTF membership only: identical for log p, perturbed, Beurling and equal weights. log p is not generic: equal weights give |M| = Q − 2π(x) (0.81 Q at 10^8), i.i.d. or shifted perturbations give |M| ≍ x/(log x)^{O(1)} (heuristic for i.i.d.), and log(p+1) (Mellin E(s)/ζ(s), RH-equivalent) is 4× farther from D_x in Hamming distance (22 %) than the planted system (5.4 %).

## Transfer (one tool family, one node)
- Tools (exact statements, O'Donnell 2105.10386): level-k inequality W^{≤k}[f] ≤ ((2e/k) ln(1/α))^k α² for k ≤ 2 ln(1/α) (level-1 sharp 2α² ln(1/α)); Peres NS_δ ≤ √(1/⌊1/δ⌋) for LTFs, W^{≥1/δ} ≤ 3 NS_δ; KKL MaxInf ≥ Var·Ω(log n/n); Friedgut ε-close to an exp(O(I/ε))-junta; Bonami–Beckner ‖T_ρ f‖_2 ≤ ‖f‖_{1+ρ²}.
- Node N10, quantity M(x) = F(P_x) (D1). N4/N5 not pursued: σ(n)/n and the Báez-Duarte distance are not Walsh coefficients of a down-set, so the tools have no object there.
- Claim tested: some inequality valid for positive-weight LTFs, monotone or small-measure f gives |F(P_x)| < Q(x). Derived (D2): hypercontractivity |M| ≤ √(2^n Q) (minimum at ρ = 1: Parseval); level-k at k = n (valid from n = 20) |M| ≤ Q((2e/n) ln(1/α))^{n/2} ≈ Q·1.94^n; Peres |M| ≤ (√3/2)2^n n^{−1/4}; KKL, Friedgut: no single-coefficient bound (a junta on J ⊊ [n] has zero top coefficient, giving only |\hat f| ≤ √α); topology (2/π²)x. Biased measure q_p = 1/(p+1) makes f_x balanced (mean → e^{−γ}; Peres keeps its bound, KKL/Friedgut lose factors in the bias 1/(x+1)), but M(x) = Π_{p≤x}(−(p+1)/√p)·\tilde f_π(P_x), amplification e^{θ(x)/2}Π(1 + 1/p) (D3, checked exactly at x ≤ 100).

## Test 1 — smallest cases, exact spectra ([fwht.c](fwht.c), [run_fwht.sh](run_fwht.sh) → [fwht_out.txt](fwht_out.txt))
FWHT on int32, 2^24 coefficients at n = 24 (peak 68.6 MB, 0.4 s); every x ≤ 96.

| x | n | Q | M | \|M\|/√Q | share of T with \|F(T)\| ≥ \|M\| | Björner | hypercontr. | level-k (k = n) | Peres |
|---|---|---|---|---|---|---|---|---|---|
| 12 | 5 | 8 | −2 | 0.71 | 0.66 | 2 | 16 | n/a | 18.5 |
| 28 | 9 | 17 | −1 | 0.24 | 1.00 | 5 | 93 | n/a | 256 |
| 52 | 15 | 32 | −2 | 0.35 | 0.81 | 10 | 1.0·10³ | n/a | 1.4·10⁴ |
| 72 | 20 | 45 | −3 | 0.45 | 0.69 | 15 | 6.9·10³ | 1.1·10⁶ | 4.3·10⁵ |
| 96 | 24 | 60 | 2 | 0.26 | 0.86 | 20 | 3.2·10⁴ | 1.7·10⁷ | 6.6·10⁶ |

- For every x ≤ 96: |M| ≤ √Q, ≥ 39 % of all 2^n coefficients are at least |M|, F(P_x) = Σ(−1)^k A_k (symmetrization), the largest |F(T)|, T ≠ ∅, sits at level 1, the level-k inequality is at most 42 % saturated at any k (28.5 % at n = 24; maximum always at k = 1), level-1 ratio W¹/(2α² ln(1/α)) = 0.39–0.74.
- Same n, x = 89: equal weights M = 8 (|M|/√Q = 1.07), i.i.d. ±0.5 perturbations |M| ≤ 3, Poisson Beurling points |M| ≤ 18. At n ≤ 24 nothing separates the systems.

## Test 2 — at scale ([scale.c](scale.c) → [scale_1e8_out.txt](scale_1e8_out.txt), [hamming_out.txt](hamming_out.txt), [fit_out.txt](fit_out.txt))
DFS over D_q(x), x ≤ 10^8 (6·10^7 nodes; ≤ 24 s, peak ≤ 156 MB). Checks: M(10^k) = 2, −23, −48, 212, 1037, 1928 and Q(10^8) = 60792694 match the literature; Chow identity b_p = Q − 2c_p holds for every p in every system.

| system (weights) | max\|log q_i − log p_i\|, p_i > 10³ | \|D Δ D_x\|/\|D_x\| | N/x | I/(2nα) | Betti/N | mean ω | M(10^8) | max\|M\|/√x, 10^7–10^8 | growth exponent |
|---|---|---|---|---|---|---|---|---|---|
| primes | 0 | 0 | 0.6079 | 1.0000 | 0.333 | 2.888 | 1928 | 0.45 | 0.55 |
| p + 1 | 0.0010 | 0.224 | 0.4717 | 1.0000 | 0.500 | 2.692 | 846 | 0.12 | 0.58 |
| smooth li-points | 0.052 | 0.030 | 0.6238 | 1.0000 | 0.333 | 2.904 | 223 | 0.08 | 0.31 |
| Poisson Beurling | 0.069 | 0.213 | 0.5254 | 1.0000 | 0.504 | 2.783 | 682 | 1.77 | 0.56 |
| i.i.d. p·e^{U[−.5,.5]} | 0.090 | 0.315 | 0.7632 | 1.0000 | 0.201 | 3.024 | 6635 | 1.34 | 0.72 |
| shift p·e^{−0.25} | 0.25 | 1.13 | 1.3357 | 1.0000 | 0.221 | 3.268 | 36003 | 3.64 | 0.82 |
| planted ρ = 0.8 + 10i | 0.065 | 0.054 | 0.5893 | 1.0000 | 0.333 | 2.860 | 43967 | 5.04 | 0.82 |
| equal weights (D6) | — | — | 0.6079 | — | — | — | 4.93·10^7 | 4.9·10³ | 1 |

- Level-1 ratio W¹/(2N² ln(1/α)) = 0.7213 in every system (→ 1/(2 ln 2): each large prime is a near-dictator on a sparse knapsack set); MaxInf/α = 2 everywhere (p > x/2).
- Predictions (D7): shift lead term 585, 4113, 30316 at 10^6, 10^7, 10^8 vs 578, 5114, 36003 (exponent 1 + e^{0.25} = 2.284; second-order ≈ 0.73 log log x/log x); i.i.d. mean 10940 at 10^8 vs realization 6635 (below √x until ≈ 10^{10}).
- Planted: decade maxima of |M|/√x 0.74, 1.22, 2.43, 5.04 (×10^{0.3} per decade = x^{0.8−0.5}); amplitude 0.0200 / 0.0201 on [10^5, 10^{6.5}] / [10^{6.5}, 10^8].

## Test 3 — sampling Walsh coefficients ([sample_out.txt](sample_out.txt))
Random T, Pr[p ∈ T] = α (same seeds across systems).

| x | system | α | samples | median F(T) | RMS (√Q) | share with \|F(T)\| ≥ \|M\| |
|---|---|---|---|---|---|---|
| 10^6 | primes | 1/2 | 2000 | 4 | 761 (780) | 0.69 |
| 10^7 | primes | 1/2 | 200 | −345 | 2078 (2466) | 0.60 |
| 10^6 | planted | 1/2 | 2000 | −1 | 756 (768) | 0.74 |
| 10^7 | planted | 1/2 | 200 | −28 | 2345 (2428) | 0.025 |
| 10^6 | primes / planted / smooth | 0.75 | 400 | −4834 / −4901 / −4722 | — | 1.00 |
| 10^6 | primes / planted / smooth | 0.9 | 400 | −1468 / −1449 / −1416 | — | 1.00 |

Selberg–Delange level profile (D4): median ≈ G_T x(log x)^{−2α}/Γ(1 − 2α) (−5490 G_T at α = 0.75, −1544 G_T at α = 0.9), zero at α = 1/2 (random signs) and α = 1 (the top): Mertens' cancellation is the second zero of 1/Γ(1 − 2α).

## Adversarial verdict (iv)
- Equal weights: Hamming ball plus a generic level-2 segment is a positive-weight LTF of every size; top = Q − 2π(x): no cancellation, so no bound in (n, |D|, LTF⁺) beats |D| − 2n.
- Generic weights: small symmetric nondegenerate i.i.d. perturbations have K(1) = E e^{−δ} ∈ (1, 2) and E M_q ≍ x/(log x)^{1+K(1)}; exact-density ones (K(1) = 1) still give x/(log x)³, since K'(1) = E[X log X] > 0 for X = e^{−δ} (Jensen) (Selberg–Delange heuristic, D7; the deterministic shift is confirmed to 20 %): the power saving is lost. log p is special to all orders in 1/log t.
- Beurling: Poisson points (exact density 1/log t) and smooth li-points keep |M| ≲ 1.8√x; the planted system has |M| ≍ x^{0.8}.
- Every Boolean statistic measured moves ≤ 3.1 % between primes and planted (the 3.1 % is the density constant N/x; the threshold x' = 1.032x removes it), and all tool bounds are functions of (n, α) only. The tools carry no RH information.

## Where it breaks
1. One coefficient at level n, inequalities averaged over levels: at the top, level-k and hypercontractivity reduce to Parseval or worse.
2. Peres, KKL and Friedgut are additive (balanced-function) statements; α = 2^{−n(1−o(1))}.
3. Balancing by the 1/(p+1)-biased measure costs the factor e^{θ(x)/2}; the uniform measure minimises the amplification (Π 1/√(q(1−q)) ≥ 2^n).
4. Noise stability on |ρ| ≤ 1 is the pair series Σ_{m,m'} t^{d(m,m')} with t ≥ 0 (no signs); M² is the same series at t = −1 (ρ → ∞); certifying the leading coefficient from [−1, 1] data via Chebyshev loses 2^{n−1} (D8, heuristic).
5. Symmetrization: the top coefficient sees only the ω-profile, so any bound is a parity bound (Selberg's parity problem) = Perron with the truncated Euler product.
6. The zero lives only in the top levels: bounded-level coefficients carry it as a zero factor (no x^ρ term) and through P(2s), at x^{β/2} < √x (D5b).

## Review
Three `reviewer` passes on D2, D3, D5, D7, D8 (verdict only, text inline). Pass 1: D5(a), D6, D7 (p+1, shift) OK; fixed D2(b) (Hamming-ball ratio 1 − 2c), D3 (KKL/Friedgut lose bias factors), D5(b) wording, D7 i.i.d. assumptions, D8 (Chebyshev is route-specific, now labelled heuristic). Pass 2: MINOR ISSUES (D5(b) "only"; integer K(1)), fixed. Pass 3: D5(b) OK; D7 fallback replaced by Jensen (K'(1) > 0), primes ≤ u_0; MINOR ISSUES, fixed.

## Prior art and novelty
Björner 1101.5704 (Betti numbers; |M| ≤ (2/π²)x + O(x^θ) is the only cancellation it gives); Pakianathan–Winfree 1104.4324 (LogPrime quota complex, χ = 1 − M, random weights → Perron of Π(1 − \hat f_j)); explore1 ([ExploreStatMech](../ExploreStatMech/REPORT.md), item C: P_x(z)); R317 (AC0/Walsh Möbius orthogonality controls no prefix sum); planted Beurling zeros: Diamond–Montgomery–Vorhauer (Math. Ann. 2006), Révész 2202.01837. New here: the spectrum/level-profile reading, the tool-by-tool top-level bounds, the Chebyshev and biased-measure obstructions, the planted-zero invisibility of bounded-level Walsh data, and the perturbation table. arXiv math.NT and math.CO listings, Oct–Nov 2026 (3134 titles), keyword scan: no Boolean-analysis approach to M(x); 2611.08289 (level-1:3 inequalities, monotone f) is low-level only. No discovery claimed.

## Suggested NOGO line
- explore33 (Orca, exploration: analysis of Boolean functions → N10) BARRIER: M(x) is the top Walsh coefficient of the squarefree knapsack down-set (spectrum = all sign-twisted squarefree sums, RMS √Q); level-k and hypercontractivity reduce to Parseval there, Peres/KKL/Friedgut are additive at measure 2^{−n(1−o(1))}, topology gives (2/π²)x; the top coefficient sees only the ω-profile (explore1's P_x(−1), Perron); a Beurling system with a planted zero 0.8 + 10i keeps bounded-level coefficients pole-free on σ > 1/2 and every Boolean statistic within 3.1 % (Hamming 5.4 %) while |M| ≍ x^{0.8} (fit 0.79, 10.00); generic perturbations of log p lose the power saving, equal weights cancel nothing. [ExploreBoolean](Riemann/Deep/ExploreBoolean/REPORT.md)

## Files
| file | what | run |
|---|---|---|
| [derivations.md](derivations.md) | D1–D8 proofs and predictions | — |
| [fwht.c](fwht.c), [run_fwht.sh](run_fwht.sh), [fwht_out.txt](fwht_out.txt) | exact spectra, n ≤ 24, 4 weight systems | 10 s, 68.6 MB |
| [scale.c](scale.c) | 7 systems to 10^8: N, M, Chow/influence, Betti, ω-profile (mode 0), sampling (1), Hamming (3) | ≤ 24 s, ≤ 156 MB per run |
| [scale_1e8_out.txt](scale_1e8_out.txt), [sample_out.txt](sample_out.txt), [hamming_out.txt](hamming_out.txt) | outputs | — |
| [fit_curves.py](fit_curves.py), [fit_out.txt](fit_out.txt) | growth exponents, planted fit; curve files (deleted) regenerate with `scale <sys> 1e8 1 <p1> <p2> 0 0 0 curve_<sys>.txt` (sys 0–6; p1, p2 = 0 0 / 0.5 0 / 0.25 0 / 0.8 10) | 2 s, 45 MB |
| [biased_check.py](biased_check.py), [sd_pred.py](sd_pred.py) | D3 identity; D7 predictions | < 1 s |

All runs single-threaded, sequential, under perl alarm; cc -O2; Python via uv (numpy) or stdlib.

RH relevance: N10 (explicit Mertens bounds); also bars Boolean/complexity routes of N12 type (R317).

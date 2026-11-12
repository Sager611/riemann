NO IDEA: Fourier-uniqueness theory reaches the tree only as exactness, and the Báez-Duarte system is already exact: the unitary Mellin multiplier (s−1)/s makes it lower-triangular on the Gauss cylinders (1/(j+1), 1/j] with Möbius coordinates (minimality: Vasyunin 1995; squarefrees suffice: Bagchi 2006), so it has no Heisenberg-type density margin (every dilation set missing a squarefree m fails unconditionally, floor (1 − log²2)/2 = 0.2598 at m = 1, ≈ 1/(2m²) at prime m, 0.02826 for prime dilations, Arb), the Hedenmalm–Montes-Rodríguez Gauss-map iteration never acts on it (one digit), its weighted form is Mayer's operator (L¹-bounded only for Re s ≥ 1, explore3), and Davenport–Heilbronn has the same floors (node N5).

# Exploration explore37: Fourier uniqueness → RH tree (12 Nov 2026, 05:59–06:55 UTC)

Field: Heisenberg uniqueness pairs (Hedenmalm–Montes-Rodríguez, Ann. Math. 2011; Canto-Martín–Hedenmalm–Montes-Rodríguez, JEMS 2014; Hedenmalm–Montes-Rodríguez 1503.04038; Radchenko–Ramos 2410.04557; Wan 2607.17159), Fourier interpolation and uniqueness pairs (Radchenko–Viazovska; Bondarenko–Radchenko–Seip 2005.02996 = BRS; Kulikov–Nazarov–Sodin 2306.14013; Lysen 2509.17600), Hardy and Beurling–Hörmander.
Sources: arXiv abstract pages; TeX of 2005.02996, 2203.05030, 2607.12084, math/0011254, math/0607733, 0812.1689, 1809.09577 (rg only); arxiv.org site search (the export API returned 429); `.scratch/lst.py` listings; Zenodo API. Nodes: N5 (main), N1, N13.

## 1. Bridges measured

| bridge | tool | node | outcome |
|---|---|---|---|
| B1 | HMR's engine: transfer operator of the Gauss-type map at the critical density | N5, N13 | edge annihilator = Gauss density; Jacobian weights give Mayer's L_s², L¹-bounded only for Re s ≥ 1; ζ only at s = ρ/2 (explore3) |
| B2 | exactness at criticality: critical interpolation pairs are complete and minimal, HMR's threshold has no slack | N5 | **chosen**: is Báez-Duarte's dilation set thinnable? No: exact on the squarefrees (Theorem A), known minimality |
| B3 | criticality of the zeta-zero pair (BRS) | N1 | known: BRS attain Kulikov's density bound; it counts zeros, not their real parts |
| B4 | Hardy, Beurling–Hörmander | N1 | vacuous on windows; the applicable uncertainty bounds relax the exact window form (explore32) |

## 2. Statements

**S0 (exact reformulation).** U = Mellin multiplier (s−1)/s is unitary on L²(0,∞) (|s − 1| = |s| on Re s = ½).
- U χ_(0,1] = e, e(x) = 1 + log x on (0,1], 0 beyond (Mellin (s−1)/s²).
- U ρ_{1/k} = −f(k·), f(x) = Σ_{i≤1/x}(1 + log ix) = ⌊1/x⌋(1 + log x) + log⌊1/x⌋! (Mellin ζ(s)(s−1)/s²; f = ½ log(2π/x) − {1/x} + O(x)).
- Hence Báez-Duarte's V1 distance (L²(0,∞), as in ../NBRate) is d_N(S)² = min_a ‖e − Σ_{k∈S, k≤N} a_k f(k·)‖²_{L²(0,1)}, with supp f(k·) = (0, 1/k].
- On the Gauss cylinder I_j = (1/(j+1), 1/j] (first continued-fraction digit j): f(kx) = n(1 + log k + log x) + log n!, n = ⌊j/k⌋. Every function of the problem is affine in log x on every cylinder; only the first digit enters.
- f is the co-Poisson transform Σ_n e(nx) − (1/x)∫e of e, and ∫_0^1 e = 0. So S0 is Báez-Duarte's general strong criterion (math/0505453) with his f := e. e is outside his smoothness hypothesis; here the criterion follows from unitarity. The V2 step-function analogue is his 2000 Lemma (math/0011254).

**S1 (triangularity).** On (1/(M+1), 1]: e = Σ_{k≤M} μ(k) f(k·) exactly (Dirichlet convolution; n ≤ 1/x ≤ M). The f(k·), k ≤ M, are independent there: f(k·) = 1 + log kx on I_k and 0 on (1/k, 1].

**Theorem A (unconditional).** For S ⊂ ℕ: d_N(S) → 0 ⟺ RH and S contains every squarefree integer.
- If a squarefree m ∉ S, then for all N and all M ≥ m: d_N(S)² ≥ η_M(S)² := dist²_{L²(1/(M+1),1]}(e, span{f(k·): k ∈ S, k ≤ M}) > 0, and η_M is non-decreasing in M.
- Proof: restrict the error to (1/(M+1), 1]; functions with k > M vanish there; by S1 the coordinates of e are μ(k), and the one at m is nonzero.
- Converse: Bagchi's remarks (math/0607733, the remarks labelled First and its sequel: RH implies the squarefree-thinned span suffices); necessity because d_N(S) ≥ d_N(ℕ).
- Reviewer (verdict only, claims inline): S0, S1, Theorem A CORRECT; overall ROUTINE given Vasyunin, Báez-Duarte and Bagchi; the prime law below is numerical; float numbers are not certificates (the Arb ones are).

**S2 (coefficient criterion, routine).** Let c*(N) be the optimal V1 coefficients for S = ℕ. Then RH ⟺ −c*_k(N) → μ(k) for every k.
- The coordinate of f(k·) on (1/(k+1), 1] is L²-continuous on the closed span and equals −c*_k(N) on the N-th optimum.
- Optima converge to the projection Pe; its coordinates are all μ(k) iff Pe = e iff RH.
- Báez-Duarte 2000 proved the pointwise analogue for V2.

**S3 (HMR edge, smallest case).** For the branch Γ₊ at αβ = 4, i.e. e^{2πimt}, e^{2πin/t} on t > 0, the annihilator is spanned by f_G = 1/(1+t) on (0,1) and −1/(t(1+t)) on (1,∞).
- Σ_{j≥0} f_G(t+j) telescopes to 0, and J f_G = −f_G for J f(u) = u^{−2} f(1/u).
- Every annihilator has h = f|_(0,1) = L²h, with L the Gauss transfer operator, so h ∝ 1/(1+x) by exactness of the Gauss map. This matches 1503.04038 (dense iff αβ < 4) and CHM-R (one-dimensional at αβ = 4).
- Weighting the inverse family by t^{2s−2} turns the equation into h = L_s²h (Mayer). Since ∫_0^1 |L_s h| ≤ ∫_0^1 u^{2Re s−2}|h|, L_s is L¹-bounded iff Re s ≥ 1.
- ζ's zeros are eigenvalue-one points of the continued L_s at Re s = ¼ (Efrat; explore3: |det(1 − L_b)| = 5.5·10⁻¹⁶ at b = ρ₁/2). No measure or L¹ uniqueness statement reaches them.
- The Báez-Duarte system uses only the first digit (S0), so no iterate of the Gauss map acts on it.

**S4 (B3, B4; argued, known).**
- BRS's pair (γ, ±(log n)/4π) is exactly critical: their §2 attains Kulikov's bound 2N(T) + e^{4πW} ≥ 4WT − C log^{2+η}(4WT) at W = log(T/2π)/4π.
- N(T) counts all zeros whatever their real parts, so criticality never sees real parts. For a Davenport–Heilbronn analogue only the conductor changes N(T), by (T/2π) log 5.
- Hardy's and Beurling–Hörmander's hypotheses fail for every nonzero window function: its transform has exponential type, so even e^{−ε|t|} decay on ℝ forces 0.
- The uncertainty bounds that do apply relax the exact window Galerkin problem (explore32: Arch > 0 only for a < 0.0856; Veltri's exact odd threshold is 0.3716).

## 3. Tests (code/, logs/)

**T1 sanity** (u1, M = 1000):
- ‖e − Σ_{k≤1000} μ(k) f(k·)‖ on (1/1001, 1] = 1.8·10⁻¹⁴.
- ⟨e, f(k·)⟩ → −b_k and ‖f‖² → G₁₁ = log 2π − γ, with tails O(log²M/M) (0.031 and 0.019 at M = 1000).
- ‖e‖² on (1/1001, 1] matches 1 − ∫_0^{1/1001}(1 + log x)² to 10 digits.

**T2 smallest case.** η₁(ℕ∖{1})² = ∫_{1/2}^1 (1 + log x)² dx = (1 − log²2)/2 = 0.2597734930 (Arb ±10⁻¹³).

**T3 single removals** (S = ℕ∖{m}, M = m; float; u1, u5):

| m | 1 | 2 | 3 | 5 | 6 | 7 | 10 | 30 | 101 | 199 | 997 | 1001 | 1999 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| η² | 0.2598 | 0.07087 | 0.04968 | 0.01918 | 0.01086 | 0.009986 | 0.003786 | 3.771e−4 | 4.901e−5 | 1.263e−5 | 5.030e−7 | 4.822e−7 | 1.251e−7 |
| m²η² | 0.260 | 0.283 | 0.447 | 0.480 | 0.391 | 0.489 | 0.379 | 0.339 | 0.49995 | 0.49999 | 0.50000 | 0.483 | 0.50000 |

- Primes: m²η² → ½ (numerical).
- Squarefree composites: 0.34–0.50 (0.339 at 30, 0.483 at 1001, 0.498 at 437); fit 0.436 m^{−1.995} on [20, 300].
- Non-squarefree m: η = 0 (≤ 10⁻²⁸).
- Arb at m = 101: 4.900966604·10⁻⁵.

**T4 families.** Floors are rigorous (Arb, u6/u6b) or float at M = 2000 (u5). The actual d_N(S)² come from the double Vasyunin Gram (u3, Gram code copied from ../NBRate/code/s1_double.py).

| S | floor, Arb (M) | floor, float M = 2000 | d_N(S)², N = 20 | N = 160 | N = 640 |
|---|---|---|---|---|---|
| ℕ∖{1} | 0.2644684568 (60) | 0.2646213 | 0.267608 | 0.266177 | 0.265828 |
| ℕ∖{2} | 0.0708870664 (120) | 0.0708874 | 0.082669 | 0.077338 | 0.075873 |
| odd | 0.0979309484 (240) | 0.0979474 | 0.109013 | 0.103409 | 0.101968 |
| {1} ∪ primes | 0.0282590219 (1000) | 0.0283334 | 0.028507 | 0.028416 | 0.028413 |
| random ½ #0 / #1 / #2 | — | 0.16580 / 0.14870 / 0.07582 | 0.16873 / 0.15425 / 0.07939 | 0.16712 / 0.15137 / 0.07739 | 0.16682 / 0.15067 / 0.07702 |
| squarefree | 0 | 0 | 0.017591 | 0.009659 | 0.0071818 |
| ℕ (reference) | 0 | 0 | 0.016504 | 0.009346 | 0.0070171 |

- Floors sit 0.3–7 % below the N = 640 values.
- Prime dilations: η_M² = η_∞² − 0.1488/M fits M = 1000, 1500, 2000 to 10⁻⁷ and gives η_∞² = 0.028408, against a stall near 0.028412. The thinned stall is a bulk quantity.
- For ℕ every floor is 0, so the whole difficulty sits near 0: the Möbius choice is exact on (1/(N+1), 1] but diverges in L² (Báez-Duarte).
- explore31's random-subset stalls (d² ≈ 0.1 to N = 640) are such floors.
- lim_N d_N(ℕ∖{1})² exists unconditionally and lies in [0.26462, 0.26583].
  - Under RH it equals lim δ²_{1,N}, δ²_{1,N} = 1/(G_N⁻¹)₁₁ = dist²(ρ₁, span{ρ_{1/k}: 2 ≤ k ≤ N}), which falls from 0.27656 (N = 20) to 0.26977 (N = 640); so RH forces it to fall by at least another 0.0039.
  - Least-squares fits L + B/log N to either sequence (N = 20–640 or 80–640) give L = 0.2641–0.2643, below the floor, so both are pre-asymptotic.

**T5 coefficients** (u3; −c*_k(N)/μ(k)):

| N | k = 1 | 2 | 3 | 5 | 6 | 10 | max\|c*\| over k = 4, 8, 9 |
|---|---|---|---|---|---|---|---|
| 20 | 0.9529 | 0.9254 | 0.9094 | 0.7679 | 0.5541 | 0.3809 | 0.130 |
| 160 | 0.9726 | 0.9541 | 0.9477 | 0.8488 | 0.7251 | 0.6314 | 0.085 |
| 640 | 0.9795 | 0.9663 | 0.9612 | 0.8923 | 0.7945 | 0.7392 | 0.058 |

1 + c*₁(N) ≈ 0.135/log N on [20, 640]. S2 makes these tend to 1 under RH, but the approach is logarithmic and cannot discriminate.

## 4. Adversarial

- **Davenport–Heilbronn (u4).** L(s)/s ↔ A(⌊1/x⌋) on (0,1], with A the partial sums of a = 1, κ, −κ, −1, 0 (mod 5), κ = 0.2840790438.
  - The triangular structure is the same: χ = Σ b_k A(⌊1/(k·)⌋) on (1/(M+1), 1], where b is the Dirichlet inverse of a (exactness 2.5·10⁻¹⁶).
  - Essential set = supp b = the n ≤ 60 coprime to 5. It contains the non-squarefree 4, 8, 9, 12, 16, … and misses the squarefree 5, 10, 15, 30, 35, 55.
  - Floors: ℕ∖{1}: ½; ℕ∖{4}: 0.03225; squarefrees: 0.0499; {1} ∪ primes: 0.0783.
  - ℕ∖{p}: p²η² ≈ |b_p|²/2 = |a_p|²/2, i.e. ½ for p ≡ ±1 mod 5 (0.4959 at 11, 0.50000 at 101 and 199) and κ²/2 = 0.040 for p ≡ ±2 (0.0400 at 3, 0.0403 at 7, 0.0347 at 2). For ζ, |μ(p)|²/2 = ½.
- Floors and essential sets are fixed by the Dirichlet coefficients and the Euler product (squarefree support of μ), not by zeros. Theorem A is the same theorem for DH, whose own d_N stays positive (explore13).
- A moved-pair model has no Dirichlet series, hence no Báez-Duarte system.
- So no RH information.

## 5. Verdict and novelty

- **Known:**
  - minimality of {ρ_{1/k}} (Vasyunin 1995, as cited in 2203.05030 l. 130);
  - a complete biorthogonal sequence with coefficients μ(k)/k in Noor's H² model (Calderaro–Manzur–Noor–Santos 2203.05030, Thm mainthm);
  - the pointwise coefficient lemma on the cylinders (Báez-Duarte math/0011254, Lemma coeffsconverge);
  - squarefree sufficiency (Bagchi math/0607733);
  - Möbius uniqueness for pointwise convergence (Pyvovarov 2607.12084);
  - BRS criticality (2005.02996 §2).
  Theorem A and S2 are immediate corollaries (reviewer: ROUTINE).
- **New here (as far as found).**
  - The floors η_M(S), Arb-certified for every N: d_N² ≥ 0.26446 without k = 1 (0.25977 = (1 − log²2)/2 from (1/2, 1] alone), ≥ 0.070887 without k = 2, ≥ 0.097930 for odd k, ≥ 0.028259 for prime dilations.
  - The prime law m²η² → ½ (numerical).
  - The bulk/layer split (T4).
  These are records and corollaries only; no node quantity moves.
- **Scans.**
  - lst.py over math.NT/CA/CV/FA/SP/DS/PR/CO, Oct–Nov 2026: 12 title hits, none on thinning or HUP ↔ NB.
  - Zenodo "Nyman-Beurling" OR "Baez-Duarte" OR "Heisenberg uniqueness": 59 records; 24054105, 23961688, 23961279, 6982337 read, none relevant.
  - arxiv.org searches: "Heisenberg uniqueness pairs", "Fourier uniqueness pairs", "Baez-Duarte", "Nyman-Beurling", "Gauss-type maps Klein-Gordon".
  - No Evidence-API call.

## Suggested NOGO line

- explore37 (Orca, exploration: Fourier uniqueness — Hedenmalm–Montes-Rodríguez Heisenberg uniqueness pairs and their Gauss-type transfer operators, Radchenko–Viazovska / Bondarenko–Radchenko–Seip interpolation, Kulikov–Nazarov–Sodin pairs, Hardy and Beurling–Hörmander → N5, N1, N13) NO IDEA: the field transfers only as exactness, and the Báez-Duarte system is already exact: the unitary Mellin multiplier (s−1)/s turns V1 into approximating 1 + log x on (0,1] by f(kx), f(x) = Σ_{i≤1/x}(1 + log ix), lower-triangular on the Gauss cylinders (1/(j+1), 1/j] with Möbius coordinates (minimality Vasyunin 1995, squarefrees suffice Bagchi 2006), so d_N(S) → 0 iff RH and S ⊇ squarefrees, and every S missing a squarefree m has an unconditional floor (Arb: 0.2645 without 1, of which (1 − log²2)/2 = 0.2598 from (1/2, 1]; 0.07089 without 2; 0.02826 for prime dilations; m²η² → ½ at primes, numerical) within 0.3–7 % of the N = 640 stalls (explore31's random stalls); RH ⟺ the optimal coefficients tend to −μ(k) (routine; 0.98 at k = 1, N = 640); the HMR Gauss-map iteration never acts (one digit), its Jacobian-weighted edge equation is Mayer's L_s² (L¹-bounded only for Re s ≥ 1; explore3); BRS's zero pair is exactly critical by counting all zeros; Hardy and Beurling–Hörmander are vacuous on windows (explore32); Davenport–Heilbronn has the same theorem with essential set supp(1/L) (the n ≤ 60 coprime to 5, 4 and 9 included) [ExploreUniq](Riemann/Deep/ExploreUniq/REPORT.md)

RH relevance: N5 (structure of the Báez-Duarte criterion: no thinning, no density margin); N1 and N13 untouched.

## Files

- **code/** (run with `uv run --no-project --with mpmath --with numpy --with scipy --with python-flint python3 <file> <args>`; single-threaded, under `perl -e 'alarm N'`):
  - u1_exact.py `1000`: S0/S1 checks, single removals, families (peak 132 MB);
  - u3_gram.py `640 1000`: actual d_N(S)², δ_{m,N}, c*_k(N); the Vasyunin Gram is copied from ../NBRate (143 MB);
  - u4_dh.py `1000`: Davenport–Heilbronn analogue (55 MB);
  - u5_floors.py `2000`: floors to M = 2000, chunked dsyrk (256 MB);
  - u6_arb.py, u6b_arb.py: Arb floors (97 MB, 191 MB).
- **logs/**: their outputs. Working notes: .scratch/explore37/STATE.md.

## Caveats

- Float numbers (T3, T4 float columns, T5) are double precision. The u3 Gram has condition number ≤ 1.6·10⁶ (../NBRate). Only the Arb floors are certificates.
- The prime law m²η² → ½ and the extrapolated limits are numerical.
- Memory: the first u5 version materialised the full matrix and fancy-index copies and peaked at 439 MB, over the 300 MB cap (a 7 s run). It was replaced by the chunked version, 256 MB with identical numbers. All other runs stayed ≤ 191 MB.

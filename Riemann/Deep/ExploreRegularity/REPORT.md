NO IDEA: Szemerédi/Frieze–Kannan regularity and dense models are blind at gap scale (every graphon limit of the gap graph is 0; a dense model misstates it by 1/(2s) + √Σ²(s)/s ≥ 0.5 at s ≤ 1) or re-encode band-limited correlations (DPP cap 1/(1+θ), LP duality); the field's one live move, majorant transference applied to Theorem U's θ < 1/2 cap, re-derives dpcbeyond's same-day Theorem U′ (λ > 7/11): independent corroboration, nothing new (node N8)

# explore43: Szemerédi regularity and the dense model theorem → N8 (Orca exploration, 12 Nov 2026, 14:21–15:25 UTC)

Derivation: [DERIVATION.md](DERIVATION.md). Code: [code/](code/). Logs: [logs/](logs/).

## 1. Transfer (bridges measured before choosing)
| bridge | tool | measurement | verdict |
|---|---|---|---|
| N8: gap graph G_s of zeros (edge if unfolded distance ≤ s), the brief's bridge | Frieze–Kannan weak regularity; regularity + counting lemma | G_s has ≈ 2sN edges: dense cut norm O(s/N) → 0. A dense model (cell counts at resolution s) misstates G_s by the atomic term 1/(2s) and the fluctuation term √Σ²(s)/s (ζ: 1.04, 0.57, 0.038 at s = 0.5, 1, 16) | blind below s ≈ 1/ε |
| N8: band-limited counts, simple/distinct zeros | counting lemma; Green–Tao/RTTV dense model; Gowers' Hahn–Banach decomposition | counts of a kernel φ, φ^ ⊂ [−θ, θ], are band-limited n-level correlations, shared by the a-discrete sine DPP (a = 1/(1+θ), no gap below a); the dense model at scale 1/θ has ‖g‖₂² − ‖g‖₁ = ∫φ^ min(\|α\|, 1): Montgomery's F; Hahn–Banach decompositions are LP duality | re-encoding: D(λ) capped at 7/11 (θ = 4/7, explore34); multiplicity constant capped by CCLM's 1.3274992 (optimal for R̂ ⊂ [−1, 1]) |
| N8: Theorem U's cap θ < 1/2 (needs ∫\|R_T\|^{2+δ} ≪ T at length T^θ > T^{1/2}) | transference (d): dominate an unbounded weight pointwise by a majorant with known moments | explicit formula: DPC's rough prime sum is an exact smoothed zero sum (error 4.1·10⁻⁹ on actual zeros); pointwise slit bound ⇒ only 2θ < κ remains | chosen; result = dpcbeyond's Theorem U′ |
| same input | Green's restriction theorem (ℓ² → L^q, q > 2, via majorant Fourier decay) | it removes a majorant's sup-loss (log N in Green's primes); here the loss is large values at 1 ≤ σ ≤ T^{θ−1/2}, where the Halász–Montgomery–Bourgain level-set argument needs σ ≥ T^{1/4} | no transfer |
| N1: Frieze–Kannan on Weil's form | – | R442 (no bounded graphon preserves Weil's form); cut norm does not control the operator norm, and the c = 47 margin is 2^{−803} | not chosen |
| N10/N11: dense model for Λ | Green–Tao majorant | transference error is the majorant's pseudorandomness, o(1) (≥ 1/log x); RH needs x^{−1/2} on interval tests | not chosen |

## 2. Statements (proofs and sketches: DERIVATION.md)
- Proposition S (§1.1). For any simple point set of unit density, a dense model of G_s at resolution s counts the N diagonal pairs (x, x): relative error exactly 1/(2s), plus cell fluctuations √Σ²(s)/s. Gap statements (s = λ < 1) are below every regularity scale.
- Proposition C (§1.1). Any counting-lemma or dense-model count of band-limited kernels is a band-limited n-level correlation sum; hence (explore34, Lagarias–Rodgers 1907.03391) no such argument proves D(λ) > 0 for λ < 1/(1+θ).
- Lemma A, Lemma B, Theorem U′ (§2): the unconditional form of DPC's Lemma 5 (prime-sum route; off-line zeros at complex argument, coupled cutoff f(Σ log q_i/L) caps the joint weight at T^{θ max d}), the pointwise slit bound T^{2θd} for the hairpins, and the claim D(λ) > 0 for every λ > 1/(1 + min(1, κ/2)): 7/11 with Conrey's κ = 8/7, 1/2 under the density hypothesis near the line, 2/3 with Jutila's κ = 1. This is dpcbeyond's Lemma Z, Lemma M, Steps 6e′, 7′ and Theorem U′ (Riemann/Deep/DPCBeyond, commit 331db28, DISCOVERIES and RECORDS), found independently here before reading it (§3).
- Observation H (§3). Global q-th moments of R_T (available for qθ < κ) with any Hölder split on the slits give the sparsity exponent E = 1 − 2/q, hence only θ < κ min(1/q, 1 − 2/q) ≤ κ/3 = 8/21 < 1/2: worse than Montgomery–Vaughan. The pointwise slit bound is necessary for any gain.

## 3. Tests
a) Smallest case, gap graph (code/partA.py, logs/partA.log; ζ: 22,462 zeros 100 < t < 2·10⁴; GUE tridiagonal bulk 9,993; Poisson 10⁴; DH on-line 304 below 420; lattice):
| s | atomic 1/(2s) | ζ | GUE | Poisson | DH | lattice |
|---|---|---|---|---|---|---|
| 0.5 | 1.000 | 1.044 | 1.061 | 1.412 | 1.034 | 1.000 |
| 1 | 0.500 | 0.567 | 0.588 | 0.986 | 0.523 | 0 |
| 4 | 0.125 | 0.155 | 0.175 | 0.506 | 0.190 | 0 |
| 16 | 0.031 | 0.038 | 0.049 | 0.257 | 0.078 | 0 |
(fluctuation column √Σ²(s)/s.) Σ²(16): ζ 0.38 (Berry saturation at T ~ 10⁴), GUE 0.62, DH 1.57, Poisson 17.1.
b) Majorant identity (code/t1.py a, logs/t1a.log): S_Λ(t) = −Σ_ρ F(L(γ − t)/2π) on the actual zeros, T = 10⁴, ρ supported in (0.05, 0.6)L and (0.05, 1.0)L (67 and 1,229 prime powers), 40 random t: max error 4.1·10⁻⁹, |S_Λ| ≈ 1.
c) Moments of DPC's rough R_T (code/t1.py b2, logs/t1b2.log; η = 0.05, 20,000 random t): normalised M₄ = 1.580 against the Montgomery–Vaughan diagonal 1.572 at T = 10⁴, θ = 0.6; 1.543 vs 1.531 at T = 10⁶, θ = 0.55; 1.578 (T = 10⁸, θ = 0.55) and 1.738 (θ = 0.6), while MV's off-diagonal factor 1 + T^{2θ−1} is 5.0–40.8. The T^{2θ} of Montgomery–Vaughan never appears (RH range: this tests only the on-line part).
d) At scale, Lemma Z's constant (code/t2.py, logs/t2.log, t2_eta01.log; θ = 0.55): max |R_T|/(1 + On)^k, On = Σ_ρ min(1, (L|γ − t|)^{−A}):
| height | T | η = 0.2 (k = 2), A = 2 / 4 | η = 0.1 (k = 5), A = 2 / 4 |
|---|---|---|---|
| t ∈ [10⁴, 1.99·10⁴], 495,000 points | 10⁴ | 0.987 / 1.208 | 1.123 / 1.852 |
| 34 Gram windows at t ≈ 10⁸, 4,080 points | 6.7·10⁷ | 0.905 / 1.134 | 0.950 / 1.737 |
| 34 Gram windows at t ≈ 9.8·10¹⁰, 4,080 points (125,854 terms) | 6.5·10¹⁰ | 0.898 / 1.163 | 0.873 / 1.916 |
No drift over seven decades of T (R454: 1.47, 1.53 at T = 2000, 6000). Weak test: |R_T| ≤ 2.5 here, and the lemma's content (off-line weights) cannot occur in the RH-verified range; R454's planted-mechanism test covers it.

## 4. Adversarial
- Propositions S and C hold for every point set: ζ ≈ GUE ≠ Poisson in Σ², DH (low height) less rigid; none of it is a node quantity or RH information.
- Theorem U′'s conclusion is RH-blind: D(7/11) = 0.18–0.21 for ζ, 0.208 for the sine process, 0.118 for DH's on-line zeros (explore34). Its inputs are ζ-specific: the Euler product (main term at Re s = 3/2) and a uniform density slope κ > 0; DH has N(σ, T) ≫ T for 1/2 < σ < 1 (Voronin), so κ_DH = 0 and the method gives nothing there. Off-line zeros only cost weight, so it carries no RH information, like every N8 gap result.

## 5. Novelty
- Field applied to ζ: arXiv API (graphon/regularity lemma/dense model/transference principle with zeta zeros): no hits; math.NT/CO/CA listings 2026-10 and 2026-11: dense-model papers (2610.39176, simplified dense model theorem for random structures; 2610.37119, fractal transference) are off-topic; gap papers are the known 2611.14486, 2610.28527, 2610.26038, 2610.22927, 2610.31744. Our files: R044 (dense model, CIRCULAR), R442 (graphon on N1), explore34 (DPP collapse).
- Theorem U′: already ours (dpcbeyond, 12 Nov 2026; R454 falsification passed; referees CORRECT/GAP-repairable). Here only an independent re-derivation of the same exponents, Observation H, and the at-scale numbers 3b–3d (new data for the DISCOVERY's 'own verification', not a new result).

## 6. Process
One process at a time, single-threaded, every run under memcap 290. Peaks: partA 26 MB, t1a 21 MB, t1b 257 MB, t1b2 207 MB, t2 122 MB, t2 (η = 0.1) 167 MB. One run (dict-based fourth-moment diagonal, T = 10⁶, θ = 0.6) reached 379 MB and was killed by the watchdog; replaced by a capped numpy diagonal (logs/t1b2.note). DPC's TeX (2610.28527v2) downloaded, read (Lemma 5, Lemma 6), deleted.

## Suggested NOGO line
- explore43 (Orca, exploration: Szemerédi regularity and the dense model theorem — Frieze–Kannan, counting lemma, Green–Tao/RTTV/Gowers dense models, transference → N8) NO IDEA: the gap graph of zeros is sparse (graphon limits 0) and a dense model misstates it by the atomic term 1/(2s) and the fluctuation √Σ²(s)/s (ζ: 0.57 at s = 1); band-limited counts and Hahn–Banach decompositions re-encode the pair LP (DPP cap 1/(1+θ), CCLM 1.3274992); Green's restriction removes only a majorant's log-loss; majorant transference on Theorem U's θ < 1/2 cap re-derives dpcbeyond's Theorem U′ (7/11) independently (explicit-formula majorant exact to 4·10⁻⁹; Lemma Z constant ≤ 1.92 up to T = 6.5·10¹⁰; global moments + sparse Hölder alone give only θ < κ/3 = 8/21) [ExploreRegularity](Riemann/Deep/ExploreRegularity/REPORT.md)

RH relevance: N8 (small gaps: corroborates the record λ > 7/11 and its mechanism; no new bound from the field); no RH information.

## Files
REPORT.md, DERIVATION.md; code/: partA.py (gap-graph dense-model errors, number variance), t1.py (a: explicit-formula identity on zeros; b, b2: moments of DPC's rough R_T), t2.py (Lemma Z constant on zeros; run from code/, arg η); logs/: partA.log, t1a.log, t1b.log (η = 0.15, near-constant R_T, uninformative), t1b2.log, t1b2.note, t2.log (η = 0.2), t2_eta01.log (η = 0.1). Inputs: Riemann/Deep/ExploreHelio/results/{zeta_zeros.npy, dh_zeros.json}, Riemann/Deep/ExploreDPP/logs/zeros_1e8.txt, zeros_1e11.txt.

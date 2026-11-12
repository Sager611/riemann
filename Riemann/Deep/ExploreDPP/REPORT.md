BARRIER: determinantal processes give no new lower bound on N8 gap quantities, but as LP duals they cap the small-gap method: Lamzouri's Thm Main2 (any θ < 4/7) cannot prove D(λ) > 0 for any λ ≤ 0.775 (Arb-checked dual; the discrete-sine DPP alone gives 7/11), so the record 0.7851 is within 0.0101 of everything the method can give, and float LPs put its floor near 0.784 (the dual floor grows like 0.784 − 0.21/X with the dual's support X); its RH version (θ = 1) cannot reach λ ≤ 0.595, against CGdL's 0.6039 (node N8)

# ExploreDPP: determinantal point processes → N8 small gaps (explore34, 12 Nov 2026)

## 1. Transfer
- Tool: HKPV (counts of a DPP are sums of independent Bernoullis), applied to the discrete sine process, the DPP on aZ with kernel
  a·sinc(x − y) (Lagarias–Rodgers 1905.12123 use a = 1/2 for the Alternative Hypothesis, AH). Randomly shifted, it is a stationary
  simple point process with Montgomery's structure factor S(α) = |α| on |α| < 1/a − 1 and no gap below a.
- Node N8, quantity D(λ) = lower density of consecutive normalised gaps ≤ λ. Record: D(λ) > 0 for λ ≥ 0.7851, unconditional
  (SmallGap: our SDP certificate for Lamzouri 2611.14486 Thm Main2; its caveat: "we have no dual bound, so the method's floor is
  unknown"). Thm Main2 is an LP: find even f ≥ 0, f^(0) = 1, f^ ≤ 0 on |u| > λ with J_θ(f) = f(0) + 2∫_0^θ f(α)α dα > 1, where
  θ < 4/7 comes from Conrey's density (θ = 1 under RH: CGdL/BGMM 0.6039).
- N13 not chosen: the CCM rungs δ_n are the HKPV Bernoulli parameters of the sine process on the window, but the ratios ε(j)/δ_4j
  are set by the primes ≤ λ² in QW_λ. HKPV knows no primes, so it only re-encodes prolate asymptotics.

## 2. Statements (proofs: DERIVATION.md)
- Collapse (§4): a lower bound on D(λ), simple-zero or gap proportions that uses the zeros only through band-limited
  (Montgomery–Rudnick–Sarnak) correlations also holds for the discrete-sine DPPs, which have every HKPV property and no gaps below
  a. HKPV, negative association, strong Rayleigh and Soshnikov add nothing on the bound side; the band-limited bounds of this type
  are LP bounds (Beurling–Selberg, CCLM, CGdL, Thm Main2). Arguments with other inputs (mollifiers, moments, resonance) are not covered.
- Lemma 1 (weak duality, §1): a measure σ = Σ_{k≠0} h(1 + η_|k|) δ_{kh} ≥ 0 with no mass in (−λ, λ) whose structure factor
  P(α) = (1 − h) + 2hΣ η_k cos(2παkh) satisfies P ≥ |α|·1[|α| ≤ θ] on [0, 1/(2h)], θ ≤ 1/(2h), kills Thm Main2 at (θ, λ):
  J_θ(f) ≤ f^(0) for every admissible f.
- Theorem A (DPP dual, §2): the 1/(1+θ)-discrete sine gives P = 1 − Σ_m Δ(α − m(1+θ)), equal to |α| on [−θ, θ] and ≥ 0 elsewhere,
  so no λ ≤ 1/(1+θ) is reachable for 0 < θ ≤ 1: 1/2 at θ = 1 (the AH obstruction), 7/11 = 0.6364 for θ < 4/7.
- Theorem B (certified lattice duals, §3): (i) for every θ ≤ 4/7 and λ ≤ 0.775 no admissible f has J_θ(f) > f^(0), so Thm Main2
  and SmallGap's Lemma T cannot reach any λ ≤ 0.775 with any θ < 4/7; (ii) θ = 1: no λ ≤ 0.595.

## 3. Tests
Smallest case: θ = 1, h = 1/2, η_k = −sinc²(k/2): P(α) = 1/2 − (4/π²)Σ_{k odd} cos(πkα)/k² = |α| on [−1, 1] exactly (AH; equality,
no slack). Finer lattices beat it. Float LP margins t* (dual_lp.py; t* > 0 means dual feasible; logs/dual_*.log):
| θ | lattice h, support of η | λ: t* |
|---|---|---|
| 1 | 0.05, (0, 10] | 0.45: 0.1412; 0.50: 0.0828; 0.55: 0.0326; 0.60: −0.0113 |
| 1 | λ/3, (0, 30] / (0, 40] | 0.58: +0.0132 / –; 0.59: +0.00426 / –; 0.595: – / +0.00116; 0.60: −0.00447 / −0.00321 |
| 4/7 − 10⁻⁵ | λ/15, (0, 10] | 0.70: 0.0716; 0.75: 0.0149; 0.76: 0.0036; 0.765: −0.0020; 0.77: −0.0075; 0.78: −0.0186 |
| 4/7 | λ/10, (0, 14] / (0, 18] | 0.77: −0.00099 / +0.00247 |
| 4/7 | λ/8, (0, 22] | 0.775: −0.00071; 0.78: −0.0062 |
| 4/7 | λ/7, (0, 26] / (0, 34] | 0.775: +0.00090 / –; 0.78: – / −0.00243 |
| 4/7 | λ/5, (0, 26] / λ/4, (0, 26] | 0.775: +0.00081 / +0.00045 |
| 4/7 | λ/5, (0, 40] | 0.78: −0.00144 |
| 4/7 | λ/4, (0, 72] | 0.78: t ≈ +0.0012 as an upper bound, cutting plane cycles under the 1,100-row cap |
| 4/7 | λ/3, (0, 70] | 0.78: +0.000224 on the 0.0002 grid (Arb check FAILED, below); 4·10⁻⁵ grid: cycles at t ≈ +0.00022 |
| 4/7 | λ/5 + kink-exact ψ_h tail, (0, 10] / (0, 20] | 0.78: −0.0184 / −0.0074 (no gain over the plain lattice) |
The support length matters more than the lattice step. Certificates (verify_dual.py: Arb, 360 bits, rational η, Taylor cells with
|P''| ≤ M2), all PASS:
| θ | λ | h | cells | min G ≥ | log |
|---|---|---|---|---|---|
| 4/7 | 0.775 | 31/280 | 90,324 | 0.000204 | logs/verify_t47exact_0775.log |
| 4/7 | 0.77 | 77/1000 | 64,937 | 0.00193 | logs/verify_t47exact_077.log |
| 4/7 (and 399993/700000) | 0.76 | 19/375 | 49,344 | 0.00317 | logs/verify_t47exact_076.log |
| 1 | 0.595 | 119/600 | 100,841 | 0.000239 | logs/verify_t1_0595.log |
| 1 | 0.59 | 59/300 | 50,848 | 0.00361 | logs/verify_t1_059.log |
| 1 | 0.55 | 1/20 | 50,000 | 0.0322 | logs/verify_t1_055.log |
Failed: θ = 4/7, λ = 0.78, h = 13/50, K = 269, 480,771 cells: min −0.0021 at α ≈ 1.4769 (logs/verify_t47exact_078.log): the
LP's 0.0002 grid missed a dip of the X = 70 polynomial. So 0.78 is not certified.
Rerun: uv run --no-project --with python-flint python3 verify_dual.py logs/dual_t47_k7X26.json 20000 4/7 (34 s, 78 MB).
So the method's floor lies in [0.775, 0.7851] for θ < 4/7 (upper end certified by SmallGap; its float SDP reaches 0.7849,
uncertified) and in [0.595, 0.6039] under RH (upper end CGdL's published 0.6039; our float SDP of the same LP: 0.6040).
Floor vs support (float LPs; λ*(X) = zero of t*(λ), slope dt*/dλ ≈ −1.11 from the λ/15 row): X = 10: 0.7632; 14: 0.7691; 18: 0.7722;
22: 0.7744; 26: 0.7758; 34: 0.7778; 40: 0.7787; 70: > 0.78 (h = λ/3, float). Fit λ*(X) ≈ 0.7839 − 0.207/X, residuals ≤ 2·10⁻⁴ (the lattice step, λ/4 to
λ/15, moves λ* by ≤ 3·10⁻⁴). So the support a certificate needs grows like 0.21/(λ_∞ − λ), λ_∞ ≈ 0.784 (extrapolated, not
certified): λ = 0.78 needs X ≈ 53, 0.782 about 105, 0.783 about 210. λ_∞ ≈ 0.784 lies just below the float primal 0.7849 (still
falling with the degree) and the certified 0.7851: consistent with no duality gap, so the method's true floor is about 0.784.
The loss is not the kink of |α| at 0: appending the h-discrete sine's exact tail (ψ_h in closed form, one LP column) leaves λ*(X)
unchanged. It is the truncation of an AH-like pseudo-lattice (doublets (λ, 2λ), (3.14λ, 4.14λ), (7.43λ, 8.43λ) at λ = 0.775).
Mixtures of a-discrete sine pair measures with a ≥ λ as the only background are worse (dual_dpp.py: t*(0.775) = −0.018 at X = 15).
θ = 1: λ*(30) = 0.5949, λ*(40) = 0.5963, two-point fit λ*(X) ≈ 0.600 − 0.17/X: the RH floor is probably near 0.600, so CGdL's
0.6039 is about 0.004 above what the method can give (float, not proved).
Dual shape (λ = 0.775): atoms of weight 0.757 (6.8 × uniform) at u = λ, 0.37 at 2λ, 0.43 at 3.14λ, 0.44 at 4.14λ, nothing on
[0.89, 1.44] and [2.77, 3.10]: an AH-like pseudo-lattice with spacing ≈ λ. It satisfies the method's pair-correlation input; the
real zeros do not look like it:

Gap laws by HKPV (Fredholm determinant of the sine kernel; Palm laws of the discrete sines) and real zeros (hkpv_gaps.py,
logs/hkpv_gaps.log; zeros: ExploreHelio's 22,491 zeros below 2·10⁴; gram_win.c -zeros at Gram index 2.48·10⁸ (t ≈ 1.0·10⁸) and
3.5·10¹¹ (t ≈ 9.8·10¹⁰), 34 windows of 62 Gram intervals each, 4 s and 91 s, 12 MB):
| configuration | D(0.5) | D(7/11) | D(0.76) | D(0.7851) | D(1) |
|---|---|---|---|---|---|
| sine process (GUE limit) | 0.1131 | 0.2078 | 0.3125 | 0.3352 | 0.5339 |
| 1/2-discrete sine (AH, θ = 1) | 0.2974 | 0.2974 | 0.2974 | 0.2974 | 0.7500 |
| 7/11-discrete sine (θ = 4/7) | 0 | 0.5046 | 0.5046 | 0.5046 | 0.5046 |
| ζ, 100 < t < 2·10⁴ (22,461 gaps) | 0.0900 | 0.1819 | 0.2937 | 0.3190 | 0.5364 |
| ζ, t ≈ 1.0·10⁸ (2,074 gaps) | 0.1148 | 0.2049 | 0.3105 | 0.3322 | 0.5154 |
| ζ, t ≈ 9.8·10¹⁰ (2,073 gaps) | 0.1143 | 0.2113 | 0.3136 | 0.3353 | 0.5408 |
| Davenport–Heilbronn, on-line zeros 20 < t < 420 (297 gaps) | 0.0505 | 0.1178 | 0.1987 | 0.2189 | 0.4949 |
(binomial s.d. ≈ 0.010 at 2,000 gaps.) At 10⁸–10¹¹ the zeros match the sine process within 1 s.d.: about 31% of the gaps are
≤ 0.76, yet the method certifies nothing at or below 0.775. HKPV detail: the sine kernel on [0, 0.76] has Bernoulli parameters
0.6554, 0.1016, 0.0030, 3·10⁻⁵; E(0; 0.76) = 0.3087; P(N ≥ 2) = 0.0685.

## 4. Adversarial test
- Theorems A and B are statements about test functions: they hold whatever the zeros are, so they carry no RH information; they
  locate where a method breaks.
- The node quantity sees ordinates only. A ζ zero set with one pair moved to β ± iγ keeps every ordinate, so D(λ) and the method's
  conclusion are unchanged (Lamzouri's E_λ already contains every off-line zero). Davenport–Heilbronn's on-line zeros have
  D(0.76) = 0.20 > 0 too (9 off-line zeros below 420 thin the line; mean normalised gap 1.06). The AH-type discrete sine and the
  dual pseudo-lattice both satisfy the pair-correlation input. So no DPP gap statement separates ζ from RH-violating or AH
  configurations, as expected for N8 gap results.

## 5. Novelty
- arXiv listings math.NT and math.PR, 2026-10 and 2026-11 (titles with gap, spacing, pair correlation, determinantal, Montgomery,
  LP): 2610.08965, 2610.22927 (moderate gaps), 2610.23467, 2610.26038 (μ ≤ 1/2 on RH, resonance), 2610.28527, 2611.14486; the DPP
  papers 2610.11554, 2610.16963, 2610.17669, 2611.01379 are off-topic. None gives a dual or a floor for a small-gap method.
- Lamzouri 2611.14486, Remark l. 199: "no serious attempt to obtain the best constant", no limit. CGdL 1810.08843: primal SDP only
  (no "dual" in the TeX). Lagarias–Rodgers 1905.12123 (AH, full Rudnick–Sarnak range) and 1907.03391 (band-limited mimicry by
  lattice processes, all n-level correlations): Theorem A is the pair-only case of their question, elementary, likely known; neither
  uses it as an LP dual. 2201.10676: limit 0.5042 for the Montgomery–Odlyzko method (infinitely many gaps, another method).
  Zenodo ("small gaps" zeta; pair correlation + dual): nothing relevant (24095120: resonance r-gaps). Our SmallGap/NOVELTY.md: floor
  unknown.
- Theorem B is new as far as checked. Gate: a barrier on the N8 record line (Thm Main2 cannot pass 0.775; record 0.7851); the
  reviewer judged the core correct and significant as a method-specific barrier (repairs below).

## 6. Process note
The first LP run (dense 10,001 × 391 rows × columns, no watchdog) reached 956 MB for 600 s before its alarm killed it: a breach of
the 300 MB cap. Every later run had a 280 MB self-abort; three runs self-aborted there (a one-shot 4,300-row LP, λ/15 on (0, 14], λ/10 on
(0, 20]) and three hit their 300–900 s alarms without output (λ/25 on (0, 12], θ = 1 with h = 0.038, λ/15 on (0, 14]). Follow-up
(dual_lp2.py, bounded active set, chunked P): the first λ/4, (0, 72] run self-aborted at 285 MB with 1,365 rows; a second dropped
binding rows and cycled (t* stable, fine-grid violations recurring), fixed by keeping all rows up to a 1,100-row cap and η ≤ 6.
Kept runs peak at 84–263 MB (LP), 59–239 MB (Arb; the failed 480,771-cell check), 31 MB (HKPV), 12 MB (gram_win).

## Suggested NOGO line
- explore34 (Orca, exploration: determinantal point processes (HKPV, Ghosh–Peres rigidity, Soshnikov, strong Rayleigh) → N8, N13) BARRIER: DPP structure adds no gap bound (the a-discrete sine DPPs share every band-limited correlation and have no gaps below a), but as LP duals they cap the small-gap method: Lamzouri's Thm Main2 (any θ < 4/7) cannot reach λ ≤ 0.775 (Arb-checked dual; the DPP alone gives 7/11), so the record 0.7851 is within 0.0101 of its floor (float: ≈ 0.784, the dual floor grows like 0.784 − 0.21/X with its support X); under RH the pair-correlation floor is ≥ 0.595 (CGdL 0.6039); N13's prolate ratios are arithmetic, HKPV only re-encodes δ_n [ExploreDPP](Riemann/Deep/ExploreDPP/REPORT.md)

RH relevance: N8 (floor of the small-gap method behind the record λ ≥ 0.7851); no RH information.

## Files
REPORT.md, DERIVATION.md; dual_lp.py (LP, cutting planes, 280 MB watchdog; args θ k0 X λ-list out.json), dual_lp2.py (bounded
active set, chunked P; args θ k0 X-list λ-list tag [-v] [-keep], env DUAL_B), verify_dual.py (Arb check;
args cert.json 1/δ [θ]), dual_dpp.py (a-discrete-sine mixtures as background; worse), hkpv_gaps.py (HKPV gap laws, real-zero statistics); logs/: dual_t47_k7X26.json (λ = 0.775 certificate),
dual_t47_k10X18.json (0.77), dual_t47_b.json (0.76), cert_t1k3_X40_lam0.595.json and cert_t1k3_X30_lam0.59.json (θ = 1),
dual_t1.json (θ = 1, 0.55), cert_t47k3_X70_lam0.78.json (float only, Arb FAIL), dual_*.log, dual2_*.log, verify_*.log/.time, hkpv_gaps.log,
zeros_1e8.txt, zeros_1e11.txt (gram_win -zeros output).

## Repairs after review (independent reviewer of DERIVATION.md: core CORRECT, SIGNIFICANT as a method-specific barrier)
1. Lemma 1 now assumes θ ≤ 1/(2h) (was θ < 1/(2h), which Theorem A's θ = 1, h = 1/2 violated). The proof only needs
   P − |α|1_θ ≥ 0 pointwise: on the base period by hypothesis, beyond it by periodicity; at θ = 1/(2h) the support of 1_θ is the closed
   base period and its shared endpoints are a null set where P ≥ 0 anyway (DERIVATION §1).
2. Theorem A: the false bound "S_a ≥ 1 − (2 − 1/a)" (S_a(1/a) = 0) is replaced by: one triangle contributes at most 1, two
   overlapping triangles sum to 2 − 1/a ≤ 1, so S_a ≥ 0 (DERIVATION §2).
3. The certificate check is spelled out (DERIVATION §3): θ is a grid point; the left cell ending at θ certifies the binding lower
   branch P(θ) − θ (G jumps up at θ); Arb outward rounding; rational η with η_k = −1 exactly for k < k0; closed cells covering
   [0, 1/(2h)].
4. Floor intervals use certified upper ends: [0.775, 0.7851] (0.7849 is an uncertified float SDP value) and [0.595, 0.6039] under RH
   (CGdL; our 0.6040 is float).
5. §2 and DERIVATION §4 no longer claim that every rigorous gap argument reduces to this LP: band-limited bounds of this type are LP
   bounds; mollifier, moment and resonance arguments are outside the statement.

Orchestrator note (12 Nov 2026): under RH, Durkan–Pearce–Crump (arXiv 2610.28527) prove D(λ) > 0 for every λ > 1/2 by another argument, so the θ = 1 floor 0.595 caps only this functional, not the RH case. Both certificates (0.775, θ = 4/7; 0.595, θ = 1) rerun by the orchestrator: PASS, identical minima.

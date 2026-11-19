AUDIT PASSED: every item of chebotarev-numeric holds as written, given Cureton–Thorner's Lean 31/32 (A) and Liu's B_2 (B) as audited; one formal slip (CT10000's binary float 0.094 sits 2.2·10⁻¹⁹ above 0.094) was closed by a decimal-ball rerun, ALL CERTIFIED.

# chebaudit: line-by-line referee of "numeric polylogarithmic least prime in every Chebotarev class" (node N11)

Scope: discoveries/chebotarev-numeric.md and Riemann/Deep/ChebExp/ (DERIVATION.md, cert.py, chains.json, direct.gp, logs), read-only on both. Chain with line numbers: [FULLPROOF.md](FULLPROOF.md). Code and logs: this directory. 19 Nov 2026, 20:20 UTC on.

## Verdict per item

| Item (page) | Verdict | Evidence (this audit) |
|---|---|---|
| (A) N𝔭 ≤ (0.49·log 3Δ_L)^32, every Galois L/K, L ≠ ℚ, every C; 𝔭 unramified, degree one | holds | implied by the certified 0.484 (monotone in c); cert rerun on a copy: ALL CERTIFIED, output identical to logs/cert_output.txt (95 s, 38 MB) |
| (A) 0.29 / 0.151 / 0.108 / 0.094 for log Δ_L ≥ 10 / 100 / 1000 / 10⁴ | holds | chains start exactly at 10, 100, 1000, 10⁴ (coverage.py, exact rationals); 0.094 certified as a decimal only after the fix (cert_dec.py, logs/cert_dec_output.txt) |
| (B) N𝔭 ≤ (0.73·log 3Δ_L)^{7.8978} given Liu's B_2; 0.595 / 0.426 / 0.366 / 0.344 | holds (input as audited) | implied by certified (0.62·log 3Δ_L)^{1/δ_2}: 1/δ_2 = 7.897769… < 7.8978 and c·log 3Δ_L ≥ 0.62·log 9 > 1 |
| A uses only H(1/32) and classical inputs; B adds only H(δ_2) (and PARI nflist for Δ_L ≤ 257) | holds | δ enters cert.py only via DELTAS; CT chains never touch B2; A's direct range is 8 Kronecker computations |
| Sharpness: ℚ(√−3), C = {1}, least split prime 7 forces c ≥ 7^{1/32}/log 9 = 0.48365; 0.49 within 1.3 % | holds | 0.4836540909 (256 bits); 0.49/0.48365 = 1.0131; (0.484·log 9)^32 = 7.162 ≥ 7 |
| Deuring reduction, weight θ ∈ [0, 1] also at ramified primes, ζ_L = ∏ L(s, χ) | holds | Lemma 3.2 and §4.1 re-derived (FULLPROOF S1, S3, S4) |
| Explicit formula, contour to −∞, trivial zeros, pole term | holds | Lemma 3.1 re-derived; orders m_j from Γ_ℝ^{r_1}Γ_ℂ^{r_2}; typo only (e^{\|s*\|U} for e^{s*U}, D l.120) |
| Hadamard identity at σ_0 > 1 with Γ-saving −n_L g(σ_0) | holds | re-derived; hadamard.gp: identity to ≤ 2.5·10⁻⁵ for ℚ(√d), d = −3, 5, −4, 12 (zeros of ζ and L(s, χ_d) separately, T = 300); g ∈ [0.8074, 1.1774] over all 405 kernels |
| Maximum principle for \|x^{s−1+δ}Φ\|/Re(1/(σ_0 − s)) (the novel step) | holds | subharmonicity re-derived; reviewer CORRECT; maxprin.py: interior/boundary max 0.014–0.86 on dense 2-D grids, six kernels |
| Ramified primes via v_𝔭(d_{L/K}) ≥ \|G\|/2 | holds | v = gf·v_𝔔(𝔇) ≥ \|G\|(e − 1)/e; three per-prime bounds, Σ min ≤ min Σ, Σ v log N𝔭 ≤ Λ |
| Prime powers: exact ≤ 2·10¹², θ(y) < 1.01624y beyond | holds | partial summation re-derived; smooth bound 7–117× the exact sum where both apply; θ(p)/p ≤ 0.99963 for p ≤ 3·10⁶ |
| n_L ≤ log Δ_L/log(π/2), affine in n_L, both ends checked | holds | Minkowski induction re-derived; negative-B endpoint harmless (reviewer CORRECT) |
| Ten chains, 405 kernels, monotone tails | holds | contiguity exact (no gap, no overlap); 46 entries (all first/last/min-margin/tail + 15 random) recomputed at 256 bits with own code: 0 flags, own margin = cert's + 0.6–1.4·10⁻³; tails positive by direct evaluation at Λ_N·f, f up to 10¹⁵ |
| Direct ranges: A Δ_L ≤ 12 (8 fields), B Δ_L ≤ 257 (563 cases) | holds | 471 C2 rows identical to own Kronecker values; 92 others re-witnessed without idealprimedec/idealfrobenius/rnfidealprimedec, all ≤ B bound; the 27 rows where the file's prime is smaller all have p \| poldisc(f) (diffcheck.gp) |
| Certificate paragraph (104 s, 38 MB, margin ≥ 0.0030, safe rounding) | holds | rerun 95 s, RSS 38 MB, min margin 0.00300; rounding directions checked line by line |
| Checks paragraph (1 013 926 cyclotomic classes, 0.4837; 508 S3 closures, 1524 classes, 0.108) | holds | Σφ(q), 3 ≤ q ≤ 2000, q ≢ 2 (4) = 1 013 926; s3_closure.gp: 508, 1524, 0.1081 |
| Before paragraph (Δ_L^{310}, 12577, ≪ Δ_L^{40}, ≪ Δ_L^{16}, FSW Cor. 1.2 non-numeric c(δ)) | holds | arXiv TeX of 2105.14181, 1508.00287, 1902.08640, 2610.18491 read (FULLPROOF S7) |

## (1) Referee notes (details in FULLPROOF.md)
- Every imported statement was checked for hypotheses and normalisation: Cureton–Thorner's Lean theorem covers every number field and 31/32 < Re s ≤ 1 (odlaudit, HPAudit; statement re-read here, Challenge/DedekindZetaChallenge.lean l.30-32); the functional equation then gives δ ≤ β ≤ 1 − δ. Liu's theorem covers every finite-order Hecke L-function over every number field, so every ζ_L (Hexagon2 l.1, l.83); only ζ_L of the top field is used (D l.16). The Γ-factors Γ_ℝ = π^{−s/2}Γ(s/2), Γ_ℂ = 2(2π)^{−s}Γ(s) enter only through the trivial-zero orders and through (3.2), both re-derived and (3.2) checked numerically.
- Directions: main term Φ(1) as a lower bound; every error term (zeros, trivial zeros, ramified, powers) as an upper bound; the zero bound uses B(σ_0) ≥ Σw ≥ 0 of the actual field, so upper-rounding q is safe even where the affine extrapolation makes B negative.
- The theorem's x is (c(Λ_a + log 3))^{1/δ} ≤ (c·log 3Δ_L)^{1/δ} because Λ ≥ Λ_a; Good(x) > 0 forces p < x.

## (2) Code audit (cert.py, chains.json)
- Outward rounding: up()/lo() give exact endpoints; every comparison is Arb's certain comparison; Φ(1), g as lower bounds, Φ(0), Φ(½), Φ(⅓), ‖k‖, S_2, S_3, A_1, A_2, P_k, τ, R_f, 1/σ_0 + 1/(σ_0 − 1), x^{−δ} as upper bounds. triv_ratio divides by the upper bound P0 but is used only as P0·ε, which is ≥ Σ x^{−j}Φ(−j): safe.
- B&B (cert l.109-144): tiles [0, G] with exact dyadic balls, accepts only certified upper bounds, tail bound decreasing; own independent sups are 0.99897–0.99902 of cert's A_1 (46 entries).
- Coverage (coverage.py): contiguity is exact, so cert's 10⁻¹⁵ slack in its own assert (cert l.340) is never used. Thresholds 2.5649 < log 13, 5.5529 < log 258, then exactly 10, 100, 1000, 10⁴ (lesson 277: no scan starts above its threshold).
- Independent recomputation (indep.py, mpmath 256 bits, own kernel, own sieve and exact prime-power sums, own digamma, no cert code): logs/indep.log, 46 entries, all margins positive and ≥ cert's (cert's are lower bounds; the gap is the 10⁻³ B&B tolerance).
- Constants: 9 of the 10 floats c lie below their decimals; CT10000's 0.094 lies 2.2·10⁻¹⁹ above. cert_dec.py (arb(repr(c))) re-certifies every chain with decimal balls: ALL CERTIFIED, same margins.

## (3) Falsification (no counterexample)
- Quadratic fields 3 ≤ |d| ≤ 10⁵ (60 786 fields, 182 358 (L, K, C) cases: K = ℚ split/inert, K = L): no violation of A (0.484, 0.49; 0.29 for log Δ ≥ 10) or B (0.62, 0.73 with 7.8978; 0.595); largest constants 0.4836541 (A) and 0.5822778 (B), both ℚ(√−3), C = {1}, p = 7; next A: 0.4710 (ℚ(√−3), K = L, p = 3), 0.4651 (inert, p = 2), 0.4232 (ℚ(i), split, p = 5) (falsify.py, logs/falsify.log).
- ℚ(ζ_q)/ℚ, 3 ≤ q ≤ 600, q ≢ 2 (4), 91 266 classes: no violation, max 0.483654 at q = 3. S_3 closures of the 508 cubic fields |d| ≤ 3000 (K = ℚ, 1524 classes): no violation, max 0.1081. D_4 quartic fields |d| ≤ 10⁴ (413 fields; K = L and K = k split/inert, 1239 cases, factormod witnesses): no violation, max 0.1813 (A), 0.2183 (B) (d4_quartic.gp, logs/d4_quartic.log).
- Small-Δ_L range where the bound is tight: A is tight only at Δ_L = 3 (7 ≤ 7.162); B's direct range rows all re-witnessed; the analytic chains start at margin 0.0030 (log Δ_L = 2.5649), where the actual least primes of all fields with 13 ≤ Δ_L ≤ 48 (69 rows of direct_range.txt) are ≤ 13 against x ≈ 9·10⁷.

## Findings (none affects a bound)
1. CT10000: binary float c = 0.094 + 2.2·10⁻¹⁹; fixed by cert_dec.py (decimal ball), all chains still certified.
2. The page states 0.49 and 0.73, weaker than the certified 0.484 and 0.62 of REPORT/DERIVATION (both true; the page could quote the stronger ones).
3. Wording: "405 kernels plus monotone tails" counts the 10 tail kernels inside the 405; D l.120 has e^{|s*|U} for e^{s*U}.

## Residual risks and what this audit did not check
- Inputs taken as audited elsewhere: Cureton–Thorner (HPAudit, odlaudit: third-party kernel logs, no local replay); Liu's B_2 (Hexagon2: conditional on Liu 2610.26835v2, unrefereed). If 31/32 fell, A falls; if B_2 fell, B falls; neither touches the other.
- Classical facts cited, not re-read from the sources (no arXiv TeX): Rosser–Schoenfeld 1962 Thm 9 (but implied by Platt–Trudgian's θ(x) < (1 + 7.5·10⁻⁷)x, arXiv 1407.1914, TeX read), Stark 1974 Lemma 1, Lagarias–Odlyzko 1977, Iwaniec–Kowalski Prop. 5.7, Minkowski, Dedekind's different bound. Each was re-derived or checked numerically where a constant enters (θ, (3.2), Minkowski); the contour estimates contribute no constant.
- PARI nflist completeness for degree ≤ 4, |disc| ≤ 257 (B only); counts match the classical tables. Hu (Hexagon 2610.00216) and Bach–Sorenson not re-read; novelty not re-audited beyond the Before items.

## (4) Dependents (unaffected: nothing fails)
DISCOVERIES.md l.11 (half-plane consequences line); TREE.md l.15 (N11 record and serrenum/koiran mentions); COROLLARIES.md l.41 (detfactor: ℓ ≤ (1.96q² log(12q²p))^{32}), l.42 (serrenum, Riemann/Deep/SerreNum), l.43 (koiran, Riemann/Deep/KoiranAM: x₀ via chebexp's certified kernel, its own constants not checked here), l.44 (harvest4: split prime ≤ (0.49 log 3|d|)^{32}, Riemann/Deep/Harvest4); Riemann/Deep/Harvest2/HARVEST.md; Riemann/Deep/LiuAudit/REPORT.md; Riemann/Tree/R399.lean.

## Review
One reviewer subagent, given the Hadamard step (3.3), the maximum principle (3.4) and the affine endpoint logic of §5.2 inline with exact formulas: (i) subharmonicity and boundary sup CORRECT; (ii) Γ-factor signs, dropping ζ'/ζ, g CORRECT; (iii) endpoint logic and upper-rounded q with possibly negative artificial-endpoint B CORRECT. Its one flag, L = ℚ (Λ = 0, n_L = 1), is excluded by the hypothesis L ≠ ℚ (D l.18) and by Λ_a ≥ 2.5649 in every interval: verdict CORRECT.

## Files
coverage.py, indep.py, gcheck.py, maxprin.py, falsify.py, compare_direct.py, direct_check.gp, diffcheck.gp, s3_closure.gp, d4_quartic.gp, hadamard.gp, cert_dec.py; logs/ (coverage, indep, gcheck, maxprin, falsify, compare_direct, direct_check, diffcheck, s3_closure, d4_quartic, hadamard, cert_rerun_output, cert_dec_output). Footprints (memcap peaks / RSS): cert rerun 38 MB RSS, cert_dec 22 MB, falsify 109 MB RSS (94 MB memcap), indep 34 MB, gcheck 14 MB, gp scripts parisize ≤ 128 MB, each < 40 s. Rerun: `uv run --no-project --with mpmath python3 indep.py 15` (≈ 2 min), falsify.py 600 (≈ 90 s); gp via a perl fork-and-kill timer.

RH relevance: N7, N11

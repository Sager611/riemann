AUDIT PASSED: Lemmas 2.5–2.9, the Riemann–Siegel input and the block bookkeeping of arXiv 2602.05614v1 are correct as our 0.568 record uses them (re-derived line by line; only misprints and omitted routine steps, all filled; Gabcke's 0.127 t^{-3/4}, t ≥ 200 matches; 1.463 needs |C₀| ≤ 0.92406, now proved in Arb; 0 failures in ≈ 1.4·10^6 Arb/exact checks). The record does not need them: with Hiary–Patel–Yang's refereed lemmas instead, the same certificate passes at 0.568 (max 0.56797).

# RevAudit — referee report on Revers, arXiv 2602.05614v1 (node N9, load-bearing input of RECORDS N9)

Sources: TeX `ZetaRevers.tex` (e-print, 2144 lines; line numbers below), Hiary–Patel–Yang (HPY) arXiv 2207.02366 TeX `main.tex`, Gabcke's thesis (2015 LaTeX edition of the 1979 dissertation, eDiss `DissGab1979.pdf` via the Wayback Machine, 133 pp.). Lemma numbers: lemma counter per section, so Lemma 2.1 = l. 180 … 2.9 = l. 455.

## 1. Statements we use

| Lemma (TeX) | Statement | Hypotheses → where used in its proof | Our use (Explicit611 `certify.py`) |
|---|---|---|---|
| 2.3 l. 272–298 | f′ monotone, continuous, z+U ≤ f′ ≤ z+V, 0<U≤V<1 ⇒ \|Σ_{a≤n≤b} e(f(n))\| ≤ (1/U+1/(1−V))/π | = HPY Lemma 2.3 (l. 179–186: ℓ+1/U_H ≤ f′ ≤ ℓ+1−1/V_H on [a,b), bound (U_H+V_H)/π) with U = 1/U_H, V = 1−1/V_H; HPY's proof (l. 650–755: summation by parts with G = 1/(1−e(g)), telescoping cot, cot y < 1/y) re-checked: correct | inside 2.6(b) |
| 2.4 l. 308–327 | Weyl differencing, any L, M ≥ 1 | = HPY Lemma 2.4 (l. 195–215) | 2.9; exact-Weyl blocks |
| 2.5 l. 329–353 | g′(L−1−m)−g′(0) ≤ g′(K−1−m)−g′(0) ≤ mKμ/W for m<L≤K | m,r,L,K ≥ 2, K ≥ K₀ > 1 → μ only | 2.6(b) l. 1273–1279 only; μ_r with K₀ = K (item 1) |
| 2.6 l. 355–398 | (a) m/W ≤ g″ ≤ mλ/W on [0,K−m]; (b) \|Σ_{n≤L−1−m} e(g(n))\| ≤ (mKμ/W−1)B + P(ε₁) + P(1−ε₂) − 2/(πΔ), B = 1+2WΔ/m+2/(πΔ) | m<L≤K (2.5), 0<Δ<1/2 (segmentation), g″ ≥ m/W (lengths); m ≥ 2 never used | 2.7(e) |
| 2.7 l. 400–427 | 1≤m≤min{M,K−1}=M, L≤K: \|Σ\| ≤ 4μK m^{1/2}/√(πW) + μKm/W + 4√(W/π) m^{−1/2} + 1 − 6/π | t ≥ 5·10^6 and φ ≤ 0.35 → (a) R ≥ 4 only; t ≥ 2012, K ≥ t^{1/3}, r ≤ R, c<3 → (b) M ≤ K−2; r ≥ 4, K ≥ t^{1/3}, m ≤ M, c<3 → (c) Δ<0.154; (d),(e) algebra | Parts C, E, F (per-m bound X_m); 1.58·10^6 ≤ t < 5·10^6 with R ≥ 4 asserted (item 8) |
| 2.8 l. 433–453 | (a) Σ(1−m/M)m^{−1/2} ≤ (4/3)M^{1/2}; (b) = (M−1)/2; (c) ≤ (4/15)M^{3/2}; (d) ≤ M²/6 | none | relaxed block bound |
| 2.9 l. 455–537 | \|Σ_{n<L} e(f(n))\|² ≤ (K/W^{1/3}+c)(αK+βW^{2/3}), L ≤ K | as 2.7; W^{1/3} ≥ 5π^{1/3} (r ≥ 4, K ≥ t^{1/3}) for β > 0 | relaxed block bound in the partial-L form of l. 502 |
| RS l. 216–234, 753–764 | \|ζ(½+it)\| ≤ 2\|Σ_{n≤n₁} n^{−½+it}\| + cos(π/8)(2π)^{1/4}t^{−1/4} + 0.127t^{−3/4} ≤ … + 1.463t^{−1/4} + 0.127t^{−3/4} | t ≥ 200 (Gabcke) | Parts B, C (t ≥ 300), E, F |
| §3 l. 765–850 | T₁ = 2Σ_{n<Kr₀} n^{−½}, blocks [Kr, K(r+1)), last block [KR, n₁]; Abel summation 1/√(Kr)·max_L | R ≥ r₀ = 4 | same split per n₁; last block with its true length L_R = n₁−KR+1 (l. 2045–2049) |

## 2. Proof check (line by line)

- **2.5** (paper's proof l. 348–353 is a pointer to HPY (5.1)–(5.2)). Own proof: g′(x) = −(tm/2π)/((Kr+x)(Kr+x+m)) is increasing (first inequality), and with ε = 1/K, (2π/(tm))K²[g′(K−1−m)−g′(0)] = 1/(r(r+mε)) − 1/((r+1−(m+1)ε)(r+1−ε)) ≤ 1/r² − 1/(r+1)² = 1/(r²(r+1)) + 1/(r(r+1)²) ≤ 1/(r²(r+1)) + 1/(r(r+1)²(1−1/K₀)) = 2μ/(r+1)³. Holds for all integers 0 ≤ m ≤ K−1, r ≥ 1 and every K₀ > 1 (K₀ ≤ K is not needed; K₀ = ∞ works). Item 1 (K₀ = K, per-block μ_r): CORRECT.
- **2.6(a)** ("standard methods", l. 1260). g″(x) = (t/π)∫_x^{x+m} dy/(Kr+y)³ with Kr+y ∈ [Kr, K(r+1)] ⇒ m/W ≤ g″ ≤ mλ/W. CORRECT for any m > 0.
- **2.6(b)** (l. 1262–1858). Checked: Estimatek (l. 1276–1278); well-definedness of x_j, y_j (g′ < 0 on (−Kr,∞), every C_j ± Δ < 0); ordering (l. 1302–1308) from 2Δ < 1 and the caps; y_j − x_j ≤ 2WΔ/m by the mean value theorem at points of [0, L−1−m] ⊂ [0, K−m]; integer counts ≤ length + 1; every Kusmin–Landau call (z, U, V) has 0 < U ≤ V < 1; the three left-boundary cases (l. 1388–1466); the right-boundary cases the paper omits (l. 1472) and the k = 1 details (l. 1577), re-derived by the reflection ε ↦ 1−ε₂; the six k = 0 cases (l. 1623–1839) are complete because ε₁ ≤ ε₂; the last step for k = 0 (l. 1852) needs mKμ/W + ε₁ − ε₂ ≥ 0, which is Estimatek with k = 0 (implicit, valid). Remark 2.1's reason ("bound independent of the length") is not a proof for closed intervals, but the claim holds: extend f linearly beyond b (f′ stays in [z+U, z+V], monotone) and apply HPY's [a, b+½) form. CORRECT.
- **2.7** (l. 1860–2027). (a) R > (√(t/2π)−1)/(t^φ+1) − 1, increasing in t, decreasing in φ: ⌈3.0115⌉ = 4. (b) 3π^{1/3}(t^{1/6}/√(2π)+1)+2 ≤ t^{1/3} for t ≥ 2012 (Arb margin +0.00253 at 2012, +0.00095 at 2011), so K ≥ cW^{1/3}+2 > M+1, i.e. M ≤ K−2. (c) Δ = √(m/(πW)) < √((c+1)/(25π^{5/3})) < 2/(5π^{5/6}) = 0.154087 < 1/6. (d) all three cases; the case-2 margin at ε = 1−Δ is Δ(1 − 1/(π(1−Δ))) > 0. (e) algebra l. 2007–2014 and the m ≥ L case (bound > 18/π + 1 > 0). CORRECT. Item 8 (Explicit611 uses the bound for 1.58·10^6 ≤ t < 5·10^6 where R ≥ 4 is asserted): CORRECT — t ≥ 5·10^6 and φ ≤ 0.35 enter only through (a).
- **2.8** (l. 451–453 cites HPY p. 215). (b) exact; (d) Σ(1−m/M)m = (M²−1)/6; (a) h(x) = x^{−½} − x^{½}/M decreases on (0, M] ⇒ Σ h(m) ≤ ∫₀^M h = (4/3)M^{½}; (c) (1−x/M)√x is concave and vanishes at 0, M ⇒ trapezoid sum ≤ ∫ = (4/15)M^{3/2}. CORRECT.
- **2.9** (l. 482–537). l. 486–496: weights, (b) used as an equality because 1−6/π < 0; l. 502–503: L ≤ K and the bound is uniform in L; l. 505–515: each of the six terms with cW^{1/3} ≤ M < cW^{1/3}+1 bounded in the right direction (6/π−1 > 0 needs the upper bound on 1/M); l. 516–519 identity; product step needs α, β > 0 (l. 521–534, βW^{1/3} ≥ 24.53). Explicit611's partial-L form (L kept in (L+M−1)(L/M+…)) and its monotonicity in L: valid. CORRECT.
- **§3 bookkeeping**. Split l. 767–777 partitions [1, n₁] when R ≥ 4; Abel l. 782–791; T₃ l. 798–813 (n₁ ≤ (R+1)K−1); identity l. 819–824; (3.1)–(3.12) directions (negative terms bounded with t ≤ t₁ and W^{1/3}/K ≥ 5π^{1/3}t^{−1/3}). The paper's own 10 regions: with full-precision coefficients from the paper's formulas (Explicit611 `sanity_authors.py`) the bound/(t^{1/6} log t) is ≤ 0.61099995 on every region (`regions_full.py`); the displayed inequalities with the printed 6-digit coefficients (a₁ rounded up) fail literally by ≤ 3·10^{−7} at three region ends (up to 0.6110003, `regions_check.py`): a display issue only.
- **Misprints (no effect)**: 2.5/2.6 assume m, L ≥ 2, 2.7 uses m = 1 and L = 1 (proofs need only m ≥ 1, L−1−m ≥ 0); l. 1316 (y_{k−1}, x_{k−1}) → (y_{k−1}, x_k); l. 1407 (g′)(C₀+Δ) → (g′)^{−1}(C₀+Δ); l. 1871, 1875, 1880 √(t₀/√(2π)) → √(t₀/(2π)) (the value 3.0115 is the correct one); l. 1922 π\frac{5}{6} → π^{5/6}; ν/υ mixed (l. 1333–1336, 1402–1405). No downward rounding: 1.463 ≥ 1.462717, 2.0895 ≤ 2.08958, 2.8805 ≤ 2.88052, 0.154088 ≥ 0.154087.

## 3. Riemann–Siegel input vs Gabcke

- Gabcke (2015 ed.), Einleitung (1), (6), (8) and **Satz 3.2.2(b)** (p. 55): Z(t) = 2Σ_{n≤N} cos(ϑ(t) − t log n)/√n + (−1)^{N−1} a^{−½} Σ_{n≤K} C_n(z) a^{−n} + R_K(t), a = √(t/2π), N = ⌊a⌋, z = 1 − 2(a−N); **for t ≥ 200, |R₀(t)| < 0.127 t^{−3/4}** (in t, not in a; crude version (a): 0.47 t^{−3/4}). C₀(z) = cos(π(z²/2+3/8))/cos(πz), C₀(±1) = cos(π/8) = 0.92388 (appendix table), so the C₀ term is ≤ cos(π/8)(2π)^{1/4} t^{−1/4} = 1.462717 t^{−1/4}.
- Paper l. 223–224, 759–762: same constant, exponent and range: **MATCH**. HPY (l. 150–163, citing Hiary 2016 Lemma 2.1) use the weaker 1.48 t^{−1/4} + 0.127 t^{−3/4}, t ≥ 200.
- One unstated step: the paper's cos(π/8) is max_{|z|≤1}|C₀(z)|, but Gabcke states only |C₀| < 0.93 (Satz 2.3.4(b)), which gives 1.4725 t^{−1/4} > 1.463 t^{−1/4}. Closed here: `c0_bound.py` proves in Arb (adaptive dyadic bisection; Cauchy mean value form at the removable singularity z = ½) that |C₀(z)| < 1.463/(2π)^{1/4} = 0.924058 on [−1, 1] (largest certified upper bound 0.923945 next to z = ±1, where C₀ = cos(π/8)). So 1.463 is valid.
- Residual (not re-derived): Gabcke's 0.127 itself rests on his numerically derived bounds for |C_n| (Satz 2.3.4(b)); it is the standard input of every explicit bound in this area, and our direct test below stays at 0.12229.
- Numbering: in the 2015 edition the bound is Satz 3.2.2; its Satz 4.2.3 is the expansion of ϑ(t). The preface says theorem references were renumbered for LaTeX; the brief's "Satz 4.2.3" is presumably the 1979 numbering. Same statement either way.
- Explicit611 applies it only for t ≥ 300, with the decreasing remainder taken at the left end of each interval: inside the range.
- Own checks (`rs_c1.py`, `falsify.py RS`): max|C₁| = 0.030597, (2π)^{3/4}·max|C₁| = 0.12143 ≤ 0.127, so the t-scale is the right one; mpmath at 856 heights in [200, 3·10^4], both sides of 64 jumps t = 2πM²: max |R₀(t)| t^{3/4} = 0.12229 (at the jump t = 2π·6²), max |C₀ term| t^{1/4} = 1.462717.

## 4. Falsification (`falsify.py`, `logs/falsify.log`; Arb 128 bits, a case fails only if the LHS ball lies above the RHS ball; overlaps: 0)

| Test | Cases | Failures | Max LHS/RHS |
|---|---|---|---|
| 2.8 (a),(c) Arb, all M ≤ 10^6; (b),(d) exact, M ≤ 3000 | 1,003,000 | 0 | margins → ζ(1/2) = −1.46 and ζ(−1/2) = −0.208 |
| 2.5, stronger form K₀ = ∞, exact rationals, K ≤ 120, r ≤ 40, every m | 290,360 | 0 | 0.99721 (m = 0, → 1 as K → ∞) |
| 2.4, random phases | 3,000 | 0 | 1.000000 (equality cases) |
| 2.6(b), random Δ ∈ [0.01, 0.49], m = 1 included | 2,994 | 0 | 0.94725 |
| 2.7 per m, stated range t ∈ [5·10^6, 10^11]: random + adversarial t-scans | 3,000 + 13,000 | 0 | 0.45759 (m = 1) |
| 2.7, item-8 range t ∈ [1.5·10^6, 5·10^6), R ≥ 4 | 3,000 + 11,750 | 0 | 0.33424 |
| 2.7(b),(c) at every sampled point | 6,000 | 0 | min K−(cW^{1/3}+2) = 90.9; min 0.154088−Δ = 0.0179 |
| 2.7(d) P-bound, grid 3000 Δ × 4001 ε + Arb recheck | 1.2·10^7 | 0 | min margin 0.0065 |
| 2.9, stated range: random + adversarial t-scans | 3,000 + 12,250 | 0 | 0.04835 |
| 2.9, item-8 range | 3,000 + 13,500 | 0 | 0.07076 |
| binding interval of 0.568 (n₁ = 605, K = 133, L = 74), 400 heights | 400 | 0 | block²/exact-Weyl bound 0.0874; \|s′_m\|/2.7 bound 0.157; \|ζ\|/full chain 0.068 |
| Gabcke R₀ and C₀ term | 856 | 0 | 0.12229/0.127; 1.462717/1.463 |

Adversarial scans step t by a fraction of the peak width of the sum (Kr/(4L) for 2.9, (Kr)²/(2mL) for s′_m) over 250 steps, with every hypothesis rechecked at both ends. Lemma 2.3 was tested in R433 (41,723 cases, max 0.966). Explicit611's certificate was not rerun (unchanged); its refereed-only variant is in §6.

## 5. Literature

- arXiv 2602.05614: v1 only (5 Feb 2026, 40 pp.), no journal reference, no erratum.
- OpenAlex W7128004822: 3 citing works, all using 0.611 (and 0.566) as input, none checking it: 2603.16023 (Tudzi, divisor problem), 2607.10053 (Tudzi, d₃), 2609.22624 (Fiori et al., zero density; "recent claimed improvements"; also uses Remark 3.1's 0.470795 t^{1/6} log t + 4.04972 t^{1/6} − 21.5437 for t > 7·10^11).
- Zenodo: "2602.05614" 0 hits; creator Revers: nothing on ζ; title (zeta AND explicit): 37 hits, none on t^{1/6} log t bounds (R433: Rafiq 22952309 re-verifies HPY, not Revers). arXiv listings Oct–Nov 2026, math.NT and math.CA (`.scratch/lst.py`): no erratum, check, improvement or formalization.
- Side finding (not our input): 2603.16023 l. 549 quotes HPY for 0.592 t^{1/6} log t on 200 ≤ t < 5.5·10^7; HPY prove 0.592 only up to 10^7 (l. 299–310) and 0.618 on [10^7, 5.5·10^7) (l. 312–337); the Lehman-type bound reaches 0.6179 at 5.5·10^7.

## 6. Dependency map

New check (`make_hpy_variant.py` → `/tmp/revaudit_hpy/certify_hpy.py`): Explicit611's certificate with Revers's 2.7/2.9 replaced by HPY's refereed Lemmas 2.5/2.6 (constant 2 − 4/π, β = 32/(3√(πc)) + (2−4/π)/W^{1/3}; 10 exact replacements, nothing else changed). HPY's lemmas have no t ≥ 5·10^6, r ≥ 4 or c < 3 hypothesis, so item 8 is moot. Run `certify_hpy.py 0.568 ABCEF 500 39894 1e10 6000` (copy under /tmp, 301 s, 47 MB; `logs/run_hpy_0568.log`): **ALL PARTS OK**; max ratios A 0.56693 (Arb, no lemma), B 0.56329, C 0.56797 at (n₁, K) = (611, 152), E (t ≥ 10^10) 0.560078, F (t ≥ 8.97·10^17) 0.482061 (same numbers in the exploratory run at 0.65, `logs/run_hpy_065.log`). Reviewer (verdict-only, formulas inline): all five swapped bounds CORRECT, "no unrefereed constant remains necessary".

| Entry | Inputs from 2602.05614 | If 2.5–2.9 failed | Refereed-only (HPY + Gabcke + Arb) |
|---|---|---|---|
| RECORDS N9: 0.568, t ≥ 3 | 2.5–2.9, RS form, item 8 | certificate void | **holds: 0.56797** |
| 0.5598, t ≥ 10^10 | 2.5–2.9 (φ = 1/3) | void | 0.5601 |
| 0.4821, t ≥ 8.97·10^17 | 2.5–2.9 (φ = 1/3) | void | **holds: 0.48207** |
| literal-hypotheses variant 0.5782 | 2.5–2.9 within stated hypotheses | void | superseded (no t ≥ 5·10^6 hypothesis) |
| Explicit611 Parts A, B (t < 1.577·10^6) | none (A: Arb ζ) / RS form (B) | hold (0.56693, 0.56329) | hold |
| 2602.05614's own 0.611 / 0.566 | all | void | HPY's 0.618 |
| HPY 0.618 (refereed), sub-Weyl 54.6 t^{27/164} (uses 0.618 on [3, e^75]) | none | hold | hold |
| citing works 2603.16023, 2607.10053, 2609.22624; Harvest3 candidate 5 | 0.611 | void | can use 0.568 |

The variant's bounds dominate the original ones pointwise (X_m gains 1 + 2/π; β_H ≥ β_R because cW^{1/3} ≥ 0.3·5π^{1/3} > 0.56 = (6/π−1)/(1+2/π)), so under the hypotheses the certificate asserts it rests on either lemma set. Suggested for the orchestrator: replace the record's caveat "rests on Lemmas 2.5–2.9 of 2602.05614 (unrefereed)" by "holds with HPY's refereed lemmas (0.56797; t ≥ 10^10: 0.5601); Revers's lemmas audited correct (RevAudit)".

## 7. Referee verdict (plain words)

The paper's new estimates (the improved second-derivative test and its consequences) are right. Every step of their proofs checks out; the paper skips some routine cases and has a handful of typos, but nothing it claims is false, and no constant is rounded the wrong way. Its Riemann–Siegel error term is exactly Gabcke's published one. Our use of the lemmas outside one stated condition (heights between 1.58 and 5 million) is justified, because that condition only serves to make the block range non-empty. One small step the paper takes for granted (the size of the first correction term) is now proved by computer. Most importantly, our record 0.568 does not actually need this unrefereed paper: redoing the certificate with the older, refereed lemmas of Hiary, Patel and Yang still passes at 0.568 (worst case 0.56797).

## Files
- `falsify.py` → `logs/falsify.log`: falsification suite (command in its header).
- `rs_c1.py` → `logs/rs_c1.log`: maxima of |C₀|, |C₁|; `c0_bound.py` → `logs/c0_bound.log`: rigorous |C₀| < 0.924058.
- `regions_check.py`, `regions_full.py` → `logs/regions_*.log`: the paper's 10 regions, printed vs full-precision coefficients.
- `make_hpy_variant.py` → `logs/run_hpy_*.log`: refereed-only variant of the Explicit611 certificate.

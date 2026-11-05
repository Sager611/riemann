RECORD: Robin verified for n <= 10^(10^38.60); Nicolas for p_k <= 1.04*10^40 (node N4)

# Robin and Nicolas: verified ranges from RH to 3·10^12 and the 7/8 half-plane

Agent robinh, 2026-11-05.

Files in this folder:
- PROOF.md: all lemmas.
- certify.py: Level-1 certificate (sharp cutoff).
- certify2.py: Level-2 certificate (smoothed cutoff).
- check_all.py: reruns every inequality in 173 s on one core.
- sanity.py, zerocheck.py: validation only.

## 1. Short answer

1. **Robin.** σ(n) < e^γ n log log n holds for every 5041 ≤ n ≤ 10^(10^38.60). Equivalently, it holds for every colossally abundant (CA) number whose largest prime is at most 9.37·10^38.
   - Inputs: RH up to height 3·10^12, no zeros with Re s > 7/8, and an explicit zero count N(T).
   - Previous best found: Polak (Zenodo 21808589, 5 Aug 2026, not refereed), n ≤ 10^(7.1·10^22) = 10^(10^22.85), largest prime 1.65·10^23.
   - Previous refereed best: Morrill–Platt, n ≤ 10^(10^13.11485), largest prime 3.0·10^13. TREE.md N4 still lists this one.
2. **Nicolas.** N_k/φ(N_k) > e^γ log log N_k holds for every primorial with p_k ≤ 1.04·10^40.
   - Previous: Nicolas 2012 (arXiv 1202.0729) computed c(N_k) for p_k ≤ 10^9 in Maple with 30 digits, without interval arithmetic.
   - No height-tied range was found anywhere.
3. **Why it works.** Both inequalities reduce to two-sided bounds on one tail integral:
   −I(x) = ∫_x^∞ (t − θ(t)) w(t) dt, with w = (1 + log t)/(t² log² t).
   - Its explicit formula has weights about 1/(ρ(ρ−1)), against a margin of about 1/(√x log x).
   - Level 1, sharp cutoff: zeros above T cost x^(β−1)/γ². With β ≤ 7/8 this reaches about 10^31. With no half-plane it reaches about 10^23, which matches Polak.
   - Level 2, cutoff smoothed over a window of log-width λ ≈ 10^(−9): zeros above T then cost x^(β−1)/(λ γ³). The window costs √x·η·λ, where η ≈ 4·10^(−11) is our own pointwise bound on abs(θ(x)/x − 1).
4. **The main lever is the smoothing, not the half-plane.** Level 2 with no half-plane already gives Robin to 10^(10^32.98) (largest prime 2.18·10^33), 10^10 times Polak's prime range. The 7/8 half-plane then adds a factor of 4·10^5 in X.
5. **Variant Re ρ ≤ 0.87338196661320752** (Bebel; its axioms are unchecked): Robin to 10^(10^38.71) (largest prime 1.19·10^39), Nicolas to p_k ≤ 1.31·10^40.

## 2. Ranges

All entries for this work are certified by check_all.py in ball arithmetic.

| method | zero input above T = 3·10^12 | Nicolas p_k ≤ | Robin: largest prime ≤ | Robin: all n ≤ |
|---|---|---|---|---|
| Level 2 (smoothed) | Re ρ ≤ 7/8 | **1.04·10^40** | **9.37·10^38** | **10^(10^38.61)** |
| Level 2 | Re ρ ≤ 0.87338… | 1.31·10^40 | 1.19·10^39 | 10^(10^38.71) |
| Level 2 | none (β < 1) | 1.42·10^34 | 2.18·10^33 | 10^(10^32.98) |
| Level 1 (sharp) | Re ρ ≤ 7/8 | 1.93·10^32 | 1.56·10^31 | 10^(10^30.83) |
| Level 1 | Re ρ ≤ 0.87338… | 2.67·10^32 | 2.13·10^31 | 10^(10^30.97) |
| Level 1 | none | 1.61·10^24 | 2.45·10^23 | 10^(10^23.03) |
| Polak 2026 (sharp, no half-plane) | – | – | 1.65·10^23 | 10^(10^22.85) |
| Morrill–Platt 2018 (Briggs' algorithm) | – | – | 3.0·10^13 | 10^(10^13.11) |
| Nicolas 2012 (Maple, floating point) | – | 10^9 | – | – |

Lower ends and overlaps:
- Nicolas: the 664,579 primorials with p_k ≤ 9,999,991 are checked directly in arb. The analytic range is 10^7 ≤ x ≤ X_N.
- Robin: the analytic range starts at largest prime 10^10. CA numbers below that have log N < 1.01·10^10 and are covered by Morrill–Platt, so there is no gap.

Tradeoffs:

| | Level 2 (this work) | Level 1 (this work) | Polak 2026 |
|---|---|---|---|
| cost of a zero above T | x^(β−1)(2/λ)/γ³ | x^(β−1)/γ² | x^(β−1)/γ² |
| extra input needed | pointwise abs(θ − t) ≤ ηt on the window (own Lemma 10_k) | none for Nicolas; Lemma 10 for Robin's second-order term | Büthe's partial-RH θ bound (x ≤ 2·10^25) |
| what limits X | window √x ηλ against zeros x^(b0−1/2) S3/λ: X ~ (S3 η)^(−1/b0) | zeros: x^(3/8) ~ T/log T | zeros: x^(1/2) ~ T/log T |
| status | reviewed (SOUND), certified | reviewed (one constant fixed), certified | Zenodo preprint |

## 3. Inputs

| input | statement | source (theorem) |
|---|---|---|
| RH to T | every zero with 0 < γ ≤ 3·10^12 has β = 1/2 | Platt–Trudgian, Bull. LMS 53 (2021), arXiv 2004.09765, main theorem |
| half-plane | no zero with Re s > 7/8 (so, by symmetry, none with Re s < 1/8) | OpenAI 2026, Lean-checked port in Riemann/Deep/QRH; variant 0.87338… (Bebel, axioms unchecked) |
| N(T) | abs(N(t) − (t/2π) log(t/2πe)) ≤ 0.1038 log t + 0.2573 log log t + 9.3675 for t ≥ e | Hasanalizade–Shen–Wong, JNT 235 (2022), arXiv 2107.06506, Cor. 1.2 |
| Σ 1/(ρ(1−ρ)) | = 2 + γ − log 4π = 0.0461914… | classical (Davenport Ch. 12) |
| ψ_1 formula | ψ_1(x) = x²/2 − Σ x^(ρ+1)/(ρ(ρ+1)) − x log 2π + c − Σ_r x^(1−2r)/(2r(2r−1)) | Ingham 1932, Thm 28 |
| θ, ψ ≤ 10^19 | abs(x − ψ) ≤ 0.94√x; 0.05√x < x − θ ≤ 1.95√x (1423 ≤ x ≤ 10^19) | Büthe, Math. Comp. 87 (2018), arXiv 1511.02032, Thm 2 |
| Chebyshev | θ(y) < y log 4 | Hardy–Wright Thm 415 |
| Mertens identity | Σ_{p≤x} 1/p = log log x + B + (θ(x)−x)/(x log x) − ∫_x^∞ (θ−t) w | Rosser–Schoenfeld 1962, (4.13)–(4.15) |
| small Robin range | Robin for 5041 ≤ n ≤ 10^(10^13.11485) | Morrill–Platt, arXiv 1809.10813, Thm "Platt" |

## 4. The explicit-formula lemma (full proofs in PROOF.md)

**Reduction** (Lemmas 1–4):
- Nicolas at p_k = x follows from −I(x) > S_2(x) ≤ 1/(2(x−1)). This needs no information about θ(x) itself, because the second-order term D(x) is ≥ 0.
- Robin at the CA number with parameter ε (x = x_1(ε), x_2 = x_2(ε)) follows from
  −I(x) + D(x) < Σ_{x_2<p≤x} p^(−2) + log(1 + log(1 + θ(x_2)/θ(x))/log θ(x)).
  Here D(x) ≤ η²/log x·(1 + …) needs a pointwise bound for θ(x).
- Under RH the normalised margins are 2 − 0.046 (Nicolas) and 2√2 − 2 − 0.046 = 0.78 (Robin).

**Lemma 5 (exact formula).** −I(x) = M(x) + Σ_ρ ρ^(−1) ∫_x^∞ t^(ρ−2) h dt + log(2π)/(x log x) + V(x), where:
- M = ∫_x^∞ (ψ − θ) w dt;
- V ∈ [−x^(−3)/log x, 0] is the trivial-zero term;
- the zero series converges absolutely. This follows from Ingham's ψ_1 formula by integration by parts and dominated convergence.

**Lemmas 6–7 (sharp cutoff).** The zero sum is at most
x^(−1/2) h(x) · 0.0461914 · (1 + κ/14.1347) + ½(x^(b0−1) + x^(−b0)) h(x)(1 + κ/T) · 2.95941·10^(−12).
The low zeros use RH below T. The high zeros use the pairing ρ ↔ 1 − ρ̄, which is maximal at β = b0 by convexity.

**Lemma 11 (ramp).** Split 1_{t≥x} = φ + 1_{[x, xe^λ]}(1 − φ), with φ linear in log t. Then:
- The window term satisfies abs(W) ≤ η H(log x) λ/2.
- M_φ ∈ [M(xe^λ), M(x)].
- Each zero satisfies abs(F_φ(ρ)) ≤ x^(β−1)(2H/λ + 3abs(H′))/abs(ρ−1)², by two integrations by parts.
- The low zeros keep the sharp bound plus λ/2.

**Lemma 10_k (pointwise).** Take k-fold one-sided box averages of ψ, where ψ is monotone and the kernel has mass h^k and mean kh/2:
abs(ψ(x) − x) ≤ kδx/2 + log 2π + x^(−2) + x^(1/2)(…) + ½(x^(b0) + x^(1−b0))(1 + kδ)^(b0+k)(2/δ)^k S_{k+1}(T).
It gives η(10^30) = 5.6·10^(−11) and η(10^38) = 4.2·10^(−11) (best k ≤ 5 and δ on a grid).

Bounds on n(x)·(−I(x)), n = √x log x, at Level 1 with 7/8 (check_all.py part 3):

| x | Level 1 | zeros above T contribute ≤ |
|---|---|---|
| 10^10 | [1.8652, 2.0563] | 8.7·10^(−9) |
| 10^20 | [1.9092, 2.0114] | 4.8·10^(−5) |
| 10^25 | [1.9134, 2.0216] | 3.6·10^(−3) |
| 10^30 | [1.6552, 2.2899] | 0.267 |

At Level 2, with Re ρ ≤ 7/8, the certified normalised margins at x = 10^35 are 1.91 (Nicolas) and 0.74 (Robin). At 10^38.5 they are 1.48 and 0.27.

## 5. Scripts

- **check_all.py**: uv run --no-project --with python-flint --with mpmath --with numpy python3 check_all.py. It prints ALL PASS in 173 s on one core. It runs:
  1. the direct Nicolas check to 10^7;
  2. Level 1 on the stated ranges for the three zero inputs;
  3. the bound table;
  4. Level 2 on the stated ranges for the three zero inputs.
- **certify.py, certify2.py**: the libraries behind check_all.py. Run directly, they search for the maximal endpoints:
  - Level 1: 10^32.287333 and 10^31.194615 (7/8); 10^32.427126 and 10^31.328757 (0.87338…); 10^24.208989 and 10^23.389562 (none).
  - Level 2: 10^40.017382 and 10^38.971833 (7/8); 10^40.117513 and 10^39.075885 (0.87338…); 10^34.154691 and 10^33.338478 (none).
  - How they work: adaptive boxes in log x (width 0.25 down to 10^(−6)), arb balls at 192 bits, x_2 enclosed by certified bisection on F(·,2), and λ, k, δ chosen per box.
- **sanity.py, zerocheck.py**: validation only (section 6).

## 6. Self-review

1. **Independent reviewers** (read-only subagents, two runs):
   - Sections 1–4: every step was OK except one wrong constant. Lemma 5 stated abs(w′) ≤ 3/(t³ log t), which is false near t = 2; the correct bound is 11/(t³ log t). The constant appears only in the dominated-convergence step and changes no number. Fixed.
   - Section 5 (Lemmas 11 and 10_k): VERDICT SOUND. The reviewer rechecked both integrations by parts, the jump terms, the bound ∫abs(G″) ≤ 3abs(H′), and the kernel's mass and mean.
   - A third pass compared the code with the formulas. Items 1–7 matched. It found two edge cases in certify2.py, both fixed:
     (a) the (1 + η_1) factor for x > 10^38 was switched on by the box, not by the box plus its window;
     (b) λ could reach 2·10^(−3), above the stated 10^(−3).
     Both are numerically negligible (about 10^(−27) and never triggered). After the fix check_all.py printed ALL PASS again with the same numbers.
2. **Per-zero bounds tested numerically** (mpmath, 48 cases: β ∈ {0.5, 0.875, 0.97}, γ from 14 to 3·10^6, log x ∈ {23, 70}, λ ∈ {10^(−2), 10^(−4)}). The worst ratios of true value to bound are 0.99999966 (Lemma 6, tight as γ → ∞), 0.951 (B1) and 0.959 (B2).
3. **Signs and constants against real zeros** (zerocheck.py). At x = 10^5, 10^6 and 3·10^6:
   - I computed the zero sum from primes alone, as −I − M − log(2π)/(x log x).
   - Summed over the first 300 zeros, the explicit formula matches it: +0.00935 vs +0.00620, +0.01117 vs +0.00988, −0.00357 vs −0.00408 (normalised). These are within the truncation error.
   - A flipped log 2π sign would shift the comparison by 0.0116 at 10^5.
4. **Bounds against truth** (sanity.py):
   - The true n(x)·(−I(x)) is 2.029, 2.018, 2.004 and 2.012 at 2.1·10^6, 3·10^6, 5·10^6 and 10^7. It lies inside both the Level-1 and the Level-2 intervals, for example [1.80, 2.32] at 10^7.
   - True Robin margins at CA numbers with largest prime 10^5–10^7 are 0.78–0.82, and the certified bounds are on the safe side of the truth.
5. **Smoothing errors.**
   - Lemma 5 is exact.
   - Lemma 11's window term is bounded with the pointwise η on [x, xe^λ]. Büthe's 1.95/√x is used only when xe^λ ≤ 10^19.
   - Lemma 10_k's bias kδx/2 is counted in full.
6. **Trivial zeros.** They give V ∈ [−x^(−3)/log x, 0], and in Lemma 10_k a term of size at most x^(−2).
7. **Zeros near T.** No sum is truncated. Zeros with abs(γ) ≤ T are on the line, and every other zero is bounded per zero with β ∈ [1/8, 7/8]. In Lemma 10_k, zeros with T < abs(γ) ≤ U are majorised by the high-zero term.

## 7. Novelty pass

- **arXiv math.NT listings** 2026-09 (812 entries), 2026-10 (854) and 2026-11 (140): titles scanned by keyword. I read the abstracts of 2609.26787, 2609.33794, 2610.12816, 2610.39120, 2610.18491 and 2609.23222 (Fiori–Jaskari). None gives a Robin or Nicolas range, and none uses the half-plane for them.
- **arXiv 2608.01498** (Mertens sums, Aug 2026) uses the same identity and an explicit formula for Σ 1/p. Its bounds are of shape (log x)^(−ℓ), and it treats neither Nicolas nor Robin.
- **arXiv search** (web_search, arXiv mode) found nothing new.
- The evidence API returned 429 (quota spent), so it was skipped.
- **Zenodo**, phrase searches for Robin, Nicolas, colossally abundant and quasi-Riemann:
  - **Polak 21808589** (Aug 2026) is our Level 1 without the half-plane: the same N(T) bound, the same pairing and a sharp cutoff, reaching 1.65·10^23. Its §10 lists "raise the height, sharpen off-line-zero estimates, or change the smoothing kernel" as routes; Level 2 and the 7/8 half-plane are the second and third of these.
  - **Fabbian 23782485** (Oct 2026) proves F(x) := log P(x) − γ − log log θ(x) ≤ 1.3064373·10^(−8) for x ≥ 3·10^13 (uniform, no range, no half-plane). Our bound at 3·10^13, 1.2·10^(−8), is consistent with it.
  - Other 2026 deposits are claimed proofs or numerics.
- **Own files.**
  - R316 missed Polak, and also said Büthe's height-tied bound alone does not give Nicolas because of the tail; the integrated formula handles the tail.
  - R294's "7/8 extends neither range" holds only for the pointwise Fiori–Jaskari route.
- **Verdict.**
  - New: the first use of the 2026 half-plane for Robin or Nicolas, the first height-tied Nicolas range, and the smoothed-cutoff method (Level 2).
  - Not new: Level 1 without the half-plane, which is essentially Polak (our 2.45·10^23 against his 1.65·10^23).

Suggested lines for the orchestrator:
- TREE.md N4 record: "Robin verified for 5041 ≤ n ≤ 10^(10^38.60), Nicolas for p_k ≤ 1.04·10^40 (RH to 3·10^12 + Re ρ ≤ 7/8; smoothed tail cutoff; robinh); without the half-plane 10^(10^32.98) and 1.42·10^34; before: Polak (Zenodo 21808589, 2026) 10^(7.1·10^22), Morrill–Platt 10^(10^13.11485)".
- DISCOVERIES.md: "Robin's inequality verified to 10^(10^38.6) and Nicolas' to p_k ≤ 10^40 from RH to 3·10^12 and the 7/8 half-plane (N4)."
- NOGO.md: none. The R294 and R316 dead-end lines are superseded for this method.

## 8. Caveats

- The half-plane is a 2026 result, Lean-checked in Riemann/Deep/QRH. The 0.87338 variant has unchecked axioms.
- Lemmas 10_k and 11 are our own and not published. Two reviewer passes and the numeric tests above support them. Lemma 10_k is needed for Robin above 10^19 and for the Level-2 window.
- The certificates trust python-flint 0.9 (arb). E_1 is computed as −Ei(−z).
- Polak's deposit is not refereed, so the refereed baseline is still Morrill–Platt.

RH relevance: N4

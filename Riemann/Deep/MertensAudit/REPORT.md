AUDIT PASSED: every step of the Mertens power-saving chain (Mertens78 PROOF.md, MertensCH DERIVATION.md, HP141's θ₁₄₁ rerun) holds as written against its fetched sources, the record's scripts certify exactly the stated inequalities with no gap, and an independent mpmath recomputation reproduces |M(x)| ≤ 9.69·10^4 x^0.9382 (96 825.485 before rounding up) and ≤ 1.52·10^4 x^0.9409 (15 192.689) for all x ≥ 1.

# MertensAudit — line-by-line referee of our record mertens-power-saving (node N10)

Agent mertensaudit, 2026-11-19 23:21 → 2026-11-20 UTC. Read-only on Riemann/Deep/Mertens78/, MertensCH/, HP141/, RECORDS.md, discoveries/. The whole chain with source line numbers: [FULLPROOF.md](FULLPROOF.md) (Parts I–V). M78 = Mertens78/PROOF.md, MCH = MertensCH/DERIVATION.md.

## Claim audited (RECORDS.md l.11)

The record: |M(x)| ≤ 9.69·10^4 x^0.9382 and ≤ 1.52·10^4 x^0.9409 for all x ≥ 1 (θ₁₄₁ since 14 Nov; was 9.77·10^4, 1.54·10^4; Mertens78's 1.84·10^6 x^0.9382 the first, 5 Nov). Certified by HP141/MertensCH/cert_ch.py --bp (log cert_ch_bp.log): Chirre–Helfgott on [1, X₀], 1817 (θ = 0.9382) / 1737 (θ = 0.9409) hybrid ranges on [X₀, X₁], final range on [X₁, ∞) with certificate ML092g8_e05_Bp; ranges use ML092g8_e05_Bp and ML093g5_e10_Bp only.

## Verdict per item of the brief

| # | item | verdict | evidence |
|---|---|---|---|
| 1 | Rademacher's convexity bound (M78 (I2)) | HOLDS | Kadiri–Lumley–Ng 2101.12263 TeX l.396–403 (quoting Rademacher Thm 4): for −1/2 ≤ −η ≤ σ ≤ 1+η ≤ 3/2, all t, \|ζ(s)\| ≤ 3\|1+s\|/\|1−s\|·(\|1+s\|/2π)^{(1−σ+η)/2}ζ(1+η) = (I2) verbatim. Used with η ∈ {0.1, 0.05, …, 0.005} and η′ = 0.1. Numerically (rademacher.py, t ≤ 10^4): ratio ≤ 0.29 everywhere, → 1 only at s = −1/2, η = 1/2, where the bound is an equality (ζ(−1/2) = −ζ(3/2)/4π) |
| 2 | harmonic majorant on the whole half-plane (A1; H1's B = 1/2 version; B′ version) | HOLDS | both σ-ranges re-derived; v = h − u ≥ 0 harmonic on all of Re s > B (FULLPROOF II) |
| 3 | Herglotz + signed Poisson certificates (A2–A4: positivity of v, linear term, negative weights) | HOLDS | re-derived; reviewer CORRECT; Lemma A2 with true ζ values: max of v(s₁) − Σλ_jv(s_j) over 160 t = −4.05 < 0; (C1) re-proved by my own mpmath.iv cover (true minimum of π(LHS − RHS) 1.5·10^−4, (C2) margin 0.39) |
| 4 | zero-free input: 7/8 (M78), θ₁₄₁ = p/q of R466 l.10–12 as B′ | HOLDS | B′ = 0.87334585610671858 ≥ θ₁₄₁ (exact Fraction: +3.2·10^−18); heights σ_j′ − B′ = σ_j − 7/8 exactly, so (C1)/(C2) are literally unchanged; κ, A recomputed with B′ |
| 5 | smoothed Perron and contour shift (B1, B2, H2) | HOLDS | no residue (1/ζ(1) = 0; 0, −1 outside); horizontal tails O(T^{a′−2}), a′ ≤ 0.759 (7/8), ≤ 0.776 (B′) |
| 6 | unsmoothing h = c_h x^θ and monotonicity (B3–B6, H0, H4, H5) | HOLDS | each factor at its worst endpoint; σ_L ≤ θ < 1, σ_H ≤ θ, θ ≥ θ* asserted in Arb |
| 7 | hybrid contour with RH to 3·10^12 (H1–H3) | HOLDS | Platt–Trudgian 2004.09765 TeX l.61–62, Thm 1: RH up to height 3 000 175 332 800 ≥ T′ + 0.52 (the Harnack discs); the abstract's rounded 3·10^12 (l.46) alone would not cover them |
| 8 | Chirre–Helfgott \|M(x)\| ≤ 3x/(π10^10) + 11.39√x on [1, X₀] | HOLDS (statement level) | cited from CHAudit2 (AUDIT PASSED; the record is robust to the loss of CH's square-free chain and to Platt's sum at 10^10+1 up to ≈ 7,690); still open there: Platt's private residue sums |
| 9 | junctions | HOLDS | M78: \|M(x)\| ≤ x on [1, x_a] ⟶ max(C, x_a^{1−θ}); MCH: CH on [1, X₀] (convex in log x ⟹ endpoints), ranges contiguous in exact exponents n/10^6, X₁ = last end, final [X₁, ∞) |
| 10 | code: outward rounding, coverage, thresholds | HOLDS | FULLPROOF IV: every inequality of the proof has its Arb check; 0 of 142 312 balls of the record's (C1) and J_L0 checks miss an endpoint of their piece; every scan starts at or below its threshold |
| 11 | constants | HOLDS | Arb maxima 96 825.4855 / 15 192.6885 rounded up to 9.69·10^4 / 1.52·10^4 and re-compared; mine (evalind.py, mpmath 40 digits, no flint, no import of audited code) 96 825.485 / 15 192.689; every M78 table constant is rounded up (msanity.log) |
| 12 | falsification of \|1/ζ(σ₁+it)\| ≤ C\|t\|^a | no counterexample | largest ratio 8.8·10^−10 (final certificate, Cor. A4, t ≤ 10^4), 1.25·10^−9 (its Φ(t), t ≤ 10^4), 4.7·10^−7 (ML093g5), 2.65·10^−5 (RH strip, Lemma H1), ≤ 6·10^−9 at t = 3·10^12; peaks at zero ordinates (14.12, 1083.2) |
| 13 | exact M(x) | HOLDS (vacuous) | own sieve to 10^7; OEIS A084237 M(10^k), k ≤ 23; Hurst 1610.08551 TeX l.479–512 extremes (\|M\| ≤ 40 371 499, q ∈ [−0.525, 0.571] to 10^16): largest \|M\|/bound 1.0·10^−5 (x = 1) |

No step fails, so nothing needs to be weakened: the record stands as written.

## Observations (no number changes)

1. **Hand-over wording.** TREE.md l.14 and discoveries/mertens-power-saving.md l.21 say "Chirre–Helfgott below 10^242.85". With θ₁₄₁ the certified run hands over at X₀ = 10^242.80 (θ = 0.9382) and 10^240.25 (θ = 0.9409) (HP141/MertensCH/logs/cert_ch_bp.log); 10^242.85 was the B′ = 0.87338… run. Suggest "below 10^242.8".
2. **θ₁₄₁ decimal.** "0.8733458561067185768" is θ₁₄₁ truncated: the exact p/q is 0.87334585610671857680118…, so "Re s > 0.8733458561067185768" (TREE.md l.502, R466 title, the brief) overstates the Lean theorem by 1.2·10^−21 (TREE.md l.11 writes "…", fine). Every script uses a value rounded up (0.87334585610671858; 0.8733458561067185768011848), so no certified number is affected. Suggest "θ₁₄₁ = 0.87334585610671857680…".
3. cert_ch.py --bp builds the B′ certificates with bprime.make() and never executes bprime.run()'s "a′ < 2" assertion (Lemma B2's horizontal tails); a′ ≤ 0.776 analytically for every B′ certificate. Harmless.
4. cert.py builds its balls as arb((lo+hi)/2, (hi−lo)/2) from floats (l.80, l.121), which could in principle miss an ulp of a piece; 0 of the 142 312 balls behind the record do (the radius is rounded up to a 30-bit mag).
5. MCH Lemma H3's hypothesis x ≥ 2·10^8 is not used in its proof (asserted anyway). Rademacher's 1959 paper itself was not fetched (paywalled); the quotation in KLN matches (I2).
6. The key lemma is loose by about 10^9 at the final certificate (largest ratio 1.3·10^−9): the harmonic-majorant step costs e^{21.5}. A sharper explicit bound for 1/ζ on Re s ≈ 0.918 above height 3·10^12 would lower c directly (not pursued).
7. COROLLARIES.md l.24 (1.84·10^6) and l.25 (9.77·10^4: Liouville and Möbius-mean corollaries) still quote earlier constants: valid, weaker, not updated to 9.69·10^4.

## Dependents (rg -n -i 'mertens78|mertensch|mertens-power' discoveries/ TREE.md COROLLARIES.md RECORDS.md DISCOVERIES.md FORMALIZATIONS.md)

RECORDS.md l.11 (the record); TREE.md l.14 (N10 record line), l.372 (R336), l.388 (R352), l.464 (R428); COROLLARIES.md l.24 (Möbius series in Hardy spaces, explicit rates), l.25 (Liouville and Möbius-mean power savings, R443); discoveries/mertens-power-saving.md l.19, l.21; rounds Riemann/Tree/R336, R352, R379, R428, R437, R443. DISCOVERIES.md, FORMALIZATIONS.md: no entry. All stand.

## Reviewer

One reviewer subagent on the most delicate step, pasted inline (Lemma A2, Theorem A3, Corollary A4 with negative weights, the (C1) tail bracket): **VERDICT: CORRECT** — "for λ < 0, λ(h − u) = −|λ|h + |λ|u ≤ −|λ|(h − L1): the inequality direction is correct"; the linear term needs only c ≥ 0 and (C2), which (C1) already implies as |w| → ∞; both tail claims re-derived by differentiation. Its one remark (κ ≥ 0 is used) holds: κ = (1 + η − B)/2 > 0.

## Method

Read every cited input at its source (KLN, Platt–Trudgian, Hurst e-prints, OEIS A084237; θ₁₄₁'s exact p/q from R466.lean l.10–12); refereed M78 and MCH line by line (FULLPROOF.md II–III); read the record's scripts line by line (IV); ran (V): falsify.py (|1/ζ| against Φ(t), C_∞t^a and C_Lt^{a_L} on t-grids to 10^4, Lemma A2 with true ζ values), mpcheck.py (the peaks in mpmath, and t = 3·10^12), c1check.py (their ball cover; my own mpmath.iv proof of (C1)), dumpplan.py + evalind.py (their float plan, my mpmath evaluation of every range), msanity.py (exact M(x), rounding of every printed constant, B′ ≥ θ₁₄₁).

## Files (Riemann/Deep/MertensAudit/)

REPORT.md, FULLPROOF.md; dumpplan.py, evalind.py, falsify.py, mpcheck.py, c1check.py, msanity.py, rademacher.py; logs/ (dumpplan.log, evalind.log, falsify.log, mpcheck.log, c1check.log, msanity.log, rademacher.log, plans.json).
Footprint (memcap peak / max RSS): dumpplan 50 / 69 MB, evalind – / 33 MB, falsify 44 / 54 MB, c1check 39 / 56 MB, mpcheck 16 / 33 MB, rademacher 18 / 37 MB, msanity – / 168 MB ("–": done before memcap's first sample). A first msanity run (unchunked numpy, 9 s) reached 511 MB max RSS before memcap could sample it — over the 290 MB cap; rewritten in 10^6-chunks and rerun (logs/msanity.log is the rerun). One process at a time, single-threaded, background QoS kept.

RH relevance: N10


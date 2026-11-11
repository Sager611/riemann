RECORD: |M(x)| ≤ 9.77·10^4 · x^0.9382 for every x ≥ 1 (record: 1.84·10^6) and |M(x)| ≤ 1.54·10^4 · x^0.9409 (was 7.91·10^5); the constant drops at every θ in [0.8992, 0.9409]. Of the 18.8× at θ = 0.9382, 6.1× comes from RH verified to 3·10^12 on a hybrid contour; the rest is routine.

# Mertens bound with RH to 3·10^12 on a hybrid contour (node N10)

Agent mertensch, 2026-11-06 and 2026-11-11 (interrupted 4.5 days by an expired login). Proofs in DERIVATION.md; certificate: sh run.sh.

## Claim

For every θ in the table and every x ≥ 1, |M(x)| ≤ c x^θ, with c from column "new". Headline: |M(x)| ≤ 9.77·10^4 x^0.9382 and |M(x)| ≤ 1.54·10^4 x^0.9409 for all x ≥ 1.

| θ | record row | published combined | routine (no RH strip) | new, 7/8 certificates | new, + B′ certificates | gain vs record | gain vs routine |
|---|---|---|---|---|---|---|---|
| 0.9409 | 10^5.90 | 10^5.26 | 1.23·10^5 | 2.25·10^4 | **1.54·10^4** | 51.4x | 7.99x |
| 0.9400 | 10^6.26 | 10^5.50 | 2.08·10^5 | 4.20·10^4 | **2.84·10^4** | 64.7x | 7.32x |
| 0.9382 | 10^6.26 | 10^5.96 | 6.00·10^5 | 1.48·10^5 | **9.77·10^4** | 18.8x | 6.14x |
| 0.9370 | 10^8.27 | 10^6.27 | 1.22·10^6 | 3.46·10^5 | **2.25·10^5** | 822x | 5.42x |
| 0.9350 | 10^8.27 | 10^6.79 | 3.95·10^6 | 1.46·10^6 | **9.15·10^5** | 202x | 4.32x |
| 0.9325 | 10^8.27 | 10^7.43 | 1.72·10^7 | 9.21·10^6 | **5.52·10^6** | 33.5x | 3.12x |
| 0.9300 | 10^8.27 | 10^8.08 | 7.49·10^7 | 6.53·10^7 | **3.58·10^7** | 5.16x | 2.09x |
| 0.9296 | 10^8.27 | 10^8.18 | 9.47·10^7 | 9.10·10^7 | **4.88·10^7** | 3.79x | 1.94x |
| 0.9277 | 10^8.88 | 10^8.79 | 3.27·10^8 | 5.25·10^8 | **2.30·10^8** | 3.27x | 1.42x |
| 0.9250 | 10^12.55 | 10^9.96 | 3.44·10^9 | 7.76·10^9 | **3.24·10^9** | 1.09e+03x | 1.06x |
| 0.9200 | 10^12.55 | 10^12.12 | 5.64·10^11 | 1.14·10^12 | **5.32·10^11** | 6.65x | 1.06x |
| 0.9192 | 10^12.55 | 10^12.47 | 1.28·10^12 | 2.63·10^12 | **1.21·10^12** | 2.92x | 1.06x |
| 0.9150 | 10^20.52 | 10^17.18 | 1.04·10^16 | 1.25·10^17 | **1.01·10^16** | 3.31e+04x | 1.03x |
| 0.9121 | 10^20.52 | 10^20.44 | 2.66·10^19 | 2.19·10^20 | **2.59·10^19** | 12.9x | 1.03x |
| 0.9100 | 10^26.09 | 10^24.47 | 2.97·10^22 | 2.15·10^24 | **2.87·10^22** | 4.27e+03x | 1.03x |
| 0.9092 | 10^26.09 | 10^26.01 | 1.11·10^24 | 7.51·10^25 | **1.07·10^24** | 114x | 1.04x |
| 0.9050 | 10^110.83 | 10^78.35 | 2.28·10^63 | 8.22·10^77 | **2.09·10^63** | 3.2e+47x | 1.09x |
| 0.9024 | 10^110.83 | 10^110.75 | 3.79·10^97 | 3.38·10^110 | **3.31·10^97** | 2.02e+13x | 1.15x |
| 0.9000 | 10^382.00 | 10^314.11 | 2.06·10^207 | 1.45·10^310 | **1.54·10^207** | 6.54e+174x | 1.34x |
| 0.8992 | 10^382.00 | 10^381.91 | 1.24·10^275 | 9.10·10^376 | **8.43·10^274** | 1.19e+107x | 1.47x |

Columns. "record row": smallest all-x constant among Mertens78 rows with θ_row ≤ θ (the record). "published combined": sup over x of the minimum of the bounds already stated (trivial, Chirre–Helfgott, every Mertens78 row), floating point. "routine (no RH strip)": the same certified pipeline with T′ = 1000, i.e. Chirre–Helfgott below X₀, Mertens78 with re-optimised smoothing, squarefree unsmoothing and the B′ certificates, but no use of RH above height 1000 (logs/cert_ch_results_bp_T1000.json, certified). "new, 7/8": hybrid contour with the nine Mertens78 certificates. "new, + B′": also the nine certificates moved to the Lean-checked half-plane B′ = 0.87338196661320752 (bold: the claim). "gain vs routine" is the part due to RH to 3·10^12.

Routine rerun with B′ (brief's second item, Mertens78's own pipeline, logs/bprime_results.json): the same nine certificates give θ smaller by 0.0013–0.0014, at a constant 1.2 times (H093) to 10^14.8 times (S089g3000) larger, because ζ is larger on the lines moved closer to 1:

| certificate | θ (7/8) | C′ all x (7/8) | θ (B′) | C′ all x (B′) | C, x ≥ x_a (B′) |
|---|---|---|---|---|---|
| H093 | 0.9409 | 10^5.90 | 0.9396 | 10^5.97 | 10^5.86 (x ≥ 10^160) |
| ML093g5_e10 | 0.9382 | 10^6.26 | 0.9368 | 10^6.35 | 10^6.26 (x ≥ 10^158) |
| ML092g8_e10 | 0.9296 | 10^8.27 | 0.9282 | 10^8.48 | 10^8.39 (x ≥ 10^169) |
| ML092g8_e05 | 0.9277 | 10^8.88 | 0.9263 | 10^9.09 | 10^9.00 (x ≥ 10^171) |
| ML091g12_e03 | 0.9192 | 10^12.55 | 0.9178 | 10^12.83 | 10^12.75 (x ≥ 10^200) |
| S090g20_e01 | 0.9121 | 10^20.52 | 0.9108 | 10^21.13 | 10^21.04 (x ≥ 10^278) |
| S090g40_e005 | 0.9092 | 10^26.09 | 0.9078 | 10^26.95 | 10^26.86 (x ≥ 10^331) |
| S089g400_e005 | 0.9024 | 10^110.83 | 0.9011 | 10^114.89 | 10^114.81 (x ≥ 10^1198) |
| S089g3000_e005 | 0.8992 | 10^382.00 | 0.8978 | 10^396.76 | 10^396.69 (x ≥ 10^3918) |

(C1)/(C2) are unchanged by the translation; rechecked anyway for H093 (analytic), ML093g5 (18 678 balls) and ML092g8 (28 634 balls).

## Before

- Mertens78 (ours, 5 Nov 2026): |M(x)| ≤ 1.84·10^6 x^0.9382 for x ≥ 1, family 7.91·10^5 x^0.9409 … 1.01·10^382 x^0.8992; the first explicit power saving.
- Chirre–Helfgott (arXiv 2511.14736v1): |M(x)| ≤ 3x/(π·10^10) + 11.39√x for x ≥ 1; smaller than the record below 10^263.5.

## Method (DERIVATION.md)

1. Where CH's linear term comes from: their weight is band-limited (transform supported on |t| ≤ T); the term ι·x is the trivial bound |Σ_{N<n≤N+H} μ(n)| ≤ H on a window of relative width π/2T around x (their remark, TeX l.5148). Its meromorphic continuation coth(δ(s − σ)) is periodic in t, so the brief's integral of |w x^{s−1}/ζ| on Re s = B′, |t| > T diverges. Zeros above T can only be used after smoothing at a scale h ≪ x.
2. Hybrid contour: Mertens78's smoothed Perron integral is moved to Re s = 0.85 for |t| ≤ T′ = 3·10^12 (no zeros there: Platt–Trudgian verified RH to 3 000 175 332 800; 1/ζ bounded by Harnack in discs of radius 0.52 tangent to the critical line, with Rademacher's majorant on Re s > 1/2), along horizontal segments at ±T′, and to a certificate line Re s = σ_H ∈ [0.888, 0.93] for |t| ≥ T′ (|1/ζ| ≤ C_∞|t|^a from the half-plane Re s > 7/8, or Re s > B′ for the moved certificates).
3. Unsmoothing by squarefree density: |M(x) − average over [x − h/2, x + h/2]| ≤ 3h/(2π²) + 6√X (elementary |R(y)| ≤ 3√y) instead of h/4 + 1.
4. Assembly: Chirre–Helfgott for x ≤ X₀; on [X₀, X₁] ranges of width 10^0.05 with their own certificate and smoothing length; Mertens78's monotone argument on [X₁, ∞).

Why the RH part matters mostly for θ ≳ 0.93: at the hand-over from Chirre–Helfgott the optimal Perron cut-off τ = 2X/h is about 10^10 for every θ (because CH's error is about x/10^10), below T′. So almost all of the integral sits on Re s = 0.85, where it is negligible, and the bound is ≈ 3h/(2π²) + (C_∞/π)X^{σ_H}τT′^{a−1}/(1 − a). This lasts while τ < T′, i.e. while x^{1−θ} grows by a factor of a few hundred (for θ = 0.9382 up to x ≈ 10^340). For θ ≳ 0.93 the constant is decided inside this window; for smaller θ it is decided later (far ranges or the final range), where the half-plane alone sets it, as in Mertens78.

| ingredient | effect on c at θ = 0.9382 |
|---|---|
| Chirre–Helfgott below X₀ instead of |M(x)| ≤ x | 2.0× (record 1.84·10^6 → published combined 9.1·10^5) |
| hybrid contour (RH to 3·10^12) | 6.1× (routine 6.00·10^5 → 9.77·10^4) |
| squarefree unsmoothing | about 1.28× inside the hybrid window (c ∝ √β there; estimate, not run separately) |
| B′ certificates | 1.5× (new 7/8 1.48·10^5 → 9.77·10^4); inside the routine column also 1.5× (published combined 9.1·10^5 → 6.0·10^5 together with re-optimisation and squarefree unsmoothing) |

## Certificate

One command, from this folder: sh run.sh (single-threaded, one process at a time, about 40 min, below 100 MB). It runs
- bprime.py: the routine B′ rerun (logs/bprime_results.json);
- cert_ch.py: the 7/8 column (logs/cert_ch_results.json, log logs/cert_ch.log);
- cert_ch.py --bp: the claimed column (logs/cert_ch_results_bp.json, logs/cert_ch_bp.log);
- MCH_TP=1000 cert_ch.py --bp: the routine baseline without the RH strip (logs/cert_ch_results_bp_T1000.json);
- table.py: the table above (logs/table.md; floating point, comparison only).
Each cert_ch line prints θ, the certified c, the float plan, the record row, X₀, X₁, the number of ranges and the final certificate.
cert_ch.py imports Riemann/Deep/Mertens78/cert.py. Floating point only chooses X₀, X₁, the ranges and smoothing lengths; every term of Lemmas H4, H5 and the CH piece is evaluated in Arb (120 bits); c is the Arb maximum rounded up to 3 digits and re-compared in Arb.

## Checks

1. Crossover structure: CH against the record row 1.84·10^6 x^0.9382 crosses at x = 10^263.5 (compare.py), as in TREE N10. The new bound hands over from CH at X₀ = 10^242.85 (B′) or 10^245.75 (7/8) at θ = 0.9382; CH is smaller below X₀, the new bound above (by construction c ≥ CH(X₀)/X₀^θ).
2. Gain at fixed θ: θ = 0.9382, 9.77·10^4 against 1.84·10^6.
3. Arb against the float plan: the Arb constant exceeds the plan by at most 0.003 dex in every row (logs/cert_ch*.log).
4. The self-checks found three bugs, all fixed before the final run: the maximum of the pieces was chosen by float, which overflows above 10^308 (caught by the Arb rounding re-check); the planner underflowed at θ ≤ 0.9 (plan 10^322.8 against Arb 10^384.0; the Arb value was a valid but poor bound; fixed with logaddexp); the range grid stopped at 10^3000 (extended to 10^100000 for θ < 0.91). logs/prelim/ keeps the runs before the final one.
5. Hypotheses asserted in cert_ch.py and bprime.py: h ≤ x on every range, x ≥ 2·10^8, σ_L ≤ θ, σ_H ≤ θ and θ ≥ θ* in the final range, contiguous ranges from X₀ to X₁, T′ + 0.52 < T_PT, A′ > 0, 0 < a < 1, σ_H ≥ 0.89 (0.888 for B′, with the B2 exponent a′ < 2 asserted).
6. Every constant of DERIVATION.md appears in the scripts (rg, see below).
7. Footprint: one process at a time, single-threaded, at most 98 MB (top).
8. Reviewer subagent pass on DERIVATION.md and certify(): FLAWED for two missing hypotheses (θ < 1 in Lemmas H4–H5; no gap between X₀ and X₁ when there are no ranges); both now stated and asserted in the code, no result changed. Lemmas H0, H1, H3 and the implementation of H4–H5 were confirmed.
9. Independent float model (explore2.py, per-x optimal h, no ranges) gives 10^5.18 at θ = 0.9382 with 7/8 certificates against the certified 10^5.17.

## Novelty

NOVELTY.md: arXiv math.NT listings of Oct and Nov 2026, the arXiv API, HAL and Zenodo (6 and 11 Nov) contain no explicit power-saving bound for M(x) other than ours; 2511.14736 and 2512.15709 are still v1. The record improved is our own. The hybrid contour (RH up to T on a line near 1/2, a zero-free region above) is standard for ψ(x); it is new for M(x) only because no power-saving partner existed.

## Caveat

- Black boxes: Chirre–Helfgott's bound (unrefereed arXiv v1) for x ≤ X₀ (10^240 at θ = 0.9409 up to 10^3840 at θ = 0.8992); Platt–Trudgian (refereed); Mertens78's certificates and lemmas (reviewer: correct, 5 Nov).
- The gain is a constant, not an exponent: θ ≥ θ* of a half-plane certificate is still needed as x → ∞. For θ ≤ 0.9277 the part due to RH to 3·10^12 is at most 1.5×; the large gains over the record row there are routine (B′ certificates, re-optimised smoothing, Chirre–Helfgott for small x).
- The real size of M(x) is about √x; all of this is far from sharp.
- Only a subagent reviewer pass so far (item 8 of Checks); the formal reviewer verdict for DISCOVERIES is the orchestrator's step.

RH relevance: N10.

## Files

Riemann/Deep/MertensCH/: REPORT.md, DERIVATION.md, NOVELTY.md, run.sh (the one command), cert_ch.py (certificate), bprime.py (B′ translation), table.py, compare.py and bptable.py (comparison tables, floating point), explore2.py (first float exploration, untrusted), logs/ (final run of 2026-11-11 12:47–13:29 UTC: cert_ch*.log, cert_ch_results*.json, bprime_results.json, bprime.log, JL0.json, cert_*_Bp.json, table.md; logs/prelim/: earlier runs, identical constants).

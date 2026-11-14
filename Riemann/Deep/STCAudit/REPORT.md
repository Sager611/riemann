**Casaña-Giner Thm 2.1: AUDIT PASSED modulo unrerun certificates (Cζ ball sum and sup of the elementary part for all six lines, recomputed in float/mpmath; (K) for ℓ150, ℓ400, ℓ6000) · Ismail Thm 1.1: AUDIT FAILED (all nine C̃3 understated by a log(1+η)-for-log ζ(1+η) slip inherited from Hasanalizade–Shen–Wong 2022; rows 7–9 also use Patel–Yang's false 66.7)** — CG's kernel inequality (K) is certified by our own interval code for ℓ5, ℓ20 (smallest margin) and ℓ60, and the HSW slip also understates Bellotti–Wong 2025's published constants.

# STCAudit — referee report (stcaudit, 12–14 Nov 2026)

Papers: [CG] Casaña-Giner, Zenodo 23866074 (1 Nov 2026); [IS] Ismail, Zenodo 23301655 (11 Oct 2026). Details: [FULLPROOF.md](FULLPROOF.md) (§D = follow-up).

## Verdicts

| Claim | Verdict | Basis |
|---|---|---|
| CG Thm 2.1, (2.1) |N−M| and (2.2) |S|, T ≥ 1 | PASSED modulo the Cζ and sup-F Arb certificates and (K) for ℓ150, ℓ400, ℓ6000 | proof re-derived; (K) certified by own arb code for ℓ5, ℓ20, ℓ60; Cζ recomputed (≥ 7·10⁻⁴ below the certified bounds); sup F sampled (≤ 1.6·10⁻⁶ below); all 18 printed constants reproduce |
| CG Thm 2.3 (asymptotic) | not audited | outside the brief |
| IS Thm 1.1, all nine rows, (1.3), Cor 1.4 | FAILED as written | C1, C2 correct (own BW transcription, ≤ 5·10⁻⁹); C̃3 understated by +0.0275…+0.1401 (slip below); rows 7–9 also k1 = 66.7 |
| IS corrected | survives | headline |S(T)| ≤ 0.134791 log T + 0.159606 log log T + 1.534881 (T ≥ e); rows 1–6 C̃3 = 1.404841, 1.534881, 1.674012, 1.891768, 2.030092, 2.199195; rows 7–9 (68.54) 2.770550, 2.938158, 3.311666 |

## The slip (new; reviewer: CORRECT)
- HSW (arXiv 2107.06506, TeX l. 693–703) bound |ζ| on −η ≤ σ ≤ 0 by interpolating ζ(1+η)(|Q+s|/2π)^{½+η} (σ = −η) and (c1/√2π)|Q+s|^{½}(log|Q+s|)^{c2} (σ = 0), then take logs as −(σ/η) log(1+η) instead of −(σ/η) log ζ(1+η). F_{c,r} (HSW l. 788, BW l. 726) and D3 (HSW l. 884, BW l. 880, IS sbound.py) inherit it. As −σ/η ∈ [0,1] there and ζ(1+η) > 1/η > 1+η, every C3, C̃3 is understated by (log ζ(1+η) − log(1+η) − c2 log 1.00212)·∫_{θ0}^{θ_{−η}}(−σ/η)dθ / (2π log(r/(c−½))).
- Sizes: IS +0.0275 … +0.0356 (rows 1–6), +0.0376, +0.1401, +0.1336 (rows 7–9). BW Table 2: +0.0068 (headline 0.10076 log T + 0.24460 log log T + 8.08344 needs 8.09025), +0.0524, +0.0189, +1.1·10⁻⁵, +0.0593. HSW Table 2: +1.5·10⁻⁵ (headline), +0.0078, +0.0090.
- No erratum found (arXiv); BW's footnotes list two other HSW typos, not this one.

## CG (Turing-type method)
- Chain (Lemmas 3.1, 3.2, Prop 4.1, (4.7)–(4.9), Lemmas 4.2–4.10) re-derived; end-to-end with Odlyzko's zeros at 5 heights × 6 kernels (identity to 8·10⁻¹⁰; Prop 4.1 holds 29/29).
- KCERT: (K) proved with python-flint arb (80 bits), own code certify_K.py: centred-form bisection on [0,½] × (0,100] and [0,½] × [−100,0) with f in the form (4.17) (covers t → 0), and Lemma 4.10's tail t²U − 4d³u − K5|u|³ > 0 on u ∈ [−1/100, 1/100] (t → ∞). ℓ5 (4 580 + 2 boxes), ℓ20 (smallest margin 1.0·10⁻⁴; 165 026 + 2 046 boxes, 23 s), ℓ60 (551 926 + 8 190 boxes): CERTIFIED. ℓ150, ℓ400, ℓ6000 not closed in the 40-min box (ℓ6000: 3.1 M boxes on t > 0 in 800 s; d = 0.85 makes b = ½ steep): there (K) rests on the paper's certificate plus our sampling (margin ≥ 1.2·10⁻⁴).
- Not rerun: the Arb certificates for Cζ (our float recomputation over all 5.76 M primes sits 0.7–3.8·10⁻³ below each certified bound, the size of their derivative correction) and for sup F (sampled on their grid; certified bounds exceed the samples by ≈ 1.6·10⁻⁶).

## IS (Backlund–Bellotti–Wong method)
- Inputs refereed and in range (HPY 0.618 t^{1/6} log t; Hoo–Teo ½ log t + 0.6633; Yang 1.546; computed |S| ≤ 1, 2, 2.5167 below T0); the Q0 = 10⁹ shift device is correct (reviewer: CORRECT). Patel–Yang's 66.7 (rows 7–9) is false (REFUTATIONS; 68.54 adds +0.0020…+0.0027; our 54.6 would give 2.753463, 2.916964, 3.289242 with the slip fixed).
- Own transcription of BW's C1, C2, C̃3 reproduces IS's raw values to ≤ 5·10⁻⁹ (κ3 over 2T0 log T0, IS's correction of BW's misprint, is right); IS's 69-check script passes. The understatement is in the shared formula, not in IS's arithmetic.

## Falsification and consistency
- Exact S(T), N(T) − M(T) on [1, 74 920] (Odlyzko's 10⁵ zeros): max|S| = 1.6287; every CG line and every IS row (printed or corrected) holds; smallest slack CG ℓ5 0.407. No counterexample: the slips are in proofs, the bounds are 2–4× the truth here.
- CG's envelope lies below IS's at every T ≥ e (by ≥ 0.47; 3.899 vs 5.256 at 3.06·10¹⁰). Leading constants differ by method: Backlund's Jensen-disc average of ζ-growth exponents (≈ 0.1 floor now) vs Turing's exact identity, whose only log T is c_{R,1}/(2π), minimised by LP. CG uses no Jensen disc in Thm 2.1, so the HSW slip does not touch it.

## For the orchestrator
- Robin2Adic/cert.py and Ingham4/cert.py use HSW's 9.3675; BW (2025, footnote, TeX l. 118) says that value is a typo for 9.4925, and the slip adds ≈ 1.5·10⁻⁵ more: re-run both with ≥ 9.4926, or with CG's 0.08382 log T + 5.326 (T ≥ 1) if CG is accepted.
- Dependents: IS rows 7–9 and BW Table 2 depend on the Patel–Yang refutation (count once); the HSW slip is a new refutation-type item (HSW → BW → IS).
- Files (Riemann/Deep/STCAudit/): certify_K.py, certify_all.sh, certify_K_*.log (KCERT); bw_constants_own.py/.log; hsw_slip_size.py/.log; earlier scripts and logs as listed in FULLPROOF.md. All processes ≤ 214 MB (memcap).

RH relevance: N7, N9

## Dependents rerun with 9.4926 (14 Nov 2026)
Copies in /tmp/stc2, 9.3675 → 9.4926 (≥ 9.4925 + 1.5·10⁻⁵) in every .py (Li200k: r(93675,10000) → r(94926,10000)); each run under memcap 290 MB; machine load ≈ 15, so wall times are 2–4× the logged ones. Logs: dep_*.log here; drivers dep_*.py, dep_rerun*.sh.

| Script (logged run) | Rerun | Result |
|---|---|---|
| Ingham3/cert.py (193 s, 199 MB) | 445 s, 186 MB; cert.log and cert_result.json differ from the logged ones only in timings | unchanged: k = 7 every n ≥ 1; x0 = e^217.45, e^288.73, e^468.84, e^1881.03 (k = 6…3) and the δ-variant rows identical |
| Li200k/verify.py (39 s, 36 MB) | 148 s, 22 MB; E(T0) = 14.213738 confirms the new constant | unchanged: N = 708,505,587,700,793,981,366,435,840 (a = ½), 717,734,280,405,823,147,676,270,592 (MTY), 952,080,909,647,480,761,027,133,440 (7/8) |
| RobinHeight/certify.py (library of check_all.py; no logged run) | record claims of check_all.py (2), (4) via dep_robinheight_record.py: 162 s, 20 MB, all PASS | unchanged: Robin to 10^(10^38.6095) (7/8) and 10^(10^38.7133) (0.87338…), Nicolas to 1.04·10⁴⁰, 1.31·10⁴⁰; Level 1 ranges PASS |
| RobinCH/certify.py (identical copy, library of certify3.py; check_all2.py 460 s, 80 MB) | record claims of check_all2.py (2), (3) via dep_robinch_record.py: 173 s, 21 MB, all PASS (first check_all2.py attempt failed: numpy missing from my uv line) | unchanged: Robin to 10^(10^46.7931), Nicolas p_k ≤ 1.78·10⁴⁸ (b0 = 0.87338…); Level 2-CH 10^(10^40.5656) |
- Not rerun: RobinHeight/check_all.py and RobinCH/check_all2.py in full (their other parts are the direct Nicolas check to 10⁷, which does not use N(T), and the b0 = 1 rows); Ingham3/check.py (float cross-check, patched but not run).

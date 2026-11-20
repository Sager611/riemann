**AUDIT PASSED (CG Thm 2.1 complete)** — our own Arb certificates close every piece stcaudit left open for all six lines: (K) for ℓ150, ℓ400, ℓ6000 (and ℓ5–ℓ60 again) by two independent methods, Cζ as an outward-rounded ball sum over all 5 761 455 primes < 10⁸, and sup F by interval bisection; every printed constant lies above our rigorous value, so |S(T)| ≤ 0.08382 log T + 4.444 and |N(T) − M(T)| ≤ 0.08382 log T + 5.326 (T ≥ 1) are audited.

# STC2 — completion of the audit of Casaña-Giner, Zenodo 23866074, Thm 2.1 (stc2, 20 Nov 2026)

Builds on [STCAudit](../STCAudit/REPORT.md) (proof chain re-derived there; this round reruns the certificates). Python-flint 0.9 (Arb), own code, one process at a time, peaks ≤ 31 MB.

## Results (summary.py → ALL CHECKS PASS)

| line | (K): min max(1,τ²)H on b = ½ / min U(b,0) − π / tail Q (certified lower bounds) | (K) 2-D | Cζ own ≤ | printed | sup F (N) own ≤ | printed | sup F (S) own ≤ | printed | A3ᴺ own → printed | A3ˢ own → printed |
|---|---|---|---|---|---|---|---|---|---|---|
| ℓ5 | 1.859e-4 / 1.987e-4 / 1.08 | ✓ | 0.70290657 | 0.70360002 | 0.66536271 | 0.66536440 | 0.13948928 | 0.13949090 | 1.36827 → 1.369 | 0.84240 → 0.844 |
| ℓ20 | 1.006e-4 / 1.992e-4 / 2.26e-2 | ✓ | 1.46670482 | 1.46807534 | 0.70206889 | 0.70207071 | −0.18278804 | −0.18278649 | 2.16877 → 2.171 | 1.28392 → 1.286 |
| ℓ60 | 1.067e-4 / 1.989e-4 / 7.08e-3 | ✓ | 2.10254945 | 2.10426949 | 0.72893981 | 0.72894156 | −0.15847731 | −0.15847575 | 2.83149 → 2.834 | 1.94407 → 1.946 |
| ℓ150 | 1.424e-4 / 1.988e-4 / 6.25e-3 | ✓ | 2.25386939 | 2.25578346 | 0.73031522 | 0.73031685 | −0.15757773 | −0.15757634 | 2.98418 → 2.987 | 2.09629 → 2.099 |
| ℓ400 | 1.298e-4 / 1.996e-4 / 1.87e-4 | ✓ | 4.46007545 | 4.46374163 | 0.73756357 | 0.73756551 | −0.14582792 | −0.14582673 | 5.19764 → 5.202 | 4.31425 → 4.318 |
| ℓ6000 | 1.211e-4 / 1.996e-4 / 1.87e-4 | ✓ | 4.58570109 | 4.58951035 | 0.73573555 | 0.73573753 | −0.14641153 | −0.14641036 | 5.32144 → 5.326 | 4.43929 → 4.444 |

Also checked: A1 = ⌈c_{R,1}/2π⌉₅ ≥ C1 (all six, Arb); the printed A3 = ⌈printed Cζ + printed sup⌉₃ (all twelve). Margins: Cζ 6.9·10⁻⁴ … 3.8·10⁻³; sup F 1.2 … 1.9·10⁻⁶ (the paper's first-order grid overshoot); A3 7.3·10⁻⁴ (ℓ5 A3ᴺ) … 4.7·10⁻³. The (K) minima reproduce stcaudit's sampled margins (ℓ20 1.006·10⁻⁴ vs 1.01·10⁻⁴; ℓ6000 1.211·10⁻⁴ vs 1.22·10⁻⁴).


## (K) for all six lines (certK2.py; logs certK2_edge_all.log, certK2_box_all.log, certK2_L20_*.log)
- **New reduction (PROOF_K.md; reviewer: CORRECT).** H_s(b,τ) = U_s − s f₊ is harmonic in (b,τ) (each term is Re/Im of a holomorphic function of d+b+iτ or an antiholomorphic one of d−b+iτ; arctan(τ/y) = arg(y+iτ)), even in b and → 0 as τ → ∞; by the minimum principle on [−½,½] × [0,T], T → ∞, (K) on the whole half-strip ⟺ (a) U(b,0) ≥ π on [0,½] and (b) H_±(½,τ) ≥ 0 for τ ≥ 0, i.e. only the lines β ∈ {0,1} and t = 0 matter (CG's remark that the constraints bind at b = ±½ is a theorem). Three 1-D certificates per line: (a), (b) on [0,100] by second-order forms v(m) − |v′(m)|r − |v″(box)|r²/2 with exact derivatives (H_bb = −H_ττ = Re ΣΨ_k″), (b) for |t| ≥ 100 by Q = t²U − 4d³u − K5|u|³ > 0 at b = ½ (Lemma 4.10, as stcaudit). Each line ≤ 27 s; accepted only within 1 % of the centre value, so the printed minima are tight lower bounds.
- **Independent 2-D check (no minimum principle).** Bisection of [0,½] × [0,100], s = ±1, with exact gradient at the centre and Hessian enclosures on the box (H(m) − |H_b|r_b − |H_τ|r_τ − |H_bb|(r_b² + r_τ²)/2 − |H_bτ|r_b r_τ > 0), and region B (|t| ≥ 100) by a second-order form in b. All six CERTIFIED: ℓ150 4 278 + 2 126 boxes, ℓ400 9 022 + 7 382, ℓ6000 9 470 + 9 106 (≤ 4 s; the paper needed 14.8 M boxes, certify_K.py did not close in 800 s). Why it closes: at large t, U − f ≈ (2dc_{R,1} − 2c_{R,2})/t² (2·10⁻⁴ for ℓ6000) is a cancellation of terms of size 0.9/t², which uniform centred forms (gradient enclosure on the box) only resolve with boxes of width ~10⁻⁴; the second-order form evaluates the gradient at a point, so the cancellation is exact.
- **ℓ20 reproduction:** CERTIFIED (2 716 + 1 960 boxes; edge mode min max(1,τ²)H ≥ 1.0062·10⁻⁴ at b = ½, τ ≈ 1.0, stcaudit sampled 1.01·10⁻⁴ there). Derivatives checked against mpmath finite differences (vg, hess, Qder; test_derivs.log).
- Inputs as the paper's certifier: c_ι the binary64 coefficients of Appendix D (kernels.json, stcaudit's extraction), d the exact rational Fraction(d).limit_denominator(10⁶) as in the paper's certifier (1, 6/5, 17/20), 80-bit balls.

## Cζ (cz_arb.py; log cz_arb_1e8.log, cz_arb_1e08.json; 289 s, 31 MB, streamed segmented sieve)
- p ≤ 1000 (168 primes): M_p ≤ min(max_x φ_p, Σ_k H_{p,k}) + (powers > 10⁸ in the tail), max φ_p by branch and bound on arcs with φ(c) + |φ′(c)|r + (Σk²H_{p,k})r²/2 (tolerance 10⁻¹³), instead of the paper's 4096-point grid + (π/M)ΣkH (that correction is most of the 0.7–3.8·10⁻³ gap).
- 1000 < p < 10⁸: one Arb ball per prime power p^k ≤ 10⁸ (5 761 287 primes, k ≤ 2), all six kernels in one pass; ψ(10⁸) ∈ [99 998 242.7966 ± 3.8·10⁻⁵].
- n > 10⁸: Lemma 4.3 with G = max(g, g̃(400)) on [10⁸, e⁴⁰⁰] (g certified decreasing there: ν = σ1(A² + B²) − (AA′ + BB′) > 0 by interval bisection, 310–3134 intervals) and G = g̃ = (Σq_iLⁱ + 3)e^{−σ1L}/π beyond (decreasing); ∫_{L0}^{400} h e^{(1−σ1)L}dL by acb.integral (finite enclosure), the rest in closed form (≤ 10⁻⁵²). Tails: ℓ6000 0.48694, ℓ400 0.45090 (they are the largest share of the margin used).
- Agreement: our values equal stcaudit's float recomputation to ≤ 10⁻⁸; the paper's own (4.11) with our tail (float estimate) is 0.8·10⁻⁶ … 1.7·10⁻³ below the printed values, the rest being their step majorant G.
- Remark: g is not monotone on all of [10⁸, ∞): for ℓ150, ν < 0 on L = log x ∈ [652.9, 653.5] (mpmath, checks_misc.log): a(L) changes sign at L ≈ 652.9 while πb ≈ 758 is a near-cancellation of terms ~10⁶, so h dips and d log g/dL = +1.7 there. Harmless for CG (their G is a running maximum and Lemma 4.4 is used only from LB = 2000/5000) and for us (G above); it corrects stcaudit's "g non-increasing from 10⁸ on", which was sampled only to L = 400.
- Also checked (checks_misc.log): kernels.json equals a fresh extraction of the six sol.json blocks of Appendix D from the PDF; g(10⁸) > g̃(400) for all six (so G(N0) = g(N0) in the tail term).


## sup F (elem_arb.py, elem_arb2.py; logs elem_arb_all.log, elem_arb2_all.log; JSON elem_*.json, elem2_*.json)
- F_u, F_ℓ (mode N) and F_u^S, F_ℓ^S (mode S) from (4.2), (4.3), (4.8), (4.16) and Lemma 4.8 (S-mode replaces M(T) by ϑ(T)/π + 1), 96-bit Arb (lgamma, polygamma). On [1, 10⁴] adaptive bisection of a ratio-1.01 partition: on [a,b] F′([a,b]) is enclosed by ball arithmetic; F′ < 0 (> 0) ⇒ bound F(a) (F(b)), else F(m) + |F′|(b − a)/2; refined until within 10⁻⁹ of the best certified point value. T ≥ 10⁴: Lemma 4.7, (4.14) with the paper's K of (4.15) (K = 3.19 … 5.38).
- Every sup is attained near T = 1: mode N by F_u at T = 1; mode S by F_ℓ at T = 1 (ℓ5) or F_u at T = 1.23 … 1.78; our upper bound equals the certified point value to ≤ 10⁻⁹; the printed bounds exceed it by 1.2 … 1.9·10⁻⁶.
- Lemma 4.7 is load-bearing only with a wide margin: its bound at T0 = 10⁴ lies ≥ 7.3·10⁻³ below each sup (smallest: ℓ6000 S −0.153704 vs −0.146412; ℓ400 S 8.2·10⁻³), i.e. ≥ 13 K/T0 (K/T0 ≤ 5.4·10⁻⁴). Cross-check elem_arb2.py: [1, 10⁸] by intervals (second-order forms F(m) + |F′(m)|r + |F″(box)|r²/2 beyond 10⁴), Lemma 4.7 only for T ≥ 10⁸ (K/T ≤ 5.4·10⁻⁸): same twelve suprema.

## Falsification
- No cell where (K) fails: both methods close with positive certified minima on every line; the smallest, ℓ20 at b = ½, τ ≈ 1.0, is 1.006·10⁻⁴ (scaled), and the t → ∞ limit 2dc_{R,1} − 2c_{R,2} ≥ 2.00·10⁻⁴ (ℓ400, ℓ6000) stays positive after Lemma 4.10's loss (Q ≥ 1.87·10⁻⁴).
- No printed constant below its rigorous value: Cζ (6), sup F (12), A1 (6), A3 (12) all hold; the tightest is sup F (≥ 1.2·10⁻⁶ slack, ℓ6000 mode S), the tightest final constant ℓ5 A3ᴺ = 1.369 (slack 7.3·10⁻⁴).

## What is now audited, what is still trusted
- CG Thm 2.1, all six lines of (2.1) and (2.2), T ≥ 1: proof chain (stcaudit) + all numerical certificates (here). In particular |S(T)| ≤ 0.08382 log T + 4.444 and |N(T) − M(T)| ≤ 0.08382 log T + 5.326 (T ≥ 1; (2.3)).
- Still taken from the literature: ψ(x) < 1.03883x (Rosser–Schoenfeld Thm 12), the Stirling remainder of Lemma 4.6 and the explicit K of Lemma 4.7 (structure re-derived by stcaudit; numerically irrelevant after the 10⁸ cross-check), and the correctness of Arb/python-flint.
- Not audited (outside the brief): CG Thm 2.3 and Sections 5–6.

## For the orchestrator
- Records: our certificates may now use CG's 0.08382 log T + 5.326 (N) / + 4.444 (S), T ≥ 1, instead of HSW with the 9.4926 patch (RECORDS.md audit note of 14 Nov); the reruns of 14 Nov already showed the HSW-based results unchanged.
- Candidate tool (not novelty-checked; gate per AGENTS.md): the harmonic reduction makes the Turing-type kernel condition (K) a pair of one-variable inequalities (t = 0 and β ∈ {0,1}) for every kernel built from Φ_X(d ± b + it)^{−j}; this would shrink CG's LP (constraints only on two lines) and any (K)-type certificate (Bellotti–Wong 2026's Turing-type bound uses kernels of this kind with ≤ 4 coefficients, per stcaudit). Second-order forms alone already cut CG's 14.8 M boxes to 1.9·10⁴.
- RH node: N7/N9 explicit S(T), N(T) inputs — moves "CG Thm 2.1" from "audited modulo certificates" to "audited".


## Files (Riemann/Deep/STC2/)
- certK2.py (modes edge, box) + certK2_edge_all.log, certK2_box_all.log, certK2_L20_edge.log, certK2_L20_box.log; certify_K_orig.py (stcaudit's code, copied unchanged); PROOF_K.md (the minimum-principle lemma); test_derivs_H.py, test_derivs_Q.py, test_derivs.log (derivative checks).
- cz_arb.py + cz_arb_1e8.log, cz_arb_1e08.json. elem_arb.py + run_elem.sh, elem_arb_all.log, elem_<line>_<N|S>.json; elem_arb2.py + run_elem2.sh, elem_arb2_all.log, elem2_*.json (cross-check to 10⁸).
- summary.py (table above; A1, A3 checks); checks_misc.py + checks_misc.log; kernels.json (copy of stcaudit's extraction).
- Resources: every run under memcap 290 MB, single-threaded, one at a time; peaks ≤ 31 MB (cz_arb 31 MB, 293 s wall; (K) ≤ 15 MB; sup F ≤ 18 MB). python-flint 0.9.0, mpmath 1.4.1 via uv. Downloads (Zenodo PDF and its pdftotext) deleted after the run; checks_misc.py needs the pdftotext in .scratch/stc2/cg_raw.txt to be rerun.

RH relevance: N7, N9

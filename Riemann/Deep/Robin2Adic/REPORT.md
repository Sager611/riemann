RECORD: Robin's inequality holds unconditionally for every n > 5040 with v_2(n) <= 72 (was v_2(n) <= 20, Axler)

# Robin2Adic — Robin for every n with v_2(n) ≤ 72, from our Robin range, the Lean half-plane and KLN zero density (node N4)

Agent robin2adic, 12 Nov 2026. Method: Gracia Llorente (Zenodo 24044336, unrefereed), every step we use re-derived (DERIVATION.md).

## Claim
**Theorem A.** σ(n) < e^γ n log log n for every n > 5040 with v_2(n) ≤ 72; in particular for every 73-free n. (KLN rows with their Table-1 δ ≈ 0.30; with δ = 1 in every row: v_2(n) ≤ 71.)
**Theorem A′ (slot form).** A counterexample n has n > 10^(10^46.79) and p_j^{a_j+1} > 1.14·10^22 for every j ≤ m (a_j = exponent of the j-th smallest prime factor, p_j = j-th prime): a_1 ≥ 73 (2^73 | n), a_2 ≥ 46, a_3 ≥ 31, a_4 ≥ 26, a_5 ≥ 21, a_6 ≥ 19, and a_j ≥ 2 whenever p_j ≤ 1.07·10^11.
**Theorem B (conditional).** If Chourasiya–Simonič's density (arXiv 2507.15184 v2, Cor. 1.1; unrefereed) holds with B_1 raised for b_4 = 0.5501: v_2(n) ≤ 75, and p_j^{a_j+1} > 1.30·10^23 in every slot (2^76 | n).
Inputs of Theorem A: RobinCH (Robin for 5041 ≤ n ≤ 10^(10^46.79)); Bebel's Lean half-plane (β ≤ 0.87338196661320752); Platt–Trudgian RH(3·10^12); Kadiri–Lumley–Ng Lemma 4.14 with RH(3·10^12) (Ingham4's rows); Hasanalizade–Shen–Wong zero count; BKLNW's θ table; Rosser–Schoenfeld θ < 1.01624x; Bertrand; Choie–Lichiardopol–Moree–Solé (odd n). No CS, no RH above 3·10^12.

## Before
- Published: v_2(n) ≤ 20 and every 21-free n (Axler, Ramanujan J. 2023); v_2(n) ≤ 19 (Hertlein 2018). Unreviewed: ν_2(n) ≤ 24, 25-free (teal-sea note, 10 Oct 2026).
- Conditional: v_2(n) ≤ 73 (Gracia Llorente, 7 Nov 2026: Re s > 7/8 zero-free, CS Table 1 with the unsupported b_4 = 0.54, finite windows to e^106). Harvest6 float estimate: K = 75 with CS; KLN untested.

## Method (DERIVATION.md)
1. GL's exact identity (Prop 2): Robin at n ⇔ B + D + H > E(x), x = p_k, E Nicolas' quantity; slot lemma (Lemma 3): B + D + H ≥ Σ_j D_{a_j}(p_j) > 2^{−v_2(n)−1} for even n. Hence Robin for every even n with 2^{−v_2(n)−1} ≥ E* := sup E(x).
2. Our range replaces GL's finite windows (his Lemma 8(1)–(3), Büthe, Morrill–Platt): n > 10^(10^46.79) forces x ≥ e^{108.57} (BKLNW + Bertrand; Arb).
3. Majorant of E(x) from GL's (5), (11), (16): zeros up to 3·10^12 via 2 + γ − log 4π; ψ − θ ≤ (1 + 1.12·10^{−9} + 1.02 y^{−1/6} log_2 y)√y; zeros above 3·10^12 by Abel summation on an exact dyadic σ-grid from ½ to B (4964 points, step 2^{−16} above 0.8), W(s) = 4∫_{H0}^∞ N(s,t)t^{−3}dt in closed form (incomplete Γ, Arb); G(d, L) ≤ d²/2·sup G_dd with |d| < 2.453·10^{−12}.
4. Every term decreases in x, so E* ≤ majorant(108.57); K = ⌊−log_2 E*⌋ − 1 (2^{−K−1} > E* asserted).
At L = 108.57 (KLN): zeros above 3·10^12: 8.738·10^{−23} (52% from the σ-increments above 0.85); ψ − θ and zeros below 3·10^12: 5.05·10^{−26}; G: 2.80·10^{−26}; log(2π)/(xL): 10^{−49}.

## Certificate
`uv run --no-project --with python-flint python3 Riemann/Deep/Robin2Adic/cert.py`
(Arb, 128 bits; 20–30 s, one thread, peak < 40 MB RSS; imports Ingham4/cert.py unchanged; writes logs/cert.log, logs/cert_result.json.)

| zero density above 3·10^12 | E* ≤ | −log_2 E* | K | 1/E* ≥ |
|---|---|---|---|---|
| **KLN Lemma 4.14, Table-1 rows, RH(3·10^12)** | **8.7462·10^{−23}** | **73.276** | **72** | 1.143·10^22 |
| KLN, δ = 1 in every row | 1.7707·10^{−22} | 72.258 | 71 | 5.65·10^21 |
| CS, b_4 corrected (conditional) | 7.6718·10^{−24} | 76.787 | 75 | 1.303·10^23 |
| CS as published, b_4 = 0.54 | 7.5573·10^{−24} | 76.808 | 75 | 1.323·10^23 |
| GL's (20) as printed at L = 108.57, β ≤ 7/8 / β ≤ B | 8.2671 / 7.8073·10^{−24} | 76.68 / 76.76 | 75 | — |
GL's (20) as printed at L = 106 (reproduction): 1.227034555·10^{−23} (GL: < 1.227035·10^{−23}).

## Checks
- Reproduction first: GL's (20) with his Table 1 at L = 106 gives 1.227034555·10^{−23} (his bound 1.227035·10^{−23}); at L = 108.57 the float model Harvest6/k_est.py (8.265·10^{−24}, 7.805·10^{−24}) is matched.
- check.py (mpmath, independent code): closed-form W(s) vs quadrature at 6 σ for both densities: agree to 4.5·10^{−9} (quadrature error at 20 digits; exact at 40); KLN closed form ≥ the exact-log(kt) quadrature (1.016 where k < 1). Continuous Stieltjes integral (no grid, no running min): cert/check = 1.0017 (CS), 1.0029 (KLN), cert ≥ check. logs/check.log.
- Monotonicity in x proved term by term (DERIVATION.md); sampled at 8 values of L from 108.57 to 2000 in cert.py (decreasing).
- Lemma 3's condition (T+U)log(T+U) < xL/0.998 and both threshold exclusions checked in Arb.
- k_of_L.py (logs/k_of_L.log): with KLN, K = 70 at L = 100, 72 at 108.57, 73 needs L ≥ 111.84 (a Robin range to about 10^(10^48.21)); with CS, K = 76 needs L ≥ 109.54.
- kln_opt.py (logs/kln_opt.log): KLN Lemma 4.14 re-optimised (μ, α, δ, d; 2520 sets, 630 admissible) at σ = 0.85–0.873: best/table ≥ 1.001, so K = 72 is the KLN limit at this range.
- Reviewer subagent (Prop 2, Lemma 3, threshold, majorant; inline): no error; its gaps (T > 1, β ≤ B, two-sided d, W̃_0 = τ/2) are answered in DERIVATION.md "Review".
- Every run < 40 MB and < 2.5 min, one process at a time under .scratch/memcap.sh.

## Novelty
NOVELTY.md: arXiv listings and API (Robin with free / valuation, Sep–Nov 2026), Zenodo ("74-free", "Robin's inequality", "Robin" AND "free", Gracia Llorente's 7 records), the teal-sea note, OpenAlex citations of Axler (4, none on valuations), our TREE/R427. No unconditional v_2 bound above 20 (published) or 24 (unreviewed) found. GL's 73 is conditional.

## Caveat
- Theorem A rests on: our RobinCH range (Chirre–Helfgott's pointwise ψ bound, arXiv 2512.15709, a preprint; Morrill–Platt below 10^10); Bebel's Lean-checked half-plane (formal, unrefereed); KLN Lemma 4.14 evaluated by us with RH(3·10^12) and Table-1 δ ≈ 0.30 (KLN's Thm 1.1 is printed with δ ≥ 1; Ingham3 explains why Lemma 4.14 holds for δ > 0; under δ = 1, K = 71). The other inputs are refereed.
- GL's note is unrefereed; we use only the steps we re-derived (Prop 2, Lemma 3, (5), (11), (12), (16), Lemma 5's summation), not his finite windows or his CS input.
- Theorem B is conditional on 2507.15184 (unrefereed; b_4 corrected via the preprint 2603.25961).
- K is set by zeros just above height 3·10^12 with real part near B; KLN allows about 11× more of them than CS. A larger K needs a sharper refereed density near σ ≈ 0.85–0.873 for T ≈ 3·10^12–10^14, RH to a greater height, or a longer Robin range.
- The gain over GL comes from inputs only (range, half-plane, density); a reviewer may call it a routine combination. Discovery gate: a new theorem (unconditional v_2 ≤ 72) whose inputs count once (RobinCH).

RH relevance: N4

## Files
- REPORT.md (this), DERIVATION.md (GL re-derived with his numbers; our inputs, threshold, majorant; review), NOVELTY.md
- cert.py (certificate), check.py (independent float check), k_of_L.py (K against the range), kln_opt.py (KLN parameter search)
- logs/cert.log, logs/cert_result.json, logs/check.log, logs/k_of_L.log, logs/kln_opt.log


**AUDIT PASSED** — Li200k's theorem (λ_n > 0 for n ≤ 78.714 T² with 0 < β < 1, 79.739 T² with MTY, 105.774 T² with 7/8, 106.2557 T² with θ₁₄₁) and DirLi's ranges (Re λ_n(χ) > 0 for n ≤ 26·t0(q)², 35·t0(q)² with 7/8, every primitive χ mod q, 3 ≤ q ≤ 400,000) hold as written; the one false input, HSW's misprinted 9.3675 still in the shipped Li200k/verify.py and REPORT, changes no certified N (9.4926 reruns: identical).

# LiRangeAudit — referee report (lirangeaudit, 19 Nov 2026)

Scope: Riemann/Deep/Li200k (REPORT l. 23–72, verify.py) and Riemann/Deep/DirLi (REPORT, PROOF.md, certify.py, results_*.txt). Chain with source line numbers: [FULLPROOF.md](FULLPROOF.md) (A = Li200k, B = DirLi, I = inputs, F = tests, R = reviewer).

## Verdicts

| item | verdict | where |
|---|---|---|
| Li200k Theorem, L l. 23–35 | holds; N(a) is defined only for the tabulated a (schema) | A0–A6 |
| λ_n = Σ_ρ[1 − (1−1/ρ)^n], rearrangement into M_n − Q_n | holds (pair sums absolutely summable) | A0 |
| per zero: β = ½ gives g(nθ) ≥ 0; \|β − ½\| ≤ a gives quartet ≥ 2g(nθ) − 2q_n (ρ and ρ* each charged once) | holds; re-derived; 18,756 + 12,240 certified corner cases | A1, B1, F1 T1, F2 D1 |
| hypotheses a/t², a²/t³ nonincreasing | holds; used only to make q_n decreasing (loss) | A1, A2 |
| zero count above T | HSW Cor. 1.2 (t ≥ e), +7/8 conservative; constant must be ≥ 9.4926: shipped verify.py has 9.3675 (false), 9.4926 reruns give identical N | I1, F4 |
| Stieltjes loss, boundary term 2q_n(T)E(T), directions | holds; w-constant 0.13 valid only for T ≥ 1.84·10⁴ (holds at T0) | A2, R |
| gain (n < 10⁴ and n ≥ 10⁴), summation by parts, piece count | holds (L l. 52 prints stale Trudgian constants 0.111/0.275; the code uses HSW's) | A3, A4 |
| n ≤ 88 via the zero at 14.13 | holds | A5 |
| chain: every n ≤ N covered, N a checked point | holds; N sits 1.07·10⁻⁴ below the exact frontier of the proven inequality (Riemann-sum slack) | A6, F3 |
| DirLi definition, Re λ_n for complex χ | holds (Re of the Taylor/Cauchy coefficient = Re λ_n; tested on a complex χ mod 5) | B0, F2 D2 |
| DirLi Lemma 1, ½M_n − ½Q_n, Lemma 2, gains A/B/C, loss, q-blocks | holds | B1–B6, R |
| Platt t0(q) | by parity of the modulus q, both signs of γ, as used | I2 |
| BMOR Thm 1.1 incl. "ℓ ≤ 1.567 ⇒ N(T, χ) = 0" | verbatim | I3 |
| DirLi headline table, RECORD line, uniform 26.12/35.73 | holds (results parsed: 749 moduli, 604 blocks covering all 299,250 moduli) | B7 |

## Defects (none changes a certified number)

1. Li200k/verify.py (l. 14 docstring, l. 31 comment, l. 33 C245 = r(93675, 10000) + 7/8) and Li200k/REPORT l. 40 still use HSW's printed C3 = 9.3675, a misprint for 9.4925 (Bellotti–Wong footnote, CG Table 1; plus the 1.5·10⁻⁵ log ζ(1+η) slip, STCAudit; BW's TeX not re-read here). The shipped file is therefore not by itself a valid certificate; Riemann/Deep/HP141/Li200k/verify.py (C3 = 9.4926) is, and certifies the same N for A, B, C (E enters the frontier only at relative order 10⁻¹⁰). Fix: one-line patch r(93675, 10000) → r(94926, 10000), or point links at the HP141 copy.
2. Li200k/REPORT l. 52 prints Trudgian's (0.111 + 0.275/log t_n); the code uses HSW's 0.1038/0.2573 (correct).
3. verify.py's w uses tE′(t) ≤ 0.13, true only for t ≥ exp(2573/262) ≈ 1.84·10⁴: fine at T0 (0.1128), invalid if the script is reused at small T (fals_zeta.py uses 0.1038 + 0.2573/log T).
4. The Theorem's N(a) is only defined for the tabulated a (the statement is a schema over the certified rows).
5. verify.py certify() returns N = 89 without error if the first block fails (89 would then be uncertified); not the case here.
6. DirLi RECORD line and TREE N3 say "every primitive χ mod q ≤ 400,000": the certificate covers 3 ≤ q (q = 1 is ζ, covered far beyond by Li200k; q = 2 has none). True, imprecise; COROLLARIES l. 20 already says 3 ≤ q.
7. DirLi REPORT l. 83–84 ("Li200k … used Trudgian 2014 … should be re-run with HSW") is outdated (done in R333). It prints the q ≤ 1000, a = 3/8 constant as 55.95 (l. 23, rounded) and 55.94 (l. 27; value 55.94929).

## Falsification (no counterexample found)

- Per-zero/pair inequality at corner zeros (β at the strip edge, γ = T + 10⁻⁶ and 1.001T for T from 14.15 to T0 = 3.0·10¹², n up to 80γ²/a, z^n computed directly): 18,756 (ζ) + 12,240 (pairs, both signs of γ) certified Arb cases; min slack/q_n 5.2·10⁻⁶.
- λ_n for n ≤ 200 by the Cauchy integral (|z| = 0.85) equals the Taylor values; ζ: λ_1 = 0.0230957089661, λ_200 = 306.6558. Dirichlet: Re λ_n(χ) for χ₋₃, a complex χ mod 5 (χ(2) = i) and the Legendre χ mod 5, all n ≤ 200, all > 0.
- The theorems at small heights (ζ: T = 100, 1000; DirLi: t0 = 50, 100 for q = 3, 5): every certified lower bound lies below the exact λ_n (n ≤ 200; max ratio 0.84). ζ against Odlyzko's 10⁵ zeros on 112 grid points n ≤ 3N_T: gain bound ≤ true M_n, loss bound ≥ Σ_{T<γ≤74920} q_n(γ), and the admissible configuration "zeros above T paired into quartets at β = ½ ± a" has λ_n ≥ our bound — 0 failures; that configuration first gives λ_n < 0 at 1.21–1.53 N_T, never below N_T.
- Independent frontier (mpmath quad of the loss with exact E′, no Riemann sums; the gain formula of V l. 84–103): every certified N lies below the frontier n* of the proven inequality, by 1.06–1.07·10⁻⁴ (the Riemann-sum slack). For a = ½, n* = 7.0858152946·10²⁶ reproduces Fabbian's continuation 708,581,529,455,821,832,027,338,082 (78.72 T²) to 10⁻¹¹: Li200k's 78.714 T² is the same inequality, rigorously discretised.


## Dependents (rg -n -i 'li200k|dirli|106\.25|105\.77' discoveries/ TREE.md COROLLARIES.md RECORDS.md DISCOVERIES.md)

- TREE.md l. 7 (N3 record line: 105.77 T², 106.2557 T², Dirichlet 26.12·t0(q)²), l. 369 (R333), l. 419 (R383), l. 436 (R400); COROLLARIES.md l. 8 (105.77 T²), l. 20 (DirLi 26.12·t0(q)², 41.30·10¹⁶/q², 35.73·t0(q)²); RECORDS.md l. 44 (9.3675 → 9.4926 audit note). discoveries/ and DISCOVERIES.md: no entry. Every one stands.
- Suggested edits (orchestrator; this audit edits neither Li200k nor DirLi): patch Li200k/verify.py C3 to 9.4926 (or link COROLLARIES l. 8 to HP141/Li200k/verify.py) and refresh Li200k/REPORT l. 40, 52; TREE l. 7 "every primitive χ mod q ≤ 400,000" → "mod q, 3 ≤ q ≤ 400,000" (as COROLLARIES l. 20).

## What survives

Everything: all printed N of Li200k (A, B, C; F = 106.245188 T² and G = 106.255746 T² with 9.4926) and of DirLi (749 per-modulus values, 604 block values, both variants, the uniform 26.12/35.73 and 41.30/55.94 constants). Variants C, G and DirLi's a = 3/8 rest on the Lean-checked half-planes (7/8, θ₁₄₁), not re-audited here; A, B and DirLi's a = ½ use refereed inputs only.

## Side observation (not a discovery)

At T = 100 and 1000 an admissible configuration (true ordinates above T, quartets at the strip edge) already makes λ_n < 0 at 1.2–1.5 N_T, so ranges of this type are within a small factor of what "RH to T + strip a" can give (consistent with TREE R400 and Fabbian's Prop. 8.1); nothing about zeros above T.


## Own reruns (this audit; one process at a time, memcap 290 MB)

| run | result | wall, peak |
|---|---|---|
| HP141/Li200k/verify.py (= Li200k with C3 = 9.4926, variant G) | E(T0) = 14.213738; part 0, 1, 3 pass; N = 708,505,587,700,793,981,366,435,840 (78.713642 T²), 717,734,280,405,823,147,676,270,592 (79.738933), 952,080,909,647,480,761,027,133,440 (105.774404), 956,413,498,941,653,124,095,410,176 (106.255746): identical to the claims | 84 s, 20 MB |
| DirLi/certify.py small | q = 3, 4, 5, 7, 8, both a: all ten N identical to the REPORT table | 5 s |
| DirLi/certify.py blocks 1.02 301948 307985 | the minimising block odd [301949, 307985]: N = 2,865,676 / 3,919,034 = 26.1273 / 35.7311 t0(q1)², identical to results_blocks.txt l. 289 | 2 s |
| fals_zeta.py, fals_dirli.py, indep_check.py, indep_g.py, parse_dirli.py | §Falsification above; all pass | 534 s / 218 MB, 87 s / 29 MB, 558 s / 32 MB, 93 s / 19 MB, < 1 s |

## Review

Reviewer subagent (verdict only, the per-zero inequality, both bookkeepings, the Stieltjes loss and the gain pasted inline): CORRECT on all four; its two "GAP as stated" remarks concerned hypotheses my paste omitted (GRH for |γ| ≤ t0 including γ = 0; the zero count valid on [t_n, ∞)), both present in the sources (Platt l. 306; HSW for t ≥ e). Its threshold for the 0.13 constant, T ≥ exp(2573/262) ≈ 1.841·10⁴, matches defect 3.

## Files (Riemann/Deep/LiRangeAudit/)

REPORT.md; FULLPROOF.md (the chain with source line numbers); fals_zeta.py, fals_dirli.py, indep_check.py, indep_g.py, parse_dirli.py, rerun.sh; logs/ (fals_zeta, fals_dirli, indep_check, indep_g, parse_dirli, rerun_verify141, rerun_dirli_small, rerun_dirli_block). Downloads (arXiv TeX of HSW, Platt, BMOR; Odlyzko's zeros1) deleted.

RH relevance: N3

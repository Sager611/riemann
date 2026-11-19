AUDIT PASSED — Wang's local pair theorem (arXiv 2609.07918 v2, thm:local-pair l.229–240), its weight removal (eq:unweighted) and prop:optimal hold as written, unconditionally, every step re-derived and falsification-tested; c(θ) is a proved lower bound, and Santibáñez-Leal's Thm 5.1 uses them only inside the proved range, so Theorem E1 and the e1arb rows stand.

# WangAudit — line-by-line referee of Wang's local pair theorem (load-bearing input of E1 and e1arb, N8)
wangaudit (Orca), 19 Nov 2026, 18:55–19:25Z. Sources: Wang v2 TeX (zeros_short_intervals.tex, 1441 l.) and v1 TeX; Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh
(BGST) arXiv 2306.04799 v1 TeX (= Acta Arith. 214; its equation numbers (2.7), (2.14) are the ones Wang cites); Santibáñez-Leal (SL) Zenodo 23134787 v0.03 PDF.
The chain with TeX line numbers: [FULLPROOF.md](FULLPROOF.md) (§ numbers below).

## Short answer
- thm:local-pair: for fixed 0 < λ < θ < 1 and fixed real even g ∈ C_c^∞, supp g ⊂ [−λ, λ], Σ_{γ,γ′∈I} ĝ(i(ρ−ρ′)L/2π)w(ρ−ρ′) = (HL/2π)(g(0) + ∫|α|g) + O_g(H + T^λL²).
  TRUE as written; no RH, no zero-free region beyond |δ| < 1/2, no zero-density or simplicity input (§2–§7, §12). Reviewer: CORRECT (§11).
- c(θ) = 2 − C_θ, C_λ = λ/2 + cot(λ/√2)/√2: PROVED lower bound. The cosine density f_λ is positive on J = [−λ/2, λ/2] (λ/√2 < π/2) and is reached inside the
  admissible class by η_ε = χ_ε√f_λ/√b_ε ∈ C_c^∞((−λ/2, λ/2)), real, even, ∫η_ε² = 1, with C(f_ε) → C_λ; then λ ↑ θ (C_λ decreasing, continuous). Minimality
  (Poincaré, coercivity 1 − 2λ²/π²) is also correct but not needed for the bound (§9).
- SL Thm 5.1 uses exactly: q_T = S_K(I)/N(I) → C(f) for each fixed admissible f (Wang §3 lemma, which SL cites as "[14, Theorem 2.2]"), then f → f_λ, then
  λ → θ; R(2 − c(θ), 1 − k) = h(θ; k) ✓. It does not use Lamzouri's Prop 2.1, (eq:both-counts) or thm:distinct-all (§10). v1 (cited by SL) and v2 have
  identical proofs (§10).
- Effect on E1's h(θ; κ) and e1arb's rows (0.0237843 at θ = 0.51, uniform 0.0092120): none.

## What was checked (MECE)
| block | source lines | method | result |
|---|---|---|---|
| lem_norm 2πF_I = ∫\|A_I\|² | Wang 268–302; BGST (2.7) | residues by hand; contour R → R + iδ crosses no pole (Im ±1, δ+δ′±1); reflection ρ′ ↦ 1 − ρ̄′; numeric with off-line zeros | ✓ (1e−31; fails without the reflection partners, as it must) |
| lem:localization O(xL³) | 305–437 | every bound re-derived (\|1+(s+iδ)²\| ≥ s² + 3/4, x^δ ≤ x^{1/2}, unit-interval counts, convexity) ; numeric T = 600 | ✓ (ratio ≤ 0.73) |
| explicit formula (eq:explicit), (eq:E) | 246–266; BGST Lemma 1 (203–227) | Landau's formula at s = 3/2 + it, −1/2 + it re-derived; jump cancellation at prime powers; trivial-zero series combined (uniform as x ↓ 1); exact identity tested with 1250 zeros | ✓ (1e−5; x\|E_x\| ≈ 1.9) |
| lem:coeff | 439–466 | partial summation by hand; sums to x = 10⁵ | ✓ |
| prop:pointwise | 468–513 | MV Cor. 3 on a length-H interval; B·D̄ by termwise IBP; Cauchy–Schwarz for E | ✓ |
| proof of thm:local-pair | 516–640 | M_0 (MVT on g), M_1, five error integrals recomputed | ✓ |
| weight removal (eq:remove-weight), (eq:unweighted) | 643–783 | 1 + π²z²/L² = 1/w(ρ−ρ′) at z = i(ρ−ρ′)L/2π; M(Q) = C(f), M(Q″) = 2∫f² − ∫f′²; conjugation invariance of the scaled zeros | ✓ |
| prop:optimal, thm:main | 790–942 | G″ = 0, C(f_λ) = C_λ, Poincaré, η_ε family; quadrature, 1200 random perturbations | ✓ |
| SL Thm 5.1 use | SL l.844–870 | exact hypotheses vs Wang's (eq:eta); order of limits; R(A, b) = h | ✓ |
| versions / status | arXiv abs, v1 vs v2 token diff | v1 7 Sep, v2 6 Oct 2026, no v3, no comments; §2–§4 proofs identical | ✓ |

## Where each hypothesis enters (thm:local-pair)
| hypothesis | used at | why |
|---|---|---|
| λ < θ | prop:pointwise (‖D_x‖² ≪ HL needs x log²x ≪ H); thm:local-pair error T^λL² = o(HL); (eq:unweighted) | the crude off-line bound x^{δ+δ′} ≤ x costs O(xL³), affordable only for x ≤ T^λ, λ < θ (BGST needed a zero-free region for x up to T; Wang does not, §12) |
| θ < 1 | (eq:A-tail) (H ≤ T, t ≤ 2T), (eq:E) (x < T), (eq:intro-rvm) (H ≤ T); prop:optimal needs λ ≤ 1 (f_λ > 0, coercivity) | |
| g fixed, real, even, C_c^∞, supp ⊂ [−λ, λ] | O_g = O(‖g‖_∞ + ‖g′‖_∞); evenness for F_I(1/x) = F_I(x); support for x = T^α ≤ T^λ | C¹ suffices for thm:local-pair itself |
| \|δ\| < 1/2 (PNT) | (eq:denominator), x^δ ≤ x^{1/2}, \|Im a\| < 1 in (2.7) | the only information on off-line zeros used |
| uniformity in T | all constants absolute or depend on λ, θ, g; "T large" only in ‖D_x‖ ≪ √(HL) | liminf statements need nothing more |

## Numerical sanity (actual zeros, rough by design; §13)
- With 1250 zeros (mpmath.zetazero, γ ≤ 1704), windows (T, θ) = (1000, 0.9), (1200, 0.75), (600, 0.95), f = f_λ, λ = 0.95θ: S_K/N = 1.408, 1.593, 1.356 vs
  C_λ = 1.451, 1.639, 1.405 when L = log(T/2π) (2–3 % low); with Wang's L = log T 16–19 % low, i.e. log 2π/log T: the O(H) term. Only rough agreement is
  possible at T ~ 10³, as expected; normalisation (HL/2π)(g(0) + ∫|α|g) is right. F(α) tracks α + L T^{−2α} for α < θ (e.g. 0.99 vs 1.00 at α = 1, T = 1000).

## Minor remarks (none affects any conclusion)
1. lem_norm applies (2.7) on the line R + iδ without stating the contour shift (no pole crossed; §2).
2. BGST Lemma 1: the two trivial-zero series diverge separately as x ↓ 1; only their combination is O(x^{−5/2}/(1+t²)) (reviewer; §11). Statement fine.
3. E_x's O(1)/x is ≈ −log 2π/x (B_x uses log(t+2) instead of log(t/2π)); harmless, but explains the finite-height discrepancy above.
4. SL cites "[14, Theorem 2.2]" (= thm:local-pair) for q_T → C(f), which is Wang's §3 lemma, a direct consequence; SL's class is C_c^∞ (the PDF layout puts
   "∞" on the next line; SantiAudit's "SL writes C_c" note is an extraction artifact).
5. v1's acknowledgment: AI-assisted drafting of the proof, author-verified; v2 drops the note. Proofs unchanged.

## Not checked (not load-bearing for E1/e1arb)
- thm:distinct-all and §5 (l.980–1300): not used by SL Thm 5.1 or E1. Lamzouri's Prop 2.1 proof (used by Wang's thm:main, not by SL/E1): citation checked
  in R481; 20000 random multisets here, no violation. Classical inputs taken as known: Montgomery–Vaughan Cor. 3, Davenport ch. 15 lemma, PNT, Landau's
  formula, Hardy–Littlewood–Pólya Thm 257 (Wirtinger). SL's parity product: audited in SantiAudit.

## Files
FULLPROOF.md (chain, §1–13); zeros.py, zeros.log, zeros1250.txt (first 1250 ordinates; scripts take its path); lemnorm_test.py/.log, int27_recheck.py/.log, bgst_lemma1.py/.log, localization_num.py/.log,
localpair_num.py/.log, finite_checks.py/.log (+ finite_checks_run1_N4e6.log). All runs single-threaded under memcap 290 MB; footprints 29–153 MB.

RH relevance: N8

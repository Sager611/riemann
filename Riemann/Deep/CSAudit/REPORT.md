AUDIT PASSED for the corrected statement — every step of Thm 2.1, of the fourth moment (Thm 2.2, Cor 2.3, Cor 2.4, §5 parts I–VI) and of Cor. 1.1 is checked and all Mathematica outputs (𝔐₁ at 9 heights, 𝒞₁ = 20.7225, M₂(10¹⁰,2·10¹⁰) ≤ 9.5079·10¹⁵) are reproduced; Cor. 1.1 holds with b₄ = 0.5501 and B₁ raised by ≤ 1.4% (rows 1–13, §D), not as printed (b₄ = 0.54 unsupported), and Lemma 5.2's second bound is false for σ₀ > 2 but unused there.

# CSAudit — referee report on Chourasiya–Simonič, "An explicit form of Ingham's zero density estimate", arXiv 2507.15184 v2 (TeX InghamPostArXiv.tex, 1920 lines)

csaudit, 12 Nov 2026, 19:31–20:22 UTC, follow-up 21:07–21:20 UTC. Line numbers = TeX lines. Proof notes: FULLPROOF.md (§1–8). Earlier audit: Riemann/Deep/Ingham4/AUDIT-2507.15184.md (R410).

## Short answer
- **Cor. 1.1 as printed (b₄ = 0.54): unproven, not refuted.** b₄ = 0.54 (footnote l. 659, "private communication") is not what Ramaré–Zuniga Alterman published (arXiv 2603.25961 v3, Thm 1.2: 0.5501 for 10¹² ≤ X < 10¹⁶, ε ≤ 1/50); printed B₁ is 0.08% (row 1) to 1.4% (row 13) too small in rows 1–13. B₂, B₃ do not involve b₄.
- **Corrected (b₄ = 0.5501): every step checked, no error.** Cor. 1.1 holds with the B₁ of §D. This includes Cor 2.4 (fourth moment): §5 parts I–VI re-derived, and 𝔐₁(10¹⁰) = 34.842, 𝒞₁(3000) = 20.7225 and the whole 𝔐₂ recursion (M₂(10¹⁰,2·10¹⁰) ≤ 9.5079·10¹⁵) recomputed from the paper's formulas (ours ≤ printed in every entry, ratios 0.99975–1).
- **Defects found (none changes a number)**: Lemma 5.2's second bound is false for σ₀ > 2 (counterexample §B.2; used only at σ₀ < 2); (eq:J22) states a per-k bound but proves the grouped one (enough); Thm 2.2's general parameter range allows Re(σ₄+c₂) > 1, outside Lemma 5.9 (all printed rows have ≤ 0.823); typos in (eq:F1F2main) and VI.
- **Still taken as published** (§C): external lemmas (CHT Σd², BBR Δ, Dona–Zuniga Alterman Γ bounds, Hiary–Patel–Yang 0.618, Preissmann/MV, Ramaré b₁–b₃, Rademacher, Gabriel); the M₂(3000) quadrature is numerical (agrees to 10⁻¹⁰), and ζ(1+ε) ≤ e^{γε}/ε (needed only for b₄ = 0.5501 exactly; b₄″ avoids it).

## A. Checked, no error (MECE by section)
| part | TeX lines | how checked | result |
|---|---|---|---|
| Cor. 1.1 from Thm 2.1 | 249–275 | re-derived (h = 1, T₀ = H₀ = H = H_RH; d ≤ 0.89 < log(3·10¹²)/32; σ₁ := σ monotonicity for row 1; near-line N(T) bound) | ✓ |
| Gabriel, Thm 3.2, Cor 3.3 | 433–577 | every step re-derived (bounds on g, integration by parts, dropping the negative u-range, F(x) ≤ (β−α)³/4, (1−1/H)² bound) | ✓ (Lemma 3.1 = Gabriel, cited to Karatsuba–Voronin; standard) |
| Littlewood setup | 612–625 | zeros of ζ ⊂ zeros of h_X; arg principal on Re = 3/2 | ✓ |
| Lemma 4.1 | 627–642 | = KLN Lemma 4.1 (KLN TeX l. 636–673), MV for real u_n, Ramaré's Σμ² bounds; ranges X ≥ 10⁹ | ✓ |
| Lemma 4.2 | 644–670 | MV split X<n<5X / n ≥ 5X (KLN (3.5), b₃ = 0.605), σ ∈ (1, 21/20] ⇔ τ ≤ 1.1; algebra to C₃, C₄ | ✓ except b₄ |
| Lemma 4.3 | 672–700 | b₆(X,1/2) = 3X^{−1/2}(log X+2+γ+7/(36X)) decreasing; chord bound for −log(1−y); MV on [H,T]; KLN (3.7) (d² ≤ d₄, Σd₄ ≤ t(log t+1)³); b₁₁(X,2) ≤ 1+3.31/log X (6/L+6/L² = 0.3035 at 10⁹) | ✓ ("continuation principle" not needed) |
| Lemma 4.4 | 702–720 | KLN Prop 4.10, Cor 4.11, (4.53) redone at η = 1/2 directly (Rademacher needs only η ≤ 1/2); true constant b₉ ≤ 38.93 + small, CS use 220 | ✓ (conservative) |
| Lemma 4.5 | 722–839 | Hölder/convexity (eq:MinHol); K₁ (Cor 2.4 at V ≥ A₀), K₂ (finer route: 2^{1/3}c[Γ(3/2)+u₀^{3/2}], c·u₀^{3/2} ≤ (4κ/3)(A₀/T₀)³T), K₃; all exponent bounds (3/4)w₁/w ∈ [e₀, e₀+3ε] etc.; L₀, L₁, L₂ | ✓ |
| Thm 2.1 proof | 841–859 | log(1+u²) ≤ u^a (a ∈ [4/3,2]); σ′ → σ exponent shifts (≤ 3d, 3d, 4d over log T); 𝔄₁…𝔄₅ | ✓ |
| Table 1 | 109–131 | independent recomputation (table1_repro.py, mpmath 30 digits) of all 16 rows at the printed (d, κ, δ, A₀) | ✓ ours/table 0.9991–1.0000 for B₁, B₂ and B₃ to the printed digits |
| §5.1–5.2 lemmas | 871–1243 | Lemma 5.1 re-derived; Lemmas 5.2–5.8 falsification grids x ≤ 2·10⁵ (divisor_falsify.py) | ✓ except Lemma 5.2b (below); 5.5b tight as σ → 1⁺ but true (slack ≥ 0 analytically) |
| Lemmas 5.9, 5.10 | 1252–1295 | numerical grids (chi_falsify.py) + Lemma 5.10 re-derived | ✓ (max ratios 0.9989 for |(z−1)χ| ≤ 6.31, 0.9934 for |χ| ≤ 0.4|y|^{1/2−x}+16) |
| Thm 2.2 parts I, II, IV (𝒥₃, 𝒥₄), V; assembly of 𝔉₁ | 1297–1345, 1568–1602 | re-derived (MV + Lemmas 5.2a, 5.5a, 5.6–5.8, χ bounds, contour shift, (A+B+C)² grouping) | ✓ |
| Cor 2.3 inputs | 371 | M₂(3000) = 703127.63, M₂(3000,6000) = 1249636.00 (Arb ζ, Gauss–Legendre, two node counts agree to 3·10⁻⁴; m2_check.py) | ✓ inside the stated brackets |
| Cor 2.4 | 346–360 | from (eq:Cor2M2) at T₀″ = 10¹⁰ ≤ A₀ (within range); a₂(T) decreasing (d a₂/du > 0, u = 1/log T) | ✓ |
| (eq:zeta2) decomposition, parts III (𝔊₁…𝔊₆, 𝔤₁…𝔤₈, 𝒥₂,₂), IV (𝒥₅, 𝒥₆), VI (𝒥₀, 𝒥₀₁, 𝒥₀₂) | 886–948, 1346–1650 | every bound re-derived (contour shifts and holomorphy, χ/Γ bounds with the right Re ranges, CBS, Lemma 5.1, Lemmas 5.2–5.8 at in-range exponents, T-powers) | ✓ (FULLPROOF §8) |
| Mathematica outputs | 186–194, 371, 396–427 | f1_repro.py: 𝔉₁ at the printed parameters for all 9 heights (10¹⁰: 34.841 ≤ 34.842), conditions (2)–(6) hold; m2rec_repro.py: 𝒞₁(3000) = 20.72246 ≤ 20.7225, recursion M₂(10ⁿ,2·10ⁿ), M₂(10ⁿ), 𝔐₂(10ⁿ) for all n (10¹⁰: 9.50787·10¹⁵ ≤ 9.5079·10¹⁵; 𝔐₂ 1.38962·10¹⁶ ≤ 1.3897·10¹⁶) | ✓ ours ≤ printed everywhere |

## B. Errors and unsupported inputs
1. **b₄ = 0.54 (Lemma 4.2, l. 653–659) unsupported** (known since R410). RZA Thm 1.2 (sagetex output re-read): S_ε(X) ≤ 0.5107 (422 ≤ X < 10⁸, ε ≤ 1/100), 0.689 (10⁸ ≤ X < 10¹², ε ≤ 1/100), 0.5936 (10¹²–10¹⁶, ε ≤ 1/25), **0.5501** (10¹²–10¹⁶, ε ≤ 1/50), 0.5103 (10¹⁶–10³³, ε ≤ 1/50), 0.4905 (X ≥ 10³³, ε ≤ 1/50). Thm 2.1 as stated (X ≥ 3·10¹⁰, τ ≤ 1.1) is not supported for X < 10¹² (published 0.689 > 0.54 there, ε ≤ 1/100 only). For Cor. 1.1: X = T ≥ 3·10¹², ε = 2δ/log X ≤ 0.0150 (rows 1–13) ⇒ b₄ = 0.5501 suffices **if** ζ(1+ε) ≤ e^{γε}/ε (true on (0, 0.1] by the Laurent series: εζ(1+ε)e^{−γε} = 1 − 0.0938ε² + O(ε³); not proved here with interval arithmetic); the rigorous elementary substitute ζ(1+ε) ≤ 1 + 1/ε gives b₄″ = 0.5501(1+(1−γ)ε), which is what Ingham4/cert.py uses.
2. **Lemma 5.2, second bound (new)**: partial summation with |Δ(u)| ≤ √u gives Σ_{n>x}d(n)n^{−σ} ≤ x^{1−σ}log x·[(σ−1) + (1+2γ(σ−1))/log x + (σ−1)²(4σ−1)/((2σ−1)√x log x)]/(σ−1)²; the printed D₂ drops the factor (σ₀−1)² in its last term. Counterexample: x₀ = 5, x = 5.99, σ = σ₀ = 4: LHS = Σ_{n≥6}d(n)n^{−4} = 0.0068135 > RHS = 0.0058938. Harmless: D₂ is used only at σ₀ = 1−σ+|c₁| ∈ (1, 3/2] (𝔤₅, l. 1510), ∈ [1/2+|c₁|, 1/2+|c₁|+c₂) ⊂ (1, 2) (𝔤₇, l. 1533) and σ₀ = 1/2+|c₁| (𝒥₀₁, l. 1622), where (σ₀−1)² < 1. Reviewer subagent: (2), (3) CORRECT; it corrected my first counterexample point (at x = x₀ = 6 the bound holds).
3. Cosmetic: (eq:MinHol) is convexity of x^{4/3} + Hölder, not Minkowski; Cor 3.3's m = 1 exceeds the pole order if M_X(1) = 0 (harmless); KLN Prop 4.10 says σ ∈ (0, 1+η] but Jensen's step needs σ ≥ 1/2 (CS use σ ≥ 1/2); KLN Lemma 4.1's proof swaps the two μ² bounds (algebra right).
4. (eq:J22) (l. 1556) states 2Σ_{k=3}^{6}|∫(J₁+J₂)J_k| ≤ 2𝒥₂,₂T log³T but proves the grouped |∫(J₁+J₂)(J₃+J₄)| + |∫(J₁+J₂)(J₅+J₆)|; the per-k form also follows (the proof splits g₁…g₄ by the triangle inequality) and only the grouped form enters |R|. Harmless.
5. Thm 2.2's hypotheses allow σ₄+c₂ (and σ+c₂ on the horizontal pieces) up to 3/2, but Lemma 5.9 is stated for Re z ∈ [−1,1]; every printed parameter row has σ₄+c₂ ≤ 0.823, so all numbers stand. Typos: (eq:F1F2main) lower limits 1+iT, iT should be σ₂+iT, σ₃+iT; VI writes χ(1/2+c₁+…) for χ(1/2+c₂+…).

## C. Taken as published (not re-proved), precisely
- External lemmas: Preissmann's constant m₀ in MV (Lemma 5.1, KLN Lemma 3.6); Ramaré's b₁, b₂, b₃; Cully-Hugill–Trudgian's Σd² formula with |D₀| ≤ 9.73x^{3/4}log x (tested to 2·10⁵, max ratio 0.19); Berkane–Bordellès–Ramaré |Δ(x)| ≤ √x; Dona–Zuniga Alterman's Γ and cos bounds (tested on grids); Hiary–Patel–Yang |ζ(1/2+it)| ≤ 0.618t^{1/6}log t; Rademacher's convexity bound; Gabriel's theorem (Lemma 3.1); RZA Thm 1.2 values; Bellotti–Wong's N(T) bound.
- Numerics: M₂(3000), M₂(3000,6000) by Gauss–Legendre on Arb values (not an interval enclosure; two node counts agree to 3·10⁻⁴, brackets have margins 28 and 364).
- ζ(1+ε) ≤ e^{γε}/ε on (0, 0.02]: series argument only; needed only to use b₄ = 0.5501 exactly — the elementary ζ(1+ε) ≤ 1+1/ε gives b₄″ (cert.py's choice, last column of §D).
- Sources of the two inputs named in the follow-up: 𝔐₁(10¹⁰) = 34.842 and M₂(10¹⁰,2·10¹⁰) ≤ 9.5079·10¹⁵ have no external source — they are the authors' Table 2, i.e. Thm 2.2 at Table 3's parameters and the recursion of l. 418–427; both now recomputed (§A). No published explicit fourth-moment bound of the shape aT log⁴T exists to cross-check; max|ζ|²·M₁ (Hiary–Patel–Yang) gives only ≈ 1.3·10¹⁷ for M₂(10¹⁰,2·10¹⁰).

## D. Corrected Cor. 1.1 rows (table1_repro.py; B₂, B₃ unchanged)
| row [σ₁, σ₂] | B₁ table | B₁ ours (0.54) | ε = 2δ/log X₀ | B₁ (b₄ = 0.5501) | B₁ (b₄″ = 0.5501(1+(1−γ)ε)) | B₂ | B₃ |
|---|---|---|---|---|---|---|---|
| 0.50000, 0.53125 | 5.069 | 5.0689 | 0.00084 | 5.0729 | 5.0730 | 9.461 | 167.8 |
| 0.53125, 0.56250 | 6.046 | 6.0454 | 0.00174 | 6.0551 | 6.0555 | 9.168 | 162.6 |
| 0.56250, 0.59375 | 7.077 | 7.0760 | 0.00265 | 7.0935 | 7.0946 | 8.875 | 157.4 |
| 0.59375, 0.62500 | 8.185 | 8.1848 | 0.00362 | 8.2124 | 8.2147 | 8.582 | 152.2 |
| 0.62500, 0.65625 | 9.387 | 9.3860 | 0.00459 | 9.4265 | 9.4307 | 8.290 | 147.0 |
| 0.65625, 0.68750 | 10.70 | 10.6905 | 0.00564 | 10.7472 | 10.7545 | 7.997 | 141.9 |
| 0.68750, 0.71875 | 12.11 | 12.1088 | 0.00675 | 12.1855 | 12.1974 | 7.704 | 136.7 |
| 0.71875, 0.75000 | 13.66 | 13.6533 | 0.00794 | 13.7547 | 13.7731 | 7.411 | 131.5 |
| 0.75000, 0.78125 | 15.35 | 15.3404 | 0.00919 | 15.4718 | 15.4995 | 7.118 | 126.3 |
| 0.78125, 0.81250 | 17.21 | 17.2024 | 0.01044 | 17.3705 | 17.4107 | 6.826 | 121.1 |
| 0.81250, 0.84375 | 19.31 | 19.3013 | 0.01183 | 19.5145 | 19.5723 | 6.533 | 115.9 |
| 0.84375, 0.87500 | 21.77 | 21.7604 | 0.01337 | 22.0299 | 22.1126 | 6.240 | 110.7 |
| 0.87500, 0.90625 | 24.83 | 24.8291 | 0.01497 | 25.1719 | 25.2897 | 5.745 | 101.9 |
Rows 14–16 (not used by us): 29.45, 35.94, 46.91 with 0.5501; row 16 has ε = 0.0211 > 1/50 for X < 1.45·10¹³ and needs 0.5936 there (B₁ 50.60).

## E. Consequences for our four conditional results (primes in (x, x+x^0.842], n⁶–(n+1)⁶, n₀(5) = e^{45.584}, Robin2Adic Thm B)
- cert.py's CSRow scales the table's B₁ by (b₄″/0.54)^{(2σ₂−1)/(2−σ₂)}; this is a valid upper bound (K₃ affine in b₄ with positive intercept, K₃ > 1 so B₁ ∝ K₃^{(2σ₂−1)/(2−σ₂)} exactly) and matches the full recomputation (row 12: 22.11). No step failed beyond b₄, so the corrected rows are unchanged and cert.py was not rerun.
- Status: our four results now rest on an audited (still unrefereed) preprint with the corrected rows that cert.py already uses; no step of CS needed by them failed, so nothing is voided. Remaining risk = the published external lemmas of §C.

## F. Files
FULLPROOF.md (step-by-step notes); table1_repro.py + logs/table1_repro.log; divisor_falsify.py + logs/divisor_falsify_{small,2e5}.log (peak 73 MB); chi_falsify.py + logs/chi_falsify.log; m2_check.py + logs/m2_check.log (peak ≈ 45 MB); sens_4thmoment.py + logs/sens_4thmoment.log; f1_repro.py + logs/f1_repro.log; m2rec_repro.py + logs/m2rec_repro.log. All runs under .scratch/memcap.sh 290; peaks ≤ 73 MB.

RH relevance: N7, N10, N4

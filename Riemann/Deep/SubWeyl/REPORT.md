RECORD: |ζ(1/2+it)| ≤ 54.6 t^{27/164} for all t ≥ 3 (before 66.7, Patel–Yang 2302.13444 = J. Number Theory 262 (2024)), using Patel–Yang's lemmas with the square-root split removed and parameters chosen per block; ball-arithmetic certificate. Patel–Yang's own computation drops a factor (2πh)^{17/328} ≈ 1.10, so their parameters give 68.54, not 66.7.

# SubWeyl: lowering Patel–Yang's explicit sub-Weyl constant 66.7 in |ζ(1/2+it)| ≤ A t^{27/164}

Source: Patel–Yang, arXiv 2302.13444v1 (TeX `main.tex`, line numbers "l." below), published as J. Number Theory 262 (2024) 301–334 (same abstract and the same lemma numbering, as cited by 2609.00537). Lemma numbers are those of v1: 1.2 Weyl differencing, 1.4 Yang's k-th derivative test (2301.03165 Lem 2.4/2.5), 2.1–2.3 Poisson / B-process, 3.1 Riemann–Siegel (Hiary), 3.2 S₁, 3.3 S₂ (4th-derivative blocks), 3.4 the ABA³B block lemma, 3.5 S₃ (sum over blocks).

## Theorem (new)
- (a) For all real t ≥ 3: |ζ(1/2+it)| ≤ 54.6 t^{27/164}. Before: 66.7, Patel–Yang Thm 1.1; 307.098, Patel 2021.
- (b) For t ≥ e^{2000}: |ζ(1/2+it)| ≤ 52.36 t^{27/164}.
- (c) The certificate also gives sharper constants on ranges, for example 35.58 on [e^{75}, e^{100}], 45.08 on [e^{175}, e^{200}] and 50.42 on [e^{275}, e^{300}] (`run_final.log`). (a) is sharper than every published Weyl-type bound for t ≥ e^{88.7}. The sharpest of those at large t are HPY's 0.478013 t^{1/6}log t + 3.853165 t^{1/6} − 2.914229 (crossing at e^{87.6}) and Revers' 0.469028 t^{1/6}log t + 3.99399 t^{1/6} − 21.4623 (crossing at e^{88.7}). (a) is sharper than our own 0.40913 t^{1/6}log t + 3.01756 t^{1/6} − 2.920 for t ≥ e^{101.3}. The range constants in (c) beat these bounds much earlier.
- Dependencies: Lemmas 1.2, 1.4, 2.1–2.3 and 3.1 as stated in the refereed Patel–Yang (1.2 and 1.4 are Yang's, from 2301.03165; 3.1 is Hiary's), the display inside the proof of Lemma 3.4 (l. 826–831; audited below), Hiary–Patel–Yang 2207.02366 Thm 1.1 (0.618 t^{1/6} log t, t ≥ 3), and arb.

## Where 66.7 binds (feasibility, done first)
- Their bound: |ζ| ≤ 2(|S₁| + |S₂| + |S₃|) + 1.48t^{-1/4} + 0.127t^{-3/4}. Here S₁ is the trivial range n ≤ θ₃t^{27/82} (Lemma 3.2). S₂ is θ₃t^{27/82} < n ≤ θ₂t^{7/17}, with geometric blocks and the 4th-derivative test (Lemma 3.3). S₃ is θ₂t^{7/17} < n ≤ √(t/2π), with geometric blocks a_k = ⌊h^{-k}√(t/2π)⌋ and Lemma 3.4, summed in Lemma 3.5.
- The constant comes from 163 intervals [e^{L}, e^{L+5}], L = 60…870 (table l. 990–1328, parameters from a stochastic optimiser, no code published). It also uses one fixed-parameter tail for t ≥ e^{875} (l. 946–953), HPY 0.618 for t < 10^{12}, and HPY §3.4 for [10^{12}, e^{60}].
- 66.70 binds at log t ∈ [540, 545] and [560, 565] (rows l. 1189, 1197). The tail gives 66.63.
- Our re-implementation (`sanity_paper.py`, log `sanity_paper.log`) reproduces all 163 rows to within −0.005…+0.115: 66.6989 at l. 1189, 66.7069 at l. 1197, 66.6299 for the tail.
- Breakdown at [540, 545]:
  - S₃ main piece (C₁, from A₅ of Yang's 5th-derivative test): 42.65.
  - S₃ secondary piece (C₂, from B₅): 18.25 as printed, 20.09 corrected (next section).
  - S₂: 5.73.
  - S₁: 0.07.
  - E-terms: below 10^{-17}.
- What each relaxation is worth at log t = 540 (float, `feasibility_proto.py`):
  - Floor-to-geometric block sums: below 10^{-50}, since all blocks have a > e^{200}.
  - R(t) ≤ R(t₁), K(t₁) with t₀: about 0.4. The mismatch between K(t₁) and t₀ costs up to (t₁/t₀)^{3/656} = 1.023 on the C₂ piece; the interval effect on C₁ is about 0.2 %.
  - Per-block θ₁, η₁ and the better of Lemma 3.4 and the 4th-derivative test per block (with the split kept): 68.54 → 65.8.
  - No square-root split in Lemma 3.4 (l. 860; the C₁ and C₂ pieces overlap over about 40 units of log a): 65.8 → 54.2. This is the main gain.
- Nothing forces the binding point. The method's own floor as t → ∞ is about 2C₁μ₃ ≈ 44.5 (float; the optimised tail at e^{2500} gives 49.4).

## An error in Patel–Yang's computation (proof gap, not a counterexample)
- Lemma 3.5, l. 876 and l. 921–923, as printed: Σ_{k≤K} a_k^{-α} < (1−ε)^{-α}(t/2π)^{-α/2}(h^{αK}−1)/(1−h^{-α}) = μ₂(α)h^{αK(t₁)}t^{-α/2}, with μ₂ = μ₁(1−ε)^{-α} and μ₁(α) = (2π)^{-α/2}(1−h^{-αK})/(h^α−1).
- The last "=" is false. Since (t/2π)^{-α/2} = (2π)^{+α/2}t^{-α/2} and (h^{αK}−1)/(1−h^{-α}) = h^α(h^{αK}−1)/(h^α−1), the true right side is (2πh)^{α} times the printed one.
- In C₄ (l. 870) the factor μ₂(17/328) multiplies the C₂ piece, which is therefore too small by (2πh)^{17/328} ≈ 1.10. (μ₂(5/246) on E₂ is negligible.)
- With the printed parameters and the correct μ₂, the rows give at most 68.54 (at [530, 535]), and 82 of the 163 rows (log t₀ from 380 to 785) exceed 66.7. So the published proof does not establish 66.7 on about [e^{380}, e^{790}].
- The tail formula C₅ (l. 883–893, μ₃, μ₄) is right. Its subtracted term √(2π)t₀^{-3/34} in μ₄ is invalid for log t above about 1.7·10^4, which is harmless and dropped here.
- Reviewer verdict on this item: CORRECT. Theorem (a) implies their Theorem 1.1, so the gap is closed.
- Typo-level slips, numerically negligible but repaired in our tail:
  - h₃ has θ₂ for θ₃ (l. 581).
  - h₀ has θ₁ for h/θ₂ (l. 864).
  - C₀ omits θ₃ in 1/(2θ₃t₀^{27/82}) (l. 560).
  - E₁ has 1/(θ₂t₀^{100/697}), but the true bound is θ₂^{-36/41}t₀^{-65/697} (l. 849).
  - q₀ has θ₂^{7/17} for θ₂^{36/41}, which is safe for θ₂ ≥ 1 only (l. 838).
  - The hypothesis θ₂t^{7/17} < a of Lemma 3.4 fails for the last block a_K ≤ θ₂t^{7/17} (l. 535 vs l. 637).
  - A₄(η₂, h₃) is written where Lemma 1.4 gives A₄(η₂, h₃⁴) (l. 581). This one is conservative, since A₄ decreases in h.

## Changed inequalities (old → new; TeX lines of 2302.13444v1)
1. **l. 860 and l. 900–941 (main gain).**
   - Old: |S_f| ≤ α(√X₁₂ + √X₃ + √X₄ + √X₅ + √X₆). Each piece is summed over geometric blocks with one θ₁, η₁ per t-interval.
   - New: per block (a, b], use the whole display l. 826–831, valid for every integer q ≥ 1:
     |S_f(a, b−a)|² ≤ (h−1+q/a)·((h−1)a²/q + (1800/2911)c₃a^{4/5}(qt)^{11/30} + (28800/54481)c₄a^{9/20}(qt)^{61/120} + (8/3)E₄a^{5/2}(tq)^{-1/2} + (9/14)E₅a^{2/3}(tq)^{1/3} + E₆a [+ items 3a].
   - Then |Σ_{a<n≤b} n^{-1/2-it}| ≤ a^{-1/2}√(that), with q = ⌊θa^{36/41}t^{-11/41}⌋ and θ, η₁ chosen per block.
   - Audit (reviewer: CORRECT): before q is chosen (l. 832), the hypothesis θ₂t^{7/17} < a is used only at l. 797–799 and in E₅ (l. 810). We keep E₅ = 0.893(5/π)((h−1)/π)^{1/3}a^{-1/3} + (6/π)h³(h−1).
   - Everything else needs only a < b ≤ ha and 1 < h ≤ 2. In particular x = πrν/t ≤ 1/4 needs r ≤ a, which holds for r < N ≤ (h−1)a; terms with r ≥ N vanish.
   - The display is uniform over prefixes (a, L] with L ≤ b, because it is stated with N ≤ (h−1)a and the same h.
2. **l. 529–552 and 576–635 (fixed split at θ₂t^{7/17}, Lemmas 3.2, 3.3).**
   - New: one block system for t ∈ [e^{L₀}, e^{L₁}]: a_j(t) = ⌊e^{ℓ_j}⌋ with ℓ_j(L) = α_j + L/2 and α_j = −min(j,J)·0.012 − max(j−J,0)·0.05 − ½log 2π, for j ≤ J₂ (down to about t^{0.245}). J and J₂ are fixed per interval.
   - Each block takes the smallest of three bounds: item 1; Lemma 1.4 with k = 4, A₄(η₂,H⁴)·H·N·λ₄^{1/14} + B₄N^{3/4}λ₄^{-1/14} with λ₄ = 3t/(π(Ha)⁴) and N ≤ (H−1)a; and the trivial 2(√b−√a).
   - The rest, n ≤ a_{J₂}, is bounded by 2√(e^{ℓ_{J₂}}+½).
3. **Proof of Lemma 3.4, l. 726–731 and l. 753 (two small gaps, repaired).**
   - (a) The partial sums S_φ(⌊α⌋, L) have up to β−α+1 terms, not β−α. Using N·λ^{1/30} ≤ (β−α)λ^{1/30} + λ^{1/30} and N^{7/8} ≤ (β−α)^{7/8} + 1, two terms are added to the second factor:
     - (1800/451)π^{1/2}h^{3/2}A₅h₅^{1/8}λ_c^{1/30}a^{9/5}t^{-19/30}q^{-19/30};
     - (1800/931)π^{1/2}h^{3/2}B₅λ_c^{-1/30}a^{6/5}t^{-11/30}q^{-11/30}.
     They come from Σ_{r<q}(1−r/q)r^s ≤ q^{1+s}/((1+s)(2+s)) at s = −19/30 and −11/30, with λ_c = 53632√2π⁴/729.
   - (b) The 5th-derivative test runs over (⌊α⌋, ⌊α⌋+L], which reaches below α. In Yang's main case (L ≥ 2337) we have β−α > 2336, and β−α ≤ 2h²(h−1)α, so α ≥ α_min := 2336/(2h²(h−1)) (about 9.4·10⁴).
     - On (α−1, β] the lower bound |φ⁽⁵⁾| ≥ λ₅ still holds, because ν^{-9/2} ≥ β^{-9/2}.
     - The upper bound holds with h₅' = h₅·(α_min/(α_min−1))^{9/2}, because ν^{-9/2} < α^{-9/2}(α/(α−1))^{9/2}.
     - x ≤ 1/4 still holds, and φ⁽⁵⁾ is monotone there (x·P'/P ≤ 0.57 < 9/2).
     - The other cases of Yang's lemma use no derivative bounds.
     - Reviewer: CORRECT.
4. **l. 832–858 (choice of q via q₀, θ₂, t₀).**
   - Q = θa^{36/41}t^{-11/41} decreases between the jumps of a_j(t) and jumps up at them, so only its lower envelope is used. That envelope, θ(e^{ℓ_j(L)}−1)^{36/41}e^{-11L/41}, increases with L (d log/dL ≥ 18/41 − 11/41 > 0), so Q ≥ Q_* := its value at L₀.
   - Q_* > 2 is asserted per block. Then q ≤ Q and q ≥ Q−1 ≥ (1−1/Q_*)Q give fixed factors g = (1−1/Q_*)^{-1}, g^{1/2}, g^{19/30}, g^{11/30}.
5. **l. 581.** We use A₄(η₂, H⁴) as Lemma 1.4 states.
6. **Tail (Lemma 3.5 second part, l. 946–953), for every t ≥ T = e^{2000}.**
   - Parameters: h = 1.006973, h₂ = 1.000094, η₁ = 1.12054, η₂ = 0.87152, θ₁ = 1.51415, θ₂ = 5.82204·10^{14}, θ₃ = 5.40246·10^{-73}.
   - All slips of the previous section are repaired. Lemma 3.4 is applied with θ₂/h₀, where h₀ = h/(1−h/(θ₂T^{7/17})). Items 3a and 3b are added in Lemma 3.4 form.
   - Every t-dependent quantity inside a constant is frozen at T in the safe direction: h₀, θ₂/h₀ < a_K(t)/t^{7/17}, q₀, α², E₁, E₂, h₃, C₀.
   - What is left in t is nonincreasing:
     - C₁μ₃(5/82) and C₂μ₄(17/328) are constants (full geometric series);
     - the S₂ terms carry t^{19/119−27/164} and t^{71/476−27/164};
     - the E-terms carry negative powers of t;
     - E₃·K(t) is bounded with the log-affine envelope K(t) ≤ (3/34·log t − log(θ₂√2π))/log h + 1. That envelope is 2.0·10⁴ at T, above (3/34)/((27/164)log h) ≈ 77 (asserted), so the product with t^{-27/164} decreases.
   - Value: 52.358 = C₁ 44.758 + C₂ 6.120 + S₂ 1.480.
7. **l. 955–960 (small t).** HPY 0.618 t^{1/6}log t is used on [3, e^{75}] only: 0.618·e^{75/492}·75 = 53.983.

## Block domain (derivation)
Line numbers refer to 2302.13444v1 `main.tex`. A block is (a, b] with a < b integers and b ≤ Ha, where H is a fixed ratio bound with 1 < H ≤ 2. In the certificate H_hi ≈ 1.01207 and H_lo ≈ 1.0513, and every block has a ≥ e^{10}. Notation: f(x) = −(t/2π)log x, N = b − a ≤ (H−1)a, S_f(a,N) = Σ_{a<n≤a+N} n^{-it} (l. 663–667).

**A. The proof of Lemma 3.4 up to the display (l. 660–831), step by step, with the hypotheses each step uses.**
1. l. 671–690: g_r(x) = f(x+r) − f(x) on [a, b−r], α = g_r'(b−r), β = g_r'(a), x_ν, φ_r. Uses only a < b and r ≥ 1.
2. l. 694–702: the formula for |φ_r⁽⁵⁾| and x = πrν/t ∈ (0, 1/4]. This needs ½r²/(a(a+r)) ≤ 1/4, i.e. r ≤ a, which the paper takes from h ≤ 2 (l. 700). Only r ≤ N−1 matter: for r ≥ N the sum S_{g_r}(a, N−r) in Lemma 1.2 is empty, and the (positive) bound used for it is still an upper bound. Since N−1 < (H−1)a ≤ a, this needs H ≤ 2 only.
3. l. 707–718: tr/(2πH²a²) ≤ α ≤ ν ≤ β ≤ tr/(2πa²), hence λ₅ ≤ |φ⁽⁵⁾| ≤ h₅λ₅ with h₅ = (76545√2/107264)H⁹. Uses only b ≤ Ha.
4. l. 722: β − α ≤ ((H−1)/π)tr/a², by the mean value theorem on [a, b−r]. Uses only b − a ≤ (H−1)a.
5. l. 725–740: Lemma 1.4 with k = 5 on S_φ, with repairs 3a and 3b. Lemma 1.4 needs only integer a and N, a monotone f⁽⁵⁾ and λ₅ ≤ |f⁽⁵⁾| ≤ h₅λ₅ (l. 129–133).
6. l. 748, 766, 773: on [a, b−r], λ₂ = tr/(πH³a³) ≤ |g_r''| ≤ H³λ₂ and λ₃ = 3tr/(πH⁴a⁴) ≤ |g_r'''| ≤ H⁴λ₃. Uses only b ≤ Ha.
7. l. 751–753: partial summation over ν, using weights |g_r''(x_ν)|^{-1/2} that are monotone in ν. No hypothesis.
8. l. 777–790: Lemma 2.3 (B-process). It needs g_r' decreasing and the bounds of step 6 (l. 467–470), and gives T₁ = 4.686λ₂^{-1/2}, T₂ = (5/π)log(β−α+2), T₃ ≤ (6/π)H³(H−1)(tr)^{1/3}a^{-1/3}, T₄ = 6.
9. l. 793–799: T₂ ≤ (5/π)(log 2 + 0.893(β−α)^{1/3}) ≤ (5/π)(log 2 + 0.893((H−1)/π)^{1/3}(tr)^{1/3}a^{-2/3}).
   - **First use of the threshold:** l. 797 replaces a^{-1/3} by θ₂^{-1/3}t₀^{-7/51}, "from a ≥ θ₂t^{7/17} ≥ θ₂t₀^{7/17}".
   - We stop one step earlier and keep E₅(a) = 0.893(5/π)((H−1)/π)^{1/3}a^{-1/3} + (6/π)H³(H−1) instead of l. 810. This holds for every a.
10. l. 802–813: the BA³B bound for |S_{g_r}(a, b−r−a)|, with E₅(a). It holds for every r ≤ N−1.
11. l. 817–822: sum over r with Σ_{r<q}(1−r/q)r^s ≤ q^{1+s}/((1+s)(2+s)) for s ∈ {11/30, 61/120, −1/2, 1/3, 0} (and −19/30, −11/30 for repair 3a). All lie in (−1, 1]. Valid for any integer q ≥ 1.
12. l. 825–831: Lemma 1.2 (any integer q > 0; a, N integers; l. 100–104) together with b − a ≤ (H−1)a gives the display |S_f(a,b−a)|² ≤ (H−1+q/a)·((H−1)a²/q + … + E₆a). It holds for every integer q ≥ 1 and every block with 1 < H ≤ 2.
- In l. 657–831 the symbols θ₂, θ₁ and t₀ occur only at l. 797, 799 and 810 (checked with rg). The display is uniform over prefixes (a, L], L ≤ b, because it is stated with N ≤ (H−1)a and the same H.

**B. The remaining uses of the threshold (all from l. 832 on, where q = ⌊θ₁a^{36/41}t^{-11/41}⌋ is chosen) and why we do not need them.**
- l. 838–840: q ≥ 1 via q₀ ≤ θ₁a^{36/41}t^{-11/41}, which needs a > θ₂t^{7/17}. Ours: Q = θa^{36/41}t^{-11/41} ≥ Q_* (lower envelope at L₀), and Q_* > 2 is asserted per block and interval, so q = ⌊Q⌋ ≥ 2.
- l. 842: a²/q ≤ θ₁^{-1}a^{46/41}t^{11/41}/(1−q₀^{-1}). Ours: the factor (1−1/Q_*)^{-1}, from q ≥ Q−1 ≥ (1−1/Q_*)Q.
- l. 845: q/a ≤ θ₁/(θ₂^{5/41}t₀^{222/697}), which needs the threshold. Ours: q/a ≤ Q/a = θa^{-5/41}t^{-11/41}, kept exactly.
- l. 848–849: (qt)^{-1/2} through θ₁ − a^{-36/41}t^{11/41} ≥ θ₁ − 1/(θ₂t₀^{100/697}), which needs the threshold (and the printed exponent is too weak). Ours: q ≥ (1−1/Q_*)Q.
- l. 860: the square-root split. Not used in Part C (see E).
- So the display, and hence our block bound, holds on every block with a < b ≤ Ha and 1 < H ≤ 2. That includes blocks with a ≤ θ₂t^{7/17} (for any θ₂) and the last block a_K of Lemma 3.5. There a_K ≤ θ₂t^{7/17} (l. 537), so the paper's own application at l. 899 violates its hypothesis at l. 637.

**C. Where the certificate applies it.**
- Part C: every block (a_j(t), a_{j−1}(t)], j ≤ J₂, with a_j = ⌊e^{ℓ_j(L)}⌋ ≥ e^{10}, ratio ≤ H < 2, and q = ⌊Q⌋ with Q_* > 2 asserted. No θ₂ exists in Part C.
- Part T (t ≥ T = e^{2000}): Lemma 3.4 is used in its stated form, with θ₂ replaced by θ₂⁻ = θ₂/h₀.
  - For every k ≤ K, a_k ≥ a_K > h^{-K}√(t/2π) − 1 > θ₂t^{7/17}/h − 1 ≥ θ₂⁻t^{7/17}. The last step holds because θ₂t^{7/17}(1/h − 1/h₀) = (t/T)^{7/17} ≥ 1. So the hypothesis θ₂⁻t^{7/17} < a holds even for a_K.
  - a_{k−1}/a_k ≤ h/(1 − h^K√(2π/t)) ≤ h₀ ≤ 2.
  - q₀ = θ₁(θ₂⁻)^{7/17}T^{65/697} > 2 and θ₂⁻ > 1, so the exponent 7/17 < 36/41 is safe.

**D. The 4th-derivative alternative.**
- Per block (a, b] with b ≤ Ha: Lemma 1.4 with k = 4, h = H⁴ and λ₄ = 3t/(π(Ha)⁴), on (a, a+N'] with N' = L − a ≤ (H−1)a for each L ≤ b, then partial summation with a^{-1/2}.
- Hypotheses: f'''' = 3t/(πx⁴) is continuous and monotone, and λ₄ ≤ |f''''| ≤ H⁴λ₄ on (a, Ha].
- This is the computation in the proof of Lemma 3.3 (l. 602–620) with h₃ := H. Lemma 3.3 uses θ₂, θ₃ and t₀ only to bound the ratio a'_{r−1}/a'_r ≤ h₃ (l. 595–600) and to sum over r (l. 627–632, R(t) ≤ R(t₁)). We use our explicit ratio bound H and add the blocks one at a time, so no threshold appears.
- In Part T the S₂ blocks keep Lemma 3.3's structure, with h₃ built from θ₃ (as proved at l. 595–600; l. 581 prints θ₂).

**E. One bound per block; the pieces are summed, not alternatives.**
- The display of A.12 is a single inequality. Its right side is (H−1+q/a) times the sum of all summands: the N/q term of Lemma 1.2 plus every term of the BA³B bound, summed over r. All summands hold simultaneously, and none alone bounds |S_f|².
- The paper's C₁ is the square root of the first two summands and its C₂ the square root of the third (l. 856–860). We take the square root of the whole sum, √(F₁·F₂) with F₂ = T1+T2 (the C₁ summands) + T3 (the C₂ summand) + T4+T5a+T5b+T6 (the E₁, E₂, E₃ summands) + T7a+T7b (repair 3a). So the C₁ and C₂ pieces are added inside one bound per block.
- Since √(x₁+…+x_n) ≤ √x₁+…+√x_n, this is at most the paper's split form with the same q and the repairs.
- The certificate's minimum over methods is taken only between three complete bounds for the same block sum: the full ABA³B bound above, the full 4th-derivative bound (both terms added), and the trivial 2(√b−√a). Each is a valid upper bound by itself, so their minimum is too. C₁ or C₂ alone is never used as a bound.

## Certificate
```
uv run --no-project --with python-flint --with numpy python3 Riemann/Deep/SubWeyl/certify.py 54.6
```
- Expected last line: `ALL PARTS OK: |zeta(1/2+it)| <= 54.6 t^(27/164) for t >= 3`. It runs in 241 s, single thread, max RSS 61 MB (log `run_final.log`).
- Part H: [3, e^{75}] by HPY, 53.983.
- Part C: [e^{75}, e^{2000}] by 325 interval evaluations. Starting widths are 25 below L = 300, then 8, 2 on [500, 950], 8, and 50 above 1300; a failing interval is halved. The maximum is 54.59967, on log t ∈ [648, 650].
- Part T: t ≥ e^{2000}, 52.358.
- Rigor:
  - python-flint arb balls at 128 bits. Decimal constants (1.48, 0.127, 4.686, 0.893, Δ = 0.012, 0.05, A) are exact rational balls; θ and η are binary floats, hence exact.
  - Every term is a monomial c·a^p t^e, bounded by e^{pα_j}·e^{(p/2+e)L} with L the ball [L₀, L₁]. This is the exact supremum of a monomial. Negative p get the fixed factor (1−e^{-ℓ_min(L₀)})^p.
  - Block bounds and the rest are multiplied by t^{-27/164} before the supremum is taken. The Riemann–Siegel error is decreasing in t and taken at L₀.
  - The final test is arb "certainly <" against the ball of A. Reported numbers are upper ends.
  - The float search only picks θ, η₁, η₂ and the method per group of 8 blocks; every value used is recomputed in arb.
- Hypotheses asserted on every interval (L₀ ≥ 75 > log 200 holds by construction, for Lemma 3.1): 1 < H < 2 for both block regions; ℓ_min(L₀) > 10; Q_* > 2 for each ABA³B block; α_min > 1 (κ finite).
- Hypotheses that hold by construction: r < N ≤ (H−1)a < a; integer a, b; |f''''| = 3t/(πx⁴) is monotone with ratio ≤ H⁴; φ⁽⁵⁾ is monotone (item 3b).
- Tail asserts: 1 < h₀ < 2, θ₂/h₀ > 1, q₀ > 2, and the K-envelope condition.
- Cross-checks:
  - `crosscheck_terms.py`: each certificate term times a·t^{27/82} equals the matching Lemma 3.4 term (squared constant over α²). The largest log ratio over 200 random cases is 4.5·10^{-13}.
  - `sanity_paper.py` reproduces the paper's table and the μ₂ correction.
  - The same code also certifies 54.55 (`certify.py 54.55`, log `run_5455.log`). That run used 475 evaluations, with intervals near the peak halved down to width 0.5; its maximum was 54.54986 on log t ∈ [638, 639]; it took 376 s and 49 MB. We claim 54.6, which passes with the default widths.

## Where 54.6 binds
- At log t ∈ [648, 650], t ≈ 10^{282}. Lemma 3.4 blocks give 51.47, 4th-derivative blocks 3.13, the rest below 10^{-10}.
- With interval width → 0 the method gives about 54.51 there (54.53 at width 0.5). The bound rises with the log-range covered by Lemma 3.4 and falls as the 4th-derivative and C₂-type terms decay like t^{-0.001}.
- Further gains need sharper lemmas. The ABA³B blocks (Weyl differencing, B-process and Yang's 5th-derivative constants A₅, B₅) carry 94 % of the bound. The 2-dimensional processes of Patel–Yang §4 are another route. The improved B-process of 2609.00537 only touches E-terms, which are negligible here.

## Novelty pass (2026-11-06)
- `python3 .scratch/lst.py 'zeta.*(explicit|estimate)|explicit.*zeta|sub-?Weyl|subconvex|exponent pair' math.NT math.CA --zenodo '"sub-Weyl"'`: 0 title hits in 2026-10/11. Zenodo returns 3 hits: Ramanujan's inequality, S(T) bounds, a short-window certificate. None improves 66.7.
- arXiv API:
  - `all:"sub-Weyl"` and `abs:"27/164"`: only 2302.13444 is explicit (others are GL(2) or Dirichlet asymptotics). 2302.13444 is still v1.
  - `ti:explicit AND ti:zeta`: the newest relevant paper is 2609.00537 (Dhiman–Kadiri–Quesada-Herrera). It improves the explicit B-process and says it "could refine" Patel–Yang, but gives no new constant.
- OpenAlex: 6 papers cite the JNT version (zero-free regions 2212.06867, Phragmén–Lindelöf 2502.13282, S(T) 2412.15470, Deuring–Heilbronn 2410.06082, log-free density 2405.12545, a 2026 book chapter). All use 66.7; none improves it, and none notes the μ₂ slip.
- HAL: "sub-Weyl" gives 0 hits; "zeta explicit critical line" gives nothing relevant.
- Repo: no earlier entry. TREE N9 tracks 66.7; R376 harvested the relaxations.

## Caveats
- Conditional on the cited lemmas. Lemma 3.4 is used through a display of its proof (item 1, audited), plus the repairs of item 3.
- The JNT version of record is CC-BY, but its host returns 403 to scripts, so all checks are against arXiv v1. The abstract (66.7) and the lemma numbering cited by 2609.00537 match v1. If the journal version had silently fixed μ₂ with new parameters, the error claim would apply to v1 only. The record (54.6 < 66.7) stands either way.
- The certificate trusts python-flint/arb.
- The reviewer checked the argument, not the code line by line. `crosscheck_terms.py` checks the term bookkeeping against Lemma 3.4.

## Review
The reviewer subagent ran three rounds, with the evidence pasted inline (`review.log`).
- Round 1: item 1 (the μ₂ error) CORRECT; item 2 (no split, per-block parameters) CORRECT, subject to the hypothesis audit; the derivative-domain repair, the uniformity in t and the tail monotonicity were GAPs.
- Round 2: the audit (P1) and the enlarged derivative domain (P2) CORRECT. Two GAPs remained: Q must be handled through its lower envelope, and K(t) through its log-affine envelope. Both were already how the code works; the specification was made explicit.
- Round 3: P3′ and P4′ CORRECT. Overall: **COMPLETE for the conditional 54.6 bound** (max{53.983, 54.59967, 52.358} < 54.6), conditional on the cited lemmas and a correct arb implementation. Two wording fixes were applied: the Q-jump description, and the 54.55 claim now rests on its own run log.

RH relevance: N9 (explicit subconvexity on the critical line; quantitative record)

## Files
- `certify.py`: the certificate (Parts H, C, T). Run with arguments `L0,L1` for single intervals or `tail`.
- `run_final.log`: certificate run with the final code (54.6).
- `run_5455.log`: the 54.55 run.
- `review.log`: the three reviewer verdicts.
- `crosscheck_terms.py`: float check of the term bookkeeping against Lemma 3.4.
- `feasibility_proto.py`: float point estimates, combined vs split square root (feasibility numbers above).
- `sanity_paper.py`, `sanity_paper.log`: re-derivation of Patel–Yang's table and of the μ₂ correction.

## Suggested lines (orchestrator)
- DISCOVERIES: |ζ(1/2+it)| ≤ 54.6 t^{27/164} for all t ≥ 3 (before 66.7, Patel–Yang); block-wise unsplit ABA³B bound, ball-arithmetic certificate.
- REFUTATIONS: Patel–Yang 2302.13444 Lemma 3.5 drops a factor (2πh)^{17/328} ≈ 1.10 in μ₂, so their parameters give 68.54, not 66.7, on [e^{380}, e^{790}]. Theorem 1.1 is restored by our 54.6.

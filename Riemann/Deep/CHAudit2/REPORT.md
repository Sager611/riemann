AUDIT PASSED: the nine technical lemmas of Chirre–Helfgott 2511.14736v1 and its square-free chain gould/qmediano → nopgik/gopnik → bataclar (c = 0.0134) → mertosimp (11.389414) re-derived line by line, no error; 1 harmless gap in a statement (qmediano), 8 new misprints; the cited Lemma 5 (Carneiro–Littmann) is used exactly as proved; still open (out of scope): Platt's private residue sums above 1.9·10^6.

# Referee audit 2 of Chirre–Helfgott arXiv 2511.14736v1 (node N10): technical lemmas and square-free propositions

Agent chaudit2, 11 Nov 2026, 17:08–17:55 UTC. Read-only. Sources: e-print `Moebius.tex` (all line numbers below are its lines), Carneiro–Littmann arXiv 1309.1527v1 (= zbMATH06384942), Mossinghoff–Oliveira e Silva–Trudgian arXiv 1912.04972v2 (= zbMATH07298460). Builds on chaudit (`Riemann/Deep/CHAudit/REPORT.md`) and R437 (random and corner-case tests of (a)); neither is repeated here. Scripts are in this directory, their output in `logs/`. No Lean, nothing committed.

## 0. Answer first

| question | answer |
|---|---|
| Any error? | No. Every inequality, constant and hypothesis of kolobrz, shiftbounded, benvenuto, paulin, pommedupe, arborio, barbar, artanor, zetinvbound, and of §8 (difintsq → andalas → gould, qmediano → nopgik, gopnik, bataclar → metamert, mertosimp) re-derived. All hold. |
| Gaps | One, harmless: the hypotheses of Prop. qmediano (l.4540–4541) do not cover v ∈ (√x, v₁], which the tail ∫₀¹ uses. The proof's (4/3)κx^{1/4} needs κ ≥ κ₋, which holds in nopgik (6.738093 > 0.570591), and the input bound there holds for all v ≥ 1 anyway. |
| Misprints | 8 new plus 2 already known (R437), none load-bearing (§3). |
| Rounding | Every rounded constant rounds in the safe direction (§2.4); none is rounded down. |
| Cited "Lemma 5" | Carneiro–Littmann, J. Fourier Anal. Appl. 19 (2013) 967–996, arXiv 1309.1527v1, Lemma 5 "Integral representation". With hypothesis λ ≥ 0 (used with λ = \|λ\|/2 > 0) it is eq. monthmay verbatim, with the same K_ν, E_ν and b. |
| Margins | c = 0.0134: nopgik gives 0.013383; 11.39 survives any c ≤ 0.013616; the stated hand-over 3.61·10^17 needs c ≤ 0.013406 (for larger c the hand-over moves, and the garadiol branch reaches 1.22·10^18). |
| Our records | Unchanged. MertensCH does not depend on the square-free chain at all (§4): with c = 1 from the elementary Lemma oldtrivial, the leading term 3x/(π·10^10) stays, and only 11.39 becomes 14.251. |

## 1. Technical lemmas (item a)

**1.1 Lemma kolobrz** (l.2567–2577, proof l.2683–2695). Statement: G meromorphic on R = (−∞,1] + i[−T,T], G(s)x₀^s bounded on ∂R ∪ L, x > x₀ ⇒ ∫ over [1−iT, 1+iT] = ∫_{C∞} + Σ_{ρ∈R} Res (sum = limit over Re ρ > σ_n, l.2578). Check: boundedness excludes poles on ∂R ∪ L; each rectangle [σ_n,1]×[−T,T] has finitely many poles. The counter-clockwise boundary gives exactly C_n (bottom edge leftward, left edge up, top edge rightward). The left edge is at most 2T·sup|Gx₀^s|·(x/x₀)^{σ_n} → 0, and the horizontal edges are dominated by (x/x₀)^{Re s}, which is integrable on (−∞,1]. ✓

**1.2 Prop. shiftbounded** (l.2697–2731). Statement eq. gruman. Check: on s = 1 − t + iξT, (s−1)/(iT) = ξ + it/T, so |Φ_λ| ≤ πt/(8T) (benvenuto) and (1/2π)(π/8)/T = 1/(16T) ✓. The hypotheses of kolobrz hold for G = Φ_λ((s−1)/iT)F, x₀′ = (1+ε)x₀ < x, L₀ = {σ_n < 1 + λT/2π}:
- on L₀, Im z = (1−σ_n)/T ≥ (1−σ_{n₀})/T > −λ/2π, and benvenuto gives a uniform bound;
- on the right edge z ∈ [−1,1] is real and Φ_λ has no real poles (Im = −λ/2π ≠ 0), so Φ_λ is continuous there;
- on the horizontal lines, linear growth times (1+ε)^{1−t} stays bounded.

The poles of Φ_λ((s−1)/iT) are s = 1 + λT/2π + 2niT; only n = 0 lies in R, and only for λ < 0 ✓. The zeros of Φ_λ((s−1)/iT) closest to R are the corners 1 ± iT, which lie on ∂R, so no pole of F is cancelled inside. Listing a non-pole adds a residue 0, which is harmless. ✓ Test: R437 end to end (F = 1/ζ, T = 20, 30), max 0.864.

**1.3 Lemma benvenuto** (l.2160–2191).
- Poles of coth(πz/2i + λ/4): z = 2n − iλ/2π, residue (sgn λ/4)(2i/π) = i sgn λ/2π ✓.
- girasol2 re-derived from girasol3 via tan αi + cot β = cos(β − αi)/(cosh α sin β) ✓. Its zeros are cos(πz/2) = 0, i.e. z odd; the denominator does not vanish there because λ ≠ 0 ✓.
- Im z ≥ c > −λ/2π ⇒ Re w ≥ πc/2 + λ/4 > 0 ⇒ coth is bounded ✓.
- Φ_λ(±1+it) = (sgn λ/4)(tanh(πt/2 + λ/4) − tanh(λ/4)) (coth(s ± iπ/2) = tanh s). Its t-derivative is (π/8)sech² ≤ π/8 and its value at t = 0 is 0 ✓.

**1.4 Lemma paulin** (l.3061–3081). For t < min(σ,1): w(t) < coth(δ(t−σ)) − tanh(δ(t−σ)) < 0, since tanh increases and t < 1, and coth x < −1 < tanh x for x < 0 ✓. w′ = δ·coth′ < 0, so |w| increases ✓. −w = coth(δ(σ−t)) − tanh(δ(σ−1)) < 1/(δ(σ−t)) + 1 + 1, because coth x − 1/x < 1 (its derivative 1/x² − 1/sinh²x is positive and the limit is 1) ✓.

**1.5 Lemma pommedupe** (l.3114–3159).
- Residue at −2n: w is holomorphic there (its only real pole is σ > −2 ≥ −2n) and the zero of ζ is simple, so Res = w(−2n)x^{−2n−1}/ζ′(−2n) ✓.
- koshmar follows from the functional equation (n = 1: ζ′(−2) = −ζ(3)/4π²) ✓.
- Consecutive ratio: |1/ζ′(−2n−2)|/|1/ζ′(−2n)| ≤ (2π)²ζ(3)/12 = 3.955. With x^{−2} ≤ 1/4 and |w(−2n)| decreasing (paulin), the terms decrease and alternate, so |Σ| ≤ the first term ≤ (1/(δ(σ+2)) + 2)(2π)²x^{−3}/ζ(3) ✓.
- So the lemma holds with O*(…) in place of "=". Its use in metamert (l.4782: σ ≥ −1, δ ≤ 1/12) gives (7/6)(2π)²/ζ(3) = 38.316 < 39 ✓.

**1.6 Lemma arborio** (l.5513–5600) and the cited Lemma 5.
- Identification: zbMATH06384942 is Carneiro–Littmann, "Entire approximations for a class of truncated and odd functions", JFAA 19 (2013), arXiv 1309.1527v1. Its theorem counter is shared (its TeX l.20–28): Thms 1–4 come before, so "Integral representation" (its l.441–457) is Lemma 5. Also, Moebius's "Theorem 1" is their two-sided Thm 1, used in carlitfou with c = 0, δ = 2: (1/δ)(1−e^{−λ/δ})/((λ/δ)(1+e^{−λ/δ})) = tanh(λ/4)/λ ✓.
- Statement of their Lemma 5: for λ ≥ 0, K_λ(x) − E_λ(x) = (sin πx/π)∫₀^∞(b(λ+w) − b(λ))e^{xw}dw for x < 0, and (sin πx/π)∫_{−∞}^0(b(λ) − b(λ+w))e^{xw}dw for x > 0, with b = 1/(1+e^w). Their K_λ = (sin πz/π)Σ_{n≥1}(−1)^n(E_λ(n)/(z−n) − e^{−λn}/z) equals eq. ombroso (E_λ(n) = e^{−λn} for n ≥ 1; Σ_n = Σ_{n≥1}, l.2079). E_λ is the same. Used with λ = \|λ\|/2 > 0: this is eq. monthmay exactly. I re-derived their x < 0 case (geometric expansion of b, Laplace integral for 1/x).
- Integration by parts monthmay → monthmayending: the boundary terms vanish ✓.
- Bound on the function: |b′| ≤ 1/4 ⇒ |K−E|(v) ≤ |sin πv|/(4πv²), and with v = 2 sgn(λ)u this gives |φ̂_λ − I_λ|(u) ≤ 1/(16πu²) for all u ≠ 0 ✓.
- Bound on the derivative: |g| ≤ 1/(4|v|), |g′| ≤ 1/(4v²) and (sinc)′ = (cos πv − sinc)/v give κ(v)/(4v²). The chain-rule factor 2 gives κ(v)/(8u²) with |v| ≥ 1, and κ < 12/11 gives 3/(22u²) ✓.
- κ = max(|cos πu|, |h|) holds with equality. On [1, 3/2] the printed chain is false (§3.5), but −1 ≤ h ≤ 2/π, so κ ≤ 1. For u ≥ 3/2, |h| ≤ (1 + 4/(3π/2)²)^{1/2} = 1.08634 < 12/11 ✓.
- Signs and monotonicity on (−½, ½): sinc > 0 and monotone; g₁ > 0, g₁′ > 0 (v < 0); g₂ < 0, g₂′ > 0 (v > 0) ✓.
- Jump −sgn λ: K is continuous and I_λ jumps by sgn λ ✓. f(±½) = 0 because K_ν(±1) = E_ν(±1) ✓.
- Tests: monthmayending (Gauss–Laguerre) against the Lerch-Φ series of K_ν agrees to 4.0·10⁻¹⁵ at 33 points (chk_artanor.py). R437: 1.2·10⁻¹¹ against the direct Fourier integral, max 16πu²|f| = 0.999918.

**1.7 Lemma barbar** (l.5602–5631).
- The substitution and the integration by parts eq. 10_32pm re-derived.
- g′(x) = 2e^x(x−β)/((β−2)x³), so the minimum is at β ✓, and the terms e^{−αu₀}/(αu₀) cancel ✓.
- At β = 6.4: β/(β−2) = 16/11 and Ei(6.4) − e^{6.4}·5.4/(6.4·4.4) = 2.52424054 ≤ 2.525 ✓.
- Test (chk_a.py): 364 pairs (α, u₀) with α ∈ [10⁻⁶, 20] and u₀ ∈ [−100, −½]. Max LHS/RHS = 0.99603 at α = 3, u₀ = −2, where −αu₀ ≈ β and the Ei bound is tight.

**1.8 Lemma artanor** (l.5633–5714; the copies at l.5716, 5837, 5892 sit in comment blocks l.5715–5889 and l.5891–5998).
- L1 bound: ∫_{1/2}^∞ ≤ e^{−α/2}/8π ✓, and ∫_{u₀}^{−1/2} ≤ (1/16π)·barbar ✓. On [−½,½], the monotone pieces give (|f(0⁻)|e^{α/2} + |f(0⁺)|)/2 ≤ e^{α/2}/2, using |f(0⁻)| + |f(0⁺)| = 1 ✓.
- TV bound: |(e^{−αu}f)′| ≤ (α/16π + 3/22)e^{−αu}/u² ✓. On [−½,½] the total is jump 1 + |f(0⁺)| + (2e^{α/2} − 1)|f(0⁻)| ≤ 2e^{α/2}. f(0) equals a one-sided limit (I_λ(0) = 1), so it adds no extra variation ✓.
- Constants: c₀, c₁, c₂ at ½ are 0.75460337, 3.39608581, 0.21281574 (< 16/21, 17/5, 2/9). At 1/6: c₀ = 0.63772139 ≤ 2/3 and c₁ = 2.83567808 ≤ 2.85. All three increase on (0, ½] (step 0.001) ✓.
- Test end to end (chk_artanor.py): 162 cases, λ ∈ ±{10⁻⁴, 0.5, 6}, α ∈ {0.02, 1/6, ½, 2}, u₀ ∈ {−40, −10, −3, −1, −½, −0.3, 0.2}. Max L1/bound = 0.489, max TV/bound = 0.877 (λ → 0, α = 0.02, u₀ = −40).

**1.9 Lemma zetinvbound** (l.6138–6211).
- Re log Γ, the identity Im z·arg z = π|Im z|/2 − |Im z|arctan(Re z/|Im z|) and the assembly to |ζ(s)| = (|1−s|/2π)^{½−σ}(2|sin|/e^{π|t|/2})e^{r₁}|ζ(1−s)| re-derived ✓.
- r₁ ≥ −Re z because arctan x ≥ 0, and Re z/(12|z|²) > |remainder| on Re z ≥ ½ ✓.
- Second bound: |sin(πs/2)| ≥ e^{π|t|/2}(1−e^{−π|t|})/2 ✓; (1−y)⁴ ≥ 1−4y; (1 + a)(1 − a/5) = 1 + a(4−a)/5 ≥ 1; t²e^{−πt} ≤ e^{−π} < 1/20 for t ≥ 1 ✓.
- Stirling input [GR 8.344]: the constant √2/180 holds. On a grid over Re z ≥ 0, sup |R₂(z)||z|³ = 0.00393 (at arg z = π/2, |z| = 0.78), below √2/180 = 0.00786. It is not load-bearing anyway: Stieltjes' sec⁴ bound 1/90 also suffices, and the needed (Re z/12|z|² + Re R₂)/(Re z/12|z|²) has minimum 0.93 > 0. R437 also tested 20,000 random s: max 0.99958.

## 2. Square-free chain (item b)

**2.1 Setup**
- **Lemma difintsq**, the convolution identity (l.3857–3904). Q(x) = Σ_{kd²≤x}μ(d) = Σ_k M(√(x/k)), and ∫₀^x M(√(x/u))du = Σ_nμ(n)x/n² = x/ζ(2) by absolute convergence ✓. The cells [k−½, k+½] for K < k ≤ x+1 tile [K+½, ⌊x⌋+3/2], and M(√(x/u)) = 0 for u > x ✓. Test: identity exact at 200 random x ≤ 10⁶.
- **Lemma sugr** (l.3925–3976). On [k−½,k] the difference is Σ_{√(x/k)<n≤√(x/u)}μ(n), whose integral is Σμ(n)(x/n² − k + ½) with each length ≤ ½; similarly on [k, k+½] ✓. Test: 72,000 cells (x ∈ [10⁸, 10¹³], k ≤ 4000), max 0.99995 (sharp: a square-free n at a cell edge).
- **Lemma trapeze** (l.4073–4092). |M(t)−M(t′)| ≤ |Q(t)−Q(t′)| ✓. f(u) = u^{−½} − (u+½)^{−½} is convex (f″ = ¾(u^{−5/2} − (u+½)^{−5/2}) > 0), so the trapezoid rule over-estimates ✓. Test: the same cells and all 7 edesmo pairs, max 0.677.
- **Prop. andalas** (l.4094–4119). Σ₁ by trapeze (telescoping); Σ₂ by sugr (disjoint cells inside (0, √(x/(K′+½)))) ✓. Test end to end, with exact R(x) and exact main terms (chk_b2.py): 156 cases (x, K, K′, pair), including the proofs' K and K′. Max |D|/bound = 0.030.
- **Lemma oldtrivial** (l.4123–4148). |M(x)| ≤ x/4 + 3 by pairing odd n with 2n ✓. Bootstrapping gives |R| ≤ 0.78342√x + 1.26135x^{1/4} + 1.5 ≤ √x from x = 2400 (it fails at 2000, so the threshold is right). Test: sup over all real x ∈ [1, 10⁷+1) of |R(x)|/√x ≤ 0.679.
- **Lemma ramando, Cor. legolas** (l.4150–4190). Midpoint convexity and tangent concavity ✓; the three side conditions hold at N = 3, 4, 38 ✓. Test, all N ≤ 10⁶: max (LHS − RHS) of synthe = −0.0296; of synthedos = 0.71078 ≤ 5/7 (at N = 3); without 5/7 the max is −0.0022 for N ≥ 37, and N = 36 gives +0.0033, so 37 is sharp.
- **Lemma edesmo** (l.4328–4351). chaudit checked it exactly over full periods. Directly on Q(t), t ≤ 10⁷ (chk_b3.py): sup(Q(t₂) − Q(t₁) − c₁(t₂−t₁)) = 1, 3/2, 8/3 (attained), 3.92, 5.51, 6.85, 8.42 against c₂ = 4.56, 7.72, 12.18, 20.02 ✓.

**2.2 Prop. gould** (l.4356–4536). Hypotheses: |M(v)| ≤ εv + κ√v for all v > 0; (c₁, c₂) short-interval pair with 0 < c₁ ≤ 1; c₃ = ½(1/2ζ(2) − c₁/4) > 0. Re-derivation:
- The K′-terms of andalas plus oldtrivial (Q(y) ≤ y/ζ(2) + √y) give αK′ + β/√(K′+½) with α = c₂ and β = 2c₃√x. The minimum of αυ + β/√υ is 3α^{1/3}(β/2)^{2/3} at υ₀ = (β/2α)^{2/3}. With K′ = ⌈υ₀ − ½⌉ the value is ≤ 3c₂^{1/3}c₃^{2/3}x^{1/3} + c₂/2, and ½(x/υ₀)^{1/4} = ½(c₂/c₃)^{1/6}x^{1/6} ✓.
- The −c₂K + c₂/2 left over is negative because K ≥ 3; the proof drops it (l.4504) ✓.
- (a) a = (8κ/3)x^{1/4}, b = (c₁/4)√x. The minimum of ay^{3/4} + b/√y is (5a/2)(2b/3a)^{3/5} at y₀ = (2b/3a)^{4/5} = x^{1/5}/(16κ/c₁)^{4/5}. This gives the x^{2/5} coefficient (5/3)c₁^{3/5}(κ/2)^{2/5} and the x^{3/5} coefficient 4ε/(16κ/c₁)^{2/5} ✓.
- Side conditions of (a): y₀ > 5/2 ⇔ x > 2¹¹5⁵(κ/c₁)⁴ = 64·10⁵(κ/c₁)⁴, which gives K ≥ 3. υ₀/y₀ = ((16κ/c₁)⁶c₃⁵c₂^{−5}x)^{2/15} ≥ 1 gives K′ ≥ K. The 5κ/7 term may be dropped when y₀ > 73/2 ⇔ x > 2¹¹·73⁵(κ/c₁)⁴, and 2¹¹·73⁵ = 4.2457·10¹² ≤ 4.25·10¹² ✓.
- (b) the minimum of 4εy + (c₁/4)/y is 2√(εc₁) at y₀² = c₁/16ε > 73/2 (because ε < c₁/584), so K ≥ 37 ✓. (8/3)y₀^{3/2} = (1/3)(c₁/ε)^{3/4} ✓. υ₀ ≥ y₀² ⇔ x ≥ (c₂/c₃)²(c₁/16ε)³ ✓.
- Test, proof-as-computation (chk_b.py Part 2): the bound obtained from the proof's own K and K′ with exact sums (andalas + oldtrivial) against the closed form. 10,012 random (a)-cases and 9,991 (b)-cases: (c₁, c₂) from edesmo or random, κ ∈ [10⁻³, 10²], ε ∈ {0} ∪ [10⁻¹², 10⁻¹] or up to 0.999999·c₁/584, x from (1+10⁻⁹)× threshold to 10¹²× threshold. K ≥ 3 (37) and K′ ≥ K never failed. Max ratio 0.99992 for (a) and 0.99994 for (b): the closed forms are sharp up to integer rounding and never exceeded.

**2.3 Prop. qmediano** (l.4538–4583).
- The K of gould (a) with κ₋ in place of κ.
- For k ≤ K and u ∈ [1, K+½], v ∈ [√(x/(K+½)), √x] ⊂ [v₀, v₁]. The v₀ condition: K+½ < y₀ + 1 ≤ x/2v₀² + x/2v₀² by (2v₀²)^{5/4}/(16κ₋/c₁) < x and x > 2v₀² ✓.
- Legolas with ε = 0, minus the [0,1]-part (4/3)κ₋x^{1/4} > (5/7)κ₋x^{1/4}, so the 5/7 term is dropped ✓. The tail ∫₀¹ ≤ 2ε√x + (4/3)κx^{1/4}, which needs the gap fix in §0.
- Test: 9,997 random cases, max 0.999994; the v₀ check never failed.

**2.4 Lemmas nopgik, gopnik, Cor. bataclar** (l.4584–4639). Exact recomputation with 40 digits (chk_b.py):

| quantity | paper | exact |
|---|---|---|
| nopgik coefficients (qmediano: κ₋ = 0.570591, v₀ = 33, v₁ = 10¹⁶, ε = π/2·10⁷, κ = 6.738093, (16/25, 114/25)) | π/10⁷, 0.772103, 0.860845, 8.98413, 0.999 | 3.14159·10⁻⁷, 0.77210204, 0.86084423, 8.984124, 0.99830 |
| nopgik lower conditions | < 10¹⁸ | 4.04·10⁶, 121.1, 1043, 2178; x ≤ v₁² = 10³² |
| \|R\|/x^{2/5} at the endpoints | ≤ 0.84491 | 0.84490422 (true sup of the bound on [10¹⁸, 10³²]: 0.84443) |
| 0.845x^{2/5} ≤ 0.0134√x | for x ≥ 10¹⁸ | holds from 9.943·10¹⁷ |
| c given by nopgik | 0.0134 | 0.0133833 |
| gopnik (a) coefficients (ε = π/2·10¹⁰, κ = 11.350514, (768/1225, 9458/1225)) | 6.51·10⁻¹¹, 2.5222, 1.042, 1.086 | 6.50752·10⁻¹¹, 2.5221914, 1.0414801, 1.0858089 |
| gopnik (a) conditions | | 6.88·10¹¹, 2.1·10⁻⁵; omitting 5κ/7 allowed from 4.56·10¹⁷ |
| 2.5222 + 1.042x^{−1/15} + 1.086x^{−7/30} at 10³² | ≤ 2.53 | 2.529853 |
| \|R\|/√x on [10³², 10⁵⁰] | ≤ 0.00161 | 0.00160283 (exact coefficients 0.0015963) |
| gopnik (b) | 1.99·10⁻⁵, 1.042, 6.01·10⁷, 1.1 | 1.98473·10⁻⁵, 1.04148, 6.00791·10⁷, 1.08581 |
| (b) conditions | ε < c₁/584; x ≥ (c₂/c₃)²(c₁/16ε)³ | yes; 1.71·10²⁹ |
| \|R\|/√x above 10⁵⁰ | ≤ 3.9·10⁻⁵ | 3.89101·10⁻⁵ |
| bataclar from MOST on [10⁸, 10¹⁸] | 1.12543x^{1/4} ≤ 0.0134√x | from 4.976·10⁷ |

Inputs:
- MOST Thm. 2 (arXiv 1912.04972v2, post-submission version, l.158–160): −1.12543x^{1/4} < R(x) < 1.11653x^{1/4} for 0 < x ≤ 10¹⁸, as quoted. Our sieve attains the upper constant 1.116522 at x = 43, and gives max |R|/x^{1/4} = 0.642 at 120 random x in [10⁷, 5·10¹³].
- Lemma hurst (refereed; chaudit: max 0.567105 on [33, 10⁹]).
- Cor. mertensimplon: garancedos (π/2·10⁷, 6.738093 = 2·3.36904620179490 rounded up) and garance2 (π/2·10¹⁰, 11.350514), both stated "for x ≥ 1" (l.3479), so "for all v > 0" (M = 0 below 1) and "for v > v₁" hold.
- No circularity (chaudit).

**2.5 → Cor. metamert, Cor. mertosimp, Cor. mertensimple.**
- Prop. mainsqfprop, eq. adanto re-derived with y₀ = −ρT/π, ρ = L/2.
- The constant 2.9: c₁ + 2π/3T ≤ 2.835678 + 0.0419 < 2.9 (artanor at α = π/T ≤ 1/6). Also 4(c₂ + (1/11 + 1/16)/T) ≤ 0.902 < 1 ✓.
- metamert: L ≥ 3 gives 1/3 + 4/9 + 1/3 = 10/9. The trivial zeros give 38.316/x² < 39/x² ≤ (39/64)x/T³ for x ≥ 4T. The final inequality holds for T ≥ 50: the 1/T² coefficient may be up to 0.9549 against 0.8727 used, slack factor 1.094 ✓.
- mertosimp: 2.9·0.0134 = 0.03886 ≤ 0.0389; 11.350514 + 0.0389 = 11.389414; the 2 = |1/ζ(0)| is absorbed for x ≥ 1.165·10⁷; e³(cT)² = 3.60656·10¹⁷ < 3.61·10¹⁷ ✓.
- Branch arithmetic of the proof of Cor. mertensimple (l.4811–4819) ✓.
- End-to-end test from chaudit's data (`CHAudit/logs/eft.log`, exact M(x) to 10⁹): max |M − P₀| = 0.452, 0.369, 0.380 times δx at T = 10³, 10⁴, 10⁵. metamert's leading term is (3/π)/(T−1)·x ≥ (6/π²)δx = 0.608δx, so the square-free gain π/2 → 3/π is consistent with the data, with 26 % margin.

## 3. Misprints and gaps found (none load-bearing)
1. pommedupe, statement l.3118–3121: "=" should be O*(…) (known from R437).
2. pommedupe l.3130: "w has no zeros to the left of σ"; the residue formula needs "no poles" (true: the only real pole is σ).
3. pommedupe l.3144: "|1/ζ′(2n)|" should be |1/ζ′(−2n)|.
4. pommedupe l.3153: the denominator ζ′(−2) < 0 should be |ζ′(−2)|.
5. arborio l.5578: the chain h(u) ≤ 2 sin πu/(πu) ≤ 4/(3π) is false on [1, 3/2], where sin πu ≤ 0 (known from R437). Correct: −1 ≤ cos πu ≤ h(u) = cos πu + 2|sin πu|/(πu) ≤ 2/π, so κ ≤ 1.
6. arborio l.5585: tan(πu)/(πu) > 1 gives monotonicity of sin(πu)/(πu) only on (0, ½); on [½, 1) it is immediate, since cos πu ≤ 0.
7. artanor l.5658: "the second integral" should be the third (∫_{1/2}^∞).
8. artanor l.5703: "continuous at ½" should be at −½ (f(−½) = 0 is what is used).
9. artanor l.5712: "Since κ_α is increasing, c₀, c₁, c₂ are increasing" is incomplete, because the e^{−α/2} terms decrease. The conclusion still holds: c₀′ ≥ e^{α/2}/4 − e^{−α/2}/16π > 0 and c₁′ ≥ e^{α/2} − (3/22 + α/16π)e^{−α/2} > 0; checked numerically as well.
10. andalas l.4107: "Assume K, K′ ≤ x" lacks its period, and the case K′ > x is not mentioned (Σ₂ is then empty; harmless).
11. Gap, qmediano l.4540–4541: as stated, the tail needs |M(v)| ≤ εv + κ√v also on (√x, v₁]; in general the term should be (4/3)max(κ, κ₋)x^{1/4}. Harmless in nopgik.

## 4. Dependency map (item 3)

| our entry | what it takes from 2511.14736 | effect of a failure in the square-free chain or the technical lemmas |
|---|---|---|
| MertensCH record \|M(x)\| ≤ 9.77·10⁴x^0.9382, 1.54·10⁴x^0.9409 (RECORDS N10, TREE N10, discoveries/mertens-power-saving.md, Riemann/Deep/MertensCH) | Cor. mertensimple as a black box on [1, X₀], X₀ = 10^{240.35} … 10^{2826.8} (logs/cert_ch_bp.log); c_CH = max(CH(1), CH(X₀)/X₀^θ) by convexity; the binding term is 3X₀^{1−θ}/(π·10¹⁰) | **None.** (i) Whole chain beyond oldtrivial lost: c = 1 (Lemma oldtrivial, elementary, checked to 10⁷) in Cor. metamert gives \|M(x)\| ≤ 3x/(π·10¹⁰) + 14.251√x for all x ≥ 1 (Cor. mertensimplon below e³T² = 2.0·10²¹), with the same leading coefficient, and CH(1) = 14.25 ≪ 1.54·10⁴. (ii) MOST lost: c = 0.019172, constant 11.41 (hand-over 7.4·10¹⁷ < 1.25·10¹⁸). (iii) Platt's sum at 10¹⁰+1 may be as large as about 7,690 (claimed 5.675256) with the record unchanged: it enters only CH(1) = 3/(π·10¹⁰) + 2Σ + 2.9c + 2 ≤ 1.54·10⁴. (iv) The square-free gain lost entirely (mainsqfprop wrong): Cor. mertensimplon's πx/(2·10¹⁰) multiplies c_CH by at most π²/6 = 1.645 at the same X₀ (estimate), θ unchanged. (v) Thm. mainthmA lost: back to Mertens78 (chaudit §5). |
| The constant 11.39 itself (quoted in TREE N10, R428, Harvest5) | gould/qmediano, nopgik/gopnik, bataclar, MOST, Hurst, Platt at 10¹⁰+1 | 11.39 holds for any c ≤ 0.013616; 11.41 without MOST; 14.251 with oldtrivial only. The technical lemmas enter only the 1/T² terms (slack factor 1.094, relative size 10⁻¹⁰), the TV constant behind 2.9 (c₁(1/6) = 2.8357 ≤ 2.85) and the contour shifts (qualitative). |
| RobinCH, OddRobin, RamRepair, Buthe78 verdict, NOGO R307/R379/R389 | nothing: rg finds no use of 2511.14736 there; they use 2512.15709 (Cor. psiexpl) only | unaffected by this audit |
| "MertensCH small-x part" | below ~10⁷⁰ even \|M(x)\| ≤ x ≤ 1.54·10⁴x^θ suffices; on [10⁷⁰, X₀] only the linear coefficient 3/(π·10¹⁰) matters | as above |

Open, unchanged: Platt's private residue sums at 10⁷+1, 10⁸+1, 10⁹ and 10¹⁰+1 (chaudit's extrapolation agrees to 0.04–0.11 %). They set 11.39 but not MertensCH (row 1, iii).

## 5. Referee verdict in plain words
The technical core of the Möbius paper holds up. I redid every estimate behind the contour shift, the extremal function, the zeta bound and the square-free-number counting, and checked each one numerically at random and extreme parameters; nothing fails. The authors cite another paper for one key formula, and that formula is exactly what the other paper proves. The only blemishes are about ten typos and one sloppily stated hypothesis, which does not matter where it is used. The constant 11.39 has a thin but real safety margin, and our own Mertens record would survive even if the whole square-free part of the paper were wrong. The one thing still on trust is D. Platt's unpublished computation of sums over zeta zeros up to height ten billion; it affects 11.39 but not our record. The paper is still an unrefereed v1 and should be cited as such.

## 6. Files
`chk_a.py` (barbar, artanor constants, Stirling/zetinvbound), `chk_artanor.py` (artanor end to end; monthmayending against the series), `chk_b.py` (nopgik/gopnik/bataclar/mertosimp numbers; proof-as-computation of gould and qmediano), `chk_b2.py` (exact μ-sieve to 10⁷: difintsq, oldtrivial, MOST spot checks, ramando, per-cell sugr/trapeze, andalas end to end), `chk_b3.py` (edesmo on the real Q), `chk_dep.py` (fallbacks and margins); outputs `logs/*.log`. Downloads (e-prints of 2511.14736v1, 1309.1527v1, 1912.04972v2) were deleted after use.

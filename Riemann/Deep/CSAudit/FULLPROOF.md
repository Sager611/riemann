# FULLPROOF — referee notes on Chourasiya–Simonič, arXiv 2507.15184 v2 (TeX InghamPostArXiv.tex, 1920 lines)

**Final verdict (follow-up, 12 Nov 2026): AUDIT PASSED for the corrected statement (b₄ = 0.5501, REPORT §D) — every step checked (§1–8); as printed Cor. 1.1 rests on the unsupported b₄ = 0.54, and Lemma 5.2's second bound is false for σ₀ > 2 (unused there).**

Line numbers = TeX lines. "✓" = step re-derived here; "✗" = error; "?" = gap.

## 1. Cor. 1.1 from Thm 2.1 (l. 228–275) — rechecked
- Thm 2.1 (l. 228–247) with T₀ = H₀ = H = H_RH, h = 1: (log hT)^{1/(2−σ)}(log T)^{(6−5σ)/(2−σ)} = (log T)^{(7−5σ)/(2−σ)};
  𝔄₂ term carries (log T)^{(1+σ)/(2−σ)} = (log T)^{(7−5σ)/(2−σ)}·(log T)^{−6(1−σ)/(2−σ)} ≤ ·(log H_RH)^{−6(1−σ₂)/(2−σ₂)} (exponent increasing in σ) ✓.
- 𝔄₃ log(HT) log T + 𝔄₄ log(hT) log T = (𝔄₃+𝔄₄) log²T + 𝔄₃ log H_RH log T ✓. N(σ, H_RH) = 0 for σ > 1/2 (RH verified to 3·10¹², Platt–Trudgian) ✓.
- σ₁ > 1/2 + d/log T: rows n ≥ 17 need d < (1/32) log(3·10¹²) = 0.8978; all d ≤ 0.89 ✓. Row [1/2, 0.53125]: Thm applied with σ₁ := σ, monotonicity of 𝔄ₙ in σ₁ (each max{K^{e(σ₁)},K^{e(σ₂)}} is non-increasing in σ₁; (3−2σ₁)/(2d) decreasing) ✓; K, L do not depend on σ₁ ✓.

## 2. §3, Gabriel (l. 433–611)
- Lemma 3.1 (Gabriel, two exponents): cited to Karatsuba–Voronin pp. 381–384; the standard form (Ivić Lemma 8.3) ✓ (not re-proved).
- Thm 3.2 proof (l. 511–545), every step re-derived:
  - |(z−1)/z| ≤ 1 ⟺ x ≥ 1/2 ✓; |e^{κ(z/T)²}| = e^{κ(x²−y²)/T²} ✓ ⇒ |g| ≤ ω₁(x)e^{−κ(y/T)²} ✓.
  - lower bound: ((x−1)²+y²)/(x²+y²) increasing in y² for x ≥ 1/2, so ≥ value at |y| = H; e^{−κy²/T²} ≥ e^{−κ} on |y| ≤ T ⇒ ω₂ ✓.
  - substitution y = T√(c/κ)v, d/dv I_f = T√(c/κ)|f|^{1/c} (both signs of v) ✓; integration by parts, boundary terms vanish by condition (3) ✓.
  - v<0 half: ∫₀^∞e^{−u}I_f(x,−T√(cu/κ),A,c)du; splitting I_f(−B,A) = I_f(−B,−A) + I_f(−A,A); on u ≤ (κ/c)(A/T)² both remaining integrals are ≤ 0, dropped ✓ ⇒ (gabrielU3) ✓.
  - F(x) = α²(β−x)+β²(x−α)−x²(β−α): F(α)=F(β)=0, concave, F((α+β)/2) = (β−α)³/4 ✓ ⇒ ω₁(α)^{w₁}ω₁(β)^{w₂}/ω₂(x) ≤ e^{κ(1+((β−α)/2T)²)}·(ratio)^{−m/2} ✓.
  - growth conditions (2),(3) for f_X = ζM_X − 1 (polynomial growth in |t|) ✓; φ = gf holomorphic (zero of order m at 1, g's pole at 0 outside the strip since α ≥ 1/2) ✓.
- Cor 3.3 (l. 558–577): m = 1; (2σ−1)/(σ²+H²) increasing on [1/2,3/2] for H > 1 (numerator of derivative −2σ²+2σ+2H² ≥ 2H²−3/2 > 0) ⇒ ratio ≥ 1−8/(9+4H²) ✓;
  1−8/(9+4H²) ≥ (1−1/H)² ⟺ 8H³−12H²+18H−9 ≥ 0 (=5 at H=1, derivative 24H²−24H+18 > 0) ✓; the factors 2^{c} from conjugate symmetry cancel since c₁ = w₁λ+w₂μ ✓.
  - remark: if M_X(1) = 0 the pole order is 0 < m = 1; the proof only needs φ holomorphic, so harmless.
- Cor 3.4 (convexity) is not used by Thm 2.1; not audited.

## 3. §4 inputs: Littlewood, Lemmas 4.1–4.4 (l. 612–721) against KLN (arXiv 2101.12263 = JMAA 465 (2018), TeX l. given as KLN l.)
- (eq:Littlewood) l. 619–624: Littlewood's lemma for h_X = 1 − f_X² on (σ′,σ″)×(H,T); ν(u) = #zeros of h_X with β > u is non-increasing, ≥ #zeros of ζ with β ≥ σ for u < σ ⇒ ∫_{σ′}^{σ″}ν ≥ (σ−σ′)(N(σ,T)−N(σ,H)) ✓. arg on Re = 3/2 principal (Re h_X > 0 there since |f_X| ≤ b₆ < 1) ✓.
- Lemma 4.1 (l. 627–642) = KLN Lemma 4.1 (KLN l. 636–673) with kH₀ = X₀ ≥ 10⁹: MV (KLN Lemma 3.6, real u_n = μ(n)n^{−1/2}) + Σμ² ≤ 0.62X (X ≥ 1700), Σμ²/n − (6/π²)log X ≤ 1.048 (X ≥ 1002) (Ramaré, Math. Comp. 2016, Lemmas 3.1, 3.4); each of C₁, C₂ decreasing in X ✓.
  (KLN's proof swaps the two labels; the algebra is right.)
- Lemma 4.2 (l. 644–670): MV ⇒ F_X(σ,0,T,1/2) ≤ πm₀Σλ²n^{1−2σ} + (T+πm₀)Σλ²n^{−2σ} ✓; λ_X(n) = 0 for n ≤ X ✓; split X<n<5X (KLN (3.5), b₃ = 0.605, X ≥ 10⁹, Ramaré Lemma 5.6) and n ≥ 5X (n^{2σ} ≥ (5X)^σ n^σ) ✓;
  Σ_{n≥1}(Σ_{d|n,d≤X}μ(d))²n^{−τ} = ζ(τ)S_{τ−1}(X) ✓ (S_ε of Ramaré–Zuniga Alterman (RZA), arXiv 2603.25961 v3 Thm 1.2, TeX l. 610–628: S_ε(X) = Σ_{d,e≤X}μ(d)μ(e)/[d,e]^{1+ε}).
  σ = 1+δ/log X: e^{2γ(σ−1)} ≤ e^{2γδ/log X₀}; X^{−2δ/log X} = e^{−2δ}; (5X)^{δ/log X} = e^δ5^{δ/log X}; (e^γ/5)^{δ/log X} ≤ 1 ⇒ C₃, C₄ ✓. Range: τ ∈ {σ, 2σ−1} ⊂ (1, 1.1] iff δ ≤ (log X)/20 ✓.
  **Only gap: b₄ = 0.54 (footnote l. 659) — see §6.**
- Lemma 4.3 (l. 672–700): |f_X(3/2+it)| ≤ Σ_{n>X}d(n)n^{−3/2} ≤ b₆(X,1/2) (KLN (4.45), Lemma 3.5 with τ = 3/2, valid for every η > 0, X ≥ 1) ✓;
  b₆(X,1/2) = 3X^{−1/2}(log X + 2 + γ + 7/(36X)) decreasing for log X + 2 + γ > 2 ✓ ⇒ b₆(X,1/2) ≤ b₆(hT₀,1/2) (= 0.0022 at 10⁹) ✓;
  −log|h_X| ≤ −log(1−|f_X|²) ≤ (−log(1−b₆²)/b₆²)|f_X|² (chord of the convex −log(1−y)) ✓;
  MV on [H,T] (KLN (3.10), constant 2πm₀): ∫_H^T|f_X(3/2+it)|² ≤ (T−H+2πm₀)Σ_{n>X}d²n^{−3} + 2πm₀Σ_{n>X}d²n^{−2};
  KLN (3.7): Σ_{n≥X}d²n^{−τ} ≤ 2τ(τ−1)^{−1}b₁₁(X,τ)log³X/X^{τ−1} (needs Σ_{n≤t}d² ≤ t(log t+1)³: true since d² ≤ d₄ and Σ_{n≤t}d₄ ≤ t(log t+1)³; ≤ 2t log³t for t ≥ 47) ✓;
  T−H+2πm₀ ≤ T (H ≥ 2πm₀), T/X² = 1/(hX), b₁₁(X,3) ≤ b₁₁(X,2) ≤ 1+3.31/log X for X ≥ 10⁹ (6/L+6/L² = 0.3035 ≤ 0.31 at L = log 10⁹) ✓; log³X/X decreasing ✓ ⇒ (eq:C5) ✓.
  "Cor. 4.11 valid for η = 1/2 by continuation": not even needed here (only b₆ is used); for Lemma 4.4 see next item.
- Lemma 4.4 (l. 702–720) from KLN Cor 4.11 (KLN l. 1262–1391) + (4.53) (KLN l. 1417–1438) at η = 1/2:
  KLN Cor 4.11 is stated for η ∈ [0.23622, 1/2); every step holds verbatim at η = 1/2: Rademacher (KLN Lemma 3.1) needs −1/2 ≤ −η ✓; Jensen circle centre 3/2, radius 2 reaches Re = −1/2 ✓; |Im w| ≥ U−2 ≥ 1000 ✓; |Re(w+1)| ≤ 2+η ✓ ⇒ direct proof, no continuation needed ✓.
  Recomputed: |arg h(τ+iT)|+|arg h(τ+iH)| ≤ b₉ + (2/log2)log(HT) + (3π/log2)log X with the same X = hT on both edges ✓;
  b₉(1/2,H) = (π/log2)log b₇ + (π/log2)log b₅ + π − (2π/log2)log(1−b₆²) + (4/log2)log(b₈/2π);
  (π/log2)·2log(3.006ζ(3/2)) = 18.68, (π/log2)log b₅(1/2) = (π/log2)log 43.59 = 17.11, π ⇒ 38.93 ≤ 220 (CS's constant is conservative by 181) ✓;
  log(1+2/(3.006ζ(3/2)u^{3/2})) ≤ log(1+u^{−3/2}) ✓; b₈(1/2,H) = √((5/2H)²+(1+2/H)²) decreasing in H ✓; integration length 3/2−σ+d/log T ≤ 3/2−σ+d/log T₀ ✓.
  KLN Prop 4.10 is stated for σ ∈ (0,1+η] but Jensen's n(u) ≥ n for u ≥ 1/2+η needs σ ≥ 1/2; CS use τ ≥ σ−d/log T ≥ 1/2 only ✓ (KLN misstatement, harmless here).

## 4. Lemma 4.5 (l. 722–792) and proof of Thm 2.1 (l. 841–859) — every exponent manipulation re-derived
- Cor 3.3 with A = A₀, α = 1/2, β = 1+ε (ε := δ/log X ≤ 1/20 ⇒ β ≤ 3/2 ✓), λ = 3/4, μ = 1/2; u-range starts at (4κ/3)(A₀/T)² ⇔ V := T√(3u/4κ) ≥ A₀ ✓.
- (eq:MinHol): |f_X|^{4/3} ≤ (|ζM_X|+1)^{4/3} ≤ 2^{1/3}(|ζM_X|^{4/3}+1) (convexity of x^{4/3}; "Minkowski" in the TeX, harmless), Hölder (3, 3/2) ⇒ ∫|ζM|^{4/3} ≤ 𝓜₂^{1/3}(∫|M|²)^{2/3} ✓.
- 𝓜₂(A₀,V) ≤ a₂(A₀)V log⁴V (Cor 2.4 needs V ≥ A₀ ≥ 10¹⁰ ✓); log V ≤ log T·(1+log(1+3u/4κ)/(2log T₀)) ✓; Lemma 4.1 at X = hT ≥ X₀ = hT₀ ✓ ⇒ first term of K₁ ✓ ((2a₂√(3/4κ))^{1/3} = (2a₂)^{1/3}(3/4κ)^{1/6}).
- F_X(1/2,0,A₀,3/4) ≤ (2a₁)^{1/3}(hC₂+A₀C₁/T₀)^{2/3}T^{2/3}log^{2/3}(hT) + 2^{1/3}A₀; T^{−1/3}log^{−4/3}T decreasing ⇒ second term of K₁ ✓.
- linear pieces: 2^{1/3}[A₀ + ∫_{u₀}^∞e^{−u}(c√u − A₀)du] (c = T√(3/4κ), A₀ = c√u₀) = 2^{1/3}c[Γ(3/2) + √u₀(1−e^{−u₀}) − ∫₀^{u₀}e^{−u}√u du] ≤ 2^{1/3}c[Γ(3/2) + u₀^{3/2}];
  c·u₀^{3/2} = (4κ/3)A₀³/T² ≤ (4κ/3)(A₀/T₀)³T ⇒ K₂ ✓ (the naive 2^{1/3}A₀ ≤ K₂T would be false; CS's finer route is correct).
- K₃: Lemma 4.2 at V′ = T√(u/2κ) (μ/κ = 1/(2κ)) with X ≥ hT₀; C₃, C₄ decreasing in X₀ ✓; πm₀/X ≤ πm₀/(hT₀) ✓.
- exponents, w = (3/4)w₁+(1/2)w₂ = (2−σ+3ε)/(2+4ε) ✓; (3/4)w₁/w = 3(1−σ+ε)/(2−σ+3ε) = 3(1−σ)/(2−σ) + 3ε(2σ−1)/((2−σ)(2−σ+3ε)) ∈ [e₀, e₀+3ε] ✓ (≤ 1, so (a+b)^x ≤ a^x+b^x ✓);
  (1/2)w₂/w = (2σ−1)/(2−σ+3ε) ∈ [(2σ−1)/(2−σ) − 3ε, (2σ−1)/(2−σ)] ✓; w₁/w = 4(1−σ+ε)/(2−σ+3ε) ≤ 4(1−σ)/(2−σ)+4ε ✓; 1/w ≤ 2/(2−σ) ⟺ σ ≥ 1/2 ✓.
  log T/log(hT) ≤ max{1, log T₀/log(hT₀)}; log log T/log(hT) ≤ max{loglogT₀/log(hT₀), loglogT₀/logT₀} for T ≥ T₀ ≥ e^e (d/dL of log L/(L+c) has numerator 1+c/L−log L < 0) ✓;
  (β−α)/(2T) = (1+2δ/log hT)/(4T) ≤ (1+2δ/log hT₀)/(4T₀); 1/(1−1/H) ≤ 1/(1−1/H₀); L₀ ≥ 1 ⇒ (·)^{1/w} ≤ L₀^{2/(2−σ)} ✓. ⇒ L₀, L₁, L₂ exactly as (eq:L0),(eq:L1),(eq:L2) ✓.
- Thm 2.1: σ′ = σ − d/log T, σ″ = 3/2; 1/(2π(σ−σ′)) = log T/(2πd) ✓. log|h_X| ≤ log(1+|f_X|²) ≤ |f_X|^{1/w}: 1/w(σ′) ∈ [4/3, 2), log(1+u²) ≤ u^a on a ∈ [4/3,2] (u ≤ 1: ≤ u²; u ≥ 1: ≤ u^{4/3}) ✓.
  σ′ → σ: exponent derivatives 3/(2−σ)², −3/(2−σ)², −4/(2−σ)², (2−σ) ≥ 1 ⇒ gaps ≤ 3d/log T, 3d/log T, 4d/log T ⇒ max{1,K₃^{−3d/logT₀}}, max{1,K₁^{3d/logT₀}}, T^{3d/log T} = e^{3d}, (log T)^{4d/log T} ≤ e^{4d·loglogT₀/logT₀} ✓;
  L₀^{2/(2−σ′)} ≤ L₀^{2/(2−σ)} (L₀ ≥ 1; CS's extra max{1,L₀^{−2d/logT₀}} = 1) ✓; endpoint maxima over σ ∈ [σ₁,σ₂] ✓; log T·(log T)^{4(1−σ)/(2−σ)} = (log T)^{(6−5σ)/(2−σ)} ✓.
  Lemma 4.3 ⇒ C₅/(2πd) ✓; Lemma 4.4 × log T/(2πd), (3/2−σ+d/logT₀)/d ≤ (3−2σ₁)/(2d)+1/logT₀ ⇒ 𝔄₃ (coefficient 2/(2π log2)), 𝔄₄ (3π/(2π log2)), C₆/(2π) ✓.
  Zeros of ζ are zeros of h_X with ≥ multiplicity ✓; pole of h_X at 1 outside the rectangle (H ≥ 1002) ✓; boundary zeros: perturbation ✓.
**§§3–4 verdict: Thm 2.1 follows from Lemmas 4.1–4.4, Cor 2.4 and b₄; no error found. The only unsupported input in §4 is b₄ = 0.54 (§6).**


## 5. Table 1 recomputed independently (table1_repro.py, logs/table1_repro.log; mpmath 30 digits, K₁ integral by two quadratures agreeing to 1e−20)
- All 16 rows from Thm 2.1/Lemma 4.5/Cor 2.4 at the printed (d, κ, δ, log₁₀A₀), b₄ = 0.54: ours/table for B₁ = 0.9991–1.0000; B₂, B₃ agree to the printed digits (ours ≤ table: table rounds up) ✓. No row fails.
- b₄ = 0.5501 (RZA Thm 1.2, 10¹² ≤ X < 10¹⁶, ε ≤ 1/50; X ≥ 10¹⁶: 0.5103) — B₁ rises by 0.08% (row 1) … 1.2% (row 12: 21.77 → 22.03) … 1.4% (row 13: 24.83 → 25.17); B₂, B₃ unchanged.
- ε = 2δ/log X ≤ 2δ/log(3·10¹²) ≤ 0.0150 ≤ 1/50 for rows 1–13 ✓ (rows 14–15 ✓; row 16: 0.0211 > 1/50 for X < 1.45·10¹³, needs RZA's 0.5936 (ε ≤ 1/25) there: B₁ 46.06 → 50.60).
- Conversion ζ(τ)S_{τ−1}(X) ≤ b₄e^{γ(τ−1)}/(τ−1) needs ζ(1+ε) ≤ e^{γε}/ε; rigorous simple substitute ζ(1+ε) ≤ 1+1/ε gives b₄″ = 0.5501(1+(1−γ)ε) (cert.py's CSRow): row 12 → 22.11, row 13 → 25.29.


## 6. §5.2 divisor-sum lemmas: falsification (divisor_falsify.py; logs/divisor_falsify_*.log)
- Grid: x₀ ∈ {2,…,10⁵}, σ (and σ₀) over each stated range, all integer x ≤ N (both one-sided limits for tails, staying ≥ x₀); tails for σ > 1 from exact ζ(σ)², ζ(σ)⁴/ζ(2σ) (mpmath) + reverse cumsums.
- **✗ Lemma 5.2, second bound (l. 963, D₂ l. 968) is false as stated for σ₀ > 2**: x₀ = 5, x = 5.99, σ = σ₀ = 4: Σ_{n>x}d(n)n^{−4} = Σ_{n≥6} = 0.0068135 > D₂(5,4)/9·x^{−3}log x = 0.0058938 (ratio 1.156; x → 6⁻: 1.161). (At x = x₀ = 6 it holds: 0.003727 < 0.005511 — reviewer correction.)
  Cause: partial summation with |Δ(u)| ≤ √u gives Σ_{n>x}d(n)n^{−σ} ≤ x^{1−σ}log x·[(σ−1) + (1+2γ(σ−1))/log x + (σ−1)²(4σ−1)/((2σ−1)√x log x)]/(σ−1)²; D₂ drops the factor (σ₀−1)² in its last term.
  Harmless: D₂ is used only at σ₀ = 1−σ+|c₁| ∈ (1, 2) (𝔤₅ l. 1510: σ ∈ [1/2,σ₄], σ₄ < |c₁|; 𝔤₇ l. 1533: σ ∈ [σ₅,1/2], σ₅ > 1/2−c₂; 𝒥₀₁ l. 1622: 1/2+|c₁| < 3/2), where (σ₀−1)² < 1 makes the printed D₂ larger than the correct one ✓.
- Lemma 5.5b (σ > 1 tail of d²): ratio → 1⁻ as σ → 1⁺ at x = x₀ (max 1.000000 at σ = 1.001, x₀ = 10⁴): the proof keeps the boundary term exactly (identity (eq:d2Aux1)); slack = (P⁺−P) + (9.73u^{3/4}log u − D₀) ≥ 0, so true given Cully-Hugill–Trudgian ✓ (coefficient check: (6/π²)(4γ−2ζ′(2)/ζ(2)) = 2.096609 vs 6D₁+2D₂⁺ = 2.096611).
- Lemmas 5.2a, 5.3 (both), 5.4 (|D₀| ≤ 9.73x^{3/4}log x: max ratio 0.19 for x ≤ N), 5.5a, 5.6, 5.7, 5.8: max LHS/RHS < 1 on the whole grid ✓ (values in the log).


## 7. §5.3 (proof of Thm 2.2) — what was checked
- Lemma 5.1 (MV, l. 871–885): Preissmann's Hilbert inequality twice, δ_n = log(1+1/n), 1/log(1+1/n) ≤ n+1/2 ✓.
- Lemma 5.9 (χ, Γ bounds, l. 1252–1273) and Lemma 5.10 (residue, l. 1275–1295): numerical falsification on grids (chi_falsify.py, logs/chi_falsify.log): all max ratios < 1
  ((eq:Gamma) 0.933, (eq:chiBound) 0.933, |χ| ≤ 0.4|y|^{1/2−x}+16: 0.9934 at x = 0, y = 10⁸ (→ 0.3989/0.4), |(z−1)χ| ≤ 6.31: 0.9989 at (1/2, 2π), |χ| ≤ 1.1 (|y| ≥ 2π): 0.909, Lemma 5.10: 0.603) ✓; Lemma 5.10's proof re-derived (|Γ(1/2+it)| ≤ √(2π)e^{−πt/2}) ✓.
- I. part (eq:J1)/(eq:cJ1): MV with a_n = d(n)n^{−1/2}, Σ(n+1/2)d²/n = Σd² + ½Σd²/n, Lemma 5.6 two-sided, (m₀/2π)log⁴T = (m₀log T/2πT)·T log³T ⇒ 𝒥₁ = 2ηD₆ + 4πm₀D₄(T₀,0) + m₀logT₀/(2πT₀) ✓.
- II. part (eq:secterm)–(eq:cJ21): J₂(t) = f(1/2+it), J₁(t) = f(1/2−it) ✓; contour shift to Re z = σ₁ ∈ [0,1/2) (f holomorphic for |Im z| > 1) ✓; |χ(z)|^{−1} ≤ (T/2π)^{x−1/2}/U⁻ for x < 1/2, |z| ≥ T ✓;
  horizontal pieces via Lemma 5.2a (σ₀ = 0) ✓; vertical piece via Lemma 5.1 + Lemma 5.5a at 2σ₁ ∈ [0,1), 2σ₁−1 ∈ [−1,0) ✓ ⇒ 𝒥₂,₁ ✓.
- Inputs of Cor 2.3 (l. 371): M₂(3000) = 703127.63 ∈ (7.031·10⁵, 7.032·10⁵) ✓, M₂(3000,6000) = 1249636.00 ∈ (1.249·10⁶, 1.25·10⁶) ✓ (m2_check.py: Arb ζ, Gauss–Legendre 6 and 10 nodes per 0.5 agree to 3·10⁻⁴; logs/m2_check.log).
- (Superseded by §8: parts III, IV (𝒥₅, 𝒥₆), VI and the Table 2/3 numerics are now checked.)


## 8. §5.3 III and VI parts, and the Mathematica outputs (follow-up, 21:07–21:30 UTC)
- (eq:zeta2) re-derived: shift Re z: 2 → c₁ crosses z = 0 (gives ζ²) and z = 1/2−it (Res); FE on Re z = c₁ (Re(1/2−it−z) = 1/2−c₁ > 1, Dirichlet series absolute); Σ_{n≤T} part shifted c₁ → c₂ < 1/2 crosses only z = 0 (gives J₁); ×χ(1/2−it) ⇒ |ζ|² = ΣJ_k ✓. Grouping (A+B+C)² ⇒ parts I–VI ✓.
- III, k ∈ {3,4} (l. 1346–1415): J₃ = g₁(1/2+it), J₄ = g₂(1/2+it), J₁ = f(1−z) ✓; f(1−z) = χ(z)Σ_{n≤T}d(n)n^{z−1} (χ(s)χ(1−s) = 1), so F₁ has no χ at all — CS's factors (2+1/T₀)^{σ₂−1/2}/U⁻² etc. are ≥ 1, conservative ✓;
  𝔊₁ (eq:fG1): CBS + Lemma 5.1 + Lemmas 5.5a/5.7/5.8 at 2−2σ₂ ∈ [0,1), 1−2σ₂ ∈ [−1,0), 2σ₂−2, 2σ₂−3, 2σ₃, 2σ₃−1, 2σ₃−2, 2σ₃−3 ∈ [−3,1) ✓, T-powers sum to T log³T ✓; |χ^{−1}(σ₃+iu)|² ≤ (2π/T)^{1−2σ₃}/U⁻² ✓;
  𝔊₂ (eq:fG2): Lemma 5.2a (σ₀ = 0, exponent 1−u ∈ [0,1/2]) × Lemma 5.3 (σ = u ∈ [0,1]) ✓; two horizontal segments ⇒ the two (≥ 1) factors ✓. Typo: (eq:F1F2main) writes lower limits 1+iT, iT for σ₂+iT, σ₃+iT.
- III, k ∈ {5,6} (l. 1416–1565): g₃(1/2+it) = J₅, g₄(1/2+it) = J₆ ✓; holomorphic for 1/2−c₂ < Re z < |c₁| (Σ_{n>T} absolute iff Re z < |c₁|), Im z > 1+B (condition (6)) ✓; |χ(z)| ≤ (2π/T)^{σ₄−1/2}/U⁻ via the lower bound at 1−z ✓;
  𝔤₁: |χ^{−1}(z)|² ≤ T^{2σ₄−1}(1/π+1/(2πT₀))^{2σ₄−1}/U⁻²; |χ(w)|⁴ (Re w = σ₄+c₁ ∈ (−1/2,0), |w| ≤ 1/2+2T+A, |Im w| ≥ T−A ≥ T₀−A₀ by (5)) ✓; weighted CBS in u with weight |Γ(c₁+iu)| ✓; Lemma 5.5b at 2−2σ₄+2|c₁|, 1−2σ₄+2|c₁| > 1 ✓; T-power 1 ✓;
  𝔤₂: Re w = σ₄+c₂ ∈ (1/2,1) for every printed row (max 0.823; the theorem's general range allows Re w up to 3/2, outside Lemma 5.9's x ∈ [−1,1] — gap in the general statement, not in the numbers), |w| ≥ T−B ✓; Lemma 5.5a at 2(1−σ₄−c₂) ∈ (0,1), 1−2(σ₄+c₂) ∈ (−1,0) ✓;
  𝔤₃, 𝔤₄ (σ₅-line, Re w ∈ (−1,0) resp. (1/2,1)) ✓; 𝔊₃, 𝔊₄ ✓; horizontal 𝔤₅–𝔤₈, 𝔊₅, 𝔊₆: χ(z)χ^{−1}(z) = 1 (conservative factors), Lemma 5.2b at σ₀ = 1−σ+|c₁| ∈ (1,2) ✓, Lemma 5.2a with σ₀ = −1/2 (𝔤₆) or 0 ✓, T-power 1 ✓.
  (eq:J22) states Σ_k|∫(J₁+J₂)J_k| but proves |∫(J₁+J₂)(J₃+J₄)| + |∫(J₁+J₂)(J₅+J₆)|; both are bounded by 𝒥₂,₂ (the proof splits g₁, g₂, g₃, g₄ by the triangle inequality), and only the grouped form is needed for |R| ✓.
- VI (l. 1603–1650): |χ(1/2+c₁+i(t+u))| ≤ 0.4(2T+|u|)^{|c₁|}+16 ≤ (2T)^{|c₁|}(0.4+16/(2T₀+A₀)^{|c₁|})(1+|u|/2T₀)^{|c₁|} (|u| ≥ A(T) ≥ A₀) ✓; |Γ(c₁+iu)| ≤ √(2π)e^{1/(2A₀²)}|u|^{−|c₁|−1/2}e^{−|u|/2}e^{−A(T)} ((π−1)/2 > 1) ✓;
  Lemma 5.2b at σ₀ = 1/2+|c₁| < 3/2 ✓; ∫_{|u|≥A} ≤ 2∫₄^∞ (A ≥ log T₀ ≥ 4) ⇒ 𝒥₀₁ ✓; 𝒥₀₂ ✓ (typo: χ(1/2+c₁+…) for χ(1/2+c₂+…); |χ| ≤ 6.31/(1/2−c₂), ∫_{|u|≥B}e^{−|u|/2} = 4e^{−B/2} ≤ 4e^{−4});
  𝒥₀: |ζ| ≤ 0.618t^{1/6}log t (Hiary–Patel–Yang Thm 1.1), 2·2^{1/3}(log 2T/log T)² ≤ 3.5T^{1/6} (T ≥ 55), T^{1−a(T)} and T^{2−b(T)} decreasing by (5), (6), (x+y+z)²/log T ≤ (3/4)(x²+y²+z²) (log T₀ ≥ 4) ⇒ 7/2 and 3/14 ✓.
- **Mathematica outputs recomputed** from the formulas above (f1_repro.py, m2rec_repro.py; logs/f1_repro.log, logs/m2rec_repro.log):
  𝔉₁ at the printed 11 parameters for T₀″ = 10⁴…10⁸, 10⁹, 10¹⁰, 10¹⁵, 10²⁰: 59.285, 48.800, 43.228, 39.714, 37.544, 36.211, **34.841**, 31.166, 29.364 ≤ printed 𝔐₁ (ratios 0.99975–0.999997) ✓; conditions (2)–(6) hold at every printed row ✓;
  𝒞₁(3000) = 20.72246 ≤ 20.7225 ✓; recursion (l. 418–427): M₂(10ⁿ,2·10ⁿ), M₂(10ⁿ), 𝔐₂(10ⁿ) reproduced for all n, e.g. **M₂(10¹⁰,2·10¹⁰) ≤ 9.50787·10¹⁵** ≤ 9.5079·10¹⁵, 𝔐₂(10¹⁰) ≥ 1.38962·10¹⁶ ≤ 1.3897·10¹⁶ ✓.
- Sources: 𝔐₁(10¹⁰) = 34.842 and M₂(10¹⁰,2·10¹⁰) ≤ 9.5079·10¹⁵ are the authors' own outputs (Table 2, l. 186–194; from Thm 2.2 at Table 3 and the recursion), no external source; no published explicit fourth-moment bound of the form a·T log⁴T exists to compare
  (max|ζ|²·M₁ via Hiary–Patel–Yang gives ≈ 1.3·10¹⁷ for M₂(10¹⁰,2·10¹⁰), 14× weaker and of the wrong shape for Cor 2.4). Cor 2.4 uses both at T₀″ = 10¹⁰ ≤ A₀, T ≥ A₀: inside their range ✓.


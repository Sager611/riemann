# WangAudit — FULLPROOF (chain of Wang arXiv 2609.07918 v2, TeX zeros_short_intervals.tex, 1441 lines; line numbers = TeX lines)

Setup (l.183–192): fixed 0 < λ < θ < 1, H = T^θ, L = log T, I = (T, T+H]; ρ = 1/2 + δ + iγ, |δ| < 1/2; ĝ(z) = ∫g(u)e^{−2πizu}du.
c(θ) = 2 − θ/2 − cot(θ/√2)/√2 (l.122); c′ = cot²(θ/√2)/2 > 0 (l.126) ✓ (differentiated by hand).

## 1. Objects (l.200–266)
- prop:finite (l.203–216) = Lamzouri 2609.02882 Prop 2.1 (checked R481).
- (eq:weighted-sum) l.218: w(z) = 4/(4 − z²), F_I(x) = Σ_{γ,γ′∈I} x^{ρ−ρ′} w(ρ−ρ′); (eq:even-F) F_I(1/x) = F_I(x): swap ρ↔ρ′, w even ✓.
- thm:local-pair (l.229–240): LHS = Σ ĝ(i(ρ−ρ′)L/2π) w(ρ−ρ′) = ∫ g(u) F_I(T^u) du (ĝ(iz L/2π) = ∫g(u)T^{uz}du) ✓.
  Main-term shape: F_I(T^α) ≈ (HL/2π)(|α| + L T^{−2|α|}); ∫ g(α) L T^{−2|α|} dα → g(0) ✓ (consistent with Montgomery).
- (eq:A) l.246: A(x,t) = Σ_ρ 2x^{δ+i(γ−t)}/(1 + ((t−γ)+iδ)²); A_I = restriction to γ ∈ I (finite sum).

## 2. lem_norm (l.268–302): 2πF_I(x) = ∫_R |A_I(x,t)|² dt — CORRECT
- |A_I|² = Σ_{ρ,ρ′} 4x^{d}/[(1+(t−γ+iδ)²)(1+(t−γ′−iδ′)²)], d = ρ + ρ̄′ − 1 = δ+δ′ + i(γ−γ′) ✓.
- u = t−γ+iδ, a = γ−γ′−i(δ+δ′) = −id ✓; ∫_R du/((1+u²)(1+(u+a)²)) = 2π/(4+a²), |Im a|<1: residues at i, −a+i ✓ (by hand).
- Contour R → R+iδ: poles Im u = ±1 and δ+δ′ ± 1; none between 0 and δ since |δ|,|δ′| < 1/2 ✓ (Wang does not state this shift; harmless).
- 4·2π/(4−d²) = 2πw(d) ✓; ρ′ ↦ 1−ρ̄′ is a multiplicity-preserving bijection of zeros fixing γ′ (FE + conjugation) ✓ ⇒ d = ρ−ρ″.
- Consequence F_I(x) ≥ 0 real ✓. No RH: only |δ| < 1/2 used.

## 3. lem:localization (l.305–437): 2πF_I(x) = ∫_I |A(x,t)|² dt + O(xL³), uniformly 1 ≤ x ≤ T^λ — CORRECT
- (number_of_zeros) l.314: N(t+1)−N(t) ≪ log(|t|+3) (Riemann–von Mangoldt) ✓.
- (eq:denominator) l.319: |1+(s+iδ)²| ≥ Re = 1+s²−δ² > s²+3/4 ✓ (only |δ|<1/2).
- (eq:A-size) l.325: |A|,|A_I| ≤ 2x^{1/2}Σ_ρ 1/|…| ≪ √x Σ 1/(1+(t−γ)²) ≪ √x L (x^δ ≤ x^{1/2} for x ≥ 1; Davenport ch.15 lemma) ✓.
- (eq:A-tail) l.330: unit-interval partition: γ ∈ (T−j−1, T−j] ⇒ t−γ ≥ a+j ✓; Σ_j 1/(1+(a+j)²) ≍ 1/(1+a) ✓;
  |γ| > 3T, t ≤ 2T ⇒ |t−γ| ≥ |γ|/3 ✓, tail ≪ L/T ✓; 1/(1+a)+1/(1+b) ≥ 2/(1+H) ≥ 1/T (convexity, a+b=H) ✓.
- (eq:A-exterior) l.417: |A_I(x,t)| ≪ √x L/(1+dist(t,I)) (I ⊂ [T,2T], O(L) zeros per unit interval) ✓.
- ∫_I ||A|²−|A_I|²| ≤ ∫_I |A−A_I|(|A|+|A_I|) ≪ xL²·2log(1+H) ≪ xL³ ✓; ∫_{R∖I}|A_I|² ≪ xL² ✓; with lem_norm ⇒ claim ✓.
- Uses only |δ| < 1/2, x ≥ 1 and H ≤ T; the constant is absolute; uniform in T ✓. No RH.

## 4. lem:coeff (l.439–466): Σa_n² = log x + O(1), Σ n a_n² ≪ x log²(2x), x ≥ 1 — CORRECT
- a_n² = nΛ²/x² (n ≤ x), x²Λ²/n³ (n > x) ✓. Partial summation with S(u) = Σ_{n≤u}Λ² = u log u + O(u):
  x^{−2}Σ_{n≤x} nΛ² = ½log x + O(1) ✓; x²Σ_{n>x}Λ²/n³ = −log x + 3x²∫_x^∞(log u/u³)du + O(1) = ½log x + O(1) ✓ (by hand).
- Second sum: x^{−2}log²x Σ_{n≤x}n² + x²Σ_{n>x}log²n/n² ≪ x log²(2x) ✓.

## 5. prop:pointwise (l.468–513): F_I(x) = (H/2π)(L²/x² + log x) + O(H + HL/x² + H√L/x + xL³ + L/√x), 1 ≤ x ≤ T^λ — CORRECT given (eq:explicit),(eq:E)
- (eq:D-mean): Montgomery–Vaughan Cor. 3 on an interval of length H (translate a_n ↦ a_n n^{−iT}); Σ n a_n² < ∞ ✓ ⇒ H log x + O(H + x log²2x) ✓.
- (eq:B-mean): log(t+2) = L + O(1) on I ⊂ (T, 2T] ✓. (eq:BD): IBP ∫_I log(t+2)n^{it}dt ≪ L/log n ✓; Σ a_n/log n ≪ √x ✓; D abs. conv. ✓.
- ‖D_x‖²_{L²(I)} ≪ H log x + H + T^λL² ≪ HL — USES λ < θ (T^λL² = o(H)) ✓. ‖E_x‖ ≪ √H/x from (eq:E) ✓.
- Cross terms by Cauchy–Schwarz: DE ≪ H√L/x, BE ≪ HL/x², EE ≪ H/x² ✓; + lem:localization O(xL³) ✓ (xL² absorbed).

## 6. Proof of thm:local-pair (l.516–640) — CORRECT
- ĝ(i(ρ−ρ′)L/2π) = ∫g(α)T^{α(ρ−ρ′)}dα ✓; finite sum ⇒ LHS = ∫_{−λ}^{λ} g F_I(T^α) dα = 2∫_0^λ g F_I(T^α) dα (eq:even-F, g even) ✓.
- M_0 = (H/2π)2L²∫_0^∞ g e^{−2Lα} = (HL/2π)g(0) + O(H‖g′‖_∞) (MVT on all α ≥ 0, g ∈ C¹) ✓.
- M_1 = (HL/π)∫_0^λ αg = (HL/2π)∫|α|g ✓.
- Errors × 2‖g‖_∞: H, H/2, H/√L, L²(T^λ−1), 2 ✓ ⇒ O_g(H + T^λL²) with O_g = O(‖g‖_∞ + ‖g′‖_∞) and λ,θ-dependent constants only via (eq:E) and ‖D‖ ✓.
- Where λ < θ enters: only ‖D_x‖ ≪ √(HL) (and to make the error o(HL)); θ < 1 enters via H ≤ T (eq:A-tail) and x < T in (eq:E).
- Fixedness of g: only through ‖g‖_∞, ‖g′‖_∞ and supp g ⊂ [−λ, λ]; smoothness beyond C¹ is not used here.

## 7. Inputs from BGST (arXiv 2306.04799 v1 = Acta Arith. 214; TeX UnconditionalPC_230606.tex; arXiv eq. numbers match Wang's citations)
- (2.7) l.242: ∫ du/((1+u²)(1+(u+a)²)) = 2π/(4+a²), |Im a| < 1 ✓ (residues by hand, §2).
- (2.14) l.281: |A(x,t)| ≪ x^{1/2}Σ_γ 1/(1+(t−γ)²) ≪ x^{1/2}log(|t|+2) ✓.
- Lemma 1 (l.203–227, Montgomery 1973 without RH): for x ≥ 1, all t, A(x,t) = −Σ_n Λ(n)n^{−1/2−it}min(n/x,x/n) + x^{−1}(log(|t|+2)+O(1))
  + O(x^{1/2}/(1+t²)) + O(x^{−5/2}/(|t|+2)). Re-derived from Landau's formula R(s) at s = 3/2+it and s = −1/2+it (σ ↦ 1−σ):
  1/(σ−½+u) − 1/(−(σ−½)+u) = (2σ−1)/((σ−½)²+((t−γ)+iδ)²), u = i(t−γ)−δ ✓; xR(3/2+it): x(Σ_{n≤x} + ζ′/ζ) = −Σ_{n>x}Λ(n)n^{−1/2−it}(x/n) ✓;
  −x^{−1}Σ_{n≤x}Λ(n)n^{1/2−it} = −Σ_{n≤x}Λ(n)n^{−1/2−it}(n/x) ✓; pole terms x^{1/2−it}[1/(½+it) + 1/(3/2−it)] ≪ x^{1/2}/τ² ✓;
  trivial-zero series ≪ x^{−5/2}/τ ✓; −x^{−1}ζ′/ζ(−½+it) = x^{−1}(log τ + O(1)) (FE + digamma; ζ′/ζ(3/2−it) = O(1)) ✓.
  s = 3/2+it, −1/2+it never hit ρ, 1, 0, −2n ✓; x = p^m by continuity in x (both sides continuous) ✓. No RH, no zero-density input.
- Wang (eq:explicit)/(eq:E) l.251–266: Λ(1) = 0 so n ≥ 2 ✓; for t ∈ I ⊂ (T, 2T]: E_x ≪ 1/x + x^{1/2}/T² + x^{−5/2}/T ✓; ≪ 1/x needs x ≤ T^{4/3} and x ≥ 1 ✓.
- MV Cor. 3 (Montgomery–Vaughan 1974): ∫_{T}^{T+H}|Σ a_n n^{−it}|² = Σ|a_n|²(H + O(n)) for Σ n|a_n|² < ∞ ✓ (standard).

## 8. Removing the weight (l.643–783) — CORRECT
- η ∈ C_c^∞((−λ/2, λ/2)) real even, ∫η² = 1; f = η², K = f̂, Q = f∗f (supp ⊂ (−λ, λ)), Q̂ = K², hat{Q″} = −4π²z²K² ✓.
- (eq:remove-weight): Q̂ − hat{Q″}/(4L²) = K²(1 + π²z²/L²); z = i(ρ−ρ′)L/2π ⇒ 1 + π²z²/L² = 1 − (ρ−ρ′)²/4 = 1/w(ρ−ρ′) ✓ (w finite: |Re(ρ−ρ′)| < 1).
- S_K(I) = W_I(Q) − W_I(Q″)/(4L²) ✓; thm:local-pair for Q, Q″ (fixed, C_c^∞, supp ⊂ [−λ, λ]) ✓ ⇒ (S_K_I) ✓.
- M(Q) = Q(0) + ∫|α|Q = ∫f² + ∬|u−v|f(u)f(v) = C(f) ✓; M(Q″) = Q″(0) + 2Q(0) = −∫f′² + 2∫f² ✓ (IBP by hand) ⇒ O(H/L).
- (eq:unweighted): S_K(I) = (HL/2π)C(f) + O_f(H + T^λL²) = (C(f) + o(1))HL/2π — needs λ < θ ✓.
- (eq:scaled-zero-set): z = i(ρ−½)L/2π = (−γ + iδ)L/2π; real ⟺ δ = 0 ✓; z̄ ↔ 1 − ρ̄ (same γ, same multiplicity) ⇒ conjugation-invariant ✓;
  K(z−s) = K(i(ρ−ρ′)L/2π) ✓; multiset sums carry multiplicity ✓. (eq:both-counts) ✓; with N(I) = (HL/2π)(1 + O(1/L)): (eq:window-bound) ✓.
- Lamzouri's η needs only L², real, even, supp ⊂ (−λ, λ), ∫η² = 1: η ∈ C_c^∞((−λ/2, λ/2)) qualifies ✓.

## 9. prop:optimal (l.790–931) and thm:main (l.934–942) — CORRECT
- f_λ = cos(√2u)/(√2 sin(λ/√2)) on J = [−λ/2, λ/2]: ∫_J f_λ = 1 ✓; λ ≤ 1 ⇒ λ/√2 < π/2 ⇒ f_λ > 0 on closed J ✓ (positivity: f = η² needs f ≥ 0).
- G = f_λ + ∫|u−v|f_λ: G″ = f_λ″ + 2f_λ = 0, even ⇒ constant = f_λ(a) + a = cot(λ/√2)/√2 + λ/2 = C_λ (zero first moment) ✓; C(f_λ) = ∫f_λG = C_λ ✓.
- Minimality: h = f − f_λ, ∫h = 0; V_h′ = 2P, P(±a) = 0 ⇒ ∬|u−v|hh = −2∫P² ✓; Wirtinger/Dirichlet–Poincaré ∫P² ≤ (λ/π)²∫P′² (HLP Thm 257) ✓
  ⇒ C(f) − C_λ ≥ (1 − 2λ²/π²)‖h‖² > 0 for λ ≤ 1 (needs only λ < π/√2) ✓. Unique minimizer ✓.
- Admissibility (the lower bound needs only this, not minimality): η_ε = χ_ε√f_λ/√b_ε ∈ C_c^∞((−λ/2, λ/2)), real, even, ∫η_ε² = 1 ✓
  (√f_λ smooth since f_λ > 0 on J); f_ε → f_λ in L¹, L² (dominated convergence) ✓; |C(f_ε) − C_λ| ≤ |‖f_ε‖² − ‖f_λ‖²| + 2λ‖f_ε − f_λ‖_{L¹} → 0 ✓.
- thm:main: for fixed λ < θ, ε: (eq:window-bound) ⇒ liminf N^s/N ≥ 2 − C(f_ε); ε → 0: ≥ 2 − C_λ; C_λ′ = −cot²(λ/√2)/2 < 0, continuous
  ⇒ sup_{λ<θ}(2 − C_λ) = 2 − C_θ = c(θ) ✓ (liminf ≥ each, hence ≥ sup; no uniformity in λ needed). Holds for every θ ∈ (0,1) (trivial if c ≤ 0).
- So c(θ) is a PROVED lower bound (explicit admissible family), not only an optimum claim.

## 10. Versions and SL's use
- v1 (TeX zeros.tex, 7 Sep) vs v2 (6 Oct), §2–§4 token diff (vdiff2.py, ratio 0.984): wording only, plus v1's distinct-zero bound
  (eq:window-bound-distinct, N^d/N ≥ 3/2 − C(f)/2) moved to v2's new §5 (thm:distinct-all). Proofs of Prop 2.1, Thm 2.2 (= thm:local-pair in both),
  lem_norm … prop:optimal identical ⇒ SL's citation of v1 is to the same proof. v1 acknowledgment: AI-assisted drafting ("ChatGPT-6 Astra"), author-verified.
- SL v0.03 Thm 5.1 (sl.txt l.844–870): uses, for η ∈ C_c^∞((−λ/2, λ/2)) real even, ∫η² = 1, 0 < λ < θ: q_T = Q_T/N → C(f) ("[14, Theorem 2.2]";
  really Wang's §3 lemma (eq:unweighted) + N(I) ~ HL/2π), then f → cosine density (C(f) → C_λ = 2 − c(λ), prop:optimal's last part), then λ → θ.
  Q_T is Wang's S_K(I) (same scaled zeros i(ρ−½)L/2π, K = f̂). SL does NOT use Lamzouri's Prop 2.1, (eq:both-counts) or thm:distinct-all
  (its parity inequality is its own [4, Prop 2.5], audited in SantiAudit). Order of limits T → ∞, f → f_λ, λ → θ matches Wang's ✓.

## 11. Reviewer (one subagent, 19:08Z; Steps lem_norm, lem:localization, BGST Lemma 1, prop:pointwise pasted inline with exact formulas): CORRECT
- Contour shift crosses no poles (Im ±1, δ+δ′±1); reflection preserves window and multiplicity; signs/normalisation of Lemma 1 right; Landau's
  zero sum taken symmetrically; the two trivial-zero series must be combined before x ↓ 1: Σ_k x^{−2k−1/2−it}·2/((2k−½+it)(2k+3/2+it))
  ≪ x^{−5/2}/(1+t²), uniformly for x ≥ 1 (each series alone ~ log 1/(1−x^{−2}) diverges as x ↓ 1; checked by hand) ✓; cross term B·D̄ by
  termwise IBP (as Wang does), not Cauchy–Schwarz ✓. Uniform for fixed 0 < λ < θ < 1 as T → ∞.

## 12. Why no zero-free region is needed (contrast with BGST Lemma 4)
- BGST's full-window localization (their Lemma 4, (2.16)) bounds x^{δ+δ′} ≤ x^{2Θ−1} ≤ x^{1−2η(T log²T)} (Vinogradov-type zero-free region) because x runs up to T there.
- Wang's lem:localization uses only x^{δ+δ′} ≤ x (|δ| < ½) and pays O(xL³); integrated this is L²(T^λ − 1) = o(HL) exactly because λ < θ. So the
  short-window range α < θ is what makes the crude bound admissible: consistent, no hidden input ✓.
- Jump cancellation at x = n (prime power) in Lemma 1: x·Λ(n)n^{−3/2−it} = x^{−1}Λ(n)n^{1/2−it} at n = x ⇒ RHS continuous in x ✓.
- SL's Thm 5.1 class: the PDF layout puts the "∞" of C_c^∞ on the next line (sl.txt l.855–856), so SL writes C_c^∞ (SantiAudit's "minor: C_c" is an extraction artifact).

## 13. Numerics (falsification; all single-threaded under memcap 290 MB; logs in this directory)
| check | script → log | result |
|---|---|---|
| BGST (2.7), 12 random a, |Im a| < 1 | int27_recheck.py | rel. err ≤ 1.0e−31 ✓ (first run 5e−3 = quadrature near the pole −a ± i, fixed by breakpoints); Im a = 1.2: fails (rel 1.3), as it must |
| lem_norm with OFF-LINE synthetic zeros (δ up to 0.45, a double zero), x = 1, 3.7, 25 | lemnorm_test.py | 2πF_I = ∫|A_I|² to 1e−31, Im F_I ≈ 1e−31 ✓; WITHOUT the partners 1−ρ̄: off by up to 37 %, F_I complex ⇒ the FE reflection is genuinely used and suffices |
| BGST Lemma 1 exact form A = xR(3/2+it) − x^{−1}R(−1/2+it), 1250 zeros (mpmath.zetazero, γ ≤ 1704.29), x = 1.5, 6.5, 30.5, t = 300, 800 | bgst_lemma1.py | agreement 1e−6–1e−5 (omitted-zero bound 2–3e−3) ✓; x|E_x| = 1.87–2.01 for all x (≈ log 2π + ζ′/ζ(3/2−it)) ⇒ E_x ≪ 1/x ✓ |
| lem:localization, T = 600, θ = 0.95 (N(I) = 336), x = 1, 3, 10, 40 | localization_num.py | (2πF_I − ∫_I|A|²)/(xL³) = −0.73, 0.003, −0.0001, −0.0007: bounded ✓ |
| prop:pointwise main term | localization_num.py | (2π/H)F_I = 23.4, 3.65, 2.23, 3.43 vs L̃²/x² + log x = 20.8, 3.41, 2.51, 3.70 (L̃ = log(T/2π)); with L = log T: 40.9, 5.6, 2.7, 3.7 — rough, as expected (errors O(H(1 + L/x²)) vs H L²/x²) |
| (eq:unweighted) S_K/N vs C_λ, f = f_λ, λ = 0.5θ, 0.8θ, 0.95θ; (T, θ) = (1000, .9), (1200, .75), (600, .95) | localpair_num.py | with L̃: S_K/N = 1.408/1.451, 1.593/1.639, 1.356/1.405 at λ = 0.95θ (2–3 % low); with Wang's L = log T: 16–19 % low = log 2π/log T, the O(H) term: ROUGH AGREEMENT ONLY |
| F(α) = 2πF_I(T^α)/(HL̃) vs α + L̃T^{−2α} | localpair_num.py | T = 1000: α = 0.5: 0.47/0.53, 0.75: 0.69/0.75, 0.9: 0.85/0.90, 1.0: 0.99/1.00; beyond θ (α = 1.2, 1.5): 1.15, 1.01 ≈ 1 (GUE-like, not claimed) |
| c(θ) values; θ₀; Wang's examples | finite_checks.py | c(1/2) = −0.1659638503650, c(0.6) = 0.1345543, c(3/4) = 0.419075 (Wang's "41.9 %") ✓, c(1) = 0.6725007037, θ₀ = 0.5501939647; D(C_{1/4}) = 0.245838 ("0.2458") ✓ |
| prop:optimal: C(f_λ) by quadrature; coercivity C(f) − C_λ ≥ (1 − 2λ²/π²)‖h‖², 300 random h per λ; η_ε approximants | finite_checks.py | C(f_λ) = C_λ to 1e−7 (λ = 0.3, 0.55, 0.8, 1) ✓; min ratio 0.985, 0.948, 0.897, 0.810 ≥ 0.982, 0.939, 0.870, 0.797 ✓; C(f_ε) ↓ C_λ as ε → 0 ✓ |
| lem:coeff, x = 1 … 10⁵ (Λ by sieve to 2·10⁶, PNT tails) | finite_checks.py | Σa_n² − log x ∈ [−0.31, 0.16] ✓; Σ n a_n²/(x log²2x) ≤ 1.68 ✓ |
| Lamzouri (eq:finite), 20000 random conjugation-invariant multisets, K = f̂_λ at complex points | finite_checks.py | min (LHS − RHS) = 0 (equality: one double real point) ✓, never negative |
- Memory: footprints 29–153 MB (/usr/bin/time -l on python); finite_checks first run (Λ sieve to 4·10⁶) had max RSS 309 MB with footprint unmeasured
  (memcap sampled 58 MB) → rerun at 2·10⁶ (footprint 153 MB), identical results; run-1 log kept as finite_checks_run1_N4e6.log.
- zeros.py stopped at 1250 zeros (CPU-starved at 18 % under background QoS; ≤ 2000 as briefed).

# FULLPROOF — referee of Robin2Adic Theorem A, step by step (r2aaudit, 12 Nov 2026)

Notation as in Robin2Adic/DERIVATION.md. "✓" = re-derived here independently; "num" = also tested numerically (falsify.py, logs/).

## 1. Prop 2 (exact identity) ✓ num
- N_k = p_1⋯p_k ≤ n < N_{k+1}, x = p_k, T = ϑ(x) = log N_k, U = log n − T ∈ [0, log p_{k+1}); n = ∏_{j≤m} q_j^{a_j}, q_1 < … < q_m.
- m ≤ k: n ≥ q_1⋯q_m ≥ N_m and n < N_{k+1}. q_j ≥ p_j. ✓
- log(σ(n)/n) = Σ_j [F(q_j) − D_{a_j}(q_j)], F(q) = −log(1 − 1/q), D_a(q) = −log(1 − q^{−a−1}). ✓
- E := Σ_{i≤k} F(p_i) − γ − log log T; B := Σ_{i≤k}F(p_i) − Σ_j F(q_j); D := Σ_j D_{a_j}(q_j); H := log log(T+U) − log log T.
  E − B − D − H = Σ_j F(q_j) − D − γ − log log(T+U) = log σ(n)/n − γ − log log log n = log(σ(n)/(e^γ n log log n)). ✓
- Needs T > 1 (log log T real): n ≥ 6 ⇒ k ≥ 2, T ≥ log 6. For n > 5040: k ≥ 5. ✓ Robin (strict) ⇔ B + D + H > E. ✓
- B ≥ 0 (F decreasing, q_j ≥ p_j, m ≤ k, F > 0); D > 0; H ≥ 0. ✓
- Numerical (Arb, 320 bits): identity residual ≤ 3·10⁻⁹² on ~10⁵ n (every n ∈ [6, 40000]; all 58778 n ≤ 10²⁰ with non-increasing exponents, incl. all 120 superabundant
  n ≤ 10²⁰; 684 colossally abundant n with largest prime ≤ 10723; 2^v·(odd part of CA), v ≤ 80; 3996 random even n with v_2 ≤ 80; 3000 primorial
  multiples N_k·2^a·r, k ≤ 2500). Robin ⇔ B + D + H > E in every case; Robin fails exactly at the 23 known n ∈ [6, 5040] (σ(n) checked as an integer).

## 2. Lemma 3 (slots) ✓ num
- H(U) = ∫_0^U du/((T+u)log(T+u)) ≥ cU for any 0 < c ≤ 1/((T+U)log(T+U)). ✓
- B + D + cU = Σ_{i≤k}(F(p_i) − c log p_i) − Σ_{j≤m} V_{a_j,c}(q_j), V_{a,c}(z) = log(1 + z^{−1} + … + z^{−a}) − ca log z (uses F − D_a = log(1+…+q^{−a}),
  U = Σ_j a_j log q_j − Σ_{i≤k} log p_i). ✓ V_{a,c} strictly decreasing on z > 1 (a ≥ 1, c ≥ 0) ⇒ V(q_j) ≤ V(p_j). ✓
- Slot j ≤ m: F(p_j) − c log p_j − V_{a_j,c}(p_j) = D_{a_j}(p_j) + c(a_j − 1)log p_j ≥ D_{a_j}(p_j) (a_j ≥ 1, incl. a_j = 1). ✓
- Unused slots m < i ≤ k: F(p_i) − c log p_i > 1/x − cL ≥ 0 iff c ≤ 1/(xL); GL's c = 0.998/(xL) gives it. ✓ (For small x the hypothesis
  (T+U)log(T+U) ≤ xL/0.998 can fail; then c' = min(1/((T+U)log(T+U)), F(x)/L) works, so B + D + H ≥ Σ_j D_{a_j}(p_j) holds for every n ≥ 6.)
- Even n: q_1 = p_1 = 2, a_1 = v_2(n), slot 1 ≥ D_{a_1}(2) = −log(1 − 2^{−a_1−1}) > 2^{−a_1−1}. ✓ Odd n: slot 1 is D_{a_1}(2) with a_1 = v_{q_1}(n);
  Theorem A uses Choie–Lichiardopol–Moree–Solé (Robin for odd n > 9) instead. ✓
- Large x: T + U < x(1 + M_0) + log(2x) (BKLNW two-sided, Bertrand) ⇒ (T+U)log(T+U) ≤ xL(1 + 2.48·10⁻¹²) < xL/0.998 at L = 108.57 (cert.log l. 9);
  the ratio decreases in x. ✓
- Numerical: B + D + H − Σ_j D_{a_j}(p_j) ≥ 0 in every tested n (equality exactly when B = H = 0, e.g. n = N_k); GL's full chain with c = 0.998/(xL)
  (B + D + cU ≥ Σ_j[D_{a_j}(p_j) + c(a_j−1)log p_j] + Σ_{m<i≤k}(F(p_i) − c log p_i), and H ≥ cU) holds in every n where its hypothesis holds (~10⁵).
- Conclusion rule: E(x) < 2^{−v_2(n)−1} ⇒ Robin; tested on every case where the premise holds (3253 n): Robin holds in all.

## 3. Threshold n > 10^(10^46.79) ⇒ x ≥ e^{108.57} ✓
- log n < log N_{k+1} = ϑ(x) + log p_{k+1}; log n > Λ = 10^{46.79} log 10, log Λ = 46.79 log 10 + log log 10 = 108.5719899… (cert.log: 108.57199). ✓
- x ≤ e^{100}: ϑ(x) + log p_{k+1} < 1.01624e^{100} + log(2e^{100}) < Λ (RS62 Thm 9 holds for all x > 0; Bertrand). ✓
- e^{100} < x < e^{108.57}: < (1 + 2.453·10⁻¹²)e^{108.57} + 108.57 + log 2 < Λ (= e^{108.57}·e^{0.00199}). ✓ Both checked in Arb (cert.log l. 8).
- RobinCH's statement is "Robin for 5041 ≤ n ≤ 10^(10^46.79)" (RobinCH/REPORT.md l. 1, 16): used exactly as stated; its proof is not re-audited here.

## 4. Integrated Mertens identity (5) and G ✓
- d/dy(1/(y log y)) = −w(y), w = (log y + 1)/(y² log² y); y w(y) = d/dy(log log y − 1/log y); partial summation with ϑ(x) = x(1+d):
  Σ_{p≤x} 1/p = log log x + B_1 + d/L + ∫_x^∞ (y − ϑ)w dy. ✓ Σ_p h(p) = γ − B_1, h(p) = F(p) − 1/p > 0. ✓
- log log ϑ(x) = log L + log(1 + log(1+d)/L) ⇒ E(x) = ∫_x^∞(y−ϑ)w − Σ_{p>x} h(p) + G(d,L), G = d/L − log(1 + log(1+d)/L). ✓ Dropping −Σ_{p>x}h: upper bound. ✓
- G(0) = G_d(0) = 0, G_dd = (L + 1 + log(1+d))/((1+d)²(L + log(1+d))²) > 0 ⇒ G ≤ (d²/2) sup_{|δ|≤|d|} G_dd; |d| < M_0 = 2.4530·10⁻¹² (BKLNW, x > e^{100}). ✓

## 5. Explicit formula, kernel bound (16), low zeros ✓
- y − ψ_0(y) = Σ_ρ y^ρ/ρ + log 2π + ½log(1 − y⁻²). Termwise integration against w on [x, ∞): via ψ_1 (Ingham Thm 28, Σ|ρ(ρ+1)|⁻¹ < ∞) and integration by parts
  (boundary terms vanish since β < 1); the resulting Σ_ρ K_ρ/ρ converges absolutely by (16). ✓ ½∫log(1−y⁻²)w < 0 dropped. ✓
- y = e^u: w dy = e^{−u}f(u)du, f = 1/u + 1/u², ∫_x^∞ w = 1/(xL). K_ρ = ∫_L^∞ e^{(ρ−1)u}f(u)du; two integrations by parts, f'' > 0, β < 1 ⇒
  |K_ρ| ≤ x^{β−1}[f(L)/|1−ρ| + 2|f'(L)|/|1−ρ|²], |f'|/f = (L+2)/(L(L+1)). ✓
- |λ| ≤ H_0 = 3·10¹²: β = ½ (Platt–Trudgian); Re 1/(ρ(1−ρ)) = (β(1−β) + λ²)/|ρ(1−ρ)|² > 0 for every zero and = 1/(|ρ||1−ρ|) on β = ½, so
  Σ_{|λ|≤H_0} 1/(|ρ||1−ρ|) ≤ Σ_ρ 1/(ρ(1−ρ)) = 2Σ_ρ 1/ρ = 2 + γ − log 4π = 0.04619142. ✓ |1−ρ| ≥ |λ| ≥ 14.13 for the ε-factor. ✓
- ψ − ϑ = Σ_{k≥2} ϑ(y^{1/k}) ≤ (1 + 1.1199·10⁻⁹)√y + (log y/log 2)·1.01624 y^{1/3} (√y > e^{54}); c_ψ(y) decreasing for log y > 6;
  ∫_x^∞ √y w = ∫_L^∞ e^{−u/2}f ≤ 2f(L)e^{−L/2}. ✓
- |λ| > H_0: x^{β−1}/(|ρ||1−ρ|) ≤ x^{β−1}/λ², factor 1 + 2(L+2)/(L(L+1)H_0). ✓
- Numerical, x ≤ 10⁷ (falsify.py ex): E(x) > 0 at all 664411 primes in [10³, 10⁷]; E ≤ 0.36·M_small (the majorant with small-x constants and the actual d);
  R = E√x/f ∈ [1.872, 2.021]; the normalised zero part Z ∈ [−0.0115, 0.0310] ⊂ [−0.0468, 0.0468] = ±S_low(1+ε_lo): the low-zero bound is respected
  with room, and the x-dependence of E is the predicted 2f(L)x^{−1/2}(1 + O(x^{−1/6})).

## 6. Zeros above H_0: Abel summation on the σ-grid ✓
- W(s) := Σ_{|λ|>H_0, β>s} λ⁻² = 2∫_{H_0}^∞ t⁻² dN(s,t) = 4∫_{H_0}^∞ N(s,t)t⁻³dt (N(s,H_0) = 0 for s > ½ by RH(H_0); t⁻²N → 0). ✓
- Pairing ρ ↔ 1 − ρ̄ (same λ, multiplicities preserved): Σ_{|λ|>H_0} x^{β−1}/λ² = Σ g_{s(ρ)}(x)/λ², g_s = (x^{s−1} + x^{−s})/2, s(ρ) = max(β, 1−β);
  g_s increasing in s ∈ [½, 1) (∂_s ∝ x^{s−1} − x^{−s} ≥ 0), decreasing in x. ✓ Tail mass μ(s(ρ) > s) = 2W(s) (s ≥ ½). ✓
- Grid ½ = s_0 < … < s_J = B_UP; s(ρ) ≤ B ≤ B_UP (half-plane + reflection). Assign g_{s_j} to s(ρ) ∈ (s_{j−1}, s_j]:
  Σ ≤ g_{s_0}·M + Σ_{j=1}^{J}(g_{s_j} − g_{s_{j−1}})·μ(s(ρ) > s_{j−1}) ≤ g_{s_0}τ + Σ_j (g_{s_j} − g_{s_{j−1}})·2W̃_{j−1}, with M = total mass ≤ τ,
  W̃ = running min of τ/2 and the upper bounds W̄(s_i), i ≤ j (valid as W is non-increasing; 2W(½) ≤ τ by the pairing). This equals Σ c_j g_{s_j}
  with c_0 = τ − 2W̃_0 = 0, c_j = 2(W̃_{j−1} − W̃_j) ≥ 0, c_J = 2W̃_{J−1}. ✓ (cert.py envelope(): identical indices; grid 4964 exact dyadics.)
- τ = (log(H_0/2π) + 1)/(πH_0) + 8((a+b)(h/2 + ¼) + c/2)/H_0² ≥ 4∫_{H_0}^∞(N(t) − N(H_0))t⁻³dt from HSW two-sided (log log t ≤ log t;
  4∫R(H_0)t⁻³ = 2R(H_0)/H_0² absorbed). ✓ = 2.959407·10⁻¹².
- W̄ closed form: (log kt)^{2s}(log t)^{5−4s} ≤ (log t)^{5−2s} (1 ≤ kt ≤ t); ∫_{H_0}^∞ t^{−q−1}(log t)^e dt = Γ(e+1, qh)/q^{e+1}, q = 2 − 8(1−s)/3 > 0;
  B(s) term with max(B, 0) and I_2 = (h²/2 + h/2 + ¼)/H_0². ✓
- Monotonicity in x: f, e^{−L/2}, c_ψ (e^{−L/6}L, L > 6), both ε, 1/(xL), sup G_dd ((L+a+1)/(L+a)², derivative −(L+a+2)/(L+a)³), and
  Σ c_j g_{s_j}(x) (c_j ≥ 0 independent of x) all decrease ⇒ sup_{x ≥ e^{108.57}} E(x) ≤ majorant(108.57). ✓

## 7. KLN Lemma 4.14 at δ < 1 and at H_0 = 3·10¹² (the K = 72 vs 71 question)
Source: KLN, arXiv 2101.12263 (TeX "…-April2018.tex"; the brief's 1708.06672 is an unrelated physics paper). Line numbers = that TeX.
- l. 194 (Thm 1.1) lists "δ ≥ 1". l. 1549 (Lemma 4.14, which proves Thm 1.1 and gives C_1, C_2): "α > 0, d > 0, δ > 0"; l. 626 (parameter conventions):
  δ > 0, σ_2 = 1 + δ/log X; l. 453 (Lemma on λ_X sums): τ > 1, δ > 0, X ≥ 10⁹. KLN's Table 1 (l. 1694ff) itself uses δ = 0.303–0.314.
- Every occurrence of δ (62 × \delta, all \del) checked; δ enters only through:
  (i) l. 465–504: Σλ_X(n)²/n^τ ≤ b_4τ²e^{γ(τ−1)}log X/(τ−1) ([Ram, Lemma 5.5], quoted for τ > 1) at τ = 1 + δ/log X, 1 + 2δ/log X; and the split of
      τ = 2 + 2δ/log X using n^{2δ/log X} ≥ e^{2δ} (n > X) and (5X)^{1+δ/log X} = 5Xe^{δ}·e^{δ log 5/log X}: valid for every δ > 0;
  (ii) σ_2 = 1 + δ/log X > 1 in the Hardy–Ingham–Pólya convexity lemma (needs only σ_1 < 1 < σ_2);
  (iii) C_5, C_6 (l. 735–738) decrease in X for every δ > 0;
  (iv) M(k,δ) (l. 1015, 1062–1075): both signs of log k + 2δ treated;
  (v) l. 1077/1099: 0 < δ < log H_0(log log H_0 − 1)/2 (26.36 at KLN's H_0; 33.9 at 3·10¹²) for the monotonicity of log log T/(log kT + 2δ) ✓ (re-derived:
      derivative ∝ 1 + c/log T − log log T, c = log k + 2δ ≤ 2δ);
  (vi) l. 1639–1643: V(α,k,δ,T) > 1 "for our choice of parameters" — a numerical condition, asserted per row in Ingham4 (`Vinf > 1`).
  No step uses δ ≥ 1. ⇒ Lemma 4.14 holds for 0 < δ < log H_0(log log H_0 − 1)/2 with V > 1 (its own conditions (v), (vi)); Thm 1.1's δ ≥ 1 only restricts.
  Verdict on this step: CORRECT for the rows used (δ = 0.303–0.314), given [Ram, Lemma 5.5] as quoted (see §9 for the reviewer's wording correction).
- [Ram, Lemma 5.5] itself (Math. Comp. 85 (2016)) not seen (not on arXiv; arXiv API: 0 hits). Independent support: Ramaré–Zuniga Alterman 2603.25961 v3
  (unrefereed) Thm 1.2: S_ε(X) ≤ 0.5936 (10¹² ≤ X < 10¹⁶, ε ≤ 1/25), ≤ 0.5103, 0.4905 beyond (ε ≤ 1/50). With Σλ²n^{−1−ε} ≤ ζ(1+ε)S_ε ≤ (1+ε)S_ε/ε
  this is ≤ 0.529(1+ε)²e^{γε}log X/ε (needs S_ε ≤ 0.529(1+ε)e^{γε}log X ≈ 14.8) at every ε used (X ≥ kH_0 ≥ 1.5·10¹²: ε ≤ 2·0.314/28.04 = 0.0224 < 1/25;
  X ≥ 10¹⁶: ε ≤ 0.0171 < 1/50): margin ≈ 25. So (3.8) at δ ≈ 0.3 holds even if Ramaré's lemma had a τ-restriction.
- b_3 = 0.605 ([Ram, Lemma 5.6], (3.7); never rechecked in this project): X·Σ_{X<n<5X}λ_X(n)²/n² = 0.45321, 0.45329, 0.45341, 0.45342, 0.45342, 0.45340, 0.45341
  for X = 10⁴, 3·10⁴, 10⁵, 3·10⁵, 10⁶, 2·10⁶, 4·10⁶ (b3check.py; exact λ_X via λ_X(n) = −Σ_{e≤4, e|n, n/e>X} μ(n/e)). Stable, below 0.605. (Stated for X ≥ 10⁹: a test, not a proof.)
- H_0 = 3·10¹²: Lemma 4.14 is printed with H_0 = 3.0610046·10¹⁰ fixed (l. 1550). Using 3·10¹² re-runs its proof (Ingham4 I4′; Cully-Hugill–Johnston, IJNT 2023,
  Lemma 7.1, refereed, did the same). Spot checks here: H_0 enters via RH(H_0), T ≥ H_0 ≥ H, X ≥ kH_0 ≥ 10⁹, log H_0, log log H_0, (v), and the per-row numerical
  conditions U(α,k,H_0) > 1, positive coefficients of J(k,T) (l. 1683, "to guarantee U decreases"), V > 1 — all re-asserted in Ingham4 at 3·10¹² (row 0.60 dropped).
  Not re-read line by line here (the δ scan above covered l. 425–1110 and 1546–1694).

## 8. K and the slot form ✓
- E* ≤ 8.7461707·10⁻²³ (cert.log; rerun by the orchestrator identical). −log₂E* = 73.2757; K = ⌊73.2757⌋ − 1 = 72; 2⁻⁷³ = 1.0588·10⁻²² > E* ✓; K = 73 would need
  E* < 2⁻⁷⁴ = 5.29·10⁻²³ ✗. With δ = 1: 1.7706578·10⁻²², −log₂ = 72.258, K = 71 (2⁻⁷² = 2.12·10⁻²² > E*) ✓. CS (b_4 corrected): 7.6717883·10⁻²⁴, K = 75 ✓.
- Rule: v_2(n) = a ≤ K ⇒ D_a(2) > 2^{−a−1} ≥ 2^{−K−1} > E* ≥ E(x) ⇒ B + D + H > E ⇒ Robin. ✓ Odd n: CLMS. n ≤ 10^(10^46.79): RobinCH. ✓
- Slot form: p_j^{a_j+1} > 1/E* ≥ 1.14336·10²²: a_1 ≥ 73, a_2 ≥ 46, a_3 ≥ 31, a_4 ≥ 26, a_5 ≥ 21, a_6 ≥ 19 ✓ (log(1/E*)/log p_j = 73.28, 46.23, 31.56,
  26.10, 21.18, 19.80). **Erratum:** "a_j ≥ 2 whenever p_j ≤ 1.07·10¹¹" is rounded the wrong way: √(1.14336·10²²) = 1.06928·10¹¹, so the clause is
  proved only for p_j ≤ 1.0692·10¹¹ (a_j = 1 is not excluded for p_j ∈ (1.06928·10¹¹, 1.07·10¹¹]). Theorem B's 1.30·10²³, a_1 ≥ 76: ✓.

## 9. Reruns, reproduction, reviewer
- cert.py rerun on a copy (/tmp, unchanged code): identical to logs/cert.log (KLN 8.7461707·10⁻²³ → K = 72; δ = 1: 1.7706578·10⁻²² → K = 71; CS corrected
  7.6717883·10⁻²⁴ → K = 75; GL (20) at L = 106: 1.227034555·10⁻²³); peak RSS 39 MB, 23 s.
- The KLN implementation used for H_0 = 3·10¹² (Ingham4 `KP`) reproduces all 20 rows of KLN Table 1 at KLN's H_0 = 3.0610046·10¹⁰ with δ as printed:
  max relative deviation 4.4·10⁻⁴ in A and B (kln_table1.py, logs/kln_table1.log; all 19 rows at 3·10¹² have Vinf > 1). At 3·10¹²: 19 rows admissible (0.60 dropped: U ≤ 1), δ ∈ [0.303, 0.307].
- Vinf(α,k,δ) = 8α(C_5 I(1,0) + C_6 I(2,0)/k) is the T → ∞ value of V; with V decreasing in T (KLN l. 1629) the per-row check Vinf > 1 gives V > 1 for all T ≥ H_0.
- Reviewer (verdict-only, the δ step pasted inline): "GAP in 'every δ > 0', not in allowing δ < 1. No quoted step requires δ ≥ 1, and δ = 0.303 satisfies the displayed
  parameter restrictions." The gap is my wording: Lemma 4.14's proof also needs δ < log H_0(log log H_0 − 1)/2 (l. 1077/1099) and V > 1 (l. 1639). Corrected
  statement: **Lemma 4.14 holds for 0 < δ < log H_0(log log H_0 − 1)/2 (33.9 at 3·10¹²) whenever V > 1; δ ≥ 1 is never used.** The table rows (δ ≈ 0.30) meet
  both, so K = 72 stands. The reviewer also notes the RZA fallback covers only X ≥ 10¹²: that is the whole range used (X = kT ≥ 0.6·3·10¹² for the admissible rows).

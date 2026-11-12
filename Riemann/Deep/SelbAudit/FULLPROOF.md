# FULLPROOF — selbaudit: line-by-line referee of the SelbergN4 chain (DERIVATION.md §1–§10)

Source: Dousselin, arXiv 2311.10285 (e-print = the TeX used by DERIVATION; same line numbers, checked: mainth l.613, eq:N l.656, min-J l.696, estiI l.711, eqIntfacile l.1586, eqintero l.1644, gainSel l.1688).
Fixed throughout: q, χ_j, c_j ≠ 0, N, θ ∈ (0,1), ϰ ∈ (0,1/8), A > 1/ϰ. Every o(1) is T → ∞ with these fixed; §10 adds P_0 (fixed), then K, D_0 → ∞ after T.

## (i) Counting: m ≥ meas(G)/H
G := {t ∈ (T,2T): 𝓕(1/2+i·) takes both signs on (t,t+H)}; G is open (two points of opposite sign persist under small shifts of t), E := (T,2T)∖G is relatively closed.
Let γ_1 < … < γ_m be the sign-change points of 𝓕 in (T, 2T+H) (zeros of odd order of F on the line; finitely many). If t ∈ G, 𝓕 changes sign at some γ_i ∈ (t, t+H), so t ∈ (γ_i − H, γ_i).
Hence G ⊂ ∪_i (γ_i − H, γ_i), meas G ≤ mH, m ≥ (T − meas E)/H = (T − meas E)log T/A. Zeros in (2T, 2T+H]: ≤ N(2T+H,F) − N(2T,F) = O(log T) by eq:N (l.656), negligible against T log T.
So N_0(2T) − N_0(T) ≥ (T − meas E)log T/A − O(log T) (N_0 counting with or without multiplicity: distinct odd-order zeros are a lower bound for both).
Comparison: Dousselin (l.736–741, 897–902) pairs abutting intervals I_1, I_2 and counts ⌊T/(2H)⌋ − n_1 with n_1H ≤ meas E: main term T/(2H). Both subtract meas(E)/H; only the main term doubles. CORRECT.

## (ii) Combination inside each dominance set
Dousselin's construction (l.811–878): S_j ⊂ (T,2T) (disjoint, j = the dominating index), S_j^* := S_j minus {t: ∫_t^{t+H}|Lη_j²|² > H loglog T}; Σ_j meas S_j^* = T + O(T/(loglogT)^{1/2−δ}) (sizeS, l.873).
(I*), (J*) (l.894, 898) hold POINTWISE for every t ∈ S_j^* (proved for each such t from eq_LL1, eq:int_Ht, mesureW, intmollif): I_j^* = c_jI_j + O(H/loglogT), J_j^* = |c_j|J_j + O(H/loglogT), O-constant depending on N only.
Let S'_j := S_j^* ∩ E. For t ∈ S'_j, 𝓕(t+u) has constant sign on (0,H) ⊃ H_t, so |I_j^*| = J_j^* (both integrals over H_t, weight |η_j²| ≥ 0). Then
  |c_j||I_j| + O(H/loglogT) ≥ |I_j^*| = J_j^* ≥ |c_j|J_j − O(H/loglogT) ≥ |c_j|(H − |M_j|) − O(H/loglogT)   (min-J, l.696),
so |I_j| + |M_j| ≥ H(1 − ε_T) on S'_j with ε_T = O(1/(|c_j| loglogT)) → 0 (c_j fixed). S'_j is measurable (S_j^* is, E is closed). Cauchy–Schwarz:
  (1 − ε_T)H·meas S'_j ≤ ∫_{S'_j}|I_j| + ∫_{S'_j}|M_j| ≤ √(meas S'_j)(‖I_j‖₂ + ‖M_j‖₂),  ‖·‖₂ on (T,2T).
With estiI/estiM (l.711–713: ‖I_j‖₂² ≲ C1 T/log²T, ‖M_j‖₂² ≲ C2 T/log²T, the same C1, C2 for every j): meas S'_j ≲ (√C1+√C2)² T/A².
E ⊂ ((T,2T)∖∪_jS_j^*) ∪ ∪_j S'_j, so meas E ≤ o(T) + N(√C1+√C2)²T/A²(1+o(1)). With (i) and eq:N (N(2T)−N(T) = (T/2π)log T(1+o(1))):
  κ_F ≥ 2π(1/A − N(√C1+√C2)²/A³).   CORRECT (it is Dousselin's own N = 1 argument, l.722–733, run on S_j^* instead of (T,2T)).
Note: Dousselin's per-j step (l.900–904) bounds |I_j| ≤ H/(2+ε) and |M_j| ≤ H/(2+ε) separately (union bound): (2+ε)²(C1+C2) ≥ 2(√C1+√C2)²; the new constant is ≤ half of it.

## (iii) M-estimate: C2 = 4(C3 − 1)/(θϰ)²
Dousselin (l.1366–1405): M(t,H) = a(t) + b(t) − c(t), a = ∫_t^{t+H}f(3/2+iu)du, b(t) = ∫_{1/2}^{3/2}f(σ+it)dσ, c(t) = b(t+H), f := Lη² − 1 (no pole: L(s,χ) entire, χ primitive mod q ≥ 3; η a polynomial).
His bounds, unchanged: ‖a‖₂² ≪ Tξ^{−θ/2+o(1)} (first-esti, l.1385; note it is ≪ T·(Σ_{n≥ξ^θ}τ_3(n)n^{−3/2})² pointwise, so o(T/log²T)); Gabriel (l.1389–1401):
‖b‖₂² ≲ Q/(θ log ξ)², where Q is ANY upper bound for ∫_T^{2T}|f(1/2+it)|²dt (his display is monotone in Q: Q^{3/2−σ}(Tξ^{−2θ}/2)^{σ−1/2}, and Q enters the σ-integral only through
(2Q/T)^{−(σ−1/2)}, which gives (1+o(1))/(θ log ξ) for every constant ratio Q/T > 0). Same for c (range (T+H, 2T+H); the extra length H adds O(∫_{2T}^{2T+H}|f|²) = O(T^{1/2+ε})).
Minkowski: ‖M‖₂ ≤ ‖a‖₂ + ‖b‖₂ + ‖c‖₂ ≲ 2√Q/(θϰ log T) (log ξ = ϰ log T), i.e. ‖M‖₂² ≲ 4Q/((θϰ)² log²T) instead of 3(‖a‖²+‖b‖²+‖c‖²) ≲ 6Q/(...).
First moment. Rectangle [1/2,3/2]×[T,2T], f entire: ∫_{1/2+iT}^{1/2+2iT} f ds = ∫_{3/2+iT}^{3/2+2iT} f ds + (horizontal sides).
 σ = 3/2: f = Σ_{n≥ξ^θ} c(n)n^{−s}, |c(n)| ≤ τ_3(n) (l.1378–1384: β = α below ξ^θ), so |∫| = |Σ c(n)n^{−3/2}(n^{−2iT}−n^{−iT})/log n| ≤ 2Σ τ_3(n)n^{−3/2}/log n = O(1).
 Sides: |η(σ+iT)| ≤ Σ_{n≤ξ}|τ_{−1/2}(n)M(n)|n^{−1/2} ≤ 2√ξ (|τ_{−1/2}| ≤ 1, 0 ≤ M ≤ 1); L(σ+iT,χ) ≪_q T^{1/4+ε} on 1/2 ≤ σ ≤ 3/2 (convexity suffices; Weyl's 1/6 is not needed).
 So each side is ≤ ∫_{1/2}^{3/2}(|L||η|² + 1)dσ ≪ T^{1/4+ϰ+ε} ≤ T^{3/8+ε} = o(T). Hence ∫_T^{2T} Lη²(1/2+it) dt = T + o(T).
Then ∫|f|² = ∫|Lη²|² − 2Re∫Lη² + T ≤ C3T − 2T + T + o(T) = (C3 − 1 + o(1))T (lemme-mollifer: ∫|Lη²|² ≲ C3T). Take Q = (C3 − 1 + o(1))T:
  C2 = 4(C3 − 1)/(θϰ)². CORRECT. (Dousselin: Q = (√C3+1)²T by Minkowski, factor 6.) C3 > 1 holds for all certified θ (C3 ≥ 18.5).

## (iv) Exact main terms from Lemma 12 (DERIVATION §4)
Inputs (Dousselin, authenticated in R416, not re-proved): Lemma 12 (l.1107–1115): for (d,q)=1, 0 < v ≤ 1/√log ξ, r = 2, 3,
  S_r(y,v) := Σ_{(n,qd)=1, n<y} τ_{−1/2}(n) log^{r−1}(y/n) n^{−1−iv} = −((r−1)!/(2√π)) Π_{qd}(v) v^{3/2−r} Δ_r(v log y) + O(∏_{p|d}(1+p^{−3/4})((log T)^{r−5/2} + (loglog T)^{r−3/2})),
  y = ξ/δ; in its proof the main term is (r−1)!·𝔓_{qd}(0)·H_r with H_r(y,v) := (1/2πi)∫_{ε−i}^{ε+i} y^s s^{−r}√(s+iv) ds (l.1124–1133).
Exact structure (l.1430–1438, 1566): for δ | d^∞, (n,qd) = 1 ⇒ (n,δ) = 1, so τ_{−1/2}(nδ) = τ_{−1/2}(n)τ_{−1/2}(δ) and
  S_d(δ,v) = τ_{−1/2}(δ)/((1−θ)log ξ)·[S_2(ξ/δ,v) − S_2(ξ^θ/δ,v)],  B_d(δ) = S_d(δ,0): the δ-dependence is EXACTLY τ_{−1/2}(δ) × (function of log δ).
4.1 (closed forms; all re-derived): φ' = e^{−ix}x^{−3/2}; ψ := Xφ − Δ_2 has ψ' = Xφ' = e^{−iX}X^{−1/2}, ψ(0) = 0, ψ(X) = 2∫_0^{√X}e^{−iy²}dy = √π e^{−iπ/4}erf(e^{iπ/4}√X);
  integration by parts in φ's definition gives φ(x) = −2x^{−1/2}e^{−ix} − 2iψ(x), hence Δ_2(X) = Xφ − ψ = −2√X e^{−iX} − (1+2iX)ψ(X). Termwise integration of φ's series:
  Δ_2(X) = −4√X + Σ_{j≥1}(−i)^jX^{j+1/2}/(j!(j²−1/4)) = √X·E(X). ψ = ψ∞ + ie^{−iX}X^{−1/2} + O(X^{−3/2}) gives Δ_2 = c1X − ψ∞ + O(X^{−1/2}), c1 = −2iψ∞ = 2√π e^{−3iπ/4};
  exactly R0(Y) := Δ_2(Y) − c1Y = −ψ(Y) − Y∫_Y^∞e^{−ix}x^{−3/2}dx, so |R0| ≤ √π + 4/√Y and |ψ| ≤ √π + 2/√Y, |ψ(Y)| ≤ 2√Y. CORRECT (mpmath check: logs/counting.log (b)).
4.2 With Π := Π_{qd}(v), X := v log ξ, u := log δ/log ξ, ℓ1 = 1 − u, ℓ2 = (θ − u)⁺: Δ_2(v log(ξ/δ)) = Δ_2(Xℓ1), so S_d(δ,v) = Π s(δ,v) + O(τ_{1/2}(δ)∏(1+p^{−3/4})(log T)^{−2/3})
  (= eqintero l.1638–1644, valid on 0 < v ≤ H since H ≤ 1/√log ξ for large T; the factor [ζ(1+iv)iv]^{−1/2} = 1 + O(v) dropped from 𝔓_{qd}(0) costs O(v·√log ξ) = O(A/√log T) = o(1)).
  s = −τ_{−1/2}(δ)√X W_X(u)/((1−θ)2√π√v log ξ) = −τ_{−1/2}(δ)W_X(u)/((1−θ)2√π√log ξ), W_X(u) := √ℓ1E(Xℓ1) − √ℓ2E(Xℓ2), |s| ≤ τ_{1/2}(δ)S_X(u)/((1−θ)2√π√log ξ)
  (|τ_{−1/2}(p^a)| = τ_{1/2}(p^a)/(2a−1) ≤ τ_{1/2}(p^a), so "=" in DERIVATION should read "≤"; harmless). CORRECT.
4.3 v → 0: S_d(δ,v) is a finite sum, continuous in v, and the error bound is uniform in v ∈ (0,H], so B_d(δ) = Π_{qd}(0)b(δ) + O(same), b(δ) = lim_{v→0}s(δ,v)
  = τ_{−1/2}(δ)(2/√π)(√ℓ1 − √ℓ2)/((1−θ)√log ξ) (E(0) = −4). B̂(u) = (√ℓ1−√ℓ2)/√(1−θ): increasing on [0,θ] (derivative (1/√(θ−u) − 1/√(1−u))/(2√(1−θ)) > 0), decreasing after, B̂(θ) = 1.
  |b(δ)| ≤ τ_{1/2}(δ)·2B̂(u)/(√π√(1−θ)√(ϰ log T)) ≤ τ_{1/2}(δ)C5'/√log T, C5' = 2/√(πϰ(1−θ)) (3.461711 at θ = 0.15, cert.log). CORRECT.
4.4 Im part (the derivative route; this is the step that replaces coro-molli2's C4). f(u) := δ^{−iu}S_d(δ,u)/Π_{qd}(u), g(u) := δ^{−iu}s(δ,u); f(0), g(0) real (B_d, b real).
  ∂_u(δ^{−iu}S_2(y,u)) = −iδ^{−iu}[log(yδ)S_2 − S_3] (log n = log y − log(y/n)); log(yδ) = log ξ (y = ξ/δ) or θ log ξ (y = ξ^θ/δ).
  Main terms: Lemma 12 gives S_2 ≈ ΠH_2, S_3 ≈ 2ΠH_3; integration by parts in s (boundary terms at ε ± i are O(1)·y^ε = O(1)): L H_2 − 2H_3 = −(1/2πi)∫y^s s^{−2}(2√(s+iu))^{−1}ds = i∂_uH_2 (L = log y).
  So main[∂_u(δ^{−iu}S_2)] = −iδ^{−iu}Π[log δ·H_2 + i∂_uH_2] = Π∂_u(δ^{−iu}H_2), error O(∏(log ξ(loglogT)^{1/2} + (log T)^{1/2} + (loglogT)^{3/2})).
  Dividing by (1−θ)log ξ and by Π (|Π|^{−1} ≤ ∏_{p|qd}(1+1/p)^{1/2} ≪ ∏(1+p^{−3/4})): f'(u) = g'(u) + O(τ_{1/2}(δ)∏(1+p^{−3/4})²(loglogT)^{1/2}) + O(|f||Π'/Π|),
  |Π'/Π| ≤ ½Σ_{p|qd}log p/(p−1) ≪ loglog T, |f| ≪ τ∏(log ξ)^{−1/2}; all uniform in 0 < u ≤ H. Hence
  |Im f(v)/v − Im g(v)/v| ≤ sup_{(0,v]}|f' − g'| = o(τ_{1/2}(δ)√log T·∏(1+p^{−3/4})²)   (g is smooth in u: W_X is entire in X).
  (DERIVATION §4.4 quotes the same error once normalised — coro-molli2's (log T)^{1/3} — and once before dividing by log ξ — log T·√loglogT against (log ξ)^{3/2}; both are o(main).)
  Exactly: Im(g(v))/v = −τ_{−1/2}(δ)√log ξ·F4_X(u)/((1−θ)2√π), F4_X(u) = Im(e^{−iXu}W_X(u))/X = cos(Xu)(ℓ1^{3/2}IE(Xℓ1) − ℓ2^{3/2}IE(Xℓ2)) − u·sinc(Xu)Re W_X(u),
  IE(Y) = Im E(Y)/Y = Σ_k(−1)^{k+1}Y^{2k}/((2k+1)!((2k+1)²−1/4)); F4_0(u) = 4u(√ℓ1−√ℓ2) − (4/3)(ℓ1^{3/2}−ℓ2^{3/2}). CORRECT.

## (v) Lemma Φ and Lemma W (DERIVATION §5)
Lemma Φ. d | δδ' ⇒ δδ' ≥ d ⇒ u_δ + u_{δ'} ≥ t_d; terms with δ or δ' ≥ ξ vanish (M(nδ) = 0), so u, u' ∈ [0,1); with f, g ≥ 0 each term has f(u)g(u') ≤ Φ_{f,g}(t_d).
  Σ_{d|δδ'|d^∞}τ_{1/2}(δ)τ_{1/2}(δ')/(δδ') = ∏_{p^k||d}Σ_{n≥k}[x^n](1−x)^{−1}·p^{−n} = ∏ 1/(p^{k−1}(p−1)) = 1/φ(d) (logs/deltasums.log (a)). Box version: if u ∈ box i, u' ∈ box j and
  u + u' ≥ t then u_i^+ + u_j^+ ≥ t, so Φ(t) ≤ max{f_i g_j : u_i^+ + u_j^+ ≥ t} with f_i, g_j box sups; Φ is non-increasing, so Φ(t_m) bounds it on [t_m, t_{m+1}]. CORRECT.
  For f = g = B̂ the sup is explicit: Φ_BB(t) = min(1, (1 − t/2)/(1−θ)) (max of (1−u)(1−u') on u+u' = t at u = u'), so ∫_0^2Φ_BB² = 2θ + 2(1−θ)/3 (0.86667 at θ = 0.15; cert step bound 0.87381).
Lemma W. d²/φ(d)³ = Σ_{kn=d}μ²(k)ψ(k)/n (local check: 1 + pψ(p) = (p/(p−1))³). For each k, partial summation of Σ_{n≤z}χ0(n)/n = (φ(q)/q)log z + c_q + O(1/z) against
  F((log k + log n)/log ξ): (φ(q)/q)log ξ∫_{t_k}^{2}F + O(sup F + Var F) ≤ (φ(q)/q)log ξ∫_0^2F + O(·) (F ≥ 0); Σ_kμ²ψχ0(k) ≤ P1. The additive O(1) is o(log ξ) for every F used
  (bounded variation uniformly in the parameters X, v: step functions with ≤ 200 steps times |sin(X/(2ϰ) + tX)|/X, X ∈ (0,1]). CORRECT.


## (vi) ∫|I|², small v (0 < v ≤ 1/log ξ; DERIVATION §7)
Representation (Dousselin, reused): ∫_T^{2T}|I|² = 4T(φ(q)/q)∫_0^H(H−v)Im[K(v)/v]dv + O(T/log³T) (eqIntfacile l.1586), K(v) = 𝒯^{iv/2}Σ_{d≤ξ²,(d,q)=1}φ_{iv}(d)A'_d(v)²,
A'_d = Σ_{d|δδ'|d^∞}S_d(δ,v)B_d(δ')/(δ^{1+iv}δ') (eqExpK l.1563, eqa' l.1566), φ_{iv}(d) = d^{1+iv}∏_{p|d}(1−p^{−1−iv}) (l.1665), 𝒯 = qT/(2π).
Identity: with Ã'_d := A'_d/(Π_{qd}(v)Π_{qd}(0)), (d,q) = 1, Π_{qd}(v)² = ∏_{p|q}(1−p^{−1−iv})^{−1}∏_{p|d}(1−p^{−1−iv})^{−1}, Π_{qd}(0)² = (q/φ(q))(d/φ(d)):
  𝒯^{iv/2}φ_{iv}(d)A'_d² = e^{ivΛ_d}(q/φ(q))²(d²/φ(d))ε_q(v)Ã'_d², Λ_d = log(√𝒯 d), ε_q(v) = (φ(q)/q)∏_{p|q}(1−p^{−1−iv})^{−1} = 1 + O_q(v). Verified term by term. CORRECT.
Bound: |Im(εe^{iα}w²)| ≤ |sin α||w|² + 2|w||Im w| + |ε−1||w|² (Im(e^{iα}w²) = sin α Re w² + cos α Im w², |Im w²| = 2|Re w||Im w|). b̃ := B_d/Π_{qd}(0) is real (B_d is a real
  sum), so Im Ã'_d = Σ b̃(δ')Im(δ^{−iv}s̃(δ,v))/(δδ'). By (iv), Lemma Φ: |Ã'_d| ≲ Φ6(t_d)/(π(1−θ)^{3/2}log ξ φ(d)), |Im Ã'_d|/v ≲ Φ4(t_d)/(π(1−θ)^{3/2}φ(d))
  (products 2/(2√π·√π) = 1/π of the two normalisations; τ_{1/2} sums give 1/φ(d)). Errors of (iv) are relative o(1) uniformly in d, δ and absolutely summable after the d-sum
  (weights (d/φ(d))∏(1+p^{−3/4})^k/d: Euler products converge). |sin(vΛ_d)| = |sin(X/(2ϰ)·(1 + log(q/2π)/log T) + Xt_d)|: function of t_d for fixed v.
  Lemma W (F = Φ6²|sin|/X and F = 2Φ6Φ4; the log ξ from Lemma W cancels one 1/log ξ): |Im(K/v)| ≲ (q/φ(q))P1/(π²(1−θ)³)·W(X), W as in DERIVATION §7. The |ε−1| term is
  O_q(v)·O(1/log ξ), giving O(TA/log³T) after ∫_0^{1/log ξ}. With v = X/log ξ, H − v = (Aϰ − X)/log ξ:
  4T(φ/q)∫_0^{1/log ξ}(H−v)Im(K/v)dv ≲ (T/log²T)(4P1/(π²ϰ(1−θ)³))∫_0^1(A − X/ϰ)W(X)dX. Constants re-derived. CORRECT.
  Comparison: Dousselin keeps |φ_{iv}(d)| ≤ σ(d) and (d/φ(d))² from |A'_d|² (l.1601–1630): local factor P2 = 79.886; the identity leaves d²/φ(d)³: P1 = 11.543.

## (vii) ∫|I|², large v (1/log ξ < v ≤ H; DERIVATION §8): integration by parts
Representation (Dousselin intKV l.1662–1669, after cancelling ∏_{p|d}(1−p^{−1−iv}) against eqADcarré's ∏_{p|qd}(1−p^{−1−iv})^{−1}): the part is
  (1/(4(1−θ)²π log²ξ))Σ_d ∫_a^H(H−v)Im[ε_q(v)(q/φ(q))𝒯^{iv/2}d^{1+iv}(Σ τ_{−1/2}(δ)B_d(δ')Δ(δ,v)/v/(δ^{1+iv}δ'))²]dv + O(A²log(Aϰ)/(log T)^{13/6}), a = 1/log ξ
  (error: Im(K/v) = O(A/(v(log T)^{7/6})) + main, l.1668, integrated with weight H − v over [a,H]). Note: Dousselin's lastone (l.1671–1676) moves ∏_{p|q}(1−p^{−1−iv})^{−1} out of the
  modulus of a v-integral, which is not an identity; it is harmless because ∏_{p|q}(...)^{−1} = (q/φ(q))(1 + O_q(v)) and the O_q(v) part costs O_q(H)∫_a^H(H−v)|D1D2|dv = O_q(A³/log T) = o(1).
  Expanding the square and taking absolute values over (δ1,δ1',δ2,δ2') (triangle inequality) leaves |J|, J := ∫_a^H(H−v)e^{ivΛ}D1D2dv, Λ = log(√𝒯d/(δ1δ2)), D_i = Δ(δ_i,v)/v = log ξ·d(x;u_i),
  x = v log ξ ∈ [1, Aϰ], d(x;u) = [Δ_2(xℓ1) − Δ_2(xℓ2)]/x, d'(x;u) = [ℓ1Δ_2'(xℓ1) − ℓ2Δ_2'(xℓ2)]/x − d/x = [ψ(xℓ1) − ψ(xℓ2)]/x² (Δ_2' = φ, Yφ − Δ_2 = ψ).
IBP. P := (H−v)D1D2, P(H) = 0, P' = −D1D2 + (H−v)(D1D2)'. J = −P(a)e^{iaΛ}/(iΛ) + (iΛ)^{−1}∫_a^HD1D2e^{ivΛ} − (iΛ)^{−1}∫_a^H(H−v)(D1D2)'e^{ivΛ}; the middle integral once more:
  |(iΛ)^{−1}∫_a^HD1D2e^{ivΛ}| ≤ Λ^{−2}[|D1D2(H)| + |D1D2(a)| + ∫_a^H|(D1D2)'|]. With |D_i(a)| = log ξ·S1(u_i), |D_i| ≤ log ξ·Dm(u_i), ∫_a^H|D_i'|dv ≤ log ξ·Dp(u_i) (dv = dx/log ξ,
  D_i' = (log ξ)²d'), (H − a)log ξ = Aϰ − 1, H log ξ = Aϰ:
  |J| ≤ (log ξ/Λ)[(Aϰ−1)S1S1 + Aϰ(Dm1Dp2 + Dm2Dp1)] + (log ξ/Λ)²[S1S1 + Dm1Dm2 + Dm1Dp2 + Dm2Dp1].   CORRECT (exact; numerically |J|/bound ≤ 0.45, logs/ibpcheck.log).
  Finiteness: d is bounded on [1,∞) (d → c1(ℓ1 − ℓ2), |R0| bounded) and |d'| ≤ 2·Psimax/x² is integrable: this is what removes Dousselin's ∫_a^H(H−v)/v dv = A log A term (gainSel, l.1688).
Λ: δ_i < ξ (else Δ-terms vanish) gives Λ ≥ ½log𝒯 + (t_d − 2)log ξ, so log ξ/Λ ≤ ϰ/(½ + (t_d−2)ϰ + O(1/log T)) = ϰr(t_d)(1+o(1)), r(t) ≤ r(0) = 1/(½ − 2ϰ) < ∞ (ϰ < 1/4).
Box bounds: u ≥ θ: Dm ≤ √(1−u)max(sup_{[0,1]}|E|, sup_{Y≥1}|Δ_2(Y)|/Y) (Y = xℓ1 ≤ 1: |Δ_2(Y)|/x = √ℓ1|E(Y)|/√x; Y ≥ 1: ℓ1|Δ_2(Y)|/Y), Dp = ℓ1∫_{ℓ1}^∞|ψ|Y^{−2};
  u < θ: d = c1(1−θ) + (R0(xℓ1) − R0(xℓ2))/x gives the x ≥ 10 tail |c1|(1−θ) + 2R0max/10, and ∫_{200}^∞|d'| ≤ 2Psimax/200. CORRECT.
Sum: |B_d(δ')| ≤ √(dq/(φ(d)φ(q)))τ_{1/2}(δ')·2B̂(u')/(√π√(1−θ)√log ξ) (4.3); the 4-fold sum factorises into (δ1,δ1') × (δ2,δ2') sums of f(u)B̂(u') with f ∈ {S1, Dm, Dp},
  bounded by Lemma Φ; Lemma W with F(t) = ϰr(t)[(Aϰ−1)Φ_{S1}² + 2AϰΦ_{Dm}Φ_{Dp}] + ϰ²r²[Φ_{S1}² + Φ_{Dm}² + 2Φ_{Dm}Φ_{Dp}] gives
  this part ≲ (T/log²T)(8P1/(π²ϰ(1−θ)³))[Aϰ(aS1+aX) − aS1 + ϰa2]  (prefactor 4T(φ/q)·(q/φ)²/(4π(1−θ)²log²ξ)·4/(π(1−θ)log ξ)·(φ/q)P1 log ξ, re-derived). CORRECT.
  C1 is therefore linear in A: Ka = pref(La + Sa), Kb = pref(Lb + Sb).


## (viii) Signed local factors at the small primes (DERIVATION §10)
Notation: P_0 = {p ≤ 10⁶, p ∤ q} (primes dividing q never occur: (d,q) = 1), P = ∏_{P_0}p; d = d_0d_1, δ = δ_0δ_1, δ' = δ_0'δ_1' (index 0: P_0-supported; index 1: coprime to P).
d | δδ' | d^∞ holds iff d_0|δ_0δ_0'|d_0^∞ and d_1|δ_1δ_1'|d_1^∞ (prime by prime); τ := τ_{−1/2} is multiplicative; u_δ = u_{δ_1} + ε_0, ε_0 = log δ_0/log ξ.
Local values: (Σ_aτ(p^a)x^a)² = 1 − x, so σ_p(k,0) := Σ_{a+b≥k}τ(p^a)τ(p^b)p^{−a−b} = −1/p (k = 1), 0 (k ≥ 2); hence Σ_{d_0|δ_0δ_0'|d_0^∞}τ(δ_0)τ(δ_0')/(δ_0δ_0') = μ²(d_0)(−1)^{ω(d_0)}/d_0
  (absolutely convergent: the |τ|-series is ≤ Σ τ_{1/2}τ_{1/2}/(δ_0δ_0') = 1/φ(d_0)). Checked numerically (logs/deltasums.log (a): −1/p and 0 to 6 digits, p ≤ 13, k ≤ 4; drift ∝ v).
Lemma S (one factor). Let F, G: [0,∞) → ℂ be continuous, = 0 on [1,∞), with modulus of continuity ω of (u,u') ↦ F(u)G(u'), and
  Σ_d := Σ_{d|δδ'|d^∞}τ(δ)τ(δ')F(u_δ)G(u_{δ'})/(δδ'). For d_0 fixed and K ≥ 1:
  |Σ_d − (μ²(d_0)(−1)^{ω(d_0)}/d_0)·I_1(d_1)| ≤ [ω(2log K/log ξ)/φ(d_0) + 2η_K(d_0)sup|F|sup|G|]/φ(d_1),  I_1(d_1) := Σ_{d_1|δ_1δ_1'|d_1^∞}τ(δ_1)τ(δ_1')F(u_{δ_1})G(u_{δ_1'})/(δ_1δ_1'),
  η_K(d_0) := Σ_{max(δ_0,δ_0')>K}τ_{1/2}(δ_0)τ_{1/2}(δ_0')/(δ_0δ_0') → 0 (K → ∞). Proof: split the (δ_0,δ_0')-sum at K; for δ_0, δ_0' ≤ K replace F(u_{δ_1}+ε_0)G(u_{δ_1'}+ε_0') by
  F(u_{δ_1})G(u_{δ_1'}) (error ω(2log K/log ξ), absolute δ-sums 1/φ(d_0)·1/φ(d_1)); the remaining δ_0-sum differs from its full value by ≤ η_K. And |I_1| ≤ Φ_{|F|,|G|}(t_{d_1})/φ(d_1) (Lemma Φ for d_1).
  Uses: §6 (F = G = B̂-profile), §7 (F = e^{−iXu}W_X(u) or the F4_X profile, G = B̂; the phase δ^{−iv} is inside F, so no separate phase; ω uniform in X ∈ [0,1] by compactness;
  the Im part is formed per δ before summing, as DERIVATION and its reviewer say).
Extension to §8 (DERIVATION §10 states every sum has the one-factor form; the large-v sum does not — J couples δ1 and δ2 — but the argument extends):
  the 4-fold sum Σ τ(δ1)τ(δ1')τ(δ2)τ(δ2')b̂(u_{δ1'})b̂(u_{δ2'})J(δ1,δ2)/(δ1δ1'δ2δ2') has J depending on δ_{i0} only through u_{δ_i} and the phase (δ_{10}δ_{20})^{−iv}·(d_0)^{iv}. For δ_{i0} ≤ K, d_0 ≤ D_0:
  |J(δ1,δ2) − J(δ_{11},δ_{21})| ≤ ∫_a^H(H−v)[|D1D2|·v·log(D_0K²) + log²ξ·2Dm·ω_d(log K/log ξ)]dv ≤ (Aϰ)³Dm²log(D_0K²)/log ξ + (Aϰ)²Dm·ω_d → 0,
  ω_d := modulus of continuity of u ↦ d(x;u), uniform on x ∈ [1,Aϰ] (compact). So J(δ1,δ2) = J(δ_{11},δ_{21}) + o(1) uniformly; the (δ_{10},δ_{10}') and (δ_{20},δ_{20}') sums then separate
  and each is Lemma S's local value; |J(δ_{11},δ_{21})| is bounded by (vii) at the d_1 level. Factor (μ²(d_0)/d_0)², as for the squares in §6–§7. CORRECT with this extension.
Σ2' (§6): Ã_d = Σ_kΛ(k)A_{[d,k]}. k = p^m, p ∉ P_0: [d,k] has P_0-part d_0, signed factor μ²(d_0)/d_0, and the d_1-part as in §6: Σ_kΛ(k)[d_1,k]/φ([d_1,k])² ≤ (d_1/φ(d_1)²)(log d_1 + 2ϰlog T + o(log T))
  (Dousselin l.1476–1484 before his step log d ≤ 2ϰ log T). k = p^m, p ∈ P_0: the signed main of A_{[d,k]} vanishes unless v_p([d,k]) = 1 (m = 1 and v_p(d) ≤ 1); these add
  ≤ (log d_0 + Σ_{p∈P_0}log p/(p−1))·(|A_d|-bound) = O_{P_0}(1)·(|A_d|-bound) = o(log T)·(…) (the extra Π_{qdp}(0)²/Π_{qd}(0)² = p/(p−1) is the 1/(p−1)). CORRECT.
d-sums. Weights become (q/φ(q))²·(μ²(d_0)/φ(d_0))·d_1²/φ(d_1)³ (§6: φ(d)(d/φ(d))²/d_0²φ(d_1)²; §7–§8: (d²/φ(d))/d_0²φ(d_1)²). Lemma W over (d_1,qP) = 1 (same proof with χ0 mod qP):
  Σ_{d_1}(d_1²/φ(d_1)³)F(t_{d_1}) ≤ (φ(qP)/(qP))∏_{p∤qP}(1+f1(p))log ξ∫_0^2F + O_{qP}(1). Σ_{d_0}μ²(d_0)/φ(d_0) ≤ ∏_{p∈P_0}p/(p−1), and (φ(qP)/(qP))∏_{P_0}p/(p−1) = φ(q)/q:
  the Euler factor of each p ∈ P_0 is exactly 1. Φ and r are evaluated at t_{d_1} ≤ t_d (Φ(t_{d_1}) ≥ Φ(t_d): the d_1-level bound is the right one; r(t_d) ≤ r(t_{d_1})).
  Remaining factor P1eff = ∏_{p>10⁶, p∤q}(1+f1(p)) ≤ exp(Σ_{p>10⁶}3.00001/(p−1)²) ≤ exp(3.00001/(10⁶−1)) (f1(p) = (3 + 1/(p(p−1)))/(p−1)²). CORRECT.
Tails. d_0 > D_0 (any K): absolute bounds of §6–§8 apply, total ≤ (Σ_{d_0>D_0, P_0-smooth}d_0²/φ(d_0)³)·O(log ξ); the full P_0-smooth series ∏_{P_0}(1 + p²/(p−1)³ + …) converges, so the
  tail → 0 as D_0 → ∞, uniformly in T. Errors of (iv): absolutely summed, relative o(1), unchanged by the split.
Order of limits. For fixed P_0, D_0, K: limsup_{T→∞}(C3- and C1-quantities) ≤ (bound with P1eff) + c(P_0)·[tail(D_0) + max_{d_0≤D_0}η_K(d_0)] (the ω- and H log K terms → 0 as T → ∞).
  The left side does not depend on K, D_0: let K → ∞, then D_0 → ∞. All constants (O_{qP}(1), c(P_0)) are fixed before T → ∞. CORRECT.
Falsification attempts (none succeeded): exact A_d/A_1 at ξ = 10⁴, 10⁵, 10⁶, q = 1 and 13 (deltasums.log (b)): d = 2: −1.28, −1.24, −1.21 → −1; d = 4, 8, 9, 12: |·| ≤ 0.017 → 0; d = 3, 5, 6, 10, 30
  approach −1/2, −1/4, +1/2, +1/4, −1/8 slowly (relative deviation ≈ c·log d/√log ξ; at ξ = 10⁶ still 20–60%). The limit is forced: G_{qd}(ξ/δ)/G_{qd}(ξ) → 1 for each fixed δ and
  dominated convergence. Aggregate: exact Σ1·log ξ (sigma1.log; transform validated against brute force) = 0.818 (θ = 0.15, ξ = 2000, q = 1; 0.880 for q = 13) against the
  §6+§10 bound 1.945 (P1eff) and 22.4 (P1): the truth is 42% of the new bound, so the Euler-factor removal is not contradicted; the old P1 bound overshoots 27-fold.
  Caveat on rate: the asymptotics are reached only for huge ξ; this is irrelevant for κ_F (a liminf) but means no finite-T check of the constants is possible.


## §9 Assembly and the negative terms
C1(A) ≤ pref[(La+Sa)A + (Lb+Sb)], pref = 8P1eff/(π²ϰ(1−θ)³): small v gives pref∫_0^1(A − X/ϰ)w̄(X)dX for ANY w̄ ≥ W/2 pointwise; since A − X/ϰ > 0 on [0,1] (A > 1/ϰ),
replacing W/2 by w̄ is monotone, and expanding gives Sa = ∫w̄, Sb = −∫Xw̄/ϰ with the SAME w̄ (code: w used in both, l.316–317). Large v: (Aϰ − 1)aS1 with Aϰ − 1 > 0 and aS1 an upper
bound, so −aS1 in Lb is legitimate. κ_F's bound 2π(1/A − N(√C1+√C2)²/A³) is decreasing in C1, C2, so upper bounds suffice; Kb ≤ 0 (asserted) gives C1 ≤ KaA and the large-N
statement N·κ_F ≥ 2π[1/α − (√(Kaα) + √(C2/N))²/α³] (A = αN), increasing in N. CORRECT.

## Code audit (selbergn4.py Part 2 against §4–§11; rerun on a copy in /tmp: identical cert.log up to timings, identical result.json; 262 s, peak 95 MB; logs/cert_rerun.log)
- Eser/IEser: 40-term series, tail ≤ ymax⁴⁰/40!·2/(40² − 1/4) added to both Re and Im (ratio ≤ 3/41 < 1/2 justifies the factor 2). psi_cf/D2_cf: closed forms (4.1). selfcheck: series = closed form
  (< 1e−30), = defining double integral (1e−8), F4 formula = Im(e^{−iXu}W)/X (< 1e−30).
- R0max = max(box sups on [0,50], √π + 4/√50) = 2.338139 (the tail bound dominates; mpmath sup on (0,50] is 2.042); DA, Emax, DAp, Psimax: box sups + valid tails; mpmath values below the
  certified ones (logs/counting.log (b)). W(s) = ∫_s^∞|ψ|Y^{−2}: backward sums of box sups × ∫Y^{−2}, |ψ| ≤ 2√Y near 0; Ipsi = sup_l l·W(l) via l ≤ g_{i+1}, W(l) ≤ W(g_i). Upward.
- g_funcs/improved_consts2: S_X, F4 sups over 128 X-boxes × (40 + 600) u-boxes (ball arithmetic in X, u, ℓ1, ℓ2, √ℓ independently: over-estimates only). phi_vec: Φ(t_m) = max_i f_i·max{g_j : u_j^+ ≥ t_m − u_i^+}
  (suffix max over boxes sorted by upper end) = max{f_ig_j : u_i^+ + u_j^+ ≥ t_m}, ×(1+1e−12). Left-end values Φ(t_m), r(t_m) on cells (both non-increasing); (6t+4) at right ends. Upward.
- Small v: |sin y|/X ≤ min(1/(2ϰ) + t_{m+1}, ms/X_k), ms = 1 if a peak (n+½)π lies in [y_lo, y_hi], else the larger endpoint value (|sin| has no interior max between peaks). The O(v) shifts of y
  (v·½log(q/2π), v·O(1)) are o(1) and not in the code (asymptotic, correct). Peak detection and Gsin's n use float midpoints: non-rigorous only if y_lo or y_hi is within ~1e−30 of (n+½)π or nπ
  (measure zero; would change a term by O(radius²)). Not a flaw of substance.
- Large v: S1 at x = 1; Dm: u < θ max(DB, S1) with DB over 8 u-boxes, x ∈ [1,10] boxes and the x ≥ 10 tail; u ≥ θ √(1−u_a)·DAp; Dp: u < θ direct sum to x = 200 + 2Psimax/200, u ≥ θ Ipsi.
  aS1, aX, a2 = ½∫r·(…) as in §8; La = ϰ(aS1+aX), Lb = −aS1 + ϰa2. C3 = C5'⁴P1x(ϰi0 + ϰ²i1) = (16P1x/(π²(1−θ)²))((1/ϰ)i0 + i1) ✓; C2 = 4(C3−1)/(θϰ)² ✓; pref ✓.
- κ_F = lower end of the Arb ball of 2π(1/A − N(√(KaA+Kb) + √C2)²/A³); (θ, A) from a float search (admissibility only). Kb ≤ 0 and κ_F(4), κ_F(10) > 10⁵ × Dousselin asserted.
- Independent check of C3's Φ-integrals: exact Φ_BB gives i0 = 0.86667, i1 = 6.20167 at θ = 0.15 (cert 0.87381, 6.29385 ≥, as upper step sums must be); C3 with exact Φ ≈ 29.47 (cert 29.80734).
Verdict on the code: it computes exactly the formulas of DERIVATION §5–§11 with upward rounding. No discrepancy.

## Numerical falsification attempts (all logs in logs/)
1. Counting (counting.log (a)): toy signals, H = 0.3…3: m ≥ meas(G)/H in every case; pairing gives up to 1.8× fewer (33 vs 59 at H = 3). Consistent with (i).
2. Local δ-sums and A_d (deltasums.log): see (viii); no contradiction; slow convergence quantified.
3. Σ1 exactly (sigma1.log): 0.818/log ξ (θ = 0.15) and 0.853/log ξ (θ = 0.05) at ξ = 2000, slowly increasing (0.787 at ξ = 100); bounds 1.945/log ξ and 1.317/log ξ (P1eff).
4. IBP (ibpcheck.log): J(B) = ∫_1^B(B−x)e^{iρx}d(x;u1)d(x;u2)dx with exact d, B ≤ 1290, ρ ∈ {2,3,6}, six (u1,u2): |J| ≤ 0.454 × IBP bound; |J|/B → constant (linear in A);
   the oscillation-free bound grows like B² (8.98e5 at B = 1290 against |J| ≈ 1.35e3).


## Reviewer on (viii) (logs/review_viii.log) and the closure of its points
Verdict (reviewer, run timed out after a full recap): "GAP: local cancellation is correct, but the stated continuity and remainder estimates do not justify the summed asymptotics":
(R3) a common modulus is needed for the NORMALISED profiles √L·b_T (L = log ξ), uniformly in T and parameters; (R4) the o(·) in the per-d bound must be summable against the weights
(Π_{qd}(0)², Lemma-12 factors ∏(1+p^{−3/4}) are unbounded); (R6) the tail must read L·Σ1,tail ≪ tail(D_0) + o(1); (R7) a T-, v-uniform modulus for D/L on x ∈ [1,Aϰ]. Closure:
(R3) The normalised main-term profiles are explicit T-INDEPENDENT functions: √L·b(δ)/τ(δ) = (2/√π)(√ℓ1 − √ℓ2)/(1−θ), √L·δ^{−iv}s(δ,v)/τ(δ) = −e^{−iXu}W_X(u)/((1−θ)2√π),
  (1/√L)·Im(δ^{−iv}s)/(vτ) = −F4_X(u)/((1−θ)2√π), with ℓ1 = 1−u, ℓ2 = (θ−u)⁺; each is continuous on [0,1]×[0,∞) in (X,u), zero for u ≥ 1 (W_X(1) = 0), hence uniformly continuous:
  one modulus ω serves all T, all d, and X ∈ [0,1]. The reviewer's counterexample L^{−1/2}cos(Lu) cannot occur. (R7) D/L = d(x;u) = [Δ_2(xℓ1) − Δ_2(xℓ2)]/x is T-independent and
  continuous on the compact [1,Aϰ]×[0,1] (x = vL ∈ [1,Aϰ] on the large-v range): uniform modulus ω_d.
(R4) Domination: write each factor as main + error, main = Π·(profile)/√L·τ, error E with |E| ≤ τ_{1/2}(δ)a(d)(log T)^{−2/3}, a(d) = ∏_{p|d}(1+p^{−3/4}) (eqintero). The freezing error of Lemma S is
  ω(2log K/L) times the OLD absolute majorant of the same term (Π_{qd}(0)²·sup²/(Lφ(d))), so after squaring and the d-sum it is ≤ ω·(old bound for Σ1, which is O(1/L) by eqprodeul):
  o(1/L). The Lemma-12 errors are exactly Dousselin's (l.1650–1658: per-d ≪ (log T)^{−7/6}a(d)²/d for A'_d; similarly for A_d), summed with Euler products Σ_d a(d)^k/d^{1+o(1)} weights
  that converge after the log ξ from Σ 1/d: they contribute o(1/L) (Dousselin's own error analysis, unchanged by §10 because errors are never summed with signs).
(R6) Normalisation: Σ1 carries the prefactor (q/φ(q))²·16/(π²(1−θ)²L²); the d_0 > D_0 part of the weighted d-sum is ≤ tail(D_0)·O(L) with tail(D_0) := Σ_{d_0>D_0,P_0-smooth}d_0²/φ(d_0)³,
  so L·Σ1,tail ≪ tail(D_0) + o(1), with an implied constant independent of T and D_0 (the old absolute bound); same for the §7/§8 sums (normalised by T/log²T). This is what the
  order of limits uses. (R5, the reviewer's own remark) The remaining product is over p ∤ q, bounded by P1eff. With (R3), (R4), (R6), (R7) the reviewer's conditions are met.


## Addenda
5. Exact finite sums vs the main-term profile of §8 (exactsums.log; q = 1, d = δ = 1, θ = 0.15, inside Lemma 12's range X ≤ √log ξ): |d_exact − d_main|/|d_main| ≤ 0.116 (median 0.083) at
   ξ = 10⁴ and ≤ 0.086 (median 0.058) at ξ = 10⁶, decreasing; J over [1, √log ξ] agrees in size (e.g. ρ = 6: −1.42+3.31i exact vs −1.65+3.05i main at ξ = 10⁶). Consistent with Lemma 12.
6. (iii) remark: Dousselin applies Gabriel's convexity theorem directly to ∫_T^{2T} (l.1389); if that step is made rigorous with an entire weight (e.g. w(t) = e^{−(t−t_0)²/Δ²} = e^{(s−1/2−it_0)²/Δ²}),
   the first-moment argument applies verbatim to ∫w(Lη² − 1) (same rectangle; w is bounded on σ ∈ [1/2,3/2]), so the replacement Q = (C3 − 1 + o(1))·(length) survives any such smoothing.
7. Characters: after Bauer's mean value theorem only χ0 mod q survives (charprin, l.1413: χ(h)χ̄(k)χ̄(h/(h,k))χ(k/(h,k)) = χ0(hk)); all δ- and d-sums are character-independent. The q = 13 runs
   (χ0 mod 13) therefore cover any choice of even primitive characters mod 13 (six even, five non-principal).
Summary: every new step (i)–(viii) and §9 holds as stated, given Dousselin's analytic inputs; DERIVATION §10 needed (and now has) explicit uniformity/domination statements and the
4-fold extension for §8; no constant changes.


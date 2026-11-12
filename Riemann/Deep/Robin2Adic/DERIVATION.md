# DERIVATION — Gracia Llorente (Zenodo 24044336, v1, 7 Nov 2026), re-derived

Equation numbers = GL's. "✓" = re-derived here; "Δ" = we change it; "!" = issue.

## Prop 2 (exact identity) ✓
Setup (2): N_k ≤ n < N_{k+1}, x = p_k, T = ϑ(x) = log N_k, U = log(n/N_k) ≥ 0; n = ∏_{j≤m} q_j^{a_j}, q_1<…<q_m.
- m ≤ k: n ≥ q_1⋯q_m ≥ p_1⋯p_m = N_m and n < N_{k+1} ⇒ m < k+1. q_j ≥ p_j (j-th smallest prime factor ≥ j-th prime).
- σ(n)/n = ∏_j (1 − q_j^{−a_j−1})/(1 − q_j^{−1}) ⇒ log(σ(n)/n) = Σ_j (F(q_j) − D_{a_j}(q_j)).
- log(e^γ n log log n) − log n = γ + log log(T+U) (log n = T+U).
- E(x) = Σ_{i≤k}F(p_i) − γ − log log T  (NB: "log ϑ(x)" inside E is log of ϑ, E subtracts log log ϑ(x)).
- H(U) = log log(T+U) − log log T. Then E − B − D − H = Σ_j F(q_j) − D − γ − log log(T+U) = log(σ(n)/(e^γ n log log n)). ✓
- Robin at n ⇔ LHS < 0 ⇔ B + D + H > E. ✓ (n > 5040 ⇒ log log n > 0, T > 1 since k ≥ 4.)
- Signs: D ≥ 0 trivially; H ≥ 0 as U ≥ 0; B ≥ 0 as F decreasing, m ≤ k, q_j ≥ p_j.

## H(U) ≥ cU (after (3)) ✓ with our x
H'(u) = 1/((T+u) log(T+u)), decreasing in u; U < log p_{k+1} (n < N_{k+1} = N_k p_{k+1}). Need (T+U)log(T+U) ≤ xL/0.998, c = 0.998/(xL).
GL: T+U < 1.0011x (Chebyshev). For x ≥ 1.4e47: T = ϑ(x) ≤ x(1+2.8e-12) [BKLNW, |d|<d*], U < log(2x) ⇒ T+U < x(1+3e-12),
log(T+U) < L + 3e-12 ⇒ (T+U)log(T+U) < xL(1+1e-11) < xL/0.998. ✓ (c could be 0.99999xL^{-1}; irrelevant: c only enters Lemma 3's nonnegativity.)

## Lemma 3 (slots), (4) ✓
- B + D + H ≥ B + D + cU, U = Σ_j a_j log q_j − Σ_{i≤k} log p_i.
- F(q) − D_a(q) = log(1+q^{-1}+…+q^{-a}) ⇒ B + D + cU = Σ_{i≤k}(F(p_i) − c log p_i) − Σ_{j≤m} V_{a_j,c}(q_j), V_{a,c}(z) = log(1+z^{-1}+…+z^{-a}) − ca log z.
- V_{a,c} strictly decreasing on z>1 (both parts), q_j ≥ p_j ⇒ V(q_j) ≤ V(p_j); regroup slot j ≤ m: F(p_j) − c log p_j − V(p_j) = D_{a_j}(p_j) + c(a_j−1) log p_j. ⇒ (4). ✓
- Nonneg.: a_j ≥ 1; for p ≤ x, F(p) > 1/p ≥ 1/x > 0.998/x = cL ≥ c log p. ✓
- n even ⇒ q_1 = p_1 = 2, slot 1 ≥ D_{a_1}(2) = −log(1−2^{−a_1−1}) > 2^{−a_1−1}. ✓
- Odd n: slot 1 is D_{a_1}(2) with a_1 = v_{q_1}(n) — not 2^{−1}-sized; GL's Thm covers odd n via v_2 = 0? (! see below: odd n need Choie et al. [11]: Robin for odd n > 9.)
- Slot form (our corollary): counterexample ⇒ E(x) ≥ Σ_{j≤m} D_{a_j}(p_j) > Σ_j p_j^{−a_j−1}; so p_j^{a_j+1} > 1/E* for every j ≤ m, a_j = exponent of the j-th smallest prime factor.
- Conclusion rule: if E(x) ≤ E* ≤ 2^{−K−1} for all x in range, every even n with v_2(n) ≤ K satisfies Robin. K = ⌊−log_2 E*⌋ − 1 (= ⌈−log_2 E*⌉ − 2 unless integer).

## (5) integrated Mertens identity ✓
- Partial summation with d/dy(1/(y log y)) = −w(y), ∫ y w = log log y − 1/log y: Σ_{p≤x}1/p − log log x − B_1 = d/L + ∫_x^∞ (y−ϑ(y))w(y)dy. ✓
- Σ_{p≤x}F = Σ_{p≤x}1/p + Σ_{p≤x}h, γ = B_1 + Σ_p h(p) (Mertens), log log ϑ(x) = log L + log(1 + log(1+d)/L).
- ⇒ E(x) = ∫_x^∞(y−ϑ)w − Σ_{p>x} h(p) + G(d,L), G = d/L − log(1 + log(1+d)/L). ✓ h(p) > 0 ⇒ dropping −Σ_{p>x}h is an upper bound. ✓ (dropped series sign OK)
- G(0)=G_d(0)=0, G_dd = (L+1+log(1+d))/((1+d)^2(L+log(1+d))^2) ✓ ⇒ 0 ≤ G ≤ (d^2/2)·max G_dd ≤ 0.52 d^2/L for L>31, |d| small. (9) ✓.


## (10)-(16) explicit formula and kernel ✓
- (10) Ingham Thm 28 (integrated ψ_1 formula, absolutely convergent); A' = y − ψ = Σ y^ρ/ρ + log 2π + ½log(1−y^{-2}). With y = e^u: w(y)dy = e^{-u}f(u)du, f = 1/u + 1/u², ∫_x^∞ w = 1/(xL). ⇒ (11) ✓; termwise limit R→∞ by dominated convergence (majorant (16) × |ρ|^{-1}, Σ 1/(|ρ||1−ρ|) < ∞). ½∫log(1−y^{-2})w < 0 dropped ✓.
- (12) Σ_ρ 1/(ρ(1−ρ)) = 2Σ_ρ 1/ρ = 2 + γ − log 4π = 0.0461914; Re(1/(ρ(1−ρ))) = (β(1−β)+λ²)/|ρ(1−ρ)|² > 0 for every zero, and = 1/(|ρ||1−ρ|) on β = ½ ⇒ Σ_{|λ|≤H0} 1/(|ρ||1−ρ|) ≤ 0.0461914 ✓.
- (13) GL: Bellotti–Wong. Ours Δ: HSW two-sided |N − M − 7/8| ≤ R; Σ_{|λ|>H0} λ^{-2} = 4∫_{H0}^∞(N(t)−N(H0))t^{-3}dt ≤ (log(H0/2π)+1)/(πH0) + 8((a+b)(h/2+¼) + c/2)/H0² = τ = 2.9594e-12. Only enters at s ≈ ½ (negligible, ~1e-35).
- (14) GL ✓ (Schoenfeld/McCurley). Ours Δ: ψ−ϑ = Σ_{k≥2}ϑ(y^{1/k}) ≤ (1+1.1199e-9)√y + 1.01624 y^{1/3} log y/log 2 (BKLNW row 50 for √y > e^50; RS62 Thm 9 for k ≥ 3) = c_ψ(y)√y, c_ψ decreasing; ∫_x^∞ √y w ≤ 2f(L)e^{−L/2}.
- (16) two integrations by parts ✓: |K_ρ| ≤ x^{β−1}[f(L)/|1−ρ| + 2|f'(L)|/|1−ρ|²], |f'|/f = (L+2)/(L(L+1)) ✓. ε = 2(L+2)/(L(L+1)H0) < 1e-12 for |λ| > H0, x ≥ 1e19 ✓. Low zeros: |1−ρ| ≥ 14.13 ⇒ factor 1 + 2(L+2)/(14.13 L(L+1)).

## Lemma 4, Lemma 5, (17)-(20) ✓
- W(s) := Σ_{|λ|>H0, β>s} λ^{-2} = 4∫_{H0}^∞ N(s,t)t^{-3}dt for s > ½ (RH to H0 ⇒ N(s,H0) = 0) ✓. GL's (log t)^p ≤ 28^{p−3}(log t)³ valid (log H0 = 28.73 > 28, p ≤ 3) ✓; I3, I2, I1 closed forms ✓. W ≤ τ/2 by ρ ↔ 1−ρ̄ ✓.
- Lemma 5 ✓: Σ x^{β−1}/λ² = Σ g_{s(ρ)}(x)/λ², g_s = (x^{s−1}+x^{−s})/2 (zero set invariant under β ↦ 1−β, g symmetric), s(ρ) = max(β,1−β); right-endpoint values + Abel summation with running-min W̃ ⇒ c_j ≥ 0, Σc_j = τ; T_J decreasing in x ✓. GL's code takes s_j for j < J, s_J = 7/8 carries c_J = 2W_{J−1} ✓.
- (20) ✓ (2.134 ≥ 2·1.02 + 0.047(1+ε')). Reproduction: our Arb run of (20) as printed (GL Table 1, J = 2048, GL's row index) at L = 106 gives 1.227034555e-23 (GL: < 1.227035e-23) ✓.

## Ours: inputs (all at x ≥ e^{108.57})
| input | statement used | source |
|---|---|---|
| range | Robin for 5041 ≤ n ≤ 10^(10^46.79) | RobinCH (ours; Chirre–Helfgott pointwise ψ, Morrill–Platt below 1e10) |
| half-plane | β ≤ B = 0.87338196661320752 (we use dyadic B_UP = B + 8.7e-14) | Bebel (Lean); Richey PR 129: β ≤ B_2 = 0.8733819666132075199 < B |
| RH height | β = ½ for 0 < |λ| ≤ 3e12 | Platt–Trudgian 2021 |
| zero count | HSW 0.1038, 0.2573, 9.3675 (two-sided) | Hasanalizade–Shen–Wong 2022 |
| density A | KLN Lemma 4.14 rows (Table 1 parameters) with RH(3e12), row 0.60 inadmissible | KLN 2018; Ingham4 build_kln_H0 |
| density B (cond.) | CS Cor. 1.1 rows with B_1 × (b_4''/0.54)^{(2s_2−1)/(2−s_2)}, b_4'' from 0.5501 | 2507.15184 v2 (unrefereed); Ingham4 CSRow |
| Chebyshev | |ϑ(x)−x| < 2.4530e-12 x (x > e^100), ϑ(t) < (1+1.1199e-9)t (t > e^50) | BKLNW Math. Comp. 2021, Table ThetaEpsilon (TeX l.5734, 5752) |
| crude | ϑ(x) < 1.01624x (x > 0); Bertrand | Rosser–Schoenfeld 1962 Thm 9 |
| odd n | Robin for odd n > 9 | Choie–Lichiardopol–Moree–Solé 2007 |


## Ours: threshold n > 10^(10^46.79) ⇒ x ≥ e^{108.57}
N_k ≤ n < N_{k+1}, x = p_k. log n < log N_{k+1} = ϑ(p_{k+1}) = ϑ(x) + log p_{k+1}, and log n > Λ := 10^{46.79} log 10 = e^{108.57199}.
- x ≤ e^100: ϑ(x) + log p_{k+1} ≤ 1.01624e^100 + log(2e^100) < Λ (RS, Bertrand): impossible.
- x > e^100, x < e^{108.57}: ϑ(x) + log p_{k+1} < (1+2.4530e-12)e^{108.57} + 108.57 + log 2 < Λ: impossible. (cert.py checks both in Arb.)

## Ours: the majorant (cert.py majorant(L)), x = e^L ≥ e^{108.57}
E(x) ≤ f e^{−L/2}[2c_ψ(x) + 0.0461915(1 + 2(L+2)/(14.13L(L+1)))] + f(1 + 2(L+2)/(L(L+1)H0))·T(x) + log(2π)/(xL) + (d*²/2) sup_{|δ|≤d*} G_dd(δ,L),
d* = 2.4530e-12, T(x) = Σ_{j=0}^{J} c_j g_{s_j}(x):
- grid: exact dyadics s_0 = ½ < … < s_J = B_UP, step 2^{-9} on [½, 0.8), 2^{-16} on [0.8, B_UP] (J+1 = 4964);
- W̃_0 = τ/2, W̃_j = min(W̃_{j−1}, W̄(s_j)) (exact upper roundings), c_0 = τ − 2W̃_0, c_j = 2(W̃_{j−1} − W̃_j), c_J = 2W̃_{J−1}; all c_j ≥ 0 asserted.
- Validity (Lemma 5 argument with the grid from ½ to B): zeros with s(ρ) ∈ (s_{j−1}, s_j] get g_{s_j}; weight with s(ρ) > s_j ≤ 2W(s_j) ≤ 2W̃_j; none with s(ρ) > B_UP.
- W̄(s) = 4∫_{H0}^∞ N̄(s,t)t^{-3}dt in closed form, Ig(q,e) := ∫_{H0}^∞ t^{−q−1}(log t)^e dt = Γ(e+1, qh)/q^{e+1} (Arb gamma_upper):
  KLN: N̄ = A(s)(log kt)^{2s}(log t)^{5−4s}t^{8(1−s)/3} + max(B(s),0)(log t)² ≤ A(log t)^{5−2s}t^{8(1−s)/3} + … (1e9/H0 ≤ k ≤ 1) ⇒ W̄ = 4[A·Ig(2 − 8(1−s)/3, 5−2s) + B·I2], min over admissible rows with s > ½ + d/log H0.
  CS: W̄ = 4[B_1'·Ig(2−α(s), p(s)) + B_2 I2 + B_3 I1], min over rows with s_1 ≤ s ≤ s_2 (exact incomplete Γ instead of GL's 28^{p−3}).
- Monotone: f, e^{−L/2}, c_ψ (e^{−L/6}L, L > 6), both ε's, g_s(e^L) (0<s<1, c_j fixed), 1/(xL), G_dd ((L+a+2)/(L+a)³ > 0) all decrease in L ⇒ E* := sup_{x ≥ e^{108.57}} E(x) ≤ majorant(108.57). Sample check at 8 values of L in cert.py.
- K := ⌊−log_2 E*⌋ − 1, asserted 2^{−K−1} > E*. Slot form: p_j^{a_j+1} > 1/E* for every j ≤ m.


## Review (reviewer subagent, 12 Nov, verdict-only on Steps Prop 2, Lemma 3, threshold, majorant; text inline)
No ERROR found. Points raised and answers:
- Prop 2 needs T > 1 (k ≥ 2): n > 5040 ≥ N_4 ⇒ k ≥ 4, T ≥ ϑ(7) > 1. ✓
- Lemma 3's c-condition was verified only for x > e^100: we use it only for x ≥ e^{108.57} (cert.py checks (T+U)log(T+U) < xL/0.998 at e^{108.57}; ratio 1 + 2.5e-12, decreasing in x). ✓
- log Λ = 108.571989946… (displayed value rounded). ✓
- β ≤ B: Bebel's Lean theorem (ζ ≠ 0 for Re s > 0.87338196661320752) and ζ(s) = 0 ⇔ ζ(1−s̄) = 0 ⇒ 1 − B ≤ β ≤ B, so s(ρ) ≤ B ≤ B_UP. ✓
- d two-sided: BKLNW table gives (1 − m_0)x < ϑ(x) < (1 + M_0)x with m_0 = M_0 = 2.4530e-12 for x > e^100. ✓
- W̃_0 ≤ τ/2: W̃_0 := τ/2 exactly, so c_0 = 0; all c_j ≥ 0 asserted in cert.py. ✓
- Density/grid data, τ, c_ψ: in "Ours" above and cert.py; checked independently by check.py (cert/check = 1.0017 CS, 1.0029 KLN; closed-form W vs quadrature agree to 1e-9, KLN closed form ≥ quadrature by design).
- K arithmetic: 2^{-73} = 1.0588e-22 > 8.75e-23 (K = 72); 2^{-76} = 1.3235e-23 > 7.67e-24 (K = 75). ✓


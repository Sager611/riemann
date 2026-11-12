# HCAudit — complete proof of Theorem 1 (hard-core augmentation of Lamzouri's Thm Main2)

Source: Y. Lamzouri, arXiv 2611.14486v1, SmallGapsPaper.tex (574 lines); "l.n" = its line numbers. f^(u) = ∫ f(α)e^{−2πiαu}dα.
J_θ(f) := f(0) + 2∫_0^θ f(α)α dα. L := log T. Unfolding x = γL/2π (fixed L, as in Lamzouri's E_λ).

## 0. Lamzouri's definitions (verbatim content)
- l.178: N(T) = #{n: 0 < γ_n ≤ T} = (T/2π)log(T/2π) − T/2π + O(log T) (S(T) = O(log T)); ordinates counted with multiplicity.
- l.183: D(λ) = liminf_{T→∞} N(T)^{−1}#{n: 0 < γ_n ≤ T, γ_{n+1} − γ_n ≤ 2πλ/L} (fixed L = log T, not log γ_n).
- l.292: w(u) = 4/(4 − u²); for on-line ρ, ρ': w(ρ − ρ') = 4/(4 + (γ − γ')²) ∈ (0, 1].
- l.298: Z(T) = multiset of nontrivial zeros with 0 < γ ≤ T.
- l.331–339 (Prop SmallSetZeros): T large, E ⊂ Z(T) invariant under ρ → 1 − ρ̄, #E = δN(T), 0 < δ < 1; 1/2 < θ' < 4/7;
  g ∈ C_c² even, g ≥ 0, supp g ⊂ [−θ', θ']. Then Σ_{ρ,ρ'∈Z∖E} g^(i(ρ−ρ')L/2π)w(ρ−ρ') = (g(0) + 2∫_0^1 αg + O_g(δ^{1/8−7θ'/32} + L^{−1/2}))N(T).
- l.405: E_λ(T) = {ρ_n: 0 < γ_n ≤ T, n ≥ 2, min(γ_{n+1} − γ_n, γ_n − γ_{n−1}) ≤ 2πλ/L}; l.411: Z_0(T) = Z(T) ∖ E_λ(T).
- External inputs: Prop SmallSetZeros rests on Conrey's density estimate (slope 8/7, Baluyot's thesis) and BGST Lemma 5;
  audited PASSED line by line by lamzaudit (Riemann/Deep/LamzAudit/REPORT.md). Not re-derived here.

## 1. Reduction to Z_0 and the lower bound (Lamzouri l.402–441, checked for our f and for the truncated F)
Claim L. Let 1/2 < θ < 4/7, f even, f ≥ 0, f ∈ L¹(R) ∩ C^∞ (Schwartz suffices) or f ∈ C_c² (the truncated F). Let δ ∈ (0,1)
and T large with #E_λ(T) < δN(T). Then, with κ := 1/16 − 7θ/64 > 0,
  S(f) := Σ_{ρ,ρ'∈Z_0(T)} f^(i(ρ−ρ')L/2π) w(ρ−ρ') ≥ (J_θ(f) − C_f δ^κ − C_f L^{−1/2}) N(T),  C_f independent of δ, T.
Proof. (a) Z_0 is on the line, simple, separated (l.443–444). If β_n ≠ 1/2, 1 − ρ̄_n = 1 − β_n + iγ_n is a distinct zero with the same
 ordinate, so in the multiplicity-counted ordering two consecutive ordinates coincide: gap 0 ≤ 2πλ/L, both lie in E_λ (their
 indices are ≥ 2: γ_1 = 14.13… is simple and on the line, γ_2 − γ_1 = 6.89 > 2πλ/L for L > 1). A multiple zero likewise gives
 gap 0. Two distinct ρ, ρ' ∈ Z_0 have |γ − γ'| > 2πλ/L: if they are consecutive zeros this is the definition (the later one has
 index ≥ 2 and is not in E_λ); if zeros lie between, |γ − γ'| ≥ γ_{a+1} − γ_a for the earlier ρ_a, which is > 2πλ/L when a ≥ 2
 (ρ_a ∉ E_λ) and = 6.89… when a = 1.
(b) E_λ is invariant under ρ → 1 − ρ̄: on-line zeros are fixed; off-line ones are paired with their mirror, which is in E_λ by (a).
 So Prop SmallSetZeros applies to E = E_λ(T) with δ' := #E_λ/N(T) < δ (if E_λ = ∅ the same proof gives the δ' = 0 case; the
 error term is increasing in δ').
(c) Cut-off (l.413–424): ε := (4/7 − θ)/2, θ' := θ + ε = (θ + 4/7)/2 ∈ (1/2, 4/7); ψ_ε ∈ C_c^∞((−θ', θ')) even, 0 ≤ ψ_ε ≤ 1,
 ψ_ε = 1 on [−θ, θ]; g := ψ_ε f. Hypotheses of Prop SmallSetZeros for g: even (both factors even) ✓; g ≥ 0 ✓; g ∈ C_c² (f ∈ C² on
 supp ψ_ε, which is all Prop SmallSetZeros needs; f Schwartz gives C_c^∞, F ∈ C_c² gives C_c²) ✓; supp g ⊂ [−θ', θ'] ✓.
 Conclusion (l.425–429): S(g) = (g(0) + 2∫_0^{θ'} αg + O_g(δ'^{1/8−7θ'/32} + L^{−1/2}))N(T), and 1/8 − 7θ'/32 = 1/16 − 7θ/64 = κ.
 Since g = f on [−θ, θ] and g ≥ 0: g(0) + 2∫_0^{θ'} αg ≥ J_θ(f).
(d) Positivity (l.430–436). For a finite multiset B of on-line zeros and h ∈ L¹ even (so h^ is real, even, bounded):
 S_B(h) := Σ_{ρ,ρ'∈B} h^((γ−γ')L/2π) · 4/(4+(γ−γ')²) = (2/π)∫_R∫_R h(α)|Σ_{ρ∈B} T^{iαγ}/(1+(γ−t)²)|² dα dt.
 Proof of the identity: expand |Σ|²; ∫_R dt/((1+(γ−t)²)(1+(γ'−t)²)) = 2π/(4+(γ−γ')²) (convolution of two Cauchy kernels;
 at γ = γ' it is ∫(1+t²)^{−2} = π/2 ✓), and ∫ h(α)T^{iα(γ−γ')}dα = h^(−(γ−γ')L/2π) = h^((γ−γ')L/2π); Fubini: the integrand is
 ≤ |h(α)|·(Σ_B (1+(γ−t)²)^{−1})², integrable. No support condition on h is used. Apply to h = f − g ≥ 0 (f ≥ g since
 0 ≤ ψ_ε ≤ 1, f ≥ 0): S(f) − S(g) = S_{Z_0}(f − g) ≥ 0. (Z_0 consists of on-line zeros by (a), so f^(i(ρ−ρ')L/2π) = f^((γ−γ')L/2π).)
(e) Combine: S(f) ≥ S(g) ≥ (J_θ(f) − C_g δ^κ − C_g L^{−1/2})N(T) = (Eq:ReductionGoodSet2) (l.438–441). □
Remark. Claim L needs neither compact support of f nor f^ ≤ 0 anywhere: Lamzouri's (Eq:ReductionGoodSet2) is a pure lower bound.
So the Schwartz f of the certificate is admissible directly; the truncation lemma (§3) is a redundant safety net.
Addendum to §1 (hidden hypotheses of Prop SmallSetZeros, l.340–390, read for this audit): its proof uses (Eq:GeneralIdentityPC)
(f even, supp ⊂ [−θ', θ'] so that T^{α(β−β')} is bounded on off-line zeros), Cauchy–Schwarz on (2/π)∫∫ f|S_Z||S_E| (f ≥ 0),
(Eq:UpperBoundKernel) |f^(i(ρ−ρ')L/2π)| ≪_f T^{θ'|β−β'|}/(1 + ((γ−γ')L)²) (two integrations by parts: f ∈ C², supp ⊂ [−θ', θ']),
the grouping under ρ → 1 − ρ̄ (E invariant), and Hölder with Lemma ZeroDensity (θ' < 4(1 − ε)²/7): δ enters only as δ^{1/(pq)},
with constants depending on f, θ' only. Nothing else: no sign of f^, no smoothness beyond C². All hold for g = ψ_εf and g = ψ_εF.

## 2. Step (iii′): the hard-core upper bound (complete)
Setting. Under #E_λ(T) < δN(T), list Z_0(T) by ordinate: γ^(1) < … < γ^(N_0) (distinct, §1(a)), N_0 = N(T) − #E_λ(T) > (1 − δ)N(T).
x_i := γ^(i)L/2π, g_i := x_{i+1} − x_i > λ (i < N_0), 0 < x_1 < x_{N_0} ≤ TL/2π. Normalisation is the SAME fixed L = log T in
E_λ (l.405) and in the pair sum (l.425, l.438): gaps > λ hold exactly in the unfolding used by f^; the local density log(γ/2π)/2π
never enters (it would only matter for the mean gap, and only the total range enters, below). So no o(1) from varying log γ.
S(f) = Σ_{i,j=1}^{N_0} f^(x_i − x_j)w_ij, w_ij = 4/(4 + (γ^(i) − γ^(j))²) ∈ (0, 1]; f^ real and even (f real, even).
(1) Split by k = |i − j| (ordered pairs: each unordered pair at k ≥ 1 occurs twice):
 k = 0: N_0 terms, f^(0)w_ii = f^(0).
 k = 1: 2Σ_{i=1}^{N_0−1} f^(g_i)w_{i,i+1};  k = 2: 2Σ_{i=1}^{N_0−2} f^(g_i + g_{i+1})w_{i,i+2}.
 k ≥ 3: |x_i − x_j| = g_i + … + g_{j−1} > 3λ ≥ U, so f^(x_i − x_j) ≤ 0 and w_ij > 0: each term ≤ 0, dropped (upper bound).
(2) Weights for k ∈ {1, 2}: 0 ≤ 1 − w = d²/(4 + d²) ≤ d²/4, d = 2π|x_i − x_j|/L, so |f^(u)w − f^(u)| ≤ π²u²|f^(u)|/L² ≤ π²M_f/L²,
 M_f := sup_u u²|f^(u)| < ∞ (f Schwartz ⇒ f^ Schwartz; F ∈ C_c² ⇒ |F^(u)| ≤ ‖F''‖_1/(4π²u²)). At most 4N_0 such ordered pairs:
 S(f) ≤ N_0 f^(0) + 2Σ_{i=1}^{N_0−1} f^(g_i) + 2Σ_{i=1}^{N_0−2} f^(g_i + g_{i+1}) + 4π²M_f N_0/L².
(3) (H) at (x, y) = (g_i, g_{i+1}), i = 1, …, N_0 − 2 (both > λ): f^(g_i) + f^(g_i + g_{i+1}) ≤ c − A(g_i) + A(g_{i+1}) + β(g_i − 1).
 Sum: ≤ (N_0 − 2)c + [A(g_{N_0−1}) − A(g_1)] + β(x_{N_0−1} − x_1 − (N_0 − 2))   (A telescopes; Σ_{i≤N_0−2} g_i = x_{N_0−1} − x_1).
 The one unmatched k = 1 term is 2f^(g_{N_0−1}) ≤ 2 sup|f^|. Hence, with K_0 := 2sup|f^| + 4|c| + 4sup|A| + 4β (A bounded),
 S(f) ≤ N_0(f^(0) + 2c) + 2β(x_{N_0−1} − x_1 − N_0) + K_0 + 4π²M_f N_0/L².
(4) Mean-gap (β) term, sign: β ≥ 0 and x_{N_0−1} − x_1 < TL/2π. l.178: N(T) = (T/2π)(L − log 2π − 1) + O(L), so
 TL/2π − N(T) = (T/2π)(log 2π + 1) + O(L), so (TL/2π − N(T))/(N(T)/L) → log 2π + 1 = 2.838 < 4: ≤ 4N(T)/L for all large T
 (reviewer's fix: the O(L) remainder is unspecified, so no explicit cutoff; at T = 2·10⁴ the actual values are 9031.6 ≤ 4N/L = 9084;
 the ratio approaches 2.838 slowly, 3.97 at L = 9.9). Thus
 x_{N_0−1} − x_1 − N_0 ≤ (N(T) − N_0) + 4N(T)/L = #E_λ(T) + 4N(T)/L < (δ + 4/L)N(T), and 2β(…) ≤ 2β(δ + 4/L)N(T).
 (β ≥ 0 makes this one-sided bound suffice: the range exceeds N_0 by ≈ #E_λ ≥ 0. A β < 0 would need range ≥ N_0 − O(δN), true by
 Riemann–von Mangoldt locally near 0 and T, but not needed: the certificates have β = 0.6496 > 0.)
(5) N_0(f^(0) + 2c) ≤ N(T)(f^(0) + 2c) + #E_λ|f^(0) + 2c| ≤ N(T)(f^(0) + 2c) + δN(T)|f^(0) + 2c| (either sign of f^(0) + 2c).
 Result: S(f) ≤ N(T)[f^(0) + 2c + δ(2β + |f^(0) + 2c|) + 8β/L + 4π²M_f/L²] + K_0.
(6) Contradiction. Claim L: S(f) ≥ N(T)(J_θ(f) − C_fδ^κ − C_fL^{−1/2}). Put η := J_θ(f) − f^(0) − 2c > 0 and fix δ ∈ (0, 1) with
 C_fδ^κ + δ(2β + |f^(0) + 2c|) < η/2. For T ≥ T_0(δ, f) the remaining terms are < η/2, a contradiction. So #E_λ(T) ≥ δN(T)
 for all T ≥ T_0 (every large T, not a subsequence).
(7) From E_λ to D(λ). Each ρ_n ∈ E_λ(T) has a gap ≤ 2πλ/L at index n (γ_n ≤ T) or at index n − 1 (γ_{n−1} ≤ γ_n ≤ T); an index m
 with γ_m ≤ T and γ_{m+1} − γ_m ≤ 2πλ/L is charged by at most ρ_m and ρ_{m+1}. So #{m: γ_m ≤ T, γ_{m+1} − γ_m ≤ 2πλ/L} ≥
 #E_λ(T)/2 ≥ δN(T)/2 for T ≥ T_0, i.e. D(λ) ≥ δ/2 > 0 (l.183 is a liminf over all T). □
Multiplicity audit. Diagonal once (N_0 f^(0)); k = 1, 2 twice (both orders; f^ even), matched to the N_0 − 2 applications of (H)
(each covering one k = 1 and one k = 2 unordered pair, times 2 for the orders) plus the leftover k = 1 term; k ≥ 3 twice, all ≤ 0. Periodic sanity check: for a cyclic gap sequence of mean 1 the bound reads
(1/N)Σ_{i,j} f^ ≤ f^(0) + 2c exactly (A-sum and Σ(g_i − 1) vanish), the form tested numerically in §6.
What (iii′) uses about f: f^ real, even, bounded; f^ ≤ 0 on [U, ∞) with U ≤ 3λ; (H); M_f < ∞. Nothing about f's support.
Edge cases. U ≤ 3λ is used once (k ≥ 3 drop) and needs only gaps > λ (strict, exact here). λ < U is not used. N_0 ≥ 3 for large T.
θ enters only through Claim L. No condition on the sign of c (the reviewer's remark is confirmed).

## 3. Truncation lemma (Lemma T′): valid, exact, and redundant
Lemma T′. Let f satisfy Theorem 1's hypotheses, η := J_θ(f) − f^(0) − 2c > 0. Then some F ∈ C_c²(R), even, F ≥ 0, has F^ < f^ on
[λ, ∞) (so F^ < 0 on [U, ∞) and (H) holds for F with the same A, β, c, with no o(1)) and J_θ(F) − F^(0) − 2c ≥ η/2.
Proof. (a) Cut-off: h ∈ C_c^∞ even, 0 ≤ h ≤ 1, h = 1 on [−1, 1]; h_R := h(·/R), R ≥ 1; E_R := f(1 − h_R) ≥ 0, Schwartz,
 supported in |x| ≥ R. Leibniz: ‖E_R^{(k)}‖_1 ≤ C_k Σ_{j≤k} ∫_{|x|≥R}|f^{(j)}| → 0 (R → ∞, k ≤ 6). With ε_R := max(‖E_R‖_1,
 (2π)^{−6}‖E_R^{(6)}‖_1): |E_R^(ξ)| ≤ ε_R ≤ ε_Rξ^{−6} for 0 < ξ ≤ 1 and |E_R^(ξ)| ≤ ‖E_R^{(6)}‖_1(2πξ)^{−6} ≤ ε_Rξ^{−6} for ξ ≥ 1.
(b) Reserve (SmallGap Lemma T(a), re-checked here): R(t) = (sin πt/πt)²/(1 − t²), R^(u) = (1 − |u| + sin(2π|u|)/2π)1_{|u|≤1}
 ≥ 0, R^ ∈ C_c² (Lamzouri l.455–462; R^' = −1 + cos 2πu, R^'' = −2π sin 2πu vanish at u = 0+, 1; numerically re-verified, §5).
 G(x) := λR^(λx) + (λ/φ)R^(λx/φ), φ = (1+√5)/2: G even, ≥ 0, C_c², G^(ξ) = R(ξ/λ) + R(φξ/λ), G^(0) = 2. For y = ξ/λ ≥ 1:
 R(y), R(φy) ≤ 0 and |R(t)| ≥ sin²(πt)/(π²t⁴), so |G^| ≥ 4M²/(π²φ⁴y⁴), M := max(‖y‖, ‖φy‖). With n, m nearest y, φy (n ≥ 1):
 ‖nφ‖ ≤ |nφ − m| ≤ φ‖y‖ + ‖φy‖ ≤ φ²M; and ‖nφ‖ ≥ 1/(√5n + 1/2) ≥ 1/(3n) (m'² − m'n − n² is a nonzero integer for m' nearest nφ);
 n ≤ 3y/2. So M ≥ 2/(9φ²y) and G^(ξ) ≤ −κ(λ/ξ)⁶, κ = 16/(81π²φ⁸) ≈ 4.26·10⁻⁴, for all ξ ≥ λ.
(c) F := fh_R + c_rG = f − E_R + c_rG, c_r := η/4, R with ε_R < c_rκλ⁶. Even ✓, ≥ 0 ✓, C_c² ✓ (so ψ_εF ∈ C_c², M_F < ∞).
 For ξ ≥ λ: F^(ξ) − f^(ξ) = −E_R^(ξ) + c_rG^(ξ) ≤ (ε_R − c_rκλ⁶)ξ^{−6} < 0. Hence for x, y ≥ λ: F^(x) + F^(x+y) < f^(x) + f^(x+y)
 and (H) passes to F verbatim; F^ < f^ ≤ 0 on [U, ∞). J_θ(F) = J_θ(fh_R) + c_rJ_θ(G) ≥ J_θ(f) (h_R = 1 on [−θ, θ], G ≥ 0).
 F^(0) = f^(0) − ∫E_R + 2c_r ≤ f^(0) + 2c_r. So J_θ(F) − F^(0) − 2c ≥ η − 2c_r = η/2. □
Status. Correct as written in DERIVATION §4 (it states the same with "F^ ≤ f^ on [λ, ∞)"). Redundant: by Claim L (§1) and §2 the
Schwartz f itself is admissible — Lamzouri's lines use compact support only for g = ψ_εf (C_c^∞ when f is Schwartz), and the
identity of §1(d) holds for any even f ∈ L¹. Either route proves Theorem 1 for Schwartz f.

## 4. Certificate scheme: why the checker's three regions bound sup (H)
Φ(x, y) := f^(x) + f^(x+y) + A(x) − A(y) + β(1 − x) on Q = [λ, ∞)². A := piecewise-linear interpolation of a_0, …, a_77 at
t_k = λ + k·dg (t_77 = U), A := a_77 on [U, ∞). Lattice x_i = λ + iδ, δ = dg/m, N = 77m (m = 5). M2 ≥ sup_{[λ,2U]}|f^''|.
Q1 = [λ, U]²: on a cell C = [x_i, x_{i+1}] × [x_j, x_{j+1}], ℓ := A(x) − A(y) + β(1 − x) is affine in x plus affine in y (m | ·:
each lattice interval lies in one node interval), so ℓ equals its bilinear interpolant. For u := f^(x) + f^(x+y),
u − I_xI_yu = (u − I_xu) + I_x(u − I_yu), I_x a convex combination: |u − Iu| ≤ (δ²/8)(sup|u_xx| + sup|u_yy|) ≤ (δ²/8)(2M2 + M2)
(u_xx = f^''(x) + f^''(x+y), u_yy = f^''(x+y), arguments in [λ, 2U]). A bilinear function peaks at a corner, so
sup_{Q1}Φ ≤ max_{i,j}Φ(x_i, x_j) + 3M2δ²/8 =: R1 (corner values with Arb upper endpoints of f^).
Q2 = [λ, U] × [U, ∞): A(y) = a_77, x + y > U ⇒ f^(x+y) < 0: Φ < f^(x) + A(x) − a_77 + β(1 − x); 1D linear interpolation error
≤ M2δ²/8: sup_{Q2}Φ ≤ max_i[f^(x_i) + A(x_i) − a_77 + β(1 − x_i)] + M2δ²/8 =: R2.
Q3 = [U, ∞) × [λ, ∞): A(x) = a_77, f^(x), f^(x+y) < 0, β(1 − x) ≤ β(1 − U) (β ≥ 0): Φ < a_77 − min_k a_k + β(1 − U) =: R3
(inf_{[λ,∞)}A = min_k a_k for a piecewise-linear A with constant tail). Q1 ∪ Q2 ∪ Q3 = Q, including x, y → ∞ and all corners.
Hence (H) holds with c = c_cert := max(R1, R2, R3) + 10⁻¹². M2 := max over a grid of step η on [λ, 2U] of |f^''| + M3η/2,
M3 ≥ sup_R|f^'''|: f^(u) = sΣ_n c_n ψ_{2n}(√(2π)su), c_n = (−1)^nB_n(2^{2n}(2n)!√π)^{1/2} (ψ_k normalised Hermite functions),
|ψ_k| ≤ 1.086435π^{−1/4} = 0.81605 (Cramér), ψ_k' = √(k/2)ψ_{k−1} − √((k+1)/2)ψ_{k+1} ⇒ |ψ_k'''| ≤ 0.81605(2(k+3))^{3/2}:
M3 = s(2πs²)^{3/2}Σ|c_n|·0.8163·(2(2n+3))^{3/2}. Positivity: f(α) = p(2πα²/s²)e^{−πα²/s²}, f^(u) = s·q(2πs²u²)e^{−πs²u²}
(FT[H_k(√(2π)α)e^{−πα²}] = (−i)^kH_k(√(2π)ξ)e^{−πξ²}; dilation by s); p(0) > 0 & no root of p in (0, ∞) ⇒ f > 0;
q(w0) < 0 & no root in (w0, ∞), w0 = 2·3.14159265·s²U² < 2πs²U² ⇒ f^ < 0 on [U, ∞).
J_θ(f) = p(0) + (s²/π)Σ_k p_k2^k k!(1 − e^{−X}Σ_{j≤k}X^j/j!), X = πθ²/s² (substitute v = πα²/s²). f^(0) = s·q(0).


## 5. Certificate re-verification (cert768; all runs on copies in /tmp, memcap 290 MB; logs in logs/)
| check | method | result |
|---|---|---|
| rerun cert_aug.py (float corners) | as shipped, EPS = 1/1000 | identical: c_cert = −0.3478616532, margin [0.001462588425 ± 3.53e-13] PASS (rerun_cert768.log) |
| rerun FULLARB=1 | corners in Arb | corner max [−0.347888828296 ± 1.6e-13], PASS (rerun_cert768_fullarb.log) |
| rerun xcheck_aug.py | mpmath quadrature, direct FT at 7 points, 10⁶ random | identical; random/edge max −0.3479745 ≤ c_cert (rerun_xcheck768.log) |
| cert770, cert775 | as shipped | margins 0.005142294469, 0.01374683343, PASS (rerun_cert770/775.log) |
| f > 0, f^ < 0 on [U, ∞) | Arb root isolation (acb complex_roots; independent of Sturm) | p: one real root −0.000218; q: real roots 3.143, 30.709, 40.357846 < w0 = 40.358051 < 2πs²U² = 40.3580513 (hc_indep768.log) |
| M3 | numerical max|f^'''| on [0, 8] | 30.2 vs checker bound 3536.5: valid, loose (hc_indep768.log) |
| sup Φ on [λ, U]² | 8× finer lattice (dg/40) in Arb + local refinement (step dg/400) | −0.3478845702 at (λ, 1.11709); c_cert above it by 2.29e-5 (hc_indep768.log) |
| sup Φ on Q2, Q3 | sampled | −0.564021 (R2 bound −0.5640124), −0.348686 (R3 bound −0.3481091) |
| second route | M2 by Arb ball cover of [λ, 2U] (no Cramér): M2 ≤ 4326.5 (dependency-inflated, rigorous); corners on dg/40 lattice in Arb | c = −0.347481068, margin [0.0007014179351 ± 4.8e-14] PASS (hc_rig768.log) |
| J | rigorous acb.integral | [0.306839281978626 ± 1.8e-16] = closed form (hc_rig768.log) |
| R^ = FT of R | mpmath quadrature at x = 0.3, 1.5, 2.7 | equal to 12 digits (hc_indep768.log) |
Code audit of cert_aug.py: Hermite recurrence, P/Q/p/q, J formula, f^ and f^'' formulas, Sturm (V(a) − V(+∞) with the
gcd-terminated sequence; a = 0, w0 are not roots), lattice indexing (k = i + j into F2; A at lattice points by exact rational
interpolation; δ = dg/5 aligned with the nodes), R2 (A(y) = A(U) = A[77] for y ≥ U), R3 (A(U) − min A + β(1 − U)), M2 grid covers
[λ, 2U], PASS requires all three: all correct. Cosmetic: l.94 comment says "+1e-15" but no slack is added there (float rounding
of A ~10⁻¹⁶ is covered by the final +10⁻¹² and by the FULLARB pass). Page erratum: discoveries/small-gap-hardcore.md displays
"c_cert ≈ −0.3478618", a mis-rounding of −0.3478616532 (to 7 decimals −0.3478617); it understates the certified bound by 1.5·10⁻⁷.
The statement is still true (sup Φ ≈ −0.3478846) but not what the checker certifies; the margin uses the full value.
Memory: xcheck_aug.py peaks at 266 MiB RSS (ru_maxrss; its last line builds 43 Hermite arrays of 400,001 floats); memcap's ~3 s
sampling reported 23 MB and reports 0 MB for runs shorter than one sample. Under the cap, but the watchdog misses short bursts.


## 6. Falsification attempts (none succeeded)
(a) Actual zeros (hc_zeros768.log; 22,491 zeros, Riemann/Deep/ExploreHelio/results/zeta_zeros.npy). T = 19998.7626 (between the
 last two zeros, so every γ_n ≤ T has its right neighbour), L = 9.90343, N(T) = 22490 (RvM main term 22489.8), Lamzouri's E_λ
 (l.405, fixed L, n ≥ 2) at λ = 0.768: #E_λ = 5886 = 0.262N, N_0 = 16604. Every step of §2 checked on this Z_0:
 - min Z_0 gap 0.768032 > λ ✓; k = 3 distances ≥ 2.54949 > U = 2.304, max f^ over all k ≥ 3 pairs = −0.0 ≤ 0, their sum −1.18 ✓;
 - (H) at all 16,602 consecutive pairs: max Φ(g_i, g_{i+1}) = −0.347999 ≤ c_cert = −0.3478617 ✓ (attained at (0.7681, 1.2335):
   real gap pairs come within 1.4·10⁻⁴ of the certified ceiling);
 - telescoping: Σ RHS(H) = 3898.3883 = closed form (N_0 − 2)c + A(g_{N_0−1}) − A(g_1) + β(x_{N_0−1} − x_1 − N_0 + 2) to 9·10⁻¹² ✓;
 - chain: S_w = 9284.45 ≤ S_{k≤2} = 9285.63 ≤ bound 24419.06 ✓; range check TL/2π − N(T) = 9031.6 ≤ 4N/L = 9084 ✓ (§2(4)).
 - Weights at finite T: |S_w − S_1| = 1130 = 0.068 per point at L = 9.9 (1 − w ≈ 0.09 at u = 1); asymptotically O(N/L²) as in
   §2(2), whose bound 4π²M_fN_0/L² = 5.19·10³ (M_f = 0.776 at u = 1.352, hc_mf768.log) covers it. Zeros escape the contradiction through #E_λ (26%): the β-term is 19,345.
(b) Synthetic hard-core sequences (hc_synth768.log), the EXACT finite inequality of §2 with w = 1:
 S = Σ_{i,j}f^(x_i − x_j) ≤ nf^(0) + 2f^(g_{n−1}) + 2[(n − 2)c + A(g_{n−1}) − A(g_1) + β(x_{n−1} − x_1 − (n − 2))], c = c_cert.
 Families: Tonks gas (λ + Exp, mean 1; 40 × 3000 gaps), all gaps = λ, runs of r λ-gaps + one long gap (r = 1..8, mean 1),
 alternating λ, 2 − λ, gaps exactly at A's nodes (±10⁻⁹), the sup-pair (λ, 1.11709) repeated, uniform [λ, 2 − λ], clusters of
 λ-gaps at both ends. max(S − bound) < 0 in every family (closest: −0.158 for the runs family). ✓
(c) Ground state of (H) (hc_synth768.log): coordinate ascent with restarts over cyclic sequences, n = 1..30, of
 E_β(g) = avg[f^(g_i) + f^(g_i + g_{i+1}) − β(g_i − 1)] (≤ c for every cyclic g ≥ λ by (H)): max −0.3479921 at g ≡ λ for every n,
 1.30·10⁻⁴ below c_cert. Full pair energy of that lattice 0.00309 ≤ f^(0) + 2c + 2β(λ − 1) = 0.00398 ✓.
(d) Sign/limit traps checked: x, y → ∞ (Q3, R3 = −0.34811 bound, sampled −0.348686); boundary x = λ (the sup sits there:
 (λ, 1.11709)); y = λ edges; the constant tail (both checkers interpolate with A(x) = A(U) beyond U); f^'s last zero 6·10⁻⁶ below U
 is harmless because k ≥ 3 distances are > 3λ = U strictly.

## 7. Verdict on Theorem 1 as written
Holds exactly as written, given Lamzouri's Prop SmallSetZeros (audited separately). Remarks, none affecting validity:
- The hypothesis "Schwartz" can be "even, ≥ 0, C² near [−θ', θ'], f ∈ L¹, sup u²|f^(u)| < ∞" (covers C_c² too); Lemma T′ is not needed.
- "λ < U" is never used; "A = β = c = 0, U = λ recovers Thm Main2" is literally outside "λ < U" but any U ∈ (λ, 3λ] works.
- The [T, 2T]/varying-log γ worry does not arise: Lamzouri's E_λ, D(λ) and the pair sum all use zeros in (0, T] and fixed L = log T.
- Corollary D(λ) > 0 for all λ ≥ 0.768: from cert768 (margin 0.00146; second route 0.00070) and monotonicity of D in λ.

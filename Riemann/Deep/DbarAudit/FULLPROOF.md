# DbarAudit — FULLPROOF: line-by-line chain of "RH ⇒ lim_{λ→0} D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|²" (dbaraudit, 20 Nov 2026)

Sources: [N] = Riemann/Deep/DbarNB/DERIVATION.md (line numbers l.), [R] = Riemann/Deep/ExploreDbar/REPAIR.md, [X] = ExploreDbar/DERIVATION.md.
Notation as in [N] l.3–5: L = log(1/λ), ℓ = L − √L, s₀ = L^{−3/4}, T₂ = t_{n₀} ∈ [L⁶, L⁶+4], ‖E‖² = ∫|E|²w, w = 1/(2π(¼+t²)).
Status marks: ✓ = re-derived by me line by line; (✓) = re-derived with a wording/constant correction that changes nothing.

## 0. External inputs, stated exactly from the source
- (TN) Friedland–Yomdin arXiv:1107.0039, Thm 1.1 (TeX l.72–78, attributed to Nazarov): p(t) = Σ_{k=0}^{m} c_k e^{λ_k t}, c_k, λ_k ∈ ℂ,
  B ⊂ ℝ an interval, Ω ⊂ B measurable: sup_B|p| ≤ e^{μ₁(B)max|Re λ_k|}(cμ₁(B)/μ₁(Ω))^{m} sup_Ω|p|, c > 0 absolute.
  So m+1 terms ↦ exponent m; real frequencies (λ_k ∈ iℝ) ↦ growth factor 1. [N] (T) l.20–23 (M terms, exponent M−1) matches.
  Jets: FY state no 1D polynomial-coefficient version (their quasipolynomial degree Σ(d_j+1), TeX l.401–405, is the multiplicity
  count). Continuity argument (mine): u^j e^{iγu} = j!·i^{−j}·lim_{ε→0}(e^{i·u})[γ, γ+ε, …, γ+jε]; so every element of
  span{u^j e^{iγ_k u}: j < m_k} is a uniform-on-compacts limit of sums of Σm_k pure exponentials with real frequencies; TN passes
  to the limit on the compact B. Exponent = (number of terms with multiplicity) − 1. ✓
- (S) Titchmarsh, Theory of the Riemann zeta-function, Thm 14.13: on RH, S(t) = O(log t/log log t). With Riemann–von Mangoldt,
  N(t+h) − N(t) = (h/2π)log(t/2π) + S(t+h) − S(t) + O(1/t) ≤ A₁ log t/log log t for h ≤ 1/log log t. [N] (Z) l.24–25 ✓
  Verbatim source: Goldston–Gonek arXiv:math/0511092, TeX l.107–122: N(t+h) − N(t) = (h/2π)log(t/2π) + S(t+h) − S(t) + O((1+h²)/t),
  "It was proved by Littlewood that on RH, S(t) ≪ log t/log log t", hence N(t+h) − N(t) − (h/2π)log(t/2π) ≪ log t/log log t for
  0 < h ≤ √t; l.127–132 sharpen it to (½+o(1))log t/log log t; l.171–176: m(γ) ≤ (½+o(1))log γ/log log γ. (Count (t, t+h]; for the
  closed [t, t+h] apply it at t − 0.)
- (C) Titchmarsh Thm 9.7: there is A such that each interval [T, T+1] contains t with |ζ(σ+it)|⁻¹ = O(t^A), −1 ≤ σ ≤ 2. [N] l.111 ✓
  (Book statement, 2nd ed. §9.7, from the standard text; no open verbatim quote found: arXiv search refused 3 queries, Ng math/0310381 has none.)
- (G) Gevrey-2 decay [N] l.26: |ζ|^N|φ̂(ζ)| ≤ ∫|φ^{(N)}| ≤ ℓBC^N(N!)²ℓ^{−N}; optimise N ≍ (|ζ|ℓ/C)^{1/2}: |φ̂| ≤ C′Bℓe^{−c(|ζ|ℓ)^{1/2}}, ζ ∈ ℝ. ✓
  The bump ψ(x) = S(4x)S(4(1−x)), S(y) = f(y)/(f(y)+f(1−y)), f = e^{−1/y}, is Gevrey-2 (used in my numerics).
- (Bu) Burnol arXiv:math/0103058, TeX l.169–175: K = L²(]0,∞[, dt), χ = 1_{]0,1]}, B_λ = finite linear combinations of
  t ↦ ρ(θ/t), λ ≤ θ ≤ 1 (no constraint Σc_kθ_k = 0); l.189: D(λ) := inf_{f∈B_λ}‖χ − f‖; l.203–206: liminf D(λ)√log(1/λ) ≥
  √(Σ_ρ m_ρ²/|ρ|²); l.197–201: if RH fails the left side is +∞ ("true but trivial"), so the bound is unconditional. Same D(λ)
  as [X] §0 / [N] (brief item (f)). ✓
- (PWS) Paley–Wiener–Schwartz: entire G with |G(s)| ≤ C(1+|s|)^N e^{a·max(0,−σ)} is ∫e^{−us}dT(u), T a distribution, supp T ⊂ [0,a].

## 1. Dual form [N] §1 l.7–18 ✓
- ∫e^{−ixt}w(t)dt = e^{−|x|/2}: ‖E‖² = ∫∫e^{−|u−v|/2}dμdμ̄; OU kernel e^{−θ|u−v|}, θ = ½, RKHS norm on [0,ℓ]:
  (2θ)⁻¹[∫(|p′|²+θ²|p|²) + θ(|p(0)|²+|p(ℓ)|²)] = N_ℓ(p). Constraint functionals E^{(k)}(γ) = ∫conj((iu)^k e^{iγu})dμ, targets (1,0,…);
  min-norm duality U = sup_a|Σa_{γ0}|²/‖p_a‖²_H and Σa_{γ0} = p_a(0). Check p = e^{iγu}: N = ℓ(γ²+¼)+1 = K(γ,γ). ✓
- Rerun from a copy (rerun/, logs/rerun_dual.log): Gram identity to 9.7e-17, L·U = 0.0403698458 / 0.0382192552 / 0.0377548398
  (L = 5, 20, 80): identical to DbarNB/logs/dual.log.
- Consequence used: ∫|p′|²χ ≤ ∫_0^ℓ|p′|² ≤ N_ℓ(p) (χ ≤ 1). ✓

## 2. Grouping [N] §3 l.28–35 ✓
- Groups = maximal runs of consecutive ordinates (sorted, with multiplicity) with gaps < s₀ (in 1-D the "gap < s₀" components are
  runs), so (i) convex hulls I_g are disjoint with gaps ≥ s₀ ✓; (ii) w_g ≤ (M_g−1)s₀/2 ✓.
- (iii) M_g ≤ M*(|γ_g|), M*(t) = A₁ log t/log log t: if M_g > M*(t), t = least |γ| of g, the first ⌊M*(t)⌋+1 ordinates lie in
  [t, t + M*(t)s₀], and M*(t)s₀ ≤ 7A₁ log L·L^{−3/4} < 1/log log t for t ≤ L⁶+4, L ≥ L₀: contradicts (S). ✓
  So e^{O(M_g)} ≤ exp(O(log|γ_g|/log log|γ_g|)) = |γ_g|^{o(1)} ≤ L^{o(1)} (RH use (b)). ✓
- Below t₀: finitely many zeros, M_g ≤ N(t₀); no group straddles ±t₀ once s₀ < dist(t₀, ordinates). ✓

## 3. Lemma A [N] l.37–45 ✓ (frequencies δ_j ∈ [−w,w], W = max(w, KM²/ℓ))
- m(δ) = e^{iδζ}χ₀(δ/W) = e^{iδζ} on |δ| ≤ 3W/2 ⊃ [−w,w], so ∫φ(s)(x+s)^k e^{iδ(x+s)}ds = (−i∂_δ)^k[e^{iδx}m(δ)] = (x+ζ)^k e^{iδ(x+ζ)}
  at every node (jets need m = e^{iδζ} on a neighbourhood: yes). φ(s) = WΦ_a(Ws), Φ_a = (e^{iηa}χ₀(η))^∨, |a| ≤ 1: Gevrey-2
  uniformly in a, |Φ_a(τ)| ≤ Ce^{−c|τ|^{1/2}}. ✓
- Tail: TN with B = hull([0,ℓ] ∪ {x+s}), Ω = [0,ℓ]: |r(x+s)| ≤ (A(1+|s|/ℓ))^{M−1}sup|r|; ∫|Φ_a(τ)|A^{M}(1+|τ|/(KM²))^{M}dτ
  = A^M CKM²∫e^{−M(c√K σ^{1/2} − log(1+σ))}dσ ≤ 4CKA^M once c√K ≥ 2 (log(1+σ) ≤ σ^{1/2} ⇔ (1−√σ)² ≥ 0). ✓ Cauchy on |ζ| = 1/W. ✓

## 4. Lemma B [N] l.46–53 ✓ (|γ_g| ≥ 14, I = [ℓ/4, 3ℓ/4])
- q_g = e^{iγ_g u}(iγ_g r + r′) ⇒ |γ_g r| ≤ |q_g| + |r′|; S = sup_I|r|² ≤ (4/ℓ)(2A)^{2M}∫_I|r|² (TN-Chebyshev set |E| ≥ |I|/2; the
  source's 8/ℓ is weaker); ∫_I|r′|² ≤ (ℓ/2)C_A^{2M}W²(2A)^{2M}S; absorb if 8(2A)^{4M}C_A^{2M}(W/γ_g)² ≤ ½; then
  |p_g(0)|² ≤ sup_{[0,ℓ]}|r|² ≤ (2A)^{2M}S ≤ 32(2A)^{4M}(ℓγ_g²)⁻¹∫_I|q_g|² ≤ B(M)(ℓγ_g²)⁻¹‖q_g‖²_χ. ✓
- Hypothesis: W_g ≤ 2M_gL^{−3/4}, γ_g ≥ 14, e^{O(M)} = L^{o(1)} ⇒ holds for all groups at once for L ≥ L₀. ✓ Uniform in the width
  (brief (b)): only W_g/γ_g enters, and W_g ≤ M_g s₀ ≪ 1.
- Numerics (audit_num.py B; logs/B.log; mpmath Newton basis + double QR, quadrature h-converged to 1e-14): the γℓ → ∞ limit
  N_B = sup|p_g(0)|²ℓγ_g²/‖p_g′‖²_χ over E_g is largest for the coalescent jet: 1.42, 9.83, 49.4, 219, 2.94e3, 2.86e4 for
  M = 1, 2, 3, 4, 6, 8 (≈ 3.4^M: e^{O(M)}, as B(M) allows); s₀-chains and random gaps (D = s₀ℓ = 2.1…31.6) are smaller; at finite
  γℓ = 280 the ratio is ≤ the limit in every case (0.58–1.00 of it): the γ² gain holds uniformly. No counterexample.

## 5. Lemma Λ [N] l.55–65 ✓
- Newton–Hermite: with D_k := Π_{l<k}(−i∂_v − δ_l), D_kν_j ≡ 0 for j < k (the multiplier Π_{l<k}(δ′−δ_l) kills the divided
  difference over δ_1..δ_j), and D_kν_j(0) = (Π_{l<k}(δ′−δ_l))[δ_1..δ_j] = δ_{jk} for j ≥ k (degree k−1). So c_k = D_k r_q(x). ✓
- |e_{k−1−n}(δ)| ≤ C(k−1,n)w^{k−1−n}, Lemma A ⇒ |c_k| ≤ eC_A^M(k−1)!W^{k−1}sup|r_q| ✓; Hermite–Genocchi ✓; Leibniz for φ_k:
  |φ_k^{(N)}| ≤ (N!)²ℓ^{−N}C^{N+1}(ℓ/2)^{k−1}Σ_j(2(k−1)/C)^j/j! ✓; (G) ⇒ |φ̂_k| ≤ C^kℓ^ke^{−c(|ζ|ℓ)^{1/2}} ✓ (real ζ: RH use (a)).
- a_g² = e^{O(M log(Wℓ+2))}ℓ‖q_g‖²_χ with M log(Wℓ) ≤ (6A₁/4 + o(1))log²L/log log L = o(log²L). ✓

## 6. Lemma D [N] l.67–73 ✓
- Windows separated by ≥ s₀ ⇒ off Ω_g, d_g ≥ s₀/2 (ξ ∉ Ω_g lies nearer another window, d_g + d_{g′} ≥ gap) ✓; per side at most one
  window can start in a shell of width s₀/2 (consecutive starts differ by > s₀) ⇒ ≤ 2 per shell ✓.
- (a) ∫_{d≥s₀/2}e^{−2c(dℓ)^{1/2}}dξ = (2/ℓ)∫_{s₀ℓ/2}^∞e^{−2c√x}dx ≍ ℓ^{−1}L^{1/8}e^{−2c(s₀ℓ/2)^{1/2}}, times a_g² ∝ ℓ: τ = e^{o(log²L)−cL^{1/8}} ✓.
- (b) Σ_g∫_{Ω_g}Σ_{g′≠g}a_{g′}²e^{−c(d_{g′}ℓ)^{1/2}} = Σ_{g′}a_{g′}²∫_{Ω_{g′}^c}e^{−c(d_{g′}ℓ)^{1/2}} ≤ τΣ‖q_{g′}‖²_χ ✓: every bound is relative to the
  group's own norm and the shell count is 2 per shell, so nothing depends on the number of groups (~L⁶log L) or on the
  conditioning inside a group (brief (c)). (c) 2|a||b| ≤ η|a|² + η⁻¹|b|² ✓. Z(T₂) is finite: all sums finite ✓.
- Numerics (audit_num.py D; logs/D.log): two groups at gap exactly s₀ (D = s₀ℓ), cos of the principal angle between E_g and E_{g′}
  in L²(χ), and K3 = sup Σ_g‖q_g‖²/‖q‖² for three groups. Coalescent jets, M = 2: cos = 0.9988, 0.847, 0.134, 1.7e-3 at D = 4.2,
  10, 31.6, 100 (L = D⁴ ≈ 3e2, 1e4, 1e6, 1e8). So the lemma is invisible below L ≈ 10⁶ and holds only asymptotically (L₀ huge,
  as the source's caveat says); the decay in D is consistent with e^{O(M log(WD))}·e^{−c√D}. Results for M = 4, 8 in §12.

## 7. Theorem R″ [N] l.75–81 ✓
- p = Σ_g p_g (unique: the u^k e^{iγu} are independent); |p(0)| ≤ Σ_g|p_g(0)| ≤ Σ_g(B_g/(ℓγ_g²))^{1/2}‖q_g‖_χ ≤ (S_L/ℓ)^{1/2}(Σ‖q_g‖²_χ)^{1/2}
  ≤ ((1+ε_L)S_L/ℓ)^{1/2}‖p′‖_χ ≤ ((1+ε_L)S_L/ℓ)^{1/2}N_ℓ(p)^{1/2} ✓ (q_g = p_g′ ∈ E_g ✓).
- S_L = Σ_g B(M_g)/γ_g² ≤ C(t₀) + 2Σ_ρ e^{βM*(|γ_ρ|)}/γ_ρ² (each group charged to its largest |γ|; |γ_g| ≥ |γ_ρ| − w_g ≥ |γ_ρ|/√2) —
  a convergent series independent of L (Σ_ρ|γ_ρ|^{−2+o(1)} < ∞ by N(T) ≪ T log T). ✓ Hence D(λ)² ≪ 1/log(1/λ) with §11.

## 8. Lemma E, Theorem R′, Corollary [N] l.83–103 ✓
- Boundary terms: |p_lo(0)|² ≤ U_{Z(H)}(ℓ)N(p_lo) ≤ C(H)ℓ⁻¹N (dual form + Claim F); |p_hi(0)|² ≤ Cℓ⁻¹N(p_hi) (R″ on high groups);
  u = ℓ by p ↦ conj(p(ℓ−·)): conj((ℓ−u)^k e^{iγ(ℓ−u)}) = e^{−iγℓ}(ℓ−u)^k e^{iγu} ∈ E_Z, N preserved ✓.
- Interior: Gram of ℓ^{−k−1/2}e′_{ρk} → blockdiag(γ_ρ²·[1/(j+k+1)]) (Hilbert matrices PD; cross blocks O(1/(ℓ·gap))) ✓; one
  integration by parts against e^{−iΩu}, |Ω| ≥ g_H, with Lemma A for r_g′ and the TN by-product: C(H)ℓ^{k−1/2}L^{1/4+o(1)}|Ω|⁻¹‖q_g‖_χ ✓;
  Σ_g 1/Ω_g² ≤ C(H) ✓; Σ|a_{ρk}|ℓ^{k−1/2} ≤ C(H)ℓ⁻¹‖p_lo′‖ ✓ ⇒ θ_L = C(H)L^{−3/4+o(1)} ✓. The ¼∫p_lo p̄_hi term is smaller ✓.
- R′: |p(0)| ≤ U_lo^{1/2}N_lo^{1/2} + V^{1/2}N_hi^{1/2} ≤ (U_lo+V)^{1/2}(N_lo+N_hi)^{1/2}; N(p) ≥ N_lo + N_hi − 2θ_L(N_loN_hi)^{1/2} ≥ (1−θ_L)(N_lo+N_hi) ✓;
  S_hi ≤ ε(H) uniformly in L (brief (d)) ✓. Corollary: U_{Z(L⁶)} ≥ U_{Z(H)} (more constraints) + Claim F at band ℓ (L/ℓ → 1) ✓.
- Note: [N] proves R′ for Z(L⁶) but §6 needs Z(T₂), T₂ ≤ L⁶ + 4: the proof only uses log T₂ ≤ 6 log L + o(1): verbatim. (✓)

## 9. §6 of [N]: Lemma C, tilings, (3v) l.110–149 ✓
- Tiling (brief (e)): t_n ∈ [3n,3n+1], t′_n ∈ [3n+3/2, 3n+5/2] from (C) (ζ ≠ 0 on the horizontal sides, −1 ≤ σ ≤ 2); every zero
  with |γ| > T₂ has 0 < β < 1 and lies in exactly one Q_n and one Q′_n; narrow and wide rectangles of equal heights hold the same
  zeros, so R_n = R^w_n identically (R_n is a residue sum, it depends only on the zeros enclosed) ✓. Partial sums of the two
  height families differ by the residues between t_N and t′_N = a contour integral → 0 for fixed s ✓: one limit R_high, no
  absolute convergence over single zeros (no bound on 1/ζ′(ρ)) ✓.
- Kernel: Res_{w=ρ}F₀(w)K(s−w), K(x) = B(x)/x = 1/x + entire, has principal part Σ_j c_j(s−ρ)^{−j} = that of F₀ at ρ ✓ (also for
  multiple poles: ((−1)^{j−1}/(j−1)!)K^{(j−1)}(x) = x^{−j} + entire).
- (i) Residue at w = s: −F₀(s)B(0) ✓. (ii) ζ(s)/(ζ(w)(s−w)) = ζ[s,w]/ζ(w) + 1/(s−w) and (2πi)⁻¹∮H(w)B(s−w)/(s−w)dw = −H(s)1_{s∈Q_n}
  ⇒ ζR_n = (2πi)⁻¹∮HBζ[s,w]/ζ(w)dw for s ≠ 1, continuous across ∂Q_n ✓. At s = 1 (divided difference has ζ's pole): not on the
  line; R_high(1) = Σ_n(I_n(1) + Ī_n(1)) directly since 1 ∉ Q_n ✓; G(1) is needed only for Nyman's constrained variant.
- (iii) on ∂Q_n: |Φ| ≤ Ce^{−c(t_n/L²)^{1/2}} ≤ Ce^{−cL²}, |Ψ| ≤ C(1+|w|)√Le^{L′}, |1/ζ| ≤ t^A (sides at good heights; σ = 3/2, −½
  bounded), |ζ[s,w]| ≤ C(2+|t|+t_{n+1})², |B| ≤ e^{L−1}min(1, 16/((L−1)dist)²): for |t| ≫ t_n the factor |t|²·|t|^{−2} stays bounded,
  so the w-integral converges and ‖ζR_high‖² ≤ e^{−c′L²} ✓ (t^K e^{−2c√t/L} decreasing for t > K²L²/c² < L⁶).
- (iv) Locally uniform convergence off the high zeros ✓; G = F₀ − R_high − R_tr entire ✓. Growth: strip
  σ ∈ [−¾, 7/4], |t| > T₂+4: the height family with sides ≥ ¼ from t, wide width, gives G = −Σ_m I^w_m − R_tr (the one rectangle
  containing s contributes I^w + F₀) ✓, |s−w| ≥ ¼, |B(s−w)| ≤ e^{(L−1)(2+max(0,−σ))} ✓, |1/ζ(−1+it)| ≪ |t|^{−3/2} ✓; compact part ✓;
  σ > 7/4 and σ < −¾ with the narrow family (s at distance ≥ ¼ from all Q_m) ✓; |1/χ(s)| ≤ C off |s+2k| < ½ for σ < −¾ ✓;
  maximum principle inside the discs ✓. Exponent L′ + 2/T₁ ≤ L − ½ ✓ ⇒ PWS: supp ⊂ [0, L−½] ✓.
- (v) 1 − ζG = (1 − Φ) + ΦE_low + ζR_tr + ζR_high on the line (ζF₀ = ΨΦ) ✓; ‖ΦE_low‖ ≤ ‖E_low‖ = U^{1/2} (|Φ| ≤ 1) ✓;
  ‖1−Φ‖² ≪ 1/T₁ = L⁻² ✓; (1+η)‖b‖² + (1+η⁻¹)‖a‖², η = L^{−1/2}: D(λ)² ≤ (1+L^{−1/2})U_{Z(T₂)}(L−√L) + O(L^{−3/2}) ✓ (RH not used).

## 10. explore59's accepted pieces (brief (g)) ✓
- [R] (3i) window: Φ = e^{−iz/T₁}F(z/T₁) has spectrum in [0, 2/T₁] (h ↦ [−½,½], h² ↦ [−1,1], shift) ✓; |Φ| = F ≤ 1 on the line ✓;
  Gevrey decay uniform in |σ − ½| ≤ 1 (φ(ξ)e^{−yξ} uniformly Gevrey) ✓; |1−Φ(t)| ≤ |t|/T₁ + Ct²/T₁² ⇒ ‖1−Φ‖² ≪ 1/T₁ ✓.
  Reproducing bound |Ψ(z)| ≤ 1 + ‖E_low‖K_{L′}(z,z)^{1/2}, ‖E_low‖ ≤ ‖1‖ = 1 (E ≡ 1 is admissible) ✓. s = 1: none needed for the
  unconstrained D (Bochner argument), and Burnol's D is the unconstrained one (§0 (Bu)) ✓; constrained variant via G − G(1) ✓.
- [R] (3ii): residue of F₀ at −2k is Ψ(−2k)Φ(−2k)/ζ′(−2k) = r_k (trivial zeros simple) ✓; |E_low(−2k)| ≤ K_{L′}(iy,iy)^{1/2}, y = 2k+½,
  K(iy,iy) = e^{2L′y}(y+½)²/(2y) + O(y) ✓; on the line each term ≤ C√k e^{−(√L−1−2/T₁)(2k+½)}/(|ζ′(−2k)||2k+½+it|) ✓.
- [X] §7 closure: T ∗ φ_δ ∈ C^∞_c([0, L−½+δ]); u ↦ ζ(s)e^{−us}/s is L²-continuous, so ∫h_δ(u)ζe^{−us}/s du is a Bochner integral in the
  closed span ✓; |φ̂_δ(½+it)| ≤ ∫φ_δe^{−u/2} ≤ 1, φ̂_δ → 1, majorant 2|ζG/s| ∈ L² ✓.
- [R] Claim F: jet Gram (¼+γ²)ΔMΔ*, M_jk = ∫_0^L u^{j+k}, [M⁻¹]₀₀ = Σ_{k<m}(2k+1)/L = m²/L (P_k(−1)² = 1) ✓; endpoint/cross terms
  relative O(1/(|γ|L)), O(1/(Ld)) ✓. Value: Σ_ρ 1/(ρ(1−ρ)) = 2 + γ − log 4π and ρ(1−ρ) = |ρ|² on the line ✓.

## 11. Assembly [N] §7 l.151–154 ✓
limsup L·D² ≤ limsup (1+L^{−1/2})L·U_{Z(T₂)}(ℓ) + O(L^{−1/2}) = Σ m²/|ρ|² (Corollary); liminf L·D² ≥ Σ m²/|ρ|² (Bu). ∎
RH is used exactly twice, both inside §§2–8: (a) real frequencies (TN without growth factor; (G) at real arguments), (b) (S).

## 12. Falsification numerics (all single-threaded, memcap peaks ≤ 227 MB; logs/ in this folder)
| test | configuration (RH-allowed worst cases) | result | verdict |
|---|---|---|---|
| Lemma B ratio (B.log) | one group, M = 1…8: coalescent jet, tight (0.01·s₀), s₀-chain (0.99·s₀), random gaps < s₀; D = s₀ℓ = 2.1…31.6 | max = jet: 1.42, 9.83, 49.4, 219, 2.94e3, 2.86e4 (≈ 3.4^M); finite γℓ = 280 ≤ limit; h-conv 1e-14; two quadratures agree to all printed digits (B_v1_vs_v2.txt) | e^{O(M)}, γ² gain uniform: no counterexample |
| Lemma D (D.log) | 2 and 3 groups at gap exactly s₀, jets and s₀-chains, M = 1, 2, 4, 8, D = 4.2, 10, 31.6, 100 | cos(principal angle) 6.7e-1 → 4.4e-4 (M = 1), jets M = 8: ≥ 0.9995 until D = 31.6, 0.109 at D = 100; K3 = sup Σ‖q_g‖²/‖q‖²: 1.0004–1.18 at D = 100, ≥ 1e10 (numerically singular) for jets M ≥ 4 at D ≤ 31.6 | decay ≈ e^{−√D} with M-dependent prefactor, as Lemma Λ allows; asymptotic only (D = L^{1/4}: L₀ ≫ 10⁸) |
| cost, rerun (rerun_clusters_m5_h0.3.log) | DbarNB clusters mode m = 5, h = 0.3 from a copy | ratios 1.06363, 1.01436, 1.00045, 1.00118, 1.00009: identical to DbarNB | reproduced |
| cost, chains (chain.log) | every 8th zero → m = 6, 8 zeros at spacing 0.99·s₀ (one proof-group); dense: every 2nd zero → 4-cluster at 0.3/L | L·U bounded (0.078–0.146), U_all/Σ_g U_g = 1.19 (L = 10) → 1.00007 (L = 160) | no counterexample |
| precision (mpcheck.log) | reduced chain config (60 zeros, m = 6 at 0.99·s₀), L = 40, 160, mpmath dps 30 vs double | ratios 0.999655810642, 0.999860709914 in both (12 digits) | the ratio < 1 is real (groups help slightly), not roundoff |
First D run killed by memcap (leggauss(4800) companion matrix, 520 MB; logs/D_killed520MB.*); rerun with composite 48-node panels.

## 13. Points the referees did not record (none breaks the chain)
1. [N] (T) cites Nazarov "Thm 1.3/1.4" (numbering not checked); the citable statement is FY Thm 1.1 (§0); the jet extension needs
   the continuity remark of §0 (absent from FY and from [N] beyond "by continuity").
2. [N] proves R″/R′ for Z(L⁶); §6/§7 use Z(T₂), T₂ ≤ L⁶ + 4: same proof (§8 note).
3. [R] (3i) "Growth" used RH bounds for 1/ζ (Titchmarsh 14.14); [N] §6(iv) replaces them by the tiling, so (3v) is unconditional.
4. Lemma B: S ≤ (4/ℓ)(2A)^{2M}∫_I|r|² suffices; the stated 8/ℓ is weaker (harmless).
5. L₀ is not just large: Lemma D's mechanism is numerically invisible below D = s₀ℓ ≈ 30–100, i.e. L ≈ 10⁶–10⁸ (§12).

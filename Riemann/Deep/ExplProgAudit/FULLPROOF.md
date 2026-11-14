# ExplProg audit — the chain, step by step (source = Riemann/Deep/ExplProg/DERIVATION.md as of 14 Nov 2026, 285 lines: "D:<n>"; cert.py 679 lines: "C:<n>")

Status marks: OK = rechecked by the auditor; GAP = missing/false step; NOTE = harmless imprecision.

## §0 Inputs (D:8-28)
- (H) B = 0.8733819666132076 > B₂ + 2.2e-254: recheck numerically below (A1).
- (FE) ρ ↦ 1 − ρ̄ is an involution of the zero multiset of a primitive L(s,χ) (FE links L(s,χ), L(1−s,χ̄); conjugation links χ̄ back to χ); fixes γ. With (H): b ≤ β ≤ B. OK.
- (BMOR) statement to be rechecked in TeX (A2).

## §1 Lemma EF (D:32-46) and orthogonality (D:48-58)
- Mellin pair W(s)=∫w t^{s−1}dt, Σ Λχ(n) w(n/x) = (1/2πi)∫_{(2)} (−L'/L)(s) x^s W(s) ds. Residues: +xW(1) (pole of ζ), −m x^ρ W(ρ), −x^{−m}W(−m) at trivial zeros; T(1)={2,4,..}, even χ≠1: {0,2,..}, odd: {1,3,..}. OK.
- W(s) = (s(s+1))^{−1}∫t^{s+1}dw'(t): boundary terms vanish because w, w' vanish outside [t₀,t₁] (w continuous, piecewise linear kernels). OK.
- Left shift: x^σ t₀^{σ+1} = t₀ (x t₀)^σ → 0 as σ → −∞ iff x t₀ > 1. OK (qualitative DAV inputs suffice).
- (1.1): x^{−m}W(−m) = ∫ w(t) t^{−1}(x t)^{−m}dt ∈ [0,(x t₀)^{−m}W(0)]; Σ_T (xt₀)^{−m} ≤ 1/(1−(xt₀)^{−2}) for each of the three T. OK.
- Orthogonality φ(q)S = Σ_χ χ̄(a)Σ Λχ w. Imprimitive difference = Σ_{p|q,p∤q*} Σ_k χ*(p^k) log p w(p^k/x); #{k: p^k∈[xt₀,xt₁]} ≤ log(t₁/t₀)/log p + 1 ⇒ (1.2). OK.
- (1.3) from EF + (1.1) + (1.2); no primitive character of conductor 2, so Ω ranges over 3 ≤ q* | q. OK.

## §2 Lemma P pairing (D:62-69)
- Pair {β+iγ, 1−β+iγ}, 1/2<β≤B, equal multiplicities, same |γ|; singleton β=1/2 bounded by min ≤ mean. OK.

## §2 Lemma Z (D:71-95), Zζ (D:96-99), Zζ' (D:101-104)
- (i) BMOR upper side: N ≤ L a(t) + b(t) + 2log(1+ℓ) using −χ(−1)/4 ≤ 1/4 and ℓ = L + c; concavity 2log(1+L+c) ≤ 2log(1+L) + 2c⁺/(1+L) (c<0: monotone). R(t;q*) increasing in q*. OK.
- (i) ℓ* ≤ 1.567 branch: needs N(T,χ)=0 from BMOR (check TeX, A2) and R(t;9) ≥ 0 on t ≥ 5/7: d/dt (t/π)log(9t/2πe) = log(9t/2π)/π > 0 iff t > 2π/9 = 0.6981 < 5/7. R(5/7;9) recomputed (A3).
- (iii) Σ_ρ G(|γ|) = ∫_0^∞ N(t)(−G'(t))dt (Fubini; real zeros γ=0 included). Upper bound uses N ≤ N⁺ on [5/7,∞), N ≤ N⁺(5/7) below, −dG ≥ 0; then integration by parts (N⁺G → 0 since G ≤ k₂/t²). dN⁺/dt = L/π + ρ₀ + e'/(1+L), e' = 2/(t+2) on t > 2π−2. Collecting terms = z_a, z_b, z_e. OK.
- (iv) antiderivatives: re-verified symbolically (A4).
- Lemma Zζ: pair sum β/|ρ|² + (1−β)/|1−ρ̄|² ≥ 1/(B²+γ²) (both moduli² ≤ B²+γ² since 1−β ≤ 1/2 < B); singleton 1/γ² = 2(1+1/(4γ²))β/|ρ|²; |γ| ≥ 14.13. Σ_ρ Re(1/ρ) = 1+γ_E/2−½log4π over ALL zeros (both signs of γ). OK; number rechecked (A3).
- Lemma Zζ': 2·(7/8+1.588) = 4.926, 2·0.137, 2·0.443. OK (Rosser Thm 19, T ≥ 2: statement-level).

## §3 Kernels (D:106-138; K1 D:117, K2 D:132)
- 3.1: W(s)=∫g(u)e^{(s−B)u}du; ∫g e^{mu} = (1/m²)∫e^{mu}dg', masses 1/a, −(1/a+1/c), 1/c. n₀=7/2, n₁=2, n₂=2(1/a+1/c). ‖w‖_∞: on [−c,0] max at −c (e^{Bc}); on [−ℓ,−c] critical point u=−ℓ+1/B inside since a>1/B: e^{Bℓ−1}/(aB). OK.
- g_L=g e^{mu}, m=1/2−B<0: ‖g_L'‖₁ ≤ ∫|g'|e^{mu}+|m|m₀; Var(g_L') ≤ ∫e^{mu}|dg'| + 2|m|∫|g'|e^{mu} + m²m₀. OK.
- K1: h=g e^{−ε(u+X)}, u+X ≥ Y on supp g; |h|,|h'|,|dh'| bounded termwise with sup_ε ε^m e^{−εY} = (m/(eY))^m. OK (sup of sum ≤ sum of sups).
- 3.2: w₋ ≤ 1_{(α,1]} ≤ w₊ needs α(1+δ) ≤ 1−δ: α ≤ 1/8, δ ≤ 1/2 ⇒ OK. W_±(1), W_±(0) formulas OK.
- K2: xt ≥ √x(1−δ) ≥ 4 ⇒ (xt)^β ≤ (xt)^B; ramp masses give ∫t^{B+1}|dw'| ≤ ((1+δ)^{B+1}+1+2α^B)/δ (w₊) and ((1+δ)^{B+1}+1)(1+α^B)/δ (w₋) ≤ A₂; ∫|w'|t^B ≤ u₁^B+v₁^B ≤ A₁; j=0 via w_± ≤ w_δ, W_δ(B). OK.

## §4 Theorem 2 (D:140-173; tail D:155-163, small q D:164, Thm 2' D:167)
- P(x) ≥ S − ‖w‖_∞(ψ−θ)(x); ψ−θ = Σ_{k≥2}θ(x^{1/k}) ≤ 1.01624(√x + x^{1/3}X/log2) (RS62 Thm 9). (4.1) = (1.3) with φ(q) ≤ q, D_q ≤ ‖w‖_∞(7ω+L), ω ≤ L/log 2, (xt₀)^{−2} = e^{−2Y}. OK.
- Tail: x^b = c₂qL ⇒ x/(q x^B L) = c₂; F(q) term by term OK. Freezing k's at Y_m is an upper bound because Z_G(L) = N⁺(5/7)n₀ + ∫N⁺'G with N⁺' = L/π + ρ₀ + … > 0 for L ≥ log 9 (min ρ₀ = ρ₀(5/7) ≈ −0.609 > −log9/π = −0.699): Z_G is increasing in G even though z_b alone is not. (The page does not state this; it holds.) NOTE.
- h(L) monotonicity criterion and the max(z_b,0), max(z_e,0) fallbacks: OK. x^{1/3−B}X decreasing for X > 1/(B−1/3) = 1.85. OK.
- Rounding 1/b = 7.897769 up to 7.8978 is safe: every base (q log q, c₂ q log q, φ(q)log q, log q) is > 1 in its range. OK.
- Thm 2': f(q) = q/(e^γ loglog q + 2.50637/loglog q) ≤ φ(q) (RS62 Thm 15, n ≥ 3); f' sign recomputed = e^γ(u − 1/log q) + 2.50637(1/u + 1/(u² log q)) > 0 for q ≥ 16. OK.

## §5 Theorem 3 (D:175-188)
- D(x) = Σ Λ(1−Re χ)w ≥ 0; if χ(n) ∈ {0,1} for n ≤ x then D ≤ Σ_{p|q} ≤ D_q. Lower bound xW(1) − Ω₁ − Ω_{χ*} − 2W(0)/(1−e^{−2Y}) − D_q. Contradiction condition has 2D_q. OK.
- After dividing by x^B L the Ω₁ term is 0.0231856(k₂ + x^{1/2−B}k₂^L)/L (no 1/q): "verbatim" needs checking in cert.py (C-check T3-Ω₁).
- q − 1 ≤ (log q)^{1/b} on 3 ≤ q ≤ e^{25.6}; minimal n prime by complete multiplicativity. OK.

## §6 Theorem 1 (D:190-224; 6.1 D:198, 6.2 D:200, 6.3 D:203, cells D:209, 6.4 D:214)
- 6.1 x<q: one n ≤ x in the class; φ(q) ≥ √(q/2) ⇒ x/φ < √(2x); log(qx) > 2log x. R(2) = max(0.19687, 0.56805) recomputed by hand. OK.
- 6.2: ψ(·;q,a) constant on [n,n+1), E linear ⇒ |E| ≤ max at piece ends, denominator ≥ value at s₀. OK. q=1,x∈[2,2.0005]: (2.0005−log2)/(2^B log²2) = 1.48533. Forced value at x=2: 1.30685/0.880182 = 1.48476. OK.
- 6.3: S₋ ≤ ψ(x)−ψ(y) ≤ S₊ in the class; xW₊(1) = x−y+(x+y)δ/2; (1.3) ⇒ |E| ≤ T₁+…+T₆ with (φ−1)Ω/φ ≤ Ω. T₅ via W_±(0) ≤ X/2+log((1+δ)/(1−δ)); T₆ via ω ≤ L/log2, ‖w_±‖_∞=1. q=1,2: Ω-term absent (no conductor ≥ 3 divides). OK.
- Cells: φ_lo = max of two valid lower bounds (√(n/2); RS62 Thm 15 with each term of the denominator bounded on [log L₀, log L₁], log L₀ > 0 since L₀ ≥ log 3). OK. x ∈ [64, e^{4.16}) lies in 6.2 (e^{4.16} = 64.07 < 66). Coverage: code check (C).
- 6.4 tail: t_a = A₁/W_B ∈ [1/N₀, K₁B] (W_B ≥ 1/B); t_b = ((1+δ)+(1+δ)^{−B})/δ ∈ [2/δ,(2+δ)/δ]; ∫G = 2k₁ + k₁log(t_b/t_a) − (5/7)n₀; ρ₀ < 0 on [5/7,1.25]; z_b, z_e, Σ_ζ bounds re-derived term by term (k₂/t_b = k₁, ∫_{t_b}^∞ 2/(t²(t+2)) ≤ 1/t_b²). M bound for left kernels OK. Limit ½(b/(4π)+b²/π) = 0.0075896 = term-by-term limit (L/(L+X)² ≤ 1/(4X)). OK.

## §7 Theorem 4 (D:226-253; EF1 D:232)
- EF1: V(s) = s^{−1}h^{−1}∫₁^{1+h}t^s dt; pole at 0 residue 1; V(−1) = −log(1+h)/h; |V(−1−m)| = (1−(1+h)^{−m})/(hm(m+1)) ≤ 1/(m+1) (Bernoulli). OK.
- G = min(1/b,1/t,((1+h)^B+1)/(ht²)) for β ∈ [1/2,B] (|s| ≥ max(1−β,|γ|), x^{β−1} ≤ x^{−b}); G_L with 2, (1+h)^{1/2}. OK.
- Prime part: Σ_{n≤y}Λ(n)/n ≤ Σ_{p≤y}log p/p + (−γ_E − E); Dusart statement: TeX check (A5).
- Tail x = (κL)^{1/b}: X − (1/b)loglog q = (1/b)log κ. OK.

## Source checks
- A2 BMOR (arXiv 2005.02989, only version v1 = Math. Comp. 90 (2021) 1455-1482; TeX CountingZeros-May5-2020.tex): Thm 1.1 at l.170-175: "χ a character with conductor q > 1, T ≥ 5/7, ℓ = log(q(T+2)/2π). If ℓ ≤ 1.567 then N(T,χ)=0. If ℓ > 1.567 then |N − ((T/π)log(qT/2πe) − χ(−1)/4)| ≤ 0.22737ℓ + 2log(1+ℓ) − 0.5." l.157-165: N counts ρ with 0<β<1, |γ| ≤ T, with multiplicity; Z(χ) = Z(χ*) for imprimitive χ. Matches D:15-18 exactly. OK.
- Low zeros at height 0 (brief question): legitimate. BMOR is applied only at T = t ≥ 5/7; for t < 5/7 only monotonicity N(t) ≤ N(5/7) is used and every zero of height < 5/7 gets the weight G(0) = n₀ = max G. No extrapolation of BMOR below 5/7. OK.
- A5 Dusart (arXiv 1002.0442 TeX, l.~699, theorem Rosser:th:lnPInvP): for x ≥ 2974, Σ_{p≤x} ln p/p − ln x − E ≤ 0.2/ln x + 0.2/ln²x, E = −γ − Σ_{n≥2}Σ_p ln p/p^n. Matches D:242-245. OK.

## Code audit cert.py (C:<line>)
- C:46-57 core: B > B₂ + 1e-17 certified; R(5/7;9) > 0, 5/7 > 2π/9, ℓ(5/7;9) > 0 certified. OK.
- C:64-157 zsum_coeffs: init za = a(5/7)n₀, zb = b(5/7)n₀, ze = 0 (= N⁺(5/7)G(0)); pieces integrated from max(lo,5/7), ze from max(u,2π−2); antiderivatives Pa, Pb, Pe = Lemma Z (iv), D:91-94 (rechecked by differentiation, A4). Unbounded piece only for m = 2 (else abort); envelope 3 vs 2 pieces decided on balls, undecided → abort. OK. NOTE: when lo > 5/7 is undecided the code integrates from 5/7: error O(radius²) ≈ 1e-110, immaterial.
- C:161-190 zsum_zeta: N⁺_ζ(14)G(14) + ∫_{14}^∞ ρ_ζ G; the 0.886/(t^{m+1} log t) parts bounded with log t ≥ log u (D:103-104). The last piece is always m = 2, so the unbounded branch is safe. OK.
- C:192-194 SUM_INV_G2 = 2(1+B²/196)(1+γ/2−½log4π). OK.
- C:199-227 Tri: n₀, n₁, n₂, W1 = I(b), W0 = I(−B), m₀ = I(1/2−B), J1, J2, m₁, m₂, ‖w‖_∞ = D:108-116; ks(Y) = Lemma K1. OK.
- C:229-238 omega_parts: L → max(L, log 9); Ω₁ ≤ x^B z1 + x^{1/2} z1L, z1 = SUM_INV_G2·k₂/2. OK.
- C:242-257 linnik_cell: one x per cell from L₀ (x ≤ (c₂ q log q)^{1/b} for q ≥ e^{L₀}), Z at L₁, q ≤ e^{L₁}: (4.1) exactly. OK.
- C:259-287 tail: k's frozen at Y_m. Simpler than D:158 says: G ≤ G_m pointwise (k's decrease in Y) and Lemma Z is applied to G_m itself, so no monotonicity of Z in G is needed. F(q) terms = D:156-157. OK.
- C:289-313 Thm 2': x_lo from f(q) ≤ φ(q); Ω₁ term divided by f(q)L; all terms at x_lo are upper bounds. OK.
- C:335-352 linnik(): each table row is certified by the TAIL at L_m = floor₆(log 10^e) ≤ log 10^e (no cells used); c₂ = 1 from q ≥ 100 (L_m = 4.605 > log 9). q < 100 → sanity.py sieve. OK.
- C:355-367 nonres_cell and C:369-382 nonres_tail: Ω₁ term is (z1 + x^{1/2−B}z1L)/L with NO 1/q (correct for Thm 3); 2·triv, 2·D_q. Check T3-Ω₁: PASSED.
- C:396-406 trivial_nonres: cells of width 1/100 (l₀ ≤ 20) / 1/1000, (1/b)log l₀ ≥ log(e^{l₁}−1), contiguous from log 3 past 25.6; tail rows start at L = 25 < 25.6 (overlap). OK.
- C:425-429 phi_lo: max(e^{L₀/2}/√2, e^{L₀}/(e^γ log L₁ + 2.50637/log L₀)); both are valid lower bounds, so whichever branch is taken is valid. OK.
- C:431-440 trap_consts = Lemma K2 with α^B ≤ al^B, al = e^{−X₀/2}. OK.
- C:442-479 pnt_cell: T1 = 1.03883 e^{(1/2−B)X₀}; T2 = δ(e^{bX₁}+e^{(1/2−B)X₀})/(2φ_lo); T3 = (Σ_ζG + e^{(1/2−B)X₀}Σ_ζG_L)/(2φ_lo); T4 = (Z+…)/2 at max(L₁, log 9); T5 = (X₁/2+lg)/(1−(e^{X₀/2}(1−δ))^{−2})/e^{BX₀}; T6 = ((L₁/log2)(X₁/2+lg)+L₁)/e^{BX₀} (q2: X₁/2+lg+log2; q1: 0); den (L₀+X₀)². δ fixed per cell (rational, ≤ 1/2). Every term is a valid cell-wide upper bound. OK.
- C:556-592 coverage: X-slabs from 4.16 (< log 66) with X₁ = X₀ + max(1/20, X₀/50) until X₀ ≥ 60; in each slab q1, q2 and L-cells from log 3 with L₁ = L₀ + max(1/10, L₀/20) while L₀ < X₁. Contiguous, no gaps; the cells also contain some (q,x) with q > x, which is harmless (6.3 never uses q ≤ x). OK.
- C:519-546 exact_small: Λ by trial division; ψ(·;q,a) accumulated; pieces of [n,n+1) for n ≥ max(2,q), n ≤ 65 (x < 66). For q = 1, a % 1 = 0 = n % 1. OK.
- C:548-554 x<q case = D:198-199. OK.
- C:481-514 pnt_tail: N₀ = (1+δ_T)^B/B ≥ W_B(δ) (mean value), K₁, t_a ≥ 1/N₀ ≥ 5/7, envelope 3-piece (A₁² ≤ W_B A₂ since W_B ≥ 1/B, A₂ ≥ 2/δ_T), ρ₀(K₁B) < 0 and ρ₀ increasing ⇒ ρ₀ < 0 on [5/7,t_a]. lam/XT = lamX, so the z_a/(4X) and log9·z_a/X² terms are the values at XT of non-increasing functions. be0 includes z_e (z_e/(1+L) ≤ z_e). λ-terms K₁λ²/π + 2K₁λ/π (z_b plus Σ_ζ). fac = (1+M x^{1/2−B})/2 with G_L ≤ M G pointwise (W_h ≤ 2(1+δ)^{1/2}, A₁^L/A₁ ≤ (1+δ)^{1/2}(1+α^{1/2}), A₂^L/A₂ ≤ ((1+δ)^{3/2}+1)(1+α^{1/2})/2). T1, T2 (φ ≥ 1, δ = θe^{−bX}), T5, T6 (L/(L+X)² ≤ 1/(4X)) OK. Matches D:214-224. OK.
- C:596-659 Thm 4: lp_kernels = D:239; lp_value: x = max(2974,(κL₀)^{1/b}) (any x valid), Z at max(L₁, log 9), prime part X + log(1+h) − γ + 0.2/X + 0.2/X², trivial 1/(x−1); cell constant = value − (1/b)log L₀ (valid since loglog q ≥ log L₀). lp_tail: every term non-increasing in L (x^{−1/2}L decreasing as 1/(2b) > 1 certified). OK.
- Claimed exponent 7.8978 > 1/b = 7.897769 and coefficient 7.8978 in Thm 4: weaker than certified since all bases > 1 / loglog q > 0 for q ≥ 3. OK.

## Independent recomputation (indep.py, mpmath 256 bits, zero sums by quadrature; logs/indep_*.txt)
- const: B₂ = 0.873381966613207519929, B − B₂ = 8.0e-17 (> 2.2e-254), 1/b = 7.8977691664599; a(5/7) = 0.4547342044, b(5/7) = −1.162576012, R(5/7;9) = 1.552099875; Σ1/γ² bound 0.0463711864607 (half 0.0231855932); W(1) = 2.424245725, W(0) = 115.9708945, ‖w‖_∞ = 38.85438136, n₂ = 1.360544218; m₀ = 12.70765942, m₁ = 12.5260612 (true ‖g_L'‖₁ = 5.489 ≤ m₁), m₂ = 12.33420029; c_∞ = 0.453310949; q = 10⁴ display row (Y = 83.2768, k₁ = 2.015461, k₂ = 1.378488, z = 2.20587, −4.675627, 0.0577325); 72.15% share. All = cert/DERIVATION.
- lemz: decomposed Z_G ≥ Stieltjes sum with the undecomposed BMOR bound at L = log 9, log 10⁴, 50 (8.368 ≥ 6.630; 31.911 ≥ 31.409; 132.827 ≥ 132.729; left kernels likewise). OK.
- linnik: every table margin reproduced to the printed digits (c₂ = 1 @10²: 0.21221; 0.79 @10²: 0.002638; 0.75 @10³: 0.0022575; 0.72 @10⁴: 0.007455; …; 0.46 @10¹⁰⁰⁰: 0.004881); φ-versions likewise (min 0.00095914 at 0.79). Direct (non-frozen) evaluation on L-grids: minimum at the left end = tail margin. OK.
- nonres: c₃ = 1 @25: 0.38871; 0.62 @25: 0.0095833; … 0.46 @10⁴: 0.0068184 (all = cert). Trivial range: min over integers 3 ≤ q ≤ 2·10⁵ of (1/b)loglog q − log(q−1) = 0.049621 (q = 3), and 0.0092 at log q = 25.6; concave in log q in between. OK.
- pnt cells: 21 cells (3 largest, first/last q1, q2, last two gen at X ≈ 60, 12 random) recomputed: indep ≤ cert in all, |rel diff| ≤ 5.3e-4 (cert's log-t bound in the Rosser part and rational δ are slightly conservative). Max cell = q = 3, x ∈ [64, e^{4.2432}]: 0.17307 (cert 0.17309) < 0.174. OK.
- ptail: my re-derivation of 6.4: 0.0142041 (X ≥ 60), 0.011127 (100), 0.0078753 (1000), 0.0075898 (10⁶) = cert. Direct non-majorised T₁…T₆ at X = 60, 61, 80, 150, 1000 over L/X ∈ [0,1]: max 0.008495 (X = 60, L/X = 0.75) < 0.0142. OK.
- lprime: first Thm-4 cell (log q ∈ [log 3, log 3 + 0.05]) C = 10.5263 (cert 10.526 ± 3e-4), later cells smaller; claimed 10.53. OK.
- small: exact ψ for q ≤ 65, x < 66, 400/8 samples per unit: max R = 1.4847782 at q = 1, x = 2 (cert piece bound 1.48535 ≥); on [64,66): 0.00845 (cert 0.00849). OK.
- kcheck.py: closed forms of W_tri, W_± agree with quadrature to 1e-31; on 4000 random (β, γ ∈ [1e-3, 3e3], x) the ratios actual/majorant are ≤ 0.9920 (K1), 0.9920 (K1 left), 0.9691 (K2), 0.9891 (K2 left), 0.9117 (Thm-4 kernel), 0.9246 (left). OK.

## Falsification (logs/falsify_output.txt, .scratch sanity_small run)
- Exact least primes q = 2..300: worst log p/log((q log q)^7.8978) = 0.4259 (q = 2, p = 3); φ-version 0.3130 (q = 3, p = 7); q = 95..105 ratios 0.109–0.158 (boundary sieve/analytic: no tension).
- Non-residues q = 3..400: tightest log-margin 0.0496 at q = 3 (n = 2 ≤ 2.10).
- Exact ψ(x;q,a), q ≤ 10, all a, x at every integer and n+1⁻ up to 10⁶: max R 1.4848 (q = 1), 0.5681 (q = 2, x = 2), 0.1190 (q = 3), 0.0775 (q = 4); x ≥ 64: ≤ 0.00935; x ≥ e¹⁰: ≤ 0.000161.
- |L'/L(1,χ)| for primitive χ mod 3, 4, 5, 7, 8: 0.368, 0.2456 (= γ + 2log2 + 3logπ − 4logΓ(1/4)), 0.181, 0.016, 0.632 ≪ 11.3–16.3.
- sanity.py at reduced size (copy in .scratch made with sed: N = 10**8 → 10**7, q-ranges 10**4 → 1500, prime powers from p ≤ 4000; footprint 19 MB; printed labels '1e8', '1e4' are hard-coded strings): failures 0; 0.4259 (q = 2), 0.3130 (q = 3), 1.4848 (q = 1, x = 2), max n_χ = 37 (q = 1364).
- Real characters: nothing in the chain treats them separately; H (β ≤ B for every χ, incl. real) is the only thing excluding a Siegel zero. A single real zero β > B would break Lemma P (weight x^β > x^B). H itself: HPAudit PASSED (statement covers every modulus and χ, ζ as q = 1).
- HP141 θ₁₄₁ rows: HP141/ExplProg/cert.py differs from ExplProg/cert.py only in B (l.42) and the core check/print (l.47, 49); its log has no CERT FAIL; margins 0.2121, 0.007355, 0.001213, 10.525. OK.

## Reviewer (most delicate step: Lemma P + Lemma Z, BMOR use, low zeros at height 0)
- Verdict CORRECT (cosmetic: R(5/7;9) ≈ 1.5521). Notes: x > 0 suffices; Tonelli justifies (iii); a γ = 0 atom contributes exactly n₀; N⁺G = O(log t/t).

## 'Before' source (1802.00085 TeX, v3 = Illinois J. Math. 62 (2018))
- Cor. uniform-cor: x ≥ exp(8√q log³q) for all q ≥ 3 (for q ≥ 58: x ≥ exp(0.03√q log³q)); main thms: c₀ = 1/160, x₀ = exp(0.03√q log³q) for q > 10⁵. The page's 'exp(8√q log³q)' quotes the all-q corollary: accurate; the sharper 0.03 (q ≥ 58) could be cited. NOTE, no effect.
- 64-bit sample (logs/indep_pnt_prec64_sample.txt): 133 cells (9 extremes, top cells with X₁ > 10 and > 20, 120 random): indep ≤ cert up to the 10-digit print rounding (≤ 3e-10), |rel diff| ≤ 9.2e-4.

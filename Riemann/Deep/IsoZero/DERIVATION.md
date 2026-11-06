# DERIVATION — explicit constants for the isolated-zero theorem of Lamzouri (arXiv 2610.26841v2)

TeX line numbers refer to IsolatedZeros2.tex of 2610.26841v2 (the latest version, 26 Oct 2026).
"cert.py §k" refers to the section comments in cert.py. Every number quoted here is printed by cert.py (logs/cert.log).

## 0. Statement and notation

N(T) = #{ρ: 0 < γ ≤ T} with multiplicity. For c > 0, B_0(T,c) = #{n ≥ 2: 0 < γ_n ≤ T, min(γ_{n+1} − γ_n, γ_n − γ_{n−1}) > 2πc/log T}
(Thm 1.1, l.211–220; ordinates listed with multiplicity, so a counted zero is simple, and it is on the line because 1 − ρ̄ has the same ordinate).
L := log T/(2π). Z_T = {−i(ρ − ½)L : 0 < γ ≤ T} (l.577–599), a conjugation-invariant multiset; real elements ↔ zeros on the line.
A real x ∈ Z_T with |x − Re z| > τ for every other z ∈ Z_T is a simple zero on the line whose ordinate is more than 2πτ/log T from every other ordinate in (0,T];
hence B_0(T,τ) ≥ |I_1| − 2 (the top zero below T and n = 1 excluded), with I_1 as in §2.

Main result (cert.py, logs/cert.log): for every c in the table of REPORT.md, first in the form

  liminf_{T→∞} B_0(T,c)/N(T) ≥ p(c) := (2 − λ − C(f) − D·α̃(c,λ)) / (1 − λ),                                   (★)

with λ = λ_C(c) (cluster loss, §2.2), C(f) the pair-correlation energy of the kernel f (§1), α̃ the off-line deficiency slope (§2.5) and D the off-line first-moment constant (§3).
The headline column H uses the hinge form (§3.7): p(c) = (2 − λ − C(f) − α_ρ B(ρ))/(1 − λ), where the off-line term α_ρ B(ρ) is below 1e-5 for small c.

## 1. Kernels

1.1 Family. f_w(x) = cos(w x)/Z_w on [−½,½], Z_w = 2 sin(w/2)/w, 0 < w < π (so f_w > 0 on [−½,½], ∫ f_w = 1). w = √2 is the Montgomery–Taylor f_0 of l.601–606.
K_w = \hat f_w (paper's K_0 for w = √2). Closed forms (cert.py §1; checked against the w = √2 formulas of l.785):
 K_w(u)  = w [sinc((w+2πu)/2) + sinc((w−2πu)/2)] / (4 sin(w/2))           (u real; 2∫_0^{1/2} f_w cos(2πux) dx),
 K_w(iy) = w² [cosh(πy) + (2πy/w) cot(w/2) sinh(πy)] / (w² + 4π²y²)        (2∫_0^{1/2} f_w cosh(2πxy) dx; ∫_0^{1/2} cos(wx)cosh(bx) dx = [w sin(w/2)cosh(b/2) + b cos(w/2)sinh(b/2)]/(w²+b²)).

1.2 Energy. C(f) := ∫ f² + ∬ |u−v| f(u) f(v) du dv is the constant of Lem:ConstructionKernel (l.670–700; = La 2609.02882, Lemma 3.2, from BGST Lemma 5):
Σ_{z,s∈Z_T} K(z−s)² = (C(f) + o(1)) N(T) for f = η², η smooth, supp η ⊂ (−½,½). (Our ZetaZeros audit re-derived C(f) from BGST; Riemann/Deep/ZetaZeros/REPORT.md item 1.)
For f_w: ∫_{−1/2}^{1/2} |u−v| cos(wv) dv = sin(w/2)/w + 2(cos(w/2) − cos(wu))/w² (integrate by parts on each side of v = u), hence
 C(f_w) = I_w (1 − 2/w²)/Z_w² + (sin(w/2)/w + 2 cos(w/2)/w²)/Z_w,   I_w = ½ + sin(w)/(2w).
At w² = 2 the first term vanishes and C(f_0) = ½ + (1/√2) cot(1/√2) = C_MT = 1.3274992963… (l.672). cert.py also integrates ∫f_0² = 1.00612719… and ∬|u−v|f_0f_0 = 0.32137210… with Arb quadrature: sum overlaps C_MT.

1.3 Smoothing (ε → 0). BGST needs η ∈ C_c^∞((−½,½)). As in l.634–668 put f_ε = ψ_ε² f_w / A_ε. Since f_w > 0 and smooth on [−½,½], η_ε = ψ_ε √f_w/√A_ε is C_c^∞. As ε → 0: f_ε → f_w in L¹ and L², so C(f_ε) → C(f_w); K_ε → K_w uniformly on compact sets of ℂ (|K_ε(ξ) − K_w(ξ)| ≤ ‖f_ε − f_w‖_1 e^{π|Im ξ|}). Every quantity in (★) is continuous under this convergence (§4), so (★) is proved for f_ε with constants → those of f_w, and the liminf bound holds with the limits.

## 2. The Hilbert-space inequality with explicit errors (replaces Pro:IsolatedSMALL, l.511–573)

Notation of the proof of Pro:IsolatedZeros (l.283–468): R_1 simple real elements, I ⊂ R_1 those > 2τ from all other elements of R_1 (|I| = d), C = R_1 \ I (clustered; split into a pairs and b triples of consecutive elements ≤ 2τ apart, |C| = 2a+3b), R_2 real elements of multiplicity ≥ 2 (r distinct, total multiplicity M_R ≥ 2r), S = {z_ℓ, z̄_ℓ} (q distinct pairs, total multiplicity M_S ≥ 2q). N = d + |C| + M_R + M_S.
The paper's estimates (Eq:MainInequalityKernel l.390, Eq:LowerAlphaU l.439, Eq:LowerAlphaV l.443, Eq:DimensionCalculations l.366) give, before any substitution (re-derived; the paper writes it after substituting N, l.451):

 Σ_{z,s} K(z−s)² ≥ 2S_0 + d + 4aκ_0 + b(16κ_0/3 + 8κ_1/3) + 4(M_R + M_S) − 4(r+q),                      (E)

S_0 = Σ_{x∈I} s_x, s_x := Σ_{j≤D_U} ⟨f_x, ψ_j⟩² ≥ 0 (squared norm of the projection of f_x onto U), κ_0 = min_{[−2τ,2τ]} K, κ_1 = min_{[−4τ,4τ]} K.
(Steps: y² ≥ 4y−4 on j ≤ D_U and y² ≥ 2y−1 on D_U < j ≤ D_V give ΣK² ≥ 2Σ_U α + 2Σ_V α − 3D_U − D_V; D_U = a+b+r+q, D_V = D_U + d; Σ_U α ≥ S_0 + a(1+κ_0) + b(1+4κ_0/3+2κ_1/3) + M_R + M_S; Σ_V α ≥ d + same.)
Subtracting 2N = 2d + 4a + 6b + 2M_R + 2M_S:

 ΣK² − 2N ≥ (2S_0 − d) − 2a·2(1−κ_0) − 3b·β_T + 2(M_R − 2r) + 2(M_S − 2q),   β_T := (6 − 16κ_0/3 − 8κ_1/3)/3.        (E')

2.1 Monotonicity. For f ≥ 0 even, supported in [−½,½]: K'(u) = −4π∫_0^{1/2} x f(x) sin(2πux) dx < 0 on (0,1] (l.705–707), so for 4τ ≤ 1: κ_00 := min_{[0,τ]}K = K(τ), κ_0 = K(2τ), κ_1 = K(4τ).

2.2 Clusters. With ε_0 = 1 − κ_0 and λ_C := max(2ε_0, β_T, 0), every λ ≥ λ_C gives −4aε_0 − 3bβ_T ≥ −λ|C|.
(Small τ: K(u) = 1 − 2π²u²μ_2 + O(u⁴), μ_2 = ∫x²f; for f_0, μ_2 ≈ 0.0775, λ_C ≈ 12.2 τ²; this is the paper's O(τ²N), made explicit.)

2.3 Split of I. I_1 = {x ∈ I: |x − Re z| > τ for all z ≠ x} (the counted set), I_2 = I \ I_1. Each x ∈ I_2 has a witness in R_2 or in S within τ (in real part); since points of I are > 2τ apart, a witness serves at most one x (injective; l.330).
I_2R: x with some R_2 witness x_m. Then s_x ≥ ⟨f_x, f_{x_m}⟩²/‖f_{x_m}‖² = K(x − x_m)² ≥ κ_00², so 2s_x − 1 ≥ 2κ_00² − 1 ≥ 0 (cert.py asserts 2K(τ)² ≥ 1 for every τ used).
I_2S: x with only complex witnesses; pick one pair z = u + iy, z̄ (|x − u| ≤ τ). g_z ∈ U, so s_x ≥ ⟨f_x, g_z⟩²/‖g_z‖² (l.549).

2.4 Explicit complex-witness bound. ⟨f_x, g_z⟩ = ∫ f(t) cos(2πδt) cosh(2πyt) dt with δ = x − u, |δ| ≤ τ (the sine part is odd). Since f ≥ 0, cosh > 0, cos θ ≥ 1 − θ²/2 and |t| ≤ ½:
⟨f_x, g_z⟩ ≥ (1 − 2π²δ²/4) K(iy) ≥ (1 − θ̄) K(iy),  θ̄ := π²τ²/2  (> 0 is asserted).  ‖g_z‖² = (1 + K(2iy))/2 (l.565). Hence
 2s_x − 1 ≥ (1−θ̄)² · 4K(iy)²/(1+K(2iy)) − 1 = −H_τ(y),   H_τ(y) := 1 − (1−θ̄)² 4K(iy)²/(1+K(2iy)).
(For τ → 0 this is the paper's (Eq:LowerBoundProjectionComplex2) without the O(τ²).) Note H_τ < 1 and H_τ is even in y.

2.5 Combination. Insert d = |I_1| + |I_2| and 2S_0 − d ≥ −|I_1| + Σ_{I_2R}(2κ_00² − 1) − Σ_{x∈I_2S} H_τ(y_{m(x)}) into (E'):
 |I_1| ≥ 2N − ΣK² − λ|C| − Σ_{x∈I_2S} H_τ(y_{m(x)}) + (nonnegative terms).
Budget: N = |I_1| + |I_2R| + |I_2S| + |C| + M_R + M_S and M_S ≥ 2q ≥ 2|I_2S| (injective witness pairs), so |C| ≤ N − |I_1| − 3|I_2S|. As λ ≥ 0:
 |I_1| ≥ 2N − ΣK² − λ(N − |I_1|) − Σ_{x∈I_2S} (H_τ(y_{m(x)}) − 3λ)
       ≥ 2N − ΣK² − λ(N − |I_1|) − Σ_{distinct pairs ℓ ≤ q} (H_τ(y_ℓ) − 3λ)^+.
Define α̃(τ,λ) := sup_{y>0} (H_τ(y) − 3λ)^+ / y. A pair ℓ comes from a distinct zero with β > ½ and |y_ℓ| = (β − ½)L, so
 Σ_ℓ (H_τ(y_ℓ) − 3λ)^+ ≤ α̃ L Σ_{0<γ≤T, β>1/2} (β − ½)   (counted with multiplicity ≥ distinct).
With Lem:ConstructionKernel (ΣK² = (C(f)+o(1))N) and Σ(β−½) ≤ (D + o(1))T (§3), using N(T) ~ (T/2π) log T = T·L:
 (1 − λ)|I_1| ≥ (2 − λ − C(f) − D α̃ + o(1)) N(T), which is (★).

2.6 Checks of (★). (i) τ → 0: λ → 0, θ̄ → 0, α̃ → α_0 := sup H_0(y)/y: (★) → 2 − C − Dα_0, the paper's l.816 constant. (ii) λ = 1/3 (κ_0 = 5/6, i.e. τ = τ_0 = 0.169641840…): α̃ = sup (H − 1)^+/y = 0 and (★) = 5/2 − 3/2 C = 0.5087510555, the paper's Pro:IsolatedZeros / (1.6). So (★) interpolates between the two results of Thm 1.1. (iii) The |C| budget is what makes the cluster loss λ(1 − p) instead of λ, and the 3λ shift in α̃ is the alternative use of the same budget by a witness pair; both are lower-order but explicit.
(iv) No other term depends on c: the Lem:ConstructionKernel constant and D are c-free; c enters only through λ_C (via K(2τ), K(4τ)), θ̄ and the assertions 2K(τ)² ≥ 1, 4τ ≤ 1.

2.7 α̃ certificate (cert.py §2). On (0, 0.1]: K(iy) ≥ ∫f = 1 and K(2iy) ≤ cosh(2πy) (l.802–807) give H_τ(y) ≤ 1 − (1−θ̄)² 4/(1+cosh(0.2π)) < 0 (asserted). On [60, ∞): (H_τ − 3λ)^+/y ≤ (1 − 3λ)^+/60 < 0.017. On [0.1, 60]: branch and bound; on each box [a,b] with midpoint m, g(y) = (H_τ(y) − 3λ)/y satisfies g(y) ≤ g(m) + |g'(m)| r + max(g''(ξ)/2, 0) r² (Lagrange remainder; g(m), g'(m) from an Arb Taylor series at the exact point m, g''(ξ)/2 enclosed by the Arb Taylor series over the box), r = max(m−a, b−m). Boxes are split until the bound is within 1e-11 of the best point value. Output: α_0 ∈ [0.214612853119, 0.214612853122] (paper: < 0.21462), maximiser y ≈ 2.42.

## 3. The off-line first moment D (replaces Lem:Conrey, l.757–767)

3.1 Input (Conrey, J. reine angew. Math. 399 (1989), Theorem 2; quoted in Bui–Conrey–Young, arXiv 1002.4127, l.338–345). Let 0 < θ < 4/7, P a real polynomial with P(0) = 0, P(1) = 1,
R > 0 fixed, σ⁻ = ½ − R/log T, y = T^θ, A(t) := Σ_{n≤y} μ(n) P(log(y/n)/log y) n^{−1/2−it}. Then (Q ≡ 1, V = ζ; ψ_1(σ⁻+it) = A(t), BCY l.233–235)
 ∫_1^T |ζ(σ⁻+it) A(t)|² dt ~ c_1 T,   c_1 = 1 + (1/θ) ∫_0^1 e^{2Rv} dv ∫_0^1 (P'(u) + θR P(u))² du.
(Q ≡ 1 satisfies the side conditions Q(0) = 1, Q(x)+Q(1−x) = const of BCY l.206.)
Independent check of this formula: our own derivation from the twisted second moment (main terms ζ(1+α+β)(h')^{−1/2−β}(k')^{−1/2−α} + (t/2π)^{−α−β} ζ(1−α−β)(h')^{−1/2+α}(k')^{−1/2+β}) gives
c(A,B) = [∫(P'+AθP)(P'+BθP) − e^{−(A+B)}∫(P'−AθP)(P'−BθP)] / ((A+B)θ) at shifts A/L, B/L; at A = B = −R this equals c_1 above, it gives 1 + 1/θ at R → 0 (Wu's value, l.1000–1003),
and it satisfies the functional-equation symmetry c(A,B) = e^{−(A+B)} c(−B,−A).

3.2 Transfer to the right of the line. Put σ⁺ = ½ + R/log T and ψ⁺(s) = Σ_{n≤y} μ(n)P(log(y/n)/log y) n^{σ⁺−1/2} n^{−s}, so ψ⁺(σ⁺+it) = A(t), coefficient of 1^{−s} equal to P(1) = 1.
Functional equation: |ζ(σ⁺+it)| = |χ(σ⁺+it)| |ζ(σ⁻+it)| and |χ(σ⁺+it)|² = (t/2π)^{−2R/log T}(1 + O(1/t)) (t ≥ 1). Fix δ > 0 with θ/(1−δ) < 4/7.
For T^{1−δ} ≤ t ≤ T the factor is ≤ e^{−2R(1−δ)}(1+o(1)); for 1 ≤ t ≤ T^{1−δ} it is ≤ (2π)^{2R/log T} = 1 + o(1). Hence
 ∫_1^T |ζψ⁺(σ⁺+it)|² dt ≤ (1+o(1)) [ e^{−2R(1−δ)} ∫_1^T |ζ(σ⁻+it)A(t)|² dt + ∫_1^{T^{1−δ}} |ζ(σ⁻+it)A(t)|² dt ].
The first integral is ~ c_1 T. The second is Conrey's theorem at height U = T^{1−δ} with the same Dirichlet polynomial A (it depends only on y and P),
mollifier exponent θ_U = θ/(1−δ) < 4/7 and shift R_U = R(1−δ) (σ⁻ = ½ − R_U/log U), so it is O(U) = o(T). Letting δ → 0:
 limsup (1/T) ∫_1^T |ζψ⁺(σ⁺+it)|² dt ≤ e^{−2R} c_1 =: c_+(P,R,θ).

3.3 Littlewood. Apply Littlewood's lemma to f = ζψ⁺ on [σ⁺, 2] × [1, T] (no zeros of ζ with γ ≤ 1): ζψ⁺(2+it) − 1 → 0 uniformly (ψ⁺(2+it) itself tends to 1/ζ(2+it), not to 1; reviewer correction) (coefficients b_m = Σ_{d|m, d≤y} μ(d)P[d] d^{R/log T} → 0 for each fixed m ≥ 2, dominated by d(m)m^{−2+o(1)}), so ∫ log|f(2+it)| = o(T) and f has no zeros with β ≥ 2;
the horizontal arguments are O(log T) (Titchmarsh §9.4 for ζ and for the Dirichlet polynomial). Zeros of f include those of ζ. With Jensen (AM–GM):
 ∫_{σ⁺}^{1} ν(σ,T) dσ ≤ (T/4π) log c_+ + o(T),   ν(σ,T) := #{ρ: β > σ, 0 < γ ≤ T}.
Layer cake: Σ_{β>1/2}(β − ½) = ∫_{1/2}^1 ν(σ,T)dσ ≤ (R/log T) ν(½,T) + ∫_{σ⁺}^1 ν(σ,T) dσ.
If liminf N_0(T)/N(T) ≥ κ (N_0: zeros on the line with multiplicity), then ν(½,T) = (N(T) − N_0(T))/2 ≤ ((1−κ)/2 + o(1)) (T/2π) log T, so
 Σ_{0<γ≤T, β>1/2} (β − ½) ≤ (T/4π) [R(1−κ) + log c_+(P,R,θ)] + o(T).                                           (D)

3.4 Optimal P. ∫_0^1(P' + bP)² = ∫P'² + b²∫P² + b (as P(0) = 0, P(1) = 1); the minimum of ∫P'² + b²∫P² is b coth b (Euler–Lagrange P'' = b²P, P = sinh(bx)/sinh b),
approached by polynomials in C¹[0,1] with the same end values. With b = θR:
 inf_P c_+ = c*(R) := e^{−2R} + ((1 − e^{−2R})/2) (1 + coth(θR)) = 1 + (1 − e^{−2R})/(e^{2θR} − 1).
So D(κ) := (1/4π) inf_{R>0} [R(1−κ) + log c*(R)] (θ → 4/7) is an admissible constant in Σ(β − ½) ≤ (D + o(1))T.
Remarks. c* ≥ 1 and c*(R) → 1 + 1/θ as R → 0+: the R → 0 limit is (1/4π) log(1 + 1/θ) → (1/4π) log(11/4) = 0.0805006428 (cert.py), which is Lamzouri's 0.0806 (Lem:Conrey; Conrey's announced Theorem 3, Bull. AMS 1989, rounded up) and BCHB's 0.0845 = (1/4π) log(1 + 17/9).
So Lem:Conrey follows from Conrey's published Theorem 2 by this argument (R → 0+); the shift R > 0 is what improves it. d/dR [R(1−κ) + log c*] at R = 0+ equals −κ, so any κ > 0 helps.

3.5 The on-line input κ, from the same inequality. Dropping S_0 ≥ 0 and the nonnegative multiplicity terms in (E') gives d + λ|C| ≥ 2N − ΣK², so for λ ≤ 1
the number of simple zeros on the line satisfies d + |C| ≥ (2 − C(f) − o(1))N (this is Lamzouri's 2609.02882 argument), and the off-line multiplicity satisfies
 M_S = N − d − |C| − M_R ≤ ΣK² − N − (1−λ)|C| ≤ (C(f) − 1 + o(1)) N,   i.e. κ = 2 − C(f) = 0.6725007 is admissible in (D).
Values (cert.py): κ = 0: D = (1/4π)log(11/4) = 0.0805006 (R → 0; no gain, since dF/dR(0+) = −κ); κ = 5/12 (Pratt–Robles–Zaharescu–Zeindler 2020): R = 0.62496, D ≤ 0.0704984, p(0.001) ≥ 0.657368;
κ = 2 − C_MT: R = 1.16435, D ≤ 0.0526935, p(0.001) ≥ 0.661189. (An earlier version of this work bootstrapped κ from the isolated proportion itself, reaching D ≤ 0.0537484 and 0.660962; superseded.)

3.6 Joint bound (column C; superseded by §3.7). The same inequality gives more: (E') with 2S_0 ≥ Σ_{x∈I_2S} 2s_x ≥ Σ_{x∈I_2S} (1 − H_τ(y_x)) (§2.4; y_x the height of x's witness pair) yields
 M_S ≤ ΣK² − N − Σ_{x∈I_2S} (1 − H_τ(y_x))   (λ ≤ 1).
Insert this in (D) written with the actual off-line fraction (ν(½,T) = M_S/2; L·(R/log T)·M_S/2 = R M_S/4π; L·(T/4π) log c_+ = (N/4π)(1+o(1)) log c_+), and use Σ_{x∈I_2S} y_x ≤ L Σ_{β>1/2}(β − ½):
 Σ_{x∈I_2S} [ y_x + (R/4π)(1 − H_τ(y_x)) ] ≤ D(R) N + o(N),   D(R) := (1/4π)[R(C(f) − 1) + log c*(R)].
Termwise (H − 3λ)^+ ≤ α̂ (y + (R/4π)(1 − H)) with α̂(τ,λ,R) := sup_{y>0} (H_τ(y) − 3λ)^+ / (y + (R/4π)(1 − H_τ(y))) (denominator > 0 as H_τ < 1). With §2.5:
 (1 − λ) |I_1| ≥ (2 − λ − C(f) − α̂(τ,λ,R) D(R) − o(1)) N,   for every R > 0;
cert.py optimises R per c (R ≈ 1.19 at small c). At c = 0.001: R = 1.1914, D(R) ≤ 0.05270413, α̂ ≤ 0.21070014, p ≥ 0.661391792 (f_0).
The witness pairs that create the off-line deficiency also lower the off-line budget, which is why α̂ < α̃ (by ≈ 1.8% at small c).

3.7 Hinge bound (the headline column H). Use §3.2–3.3 at σ⁺ = ½ + R/log T, R = 2πρ for a fixed ρ > 0, but keep only the part of the layer cake above σ⁺:
 Σ_{β>σ⁺} (β − σ⁺) = ∫_{σ⁺}^1 ν(σ,T) dσ ≤ (T/4π) log c*(R) + o(T).
In the normalised height y = (β − ½)L (L = log T/2π) one has β − σ⁺ = (y − ρ)/L and T·L ~ N(T), so
 Σ_{0<γ≤T, β>1/2} (y − ρ)^+ ≤ (B(ρ) + o(1)) N(T),   B(ρ) := (1/4π) log c*(2πρ) = (1/4π) log(1 + (1 − e^{−4πρ})/(e^{(16π/7)ρ} − 1)).
No information on the proportion of zeros on the line is needed. B(ρ) ≈ e^{−7.18ρ}/(4π): B(0.5) ≤ 2.223e-3, B(1) ≤ 6.059e-5, B(1.2) ≤ 1.441e-5, B(2) ≤ 4.610e-8 (cert.py prints 2.222096e-03, 6.058700e-05, 1.440588e-05, 4.609366e-08).
(In words: zeros farther than 2πρ/log T to the right of the line are exponentially rare; this is the integrated form of the density estimate behind Conrey's N(σ,T) ≪ T^{1−8(σ−1/2)/7} log T.)
Deficiency. h(y) := (H_τ(y) − 3λ)^+ vanishes for 0 < y ≤ ρ_max(τ) (certified: H_τ − 3λ < 0 on (0, ρ_max]; ρ_max ≈ 1.22 for small τ, larger for larger τ).
For 0 < ρ ≤ ρ_max put α_ρ := sup_{y>ρ} h(y)/(y − ρ) (certified branch and bound; finite because h(ρ) = 0 is approached from below). Termwise h(y) ≤ α_ρ (y − ρ)^+, so
 Σ_{x∈I_2S} (H_τ(y_x) − 3λ) ≤ Σ_{distinct pairs ℓ} h(y_ℓ) ≤ α_ρ (B(ρ) + o(1)) N,
and with §2.5:
 (1 − λ)|I_1| ≥ (2 − λ − C(f) − α_ρ B(ρ) − o(1)) N.                                                         (H)
At c = 0.001: ρ = 1.2036, α_ρ ≤ 0.670124, B(ρ) ≤ 1.404e-5, off-line term ≤ 9.5e-6, p ≥ 0.672487 — within 1.4e-5 of 2 − C_MT = 0.672501, the simple-and-on-line proportion of the same Hilbert-space method (La 2609.02882).
Reading: an off-line zero at normalised height y < ρ_max is "rigid" (H_τ(y) < 0: it is charged like a double zero and cannot spoil an isolated point), and higher off-line zeros are exponentially rare. The first-moment bounds of §3.4–3.6 are the special case ρ → 0 of this family plus a count of off-line zeros; they are superseded.

## 4. Limits and quantifiers

For fixed c the inequality (★) is proved for every admissible (ε, θ, P): ε ∈ (0, 1/10) (kernel smoothing), θ < 4/7, polynomial P. For each such choice,
liminf_T B_0(T,c)/N(T) ≥ p(c; ε, θ, P). The right side is continuous and its limit as ε → 0, θ → 4/7, P → sinh(bx)/sinh b is the certified p(c);
since the liminf is at least the sup over admissible choices, it is at least the limit. Continuity: C(f_ε) → C(f_w); K_ε(u) → K_w(u) on [0,4τ] (so λ_C converges and stays ≥ 0);
α̃_ε → α̃ because H_{τ,ε} → H_τ uniformly on [0.1, 60], H_{τ,ε} < 0 on (0, 0.1] for every ε (bounds K_ε(iy) ≥ 1, K_ε(2iy) ≤ cosh 2πy hold for every ε), and the tail bound (1−3λ)^+/60 is ε-free; c_+ → c* and D(θ) → D(4/7). Hinge form: H_{τ,ε} − 3λ_ε → H_τ − 3λ uniformly on [0.1, ρ], where the limit is ≤ −η < 0 (certified), so the sign condition holds for all small ε; α_{ρ,ε} → α_ρ for the same reasons as α̃, and B(ρ; θ) → B(ρ; 4/7) as θ → 4/7.
The only c-dependence is listed in §2.6(iv). Finite T: all statements are liminf statements; no claim is made for a specific T.

## 5. Kernel optimisation (cert.py "kernel f_w")

The paper's Remark 1.2 (l.222–227) suggests optimising the kernel. Within f_w = cos(wx)/Z_w (per-c golden section on w², certified at the chosen w):
hinge bound (H): at small c the off-line term is negligible, so f_0 (which minimises C(f)) is optimal: w² → 2.00 at c = 0.001 and 0.01 (gain < 1e-7). Larger c favours peaked kernels (smaller μ_2, smaller λ): c = 0.05, w² = 2.1257: 0.662237 (f_0: 0.662211); c = 0.1, w² = 2.5315: 0.628354 (f_0: 0.627814); c = 1/6, w² = 3.75: 0.524702 (f_0: 0.516706; paper: 0.50875).
Paper inputs (D = 0.0806): c = 0.001, w² = 1.7547: 0.655296 (f_0: 0.655200), the kind of gain Lee's forty-mode kernel exploits. Each optimised value is certified like the f_0 values (§2 with K_w, C(f_w)).

## 6. Rounding (cert.py) and an independent float check

sanity_mp.py (mpmath, K by quadrature, sup on a grid of step 5e-4) reproduces the c = 0.001 hinge numbers: ρ_max = 1.22365, sup h/(y − ρ) = 0.6701238 at ρ = 1.2036, B = 1.40382e-5, p = 0.6724873 (logs/sanity.log).


Arb balls at 160 bits throughout. Float conversions use one-ulp directed rounding (flo, fhi). λ ≥ λ_C (upper endpoint). α̃: upper endpoint of the B&B (with 1e-14 relative padding of float bound arithmetic). D: κ rounded down, F rounded up.
p: lower endpoint of the Arb expression (2 − λ − C − D·α̃)/(1 − λ). D(R), B(ρ): log c* rounded up. α̂, α_ρ: the same branch and bound (denominators y + (R/4π)(1 − H), y − ρ; boxes with H − 3λ ≤ 0 are discarded first, so y − ρ is never 0 on a used box). ρ_max: bisection, each candidate certified by interval enclosures of H on boxes of width ≤ 1e-7. C(f_w): exact closed form in Arb.

## 7. Remarks

7.1 log T versus log γ_n. Thm 1.1 uses 2πc/log T; the abstract uses 2πc/log γ_n. For γ_n ∈ (T^{1−δ}, T], 2πc/log γ_n ≤ 2πc/((1−δ) log T), so B_0^{(log γ)}(T,c) ≥ B_0(T, c/(1−δ)) − N(T^{1−δ}).
Since B_0 decreases in c and N(T^{1−δ}) = o(N(T)), the table value at c holds for the log γ_n definition at every c' < c (take δ with c'/(1−δ) ≤ c).
7.2 Multiplicity. N(T), Lem:ConstructionKernel and Σ(β−½) count with multiplicity; the counted isolated zeros are simple; a non-real or multiple witness is charged by its multiplicity through M_R, M_S ≥ 2r, 2q.
7.3 What is not used: RH, GRH, Montgomery's conjecture, the well-spacing hypothesis. Inputs: BGST (Acta Arith. 2024, refereed), La 2609.02882 Lemma 3.2 = Lem:ConstructionKernel (2026 preprint, also used by the paper), Conrey 1989 Theorem 2 (refereed), classical Littlewood lemma and the functional equation.

## 8. Side result: Dirichlet L-functions and Theorem 1.4 (same footing as the paper's Theorem 1.4)

The paper extends (1.5) to a fixed primitive χ (l.966–1032) with Lem:GaPa (dyadic pair correlation; same constant C(f)) and Lem:AnalogueConreyL (0.0806, via Wu's mean value on the line).
The hinge lemma transfers: Wu, Math. Z. 293 (2019), arXiv 1802.09704, Theorem 1.1(B), gives I_R(Q,χ) = ∫_T^{2T} |V(½ + it + R/L_χ, χ) B(½+it, χ)|² dt asymptotically for θ < 4/7,
fixed shifts a = b = −R (a, b ∈ ℂ, a, b ≪ 1), mollifier coefficients μ(n)F_0(n) with F_0(x) = P(log(y/x)/log y) admissible; with Q = 1 its main term has the structure of §3.1
(L(1+α+β, χ_0) and (qT/2π)^{−α−β}-weighted L(1−α−β, χ_0) terms, the q/φ(q) factors cancel against the Euler products), so c_1 is unchanged.
(Wu's displayed (1.10) shows (h,k)^{α+β} where (1.12) and the α = β = 0 case need (h,k)^{1+α+β}; we read it as a typo; this is the only check made.)
The functional equation |L(½+a+it, χ)| = (qt/2π)^{−a} |L(½−a+it, χ)| (1 + O(1/t)) moves it to the right of the line; on dyadic blocks t ∈ [T, 2T] the factor is ≤ e^{−2R}(1+o(1)) directly.
Hence Σ_{T<γ≤2T, β>1/2} (y − ρ)^+ ≤ (B(ρ) + o(1)) (T/2π) log T for L(s,χ), and §2 applies verbatim to the dyadic multisets Z_T^{(j)} (l.1018–1025).
So each L(s,χ_j) has ≥ 0.672487 of its zeros in [T,2T] 0.001-isolated, and Remark 1.5 / Pro:BombieriHejhal (l.251–252, 868–870, not re-audited here) give
 ≥ 2·0.672487 − 1 = 0.344974 of the zeros of b_1e^{iθ_1}L(s,χ_1) + b_2e^{iθ_2}L(s,χ_2) on the line with odd multiplicity (paper: 31%; Bombieri–Hejhal under GRH: 34.5%).
Status: derived, not separately certified or reviewed; depends on Wu's theorem as stated and on the paper's Pro:BombieriHejhal.

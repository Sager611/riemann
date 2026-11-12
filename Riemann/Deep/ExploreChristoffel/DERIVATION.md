# explore38 derivation: Burnol's zero vectors at finite λ are Christoffel–Darboux kernels of a de Branges-type space

Notation: K = L²(0,∞), Mellin on Re s = 1/2 with dτ/2π; (f,g) linear in f. χ = 1_(0,1], χ₁ = (1 + log t)χ (Mellin 1/s − 1/s²), M f(t) = t⁻¹∫₀ᵗ f (multiplier 1/(1−s)),
U = Báez-Duarte operator (multiplier (s/(1−s))·ζ(1−s)/ζ(s)), V = (1 − M)²U, multiplier V(s) = (s/(1−s))³/χ(s) with ζ(s) = χ(s)ζ(1−s); |V| = 1 on the line,
so V is unitary. ψ_w = t^{−w}χ (in K for Re w < 1/2), Q_λ = projection onto L²[λ,∞), L = log(1/λ). A(t) = inverse Mellin of Z(s) = (s−1)ζ(s)/s² (Burnol Prop.),
C_λ = span{D_θ A : λ ≤ θ ≤ 1}; D(λ) = dist(χ₁, C_λ) is BBLS's continuous Nyman–Beurling distance (Burnol §3), and d_N ≥ D(1/N) (dilations 1/k ∈ [1/N, 1]).

## 1. Burnol (math/0103058, Thm-Def + Thm theo + Cor)
Y_ρ := lim_{w→ρ, Re w<1/2} V⁻¹Q_λVψ_w exists for Re ρ = 1/2, and (D_θA, Y_ρ) = θ^{ρ−1/2}Z(ρ) for θ ∈ [λ,1]; so Y_ρ ⊥ C_λ iff ζ(ρ) = 0.
Hence for any finite set R of on-line zeros: D(λ)² ≥ ‖P_{span Y_R} χ₁‖² = cᴴH⁻¹c, H_ij = (Y_j, Y_i), c_i = (χ₁, Y_i). Burnol takes λ → 0 only
(Gram → identity after rescaling, c_ρ → (ρ−1)/ρ²), giving liminf D(λ)² L ≥ Σ m_ρ²/|ρ|². Here: everything at finite λ.

## 2. Pointwise form of Vψ_w (boundary term!)
Burnol's φ_{w,0} = Uψ_w is ∫_δ¹ v^{−w} (d/dv) k(t/v) dv with k(x) = sin(2πx)/(πx). Integrating by parts gives
φ_{w,0}(t) = k(t) + w∫₁^∞ u^{w−1}k(tu)du. (Burnol's Thm pretheo is stated for k ≥ 1, where the boundary term vanishes; the k = 0 case needs
the extra k(t). His final series extends correctly to k = 0. Reviewer: correct as a k = 0 extension, not an erratum.) Mellin check: k̂(s) = γ₊(s)/(1−s) and k̂(s)(1 + w/(s−w)) = U(s)/(s−w). Hence for every t > 0
  Vψ_w(t) = V(w)t^{−w} + F_w(t),  F_w(t) = Σ_{j≥1} a_j(w)(2πt)^{2j},  a_j(w) = 4j(−1)^j(2j/(2j+1))² / ((2j+1)!(w+2j)),
i.e. the residues of V(s)t^{−s}/(s−w) at the poles s = −2j of V (Res_{−2j}V = −2(−1)^j(2π)^{2j}(2j/(2j+1))³/(2j)!; checked equal).
At w = −2j the two displayed pieces diverge separately; the combined expression is read as its removable limit.
Similarly Vχ₁ = (1−M)k = Σ_{j≥1} b_j(2πt)^{2j}, b_j = 2(−1)^j(2j/(2j+1))/(2j+1)!.
(Without the k(t) term the Gram matrix below is not positive definite at N = 2: eigenvalue −4.65; with it λ_min(H) ≥ 1.04 for every N ≥ 2 tested.)

## 3. Exact Gram at finite λ
V unitary ⇒ (Q_λVψ_w, Q_λVψ_w') = (ψ_w,ψ_w') − ∫₀^λ Vψ_w·conj(Vψ_w') dt, (ψ_w,ψ_w') = 1/(1 − w − w̄'). All integrals are of monomials:
  (Y_w,Y_w') = [1 − V(w)conj V(w') λ^{1−w−w̄'}]/(1 − w − w̄') − V(w)I₁ − conj V(w') I₂ − I₃,
  I₁ = Σ_j conj(a_j(w'))(2π)^{2j} λ^{2j+1−w}/(2j+1−w),  I₂ = Σ_j a_j(w)(2π)^{2j} λ^{2j+1−w̄'}/(2j+1−w̄'),
  I₃ = Σ_{j,k} a_j(w)conj(a_k(w'))(2π)^{2j+2k} λ^{2j+2k+1}/(2j+2k+1).
On-line limits (w → ρ = ½+iγ, w' → ρ' = ½+iγ'): 1 − ρ − ρ̄' = i(γ'−γ), λ^{i(γ'−γ)} = N^{i(γ−γ')}; V(½+it) = ((½+it)/(½−it))³e^{2iθ(t)}
(θ = Riemann–Siegel theta; χ(½+it) = e^{−2iθ(t)}). Diagonal: with w = ρ − ε/2, |V(w)|² = 1 − ε Re(V'/V)(ρ) + O(ε²), λ^ε = 1 − εL, so
  (Y_ρ,Y_ρ) = L + Re(V'/V)(ρ) − 2Re(V(ρ)I₁) − I₃,  Re(V'/V)(ρ) = 3/|ρ|² − χ'/χ(ρ) = 3/|ρ|² + 2θ'(γ).
Writing Φ_N(t) := t log N + 2θ(t) + 6 arctan(2t) (so e^{iΦ_N(t)} = V(½+it)N^{it} and Φ_N' = L + 2θ' + 3/(¼+t²)):
  G_{ρρ'} := (Y_ρ, Y_ρ') = (1 − e^{i(Φ_N(γ)−Φ_N(γ'))})/(i(γ'−γ)) + O(N^{−5/2}),  G_{ρρ} = Φ_N'(γ) + O(N^{−5/2}),
and H = Gᵀ = conj(G) in the convention H_ij = (Y_j, Y_i) of §1. G has the form of a model-space reproducing kernel with phase Φ_N sampled at
the zeros, exact up to the explicit pole terms (poles of V at −2, −4, …); 'de Branges-type' is an analogy (Θ = e^{iΦ_N} is not inner: V has poles).
The single-zero Christoffel-type quantity is |c_ρ|²/G_ρρ = 1/(|ρ|²Φ_N'(γ)), Φ_N'(γ) = log N + log(γ/2π) + O(1/γ²).
Error constants depend on the fixed zero set R; nothing here is uniform over growing R.
Scalar products: c_ρ = (χ₁,Y_ρ) = (χ₁,ψ_ρ) − ∫₀^λ Vχ₁ conj(Vψ_ρ) = (ρ−1)/ρ² − Σ_j b_j(2π)^{2j}[conj V(ρ) λ^{2j+1−ρ̄}/(2j+1−ρ̄)
  + Σ_k conj(a_k(ρ))(2π)^{2k}λ^{2j+2k+1}/(2j+2k+1)];  |(ρ−1)/ρ²|² = 1/|ρ|².

## 4. Theorem (finite-λ Christoffel bound; unconditional)
For every N ≥ 2 and every finite conjugation-closed set R of simple zeros on the critical line,
  d_N² ≥ D(1/N)² ≥ B_R(N) := c_Rᴴ H_R(N)⁻¹ c_R,
with H_R, c_R the explicit expressions of §3. (Proof: §1–3; only Burnol's existence/orthogonality theorem, unitarity of V and termwise
integration of convergent power series are used.)

## 5. Expansion in 1/L (fixed R)
H = L·I + A + E(L), A = diag(Φ₀'(γ)), Φ₀' = 2θ' + 3/|ρ|², E_{ρρ'} = conj[(1 − e^{i(Φ_N(γ)−Φ_N(γ'))})/(i(γ'−γ))] (convention H_ij = (Y_j,Y_i)).
  B_R = |c|²/L − cᴴ(A + E)c/L² + O(L⁻³).
Mean part (E's non-oscillating piece i/(γ'−γ)): c₂^B = −(S₁ + P₀), S₁ = Σ_ρ Φ₀'(γ)/|ρ|², P₀ = Re Σ_{ρ≠ρ'} conj(c_ρ)c_ρ' i/(γ'−γ).
All zeros (first 3000 pairs exact + file to 2·10⁴ + density tail): S₁ = 0.100516, P₀ = +0.0030 (K = 3000; +0.00254 at K = 200), c₂^B(∞) = −0.1035.
Oscillating part: −Re Σ conj(c_ρ)c_ρ' (−i e^{−i(Φ_N(γ)−Φ_N(γ'))})/(γ'−γ): almost periodic in log N with the zero-difference frequencies γ_j − γ_k.
Finite-L values of L²B − CL include large third-order terms (+Σ Φ₀'²/|ρ|²/L³ …): −0.073 (N = 1280), −0.085 (10⁵), −0.090 (10¹⁰), −0.095 (10³⁰).

## 6. Davenport–Heilbronn and L(s, χ₋₄) (gamma-factor case only)
Scope: f with real coefficients, f(s) = X(s)f(1−s) with X(s) = (q/π)^{1/2−s}Γ((2−s)/2)/Γ((1+s)/2) (odd character type: DH with q = 5,
L(s, χ₋₄) with q = 4). Reviewer: the residue-only expansion is false in general (f = (1 + 2^{−s})² gives double poles of V_f and log t terms);
it is used here only for this X, whose 1/X has simple poles at p_k = −(2k+1) and Stirling decay to the left (the same contour argument
that reproduces Burnol's explicit series for ζ; numerically the resulting Gram is positive definite and the code reproduces ζ to 10 digits).
F = f(s)/s ∈ H² (NB generator A(1/t), A(x) = Σ_{n≤x}a_n), target χ, V_f = conj(F)/F = (s/(1−s))/X(s), V_f(0) = 0 (so V_fχ has no t⁰ term);
residues R_k = (p_k/(1−p_k))(q/π)^{p_k−1/2}·2(−1)^k/(k!Γ((2−p_k)/2)).
Vψ_w(t) = V_f(w)t^{−w} + Σ_k R_k t^{−p_k}/(p_k − w) (t < 1), V_fχ = Σ_k (R_k/p_k) t^{−p_k}; Gram/scalar products as in §3 with these terms
(s7_dh.py, function bound(); on ζ it reproduces s3 to 10 digits). Off-line zeros ρ₀ = β + iγ (β > 1/2): w = 1 − ρ̄₀ has Re w < 1/2, ψ_w ∈ K
and (ψ_w, D_θ a) ∝ conj F(ρ₀) = 0 for every θ; since V_f a = Ja and Q_λD_θJa = D_θJa for θ ≥ λ, the truncated V⁻¹Q_λVψ_w (no limit needed)
is ⊥ C_λ. Its norm² (1 − |V_f(w)|²N^{1−2β})/(2β−1) − (pole terms) interpolates between L + Re(V_f'/V_f) (for (2β−1)log N ≪ 1, the on-line
value) and 1/(2β−1) (Beurling's floor vector ψ_w) as N → ∞; the single-zero quantity |c|²/norm² is its reciprocal times 1/|ρ₀|².
On the line Re(V_f'/V_f)(ρ) = 1/|ρ|² + log(q/π) + Re ψ(3/4 + iγ/2) = log(qγ/2π) + O(1/γ²).

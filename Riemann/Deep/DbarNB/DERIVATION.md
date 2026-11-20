# DbarNB — DERIVATION (dbarnb, 20 Nov 2026). Lemmas R″ and R′ under RH; then Lemma C and (3v)

Notation (explore59): V_L = {E = ∫_[0,L] e^{−iut}dμ(u)}, ‖E‖² = ∫∫e^{−|u−v|/2}dμdμ̄ = ∫|E|²w, w = 1/(2π(¼+t²)); Z = Z(T₂),
T₂ = L⁶: ordinates γ (both signs, real by RH) with multiplicity m_γ; U_Z(ℓ) = min{‖E‖² : E ∈ V_ℓ, E = 1 on Z with mult.};
ℓ := L′ = L − √L throughout (interval/band length). A, c, C: absolute constants, changing from line to line.

## 1. Dual form (exact): U_Z(ℓ) = sup_{p∈E_Z} |p(0)|²/N_ℓ(p)
E_Z = span{u^k e^{iγu} : γ ∈ Z, k < m_γ} (functions of u ∈ [0, ℓ]); N_ℓ(p) = ∫_0^ℓ (|p′|² + ¼|p|²)du + ½|p(0)|² + ½|p(ℓ)|².
Proof. The constraints are E^{(k)}(γ) = ∫(−iu)^k e^{−iuγ}dμ = ∫ conj((iu)^k e^{iγu}) dμ, with targets 1 (k = 0), 0 (k ≥ 1); so
Σ ā_{γk}E^{(k)}(γ) = ∫ p̄_a dμ with p_a(u) = Σ a_{γk}(iu)^k e^{iγu}. For the covariance norm ‖μ‖² = ∫∫e^{−|u−v|/2}dμdμ̄ on [0,ℓ],
sup_μ |∫p̄ dμ|²/‖μ‖² = ‖p‖²_H, H = RKHS of the Ornstein–Uhlenbeck kernel e^{−θ|u−v|}, θ = ½, whose norm is
(2θ)⁻¹[∫(|p′|² + θ²|p|²) + θ(|p(0)|² + |p(ℓ)|²)] = N_ℓ(p); the representer of p is the finite measure ν_p = (¼p − p″)du +
(½p(0) − p′(0))δ₀ + (½p(ℓ) + p′(ℓ))δ_ℓ, so the minimum is attained in V_ℓ. Lagrange duality for min{‖μ‖² : jet constraints = targets} gives
U = sup_a |Σ_j a_j|²/‖p_a‖²_H, and Σ_j a_j = p_a(0) (jets have value 0 at u = 0 for k ≥ 1). ∎
Check: p = e^{iγu}: N = ℓ(γ² + ¼) + 1 = K(γ,γ) (explore59's kernel diagonal). Numerically checked (dual_check, §6).
Consequence: Lemma R″ ⇔ L|p(0)|² ≤ C·N_ℓ(p) for all p ∈ E_{Z(L⁶)}; it suffices to bound |p(0)|² by ∫_0^ℓ|p′|²χ, where
χ = ψ², ψ ≥ 0 Gevrey-2, supp ψ ⊂ [0, ℓ], ψ = 1 on [ℓ/4, 3ℓ/4], |ψ^{(N)}| ≤ C^{N+1}(N!)²ℓ^{−N} (since ∫|p′|²χ ≤ N_ℓ(p)).

## 2. Inputs
(T) Turán–Nazarov lemma (Nazarov 1993, St. Petersburg Math. J. 5 (1994) 663–717, Thm 1.3/1.4): for f(u) = Σ_{j=1}^{M} b_j e^{iλ_j u},
λ_j ∈ ℝ (jets allowed, by continuity), I an interval, E ⊂ I measurable: sup_I|f| ≤ (A|I|/|E|)^{M−1} sup_E|f|.
  (T′) E = {u ∈ I : |f|² ≤ 2 avg_I|f|²} has |E| ≥ |I|/2, so sup_I|f|² ≤ 2(2A)^{2M} |I|⁻¹∫_I|f|² ("Nikolskii").
  (T″) sup_{[0,ℓ]}|f| ≤ (2A)^{M} sup_{[ℓ/4,3ℓ/4]}|f|.  Write e^{O(M)} for products of such factors.
(Z) RH counting (Littlewood; Titchmarsh 14.13): S(t) ≪ log t/log log t, so for t ≥ t₀ every interval [t, t+h], h ≤ 1/log log t,
contains ≤ A₁ log t/log log t ordinates (with multiplicity). This and "all zeros real in z" are the only uses of RH in §1–5.
(G) Gevrey-2 Fourier decay: if φ ∈ C_c, supp φ ⊂ [0,ℓ], |φ^{(N)}| ≤ B C^N (N!)² ℓ^{−N}, then |φ̂(ζ)| ≤ C′ B ℓ e^{−c(|ζ|ℓ)^{1/2}}.

## 3. Grouping (scale s₀ = Δ/L, Δ = L^{1/4})
Groups g = connected components of Z under "|γ − γ′| < s₀". Window I_g := [min_g γ, max_g γ] (convex hull), γ_g := its midpoint,
w_g := its half-length. Then (i) the windows of g ≠ g′ are disjoint with gaps ≥ s₀; (ii) w_g ≤ M_g s₀/2 (M_g = #g with
multiplicity); (iii) M_g ≤ M*(|γ_g|), M*(t) := A₁ log t/log log t, for groups above t₀ (t₀ large, not an ordinate; L ≥ L₀ so that
no group straddles ±t₀); groups below t₀ have M_g ≤ N(t₀). Proof of (iii) (reflect negative groups): list g's ordinates in increasing
order with multiplicity, t := the smallest; if M_g > M*(t), the first ⌊M*(t)⌋+1 lie in [t, t + h], h ≤ (M*(t)+1)s₀ ≤ 7A₁ log L·L^{−3/4}
< 1/log log t (L large): contradicts (Z); and M*(t) ≤ M*(|γ_g|) (M* increasing). ∎ Hence M_g ≪ log L/log log L, e^{O(M_g)} = L^{o(1)}.

## 4. Local lemmas (one group; f = e^{iγ_g u} r(u), r with M terms and frequencies δ_j = γ_j − γ_g ∈ [−w, w])
Lemma A (complex extension). Let W := max(w, K M²/ℓ), K absolute. For x ∈ [0,ℓ], ζ ∈ ℂ, |ζ| ≤ 1/W:
|r(x + ζ)| ≤ C_A^{M} sup_{[0,ℓ]}|r|; hence (Cauchy) |r^{(n)}(x)| ≤ C_A^{M} n! W^n sup_{[0,ℓ]}|r| for x ∈ [0,ℓ], n ≥ 0.
Proof. Fix χ₀ ∈ G²_c, χ₀ = 1 on [−3/2, 3/2], supp ⊂ [−2,2]. m(δ) := e^{iδζ}χ₀(δ/W) equals e^{iδζ} near [−w,w] and φ := m^∨,
φ(s) = WΦ_a(Ws), a = Wζ, |a| ≤ 1, |Φ_a(τ)| ≤ C e^{−c|τ|^{1/2}} uniformly (G). Since ∫φ(s)(x+s)^k e^{iδ(x+s)}ds =
(−i∂_δ)^k[e^{iδx}m(δ)] = (x+ζ)^k e^{iδ(x+ζ)} for δ near [−w,w]: r(x+ζ) = ∫φ(s) r(x+s) ds. Where x+s ∈ [0,ℓ]: ≤ ‖Φ‖₁ sup|r|.
Elsewhere (T) with I ⊃ [0,ℓ] ∪ {x+s}, |I| ≤ ℓ + |s|: |r(x+s)| ≤ (A(1 + |s|/ℓ))^{M−1} sup|r|; the tail is
≤ A^{M}∫|Φ(τ)|(1 + |τ|/(KM²))^{M}dτ ≤ A^{M} C KM² ∫_0^∞ e^{−M(c√K σ^{1/2} − log(1+σ))}dσ ≤ 4CK A^{M} once c√K ≥ 2
(log(1+σ) ≤ σ^{1/2}; substitution τ = KM²σ). Cauchy's estimate on |ζ| = 1/W gives the derivative bound. ∎

Lemma B (endpoint value with the γ² gain). Let |γ_g| ≥ 14, f = p_g ∈ E_g, q_g = p_g′, I = [ℓ/4, 3ℓ/4]. If
8(2A)^{4M}C_A^{2M}(W/γ_g)² ≤ ½, then |p_g(0)|² ≤ B(M)·(ℓγ_g²)^{−1}∫|q_g|²χ, B(M) := 32(2A)^{4M}.
Proof. p_g = e^{iγ_g u}r, q_g = e^{iγ_g u}(iγ_g r + r′), so |γ_g r| ≤ |q_g| + |r′|. S := sup_I|r|² ≤ (8/ℓ)(2A)^{2M}∫_I|r|² (T′),
∫_I|r|² ≤ 2γ_g^{−2}(∫_I|q_g|² + ∫_I|r′|²), ∫_I|r′|² ≤ (ℓ/2)C_A^{2M}W² sup_{[0,ℓ]}|r|² ≤ (ℓ/2)C_A^{2M}W²(2A)^{2M}S (Lemma A, T″).
So S ≤ 16(ℓγ_g²)^{−1}(2A)^{2M}∫_I|q_g|² + 8(2A)^{4M}C_A^{2M}(W/γ_g)²S; absorb (factor 2); |p_g(0)|² = |r(0)|² ≤ (2A)^{2M}S; χ = 1 on I. ∎
The hypothesis holds for every group once L ≥ L₀: W_g ≤ 2M_g L^{−3/4} and (2A)^{2M}C_A^M = L^{o(1)} by §3(iii).
By-product (used in Lemma E): sup_{[0,ℓ]}|r|² ≤ B(M)(ℓγ_g²)^{−1}‖q_g‖²_χ and sup_{[0,ℓ]}|r_q|² ≤ e^{O(M)}ℓ^{−1}‖q_g‖²_χ, where
q_g = e^{iγ_g u} r_q (T′, T″); ‖·‖_χ := (∫|·|²χ)^{1/2}.

Lemma Λ (leakage of one group in frequency). F_g(ξ) := ∫ q_g ψ e^{−iξu}du, d(ξ) := dist(ξ, [γ_g − w, γ_g + w]). Then
|F_g(ξ)| ≤ e·C_A^{M}·C M (CWℓ + 1)^{M}·ℓ sup_{[0,ℓ]}|r_q|·e^{−c(d(ξ)ℓ)^{1/2}}.
Proof. Newton expansion of δ ↦ e^{iδv} at the nodes δ_1..δ_M (repeated for jets), centre x = ℓ/2:
r_q(u) = Σ_{k≤M} c_k ν_k(u − x), ν_k(v) = (e^{i·v})[δ_1..δ_k] = ∫_{S_{k−1}}(iv)^{k−1}e^{iθ_t v}dt (Hermite–Genocchi, θ_t = Σt_lδ_l ∈ [−w,w],
|S_{k−1}| = 1/(k−1)!), c_k = [Π_{l<k}(−i∂ − δ_l) r_q](x) (exact on exponentials and jets, since the Newton–Hermite interpolant
is exact at the nodes). |c_k| ≤ Σ_n binom(k−1,n) w^{k−1−n}|r_q^{(n)}(x)| ≤ e C_A^M (k−1)! W^{k−1} sup|r_q| (Lemma A).
F_g(ξ) = Σ_k c_k ∫_{S_{k−1}} e^{−iθ_t x} φ̂_k(ξ − γ_g − θ_t) dt, φ_k(u) := (i(u − x))^{k−1}ψ(u); φ_k is Gevrey-2 with
B = (ℓ/2)^{k−1}e^{2k}C (Leibniz; binom(N,j)((N−j)!)² ≤ (N!)²/j!), so |φ̂_k(ζ)| ≤ C^k ℓ^k e^{−c(|ζ|ℓ)^{1/2}} (G), and
|ξ − γ_g − θ_t| ≥ d(ξ). Summing: Σ_k e C_A^M (k−1)!W^{k−1} sup|r_q| ·C^kℓ^k/(k−1)! · e^{−c(dℓ)^{1/2}}. ∎
With the by-product: |F_g(ξ)|² ≤ a_g² e^{−2c(d(ξ)ℓ)^{1/2}}, a_g² := e^{O(M_g log(W_gℓ + 2))}ℓ‖q_g‖²_χ = e^{o(log²L)}ℓ‖q_g‖²_χ.

## 5. Assembly
Lemma D (almost orthogonality of groups). For p ∈ E_Z, q = p′ = Σ_g q_g: Σ_g ‖q_g‖²_χ ≤ (1 + ε_L)‖q‖²_χ, ε_L ≤ e^{−cL^{1/8}}.
Proof. Partition ℝ_ξ into cells Ω_g (points nearer to g's window than to any other); windows have gaps ≥ s₀, so d_g ≥ s₀/2 off
Ω_g and #{g′ : d_{g′}(ξ) ∈ [js₀/2, (j+1)s₀/2)} ≤ 2 (j ≥ 1). F := (qψ)^ = Σ_g F_g, ∫|F|² = 2π‖q‖²_χ, ∫|F_g|² = 2π‖q_g‖²_χ.
(a) ∫_{Ω_g^c}|F_g|² ≤ a_g²∫_{d≥s₀/2}e^{−2c(dℓ)^{1/2}}dξ ≤ τ‖q_g‖²_χ, τ := e^{o(log²L) − c(s₀ℓ/2)^{1/2}} ≤ e^{−c′L^{1/8}} (s₀ℓ ≍ L^{1/4}).
(b) On Ω_g, by Cauchy–Schwarz, |Σ_{g′≠g}F_{g′}|² ≤ (Σ_{g′≠g}e^{−c(d_{g′}ℓ)^{1/2}})(Σ_{g′≠g}a_{g′}²e^{−c(d_{g′}ℓ)^{1/2}}), first factor
≤ 4Σ_{j≥1}e^{−c(js₀ℓ/2)^{1/2}} ≤ τ; integrating over Ω_g and summing over g: Σ_g∫_{Ω_g}|Σ_{g′≠g}F_{g′}|² ≤ τ²Σ_{g′}‖q_{g′}‖²_χ.
(c) |a+b|² ≥ (1−η)|a|² − η⁻¹|b|² on each Ω_g with η = τ^{1/2}: 2π‖q‖²_χ ≥ (1−η)(2π − τ)Σ‖q_g‖²_χ − η⁻¹τ²Σ‖q_g‖²_χ. ∎

Theorem R″ (RH). For L ≥ L₀: L·U_{Z(L⁶)}(L − √L) ≤ (1 + o(1))S, S := sup_L Σ_g B(M_g)/γ_g² < ∞.
Proof. |p(0)| ≤ Σ_g|p_g(0)| ≤ Σ_g (B(M_g)/(ℓγ_g²))^{1/2}‖q_g‖_χ ≤ (S/ℓ)^{1/2}(Σ_g‖q_g‖²_χ)^{1/2} ≤ ((1+ε_L)S/ℓ)^{1/2}‖p′‖_χ, and
‖p′‖²_χ ≤ ∫_0^ℓ|p′|² ≤ N_ℓ(p); then §1 and L/ℓ → 1. Finiteness of S: charge each group to its largest-|γ| ordinate ρ (distinct
representatives, |γ_g| ≤ |γ_ρ| ≤ √2|γ_g|);
by §3(iii) S ≤ C(t₀) + 2Σ_ρ exp(β A₁ log|γ_ρ|/log log|γ_ρ|)/|γ_ρ|² < ∞ (B(M) = e^{βM}; Σ_ρ|γ|^{−2+o(1)} converges). ∎
So, with explore59's reduction (3v): RH (+ Lemma C) ⇒ D(λ)² ≪ 1/log(1/λ).


Lemma E (low/high decoupling; H fixed, H ≠ |γ| for all γ). Z_lo = Z(H), Z_hi = Z(L⁶)∖Z(H); for L ≥ L₀(H) the groups do not
straddle ±H (s₀ < g_H := dist(Z_lo, Z_hi) > 0). For p_lo ∈ E_{Z_lo}, p_hi ∈ E_{Z_hi}:
|⟨p_lo, p_hi⟩_N| ≤ θ_L (N(p_lo)N(p_hi))^{1/2}, θ_L = C(H)L^{−3/4+o(1)}, where ⟨p,q⟩_N = ∫_0^ℓ(p′q̄′ + ¼pq̄) + ½p(0)q̄(0) + ½p(ℓ)q̄(ℓ).
Proof. Boundary: |p_lo(0)|² ≤ U_{Z_lo}(ℓ)N(p_lo) ≤ C(H)ℓ⁻¹N(p_lo) (§1; Claim F) and |p_hi(0)|² ≤ Cℓ⁻¹N(p_hi) (proof of R″ applied
to the high groups); the same at u = ℓ, since p ↦ conj(p(ℓ − ·)) maps E_Z onto E_Z (same frequencies; jets stay jets) and preserves N.
Interior: p_lo = Σ a_{ρk}e_{ρk}, e_{ρk} = u^k e^{iγ_ρ u}, |γ_ρ| < H, k < m_ρ: n_H fixed functions; as ℓ → ∞ the Gram matrix of
ℓ^{−k−1/2}e′_{ρk} on [0,ℓ] tends to blockdiag(γ_ρ²·Hilb_{m_ρ}) (cross blocks O(1/(ℓ·min gap of Z(H)))), so ∫_0^ℓ|p_lo′|² ≥ c(H)Σ|a_{ρk}|²ℓ^{2k+1}.
For a high group g, q_g = e^{iγ_g u}r_g, Ω := γ_g − γ_ρ (|Ω| ≥ g_H, and |Ω| ≥ |γ_g|/2 for |γ_g| ≥ 2H), v := (ku^{k−1} + iγ_ρu^k):
one integration by parts and Lemma A give |∫_0^ℓ v e^{−iΩu} r̄_g| ≤ C(H)ℓ^k|Ω|⁻¹(3 + ℓC_A^{M}W_g) sup|r_g| ≤ C(H)ℓ^{k−1/2}L^{1/4+o(1)}|Ω|⁻¹‖q_g‖_χ
(Lemma B by-product; ℓW_g ≤ 2M_g L^{1/4}). Summing over g by Cauchy–Schwarz with Σ_g e^{O(M_g)}/Ω_g² ≤ C(H) and Lemma D:
|∫_0^ℓ p_lo′ p̄_hi′| ≤ C(H)L^{1/4+o(1)}Σ_{ρ,k}|a_{ρk}|ℓ^{k−1/2}‖p_hi′‖_χ ≤ C(H)L^{−3/4+o(1)}‖p_lo′‖‖p_hi′‖_χ. The ¼∫p_lo p̄_hi term is the
same with v = u^k and r_{p,g} (sup|r_{p,g}|² ≤ B(M)(ℓγ_g²)⁻¹‖q_g‖²_χ). Finally ‖p_lo′‖², ¼‖p_lo‖², ‖p_hi′‖²_χ ≤ N. ∎

Theorem R′ (RH). For each fixed H (not an ordinate), with ε(H) := 2Σ_{|γ_ρ|>H−1} B(M*(γ_ρ))/|γ_ρ|² (+ the finite part below t₀
when H < t₀) → 0 as H → ∞: for L ≥ L₀(H),
L·U_{Z(L⁶)}(L − √L) ≤ (1 − θ_L)⁻¹ [L·U_{Z(H)}(L − √L) + (1 + ε_L)(L/ℓ) ε(H)].
Proof. p = p_lo + p_hi. |p_lo(0)|² ≤ U_{Z(H)}N(p_lo); |p_hi(0)|² ≤ (1+ε_L)(S_hi/ℓ)‖p_hi′‖²_χ ≤ (1+ε_L)(S_hi/ℓ)N(p_hi), S_hi :=
Σ_{g high}B(M_g)/γ_g² ≤ ε(H). Cauchy–Schwarz: |p(0)|² ≤ (U_{Z(H)} + (1+ε_L)S_hi/ℓ)(N(p_lo) + N(p_hi)), and by Lemma E
N(p) ≥ (1 − θ_L)(N(p_lo) + N(p_hi)). Then §1. ∎
Corollary (RH). lim_{L→∞} L·U_{Z(L⁶)}(L − √L) = Σ_ρ m_ρ²/|ρ|². Proof: ≥ L·U_{Z(H)}(ℓ) for every H (fewer constraints), and
Claim F (explore59, REPAIR (2)) gives L·U_{Z(H)}(ℓ) → Σ_{|γ|<H}m²/|ρ|²; R′ gives limsup ≤ Σ_{|γ|<H}m²/|ρ|² + ε(H); let H → ∞. ∎
Remark. R′ is proved in the form needed (H fixed, L → ∞); the brief's uniformity in H ≤ L⁶ is not needed and not claimed.
Remark (constants). Only (T) with its exponential constant is used; Erdélyi (arXiv:1602.02315, eq. (1.4)–(1.5)) gives the
polynomial endpoint bound |f(0)| ≤ C n‖f‖_{L²[0,1]} for n real frequencies, consistent with the orchestrator's numerics
(coalescent jet worst, C(x)/m² ≤ 1, R508); under RH the exponential constants cost only |γ|^{o(1)} (§3(iii)).


## 6. Lemma C and the reduction (3v), written out (closes the referee's (3iii) gap; Cartan curves no longer used)
Lemma C (Titchmarsh, Theory of the Riemann zeta-function, Thm 9.7; unconditional). There is A such that every interval
[T, T+1], T ≥ 2, contains t with |ζ(σ+it)|⁻¹ ≤ t^A for −1 ≤ σ ≤ 2. Also (convexity, Cauchy): |ζ(w)|, |ζ′(w)| ≤ C(2+|Im w|)²
for −1 ≤ Re w ≤ 2, |Im w| ≥ 1; |1/ζ(w)| ≤ C for Re w = 3/2 and, by the functional equation, for Re w = −½, |Im w| ≥ 1.
Tiling. Pick good heights t_n ∈ [3n, 3n+1] and t′_n ∈ [3n + 3/2, 3n + 5/2] (Lemma C); n₀ := ⌈L⁶/3⌉, T₂ := t_{n₀} ∈ [L⁶, L⁶+4].
Low set Z(T₂): all zeros with |γ| ≤ T₂ as complex nodes z_ρ = γ − i(β − ½), Hermite conditions (real under RH; §1–5 verbatim).
Families: Q_n := [−½, 3/2] × [t_n, t_{n+1}] (n ≥ n₀); Q′_{n₀} := [−½, 3/2] × [T₂, t′_{n₀}], Q′_n := [−½, 3/2] × [t′_{n−1}, t′_n] (n > n₀);
wide versions Q^w, Q′^w with σ ∈ [−1, 2]; mirror images of all. Every zero with |γ| > T₂ (0 < β < 1 unconditionally) lies in exactly
one rectangle of each family (ζ ≠ 0 on the horizontal sides by Lemma C): no omission, no double count.
Corrections. H := ΨΦ (entire), F₀ := H/ζ, B(x) = ((1 − e^{−ax})/(ax))², a = (L−1)/2 (|B(x)| ≤ min(1, 4/(a|x|)²)e^{(L−1)max(0,−Re x)}).
R_n(s) := Σ_{ρ∈Q_n} Res_{w=ρ}[F₀(w)B(s−w)/(s−w)] (all multiplicities), R_high := Σ_n (R_n + R̄_n).
Family independence: partial sums of two height families up to t_N and t′_N differ by the residues of the zeros between these
good heights, = (2πi)⁻¹∮ over that rectangle of F₀(w)B(s−w)/(s−w)dw = O_s(t_N^{A+2}√L e^{3L}e^{−c(t_N/L²)^{1/2}}) → 0; narrow and wide
rectangles of equal heights hold the same zeros. So R_high is the same limit for all four families (no absolute convergence over
single zeros, i.e. no bound on 1/ζ′(ρ), is needed).
(i) Residue theorem (the only other pole of the integrand is w = s, residue −F₀(s)B(0) = −F₀(s)):
R_n(s) = I_n(s) := (2πi)⁻¹∮_{∂Q_n}F₀(w)B(s−w)/(s−w)dw for s ∉ Q_n, R_n(s) = I_n(s) + F₀(s) for s ∈ Q_n° off the zeros.
(ii) Divided-difference form (referee's repair): with ζ[s,w] := (ζ(s) − ζ(w))/(s − w) (holomorphic for s, w ≠ 1, removable at s = w), for every s ≠ 1
ζ(s)R_n(s) = (2πi)⁻¹∮_{∂Q_n} H(w)B(s−w)ζ[s,w]/ζ(w) dw,
because (2πi)⁻¹∮H(w)B(s−w)/(s−w)dw = −H(s)·1_{s∈Q_n} (H, B entire) and ζF₀ = H. No singularity where Re s = ½ crosses ∂Q_n.
(iii) Bound on the line s = ½ + it: on ∂Q_n, |1/ζ| ≤ t_{n+1}^A, |B(s−w)| ≤ e^{L−1}min(1, 16/((L−1)dist(s,∂Q_n))²),
|ζ[s,w]| ≤ C(2 + |t| + t_{n+1})² (|s − w| ≤ 1: sup|ζ′| on the segment; else |ζ(s)| + |ζ(w)|), |Φ(w)| ≤ Ce^{−c(t_n/T₁)^{1/2}} ≤ Ce^{−cL²}
(t_n ≥ L⁶, T₁ = L²), |Ψ(w)| ≤ C(1+|w|)√L e^{L′} (REPAIR (3i)). len ∂Q_n ≤ 12. Hence
‖ζR_high‖²_w ≤ C Σ_{n≥n₀} t_{n+1}^{2A+8} L e^{4L} e^{−2c(t_n/L²)^{1/2}} ≤ e^{−c′L²}, and |R_high(1)| = |Σ_n (I_n(1) + Ī_n(1))| ≤ e^{−c′L²} directly (1 lies outside all rectangles).
(iv) Entire + growth. Σ_n R_n converges locally uniformly off the high zeros (for s in a compact set, R_n = I_n for n large,
|I_n(s)| ≤ C t_{n+1}^{A+2}√L e^{3L}e^{−c(t_n/L²)^{1/2}}); F₀ − R_high has no poles at high zeros (principal parts cancel), and
R_tr removes the trivial ones (REPAIR (3ii)), so G := F₀ − R_high − R_tr is entire. Growth: R_high is a sum over all high
poles, so any of four tilings may represent it: widths [−½, 3/2] (narrow) or [−1, 2] (wide, Lemma C covers −1 ≤ σ ≤ 2), heights
t_n or t′_n (horizontal sides ≥ ½ apart, except the shared bottom T₂, excluded by |t| > T₂ + 4). For σ ∈ [−¾, 7/4], |t| > T₂ + 4: wide family whose horizontal
sides are ≥ ¼ from t; G(s) = −Σ_m I^w_m(s) − R_tr(s), |I^w_m(s)| ≤ C t_{m+1}^{A+1}√L e^{2L}e^{−c(t_m/L²)^{1/2}}e^{(L−1)(2+max(0,−σ))}
·min(1, C(L|t − t_m|)⁻²) (dist(s, ∂Q^w_m) ≥ ¼), summable with polynomial growth in |t|. For |t| ≤ T₂ + 4, σ ∈ [−¾, 7/4]: compact, G
bounded by C_L. For σ ∉ [−¾, 7/4]: narrow family (s outside all Q_m, distance ≥ ¼), G = F₀ − Σ_m I_m − R_tr with σ > 7/4:
|F₀| ≤ C(1+|s|)√L; σ < −¾: |F₀| ≤ C_L(1+|s|)^C e^{(L′+2/T₁)|σ|} off the discs |s+2k| < ½ (|1/χ(s)| ≤ C, |1/ζ(1−s)| ≤ C); inside the discs: bound F₀,
R_high, R_tr separately on |s+2k| = ½, then the maximum principle for the entire G, |I_m(s)| ≤ C_L e^{(L−1)(|σ|+3/2)}. So |G(s)| ≤ C_L(1+|s|)^C e^{(L−½)max(0,−σ)} (L′ + 2/T₁ ≤ L − ½), PWS applies and
DERIVATION §7 of explore59 (closure) puts ζG/s in the closed span. No use of RH in (i)–(iv).
(v) On the line 1 − ζG = (1 − Φ) + ΦE_low + ζR_tr + ζR_high (ζF₀ = ΨΦ, Ψ = 1 − E_low), ‖ΦE_low‖ ≤ ‖E_low‖ = U^{1/2}
(|Φ| ≤ 1), ‖1 − Φ‖² ≪ L⁻², ‖ζR_tr‖² ≪ e^{−5(√L−2)} (REPAIR (3ii)), ‖ζR_high‖² ≤ e^{−c′L²} by (iii). With ‖a+b‖² ≤ (1+η)‖b‖² +
(1+η⁻¹)‖a‖², η = L^{−1/2}: D(λ)² ≤ (1 + L^{−1/2}) U_{Z(T₂)}(L − √L) + O(L^{−3/2}) — unconditionally (complex nodes).
Where RH is used in the whole chain: only §1–5 (frequencies real; cluster bound §3(iii)). Off-line zeros below T₂ would make
U ≥ floor ≈ 2b/|ρ|² (explore59 §4), so the reduction is consistent.

## 7. Theorem
RH ⇒ lim_{λ→0} D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|²  (= 2 + γ − log 4π if all zeros are simple).
Proof. Upper: §6(v) + Corollary of R′ (§5): limsup L·D² ≤ lim L·U_{Z(T₂)}(L − √L) = Σ m²/|ρ|². Lower: Burnol (math/0103058),
liminf L·D² ≥ Σ m²/|ρ|² unconditionally. ∎   (R″ alone gives D(λ)² ≪ 1/log(1/λ).)


# NBDir — DERIVATION (nbdir, 20 Nov 2026). GRH ⇒ lim L·D_χ(λ)² = Σ_ρ m_ρ²/|ρ|² for primitive χ

Notation. χ primitive mod q ≥ 3 (hence nonprincipal, conductor q; induced characters have extra zeros on Re s = 0), parity a ∈ {0,1} (χ(−1) = (−1)^a); §§1–5 hold for every primitive χ; §§6–7 are written for real χ, §8 extends them. K = L²(]0,∞[, dt), (f,g) = ∫fḡ;
Mellin f̂(s) = ∫_0^∞ f t^{s−1}dt, Re s = ½. a_χ(t) := A_χ(1/t) (= 0 for t > 1, |a_χ| ≤ q), â_χ(s) = L(s,χ)/s (Re s > 0); 1 := 1_{]0,1]},
1̂ = 1/s; D_θf(t) = θ^{−½}f(t/θ) (multiplier θ^{s−½}); A_χ(θ/t) = θ^{½}D_θa_χ(t); B_λ := span{D_θa_χ : λ ≤ θ ≤ 1}; D_χ(λ) := dist(1, B_λ);
L := log(1/λ). Q_λ := multiplication by 1_{[λ,∞)}. ψ_{w,k}(t) := (log 1/t)^k t^{−w}1(t) ∈ K for Re w < ½.

## 1. Phase operator (Burnol math/0103058 l.296–364, with ζ replaced by L(s,χ))
γ_a(s) := (q/π)^{s−½}Γ((s+a)/2)/Γ((1−s+a)/2). Functional equation (Davenport, Multiplicative NT, ch. 9): ξ(s,χ) := (q/π)^{(s+a)/2}Γ((s+a)/2)L(s,χ)
= W(χ)ξ(1−s,χ̄), W(χ) = τ(χ)/(i^a√q), |W| = 1, W = 1 for real χ. So L(1−s,χ̄) = W(χ)^{−1}γ_a(s)L(s,χ). On Re s = ½, conj L(s,χ) = L(1−s,χ̄),
so the phase conj(â_χ)/â_χ equals V_χ(s) := W(χ)^{−1}(s/(1−s))γ_a(s) off the zeros; |V_χ| = 1 on the line ((1−s+a)/2 = conj((s+a)/2), |s| = |1−s|).
Operators: M f(t) = t^{−1}∫_0^t f (multiplier 1/(1−s); M − 1 ↦ s/(1−s), unitary); I f(t) = t^{−1}f(1/t) ((If)^(s) = f̂(1−s)); F₊f = 2∫_0^∞f(u)cos(2πtu)du,
F₋f = 2∫_0^∞f(u)sin(2πtu)du (unitary on L²(0,∞)). A Hankel operator ∫k(tu)f(u)du has Mellin k̂(s)f̂(1−s), so F_±I ↦ γ_±(s) := Mellin of 2cos/2sin(2πx):
γ₊ = 2(2π)^{−s}cos(πs/2)Γ(s) = π^{½−s}Γ(s/2)/Γ((1−s)/2) (Burnol l.315), γ₋ = 2(2π)^{−s}sin(πs/2)Γ(s) = π^{½−s}Γ((1+s)/2)/Γ(1−s/2) (duplication +
reflection; checked numerically to 1e-21, part A of logs/burnol_m4_AB.log and logs/burnol_5_ABCDE.log). Since γ_a = q^{s−½}γ_± (+ for a = 0, − for a = 1):
  V_χ = W^{−1}(M − 1)D_qF_±I  (a dilated cosine (a = 0) or sine (a = 1) transform); unitary, commutes with all D_θ.
(P1) V_χa_χ = Ja_χ, J g := conj(Ig) (Mellin of Jg on the line = conj ĝ); Ja_χ(t) = t^{−1}A_χ̄(t) vanishes for t < 1, so for θ ∈ [λ,1]
D_θJa_χ is supported in [θ,∞) ⊂ [λ,∞): Q_λD_θJa_χ = D_θJa_χ. (For ζ Burnol needs Z = ((s−1)/s)ζ(s)/s and V = (s/(1−s))³γ₊; here â_χ ∈ H² already.)

## 2. Pointwise kernel lemma (the riskiest step; a = 1 is new, a = 0 is Burnol's Thms theo2/pretheo, l.477–531)
g₀(x) := sin(2πx)/(πx), g₁(x) := (1 − cos 2πx)/(πx); U_a(w) := (w/(1−w))γ_±(w). Kernel: (M−1)[2cos(2πt/v)/v] = sin(2πt/v)/(πt) − 2cos(2πt/v)/v
= ∂_v g₀(t/v), (M−1)[2sin(2πt/v)/v] = (1−cos(2πt/v))/(πt) − 2sin(2πt/v)/v = ∂_v g₁(t/v) (g_a(t/v) = (v/πt)(…)).
Lemma Bχ. For Re w < 1, k ≥ 0, t > 0, φ_{w,k}(t) := lim_{δ→0}∫_δ^1(log 1/v)^k v^{−w}∂_v g_a(t/v)dv exists, is holomorphic in w, φ_{w,k} = ∂_w^kφ_{w,0}, and
 (i) |φ_{w,k}(t)| ≤ C/t on [1,∞), uniformly for w in compacts of {Re w < 1};
 (ii) φ_{w,k}(t) = ∂_w^k(U_a(w)t^{−w}) + O(1) on ]0,1], uniformly for 0 < Re w < 1, |Im w| ≤ B (O(t) for a = 1, k ≥ 0).
Proof. Integrate by parts, u = 1/v (boundary terms at v = δ are O(δ^{1−Re w}(log 1/δ)^k) → 0; at v = 1 they vanish for k ≥ 1):
 φ_{w,0}(t) = g_a(t) + (w/(πt))∫_1^∞u^{w−2}c_a(2πtu)du, φ_{w,k}(t) = (πt)^{−1}∫_1^∞(k + w log u)(log u)^{k−1}u^{w−2}c_a(2πtu)du (k ≥ 1),
 c₀ = sin, c₁ = 1 − cos. Absolutely convergent for Re w < 1 (|c_a| ≤ 2): existence, holomorphy, ∂_w^k, and (i) (|g_a(t)| ≤ 2/(πt)).
 For 0 < Re w < 1 split ∫_1^∞ = ∫_0^∞ − ∫_0^1. With μ = w − 1: ∫_0^∞x^{μ−1}sin x dx = Γ(μ)sin(πμ/2) (−1 < Re μ < 1) and ∫_0^∞x^{μ−1}(1 − cos x)dx
 = −Γ(μ)cos(πμ/2) (−2 < Re μ < 0) give (w/(πt))∫_0^∞u^{w−2}c_a(2πtu)du = (w/(1−w))·2(2π)^{−w}Γ(w)·{cos, sin}(πw/2)·t^{−w} = U_a(w)t^{−w}.
 Remainders: (w/(πt))∫_0^1u^{w−2}sin(2πtu)du = 2 + w·Σ_{j≥1}(−1)^j 2(2πt)^{2j}/((2j+1)!(w+2j)) and (w/(πt))∫_0^1u^{w−2}(1 − cos 2πtu)du
 = w·Σ_{j≥1}(−1)^{j+1}(2π)^{2j}t^{2j−1}/(π(2j)!(w+2j−1)): |w + 2j − 1| ≥ 1, so they and their w-derivatives are O(1) (a = 0), O(t) (a = 1)
 uniformly; g₀ = O(1), g₁ = O(t) on ]0,1]. ∎   (No hypothesis on L(s,χ) enters; the lemma is about the cosine/sine kernels only.)

## 3. The vectors Y^λ_{s,k} (Burnol Thm-Def l.400–466 and Thm theo l.555–570, for χ)
Lemma Yχ. (a) For Re w < ½: V_χψ_{w,k} = Φ_{w,k} a.e., Φ_{w,k}(t) := W^{−1}q^{−½}φ_{w,k}(t/q).
(b) For Re s = ½ the L²-limit Y^λ_{s,k} := lim_{w→s, Re w<½}V_χ^{−1}Q_λV_χψ_{w,k} exists, V_χY^λ_{s,k} = 1_{[λ,∞)}Φ_{s,k}, and
 Φ_{s,k}(t) = ∂_s^k(V_χ(s)t^{−s}) + O(1) on ]0,1], = O(1) on [1,q], = O(1/t) on [q,∞); constants independent of λ, uniform for |Im s| ≤ B.
Proof. (a) ψ^δ := ψ_{w,k}1_{[δ,1]} → ψ_{w,k} in K; Iψ^δ is bounded with compact support, so F_±Iψ^δ(t) = 2∫_δ^1ψ(v){cos,sin}(2πt/v)dv/v is an
absolutely convergent integral, and M − 1 passes under it (Fubini): V_χψ^δ = W^{−1}D_q[∫_δ^1ψ(v)∂_vg_a(·/v)dv] → Φ_{w,k} pointwise (Lemma Bχ);
L²- and pointwise limits agree a.e. (b) Q_λΦ_{w,k} → Q_λΦ_{s,k} pointwise (holomorphy); domination on [λ,q] by C_λ (Lemma Bχ (ii): the main
term is ≤ C(1 + log(q/t))^k(t/q)^{−Re w} with U_a holomorphic on 0 < Re w < 1) and on [q,∞) by C/t (Lemma Bχ (i)); dominated convergence,
V_χ^{−1} unitary. Estimates: Lemma Bχ at t/q, with W^{−1}q^{−½}q^{s}U_a(s) = W^{−1}(s/(1−s))q^{s−½}γ_±(s) = V_χ(s). ∎
(Simpler than Burnol: V_χ contains (M − 1) once, so no M^L-iterates; even χ: V_χ = D_q U with U Báez-Duarte's operator itself.)

## 4. Orthogonality, Gram, target (Burnol l.418–466, theogram l.600–625, theoscal l.627–650)
(O) θ ∈ [λ,1], Re w < ½: (V_χ^{−1}Q_λV_χψ_{w,k}, D_θa_χ) = (Q_λV_χψ, D_θJa_χ) = (V_χψ, D_θJa_χ) (P1) = (ψ_{w,k}, D_θa_χ)
 = ∂_w^k θ^{½−w}∫_0^1 t^{−w}conj(a_χ)dt = ∂_w^k[θ^{½−w}L(1−w,χ̄)/(1−w)] (holomorphic for Re w < 1; a pole at w = 1 if L(0,χ̄) ≠ 0). Let w → s and conjugate (conj(L(1−s,χ̄)/(1−s)) = L(s,χ)/s):
 (D_θa_χ, Y^λ_{s,k}) = (−∂_s)^k[θ^{s−½}L(s,χ)/s]. Vanishing for all θ in [λ,1] ⟺ (L(s,χ)/s)^{(j)} = 0 for j ≤ k. So Y^λ_{ρ,k} ⊥ B_λ for each zero ρ
 of L(s,χ) on the line and k < m_ρ.
(G) X^λ_{ρ,k} := L^{−½−k}Y^λ_{ρ,k}: (X^λ_{ρ,k}, X^λ_{ρ′,l}) → δ_{ρρ′}/(k+l+1). Indeed (Y,Y′) = (V_χY, V_χY′) = ∫_λ^∞Φ_{ρ,k}conj(Φ_{ρ′,l}); with
 ∂_s^k(V(s)t^{−s}) = V(s)(log 1/t)^k t^{−s} + (lower log powers)·t^{−s} and |V_χ(ρ)| = 1: ∫_λ^1(log 1/t)^{k+l}t^{−1}dt = L^{k+l+1}/(k+l+1) for ρ = ρ′,
 O(L^{k+l}) for γ ≠ γ′ (one integration by parts of t^{−1−i(γ−γ′)}); main × O(1) terms: ∫_0^1 t^{−½}(log 1/t)^k dt < ∞; [1,∞): bounded.
(T) √L(1, X^λ_{ρ,k}) → δ_{k0}/ρ. Mellin of V_χ1 is V_χ(s)/s = W^{−1}γ_a(s)/(1−s) = W^{−1}γ_a(s)(I1)^(s), so V_χ1 = W^{−1}D_qF_±(1) = W^{−1}q^{−½}g_a(t/q):
 bounded on ]0,q], O(1/t) beyond. (1, Y^λ_{ρ,k}) = (V_χ1, Q_λΦ_{ρ,k}) has a finite limit as λ → 0 (Φ_{ρ,k} = O(t^{−½}(1+log 1/t)^k) on ]0,1]), so
 L^{−k}(1,Y) → 0 for k ≥ 1. k = 0: for Re w < ½, ∫_0^∞V_χ1·conj Φ_{w,0} = (1, ψ_{w,0}) = 1/(1−w̄), and |∫_0^λ V_χ1·conj Φ_{w,0}| ≤ C∫_0^λ(t^{−Re w} + 1)dt
 =: ε(λ) → 0 uniformly near ρ (Lemma Bχ (ii)); so |∫_λ^∞V_χ1 conj Φ_{ρ,0} − 1/ρ| ≤ ε(λ) (w → ρ, 1/(1−ρ̄) = 1/ρ), i.e. (1, Y^λ_{ρ,0}) → 1/ρ.

## 5. Theorem Bχ (Burnol's liminf for χ; unconditional, every primitive χ)
liminf_{λ→0} L·D_χ(λ)² ≥ Σ_{ρ: Re ρ = ½} m_ρ²/|ρ|² (distinct zeros of L(s,χ) on the line). If L(s,χ) has a nontrivial zero off the line, D_χ(λ) ≥ c_χ > 0.
Proof. R a finite set of zeros on the line: by (O), D_χ(λ)² ≥ ‖P_R1‖², P_R = projection onto span{X^λ_{ρ,k}: ρ ∈ R, k < m_ρ}; by (G), (T) and
[Hilb_m^{−1}]₀₀ = m² (Cauchy; Hilb_{jk} = 1/(j+k+1)), L‖P_R1‖² → Σ_{ρ∈R}m_ρ²/|ρ|²; take sup over R. Off the line: the nontrivial zeros of L(s,χ) are symmetric
under ρ ↦ 1 − ρ̄ (functional equation and conj; not the trivial ones), so some nontrivial zero has β₀ > ½ (and β₀ < 1). Every 1 − f, f ∈ B_λ, lies in L²(]0,1]), its Mellin is (1 − L(s,χ)G(s))/s,
and |∫_0^1 g t^{ρ₀−1}dt| ≤ ‖g‖(2β₀ − 1)^{−½}: ‖1 − f‖ ≥ (2β₀−1)^{½}/|ρ₀|. ∎

## 6. Upper bound under GRH: DbarNB DERIVATION §1–7 + ExploreDbar REPAIR (3i), (3ii), §7, Claim F with ζ → L(s,χ)
Mellin form: D_χ(λ)² = inf_G ∫|1 − L(½+it,χ)G(½+it)|²w dt over G(s) = ∫_{[0,L]}e^{−us}dν(u) (closure of B_λ), w = 1/(2π(¼+t²)), exactly as for ζ
(Mellin of 1 − Σc_jD_{θ_j}a_χ is (1 − L(s,χ)Σc'_jθ_j^s)/s; no sign, no pole). Z_χ(T) := zeros of L(s,χ) with |γ| ≤ T, nodes z_ρ = γ − i(β − ½).
Inputs, q fixed (constants depend on q):
 (Zχ) GRH: S(t,χ) ≪ log q(|t|+2)/log log q(|t|+3) (Titchmarsh 14.13's proof verbatim: it uses only the Hadamard/explicit formula for
   log L and the Dirichlet series). One-sided count: N⁺(T,χ) := #{0 < γ ≤ T} = (T/2π)log(qT/2πe) + S(T,χ) + O(1) (Davenport ch. 16; the
   zeros with −T ≤ γ < 0 are the conjugates of those of L(s,χ̄), counted by N⁺(T,χ̄)); so every [t, t+h] or [−t−h, −t], h ≤ 1/log log t, holds
   ≤ A₁ log qt/log log qt ordinates. Replaces DbarNB (Z); §3 grouping then gives M_g ≤ M*(|γ_g|) := A₁ log(q|γ_g|)/log log(q|γ_g|).
 (Cχ) Lemma C for χ (unconditional): every [T, T+1] and [−T−1, −T], T ≥ 2, contains t with |L(σ+it,χ)|^{−1} ≤ (qT)^A on −1 ≤ σ ≤ 2.
   Titchmarsh 9.7's proof verbatim: log|L(s,χ)| = Σ_{|t−γ|≤1}m_ρ log|s−ρ| + O(log qT) on −1 ≤ σ ≤ 2 (from L′/L = Σ_{|t−γ|<1}m_ρ(s−ρ)^{−1} + O(log q(|t|+2)), distinct ρ,
   Davenport ch. 16); log|s−ρ| ≥ log|t−γ|; the average over t ∈ [T,T+1] of Σ_{|t−γ|≤1}m_ρ log|t−γ| + Σ_{|t+γ|≤1}m_ρ log|t+γ| is ≥ −A log qT
   (O(log qT) zeros in [T−1, T+2] ∪ [−T−2, −T+1]), and both sums are ≤ 0: one t is good at +t and at −t simultaneously (symmetric cutoff T₂).
 (Fχ) Convexity/functional equation: |L(½+it,χ)| ≪_ε (1+|t|)^{¼+ε} (so L(s,χ)/s ∈ L² of the line and ‖L R_tr‖_w is finite, as |ζ| ≪ t^{1/6}
   for ζ); |L(w,χ)|, |L′(w,χ)| ≤ C(2+|Im w|)² on −1 ≤ Re w ≤ 2; |1/L| ≤ ζ(3/2)/ζ(3) on Re w = 3/2; on Re w = −½,
   |Im w| ≥ 1: L(w) = γ_a(1−w)W^{±1}L(1−w, χ̄), |γ_a(3/2−it)| ≍ q|t| so |1/L| ≤ C; for σ < −¾, |1/L(s,χ)| ≤ C off |s + a + 2k| < ½.
Steps verbatim (L(s,χ) enters only through real ordinates, (Zχ), (Cχ), (Fχ)): dual form §1; Lemmas A, B, Λ, D; Theorems R″, R′ and the
Corollary; Lemma E (with (c6)); §6 tilings, family independence, (i)–(v); REPAIR (3i) window Φ, Ψ = 1 − E_low; §7 closure (PWS + Bochner
integral; u ↦ e^{−us}L(s,χ)/s is L²-continuous since θ ↦ D_θa_χ is); Claim F (any finite set of real nodes). Changes:
 (c1) No pole: L[s,w] := (L(s,χ) − L(w,χ))/(s − w) is entire in (s,w); §6(ii)–(iii) hold for every s, and the s = 1 exceptions (R_high(1), G(1))
   disappear (G(1) mattered only for Nyman's constrained variant).
 (c2) Trivial zeros s_k = −a − 2k, k ≥ 0 (simple; s₀ = 0 for even χ): R_tr(s) := Σ_{k≥0} r_k e^{−(L−1)(s−s_k)}/(s−s_k), r_k := Ψ(s_k)Φ(s_k)/L′(s_k,χ).
   REPAIR (3ii) with y_k := a + 2k + ½ for 2k + ½: |Ψ(s_k)| ≤ 1 + C√(k+1)e^{L′y_k}, |Φ(s_k)| ≤ Ce^{(2/T₁)y_k}, |L′(s_k,χ)| = ½(q/π)^{y_k}Γ(k+a+½)k!|L(1+a+2k,χ)|
   (from L(s) = γ_a(1−s)L(1−s), 1/Γ((s+a)/2) ≈ (−1)^k k!(s−s_k)/2): on the line each term ≤ C√(k+1)e^{−(√L−1−2/T₁)y_k}/(|L′(s_k)||s−s_k|), so
   ‖L R_tr‖²_w ≪ e^{−(2a+1)(√L−2)} = o(L^{−2}) (worst: even χ, the zero at s = 0). R_tr is meromorphic (no spectrum claimed); off the discs
   |s − s_k| < ½ its terms grow at most like e^{(L−1)max(0,−σ)}, which is all (iv) uses; PWS is applied only to the entire G.
 (c3) Growth (iv): σ > 7/4: |F₀| = |ΨΦ/L| ≤ C(1+|s|)√L; σ < −¾ off the discs: (Fχ); inside the discs |s − s_k| < ½: bound F₀, R_high, R_tr on
   the circles, then the maximum principle for the entire G (not for its meromorphic summands), as for ζ.
 (c4) Zeros above T₂ satisfy 0 < β < 1 (L ≠ 0 on Re s = 1; on Re s = 0 only s₀ = 0 for even χ, not high): tilings hold.
 (c5) Low zeros below 14: Lemma B (needs |γ_g| ≥ 14) is used only for high groups (|γ| > H ≥ 14) in R′ and Lemma E; low zeros enter only
   through Claim F. (DbarNB's R″ as stated covers all groups; R′ + Claim F suffice for the limit; R″ for the high part is unchanged.)
 (c6) Central zero ρ = ½ (γ = 0; GRH allows it; for real χ conjecturally absent, Chowla): Claim F holds (Gram (¼+γ²)ΔMΔ*, relative error O(1/L));
   in Lemma E replace ‖p_lo′‖² ≥ c(H)Σ|a_{ρk}|²ℓ^{2k+1} (false for e_{ρ0} ≡ 1) by N(p_lo) ≥ c(H)Σ|a_{ρk}|²ℓ^{2k+1} (the ⟨·,·⟩_N-Gram of ℓ^{−k−½}e_{ρk}
   tends to blockdiag((¼+γ_ρ²)Hilb_{m_ρ}), nondegenerate at γ = 0); the cross terms ∫p_lo′p̄_hi′ and ¼∫p_lop̄_hi are bounded as there (|Ω| ≥ g_H > 0):
   θ_L = C(H)L^{−3/4+o(1)} unchanged.
Result (3v)χ, unconditional: D_χ(λ)² ≤ (1 + L^{−½})U_{Z_χ(T₂)}(L − √L) + O(L^{−3/2}); under GRH (Corollary of R′ + Claim F):
lim L·U_{Z_χ(T₂)}(L − √L) = Σ_ρ m_ρ²/|ρ|². Hence limsup L·D_χ² ≤ Σ_ρ m_ρ²/|ρ|².

## 7. Theorem and equivalent form
Theorem. χ primitive real, GRH for L(s,χ) ⇒ lim_{λ→0} D_χ(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|² (distinct zeros, multiplicities m_ρ; no condition on L′(ρ,χ)).
Proof. limsup ≤: §6. liminf ≥: Theorem Bχ (§5; under GRH every zero is on the line). ∎
Constant. Davenport ch. 12 (sums over distinct zeros with multiplicity weights): ξ′/ξ(s,χ) = B(χ) + Σ_ρ m_ρ(1/(s−ρ) + 1/ρ),
Re B(χ) = −Σ_ρ m_ρ Re(1/ρ) (unconditional; the real-part sum converges absolutely, Σ 1/ρ need not); B(χ) = ξ′/ξ(0,χ)
= −ξ′/ξ(1,χ̄) (functional equation), conj B(χ) = B(χ̄). So Σ_ρ m_ρ Re(1/ρ) = Re ξ′/ξ(1,χ) = ½log(q/π) + ½ψ((1+a)/2) + Re L′/L(1,χ). Under GRH
Re(1/ρ) = 1/(2|ρ|²): Σ_ρ m_ρ/|ρ|² = C_χ := 2Re L′/L(1,χ) + log(q/π) + ψ((1+a)/2) (ψ(½) = −γ − 2log 2, ψ(1) = −γ). χ₋₄: 0.1555679799 (brief: 0.15556798),
χ₅: 0.1565569536 (logs/burnol_*_CDE*.log, part D; numerics in REPORT Evidence).
Corollary (equivalent form). lim_{λ→0} D_χ(λ)² log(1/λ) = C_χ ⟺ GRH for L(s,χ) and all its zeros simple.
Proof. ⇐: Theorem with m_ρ = 1. ⇒: if GRH fails, D_χ(λ) ≥ c_χ > 0 (Theorem Bχ), so L·D_χ² → ∞; if GRH holds and some m_ρ ≥ 2, the limit is
Σ m_ρ²/|ρ|² > Σ m_ρ/|ρ|² = C_χ. ∎
Uniformity in q: not obtained. L₀ depends on χ through Claim F (the zeros below H, their gaps g_H) and the constants of (Zχ), (Cχ); the
theorem is a limit for each fixed χ.

## 8. Complex primitive χ (brief item (4)): no step fails; no (χ, χ̄)-pairing needed
Lower bound: §§1–5 were written for all primitive χ: the phase is V_χ = W(χ)^{−1}(s/(1−s))γ_a(s), explicit because conj L(s,χ) = L(1−s,χ̄) on
the line and L(1−s,χ̄) = W(χ)^{−1}γ_a(s)L(s,χ); W(χ)^{−1} is a unimodular constant; (O) produces the zeros of L(s,χ) itself (not of L(s,χ̄)).
Upper bound: the ordinates of L(s,χ) are no longer symmetric and conj L(s̄,χ) = L(s,χ̄) ≠ L(s,χ). DbarNB uses real coefficients only in
R_high = Σ_n(R_n + R̄_n), R̄_n(s) := conj R_n(s̄) = the residue sum of F₀ over the mirror rectangle. Repair: (Cχ) gives heights t_n, t′_n good at
+t and −t simultaneously; Q_n^± := [−½, 3/2] × ±[t_n, t_{n+1}] (and the t′-family, narrow/wide); R_n^± := Σ_{ρ∈Q_n^±}Res_{w=ρ}[F₀(w)B(s−w)/(s−w)]
(actual residues of F₀ = ΨΦ/L(·,χ)); R_high := Σ_n(R_n^+ + R_n^−). Every step of §6(i)–(v) is a per-rectangle estimate summed in absolute value:
(i)–(ii) per rectangle; (iii) ‖L R_high‖²_w ≤ C Σ_±Σ_n t_{n+1}^{2A+8}Le^{4L}e^{−2c(t_n/L²)^{1/2}} ≤ e^{−c′L²} (T₂ = t_{n₀} ∈ [L⁶, L⁶+4] as for ζ;
the same bound for each half: |1/L| ≤ (qt)^A on the sides of Q_n^±, |Φ(w)| ≤ Ce^{−c(|Im w|/T₁)^{1/2}} in both half-planes, |L[s,w]| by convexity);
(iv) for σ ∈ [−¾, 7/4] and ±t > T₂ + 4 use the wide family of that half-plane whose horizontal sides are ≥ ¼ from t: G = −Σ_m I^{w,±}_m − Σ_m I^{w,∓}_m − R_tr
(s lies in no rectangle of the opposite half, so there R_m^∓ = I^{w,∓}_m by §6(i)), with |I^{w,±}_m(s)| ≤ C t_{m+1}^{A+1}√L e^{2L}e^{−c(t_m/L²)^{1/2}}
·e^{(L−1)(2+max(0,−σ))}min(1, C(L|t ∓ t_m|)^{−2}) (DbarNB §6(iv) verbatim, per rectangle) and |I^{w,∓}_m(s)| ≤ the same with min(1, C(L(|t| + t_m))^{−2})
(separation ≥ |t| + t_m − 1; same boundary bounds on ∂Q_m^∓); |t| ≤ T₂ + 4: compact; σ ∉ [−¾, 7/4]: the narrow families of both halves, (Fχ) and the circles |s − s_k| = ½. Hence for all s
|G(s)| ≤ C_L(1+|s|)^C e^{(L−½)max(0,−σ)}, and PWS gives supp ⊂ [0, L − ½]. No cancellation between R_n and R̄_n is used for ζ either, and the high part enters (v) only as an
O(e^{−c′L²}) error: the constant Σm²/|ρ|² comes from U_{Z_χ(T₂)} alone (§§1–5 of DbarNB, any real node set). Lemma E's reflection p ↦ conj p(ℓ − ·) keeps each real frequency (no symmetry used). Claim F, R′, D: any real
node set. GRH for χ ⟺ GRH for χ̄. Theorem and Corollary hold for every primitive χ (q ≥ 3) with C_χ = 2Re L′/L(1,χ) + log(q/π) + ψ((1+a)/2).
Repair (orchestrator's §8 referee, 20 Nov: (iv) at height ±t kept only that half's contour sum): R_high = R_high^+ + R_high^−, so G = F₀ − R_high − R_tr
at ±t > T₂ + 4 equals −Σ_m I^{w,±}_m − Σ_m I^{w,∓}_m − R_tr; the opposite-half terms obey the same estimates with separation |t| + t_m, so the
growth bound and PWS stand. ROUTINE; no (χ, χ̄)-pairing obstruction.

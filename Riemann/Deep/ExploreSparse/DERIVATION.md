# explore61 — derivations (exact statements; numerics in REPORT.md)

H = L²(0,∞), ρ_θ(x) = {θ/x}, χ = 1_(0,1). V_N = span{ρ_{1/k} : k ≤ N}, C_λ = closed span{ρ_θ : λ ≤ θ ≤ 1}, λ = 1/N, so V_N ⊂ C_λ.
d_N = dist(χ, V_N), D = D(λ) = dist(χ, C_λ). h* = P_C χ, r* = χ − h* ⊥ C_λ; h_N = P_V χ, r_N = χ − h_N ⊥ V_N. P := d_N² − D².
Harmonic intervals I_m = (1/(m+1), 1/m), m ≥ 1.

## 1. Bagchi's space W1 (known: Bagchi, math/0607733, §§1, 4; restated in L²(0,∞))
W1 := {f ∈ H : f = a_m + b/x on each I_m, f = b/x on (1,∞), one b for all m}. Closed (two-dimensional on each piece, b shared).
- ρ_{1/k} ∈ W1: on I_m, 1/(kx) ∈ (m/k, (m+1)/k) contains no integer, so {1/(kx)} = (1/k)(1/x) − ⌊m/k⌋; on (1,∞) it is (1/k)/x.
  Equivalently ρ_{1/k} = (1/k)ρ_1 + Σ_m {m/k} 1_{I_m} (Bagchi's g_k ∈ M, γ_k = ({m/k})_m). χ ∈ W1 (a_m = 1, b = 0).
- Hence V_N ⊂ W1 for every N. W1 = ℂρ_1 + PC, a direct sum that is NOT orthogonal (⟨ρ_1, 1_{I_m}⟩ = log(1 + 1/m) − 1/(m+1) > 0),
  PC = closed span{1_(0,1/n)} = piecewise constants on the I_m.
- RH ⟺ V_∞ := closure(∪V_N) = W1. (⇒) Báez-Duarte: χ ∈ V_∞; the contractions f ↦ f(k·) (norm k^{-1/2}) map V_∞ into itself,
  so 1_(0,1/k) ∈ V_∞, PC ⊂ V_∞, and ρ_1 ∈ V_∞. (⇐) χ ∈ W1 = V_∞ ⊂ closure ∪C_λ, so RH (Nyman–Beurling). Bagchi's form: γ_k total.

## 2. Sampling defect (barrier for MZ/frame transfer; elementary)
For θ ∉ {1/k}, ⌊θ/x⌋ jumps by 1 at x = θ (inside some I_m), so ρ_θ ∉ W1 and ω(θ) := dist(ρ_θ, W1) > 0. Since V_N ⊂ W1:
  dist(ρ_θ, V_N) ≥ ω(θ) for every N;  S_N ⪰ Ω  (S_N = Gram of P_{V_N}^⊥ρ_θj, Ω = Gram of P_{W1}^⊥ρ_θj) for every finite θ-set.
Closed form (s1_split.w1_perp_gram, s9_cert): ρ_θ = θρ_1 + [θ⌊1/x⌋ − ⌊θ/x⌋]; on I_m the bracket is constant minus Σ 1[x < θ/n]
(jumps at θ/n ∈ I_m, fraction q = (θ/n − 1/(m+1))/|I_m|). With β the ρ_1-coefficient, a_m the interval means:
  ω(θ)² = A − B²/C′, A = Σ |I_m|(q − q²), B = −Σ (log(x_j/l) − q log(r/l)), C′ = 1 + Σ_m (1 − m(m+1) log²(1 + 1/m)),
summed over jumps x_j ∈ I_m = (l, r) and all m; dropping m > M lowers it (every term of the β-quadratic is ≥ 0): certified lower bounds.

## 3. Exact split and Lemma B (new as far as found)
(a) P = ‖h* − h_N‖² (Pythagoras: r_N = r* + (h* − h_N), r* ⊥ C ∋ h* − h_N), and h_N = P_V h* (V ⊂ C).
(b) Split by W1 (h_N ∈ V ⊂ W1): P = R + Q, R := ‖P_{W1}^⊥ h*‖² = dist(h*, W1)², Q := ‖P_{W1}h* − h_N‖² = dist(P_{W1}h*, V_N)²
    (the nearest point of V_N to P_W h* is P_V P_W h* = P_V h* = h_N).
(c) ⟨P_W r*, P_W h* − h_N⟩ = ⟨r*, P_W h*⟩ (r* ⊥ h_N) = −⟨r*, P_W^⊥ h*⟩ (r* ⊥ h*) = −⟨P_W^⊥ r*, P_W^⊥ h*⟩ = ‖P_W^⊥ h*‖² = R,
    since P_W^⊥ r* = −P_W^⊥ h* (χ ∈ W1). ‖P_W r*‖² = D² − R. Cauchy–Schwarz: **R² ≤ (D² − R)·Q** (Lemma B).
(d) Adding RQ to both sides of (c): R(R + Q) ≤ D²Q, i.e. **R·P ≤ D²·Q** (always). If R > 0: P ≤ (Q/R)·D², i.e.
    d_N² ≤ (1 + Q/R)·D(1/N)². R = 0 happens (e.g. N = 1, where C = V_1), so the form used everywhere below is:
    **if Q ≤ ε·R for some ε ≥ 0, then P ≤ ε·D²** (R > 0: R·P ≤ D²Q ≤ εRD²; R = 0: then Q = 0 and P = R + Q = 0).
    Along N: Q_N ≤ ε_N R_N with ε_N → 0 (Q_N = 0 where R_N = 0) gives P_N ≤ ε_N D_N². (Reviewer's and referee's GAP, fixed.)
    The same holds for every target g ∈ W1 in place of χ, and for every closed C with V ⊂ C (e.g. the head-hybrid spans of s1).

## 4. Theorem C (discretization through W1; the MZ power-function bound in the node space)
For every h ∈ H: d_N ≤ ‖χ − P_W h‖ + dist(P_W h, V_N) ≤ ‖χ − h‖ + dist(P_W h, V_N) (P_W contracts toward χ ∈ W1).
For h = Σ_{k≤N} c_kρ_{1/k} + ∫ρ_θ dν(θ) (ν a finite Borel measure on [λ,1]; θ ↦ ρ_θ is continuous into H with ‖ρ_θ‖² = θ‖ρ_1‖²,
so the integral is a Bochner integral and the bounded maps P_W, I − P_V pass inside it): dist(P_W h, V_N) ≤ ∫ π^W_N(θ) d|ν|(θ), π^W_N(θ) := dist(P_{W1}ρ_θ, V_N) = (π_N(θ)² − ω(θ)²)^{1/2}.
Contrast: the full-space power function π_N(θ) = dist(ρ_θ, V_N) ≥ ω(θ) never decays; π^W_N → 0 under RH (V_∞ = W1).
Corollary (RH, with dbarnb's THEOREM D² log(1/λ) → C): if continuous candidates h_N with ‖χ − h_N‖² ≤ (C + o(1))/L have non-integer
measures ν_N with ∫π^W_N d|ν_N| = o(L^{−1/2}) (or ∫π_N d|ν_N| = o(L^{−1/2})), then d_N² log N → C (BBLS's discrete conjecture).

## 5. Zero formula (numerical, not proved)
Expected from Burnol's zero vectors (⟨g, Y_ρ⟩ ≈ G(ρ)(ρ−1)/ρ, ‖Y_ρ‖² ≈ L; not checked here beyond the numerics), for a target g with
Mellin transform G(s) = ∫g x^{s−1}dx: liminf L·dist(g, C_λ)² ≥ Σ_ρ|G(ρ)|² (for χ: Σ1/|ρ|² = C, Burnol's theorem). For g = P_W ρ_θ ∈ W1: G(ρ) = 𝓜[ρ_θ](ρ) − 𝓜[P_W^⊥ρ_θ](ρ) = −𝓜[P_W^⊥ρ_θ](ρ) (ζ(ρ) = 0). Observed (s5, s6, s10):
  L·dist(P_W ρ_θ, V_N)² = Σ_ρ |𝓜[P_W^⊥ρ_θ](ρ)|² (1 ± .06), and along the optimizer  L·Q = Σ_ρ |𝓜[P_W^⊥h*](ρ)|² (1 ± .03).
𝓜 of the roughness, per jump x_j ∈ (l, r): −[(x_j^s − l^s) − q(r^s − l^s)]/s; β-part −β[Σ_m((r^{s−1} − l^{s−1})/(s−1) − avg_m(r^s − l^s)/s) + 1/(1−s)].
Consequence: Q/R = Z/L with Z := Σ_ρ|𝓜g(ρ)|²/‖g‖² (zero visibility of the roughness g = P_W^⊥h*), and Lemma B with cos := R/√((D²−R)Q):
  R = cos²·Z·(D² − R)/L, so L²P → cos²·Z·C if Z, cos converge (measured Z = 1.36–1.70, cos = .92–.96).

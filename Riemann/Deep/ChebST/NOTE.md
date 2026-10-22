# ChebST: Sato–Tate in every ideal class, for almost all discriminants (proofs)

Date 2026-10-21. Brief: `.scratch/orca/chebst-brief.md`. Scratch: `.scratch/chebst/` (deleted at the end of the run; the outputs
quoted here are reproduced in the text): `cms.py` (the polynomials of §1 built by Hermite interpolation and checked on grids),
`cms_table.py` (degrees for sample intervals; the uniform bound of Corollary M5 on 287 intervals), `params.py` (every
parameter inequality of §2 for 2 ≤ m ≤ 80, routes F, H, H′, at A = threshold + 10⁻⁴, 10⁻², 0.1, 0.25: 9908 checks, 0 failures),
novelty scripts ([NOVELTY.md](NOVELTY.md)).
References: [RS] = [RSFam NOTE](../RSFam/NOTE.md) (§0 notation, (R1)–(R4), (R1′), Lemmas A, B, Theorem C, Proposition F, §5, §5.6, §6,
§7); [C] = [Combos NOTE](../Combos/NOTE.md) (Lemmas 2, 3, 5; (F1)–(F4)); [R] = [RCApps NOTE](../RCApps/NOTE.md) ((I4)–(I6), §2,
Theorem 3); [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md) (Theorems 2′, 3′, 4, Corollary 5); [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md)
(§2: the squarefree zero detector); [D2] = [Density2 NOTE](../Density2/NOTE.md) (Lemma 5, Theorem A, (2.1)).
HT = Humphries–Thorner, Compos. Math. 160 (2024), arXiv 2103.05634v2, Theorem 4.1 (TeX label thm:main_theorem_1, read in the
source); PT = Pascadi–Thorner, *Large sieves for GL_n and applications*, arXiv 2508.14888v2 (2025, preprint), Theorem 1.1 and
Corollary 1.2 (labels thm:large_sieve2, cor:large_sieve2, read in the source). Paper proofs, not refereed, no Lean, no 2026
zero-free half-planes.

## Verdicts

| Item | Result | Status |
|---|---|---|
| minorants (§1) | for every interval J ⊂ [−1,1] of the variable x = cos θ and every n ≥ 1: polynomials P⁻ ≤ 1_J ≤ P⁺ on [−1,1] of degree ≤ 2n with ∣∫P^± dμ_ST − μ_ST(J)∣ ≤ λ_n(α) + λ_n(β) (α, β the endpoints inside (−1,1)), λ_n(cos θ) = 4sin³θ/((2n+3)sin θ − sin(2n+3)θ) ≤ 2/(n+1); positive mean once n ≥ ⌈37.7/L⌉, L = length of the θ-interval | classical (Chebyshev–Markov–Stieltjes), written out with explicit constants |
| Proposition T_m, every m ≥ 2 (§2) | for A > 𝒜(m): ∣h S_C^{(m)}∣ ≤ η₀ν_D H for all classes C, x ≥ ∣D∣^A, outside ≤ X^{1−δ/8} discriminants | route F (A > (m+7)/3) is [RS] §5 for general m ([RS] §7 sketched it); routes H (A > max{5, (m+3)/3}, published inputs only) and H′ (A > max{11/3, (m+3)/3}, uses PT Cor 1.2) are new |
| Proposition F′ (§3) | Σ_{∣D∣~X} Σ_{χ∈Ĝ_c} N(σ; t, Sym^m f ⊗ θ_χ) ≪ X^{δ′}[X^{1+3(1−σ)} + X^{3/4+(11/2)(1−σ)} + (X^{m+6}T^{2m+2})^{(1−σ)/(3−2σ)}] (unit window); with HT alone X^{1+5(1−σ)} replaces the first two terms; Proposition F″: (X^{m+12}T^{2m+2})^{(1−σ)/(3−2σ)} with the GL₂ sieve only (used near σ = 1) | new; against X^{(m+10)(1−σ)/(3−2σ)} ([RS] Prop F): per-D mollifiers (family size X^{3/2}) in class II; class I by the per-D sieve (N ≤ X^{3/2}), PT's Hölder sieve (X^{3/2} ≤ N ≤ X^{11/4}) and HT/PT above |
| Theorem 1 (§4) | every interval J with nonempty interior: for A > A(J) := 𝒜(M*(J)), almost all D (both signs) have primes 𝔭 of degree one with cos θ_p ∈ J in every (narrow) class and every (x, x + x∣D∣^{−κ}], x ≥ ∣D∣^A; A(J) ≤ max{11/3, 25.2/L + 1.7} | **new** statement (no source found); for M*(J) ≤ 4 the threshold is that of [RS] Theorem D_M |
| Theorem 2 (§5) | Sato–Tate inside every class: discrepancy ≤ 5/(n+1) (intervals), ≤ 3/(n+1) (half-lines) for A > 𝒜(2n); for A > 11/3: ≤ 10/(3A−3) | **new** |
| Theorem 3 (§6) | sectors for D < 0 (Hecke angle in arcs of length ≥ ∣D∣^{−κ}) with the same thresholds (8/3 instead of 9/4 for degree one) | new; proof by listed changes |
| uniformity in m (§7) | non-explicit in m: HT's constant for π′ = Sym^m f (route F only), BTZ's B_m, Lemma B's C(ε, m); m-free: HT/PT for π′ = 𝟙 (routes H, H′); explicit: Thorner's zero-free region, conductors, coefficient bounds, c_m = 10⁷(m+1)⁴ | — |

Thresholds (minorant degree M ↦ 𝒜(M); the published densities would give (M+8)/3, [RS] §7):

| M | 0 | 1 | 2 | 3 | 4–8 | ≥ 8 |
|---|---|---|---|---|---|---|
| 𝒜(M) | 3/2 | 9/4 | 3 | 10/3 | 11/3 | (M+3)/3 |
| route | [R2] | [RS] (Thm C) | F | F | F (M = 4), H′ | H′ |
| without PT (HT only) | 3/2 | 9/4 | 3 | 10/3 | (M+7)/3 (M ≤ 8) | 5 (8 ≤ M ≤ 12), (M+3)/3 (M ≥ 12) |

## 0. Setup

- Notation of [RS] §0, §5: f ∈ S_k(Γ₀(N)) a newform without CM, k ≥ 2, trivial nebentypus; λ_f(p) = a_f(p)p^{−(k−1)/2} = 2cos θ_p,
  0 ≤ θ_p ≤ π (p ∤ N); **x_p := cos θ_p**. μ := μ_ST, (2/π)sin²θ dθ on [0, π], i.e. (2/π)√(1−x²)dx on [−1, 1].
  U_m(cos θ) = sin((m+1)θ)/sin θ: orthonormal for μ, ∣U_m∣ ≤ m+1 on [−1,1], ∫U_m dμ = δ_{m,0}.
- 𝔇, 𝔇(X), K = Q(√D), G, h, Ĝ, Ĝ₂, Ĝ_c, 𝔉(D), ε_d(C), ν_D, b_m(𝔞), L_m(s,χ) = L(s, Sym^m f ⊗ θ_χ), N_m(σ;t,χ),
  N_D^{c,m}(σ,T), S_C^{(m)} := Σ_{𝔞∈C}Λ_K(𝔞)b_m(𝔞)w(N𝔞), w(u) := w₀((u−x)/H): as in [RS] §0, §5 (w₀ ∈ C_c^∞((0,1)), 0 ≤ w₀ ≤ 1,
  w₀ = 1 on [¼, ¾]). For D > 0: narrow classes, h⁺, Ĝ⁺ ([RS] §6, [R2] §3). b_m(𝔭) = U_m(x_p) for degree-one 𝔭, p = N𝔭 ∤ N ((R1)).
- Π_m := Sym^m π_f is cuspidal automorphic on GL_{m+1}(A_Q) for **every** m ≥ 1: Newton–Thorne, Publ. Math. IHÉS 134 (2021),
  I (level one and a class including semistable curves) and II ("the existence of the symmetric power liftings of all holomorphic
  Hecke eigenforms"); cuspidal because f has no CM; tempered everywhere and self-dual ([RS] §0). [RS] states (R1)–(R4), Lemmas A, B,
  Theorem C and Proposition F for every m ≥ 1, and these hold as stated; only their constants depend on m (§7).
- Intervals: for a θ-interval I = [a, b] ⊂ [0, π] put J := cos(I) = [cos b, cos a], L := b − a, ℓ := ∣J∣ ≤ L. Endpoints may be open or
  closed; every minorant below is ≤ 1 on the open interval, so the type of an endpoint never matters.
- ∣x_p∣ < 1 for all but finitely many p ∤ N: k is even (trivial nebentypus), so a_f(p)² = 4p^{k−1} would put √p in the coefficient
  field of f, which has finitely many quadratic subfields. Those finitely many primes lie below x once ∣D∣ ≥ D₀, so only the values of
  a test polynomial on (−1, 1) matter.
- **M\*(J)**, for an interval J ⊂ [−1,1] with nonempty interior: the least M ≥ 0 such that some real polynomial P of degree ≤ M has
  P ≤ 1_J on (−1, 1) and ∫P dμ > 0 (finite by Corollary M5; M\*(J) = 0 iff J ⊃ (−1, 1)).
- **𝒜(M)**: 𝒜(0) := 3/2, 𝒜(1) := 9/4, 𝒜(M) := min{(M+7)/3, max{11/3, (M+3)/3}} for M ≥ 2 (using PT Corollary 1.2; with
  published inputs only, replace 11/3 by 5). 𝒜 is nondecreasing (`params.py`).
- e(σ) := (1−σ)/(3−2σ). C(·), c(·): effective constants depending only on f and the displayed parameters.

## 1. Minorants and majorants for the Sato–Tate measure

K_n(s,t) := Σ_{k=0}^{n}U_k(s)U_k(t) and λ_n(t) := 1/K_n(t,t) (the Christoffel function of μ). From Σ_{j=1}^{n+1}sin²(jθ) =
(n+1)/2 + ¼ − sin((2n+3)θ)/(4sin θ) and Σ_{j≤n+1}j² = (n+1)(n+2)(2n+3)/6:

  (1.1)  λ_n(cos θ) = 4sin³θ/((2n+3)sin θ − sin((2n+3)θ))  (0 < θ < π),   λ_n(±1) = 6/((n+1)(n+2)(2n+3)).

**Lemma M1 (quadrature with a prescribed node).** Let n ≥ 1, t ∈ R, ω_t(x) := U_n(t)U_{n+1}(x) − U_{n+1}(t)U_n(x), and let
y_0 < ⋯ < y_r be the real zeros of ω_t. Then r = n if U_n(t) ≠ 0 and r = n−1 if U_n(t) = 0; t is one of the y_j; λ_n(y_j) > 0;
and ∫q dμ = Σ_j λ_n(y_j)q(y_j) for every polynomial q of degree ≤ 2r.

*Proof.* If U_n(t) ≠ 0, ω_t is a nonzero multiple of U_{n+1} − cU_n, c := U_{n+1}(t)/U_n(t); such a quasi-orthogonal polynomial has n+1
real simple zeros (they interlace with those of U_n; Szegő, *Orthogonal polynomials*, §3.3), and ω_t(t) = 0. If U_n(t) = 0 then
U_{n+1}(t) ≠ 0 (consecutive orthogonal polynomials have no common zero) and ω_t = −U_{n+1}(t)U_n: its zeros are the n zeros of U_n,
t among them. *Exactness.* Write q = ω_t g + ρ with deg ρ ≤ r and deg g ≤ 2r − r − 1 ≤ n − 1; ω_t is orthogonal to polynomials of
degree ≤ n−1 (it is a combination of U_n and U_{n+1}), so ∫q = ∫ρ, and the interpolatory rule on r+1 nodes integrates ρ exactly.
*Weights.* The weight at y_j is ∫ℓ_j²dμ > 0 (ℓ_j the Lagrange basis, deg ℓ_j² = 2r). Christoffel–Darboux,
K_n(x,y) = (U_{n+1}(x)U_n(y) − U_n(x)U_{n+1}(y))/(2(x−y)), gives K_n(y_i, y_j) = 0 for i ≠ j: if U_n(t) ≠ 0 both satisfy
U_{n+1}(y) = cU_n(y); if U_n(t) = 0 both are zeros of U_n. Exactness applied to K_n(·, y_j)² (degree 2n, or 2n−2 = 2r when U_n(y_j) = 0,
since then K_n(·,y_j) = K_{n−1}(·,y_j)) gives K_n(y_j,y_j) = Σ_i w_iK_n(y_i,y_j)² = w_jK_n(y_j,y_j)², so w_j = λ_n(y_j). □

**Lemma M2 (one-sided polynomials; Chebyshev–Markov–Stieltjes).** For n ≥ 1 and t ∈ R there are real polynomials P_t⁻, P_t⁺ of degree
≤ 2n with P_t⁻ ≤ 1_{(−∞,t)} and P_t⁺ ≥ 1_{(−∞,t]} on all of R, ∫P_t⁻dμ ≥ 0, and ∫P_t⁺dμ − ∫P_t⁻dμ = λ_n(t). Consequently
  μ([−1,t]) − λ_n(t) ≤ ∫P_t⁻dμ ≤ μ([−1,t)),   μ([−1,t]) ≤ ∫P_t⁺dμ ≤ μ([−1,t)) + λ_n(t).

*Proof.* Nodes y_0 < ⋯ < y_r of Lemma M1, y_ν = t. Let P⁻ be the Hermite interpolant of degree ≤ 2r with P⁻(y_j) = 1 for j < ν,
P⁻(y_j) = 0 for j ≥ ν, (P⁻)′(y_j) = 0 for j ≠ ν (2r+1 conditions). If ν = 0, P⁻ ≡ 0. Otherwise P := P⁻ is not constant and P′ has
degree ≤ 2r − 1 and at least 2r − 1 distinct zeros: the r points y_j (j ≠ ν) and, by Rolle, one point strictly inside each of the
r − 1 intervals (y_j, y_{j+1}) whose end values agree (all except (y_{ν−1}, y_ν)). These are therefore all its zeros, each simple, and
P′ changes sign at each. P′ has no zero in (y_{ν−1}, y_ν) and P goes from 1 to 0 there, so P′ < 0 on it; following the sign changes
outwards, every y_j (j ≠ ν) is a strict local maximum, every Rolle point a strict local minimum, P′ > 0 on (−∞, y_0) and P′ < 0 on
(y_r, ∞) (on (t, ∞) if ν = r). Hence P ≤ 1 on (−∞, t) and P ≤ 0 on [t, ∞). P⁺ is defined by P⁺(y_j) = 1 (j ≤ ν), 0 (j > ν),
(P⁺)′(y_j) = 0 (j ≠ ν) (P⁺ ≡ 1 if ν = r); the same count, with P′ < 0 on (y_ν, y_{ν+1}), makes the y_j (j ≠ ν) local minima, so
P⁺ ≥ 1 on (−∞, t] and P⁺ ≥ 0 on (t, ∞). By Lemma M1, ∫P⁻ = Σ_{j<ν}λ_n(y_j) ≥ 0 and ∫P⁺ = Σ_{j≤ν}λ_n(y_j), whose difference is
λ_n(t). The displayed inequalities follow by integrating P⁻ ≤ 1_{(−∞,t)} ≤ 1_{(−∞,t]} ≤ P⁺ against μ. □

`cms.py` builds P_t^± for n ∈ {1,2,3,5,8,12} and 8 values of t (96 polynomials) and confirms P_t⁻ ≤ 1_{(−∞,t)}, P_t⁺ ≥ 1_{(−∞,t]}
on grids of [−2.5, 2.5] (worst violation 1.1·10⁻¹⁰, rounding), Σ_j λ_n(y_j) = 1, and ∫P_t^±dμ = Σ λ_n(y_j) to 10⁻⁹ except in two
ill-conditioned cases (n = 12, error 2·10⁻⁴ in the double-precision interpolation; the lemma does not depend on this check).

**Lemma M3 (intervals).** Let −1 ≤ α < β ≤ 1, J an interval with endpoints α, β, n ≥ 1, and ϵ_α := λ_n(α) if α > −1, ϵ_α := 0 if
α = −1; ϵ_β := λ_n(β) if β < 1, ϵ_β := 0 if β = 1. Put P_J⁻ := P_β⁻ − P_α⁺ and P_J⁺ := P_β⁺ − P_α⁻, where the α-term is omitted if
α = −1 and P_β^± is replaced by 1 if β = 1. Then on [−1, 1]: P_J⁻ ≤ 1_{(α,β)} ≤ 1_J ≤ 1_{[α,β]} ≤ P_J⁺; deg P_J^± ≤ 2n;
∣∫P_J^±dμ − μ(J)∣ ≤ ϵ_α + ϵ_β; ∫∣P_J^±∣dμ ≤ 4; and the coefficients c_m^± := ∫P_J^±U_m dμ in P_J^± = Σ_{m≤2n}c_m^±U_m satisfy
c_0^± = ∫P_J^±dμ, ∣c_m^±∣ ≤ 4(m+1), Σ_{m=1}^{2n}∣c_m^±∣ ≤ 4n(2n+3).

*Proof.* On [−1,1], P_β⁻ − P_α⁺ ≤ 1_{(−∞,β)} − 1_{(−∞,α]} = 1_{(α,β)}, and 1 − P_α⁺ ≤ 1_{(α,1]}; similarly for P_J⁺. By Lemma M2,
∫P_J⁻ ≥ (μ([−1,β]) − λ_n(β)) − (μ([−1,α)) + λ_n(α)) = μ(J) − ϵ_α − ϵ_β (μ has no atoms), and ∫P_J⁻ ≤ μ(J) since P_J⁻ ≤ 1_J; the
same for P_J⁺. Norms: P_β⁻ ≤ 1 with ∫P_β⁻ ≥ 0 gives ∫∣P_β⁻∣ = ∫(2(P_β⁻)⁺ − P_β⁻) ≤ 2; P_α⁺ ≥ 0 with ∫P_α⁺ ≤ 1 + λ_n(α) ≤ 2
(Lemma M4); so ∫∣P_J^±∣ ≤ 4. Finally ∣c_m∣ ≤ (m+1)∫∣P∣dμ and Σ_{m=1}^{2n}(m+1) = n(2n+3). □

**Lemma M4 (size of λ_n).** For n ≥ 1, N := 2n+3 and 0 < θ < π: λ_n(cos θ) ≤ 8sin²θ/N + 3π²/(4(n+1)³), and λ_n(t) ≤ 2/(n+1) for
every t ∈ [−1, 1].

*Proof.* If N sin θ ≥ 2, (1.1) and ∣sin Nθ∣ ≤ 1 give λ_n ≤ 4sin³θ/(N sin θ − 1) ≤ 8sin²θ/N; as s ↦ s³/(Ns − 1) increases for
s > 3/(2N), also λ_n ≤ 4/(N−1) = 2/(n+1). If N sin θ < 2, by the symmetry θ ↦ π − θ take θ ≤ π/2; then θ ≤ (π/2)sin θ < π/N, so
(k+1)θ < π/2 for k ≤ n and U_k(cos θ) ≥ (2/π)(k+1)θ/sin θ ≥ (2/π)(k+1); hence K_n ≥ (4/π²)(n+1)³/3 and λ_n ≤ 3π²/(4(n+1)³) ≤ 2/(n+1).
At t = ±1, (1.1) gives λ_n(±1) ≤ 2/(n+1). □

**Corollary M5 (degrees).** Let I = [a, b] ⊂ [0, π], L := b − a > 0, J := cos(I), n ≥ ⌈37.7/L⌉. Then ∫P_J⁻dμ ≥ μ(J)/3 > 0. Hence
  M\*(J) ≤ 2⌈37.7/L⌉ ≤ 75.4/L + 2 ≤ 75.4/ℓ + 2.

*Proof.* (i) sin is concave on [0, π], so sin ≥ its linear interpolant ℓ_{ab} ≥ 0 on [a, b] and
μ(J) = (2/π)∫_a^b sin² ≥ (2/π)∫_a^b ℓ_{ab}² = (2L/(3π))(sin²a + sin a sin b + sin²b) ≥ (2L/(3π))(sin²a + sin²b).
(ii) Among θ-intervals of length L the integral of the symmetric unimodal function sin² is smallest at an end, so
μ(J) ≥ μ(cos[0, L]) = (L − sin L cos L)/π ≥ (2/(3π))sin³L ≥ 16L³/(3π⁴) for L ≤ π/2 (both sides vanish at 0 and the derivatives are
2sin²L ≥ 2sin²L cos L), and μ(J) ≥ ½ for L ≥ π/2. (iii) N = 2n+3 ≥ 75.4/L > 24π/L. Lemma M4 at the endpoints inside (0, π) gives
ϵ_α + ϵ_β ≤ (8/N)(sin²a + sin²b) + 3π²/(2(n+1)³), where (8/N)(sin²a + sin²b) ≤ μ(J)/2 by (i), and 3π²/(2(n+1)³) ≤ μ(J)/50 by
(ii): for L ≤ π/2, n+1 > 37.7/L and (3π²/2)/37.7³ = 2.76·10⁻⁴ ≤ (1/50)·16/(3π⁴); for L > π/2, n+1 ≥ 14 and 3π²/(2·14³) ≤ 0.0027 ≤ 1/100.
So ∫P_J⁻ ≥ (1 − 0.52)μ(J) ≥ μ(J)/3 by Lemma M3. □

Sharpness and examples (`cms_table.out`). "CMS degree" is 2n for the least n with μ(J) > ϵ_α + ϵ_β (Lemma M3 then gives a
positive mean). "M\* ≥" is a certificate: if the closure of J contains no node of a positive quadrature on [−1, 1] that is exact for
degree d, then every P ≤ 1_J has ∫P = Σw_jP(y_j) ≤ 0, so M\*(J) ≥ d+1; the column uses the Gauss rules (nodes cos(jπ/(n+1)), exact
for 2n−1) and the rules of Lemma M1 with all nodes in [−1, 1].

| θ-interval | L | μ_ST | CMS degree | M\* ≥ | 𝒜(CMS degree) | known better |
|---|---|---|---|---|---|---|
| [0, π/2] (λ_f(p) > 0) | 1.571 | 0.5 | 8 | 2ᵃ | 11/3 | degree 2 ([RS] Cor D1): M\* = 2, A > 3 |
| [π/2, π] (λ_f(p) < 0) | 1.571 | 0.5 | 4 | 2ᵃ | 11/3 | degree 2: M\* = 2, A > 3 |
| [π/3, 2π/3] (∣x_p∣ < ½) | 1.047 | 0.609 | 8 | 4ᵇ | 11/3 | — |
| [π/4, π/2] | 0.785 | 0.409 | 12 | 3 | 5 | — |
| [0, π/4] | 0.785 | 0.0909 | 20 | 5 | 23/3 | — |
| [0.4, 0.6] | 0.2 | 0.0295 | 64 | 27 | 67/3 | — |
| [π/2 − 0.05, π/2 + 0.05] | 0.1 | 0.0636 | 124 | 60 | 127/3 | — |
| [0.05, 0.15] | 0.1 | 6.9·10⁻⁴ | 144 | 59 | 49 | — |
| [0, 0.1] | 0.1 | 2.1·10⁻⁴ | 182 | 60 | 185/3 | — |

ᵃ A degree-one P ≤ 1_J has P(0) ≤ 0 (0 ∉ J, or by continuity), and ∫P dμ = P(0) for linear P ([RS] Remark after Corollary D1).
ᵇ Symmetrizing P(x) ↦ (P(x) + P(−x))/2 keeps the mean, so degree ≤ 3 reduces to c₀ + c₂U₂ with U₂(±½) = 0, forcing c₀ = P(½) ≤ 0.

Asymptotically (L = 0.2 … 0.01): CMS degree × L → 4π for intervals centred at π/2, → 6π for intervals [0, L], ≈ 14.7 for
[L/2, 3L/2]; the certificate is ≈ 0.96·2π/L for centred intervals (60 at L = 0.1). In general, if I lies strictly between jπ/(n+1) and (j+1)π/(n+1), then
M\*(J) ≥ 2n; so for every L there are intervals of length L with M\* ≥ 2⌊π/L⌋ − 2, and A(J) ≍ 1/L is the true order of the method.
The uniform bound n ≤ ⌈37.7/L⌉ of Corollary M5 was checked on 287 intervals (7 lengths × 41 positions): all OK.

Remarks. (1) The CMS polynomials solve one-sided extremal problems on R; for intervals touching ±1 the extremal problem on [−1,1]
allows nodes at ±1 and lower degrees (degree 2 for the sign, against 8 above). A linear program gives M\*(J) exactly for any given J;
Theorem 1 is stated with M\*(J). (2) A Selberg-type trigonometric minorant in θ has Sato–Tate mean ≥ μ(J) − 4/(M+1), i.e. degree
≈ 2π/L for intervals centred at π/2, but it fails near the ends of [0, π], where μ(J) ≍ L³; not used.

## 2. Proposition T_m for every m ≥ 2

**Proposition T_m.** Let m ≥ 2, A > 𝒜(m), η₀ ∈ (0, 1]. There is an effective δ_m = δ_m(f, m, A, η₀) > 0 such that for every
δ ∈ (0, δ_m], with κ := δ/4, and X ≥ X₀(f, m, A, η₀, δ), there is E^{(m)}(X) ⊂ 𝔇(X) with #E^{(m)}(X) ≤ X^{1−δ/8} such that for every
D ∈ 𝔇(X) ∖ E^{(m)}(X), C ∈ G, x ≥ ∣D∣^A and H ∈ [x∣D∣^{−κ}, x]:  ∣h S_C^{(m)}∣ ≤ η₀ν_D H.
(m = 1, A > 9/4: [RS] Proposition T.) The same holds for D > 0 with narrow classes (§6).

[RS] proved this for m = 2 (A > 3) and sketched route F for m ≥ 3 in its §7. Here the routes are written for general m with every
m-dependence displayed. The statement for A implies it for every A′ > A (same δ, κ, exceptional set), so fix a route with threshold
A^♭ < A ≤ A^♭ + ¼: **route F**, A^♭ := (m+7)/3 (inputs as in [RS]); **route H′**, A^♭ := max{11/3, (m+3)/3} (Proposition F′ with
PT Corollary 1.2), used when it is the smaller one (m ≥ 5); **route H**, A^♭ := max{5, (m+3)/3} (Proposition F′ without PT, so
with published inputs only), smaller than (m+7)/3 for m ≥ 9. 𝒜(m) is the minimum of the thresholds of F and H′.

### 2.1 Parameters (all checked in `params.py` for 2 ≤ m ≤ 80, every route)
- σ₀ := 1 − 1/(2A^♭); s₁ := A(1−σ₀) − ½ = (A − A^♭)/(2A^♭) > 0.
- Route F: ε′ := s₁/4. With φ_F(σ) := ((m+10)/(3−2σ) − A)(1−σ): φ_F(σ₀) = 1 − s₁ = 1 − 4ε′, because (m+10)/(1 + 1/A^♭) = 3A^♭;
  φ_F′(σ) = A − (m+10)/(3−2σ)² is decreasing and φ_F′(σ₀) = A − 3(A^♭)²/(A^♭+1) < 0 (A^♭ ≥ 3, A ≤ A^♭ + ¼); so φ_F ≤ 1 − 4ε′ on
  [σ₀, 1]. Also (2m+10)e(σ) ≤ (2m+10)e(σ₀) = 3(m+5)/(m+10) ≤ 3 on [σ₀, 1].
- Routes H′ and H: σ₂ := 1 − 1/(4(m+10)), φ_H(σ) := ((m+6)/(3−2σ) − A)(1−σ), and the "mid-range" exponents of Proposition F′:
  g(σ) := max{1 + 3(1−σ), ¾ + (11/2)(1−σ)} (H′), g(σ) := 1 + 5(1−σ) (H). Put
  ε′ := min{s₁/4, (A−3)/(16(m+10)), (2A + A^♭ − 11)/(16A^♭), (1 − φ_H(σ₀))/4, 3/16} (H′), resp.
  ε′ := min{s₁/4, (A−5)/(16(m+10)), (1 − φ_H(σ₀))/4, 3/16} (H). On [σ₀, σ₂], g(σ) − A(1−σ) ≤ 1 − 4ε′: it is linear in 1−σ on each
  piece, ≤ 1 − (A−3)(1−σ₂) resp. 1 − (A−5)(1−σ₂) at σ₂, and at σ₀ (1−σ₀ = 1/(2A^♭)) the H′ term ¾ + (11/2 − A)/(2A^♭) equals
  1 − (2A + A^♭ − 11)/(4A^♭) with 2A + A^♭ − 11 > 3A^♭ − 11 ≥ 0. Further φ_H ≤ φ_H(σ₀) ≤ 1 − 4ε′ on [σ₀, 1], since φ_H′ is
  decreasing with φ_H′(σ₀) = A − (m+6)(A^♭)²/(A^♭+1)² < 0 (if A^♭ = (m+3)/3 then m+6 = 3(A^♭+1) and the subtracted term is
  3(A^♭)²/(A^♭+1) > A^♭ + ¼; if A^♭ = 11/3, m ≥ 5, it is ≥ 11·121/196 > 3.92; if A^♭ = 5, m ≥ 9, it is ≥ 15·25/36); and
  1 − φ_H(σ₀) > 0 (= s₁ if A^♭ = (m+3)/3, and a direct computation otherwise). On [σ₂, 1]: (m+12)e(σ) ≤ (m+12)(1−σ) ≤ ⅓
  and (2m+2)e(σ) ≤ ½. On [σ₀, σ₂]: (2m+2)e(σ) ≤ (m+1)/(A^♭+1) ≤ 3.
- Common: c_m := 10⁷(m+1)⁴ and B_m from (J4) of [RS]; c₁′ := c/(2m²), c the absolute constant in Thorner's zero-free region ((J5) of
  [RS]: L(σ+it, Sym^m f ⊗ χ) ≠ 0 for σ ≥ 1 − c/(m²[m^y log(mC(f)) + log(C(χ)(3+∣t∣))]), uniform in χ); σ_g := 1 − 1/(2c_m);
  δ_m := min{1/2840, s₁/10, ε′/10, A/(8c_m), δ_e}, δ_e := A/(1.3(c_m + c₁′^{−1}log(64B_m/η₀))).
  For δ ≤ δ_m: κ := δ/4, ε₁ := δ/100, ε₀ := δ/(100(m+1)), δ′ := δ/100 (Proposition F′), σ_j := σ₀ + jε′/A (0 ≤ j ≤ J, J minimal
  with σ_J ≥ 1 − ε′/A), T_X := 2X^{κ+ε₁}, Q₀(X) := (X/2)^δ; for a given D: Q₀ := ∣D∣^δ, T₁ := 2∣D∣^{κ+ε₁} ≤ T_X, 𝔉_s := {d ∈ 𝔉(D) :
  ∣d∣ ≤ Q₀}, 𝔉_l := 𝔉(D) ∖ 𝔉_s.

### 2.2 Exceptional sets (X a power of 2)
- E₁^{(m)}(X) := {D ∈ 𝔇(X) : N_D^{c,m}(σ_j, T_X) > X^{A(1−σ_j)−2ε′} for some j}. By Markov, #E₁ ≤ Σ_jX^{2ε′−A(1−σ_j)}Σ_DN_D^{c,m}(σ_j, T_X).
  Route F: [RS] Proposition F (ε = ε₀) bounds the inner sum by C(X^{m+10}T_X^{2m+10})^{e(σ_j)+ε₀}, so the X-exponent of the j-th
  term is ≤ 2ε′ + φ_F(σ_j) + (m+10)ε₀ + (2m+10)(κ+ε₁)(e(σ_j)+ε₀) + o(1) ≤ 1 − 2ε′ + 0.11δ + 3.1(κ+ε₁) ≤ 1 − ε′.
  Routes H, H′, σ_j ≥ σ₂: Proposition F″ (§3) on the ≤ 3T_X windows bounds the inner sum by CT_XX^{δ′}(X^{m+12}(2T_X)^{2m+2})^{e(σ_j)};
  the exponent is ≤ 2ε′ + 2(κ+ε₁) + δ′ + ⅓ ≤ 1 − ε′. (So routes H, H′ never use a Rankin–Selberg large sieve with π′ = Sym^m f.)
  Routes H, H′, σ_j < σ₂: Proposition F′ (§3) on each of the ≤ 3T_X unit windows ∣γ − k∣ ≤ 1 (∣k∣ ≤ T_X, ∣k∣+3 ≤ 2T_X) bounds
  the inner sum by CT_XX^{δ′}[X^{g(σ_j)} + (X^{m+6}(2T_X)^{2m+2})^{e(σ_j)}]; the exponent is
  ≤ 2ε′ + (κ+ε₁) + δ′ + max{g(σ_j) − A(1−σ_j), φ_H(σ_j) + 3(κ+ε₁)} ≤ 1 − 2ε′ + 4(κ+ε₁) + δ′ ≤ 1 − ε′.
  Hence #E₁ ≤ CJX^{1−ε′} (δ ≤ ε′/10).
- E₂^{(m)}(X): D such that some d ∈ 𝔉(D) with ∣d∣ > Q₀(X) has a zero of L(s, Π_m ⊗ χ_d) with β > σ_g, ∣γ∣ ≤ T_X. As [RS] §5.2
  ((J4) with ν^♯ = 1 at σ_g; each d with ∣d∣ ∈ (Z, 2Z] divides ≤ 2X/Z elements of 𝔇(X)): #E₂ ≤ C(m)XQ₀(X)^{−1/2}T_X^{1/2}log X ≤ X^{1−δ/4}.
- E^{(m)}(X) := E₁ ∪ E₂, #E^{(m)}(X) ≤ X^{1−δ/8} for X ≥ X₀.

### 2.3 Decomposition
Verbatim [RS] §5.3: for D ∈ 𝔇(X) ∖ E^{(m)}(X), ∣D∣ ≥ D₀, C, x, H: hS_C^{(m)} = Σ_χχ̄(C)S_χ^{(m)}, the genus part is
Σ_{d∈𝔉(D)}ε_d(C)ψ_w(Π_m ⊗ χ_d), (J6) is applied to every term (no main terms, m ≥ 1), and E_a, …, E_f bound: (a) zeros of L_m(·,χ),
χ ∈ Ĝ_c, β ≥ σ₀, ∣γ∣ ≤ T₁; (b) zeros of L_m(·,χ), all χ ∈ Ĝ, β < σ₀, ∣γ∣ ≤ T₁; (c) all zeros with ∣γ∣ > T₁; (d), (e) zeros of
L(s, Π_m ⊗ χ_d) for d ∈ 𝔉_l, resp. d ∈ 𝔉_s, ∣γ∣ ≤ T₁; (f) remainders.

### 2.4 Bounds (∣D∣ ≥ D₀(f, m, A, η₀, δ); x = ∣D∣^A y, y ≥ 1)
- (a) D ∉ E₁ and x^{σ−1} ≤ 2^AX^{−A(1−σ)}: E_a ≤ HΣ_{j<J}N_D^{c,m}(σ_j,T_X)x^{σ_{j+1}−1} ≤ J2^AHX^{−ε′}; N_D^{c,m}(σ_J,T_X) < 1.
- (b) Only the trivial count (Theorem C is not needed for m ≥ 2). In a window ∣γ − k∣ ≤ 1, N(σ) := N_D^{(m)}(σ;k) ≤ N(0) ≤
  C(m)h log(∣D∣(∣k∣+3)) ≤ ∣D∣^{1/2+o(1)}(∣k∣+3)^{o(1)} ((R4), [C] (F3)), and Σ_{β<σ₀}x^{β−1} ≤ x^{−1/2}N(0) + log x∫_{1/2}^{σ₀}x^{σ−1}N(σ)dσ
  with x^{σ−1}log x ≤ ∣D∣^{−A(1−σ)}(A log∣D∣ + 2A^♭) for σ ≤ σ₀ (y^{−u}log y ≤ 1/u, u = 1−σ ≥ 1/(2A^♭)). Each window is
  ≤ ∣D∣^{1/2−A(1−σ₀)+o(1)} = ∣D∣^{−s₁+o(1)}; over ≤ 3T₁ windows, E_b ≤ H∣D∣^{−s₁+κ+ε₁+o(1)}.
- (c) Each window holds ≤ C(m)h log(∣D∣(∣k∣+3)) zeros; ∣W(ρ)∣ ≤ C_kHx^{β−1}(2x/(H∣ρ∣))^k with k ≥ 2/ε₁ and x/H ≤ ∣D∣^κ: E_c ≤ H∣D∣^{−1}.
- (d) D ∉ E₂: these zeros have β ≤ σ_g and number ≤ 2^{ω(D)}C(m)T₁log(∣D∣T₁): E_d ≤ H∣D∣^{−A/(2c_m)+κ+ε₁+o(1)}.
- (e) As [RS] §5.4(e): L := log(3C(Π_m)Q₀T₁) ≤ 1.3δ log∣D∣ (∣D∣ ≥ D₀(m)), ℓ := log x ≥ A log∣D∣ ≥ 2c_mL (δ ≤ A/(8c_m)); no zero with
  β > 1 − c₁′/L ((J5), Thorner, for every m); (J4) with the Landau–Siegel repulsion ν^♯ ≤ 2ν_D; with δ ≤ δ_e this gives
  E_e ≤ η₀ν_DH/8.
- (f) Σ∣R∣ ≤ C(m)h(H/x + x^{−1/2})log(∣D∣x) ≤ H∣D∣^{−1} (x ≥ ∣D∣²).

### 2.5 Conclusion
ε′ > δ/2, s₁ − κ − ε₁ > δ/2 and A/(2c_m) − κ − ε₁ > δ/2 (δ ≤ δ_m) give E_a + E_b + E_c + E_d + E_f ≤ H∣D∣^{−δ/2−η₁} for some
η₁ > 0, which is ≤ η₀ν_DH/8 for ∣D∣ ≥ D₀, since ν_D ≥ c′∣D∣^{−δ/2}/log∣D∣ ([R] (I5)). With (e), ∣hS_C^{(m)}∣ ≤ η₀ν_DH/4. □

Exceptional zeros are handled as in [RS] §5.6: those of the non-genus twists (χ ∈ Ĝ_c) by discarding E₁ (for D ∉ E₁ there is no
zero with β ≥ σ_J, ∣γ∣ ≤ T_X), those of the genus twists L(s, Sym^m f ⊗ χ_d) by Thorner's theorem (no Landau–Siegel zero, every m) and
the repulsion of (J4).

## 3. Propositions F′ and F″: per-D mollifiers, three large sieves

**Proposition F′.** Let m ≥ 1, 0 < δ′ ≤ 1/100, X ≥ X₀(f, m, δ′), ½ + 3δ′ ≤ σ ≤ 1, t ∈ R, T := ∣t∣ + 3 ≤ X. Then
  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ_c(D)}N_m(σ; t, χ) ≤ C(δ′, m)X^{δ′}[X^{1+3(1−σ)} + X^{3/4+(11/2)(1−σ)} + (X^{m+6}T^{2m+2})^{(1−σ)/(3−2σ)}],
and, without PT Corollary 1.2 (published inputs only), the same with X^{1+5(1−σ)} in place of the first two terms.

Compare [RS] Proposition F, (X^{m+10}T^{2m+10})^{e(σ)}. There, the Humphries–Thorner detector uses one mollifier for the whole family,
of length equal to the family large-sieve threshold Q²∣𝒮∣ ≍ X^{7/2}; here every D gets its own mollifier of length ∣D∣^{1/2}, so class
II costs the family size X^{3/2} instead of X^{7/2}, which replaces m+10 by m+6. The price is in class I, whose sums of length N
run far below the family threshold X^{7/2}: the per-D sieve costs X·N, PT's Hölder form (k = 2) costs N·X^{3/4} + X^{5/2}, and
HT costs N + X^{7/2}; the worst lengths are N = X^{3/2} (per-D meets PT, giving X^{1+3(1−σ)}) and N = X^{11/4} (PT meets HT,
X^{3/4+(11/2)(1−σ)}); without PT, N = X^{5/2} (per-D meets HT, X^{1+5(1−σ)}). The family sieves need class I coefficients
independent of D, which the squarefree detector of [R3] §2 provides.

*Setup.* δ″ := δ′/(4m+20); P_m := 64(m+1)². n is *admissible* if it is squarefree and every prime p ∣ n has p > P_m and p ∤ N.
a_V(n) := Σ_{k∣n, k≤V}μ(k). For squarefree n, λ_m(n) := Π_{p∣n}U_m(x_p) (= λ_{Sym^m f}(n)). For D ∈ 𝔇(X), χ ∈ Ĝ(D):
λ_χ(n) := Σ_{N𝔞=n}χ(𝔞) (multiplicative; at p: χ(𝔭) + χ(𝔭̄) split, 0 inert, χ(𝔭) ramified) and a(n) := λ_m(n)λ_χ(n).
(F0) For admissible p, a(p) is the coefficient of p^{−s} in the local factor L_{m,p}(s,χ) := Π_{𝔭∣p}Π_i(1 − α_{m,i}(𝔭)χ(𝔭)N𝔭^{−s})^{−1}
of L_m(s,χ): by (R1) and (R1′), at split p the two factors have parameters (e^{i(m−2j)θ_p})_j times χ(𝔭), resp. χ(𝔭̄); at inert p,
N𝔭 = p²; at ramified p, λ_m(𝔭) = U_m(x_p). Also ∣a(p)∣ ≤ 2(m+1), and L_{m,p}^{−1} is a polynomial in p^{−s} of degree ≤ 2m+2 whose
k-th coefficient is ≤ binom(2m+2, k) in modulus. (So a(n) is the n-th coefficient of L_m for admissible n, and only the factor λ_χ(n)
depends on D.)

**Lemma F1.** Let 𝒢(s,χ) := Σ_{n adm}a_V(n)a(n)n^{−s}.
(a) For Re s > 1, 𝒢 = L_m(s,χ)·E·H·M̃ with E := Π_{p≤P_m or p∣N}L_{m,p}(s,χ)^{−1}, H := Π_{p adm}(1 + a(p)p^{−s})L_{m,p}(s,χ)^{−1},
M̃ := Σ_{k≤V adm}μ(k)k^{−s}Π_{p∣k}ψ_p(a(p)), ψ_p(u) := u/(1 + up^{−s}). E is entire and ∣E∣ ≤ C(f, m) on Re s ≥ 0; H and M̃ are
holomorphic on Re s > ½ and ∣H∣ ≤ C(δ″, m) on Re s ≥ ½ + δ″. So 𝒢 is holomorphic on Re s > ½, of polynomial growth in vertical
strips there, and vanishes at every zero of L_m(·,χ) with real part > ½.
(b) The coefficient of 𝒢 at n = 1 is 1, and at every 1 < n ≤ V it is 0.
(c) On Re s = ½ + δ″, for V ≤ ∣D∣^{1/2}: Σ_{χ∈Ĝ}∣M̃(s,χ)∣² ≤ C(δ″, m)∣D∣^{1/2}(log 3∣D∣)^{14(m+1)²+2}.
(d) For v ∈ R: Σ_{χ∈Ĝ}∣L_m(½+δ″+iv, χ)∣² ≤ C(ε, m)∣D∣^{(m+1)/2+ε}(∣v∣+3)^{m+1+ε}.

*Proof.* (a) Since a is multiplicative on coprime squarefree arguments, Σ_{n adm}a_V(n)a(n)n^{−s} = Σ_{k≤V adm}μ(k)a(k)k^{−s}
Π_{p adm, p∤k}(1 + a(p)p^{−s}) = M̃·Π_{p adm}(1 + a(p)p^{−s}), and Π_{p adm}(1 + a(p)p^{−s}) = L_m·E·H. For p > P_m and Re s ≥ ½,
∣a(p)p^{−s}∣ ≤ 2(m+1)p^{−1/2} < ¼, so ψ_p is holomorphic there. By (F0), (1 + a(p)p^{−s})L_{m,p}^{−1} = 1 + Σ_{k≥2}c_kp^{−ks} with
∣c_k∣ ≤ binom(2m+2,k) + 2(m+1)binom(2m+2,k−1) ≤ 8^{m+2}; on Re s ≥ ½ + δ″ the sum is ≤ 2·8^{m+2}p^{−1−2δ″}, so H converges
absolutely and ∣H∣ ≤ exp(8^{m+3}ζ(1+2δ″)). E is a finite product of polynomials in p^{−s} bounded by 4^{m+1} on Re s ≥ 0. Growth: (R2).
(b) For admissible n ≤ V every divisor is ≤ V, so a_V(n) = Σ_{k∣n}μ(k).
(c) As [R3] Lemma 3(c). With θ_p(u) := up^{−s}/(1 + up^{−s}), ψ_p(u) = u(1 − θ_p(u)) and ∣θ_p(a(p))∣ ≤ (8/3)(m+1)p^{−1/2}. Expanding
Π_{p∣k}(1 − θ_p) and writing k = tk′: M̃ = Σ_{t≤V adm}Θ_t(χ)M_t(χ), Θ_t := μ(t)(−1)^{ω(t)}a(t)t^{−s}Π_{p∣t}θ_p(a(p)),
M_t(χ) := Σ_{k′≤V/t adm, (k′,t)=1}μ(k′)a(k′)k′^{−s}, and ∣Θ_t∣ ≤ (6(m+1)²)^{ω(t)}t^{−1}. Cauchy–Schwarz over t gives
Σ_χ∣M̃∣² ≤ (Σ_{t≤V}(6(m+1)²)^{ω(t)}t^{−1})·Σ_t(6(m+1)²)^{ω(t)}t^{−1}Σ_χ∣M_t∣², and the first factor is ≤ C(log 3V)^{6(m+1)²}. Now
M_t(χ) = Σ_𝔞c_𝔞χ(𝔞) with c_𝔞 := μ(N𝔞)λ_m(N𝔞)N𝔞^{−s} on the (primitive) ideals of admissible norm ≤ V/t prime to t; [C] Lemma 2 gives
Σ_χ∣M_t∣² ≤ C log²(3∣D∣)(V/t + ∣D∣^{1/2})Σ_{k′≤V}d(k′)(m+1)^{2ω(k′)}k′^{−1} ≤ C(m)(log 3∣D∣)^{2(m+1)²+2}∣D∣^{1/2}; with the first factor squared
(Σ_tw_t ≤ C(log 3V)^{6(m+1)²}) this gives (c).
(d) The proof of [RS] Lemma B (V = 1) on the line Re s = ½ + δ″: in the approximate functional equation L_m = S₁ + ε_m(χ,s)S₂ the dual
sum S₂ has coefficients of size ∣a∣n^{−1/2+δ″}, so its large-sieve bound gains a factor P^{2δ″} (P ≤ L₀𝔏^ε, L₀ := √q_m(∣v∣+3)^{m+1}),
while ∣ε_m(χ,s)∣ ≍ (q_m·Π_j(∣v + Im μ_j∣ + 3))^{−δ″} ≪ L₀^{−2δ″}(∣D∣(∣v∣+3))^{O(m)ε} by Stirling; S₁ is as in Lemma B. □

**Lemma F2 (detector).** Let ∣D∣ ≥ D₀(δ″, m), ½ + 2δ″ ≤ σ ≤ 1, Y ≥ V ≥ 2, Z := 2Y log²Y, χ ∈ Ĝ(D), and ρ = β + iγ a zero of
L_m(·,χ) with σ ≤ β ≤ 1. Then (I) ∣Σ_jF_j(χ,ρ)∣ ≥ ¼, where F_j(χ,ρ) := Σ_{n adm}a_V(n)a(n)φ(n/2^j)n^{−ρ}e^{−n/Y} and j runs over
2^{j+1} > V, 2^{j−1} < Z (φ the dyadic partition of unity of [C] Theorem 1), or (II) ∫∣𝒢(½+δ″+i(γ+u), χ)∣e^{−∣u∣}du ≥ c₁Y^{β−1/2−δ″},
c₁ = c₁(δ″) > 0.

*Proof.* [R3] Lemma 4 without its pole-killing factor (L_m has no pole for m ≥ 1, also for χ = 1). (2πi)^{−1}∫_{(2)}𝒢(ρ+w)Γ(w)Y^wdw
= Σ_n c_nn^{−ρ}e^{−n/Y} = e^{−1/Y} + Σ_jF_j + (n > Z), by Lemma F1(b); ∣c_n∣ ≤ d(n)(2(m+1))^{ω(n)} ≤ C(m)n^{1/4}, so the tail is ≤ 1/16
for ∣D∣ ≥ D₀. Move the line to Re w = ½ + δ″ − β ∈ [−½ + δ″, −δ″]: only w = 0 is crossed, with residue 𝒢(ρ) = 0 (Lemma F1(a)), and
∣Γ(½+δ″−β+iu)∣ ≤ Cδ″^{−1}e^{−∣u∣}. If (I) fails, the new integral has modulus ≥ ½, which is (II) with c₁ := πδ″/C. □

*Class I coefficients.* Removing β as in [C] Theorem 1 / [R3] §2 (Φ_β(v) := φ(v)v^{σ−β}, ∣Φ̂_β(τ)∣ ≤ C_A(1+∣τ∣)^{−A} uniformly,
φ₂ ∈ C_c^∞((1/4,4)), φ₂ = 1 on supp φ): ∣F_j(χ,ρ)∣ ≤ (2π)^{−1}∫∣Φ̂_β(τ)∣∣P_j(χ, γ+τ)∣dτ with
  P_j(χ,t′) := Σ_nb_{j,t′}(n)λ_χ(n),  b_{j,t′}(n) := 1_{adm}(n)a_V(n)λ_m(n)n^{−σ−it′}e^{−n/Y}φ₂(n/2^j).
**b_{j,t′} depends on neither D nor χ.** ∣b_{j,t′}(n)∣ ≤ d(n)(m+1)^{ω(n)}n^{−σ}, ∂_{t′}b_{j,t′} = −i(log n)b_{j,t′}, supports in (2^j/4, 2^{j+2}).

**Lemma F3 (three large sieves).** Let N := 2^j, b := b_{j,t′} or ∂_{t′}b_{j,t′}, ‖b‖_♯² := Σ_nτ(n)³∣b(n)∣², ‖b‖_∞ := max_n∣b(n)∣.
(LS1) For each D ∈ 𝔇(X): Σ_{χ∈Ĝ(D)}∣Σ_nb(n)λ_χ(n)∣² ≤ C log²(3∣D∣)(4N + ∣D∣^{1/2})‖b‖_♯².
(LS2) Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ_c(D)}∣Σ_nb(n)λ_χ(n)∣² ≤ C(ε)X^ε(N + X^{7/2})‖b‖_♯².
(LS3) Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ_c(D)}∣Σ_nb(n)λ_χ(n)∣² ≤ C(ε)X^ε(N X^{3/4} + X^{5/2})·N‖b‖_∞².

*Proof.* (LS1): Σ_nb(n)λ_χ(n) = Σ_𝔞b(N𝔞)χ(𝔞); ideals of squarefree norm are primitive and there are ≤ τ(n) of norm n; [C] Lemma 2.
(LS2): for squarefree n write n = n₁n₂ with n₁ := (n, D); then λ_χ(n) = λ_χ(n₁)λ_χ(n₂) with ∣λ_χ(n₁)∣ = 1, and by Cauchy–Schwarz over
the τ(D) ≤ X^ε divisors n₁ ∣ D, ∣Σ_nb(n)λ_χ(n)∣² ≤ X^εΣ_{n₁∣D}∣Σ_{(n₂,D)=1}b(n₁n₂)λ_χ(n₂)∣². For fixed n₁ apply HT Theorem 4.1 with
n = 2, F = Q, π′ = 𝟙, T = 1, θ₂ = 0 (every π_χ satisfies GRC) to 𝒮_{n₁} := {π_χ : D ∈ 𝔇(X), n₁ ∣ D, χ ∈ Ĝ_c(D)}: these are cuspidal
with unitary central character χ_D and analytic conductor ≤ 12X =: Q, each π_χ = π_χ̄ occurs for two characters, and HT's sum over
N𝔫 ∈ (y, ey] with gcd(𝔫, 𝔮_π) = 1 is exactly our sum over n₂ (the support of b(n₁·) is covered by three such intervals). HT gives
Σ_{π∈𝒮_{n₁}}∣Σ_{(n₂,D)=1}b(n₁n₂)λ_π(n₂)∣² ≤ C(ε)X^ε(4N/n₁ + Q²∣𝒮_{n₁}∣)Σ_{n₂}∣b(n₁n₂)∣², with ∣𝒮_{n₁}∣ ≤ Σ_{∣D∣≤X}h(D) ≤ CX^{3/2}log X.
Sum over n₁ (each n is n₁n₂ for ≤ τ(n) choices). (PT Theorem 1.1 gives (LS2) directly, without the coprimality condition.)
(LS3): PT Corollary 1.2 with n = 2, π₀ = 𝟙, k = 2, applied to 𝒮 := {π_χ : D ∈ 𝔇(X), χ ∈ Ĝ_c(D)} (∣𝒮∣ ≤ CX^{3/2}log X, Q ≤ 12X):
Σ_{π∈𝒮}∣Σ_{n≤4N}λ_π(n)a(n)∣² ≪ (NQ)^ε(N∣𝒮∣^{1/2} + Q∣𝒮∣)·N‖a‖_∞², all coefficients, no coprimality condition and no Ramanujan
hypothesis; λ_π(n) = λ_χ(n) for π = π_χ (Hecke eigenvalues of θ_χ), and each π occurs for two χ. Its implied constant depends on
k = 2 and n = 2 only. □

*Proof of Proposition F′.* Let 𝒳(D) be the set of χ ∈ Ĝ_c(D) for which L_m(·,χ) has a zero in [σ,1] × [t−1, t+1]; fix one, ρ_χ, for
each. By (R4), Σ_{χ∈Ĝ_c}N_m(σ;t,χ) ≤ C(m)log(XT)∣𝒳(D)∣. Take V := X^{1/2}/2 ≤ ∣D∣^{1/2}, Y := (X^{m+6}T^{2m+2})^{1/(2(3−2σ))}
(≥ X^{7/4} ≥ V), Z := 2Y log²Y. Lemma F2 splits 𝒳(D) = 𝒳_I ∪ 𝒳_II.
*Class II.* Sum (II) over 𝒳_II(D), substitute v = γ_χ + u (e^{−∣u∣} ≤ e^{1−∣v−t∣}) and apply Cauchy–Schwarz over all of Ĝ at each v,
with Lemma F1(a), (c), (d):
  c₁Y^{σ−1/2−δ″}∣𝒳_II(D)∣ ≤ C(δ″,m)∫e^{−∣v−t∣}(Σ_χ∣L_m∣²)^{1/2}(Σ_χ∣M̃∣²)^{1/2}dv ≤ C(ε,δ″,m)∣D∣^{(m+2)/4+ε}T^{(m+1)/2+ε}.
Summing over D: Σ_D∣𝒳_II(D)∣ ≤ CX^{(m+6)/4+ε}T^{(m+1)/2+ε}Y^{1/2+δ″−σ} = CX^εT^εY^{δ″}(X^{m+6}T^{2m+2})^{e(σ)}, because
(½ − σ)/(2(3−2σ)) + ¼ = e(σ); and Y^{δ″} ≤ (X^{3m+8})^{δ″/2} ≤ X^{δ′/2} (T ≤ X).
*Class I.* As in [C] Theorem 1 and [R3] Theorem 2: each χ ∈ 𝒳_I(D) has ∣F_j(χ,ρ_χ)∣ ≥ 1/(4J) for some j (J ≤ 2 log Z blocks), and by
Cauchy–Schwarz in τ and Sobolev on I_τ := [t+τ−1, t+τ+1],
  Σ_D∣𝒳_I(D)∣ ≤ CJ²Σ_j∫(1+∣τ∣)^{−A}∫_{I_τ}Σ_DΣ_{χ∈Ĝ_c}(∣P_j(χ,t′)∣² + ∣∂_{t′}P_j(χ,t′)∣²)dt′dτ.
Here ‖b‖_♯² ≤ C(ε,m)N^{1−2σ+ε}log²(2N) and N‖b‖_∞² ≤ C(ε,m)N^{1−2σ+ε} for both sequences (∣b(n)∣ ≤ τ(n)(m+1)^{ω(n)}n^{−σ}),
and N = 2^j ∈ (V/2, 2Z); write η := 1−σ, so N^{1−2σ} = N^{2η−1}.
- N ≤ X^{3/2}, (LS1) summed over the ≤ X discriminants: ≤ CX^{1+ε}(4N + X^{1/2})N^{2η−1} ≤ 12X^{1+ε}N^{2η} ≤ 12X^{1+3η+ε}.
- X^{3/2} < N ≤ X^{11/4}, (LS3): ≤ CX^ε(X^{3/4}N^{2η} + X^{5/2}N^{2η−1}) ≤ CX^ε(X^{3/4+(11/2)η} + X^{1+3η}).
- X^{11/4} < N ≤ X^{7/2}, (LS2): ≤ CX^ε(N^{2η} + X^{7/2}N^{2η−1}) ≤ CX^ε(X^{7η} + X^{7/2+(11/4)(2η−1)}) ≤ 2CX^{3/4+(11/2)η+ε}.
- N > X^{7/2}, (LS2): ≤ CX^ε(N^{2η} + X^{7/2}N^{2η−1}) ≤ 2CX^εZ^{2η}.
(Without PT: (LS1) for N ≤ X^{5/2}, (LS2) above, giving X^{1+5η} in place of the first three lines.) Finally
Z^{2−2σ} ≤ C(log Y)⁴(X^{m+6}T^{2m+2})^{e(σ)}. Adding the classes, with ε := δ′/10, proves the proposition. □

**Proposition F″ (near σ = 1).** For m ≥ 2 and the same ranges, Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ_c(D)}N_m(σ;t,χ) ≤ C(δ′,m)X^{δ′}(X^{m+12}T^{2m+2})^{e(σ)}.
*Proof.* The proof of Proposition F′ with V := X^{7/2}, Y := (X^{m+12}T^{2m+2})^{1/(2(3−2σ))} (≥ X^{(m+12)/4} ≥ V for m ≥ 2). Lemma F1(c)
holds for every V ≥ 1 in the form Σ_χ∣M̃∣² ≤ C(δ″,m)(V + ∣D∣^{1/2})(log 3V∣D∣)^{14(m+1)²+2} (same proof), so class II is
≤ CX^εY^{1/2+δ″−σ}X^{1+(m+1)/4}V^{1/2}T^{(m+1)/2} = CX^εY^{δ″}(X^{m+12}T^{2m+2})^{e(σ)}; every class I block has N > V/2 and (LS2) gives
≤ CX^ε(N + X^{7/2})N^{1−2σ} ≤ 3CX^εZ^{2−2σ}. □ This is [RS] Proposition F with exponent m+12 instead of m+10, but its only family large
sieve is the GL₂ one (π′ = 𝟙), whose constant does not depend on m.

*D > 0.* The same statements hold for Σ_{D∈𝔇⁺(X)}Σ_{χ∈Ĝ⁺_c(D)}, with the dictionary of [RS] §6: (LS1) by [D2] Lemma 5 (𝔔 = 9D,
M = 0); (LS2), (LS3) by HT and PT for 𝒮⁺ := {π_χ}, Maass forms of eigenvalue ¼ satisfying GRC, C(π_χ) ≤ 12D, ∣𝒮⁺∣ ≤ CX^{3/2}
([R2] Prop 2);
Lemma F1(c) with [D2] Lemma 5; Lemma F1(d) by Lemma B for each archimedean type a ∈ {0, 1}; (R4) with h⁺ ≤ 5D^{1/2}log²(3D).

## 4. Theorem 1: every interval

**Theorem 1.** Let J ⊂ [−1,1] be an interval with nonempty interior, M := M\*(J), and A > 𝒜(M). There are effective c, κ > 0 and
X₀ (depending on f, J, A) such that for X ≥ X₀ all fundamental discriminants D with ∣D∣ ≤ X, negative or positive, outside a set
of size ≤ X^{1−c}, satisfy: for every ideal class C of Q(√D) (narrow class if D > 0), every x ≥ ∣D∣^A and every H ∈ [x∣D∣^{−κ}, x],
  Σ_{𝔭∈C of degree one, N𝔭 = p ∤ ND, x < p ≤ x+H, cos θ_p ∈ J} log p ≥ c₀ν_D H/(32h^{(+)}(D)),
where c₀ > 0 is the Sato–Tate mean of a degree-M minorant of 1_J and ν_D ≥ c′∣D∣^{−δ/2}/log∣D∣ is as in [R] Theorem 2 (ν_D = 1 unless a
genus factor χ_d, ∣d∣ ≤ ∣D∣^δ, has a Landau–Siegel zero). For a θ-interval I of length L and J = cos(I), Corollary M5 gives
M\*(J) ≤ 2⌈37.7/L⌉, so every A > max{11/3, 25.2/L + 1.7} is admissible, with c₀ ≥ μ_ST(J)/3 for the CMS minorant (with
published inputs only: A > max{5, 25.2/L + 1.7}).

*Proof.* If M = 0 then J ⊃ (−1,1) and the claim is [R2] Theorems 2′ and 4. Let M ≥ 1 and P = Σ_{m≤M}c_mU_m with P ≤ 1_J on (−1,1),
c₀ = ∫P dμ > 0. For each 1 ≤ m ≤ M with c_m ≠ 0, Proposition T_m holds ([RS] Proposition T for m = 1, §2 for m ≥ 2), since
A > 𝒜(M) ≥ 𝒜(m); use η₀ := c₀/(16Σ_{m≥1}∣c_m∣). For m = 0, [R2] Theorem 2′ (D < 0) and Theorem 4 (D > 0) give hS_C^{(0)} ≥ ν_DH/8
outside their exceptional sets for every δ ≤ 1/2840 with κ = δ/4 ([RS], proof of Theorem D). Take δ := min{1/2840, δ_1, …, δ_M},
κ := δ/4, and remove the union over dyadic X′ ≤ X of all these exceptional sets (≤ X^{1−δ/9}). The proof of [RS] Theorem D then applies
verbatim: with P_𝔞 := Σ_mc_mb_m(𝔞), hΣ_{𝔞∈C}Λ_K(𝔞)P_𝔞w(N𝔞) ≥ c₀ν_DH/8 − Σ_{m≥1}∣c_m∣η₀ν_DH = c₀ν_DH/16; P_𝔭 = P(cos θ_p) for degree-one
𝔭 with p ∤ ND; the other 𝔞 have ∣P_𝔞∣ ≤ Σ_m(m+1)∣c_m∣ and number ≤ C(H/√x + ω(ND) + log x)log x, so they contribute ≤ c₀ν_DH/(32h)
for ∣D∣ ≥ D₀ (x ≥ ∣D∣², h ≤ ∣D∣^{1/2+o(1)}); and P(x_p) ≤ 1_J(x_p) since ∣x_p∣ < 1 for p > x ≥ ∣D∣^A (§0), 0 ≤ w ≤ 1_{(x,x+H]}. □

**Corollary 1.1.** Let D be outside the exceptional set of Theorem 1. Every primitive binary quadratic form of discriminant D (positive
definite if D < 0, indefinite if D > 0) represents ≥ c₀ν_DH/(64h log 2x) primes p ∤ ND in (x, x+H] with cos θ_p ∈ J; in particular
one with p ≤ 2∣D∣^A. For E/Q without CM (f = f_E): a_E(p)/(2√p) ∈ J. For all but O(N^{1−c}) squarefree n ≤ N with n ≡ 1, 2 (mod 4)
there is a prime p = a² + nb² ≤ 2(4n)^A with a_E(p)/(2√p) ∈ J (D = −4n, principal class); for n ≡ 2, 3 (mod 4), p = a² − nb²
(D = 4n, principal narrow class). *Proof.* [RS] Corollary D1 ([C] (F1), (F2) for D < 0; [R2] Corollary 5 for D > 0); each p lies
under ≤ 2 prime ideals. □

Examples (exact thresholds): λ_f(p) > 0, < 0, ≠ 0: A > 3 ([RS]); λ_f(p) > −α, α ∈ (0, 2): A > 9/4 ([RS]); θ_p ∈ [π/4, π/2]: A > 5;
∣θ_p − π/2∣ ≤ 0.05: A > 127/3; θ_p ∈ [0.4, 0.6]: A > 67/3 (CMS degrees of §1; the true M\* may be smaller).

## 5. Theorem 2: Sato–Tate equidistribution inside every class

**Theorem 2.** Let n ≥ 1 and A > 𝒜(2n). There are effective c, κ > 0 and X₀ (depending on f, n, A) such that for X ≥ X₀ all
fundamental D with ∣D∣ ≤ X (both signs) outside a set of size ≤ X^{1−c} satisfy, for every (narrow) class C, x ≥ ∣D∣^A and
H ∈ [x∣D∣^{−κ}, x], with Π_C(J) := Σ_{𝔭∈C of degree one, p∤ND, cos θ_p∈J}log p·w₀((p−x)/H):
 (i) Π_C([−1,1]) ≥ ν_DH/(16h^{(+)});
 (ii) for every interval J ⊂ [−1,1] with endpoints α ≤ β:
   ∣Π_C(J) − μ_ST(J)Π_C([−1,1])∣ ≤ (ϵ_α + ϵ_β + 1/(n+1))Π_C([−1,1]),  ϵ_α = λ_n(α)1[α > −1], ϵ_β = λ_n(β)1[β < 1];
 (iii) consequently sup_J∣Π_C(J)/Π_C([−1,1]) − μ_ST(J)∣ ≤ 5/(n+1), and sup_t∣Π_C([−1,t])/Π_C([−1,1]) − μ_ST([−1,t])∣ ≤ 3/(n+1).
In terms of A: for A > 11/3 take n := ⌈(3A−3)/2⌉ − 1 (then n ≥ 4, so 𝒜(2n) = (2n+3)/3 < A, and n+1 ≥ (3A−3)/2): the discrepancy is
≤ 10/(3A−3) (e.g. ≤ 1/3 for A ≥ 11, ≤ 0.1 for A ≥ 35); for 3 < A ≤ 11/3 the bound 10/(3A−7) is trivial. The bound (ii) is finer near ±1, where λ_n(cos θ) ≈ 2sin²θ/n.

*Proof.* Let η := 1/(64n(2n+3)(n+1)). Apply Proposition T_m with η₀ := η for every 1 ≤ m ≤ 2n (A > 𝒜(2n) ≥ 𝒜(m)) and the m = 0 lower
bound; take δ as in Theorem 1 and remove the exceptional sets. Fix D outside them, C, x, H; write S^{(m)} := S_C^{(m)} and let R be the
part of S^{(0)} from the 𝔞 ∈ C that are not degree-one primes with p ∤ ND. Then hS^{(0)} ≥ ν_DH/8, ∣hS^{(m)}∣ ≤ ην_DH ≤ 8ηhS^{(0)}
(1 ≤ m ≤ 2n), and R ≤ C(H/√x + ω(ND) + log x)log x ≤ ∣D∣^{−1/2}S^{(0)} for ∣D∣ ≥ D₀ (x ≥ ∣D∣³, H ≥ x∣D∣^{−κ}, h ≤ ∣D∣^{1/2+o(1)},
ν_D ≥ ∣D∣^{−δ}); so Π_C([−1,1]) = S^{(0)} − R gives (i).
For P = Σ_{m≤2n}c_mU_m, Σ_{𝔭}P(cos θ_p)log p·w₀((p−x)/H) = Σ_mc_mS^{(m)} − Σ_{other 𝔞}P_𝔞Λ_K(𝔞)w(N𝔞)
= c₀S^{(0)} + O(8ηΣ_{m≥1}∣c_m∣·S^{(0)}) + O(Σ_m(m+1)∣c_m∣·R). For P = P_J^± (Lemma M3), 8ηΣ_{m≥1}∣c_m∣ ≤ 1/(2(n+1)) and
Σ_m(m+1)∣c_m∣R ≤ 4(2n+1)³∣D∣^{−1/2}S^{(0)}. Since P_J⁻ ≤ 1_J ≤ P_J⁺ on [−1,1], w₀ ≥ 0 and ∣c_0^± − μ(J)∣ ≤ ϵ_α + ϵ_β:
  (μ(J) − ϵ_α − ϵ_β − 1/(2(n+1)) − 4(2n+1)³∣D∣^{−1/2})S^{(0)} ≤ Π_C(J) ≤ (μ(J) + ϵ_α + ϵ_β + 1/(2(n+1)) + 4(2n+1)³∣D∣^{−1/2})S^{(0)},
and S^{(0)}(1 − ∣D∣^{−1/2}) ≤ Π_C([−1,1]) ≤ S^{(0)}. For ∣D∣ ≥ D₀(n) this gives (ii). (iii) follows from Lemma M4 (λ_n ≤ 2/(n+1)). □

**Corollary 2.1 (A → ∞).** There are functions A_f(X) → ∞, κ_f(X) > 0, ε_f(X) → 0, effective in principle, and a set 𝔈_f of fundamental
discriminants with #{D ∈ 𝔈_f : ∣D∣ ≤ X} = o(X), such that for D ∉ 𝔈_f, every (narrow) class C, x ≥ ∣D∣^{A_f(∣D∣)} and
H ∈ [x∣D∣^{−κ_f(∣D∣)}, x], the Sato–Tate discrepancy of (cos θ_p) over the degree-one primes of C in (x, x+H] (weights log p·w₀) is
≤ ε_f(∣D∣). *Proof.* For each n let X_n be increasing with X ≥ X₀(n, 𝒜(2n) + 1/n) and X^{−c_n} ≤ 1/n for X ≥ X_n; for dyadic X let
n(X) := max{n : X_n ≤ X} and 𝔈_f := ∪_X E_{n(X)}(X). Then #{D ∈ 𝔈_f : ∣D∣ ≤ Y} ≤ Σ_{X≤Y}X/n(X) = o(Y), and Theorem 2 gives the rest with
A_f(X) := 𝒜(2n(X)) + 1/n(X), ε_f(X) := 5/(n(X)+1). □ The rates are not explicit because δ_m and X₀ depend on BTZ's constant B_m and,
in route F, on HT's constant for π′ = Sym^m f (§7).

## 6. Positive discriminants and sectors

*D > 0.* Theorems 1 and 2 were stated for both signs. For D > 0 every ingredient has its narrow-class counterpart: m = 0, [R2]
Theorem 4; Proposition T_m for m ≥ 2, route F with [RS] Proposition F⁺ ([RS] §6) and routes H, H′ with the D > 0 form of Proposition F′
(end of §3); m = 1, [RS] Theorem D⁺ (Theorem C for narrow class characters); genus theory for narrow classes and (J4), (J5) concern
the same Dirichlet twists ([R2] Theorem 4, [RS] §6); the representation step is [R2] Corollary 5.

**Theorem 3 (sectors, D < 0).** Let J, M := M\*(J), and A > 𝒜_sec(M), where 𝒜_sec(M) := 𝒜(M) for M ≠ 1 and 𝒜_sec(1) := 8/3. There
are effective c, κ > 0 such that all but O(X^{1−c}) D ∈ 𝔇 with ∣D∣ ≤ X satisfy: for every C ∈ G, every integral 𝔠 ∈ C^{−1}, every
arc I ⊂ R/(2π/w)Z of length ϑ ≥ ∣D∣^{−κ}, every x ≥ ∣D∣^A and H ∈ [x∣D∣^{−κ}, x],
  Σ_{𝔭∈C of degree one, p∤ND, x<p≤x+H, cos θ_p∈J, arg β_𝔭∈I}log p ≥ c₀wϑν_DH/(256πh),  where 𝔭𝔠 = (β_𝔭).
In particular, for almost all squarefree n ≡ 1, 2 (mod 4): a prime p = a² + nb² ≤ 2(4n)^A with arg(a + ib√n) in any arc of length
≥ (4n)^{−κ} and a_E(p)/(2√p) ∈ J.

*Proof (changes).* (i) The angular expansion of [R] Theorem 3 (ii)–(iii) with weight g supported in I (ĝ(0) ≥ wϑ/(4π),
∣ĝ(k)∣ ≤ ĝ(0), ∣ĝ(k)∣ ≤ C_rϑ(ϑ∣k∣)^{−r}), applied to Λ_K(𝔞)b_m(𝔞)w(N𝔞) for each m: the sector sum of degree m equals
ĝ(0)S_C^{(m)} + Σ_{k≠0}ĝ(k)h^{−1}Σ_{ψ∈χ_{wk}Ĝ}ψ(𝔠)S_ψ^{(m)}, S_ψ^{(m)} := Σ_𝔞Λ_K(𝔞)b_m(𝔞)ψ(𝔞)w(N𝔞) = ψ_w(Sym^m f ⊗ θ_ψ); frequencies
∣k∣ > K := ϑ^{−1}∣D∣^{ε₂} (ε₂ := δ/100) are trivial (∣b_m∣ ≤ m+1, Brun–Titchmarsh). (ii) m = 0: [R2] Theorem 3′ gives
≥ ĝ(0)ν_DH/(16h). (iii) m ≥ 1, ψ of ∞-type wk ≠ 0: θ_ψ is a cuspidal CM form of weight w∣k∣ + 1 (ψ ≠ ψ^c), GRC holds,
L(s, Sym^m f × θ_ψ) is entire with no Dirichlet (genus) factors, so only the analogues of (a), (b), (c), (f) of §2.3 occur. (b): the
trivial count over the ≤ (2K+1)h characters, ≤ ∣D∣^{1/2+κ+ε₂+o(1)} zeros per window, gives H∣D∣^{−s₁+2κ+ε₁+ε₂+o(1)}. (a): an
exceptional set E₁′ defined by Markov from the angular family densities: route F, HT with 𝒮′ := {π_ψ : D ∈ 𝔇(X), ψ ∈ 𝒜(D,M_X) ∖ Ĝ},
∣𝒮′∣ ≤ CX^{3/2}M_X log X, C(π_ψ) ≤ 12X(M_X+4)², and I₂′ ≤ C(ε)X^{(m+3)/2+ε}M_X^{m+2+ε}T^{m+2+ε} from [RS] Lemma B for each ∞-type
(absorb the fixed character χ_ℓ into the coefficients; the gamma factor of type ℓ raises the length to √q_m(∣ℓ∣ + km + ∣v∣ + 3)^{m+1}):
density (X^{m+10}M_X^{2m+14}T^{2m+10})^{e(σ)+ε}; routes H, H′, Proposition F′ with λ_ψ(n) := Σ_{N𝔞=n}ψ(𝔞) (still D-independent class I
coefficients), (LS1) by the weighted large sieve for 𝒜(D, M) behind [D2] Theorem A (threshold ≍ ∣D∣^{1/2}(M+3)), (LS2), (LS3) by
HT and PT for 𝒮′ (thresholds X^{7/2}M_X^5, X^{5/2}M_X^3, and N∣𝒮′∣^{1/2} ≤ NX^{3/4}M_X^{1/2}). Here M_X := 6X^{κ+ε₂}, so all extra factors are X^{O(m(κ+ε₂))}, absorbed by taking δ smaller by a factor C·m.
(c), (f) as before. (iv) Positivity: with η₀ := c₀/(64Σ_{m≥1}∣c_m∣), the m ≥ 1 sector sums are ≤ (9/8)η₀ĝ(0)ν_DH/h each, so the total is
≥ c₀ĝ(0)ν_DH/(16h) − (9/8)(c₀/64)ĝ(0)ν_DH/h ≥ c₀ĝ(0)ν_DH/(32h); remove non-degree-one 𝔞 (≤ c₀ĝ(0)ν_DH/(64h)). For M = 1 the per-D
input at A ≤ 3 would be Theorem C for angular characters, not written here; route F's trivial-count argument gives A > 8/3. □

## 7. Uniformity in m

| input | where used | dependence on m | explicit? |
|---|---|---|---|
| Newton–Thorne (Sym^m f cuspidal), Deligne + local–global compatibility (GRC for Sym^m f) | everywhere | qualitative, every m | — |
| conductor q_m(D) ≍ ∣D∣^{m+1} (constants N^{O(m)}); archimedean parameters ≈ (k−1)(m−2j)/2, analytic conductor ≤ ∣D∣^{m+1}N^{O(m)}(km + ∣t∣ + 3)^{2m+2} | (R2)–(R4), (J6), Lemma B, (e) | exponents linear in m | yes |
| coefficients: ∣λ_m(𝔞)∣ ≤ d_{m+1}(𝔞), ∣μ_m(𝔞)∣ ≤ d_{m+1}(𝔞), ∣b_m∣ ≤ m+1, ∣λ_m(p)∣ ≤ m+1, Euler factors of degree ≤ 2m+2 | Lemmas A, B, F1 | through divisor sums: losses (log X)^{O(m²)} and constants 8^{O(m)} | yes |
| per-D weighted large sieve ([C] Lemma 2; [D2] Lemma 5) | Lemma A, (LS1), F1(c) | none (coefficients absorbed) | yes (30/ε) |
| [RS] Lemma B (second moment) | class II, Prop F, I₂ | C(ε, m) from divisor sums of d_{2m+2} | effective, not computed |
| HT Theorem 4.1 with π′ = Sym^m f | [RS] Prop F (route F) | ≪_{ε,n′}, n′ = m+1 | **not explicit in m** |
| HT Theorem 4.1 / PT Theorem 1.1 and Corollary 1.2 with π′ = 𝟙 | (LS2), (LS3) (routes H, H′) | none: GL₂ family, k = 2 | constant not explicit but m-free |
| BTZ (J4): exponent c_m = 10⁷(m+1)⁴, constant B_m | E₂, (e) | c_m explicit; B_m effective, not explicit | **B_m not explicit** (enters δ_e through log B_m) |
| Thorner, IMRN 2025 (J5): no Landau–Siegel zero, region c/(m²[m^y log(mC(f)) + log(C(χ)(3+∣t∣))]) | (e) | explicit, c absolute | yes |
| mollifier lengths | Prop F: the HT length (family threshold); F′: V = X^{1/2}/2 per D; admissible primes p > 64(m+1)² | explicit | yes |
| δ_m, κ, exceptional exponent | §2.1 | δ_m = min{1/2840, s₁/10, ε′/10, A/(8c_m), δ_e} with s₁, ε′ ≫ (A − 𝒜(m))/m explicit, A/(8c_m) ≍ m⁻⁴, δ_e ≍ 1/(c_m + m² log(B_m/η₀)) | explicit except B_m |
| X₀ | everywhere | all of the above | effective, not explicit |

So for each fixed M every constant is effective; δ, κ and c are explicit up to log B_m; nothing is uniform enough to let M grow
with X at an explicit rate, which is why Corollary 2.1 has no rate. Routes H, H′ use the family large sieves only with π′ = 𝟙, so
the only non-explicit m-dependence there is B_m (genus twists) and Lemma B's C(ε, m).

## 8. What fixes the thresholds; limits

- *Binding point.* For the top twist m the trivial per-D count handles β < 1 − 1/(2A), and Markov over D needs the family density
  below X^{1+A(1−σ)} beyond. Route F: ((m+10)e(σ) − A(1−σ)) < 1 at σ = 1 − 1/(2A) ⟺ A > (m+7)/3. Route H′: the second-moment term
  gives A > (m+3)/3 and the mid-range terms give A > 3 and A > 11/3 (at σ = 1 − 1/(2A), ¾ + (11/2 − A)/(2A) < 1 ⟺ A > 11/3).
- *Theorem C is not needed for m ≥ 2.* Its per-D exponent (m+2)e(σ) at σ = 1 − 1/(2A) is (m+2)/(2(A+1)), below the trivial ½ only
  if A > m+1; at every threshold above A < m+1 for m ≥ 2. The brief's idea (Theorem C for every m via Newton–Thorne) is already in
  [RS], where Theorem C is stated for all m; it matters only for m = 1 (9/4).
- *The slope ⅓ per degree* comes from the second moment Σ_χ∣L(½, Sym^m f ⊗ θ_χ)∣² ≍ ∣D∣^{(m+1)/2}: the approximate-functional-equation
  length ∣D∣^{(m+1)/2} far exceeds the family size h ≍ ∣D∣^{1/2}, and the large sieve is at its diagonal. Lindelöf on average over the
  family (I₂ ≍ X^{3/2}) would remove the slope (route H′: A > 11/3 for every M; route F: A > 7/3). Out of reach (same obstruction as
  [RS] §9 for m = 2).
- *The constant 11/3 (5 without PT)* is the cost of class I lengths X^{3/2} ≤ N ≤ X^{7/2}, below the family threshold. A family large
  sieve at the family size X^{3/2} for these lengths (the "beyond the conductor" problem of [R3] §6) would remove it; [R3] Theorem 1
  reaches only N ≤ X^{1−4δ}.
- *Minorants.* The CMS polynomials are within a factor ≈ 2 (centred intervals) to ≈ 3 (intervals at an end) of the certified lower
  bound; for any given J a linear program gives M\*(J), and Theorem 1 is stated with M\*(J).
- *Equidistribution needs A → ∞.* With A fixed, only twists up to degree ≈ 3A are controlled, so the discrepancy bound stays ≍ 1/A;
  Corollary 2.1 lets A grow with ∣D∣ at a non-explicit rate.
- *Per-D statements* (every D, not almost all) would need a log-free per-D density near σ = 1 ([RS] §9 lead); not attempted. A
  per-D bound ∣D∣^{C} with C ≈ 10^{10} by Linnik's method is plausible from published inputs ([RS] NOVELTY caveat), for every interval.
- *Sectors with degree-one minorants* at A ≤ 3 need Theorem C for angular characters (not written); Theorem 3 uses 8/3 there.

## 9. Novelty (details in [NOVELTY.md](NOVELTY.md))
No source found states, for varying quadratic fields, a Sato–Tate condition on primes in every ideal class with any power of ∣D∣,
per D or on average, nor any joint (ideal class, Sato–Tate) equidistribution in the range x ≥ ∣D∣^{O(1)}. Closest: Lemke Oliver–Thorner
(IMRN 2019; Hoheisel-type short intervals and the least prime with θ_p ∈ I, for a fixed form over a fixed field), Thorner (Res. Math.
Sci. 2021; IMRN 2025), Wong (JNT 2019) and later Chebotarev–Sato–Tate papers (Amri; Malik–Prabhu; Lee; Yang; Kumar–Kumari–Mishra):
a fixed Galois extension, where ring class fields give bounds exponential in h(D); Pierce–Turnage-Butterbaugh–Wood (Invent. 2020):
Chebotarev for almost all fields of a fixed Galois group, no Sato–Tate, dihedral groups of growing order not covered; Thorner–Zaman,
Chen–Park–Swaminathan: no varying class groups with a Sato–Tate condition. [RS] (this project): Theorems D, D_M (degree ≤ 2 in
detail, all M sketched, threshold (M+7)/3).

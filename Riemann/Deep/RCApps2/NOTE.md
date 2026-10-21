# RCApps2: almost-all exponent 3/2 for primes in ring classes, and indefinite forms: proofs

Date 2026-10-21. Brief: `.scratch/orca/rcapps2-brief.md`. Files in `.scratch/rcapps2/`: `params2.py` → `params2.out`
(thresholds, explicit parameters, conditional exponents), `hyb.py` → `hyb.out` (k-th powers in the family detector),
`hyp.py` → `hyp.out` (what hypothetical family bounds would give), `axq.py`, `axq2.py` → `axq.out`, `axq2.out` (arXiv),
`gxl.py`, `gxl2.py` → `gxl*.out`, `gxl*_*.json` (evidence API). References: [C] = [Combos NOTE](../Combos/NOTE.md)
(Lemmas 2–5, Theorem 1), [D2] = [Density2 NOTE](../Density2/NOTE.md) (Lemmas 4–5, Proposition 6, Theorems A–C, Cor B′),
[R] = [RCApps NOTE](../RCApps/NOTE.md) (inputs (I1)–(I6), Proposition 1, Theorems 2–3, Corollary 4). Only §5 uses a 2026
half-plane, and says so.

## Verdicts

| Item | Result | Status |
|---|---|---|
| tool: per-D density | every family of [D2] Prop 6: Σ_{χ∈ℱ}N(σ;t,χ*) ≤ C(ε)(𝔔T)^{2(1−σ)/(3−2σ)+ε}, 0.6 ≤ σ ≤ 1 (Proposition 1). Covers ring class characters of every imaginary quadratic order, anticyclotomic characters, narrow ring class characters of real quadratic orders, characters of bounded modulus (Theorem 1′) | new; the exponent 2σ(1−σ) of [C] Thm 1 and [D2] Thms A–C drops for every σ ∈ (½,1) (by a factor up to 9/8, at σ = ¾); the T³ per window becomes T^{2(1−σ)/(3−2σ)} |
| target 1, D < 0 | A₀ = 3/2: [R] Theorems 2, 3 and Corollary 4 hold for every A > 3/2 (Theorems 2′, 3′, Corollary 4′) | new record; was A₀ = 1.653595 ([R]); literature: \|D\|^{20/3+ε} outside a density-zero set (Ditchen 1312.1502) |
| target 2, D > 0 | Theorem 4, Corollary 5: for every A > 3/2, all but O(X^{1−c}) fundamental D ∈ (0, X]: every narrow class has ≥ ν_D H/(16h⁺log 2x) prime ideals of norm in (x, x+H] for x ≥ D^A, H ≥ xD^{−κ}; every primitive indefinite form of discriminant D represents a prime ≤ 2D^A | new; for every D only D^{O(1)} (Weiss), explicitly about D^{455} (Zaman); no average result found |
| A₀ below 3/2 | needs a family bound that beats the sum over D of the per-D bounds by some X^η near σ = 5/6; the Humphries–Thorner sieve on k-th powers and intermediate mollifier lengths give nothing (§4) | NOGO with current tools; lead (§4) |
| conditional on (L) or (Hu) | least prime ≪ \|D\|^{2/(3−2B♭)+ε} = \|D\|^{1.59971} (L), \|D\|^{8/5+ε} (Hu); ring class PNT for \|D\| ≤ x^{0.62511} (L); same exponent for indefinite forms (§5) | routine; supersedes [C] Thm 2 (1.7498), Cor 3 (0.5715), [D2] Cor B′ |

## 0. Conventions
Notation of [C] §0, [D2] §0 and [R] §0. 𝔇 (𝔇⁺) is the set of negative (positive) fundamental discriminants, and 𝔇(X)
(𝔇⁺(X)) its elements with X/2 < |D| ≤ X. Put e(σ) := 2(1−σ)/(3−2σ). For D ∈ 𝔇⁺: K = Q(√D), G⁺ = Cl⁺(K) (narrow class
group), h⁺ = |G⁺|, Ĝ⁺₂ = characters of order ≤ 2 (genus characters), Ĝ⁺_c = Ĝ⁺ ∖ Ĝ⁺₂.

## 1. Per-D density with a long mollifier

**Proposition 1.** Assume the hypotheses of [D2] Proposition 6 (K quadratic, 𝒵 a Z-trivial family with the weighted large
sieve (LS), B ≤ 50log²(3𝔔), Λ ≤ 2𝔔^{1/2}, ψ₀, ℱ = ψ₀𝒵, 𝔔 = |D|(M₀+3)², gamma shifts bounded by M₀+1). For ε > 0,
t ∈ R, T = |t|+3 and 0.6 ≤ σ ≤ 1,

  Σ_{χ∈ℱ} N(σ;t,χ*) ≤ C(ε)(𝔔T)^{e(σ)+ε},  e(σ) = 2(1−σ)/(3−2σ).

*Proof.* (LS) with the single coefficient c_{O_K} = 1 gives |ℱ| ≤ (2B)(1+Λ) ≤ C𝔔^{1/2}log²(3𝔔). Each χ has ≤ c₀log(16𝔔T²)
zeros in the window ([D2] Prop 6, step 4), so we may assume 𝔔 ≥ 𝔔₀(ε). Put ε₀ := ε/10, V := 𝔔^{1/2},
Y := (𝔔T)^{1/(3−2σ)} ≥ V (as 3−2σ ≤ 2) and Z := 2Y log²Y. For each χ ∈ ℱ ∖ {1} with a zero in [σ,1] × [t−1,t+1], fix one
such zero ρ_χ = β_χ + iγ_χ; let 𝒳 be the set of these χ. [C] Lemma 5 with mollifier length V ≤ Y (valid for ℱ by [D2] Prop 6,
step 3) splits 𝒳 = 𝒳_I ∪ 𝒳_II.

*Class II, first moment.* For χ ∈ 𝒳_II, c₁Y^{σ−1/2} ≤ e∫e^{−|v−t|}|L(½+iv,χ)M_V(½+iv,χ)|dv, because |γ_χ − t| ≤ 1. Sum over
𝒳_II and apply Cauchy–Schwarz at each v over all of ℱ (not only over 𝒳_II):

  c₁Y^{σ−1/2}|𝒳_II| ≤ e∫e^{−|v−t|}(Σ_{χ∈ℱ, χ≠1}|L(½+iv,χ)|²)^{1/2}(Σ_{χ∈ℱ}|M_V(½+iv,χ)|²)^{1/2}dv.

By [D2] (2.1) with mollifier length 1, the first sum is ≤ C(ε)𝔔^{1/2+ε₀}(|v|+3)^{1+ε₀}. By (LS) ([D2] Prop 6, step 1) for
c_𝔡 = μ_K(𝔡)ψ₀(𝔡)N𝔡^{−1/2−iv} (N𝔡 ≤ V, (𝔡,f) = 1) and [C] Lemma 3, the second is ≤ (B/ε₀)(V+Λ)·C(ε₀)V^{2ε₀}
≤ C(ε)𝔔^{1/2+ε₀}. Hence |𝒳_II| ≤ C(ε)Y^{1/2−σ}𝔔^{1/2+ε₀}T^{1/2+ε₀}.

*Class I, no k-th powers.* Exactly as in the proof of [C] Theorem 1 (tail beyond Z, dyadic blocks j with 2^{j+1} > V and
2^{j−1} < Z, J ≤ 2log Z of them, removal of β), |𝒳_{I,j}| ≤ CJ²∫(1+|τ|)^{−2}Σ_{χ∈ℱ}sup_{t′∈I_τ}|P_j(χ,t′)|²dτ, where P_j has
coefficients b_𝔞 = a_𝔞e^{−N𝔞/Y}N𝔞^{−σ−it′} on 2^{j−1} < N𝔞 < 2^{j+1} and |b_𝔞| ≤ d(𝔞)2^{−(j−1)σ}. Now M := 2^{j+1} ≥ V
= 𝔔^{1/2} ≥ Λ/2, so (LS) has constant ≤ 3(B/ε₀)M, and Sobolev with [C] Lemma 3 gives
Σ_χ sup_{I_τ}|P_j|² ≤ C(ε)log⁴(3𝔔Y)·M^{2−2σ+2ε₀}. So |𝒳_I| ≤ C(ε)(𝔔T)^{ε₀}Z^{2−2σ+2ε₀}.

*Conclusion.* Σ_χ N(σ;t,χ*) ≤ c₀log(16𝔔T²)(|𝒳_I| + |𝒳_II| + 2), and Y^{2−2σ} = Y^{1/2−σ}(𝔔T)^{1/2} = (𝔔T)^{e(σ)} by the
choice of Y. Rename ε. □

**Theorem 1′.** With effective C(ε), for 0.6 ≤ σ ≤ 1, t ∈ R, T = |t|+3:
- (a) ring class characters of the order of discriminant D < 0 (the family of [C] Thm 1; (LS) = [D2] Lemma 4 with M = 0):
  Σ_{χ∈Ĝ}N(σ;t,χ*) ≤ C(ε)(|D|T)^{e(σ)+ε};
- (b) the anticyclotomic family 𝒜(D,M) of [D2] Thm A: ≤ C(ε)(|D|(M+3)²T)^{e(σ)+ε};
- (c) narrow ring class characters of the real quadratic order of discriminant D > 0 ([D2] Thm B, (LS) = [D2] Lemma 5):
  ≤ C(ε)(DT)^{e(σ)+ε};
- (d) all finite-order Hecke characters of conductor dividing 𝔪, q = min(𝔪 ∩ N) ([D2] Thm C, φ(q) cosets):
  ≤ C(ε)φ(q)(q²|d_K|T)^{e(σ)+ε}.

*Proof.* Proposition 1 with the inputs listed in the proofs of [D2] Theorems A, B, C (Thm 1 of [C] is Thm A with M = 0). □

Comparison. e(σ) < 2σ(1−σ) on (½,1), since the ratio is σ(3−2σ) > 1 there; the ratio is at most 9/8 (at σ = ¾). In terms of
the family size |ℱ| ≍ 𝔔^{1/2±ε} the bound reads |ℱ|^{4(1−σ)/(3−2σ)}: the exponent of the Humphries–Thorner detector when the
large sieve threshold equals the family size, which the per-D sieves of [C] and [D2] achieve. The Carlson detector of [C]
needed k-th powers to reach that threshold and lost the factor 𝔔^{(1/2)(2−2σ)} on the way. For Gaussian angular characters
([D2] Cor A1) Theorem 1′(b) gives K₀^{4(1−σ)/(3−2σ)+ε}, below Ricci's K₀^{(10/3)(1−σ)} for σ < 9/10 (was σ < 5/6); the
sector exponents 3/10 and 3/5 do not change, because Ricci's bound still decides the range σ ≥ 9/10.

## 2. Negative discriminants: A₀ = 3/2

**Theorem 2′.** [R] Theorem 2 holds with A₀ = 3/2: for every A > 3/2 there are effective c, κ, δ > 0 and X₀ with
#E_A(X) ≪_A X^{1−c} (E_A(X), ν_D and the count (2.1) as in [R]). For A ≥ 1.55 one may take δ = 1/2840, κ = 8·10^{−5},
c = 1.7·10^{−4}.

*Proof.* [R] §2 with three changes; the statement for A implies it for every A′ > A, so take 3/2 < A ≤ 1.6.
1. *Parameters ([R] §2.1).* Let σ₁(A) := 3/2 − 1/A, so that 2/(3−2σ) < A ⟺ σ < σ₁(A). Keep φ(σ) := (10/(3−2σ) − A)(1−σ).
   Since 3 − 2σ₁ = 2/A and 1 − σ₁ = (2−A)/(2A), φ(σ₁) = (5A − A)(2−A)/(2A) = 2(2−A), which is < 1 iff A > 3/2. Fix
   σ₀ ∈ [0.8, σ₁) with φ(σ₀) < 1 (σ₀ = 0.84 at A = 1.55) and put s₁ := (1−σ₀)(A − 2/(3−2σ₀)) > 0,
   ε′ := (1−φ(σ₀))/4, δ := min{A/4400, s₁/1.2, ε′}, κ := δ/4, ε₀ := ε₁ := δ/100; σ_g, σ_j, J, T_X, Q₀ as in [R].
   φ is decreasing on [0.8, 1]: with u = 3−2σ ∈ [1, 1.4], φ′ = A − 10/u² < 0.
2. *Zeros with β < σ₀ ([R] §2.4(b)).* Replace (I1) by Theorem 1′(a). For σ ∈ [0.6, σ₀] and a window |γ − k| ≤ 1,
   x^{σ−1}log x·N(σ) ≤ C(ε₀)|D|^{−g(σ)+ε₀}(|k|+3)^{e(σ)+ε₀}(A log|D| + 10) with g(σ) := A(1−σ) − e(σ)
   = (1−σ)(A − 2/(3−2σ)). On [0.6, σ₀] both factors of g are positive and decreasing, so g ≥ s₁. The part σ ∈ [½, 0.6] and the
   term x^{−1/2}N(0) are bounded with [C] Theorem 1(b) as in [R]: ≤ |D|^{1/2−0.4A+o(1)} ≤ |D|^{−0.1+o(1)}. Summing over
   |k| ≤ T₁ (e(σ) ≤ 4/9): E_b ≤ H|D|^{−s₁+ε₀+1.45(κ+ε₁)+o(1)} + H|D|^{−0.1+1.1(κ+ε₁)+o(1)}.
3. *Conclusion ([R] §2.5).* Use s₁ − ε₀ − 1.45(κ+ε₁) > δ/2 and 0.1 − 1.1(κ+ε₁) > δ/2 in place of [R]'s
   s₁ − ε₀ − 4(κ+ε₁) > δ/2.

Everything else is unchanged: E₁ uses [R] Proposition 1 for σ ≥ σ₀ with the same Markov step (φ(σ_j) ≤ φ(σ₀) = 1 − 4ε′);
E₂, §2.3 and §2.4 (a), (c)–(f) do not involve σ₀. `params2.py` checks every inequality of [R] §2.1 with the new σ₀, s₁, δ
(A = 1.55, σ₀ = 0.84: φ(σ₀) = 0.96412, ε′ = 0.00897, s₁ = 0.00558, δ = A/4400 = 1/2839, c = δ/2 = 1.76·10^{−4}). □

**Theorem 3′ (sectors).** [R] Theorem 3 holds for every A > 3/2 (for A ≥ 1.55 with the numbers of Theorem 2′). *Proof.* The
changes of Theorem 2′, with Theorem 1′(b) in place of [D2] Thm A in step (iv) of [R] Thm 3: since 𝔔 = |D|(M+3)² ≤ C|D|^{1+2μ}
for M ≤ 6|D|^{κ+ε₂}, μ := κ + ε₂ + o(1), the quantity s₁ becomes s₁′ := (1−σ₀)(A − 2(1+2μ)/(3−2σ₀)); φ′ is as in [R].
`params2.py`: A = 1.55, μ = 9.2·10^{−5}: φ′(σ₀) = 0.96426, s₁′ = 0.00553, all checks pass. □

**Corollary 4′.** For every A > 3/2 (for A ≥ 1.55: κ = 8·10^{−5}, c = 1.7·10^{−4}), all but O(N^{1−c}) squarefree n ≤ N
with n ≡ 1, 2 (mod 4) satisfy: for all x ≥ (4n)^A, H ≥ x(4n)^{−κ} and every arc I ⊂ R/πZ of length ≥ (4n)^{−κ}, there is a
prime p = a² + nb² in (x, x+H] with arg(a + ib√n) ∈ I. In particular the least prime a² + nb² is ≤ 2(4n)^A. (Proof of [R]
Corollary 4 with Theorem 3′.)

## 3. Positive discriminants: indefinite forms

**Proposition 2 (family density).** For X, T ≥ 3, 0 ≤ σ ≤ 1 and ε > 0:
Σ_{D∈𝔇⁺(X)}Σ_{χ∈Ĝ⁺_c(D)}N(σ,T;χ) ≤ C(ε)(XT)^{10(1−σ)/(3−2σ)+ε}.

*Proof.* [R] Proposition 1 with the following inputs. For χ ∈ Ĝ⁺_c let π_χ := AI_{K/Q}(χ). Since 𝔞𝔞^c = (N𝔞) with N𝔞 > 0
totally positive, χ^c = χ̄ ≠ χ, so π_χ is cuspidal on GL₂(A_Q), and π_χ = π_{χ′} iff χ′ ∈ {χ, χ̄} (π_χ determines its
conductor D, hence K). As χ((α)) = sgn(Nα)^a for α ≡ 1 (mod^× f) ([D2] Thm B proof), the archimedean component is the
principal series sgn^a ⊞ sgn^a: π_χ is a Maass form of eigenvalue ¼, tempered at ∞. The Satake parameters are χ(𝔭), χ(𝔭^c)
at split p and ±1 at inert p (χ((p)) = 1), so GRC holds and θ₂ = 0 is allowed in (I3). C(π_χ) ≤ 12D. Next,
|𝒮| ≤ Σ_{D≤X}h⁺(D) ≤ CX^{3/2}: h⁺ ≤ 2h and h·log ε_D = ½√D L(1,χ_D) with log ε_D ≥ log((1+√D)/2). Finally
|L^{ur}(½+iv,π_χ)| ≤ 2^{ω(D)}|L(½+iv,χ)|, and [D2] (2.1) for 𝒵 = Ĝ⁺ (𝔔 = 9D, mollifier length 1) gives
I₂ ≤ C(ε)X^{3/2+2ε}T^{2+ε}. The rest of [R] Proposition 1 is unchanged. □

**Theorem 4.** For every A > 3/2 there are effective c, κ, δ > 0 (for A ≥ 1.55: δ = 1/2840, κ = 8·10^{−5}, c = 1.7·10^{−4})
and X₀ such that for X ≥ X₀ all D ∈ 𝔇⁺ with D ≤ X outside a set of size ≪_A X^{1−c} satisfy: for every C ∈ G⁺, every
x ≥ D^A and every H ∈ [xD^{−κ}, x],

  #{𝔭 ∈ C prime of degree one : x < N𝔭 ≤ x+H} ≥ ν_D H/(16 h⁺(D) log 2x),

where ν_D is defined as in [R] Theorem 2 (so ν_D ≥ cD^{−δ/2}log^{−1}D, and ν_D = 1 unless some χ_d, d | D, |d| ≤ D^δ, has a
Siegel zero).

*Proof.* [R] §2 together with the changes of Theorem 2′, under this dictionary.
- Ĝ, h, Ĝ_c, Ĝ₂ become Ĝ⁺, h⁺, Ĝ⁺_c, Ĝ⁺₂. Genus theory for narrow classes of a fundamental D > 0 (Gauss; e.g. Hecke,
  *Vorlesungen*, §48; Buell, *Binary quadratic forms*, Ch. 4): Ĝ⁺₂ = {χ_{d,D/d}}, one character per unordered pair
  {d, D/d} with d, D/d ∈ {1} ∪ 𝔇 ∪ 𝔇⁺ (necessarily of equal sign), χ_{d,D/d}(𝔞) = χ_d(N𝔞) when (N𝔞, d) = 1, and
  L(s,χ_{d,D/d}) = L(s,χ_d)L(s,χ_{D/d}) (Kronecker). So (G) and (2.2) of [R] hold verbatim.
- (I1): Theorem 1′(c) for 0.6 ≤ σ ≤ 1; for every σ, Σ_{χ∈Ĝ⁺}N(σ;t,χ) ≤ c₀h⁺log(144DT²), and
  h⁺ ≤ 1.04D^{1/2}L_D ≤ 5D^{1/2}log²(3D) ([D2] §0), in place of [C] Theorem 1(b) and (F3).
- [R] Proposition 1: Proposition 2. The Markov step and φ are literally the same.
- (I4), (I5) concern the Dirichlet characters χ_d and do not change; in E₂, each d divides ≤ X/|d| + 1 elements of 𝔇⁺(X).
- (I6) holds for Hecke characters of K (Iwaniec–Kowalski Thm 5.11); [R] §2.4(f) becomes
  Σ|R| ≤ C(h⁺ + 2^{ω(D)})(H/x + x^{−1/2})log(Dx) ≤ HD^{−1}; in §2.4(c) each window holds ≤ Ch⁺log(D(|k|+3)) zeros.
- [R] §2.5: prime powers, inert and ramified primes contribute ≤ C(H/√x + 1)log²x ≤ ν_D H/(16h⁺), since √x ≥ D^{3/4}.
Every exponent condition involves only A, σ₀, δ, κ, ε₀, ε₁, so the numbers of Theorem 2′ work. □

**Corollary 5.** Let A > 3/2 and let D ∈ 𝔇⁺ be outside the exceptional set of Theorem 4. Every primitive indefinite binary
quadratic form Q of discriminant D represents ≥ ν_D H/(32h⁺ log 2x) primes p ∤ D in (x, x+H] for all x ≥ D^A and
H ∈ [xD^{−κ}, x]; in particular a prime p ≤ 2D^A. Hence, for every A > 3/2 (A ≥ 1.55: c = 1.7·10^{−4}), all but
O(N^{1−c}) squarefree n ≤ N with n ≡ 2, 3 (mod 4) have a prime p = a² − nb² in every such interval with x ≥ (4n)^A, and
the least one is ≤ 2(4n)^A.

*Proof.* The dictionary of [D2] Cor B′: replace Q by an equivalent form (a,b,c) with a > 0; 𝔞_Q := Za + Z(−b+√D)/2 is an
invertible ideal of norm a. If 𝔭 lies in the narrow class of 𝔞_Q^{−1}, then 𝔭𝔞_Q = (α) with α totally positive, α ∈ 𝔞_Q,
α = ua + v(−b+√D)/2 with u, v ∈ Z, and a·Q(u,−v) = N(α) = a·N𝔭 > 0, so Q(u,−v) = N𝔭 = p. Each p lies under ≤ 2 prime ideals. For x² − ny²
take D = 4n and the principal narrow class. □

## 4. Why 3/2, and what would lower it
- *The binding point.* At A = 3/2 and σ = 5/6 three exponents meet: the sum over D of Theorem 1′(a), 1 + e(5/6) = 5/4; the
  family bound of [R] Prop 1, 10(1/6)/(4/3) = 5/4; and the Markov bar 1 + A(1−σ) = 5/4. For A < 3/2 neither bound works on
  [σ₁(A), σ₂(A)], where φ(σ₂) = 1 (A = 1.45: [0.810, 0.835]); there the sum of the per-D bounds misses the bar by at most
  X^{0.009}.
- *Class II cannot improve.* Σ_{χ∈Ĝ}|L(½,χ)M_V(½,χ)| ≍ h in mean, so its Cauchy–Schwarz bound is sharp up to |D|^ε. A
  saving must come from class I: a large sieve over all D ≤ X, Σ_DΣ_{χ∈Ĝ_c(D)}|Σ_{N𝔞∼N}c_𝔞χ(𝔞)|² ≤ X^{1−η}N‖c‖², at
  N ≈ X^{3/4}, for coefficients that depend only on N𝔞 (true for the mollified coefficients on primitive ideals).
- *Tried, no gain.* (i) Humphries–Thorner's sieve (threshold X^{7/2}) on k-th powers of the class I blocks inside the family
  detector: at the handover block N = X^{A/2} the best k costs X^{1.375} against the bar X^{1.275} (A = 1.45; `hyb.py`).
  (ii) Mollifier lengths between X^{1/2} and X^{7/2} in the Humphries–Thorner detector: class I costs a factor X for tail
  lengths below X^{5/2}, so the two endpoints (sum of per-D bounds; [R] Prop 1) are optimal. (iii) Pointwise or averaged
  Rankin–Selberg bounds for the D ≠ D′ terms of the dual sieve give thresholds ≥ X^{7/2}.
- *Payoffs* (`hyp.py`, exponents in log X): a class I saving X^{0.01} gives A₀ ≈ 1.472; a family fourth moment
  Σ_DΣ_χ|L(½,χ)|⁴ ≪ X^{3/2+θ} with θ < 1/4 gives A₀ < 3/2 (θ = 0: 1.294); a family sieve constant N² + X^{3/2} gives 1.208;
  the optimal constant X^{3/2+ε} for N ≤ X^{1−η} gives A₀ = 1, the GRH exponent.
- *Lead (trace formula).* Two primitive ideals of norms n, n′ lie in one class of O_D iff t² + |D|y² = 4nn′ has a compatible
  solution with y ≠ 0. Summed over D with the weight h(D), the class I second moment becomes Σ_{n,n′}c_n c̄_{n′}𝓦(nn′), with
  𝓦(m) a partial Hurwitz class number sum Σ_y Σ_t h((t² − 4m)/y²) over (4m − t²)/y² ∈ 𝔇(X). Zagier's form of the
  Eichler–Selberg trace formula expands such sums into Hecke eigenvalues λ_j(m) of level-one Maass forms, with spectral cutoff
  ≈ N/X^{1/2}; the Deshouillers–Iwaniec spectral large sieve then bounds the cuspidal part by about X^{3/2}‖c‖². Open step:
  identify the identity and Eisenstein terms with the genus-character part (a Siegel–Weil-type identity), uniformly, with the
  cutoffs to fundamental D and to D ∈ [X/2, X]. Not attempted.

## 5. Conditional corollaries (assume (L): finite-order Hecke L-functions do not vanish for Re s > B♭ = 0.874887; or (Hu): 7/8)
- [C] Theorem 2 with Theorem 1′(a): for σ ∈ [0.6, B♭] and x = |D|^A y, x^σN_D(σ;k) ≤ C x|D|^{−(1−σ)(A−2/(3−2σ))+ε′}
  y^{−(1−σ)}(|k|+3)^{1/2}, and 2/(3−2σ) ≤ 2/(3−2B♭) on that range. So every primitive positive definite form of discriminant
  D represents a prime ≤ C(ε)|D|^{2/(3−2B♭)+ε}, 2/(3−2B♭) = 1.599710 (was 2B♭ = 1.749773); under (Hu), 8/5 (was 7/4). Same
  for x² + ny² and dihedral Chebotarev ([C] Cor 4).
- [C] Corollary 3: π_C(x) = (Li(x)/h)(1 + O(x^{−c})) uniformly for |D| ≤ x^{(3−2B♭)/2−ε} = x^{0.625113−ε} (was x^{0.571503}),
  since |D|^{e(σ)} ≤ x^{(1−σ)−c} for σ ≤ B♭ when |D| ≤ x^{(3−2B♭)/2−ε}; under (Hu), x^{5/8−ε}.
- [D2] Cor B′ (indefinite forms) and Cor C′ (bounded q): exponent 2/(3−2B♭) in place of 2B♭.
These rest on unrefereed 2026 claims and are routine: COROLLARIES.md, not DISCOVERIES.md.

## 6. Novelty audit (2026-10-21)
- arXiv API, 34 queries (`axq.out`, `axq2.out`): indefinite (binary) quadratic forms × least prime / primes / represented;
  narrow class (group) × prime ideal; real quadratic × least prime (ideal) / Chebotarev / almost all / Linnik / short
  intervals; least prime × ideal class; Chebotarev × dihedral / ring class field / almost all; ring class / class group
  (characters, L-functions) × zeros / zero density; zero density × dihedral / real quadratic / Hecke × imaginary quadratic;
  Humphries–Thorner; au:Thorner × zero density; mollifier × zero density × family. Relevant hits: Ditchen 1312.1502
  (positive definite forms only, on average), Zaman 1710.08914 (upper bounds) and 1506.01635 (exceptional-zero case:
  D^{9.5}, D^6), Sardari 1802.06193 and 1803.03218 (GRH), Gaudet 2408.01610 (sieve), Pascadi–Thorner 2508.14888, Blomer–Thorner 2410.17158
  (relative trace formula; GL(m) families, not class groups), Hashimoto 0807.0056 / 1003.3716 (class number sums of indefinite
  forms, no primes).
- GXL evidence API, 6 claims (`gxl.out`, `gxl2.out`): only Ditchen (20/3, density-zero exceptions, definite forms), Zaman
  (exceptional case) and Sardari (GRH) support anything; nothing on indefinite forms on average over D, and no per-D density
  for class group families below 2σ(1−σ).
- The tool (Proposition 1) is a known detector (Humphries–Thorner, Compos. Math. 2024, §5) fed with the project's per-D
  large sieves; the statements of Theorem 1′, Theorems 2′–4 and Corollaries 4′–5 were not found.

## 7. Dead ends (NOGO candidates)
- Below A₀ = 3/2 with the present inputs: k-th powers with the Humphries–Thorner sieve, intermediate mollifier lengths, and
  Rankin–Selberg bounds for the D ≠ D′ terms all fail (§4).
- Per-D beyond e(σ) with second moments only: class II is sharp, the class I count Y^{2−2σ} is the diagonal of the sieve,
  Halász–Montgomery off-diagonals are ≈ h, and subconvexity for class group L-functions (Duke–Friedlander–Iwaniec) does not
  enter (in the extremal configuration |L(½,χ)| ≈ Y^{(2σ−1)/4}, about |D|^{0.125} at σ = 5/6, far below |D|^{1/4}). The
  per-D fourth moment is only |D|^{1+ε}.
- Sectors for real quadratic fields (hyperbolic angles): the Grössencharacter large sieve has constant N + D^{1/2}(M+1)log ε_D
  ([D2] §8), so only small regulators are covered. Not pursued.

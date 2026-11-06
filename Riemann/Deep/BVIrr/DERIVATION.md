# DERIVATION — Breuillard–Varjú's random-polynomial theorem without the Riemann hypothesis

Agent bvirr, 2026-11-06. Every step and estimate behind REPORT.md. TeX line numbers refer to the arXiv sources
(BV = 1810.13360, file `irreducibility-via-large-primes.tex`; BSHKP = 2308.04878 `main.tex`; BBK = 2609.40123
`rand_pols_final.tex`; VX = 2511.04240 `Fekete_Polynomial.tex`; BGS = 2502.17218 `tridiag009__1_.tex`;
Eb = 2008.01223 `main.tex`; FJSS = 2106.04049 `main.tex`; DKN = 2304.12075 `main.tex`).
BV's numbering (one counter for all environments): Thm 1 `th:main` l.226, Thm 2 `th:GRH` l.283, Cor 3 `cr:00` l.333,
Cor 4 `cr:01` l.343, Thm 5 `th:mild-hyp1` l.378, Thm 6 `th:mild-hyp2` l.419, Thm 7 `th:exp-zero` l.438,
Prop 9 `pr:prime-sum-GRH` l.650, Thm 10 `th:explicit` l.704, Lemma 12 `lm:Stark` l.851, Prop 13 `pr:prime-sum` l.904,
Prop 16 `pr:ApBp` l.1072, Prop 19 `pr:Bp-sum-GRH` l.1125, Lemma 21 `lm:discriminant` l.1154.

Notation. K a number field, n_K = [K:Q], Δ_K ≥ 1 the absolute value of its discriminant, r_1, r_2 its signature,
Λ_K(N) = Σ_{N𝔭^j = N} log N𝔭 (N ≥ 2), ψ_K(x) = Σ_{N ≤ x} Λ_K(N), A_K(p) = #{𝔭 : N𝔭 = p}. x = e^X.
ℓ_K(t) := log Δ_K + n_K log(t + 2) for t ≥ 0. C_0, C_1, ... are absolute constants.

## 0. Inputs

- (CT) Cureton–Thorner (Hexagon 2610.00200v2; Lean theorem `dedekindZeta_ne_zero`, github.com/AxiomMath/DedekindZetaZeroFree,
  Main.lean l.9–14, no sorry, checked in our R399): for every number field K, ζ_K(s) ≠ 0 whenever Re s > 31/32, s ≠ 1.
- (FE) The completed function ξ_K(s) = Δ_K^{s/2} Γ_R(s)^{r_1} Γ_C(s)^{r_2} ζ_K(s) is meromorphic with simple poles at
  s = 0, 1 only and ξ_K(s) = ξ_K(1 − s), ξ_K(s̄) = conj ξ_K(s). Its zeros are the nontrivial zeros ρ = β + iγ of ζ_K,
  and 0 < β < 1. The other zeros of ζ_K (trivial zeros, from the poles of the Γ-factors) are s = 0 of order r_1 + r_2 − 1,
  s = −2j of order r_1 + r_2 and s = −(2j − 1) of order r_2 (j ≥ 1); every multiplicity is ≤ n_K.
  With (CT): ρ nontrivial ⇒ β ≤ 31/32 and (applying CT to 1 − ρ̄, also a zero) β ≥ 1/32. Hence |ρ| ≥ 1/32 and
  |1 − ρ| ≥ 1/32 for every nontrivial zero; trivial zeros have |1 − ρ| ≥ 1.
- (EF) Explicit formula, BV Thm 10 (l.704–743): if g ∈ C²(R) has compact support in (0, ∞) and G(s) := ∫ e^{su} g(u) du,
  then Σ_{N ≥ 2} Λ_K(N) g(log N) = G(1) − Σ_ρ G(ρ), the sum running over all zeros of ζ_K, trivial ones included, with
  multiplicity (absolutely convergent: |G(σ + it)| = O_σ(|t|^{−2})). BV prove it by Perron inversion and contour shift;
  it is the classical Weil explicit formula for ζ_K.
- (ZC) Zero count, Lagarias–Odlyzko 1977 Lemma 5.4, quoted by BV l.892–897: for every t ∈ R,
  #{ρ nontrivial : |γ − t| ≤ 1} ≤ C_0 ℓ_K(|t|). (For K = Q it reads C_0 log(|t| + 2).)
- (BV-free) Everything in BV outside Prop 9, Prop 19, the proof of Thm 2 (§8.1, l.2411–2470) and the statements that quote
  them is free of any hypothesis on zeros: checked by reading every RH/hypothesis occurrence (REPORT, Checks). In particular
  Prop 16 (l.1072), Lemma 21 (l.1154), Prop `pr:exp-root-2` / Cor `cor:exp-root-2` (l.1978–2010, random walks),
  Prop `pr:small-Mahler` (l.2142), Prop `pr:powers` (l.2261; its PNT in progressions, l.2314, is for moduli
  k ≤ d^{1/5} ≤ (log x)^{O(1)}, effective and unconditional) and Lemma `lm:cycl-div` (l.2480).

## 1. Prop 9′: the dyadic prime ideal theorem from the strip

**Proposition 9′.** There is an absolute C_2 such that for every number field K and every X ≥ 2,
  |ψ_K(e^X) − ψ_K(e^X/2) − e^X/2| ≤ C_2 x^{31/32} (1 + X)(log Δ_K + n_K(1 + X)).
In the brief's form: |ψ_K(x) − ψ_K(x/2) − x/2| ≤ C x^{31/32}(log x · log Δ_K + n_K log² x) (x ≥ e²). Summing dyadic
pieces gives the same bound for |ψ_K(x) − x|.

*Proof.* Put a = X − log 2, b = X, S := Σ_{a < log N ≤ b} Λ_K(N) = ψ_K(x) − ψ_K(x/2), and fix η ∈ (0, 1/4].

Step 1 (kernel). φ_η := I_{η/4} * I_{η/4} * I_{η/4} * I_{η/4}, where I_r = r^{−1} 1_{[−r/2, r/2]}. Then φ_η ≥ 0,
∫ φ_η = 1, supp φ_η = [−η/2, η/2], φ_η ∈ C² (piecewise cubic), φ_η(±η/2) = 0, and
  Φ_η(s) := ∫ e^{sv} φ_η(v) dv = (sinh z / z)^4, z = sη/8,
since ∫ e^{sv} I_r(v) dv = sinh(sr/2)/(sr/2).

Step 2 (two bounds). For z ∈ C: |sinh z| ≤ cosh(Re z) ≤ e^{|Re z|}, and sinh z / z = ∫_0^1 cosh(tz) dt gives
|sinh z / z| ≤ e^{|Re z|}. Hence, with σ = Re s,
  |Φ_η(s)| ≤ e^{|σ|η/2} · min(1, (8/(η|s|))^4).

Step 3 (sandwich). J_+ = [a − η/2, b + η/2], J_− = [a + η/2, b − η/2] (nonempty: b − a = log 2 > η), g_± := 1_{J_±} * φ_η.
Both are C², supported in [a − η, b + η] ⊂ (1, ∞) (X ≥ 2), and
  0 ≤ g_− ≤ 1_{(a, b)} ≤ 1_{(a, b]} ≤ 1_{[a, b]} ≤ g_+ .
(For u ∈ [a, b] and |v| ≤ η/2, u − v ∈ J_+, so g_+(u) = 1; g_−(u) > 0 forces u − v ∈ J_− for some |v| < η/2, so a < u < b;
g_− ≤ 1 everywhere.) Since Λ_K ≥ 0: Σ Λ_K(N) g_−(log N) ≤ S ≤ Σ Λ_K(N) g_+(log N).

Step 4 (transforms). G_±(s) = E_±(s) Φ_η(s) with E_±(s) = ∫_{J_±} e^{su} du = (e^{s b_±} − e^{s a_±})/s (s ≠ 0),
where a_± = a ∓ η/2, b_± = b ± η/2.

Step 5 (main term). e^{b_+} − e^{a_+} − x/2 = x(e^{η/2} − 1) + (x/2)(1 − e^{−η/2}) ∈ [0, 0.82ηx] (e^{η/2} − 1 ≤ 0.57η) and
e^{b_−} − e^{a_−} − x/2 = −x(1 − e^{−η/2}) − (x/2)(e^{η/2} − 1) ∈ [−0.82ηx, 0]; Φ_η(1) = (sinh(η/8)/(η/8))^4 ∈ [1, 1 + η²/64].
Hence |G_±(1) − x/2| ≤ ηx.

Step 6 (nontrivial zeros). For ρ = β + iγ, 1/32 ≤ β ≤ 31/32 (FE + CT): |E_±(ρ)| ≤ (e^{β b_±} + e^{β a_±})/|ρ|
≤ 2e^{η/2} x^{31/32}/|ρ| ≤ 2.28 x^{31/32}/|ρ|, and |Φ_η(ρ)| ≤ e^{η/2} min(1, (8/(η|ρ|))^4) by Step 2. So
  |G_±(ρ)| ≤ 2.6 x^{31/32} |ρ|^{−1} min(1, (8/(η|ρ|))^4).
Let Z_j = #{ρ : j ≤ |γ| < j + 1}; these zeros lie in {|γ − (j + ½)| ≤ 1} ∪ {|γ + (j + ½)| ≤ 1}, so Z_j ≤ 2C_0 ℓ_K(j + 1)
by (ZC). Put J = ⌈8/η⌉. Then
  Σ_ρ |ρ|^{−1} min(1, (8/(η|ρ|))^4) ≤ 32 Z_0 + Σ_{j=1}^{J} Z_j / j + (8/η)^4 Σ_{j>J} Z_j j^{−5},
using |ρ| ≥ 1/32 for j = 0 and |ρ| ≥ |γ| ≥ j for j ≥ 1. The middle sum is ≤ 2C_0 ℓ_K(J + 1)(1 + log J).
For j ≥ J, (j + 3)/(J + 3) ≤ j/J, so ℓ_K(j + 1) ≤ ℓ_K(J + 1) + n_K j/J; with Σ_{j>J} j^{−5} ≤ 1/(4J⁴),
Σ_{j>J} j^{−4} ≤ 1/(3J³) and (8/η)^4 ≤ J⁴, the last sum is ≤ 2C_0(ℓ_K(J + 1)/4 + n_K/3) ≤ 2C_0 ℓ_K(J + 1).
Since J ≤ 9/η and J + 3 ≤ 12/η:
  Σ_{ρ nontrivial} |G_±(ρ)| ≤ 5.2 C_0 x^{31/32} (34 + log(9/η)) (log Δ_K + n_K log(12/η)).

Step 7 (trivial zeros). ∫ g_± = |J_±| ≤ log 2 + η ≤ 0.95, and a − η ≥ 1.05. So |G_±(0)| ≤ 0.95 and
|G_±(−j)| ≤ e^{−j(a − η)} · 0.95 ≤ 0.95 e^{−j}; with multiplicities ≤ n_K the trivial zeros contribute ≤ 1.6 n_K.

Step 8 (conclusion). By (EF) and Steps 3, 5–7,
  |S − x/2| ≤ ηx + 5.2 C_0 x^{31/32}(34 + log(9/η))(log Δ_K + n_K log(12/η)) + 1.6 n_K.
Take η = x^{−1/32}/4 (≤ 1/4): ηx = x^{31/32}/4, log(9/η) = X/32 + log 36, log(12/η) = X/32 + log 48. This is
≤ C_2 x^{31/32}(1 + X)(log Δ_K + n_K(1 + X)). ∎

**Remark 1.1 (any strip).** If ζ_K has no zero with Re s > 1 − δ (0 < δ ≤ 1/2), the same proof (with |ρ| ≥ δ in the j = 0
term and η = x^{−δ}/4) gives |S − x/2| ≤ C x^{1−δ} (δ^{−1} + X)(log Δ_K + n_K(1 + X)). This matches Frei–Sofos–Widmer
2610.18491 Thm 1.1 (l.194–206: |π_C − (#C/#G)Li| ≤ c_1 (#C/#G) x^{1−δ} log(Δ_L^{1/δ} x^{n_L}); take L = K, as
Baldwin–Sivek 2611.06705 Cor `cor:landau`, l.379–387, do). Our derivation uses only BV's own Thm 10 and (ZC).

**Remark 1.2 (real Dirichlet characters).** The sandwich of Step 3 needs nonnegative coefficients, so it is not applied to
Λχ directly. For a fundamental discriminant D ≠ 1 and K = Q(√D): Λ_K = Λ·(1 + χ_D) exactly (p split: two primes of norm p;
inert: one of norm p²; ramified: one of norm p), hence ψ(x, χ_D) = ψ_K(x) − ψ(x), and ζ_K = ζ·L(s, χ_D). If ζ and L(s, χ_D)
both have no zero with Re s > T = 1 − δ, Remark 1.1 for K and for Q gives
|ψ(x, χ_D) − ψ(x/2, χ_D)| ≤ C x^{T}(δ^{−1} + X)(log|D| + 1 + X). This is all §6.3 uses (only quadratic characters occur).

## 2. Prop 9′ in BV's normalisation

BV's weight: h_X(u) = 2e^{−X} 1_{(X − log 2, X]}(u). For a rational prime p, Λ_K(p) = A_K(p) log p. The other N in
(x/2, x] with Λ_K(N) ≠ 0 are prime powers p^e, e ≥ 2; there are ≤ √x + (X/log 2) x^{1/3} ≤ 2X√x of them (X ≥ 2), and
Λ_K(p^e) = Σ_{𝔭 | p, f_𝔭 | e} f_𝔭 log p ≤ n_K log p ≤ n_K X. Hence S − Σ_{x/2<p≤x} A_K(p) log p ∈ [0, 2n_K X² √x] and,
by Prop 9′,

**(9′)** Σ_p A_K(p) log(p) h_X(log p) = 1 + O( e^{−X/32} (1 + X)(log Δ_K + n_K(1 + X)) )   (X ≥ 2, every K).

For K ≠ Q, Minkowski (n_K ≤ C log Δ_K) turns this into BV's shape 1 + O(X² log Δ_K e^{−X/32}): BV Prop 9 (l.659–661) with
e^{−X/32} for e^{−X/2}. For K = Q: w := Σ_p log(p) h_X(log p) = 1 + O(X² e^{−X/32}).
(Slip in BV Prop 9 as printed: for K = Q, log Δ = 0 and the stated error vanishes; the Grenié–Molteni input l.670–675 has
the term d log² x, so the right form is X²(log Δ_K + n_K)e^{−X/2}. VX's restatement, 2511.04240 l.384–394, already writes
max(log|Δ|, 1). Harmless.)

## 3. Prop 19′ (BV Prop 19 without RH)

**Proposition 19′.** There is an absolute D_0 such that for d ≥ D_0, m ≥ 1, X ≥ m d^{1/10}, every P ∈ Z[x] of degree ≤ d with
coefficients in [−H, H], H = exp(d^{1/10}), and every κ ∈ (0, 1/100) (admissibility as in BV §4 with parameters (κ/m, mX)):
  Σ_p B_P(p)^m log(p) h_X(log p) = |Ω^m/G| + O(e^{−X/40}).

*Proof.* BV's proof (l.1198–1231) with (9′) in place of Prop 9.
(a) Main sum. For each orbit O ∈ Ω^m/G, K_O = Q(x_1, …, x_m) has n_{K_O} = |O| ≤ d^m (the stabiliser of a representative
fixes K_O) and log Δ_{K_O} ≤ 2m d^m log(Hd) ≤ 4m d^{m + 1/10} (Lemma 21, l.1166–1168; log(Hd) ≤ 2d^{1/10}). There are
≤ |Ω^m| ≤ d^m orbits. By (9′),
  Σ_O Σ_p A_{K_O}(p) log(p) h_X(log p) = |Ω^m/G| + E_1,  |E_1| ≤ C_3 e^{−X/32} d^m (1 + X)(4m d^{m+1/10} + d^m(1 + X)).
As 4m d^{m+1/10} ≤ 4 d^m X (since X ≥ m d^{1/10}) and X ≥ 2: |E_1| ≤ 9C_3 d^{2m} X² e^{−X/32}.
Derivation of the discriminant bound (BV Lemma 21, RH-free): for Q | P, M(Q) ≤ M(P) ≤ ‖P‖₂ ≤ H(d + 1)^{1/2} (Mahler measure is
multiplicative and ≥ 1 on Z[x] \ {0}; Landau), so Mahler's |Δ_Q| ≤ (deg Q)^{deg Q} M(Q)^{2 deg Q − 2} gives |Δ_Q| ≤ (Hd)^{2d};
for a root α of an irreducible factor Q, |Δ_{Q(α)}| ≤ |Δ_Q| (the invariant order of a primitive Q has discriminant Δ_Q:
Simon 2001, Del Corso–Dvornicich–Simon 2005; even the crude route through the monic polynomial of aα, a the leading
coefficient, costs only a factor H^{d²}, i.e. d^{2.1} in log Δ, harmless below); Toyama's |Δ_{L_1L_2}| ≤ |Δ_{L_1}|^{[L_1L_2:L_1]} |Δ_{L_2}|^{[L_1L_2:L_2]} by induction over Q(α_1), …, Q(α_m)
gives |Δ_{K_O}| ≤ ∏_i |Δ_{Q(α_i)}|^{[K_O:Q(α_i)]} ≤ (Hd)^{2d·m d^{m−1}}.
(b) Discrepancy (RH-free, BV l.1205–1231). By Prop 16, B_P(p)^m = Σ_O A_{K_O}(p) unless p divides Δ_{P̃} or Res(P̃, R) for an
exceptional R (irreducible, deg R ≤ 10mX, M(R) ≤ e^{κ/m}). There are ≤ e^{X/10} such R (Dubickas–Konyagin, as BV state at
l.2041: polynomials of degree ≤ 10mX and Mahler measure ≤ e^{κ/m}, κ < 1/100). From Res(A, B) = a^{deg B} ∏ B(α_i) and
|B(α)| ≤ ‖B‖₁ max(1, |α|)^{deg B}, ‖B‖₁ ≤ 2^{deg B} M(B): |Res(A, B)| ≤ 2^{deg A·deg B} M(B)^{deg A} M(A)^{deg B}; with A = P̃,
B = R: log|Res(P̃, R)| ≤ 10mX(d log 2 + d^{1/10} + log d) + d ≤ 11 m d X (d ≥ D_0). A prime p > e^X/2 dividing an integer
N ≠ 0 uses ≥ X − log 2 ≥ X/2 of log|N|, and log|Δ_{P̃}| ≤ 2d log(Hd) ≤ 4d^{11/10}; so at most
(2/X)(4d^{11/10} + e^{X/10}·11 m d X) ≤ 8d + 22 m d e^{X/10} ≤ 23 m d e^{X/10} primes contribute (X ≥ d^{1/10}); m is kept.
Each contributes ≤ d^m · X · 2e^{−X} (0 ≤ B_P(p)^m ≤ d^m, 0 ≤ Σ_O A_{K_O}(p) ≤ |Ω^m| ≤ d^m). So
|E_2| ≤ 46 m d^{m+1} X e^{−0.9X}. No zeros enter. (The cruder 40 m d² X resultant bound gives ≪ m d^{m+2} X e^{−0.9X}; same use.)
(c) Margin. 9C_3 d^{2m} X² e^{−X/32} ≤ ½e^{−X/40} ⟺ X/160 ≥ 2m log d + 2 log X + log(18C_3). Since X ≥ m d^{1/10},
2m log d ≤ X/320 once d^{1/10} ≥ 640 log d (true for d ≥ 10^{48.6}), and 2 log X + log(18C_3) ≤ X/320 once
X ≥ X_0(C_3) (true as X ≥ d^{1/10}). And 46 m d^{m+1} X e^{−0.9X} ≤ ½e^{−X/40} ⟺ log(92m) + (m + 1) log d + log X ≤ 0.875X:
here (m + 1) log d ≤ 2m log d ≤ X/4 once d^{1/10} ≥ 8 log d, and log(92m) + log X ≤ 2 log X + 5 ≤ X/4 (m ≤ X), uniformly in m. ∎
(If X also satisfies BV's window (condiX), l.1980–1983, X > 100 m² d^{1/10}, the margin needs only
d^{1/10} ≥ 7 log d + O(1), i.e. d ≥ about 10^{26}; see params_output.txt.)

## 4. BV Theorem 2 without RH (redo of §8.1, l.2411–2470)

**Theorem 2′.** There are absolute effective c, C > 0 such that the following holds. Let P = A_d x^d + … + A_0 have independent
coefficients, A_1, …, A_{d−1} with common law μ, all coefficients bounded by exp(d^{1/10}) a.s., ‖μ‖₂² < 1 − τ.
(1) If dτ⁴ ≥ C(log d)², then with probability ≥ 1 − exp(−cτd^{1/2}/log d): P = Φ P̃, P̃ irreducible, Φ a product of cyclotomic
polynomials and x^m, deg Φ ≤ (C/τ)√d. (2) If dτ⁴ ≥ C(log d)^{12}, then with probability ≥ 1 − exp(−cτd^{1/2}/(log d)²) moreover
Gal(P̃) ⊇ Alt(deg P̃). In particular both hold for every fixed τ and d ≥ d_0(τ). No hypothesis on zeros.

*Size hypothesis (RH-independent; it affects BV's own proof equally).* BV state Thm 2 for d ≥ C/τ⁴ (l.288) and use "τ⁴d large"
(l.2430, l.2463). But the lower end of (condiX), 100 m² max{κ^{−1}, τ^{−1}, d^{1/10}} with κ = τ/100, forces X > 10⁴m²/τ. With
X = (c_0^{1/2}/10) τ d^{1/2}/(m^{1/2} log(md)) (BV's choices; m^{1/2} log(md) = log d for m = 1) this reads
  dτ⁴ > (10^{10}/c_0) m⁵ (log md)²,  i.e. dτ⁴ ≥ C(log d)² (m = 1) and dτ⁴ ≥ C(log d)^{12} (m = ⌈30(log d)²⌉).
For τ = (C/d)^{1/4} one has X/(10⁴m²/τ) → 0, so "d ≥ C/τ⁴" does not suffice. Under the stated hypotheses the rest of (condiX)
holds for d ≥ D: X > 100 m² d^{1/10} (since τ ≥ (C/d)^{1/4}(log d)^{1/2} resp. (log d)³ gives X ≫ C^{1/4} d^{1/4}(log d)^{−1/2}
resp. C^{1/4} d^{1/4} log d), and X < τ²d/(100 C_0 (m log md)³) ⟺ τd^{1/2} > 10 C_0 c_0^{1/2} m^{5/2}(log md)² ≍ (log d)^{O(1)}.
Also X > 10⁴m²/τ ≥ 100 m² gives B_m ≤ 2^{m²} ≤ e^{X/100} directly, and τ < 1 gives d ≥ C ≥ D_0 for Prop 19′.

*Proof of (1).* BV's choices: κ = τ/100, m = 1, X = τ(c_0 d/100)^{1/2}/log d (l.2430), g = h_X.
- Cor `cor:exp-root-2` (RH-free): with probability ≥ 1 − 2exp(−X/6), |Σ_p B_P(p) log(p) h_X(log p) − w| < ½; with w → 1
  (PNT) this is (good-poly), l.2424–2427: |Σ_p B_P(p) log(p) h_X(log p) − 1| ≤ 2/3.
- Prop 19′ applies: X ≥ d^{1/10} and d ≥ D_0 hold when dτ⁴ ≥ C(log d)² with C large (size hypothesis above; X lies in
  (condiX)). It gives Σ_p B_P(p) log(p) h_X(log p) = |Ω/G| + O(e^{−X/40}).
- Hence | |Ω/G| − 1 | ≤ 2/3 + O(e^{−X/40}) < 1, so |Ω/G| = 1: P̃ (the product of the admissible irreducible factors)
  is irreducible. This replaces l.2433–2441, the only place RH enters part (1).
- The rest is RH-free and unchanged (l.2443–2451): Prop `pr:small-Mahler` (exceptional factors cyclotomic with probability
  ≥ 1 − 2exp(−cτd^{4/5})), Prop `pr:powers` (no proper power, ≥ 1 − 2exp(−c(τd)^{1/2})), Boyd–Montgomery + Lemma `lm:indi1`
  (deg Φ ≤ (C/τ)√d, ≥ 1 − exp(−Cd^{1/2}/4)). The union of the failure events has probability ≤ exp(−c′τd^{1/2}/log d).

*Proof of (2).* BV's choices (l.2463): κ = τ/100, m = ⌈30(log d)²⌉, X² = c_0τ²d/(100 m (log md)²), g = h_X. (m ≥ 6 log d,
Wielandt's bound, Dixon–Mortimer Thm 5.5.B, cited by BV at l.2453; their "(30 log d)²" for Bochert–Jordan in the same
sentence is a harmless misquote; a larger m such as ⌈(30 log d)²⌉ would change only constants.) So X ≍ τd^{1/2}/(log d)². With dτ⁴ ≥ C(log d)^{12} (size hypothesis above; BV's "τ⁴d large enough" at l.2463 understates
it), X lies in (condiX), in particular X > 100 m² d^{1/10} ≥ m d^{1/10}. The failure probability 2exp(−X/6) is
exp(−cτd^{1/2}/(log d)²) with a smaller c.
- (m-bd), l.2458 (RH-free): with probability ≥ 1 − 2exp(−X/6), |Σ_p B_P(p)^m log(p) h_X(log p) − B_m w| < ½.
- Prop 19′: Σ_p B_P(p)^m log(p) h_X(log p) = |Ω^m/G| + O(e^{−X/40}).  [replaces l.2464–2466]
- §2 with K = Q: |w − 1| = O(X² e^{−X/32}).  [replaces l.2467, which used RH for ζ_Q]
- B_m ≤ m^m ≤ 2^{m²} ≤ e^{X/100} (BV l.2467; m² ≍ (log d)⁴ ≪ X), so
  |B_m − |Ω^m/G|| ≤ ½ + B_m|1 − w| + O(e^{−X/40}) ≤ ½ + O(X² e^{X/100 − X/32}) + O(e^{−X/40}) < 1 for d ≥ D_1;
  thus |Ω^m/G| = B_m, G is m-transitive (Lemma `lm:perm`, which needs |Ω| = deg P̃ ≥ d − (C/τ)√d ≥ m), and
  G ⊇ Alt(deg P̃) because m ≥ 6 log d ≥ 6 log deg P̃ (Wielandt, as cited by BV, l.2453).
- Finally P = Φ P̃ as in part (1). ∎

The strip enters only through e^{−X/32} versus e^{−X/2}; at BV's X (≍ d^{1/2} up to logs) the loss is a constant in c.

## 5. BV's other statements

- Thm 1 (l.226) and Cor 4 (l.343, the 0/1 model: P(irreducible) = 1 − √(2/(πd)) + O(1/d)) are now unconditional; both were
  already proved unconditionally by Kogler 2609.26771 (Thm 1.1, by twin factorizations).
- Cor 3 (l.333): for every μ with finite second moment, not a point mass, and every N,
  P(P_d irreducible) = 1 − P(∃ω ∈ U_N : P_d(ω) = 0) + O_{μ,N}(d^{−N/2}). Its proof (l.2494–2521: Lemma `lm:cycl-div`,
  truncation at e^{d^{1/10}} by Chebyshev, Thm 2(1)) uses RH only through Thm 2(1). Unconditional now.
- Thm 5 (l.378): its hypothesis "no ρ with |1 − ρ| < (log d)^α/d" holds for every K by (FE)+(CT) once (log d)^α/d < 1/32;
  Thm 6 (l.419) likewise (no zero at all within 1/32 of 1). Both become unconditional for d ≥ d_0(α), and are implied by
  Thm 2′ when e^{d^{1/10}} ≥ d^{1/τ}. Thm 7 (l.438) is superseded: every ζ_K is zero-free in |1 − s| < 1/32.

## 6. Descendants

### 6.1 Littlewood polynomials (BSHKP 2308.04878, l.113, l.128)
Model: f = Σ_{i=0}^{n} a_i x^i, a_i i.i.d. uniform on {±1}. Cor 3 with μ uniform on {±1} (finite variance, ‖μ‖₂² = ½) and
N = 2, U_2 = {0, 1, −1}; f(0) = a_0 ≠ 0. So P(f irreducible) = 1 − P(f(1) = 0 or f(−1) = 0) + O(1/n).
- n even: f(±1) is a sum of n + 1 (odd) terms ±1, never 0: P(f irreducible) = 1 − O(1/n).
- n odd: f(1) and f(−1) = Σ (−1)^i a_i are each a sum of n + 1 i.i.d. signs: P(f(1) = 0) = C(n+1, (n+1)/2) 2^{−n−1}
  = √(2/(πn)) + O(n^{−3/2}), same for f(−1); P(both) = P(Σ_{i even} a_i = 0) P(Σ_{i odd} a_i = 0) = O(1/n).
  So P(f irreducible) = 1 − 2√(2/(πn)) + O(1/n).
Hence lim_n P(f irreducible) = 1 (the folklore conjecture of BSHKP's abstract; they proved limsup = 1). With Thm 2′(2),
Gal(f) ⊇ A_n with probability → 1 as well. Exhaustive counts n ≤ 20 (littlewood_output.txt): the reducible f with a root ±1 are
exactly the predicted count, and n·(other reducible)/2^n stays below 2.7.

### 6.2 Blachar–Breuillard–Kozma 2609.40123 (l.189)
Model: R = x^n + ω_{n−1}x^{n−1} + … + ω_0, independent, ω_1..ω_{n−1} ~ μ, |ω_j| ≤ H. Their Thm 1.1 (l.165, unconditional):
H_2(μ) > 12 ⇒ P(disc R is a square) ≤ C(μ)/log n. Their l.189: for arbitrary μ with H_2(μ) > 12 the Sym(n) corollary needs GRH
"from [BV, Thm 2]". Now: condition on ω_0 ≠ 0 (a new law for A_0; BV allow it). Thm 2′(2) (|ω_j| ≤ H ≤ e^{n^{1/10}},
‖μ‖₂² = e^{−H_2} < 1 − ½): with probability ≥ 1 − exp(−cn^{1/2}/(log n)²), R = Φ R̃, Gal(R̃) ⊇ Alt. x ∤ R; P(R(±1) = 0) = O(n^{−1/2})
(local limit theorem; Lemma `lm:cycl-div` with n = 1, 2); Σ_{k≥3} P(ω_k | R) = O(1/n) (BV l.2496–2513). So
P(R irreducible, Gal(R) ⊇ A_n | ω_0 ≠ 0) ≥ 1 − O(n^{−1/2}), and with Thm 1.1 (Gal ⊆ A_n ⇒ disc square):
**P(R irreducible with Gal(R) = S_n | ω_0 ≠ 0) ≥ 1 − C′(μ)/log n for every μ with H_2(μ) > 12**, unconditionally.

### 6.3 Varjú–Xu 2511.04240 (random multiplicative ±1 coefficients)
RH enters only via Prop `pr:prime-number-theorem` (l.384–394, = BV Prop 9) inside Prop `pr:number-distinct-factors`
(l.361–375, proof l.416–458), in §6 (l.972–1031), and, for Cor 1.3 only, GS's Prop 9.1 under RH for Dirichlet L (l.1038–1050).
(a) Prop `pr:number-distinct-factors` with (9′): each K_Q has n_{K_Q} ≤ d and log Δ_{K_Q} ≤ log((2d)^{4d} M^{2d}) (l.420–424);
≤ d factors; so the error is O(d² X² log(dM) e^{−X/32}) (the q | Δ term C d² log(dM) e^{−X} is unchanged).
(b) §6: replace X = 5 log d (l.976) by X = 81 log d and l > 5/ε (l.974) by l > 81/ε. Prop `pr:number-roots-random-poly`
(l.821–830) needs only q > d, l ≥ ε^{−1} log q/log d (log q ≤ X = 81 log d) and d, q large in terms of l. Then (eq:weights):
Σ_q log(q) h_X(log q) = 1 + O(X² d^{−81/32}); E[#distinct irreducible factors of R] = 1 + O(d^{−1/2+ε}) +
O(d² X² log d · d^{−81/32}) = 1 + O(d^{−1/2+ε}) since 2 − 81/32 = −0.53; the lower bound "≥ ½ for every realization" and the
proper-power step (l.1004–1027) go through verbatim. **Thm 1.1: P(P_d(x)/x irreducible) ≥ 1 − C(ε)d^{−1/2+ε}, unconditionally.**
(c) Cor 1.3 (Legendre-symbol model F_{d,f}) needs a count of primes p ≤ x with (r/p) = ω_r for all primes r ≤ d.
Write the indicator as 2^{−π(d)} Σ_{m | Π_{r≤d} r} ω_m (m/p) (p > d). For squarefree m > 1, (m/p) = χ_D(p) for odd p ∤ m,
D the discriminant of Q(√m) (|D| ≤ 4m), χ_D primitive. For K = Q(√D), Λ_K = Λ·(1 + χ_D) exactly (split, inert, ramified),
so ψ(x, χ_D) = ψ_K(x) − ψ(x), and Prop 9′ (dyadic sum) gives |ψ(x, χ_D)| ≤ C x^{31/32}(1 + log x)(log|D| + 1 + log x).
With log m ≤ θ(d) ≤ 2d, θ − ψ = O(√x log x) and the primes p ≤ d (O(d log d)):
  Σ_{p ≤ x, (r/p) = ω_r ∀ r ≤ d} log p = x 2^{−π(d)} + O(x^{31/32}(d + log x)²),
GS Prop 9.1 (l.1038–1050) with x^{31/32} for x^{1/2}. VX's proof (l.1051–1090) needs the error ≤ C f(d)2^{−π(d)} at x = f(d):
true when f(d)^{1/32} ≥ C 2^{π(d)}(d + log f(d))², e.g. f(d) ≥ 2^{32π(d)} d^{65} for d ≥ d_0 (then f^{1/32} = 2^{π(d)}d^{2.03}
and d + log f ≤ 3d). With the Lean-checked half-plane for
ζ and every Dirichlet L-function (Bebel, T* = B_2 + 2.2·10^{−254} = 0.8733819666…, TREE.md N7) and Remark 1.2 instead:
f(d) ≥ 2^{a π(d)} with any a > 1/(1 − T*) = 7.8978 suffices for d ≥ d_0(a). (GRH gives VX's f(d) > 2^{2π(d)}d⁴.)
So Cor 1.3 is unconditional in the range f(d) > 2^{7.9π(d)} (d large); the GRH range 2^{2π(d)}d⁴ needs strip width ½.

### 6.4 Bary-Soroker–Garzoni–Sodin 2502.17218 (random tridiagonal matrices)
RH enters only via (eq:bv) (l.936–947, "from BV Prop 9 and Lemma 21"), used at l.1149, l.1241 and l.2043.
(a) (eq:bv′). For P ∈ Z[x] of degree ≤ n, height H, k ≥ 1, x ≥ 2: for p ∤ Δ_{P̃}, #(Z(P,p))_k = Σ_{O ∈ (Z(P))_k/Gal} A_{K_O}(p)
(Frobenius count as in BV Prop 16, distinct tuples). Prop 9′ on (x, 2x] for each K_O (≤ n^k orbits, n_{K_O} ≤ n^k,
log Δ_{K_O} ≤ 2k n^k log(Hn)), plus ≤ 2n log(Hn)/log x primes p | Δ_{P̃} contributing ≤ n^k log(2x)/x each:
  |x^{−1} Σ_{x<p≤2x} log(p) #(Z(P,p))_k − #((Z(P))_k/Gal)| ≤ C k n^{2k} log(Hn) log²(2x) x^{−1/32}.
(b) Stretched-exponential version (l.1145–1160): x = exp((n/k⁸)^{1/(2C)}), k ≤ n^c, log H ≤ exp(n^ε) with probability
≥ 1 − C_ε exp(−(b/2)n^ε) (Cor `cor:height`). The right side of (eq:bv′) is ≤ exp(−n^c) once n^ε + 2k log n ≪ (log x)/32,
which holds for c, ε small, exactly as with x^{−1/2}.
(c) Exponential version (l.1239–1245): x = exp(c′n) with c′ ≤ 1/C_k (Lemma `cor:mix-upgr` needs C_k log x ≤ n ≤ c_k x^{1/(2k)});
on log H ≤ e^{εn} (probability ≥ 1 − C_ε e^{−(b/2)εn}) with ε < c′/64, (eq:bv′) ≤ C k n^{2k}(e^{εn} + log n)(c′n)² e^{−c′n/32}
≤ e^{−c′n/64}. So #orbits = 1 for k = 1 and k = 6 with probability ≥ 1 − C e^{−cn}: **Theorem 1 (l.357: irreducible, Gal ∈ {S_n, A_n},
probability ≥ 1 − Ce^{−cn}) and Theorem 2 (l.821: Gal ≥ U_m ⋊ A_m) unconditional**, and the minors corollary they draw (l.685).
(d) Scope slip: their (eq:rh), l.318–322, assumes RH only for Q(a), a a single root, but (eq:bv) for k ≥ 2 needs ζ_{K_O},
K_O = Q(a_1, …, a_k) (BV's Thm 2(2) assumes RH for these). Moot now.

### 6.5 Further ERH-conditional descendants found during the search
- Eberhard 2008.01223, Thm 1.3 (label `main`, l.154–171; 2304.12075 cites it as Theorem 1.3): under ERH, P(φ irreducible) = 1 − O(e^{−cn}), P(Gal(φ) ≥ A_n) = 1 − O(e^{−cn^{1/2}})
  for φ = det(t − M), M with i.i.d. entries of law μ. RH enters via BV Prop 9 (l.605–612) in Prop `local-global` (l.693–700):
  error O(m n^{2m+1} X² e^{−X/2} log(Hn)), which becomes O(m n^{2m+1} X² e^{−X/32} log(Hn)) by (9′) (log Δ_{Q(ω)} ≤
  m n^{m−1} log|Δ_φ̃| and |Δ_φ̃| ≤ H^{n²} n^{n²}). m = 1 (l.722–743): E R_φ(p) = 1 + O(1/p + p e^{−cn}) for all p > 2H + 1, so
  with X = 32cn/33 (for cn/2) both errors are O(n⁵ log(Hn) e^{−cn/33}): P(irreducible) = 1 − O(e^{−c″n}) (plus their
  RH-free Prop `one-simple-eig`). m ≥ 2 (l.788–842): X = c n^{1/2}/m^{5/2} (l.811) is unchanged; the extra factor needs
  X/32 ≫ (2m + 1) log n, true for fixed m = 6. **Eberhard's Thm 1.3 holds unconditionally, same rates, smaller c.**
- FJSS 2106.04049, Thm 1.1 (l.157–160): under ERH, the characteristic polynomial of a random symmetric ±1 matrix is irreducible
  with probability ≥ 1 − 2exp(−cn^{1/4}). RH enters via BV Prop 9 (l.481–484) and Eb Prop 3.5 (l.486–489) at l.514–521, with
  X = c n^{1/4} (l.518, the top of the local range p ≤ exp(cn^{1/4})): error n³X² e^{−X/32} log n = exp(−c n^{1/4}/32 + O(log n)).
  **Unconditional with probability ≥ 1 − 2exp(−c′n^{1/4}).**
- 2304.12075 (eigenvector minors): Cor 1.2 (l.145) rests on Eberhard's Thm 1.3 and becomes unconditional with it;
  Cor 1.3 (l.149, G(n,p) adjacency matrices) rests on FJSS (l.237); we did not check that model transfer.

## 7. How much strip is needed (barrier bookkeeping)

With Remark 1.1 (strip width δ), the applications need δX ≥ (exponent of the polynomial prefactor) + O(log X + log(1/δ)):
| use | X used | prefactor | δ that suffices |
|---|---|---|---|
| BV Thm 2(1), Prop 19′ m = 1 | ≍ τd^{1/2}/log d | d^{2.1}X² | δ ≥ C(log d)²/(τd^{1/2}) |
| BV Thm 2(2), m ≍ (log d)² | ≍ τd^{1/2}/(log d)² | d^{2m}X² (Prop 19′); B_m ≤ m^m (the w term, ζ only) | δ ≥ C(log d)⁵/(τd^{1/2}) |
| VX Thm 1.1 | free: X = (2.6/δ) log d | d²X² log d vs d^{−1/2} | any fixed δ > 0 (l > 2.6/(δε)) |
| VX Cor 1.3 | x = f(d) | 2^{π(d)} | range f ≥ 2^{π(d)/δ}d^{O(1/δ)}; GRH range needs δ = ½ |
| BGS, Eberhard, FJSS | x = e^{c′n}, e^{cn}, e^{cn^{1/4}} | n^{O(1)}e^{εn} etc. | any fixed δ > 0 |
In Thm 2(2) the term B_m|1 − w| needs δ_ζ X ≥ log B_m + O(log X) for the strip width δ_ζ of ζ alone, i.e.
δ_ζ ≥ C(log d)⁴ log log d/(τd^{1/2}) (B_m ≤ m^m). So every use tolerates δ → 0 polynomially in d; δ = 1/32 is far more
than enough once d ≥ D_0 (params_output.txt: delta_min columns). The classical region δ ≍ 1/log Δ_K ≍ 1/(d^{1.1}) does not suffice anywhere, which is
why BV needed RH (their l.548–551). No step fails with δ = 1/32.

## 8. Not re-checked here
BV's RH-free parts (random walks §5–6, Konyagin's estimates, Mahler measure §7, proper powers §8) and the RH-free parts of
the descendants (equidistribution, Green–Tao input in VX, SL_2 expansion in BGS, local problems in Eb and FJSS) are used as
published; we only verified that they contain no hypothesis on zeros. BV Thm 10's contour-shift details are left to the
reader there; it is the classical explicit formula.

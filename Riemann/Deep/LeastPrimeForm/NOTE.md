# Least prime represented by a binary quadratic form: an explicit unconditional exponent? (2026-10-20)

Brief: `.scratch/orca/lpqf-brief.md`. Notation, the inputs (F1)–(F6), Lemma 2 (weighted large sieve) and
Theorem 1 (ring class zero density) are those of [Combos NOTE](../Combos/NOTE.md). Here L := log|D|, a zero
is written ρ = β + iγ = 1 − λ/L + iγ, and N(λ) is the number of characters χ of the ring class group whose
L(s,χ*) has a zero with λ_ρ ≤ λ and |γ| ≤ 1. Nothing below uses the 2026 half-planes (L), (O), (Hu).

## Verdicts

| Question | Answer |
|---|---|
| Explicit unconditional A below 455 for **all** D, from Theorem 1 + zero-free region + Deuring–Heilbronn + log-free density | **No** (§1). The exponent is fixed by an explicit count of zeros with λ ≈ 0.6; Theorem 1 gives no information there. |
| Is Theorem 1 log-free near σ = 1? | No. At σ = 1 − λ/L it gives C(ε)e^{2λ}\|D\|^ε (§1.2). A log-free version with an inexplicit constant would not give an explicit A either. |
| Side result: density-one D | **New candidate** (Theorem B, §2): for every A > A* = 1.72419, all fundamental D < 0 outside a density-zero set have p(Q) ≤ \|D\|^A for every form Q of discriminant D. Previous: \|D\|^{20/3+ε} (Ditchen 2013). |

## 0. Records (checked 2026-10-20, §3)

- All D, unconditional: |D|^455 (Zaman, PhD thesis, Toronto 2017, Thm 7.4.1 = D_K^455 Q^342 and Cor 1.4.1);
  |D|^694 (Thorner–Zaman, ANT 11 (2017), arXiv 1604.01750). Exceptional case only (an η-Siegel zero with
  η ≥ η(δ)): |D|^{6+δ} or |D|^{9.5+δ} (Zaman, Cor 1.4.2). Inexplicit: Fogels, Weiss, Kowalski–Michel 2002
  Cor 7, Gaudet arXiv 2408.01610 (prime D, sieve). Thorner–Zhang arXiv 2412.01802 improve Weiss's bound
  group-theoretically (abelian H replaced by AHC subgroups), with unspecified exponents; for forms H is
  already abelian, so nothing changes.
- Density-one D (all classes): |D|^{20/3+ε} (Ditchen arXiv 1312.1502, Cor 1.3, D ≢ 0 mod 8); |D|^{3+ε}
  for density-one D and almost all classes. Gaudet (arXiv 2408.01610, §1) still quotes 20/3 as the record. Zaman
  arXiv 1710.08914 treats upper bounds only. Sardari arXiv 1802.06193, 1803.03218: GRH, or a positive
  proportion of classes.

## 1. Why Theorem 1 gives no explicit exponent for all D

### 1.1 What fixes the exponent
In the explicit-formula method (Combos NOTE (1)), hS(x) = xW(1) − Σ_χ χ̄(C)Σ_ρ x^ρW(ρ) + O(h log x), and one
must show Σ_{χ,ρ} x^{β−1}|W(ρ)| < W(1). With x = |D|^A a zero with parameter λ contributes ≍ e^{−Aλ}. If
K zeros with parameter ≈ λ₀ cannot be excluded, the termwise bound only certifies A ≥ (log K)/λ₀. So an
explicit A needs an **explicit** bound for N(λ) at every λ where the method is not yet free of zeros.

Zaman's 455 is exactly this quotient. In the proof of Thm 7.4.1 (case λ₁ ≥ 10^{−2}):
- positivity (explicit inequality of Heath-Brown 1992 §12; Zaman Cor 4.5.3) gives N(0.569) ≤ 3365;
- beyond 0.569 only the power-sum bound is available (Zaman Table 5.1): log N(0.6) ≤ 258.0;
- with B' = 453.5 the first term of his Z₁ is e^{−453.5·0.569}e^{258.0} ≈ 0.96, and his total is
  Z₁ ≤ 0.9649. So B' ≥ 258.0/0.569 = 453.4 is forced; hence his B = 454 and the range x ≥ D_K^455 Q^342.
The tail λ ≥ 1 (log-free density e^{162λ+188}) contributes Z₂ ≤ e^{−100}: it is not binding. In the
sub-case 10^{−3} ≤ λ₁ < 10^{−2} the binding term is again a low-lying count (e^{162λ+188} at λ = 0.968).

### 1.2 Theorem 1 near σ = 1
At σ = 1 − λ/L, Theorem 1(a) reads N_D(σ;t) ≤ C(ε)e^{2λ(1−λ/L)}|D|^ε T³. For fixed λ this grows with |D|,
so it says nothing about N(λ). The loss |D|^ε is built into the proof: the mollifier has length
X = |D|^{ε₀}, its coefficients are bounded only by d(𝔞), the Sobolev step costs log²M, the dyadic
choice costs J^{2k}, and Lemma 4 loses |D|^ε.

A log-free version N(1 − λ/L) ≤ C₀e^{cλ} with an unspecified C₀ (Kowalski–Michel, Pacific J. Math. 207
(2002), Thm 2 gives this for the family, any c > 5.5) also yields no explicit A: by §1.1, K = C₀e^{cλ₀}
unknown means A unknown. This is why Kowalski–Michel's Cor 7 has an unspecified exponent.

Theorem 1 does improve the part of Zaman's argument that handles zeros with λ ≥ εL, but that part
is not binding (§1.1).

### 1.3 The sharper large sieve does not reach the power-sum bounds
Zaman, Thm 5.3.3: Σ_χ N(σ,T,χ) ≤ (4ξ/√(ξ²−1))(C₄λ⁴+C₃λ³+C₁λ+C₀)e^{B₁λ+B₂}(1−J(ξλ))^{−2}, with
B₁ = 4φAξ log(4eα^{−1}(1+α)^{2(1+α)/α}) and B₂ = 16 log(4eα^{−1}(1+α)^{2(1+α)/α}). B₁, B₂ come from the
Kolesnik–Straus power-sum inequality and do not involve the large sieve. The large sieve (his Thm 5.2.2)
enters only through the polynomial prefactor C₄ (its factor 5π and the υ-terms) and through the side condition
(5.10), y ≥ C h_H D_K^{3/2+υ}Q^{1/2}T^{…}. For a ring class family of fundamental D this needs y = e^{YL}
with Y ≥ 2 + υ; the power-sum parameters already force Y_λ ≈ 5.3–5.8 (Table 5.1). So replacing his large
sieve by Lemma 2 (constant N + |D|^{1/2}) changes log N(0.6) = 258 by O(1), and 455 by a few units at
most.

### 1.4 Positivity is capped by the conductor term
Zaman Thm 4.5.1 and Prop 4.5.2 (after Heath-Brown §12) count characters with a zero in λ ≤ λ₀ using
L^{−1}K(s,χ;f) ≤ f(0)φL_χ/L − F(λ̃ − b) + ε, φ = 1/4 (convexity). They need F(λ − b) > (4/3)f(0)φ; the
off-diagonal terms add the N² condition. The term f(0)φL_χ/L is the conductor term of the explicit
inequality; it needs f(0) > 0, i.e. it comes from the small primes, where no large sieve acts. A
nonnegative weight vanishing near 0 cannot have Re F ≥ 0 on Re z ≥ 0 (a positive-definite function
peaks at 0; so Condition 2 even forces f(0) = max f), and Lemma 2 cannot remove the cap. This is why positivity stops at λ = 0.569.

### 1.5 The mollifier route needs Möbius cancellation in Q(√D) below √|D|
Heath-Brown's (and Graham's) mollifier method gives the good constants for Dirichlet characters. For the
ring class family the zero detector Σ_𝔞 Δ(𝔞)χ(𝔞)N𝔞^{−ρ}, Δ = 1 ∗ λ, must have Δ(𝔞) = 0 for
1 < N𝔞 ≤ w with w ≥ |D|^{1/2}: below √|D| the class-group large sieve constant is ≍ h, not N (Combos NOTE,
proof of Lemma 2). Möbius inversion then forces λ_𝔡 = μ_K(𝔡) for N𝔡 ≤ w, and Graham's mean value of Δ²
needs cancellation in Σ_{N𝔞≤u} μ_K(𝔞)N𝔞^{−1} at u ≈ |D|^{1/2}, i.e. control of 1/ζ_K = 1/(ζ·L(·,χ_{d_K}))
near s = 1, including residues 1/L'(ρ). No d_K-uniform bound of that kind is known (Zaman §5.1.1). An
integer mollifier with Selberg's pseudo-characters (Kowalski–Michel §3–4) avoids this but only with
unspecified constants.

### 1.6 Conclusion
No explicit unconditional A below 455 follows. What would: any explicit N(λ) for λ ∈ [0.569, 1] better than
Zaman's Table 5.1; since A ≈ log N(0.6)/0.569, halving log N(0.6) roughly halves A.

## 2. Theorem B: density-one discriminants

### Inputs
- (T1) Theorem 1 of the Combos NOTE (for fundamental D: f = 1, χ* = χ).
- (KM) Kowalski–Michel, Pacific J. Math. 207 (2002), Thm 2: for a family S(q) of cuspidal automorphic
  representations of GL(n)/ℚ that satisfy Ramanujan, have conductor ≤ q^A, |S(q)| ≪ q^d and a common
  component at ∞: Σ_{π∈S(q)} N(π;α,T) ≪ T^b q^{c₀(1−α)/(2α−1)} (some b depending on the family) for α ≥ 3/4, T ≥ 2, any c₀ > 5nA/2 + d.
  Applied to:
  - S₂(q) = {θ_χ : D fundamental, −q ≤ D < 0, χ ∈ Cl(D)^, χ² ≠ 1}. θ_χ is a weight-one cuspidal newform of
    level |D| with L(s,θ_χ) = L(s,χ) (Hecke; Iwaniec, Topics in classical automorphic forms, Ch. 12);
    θ_χ = θ_ψ iff ψ ∈ {χ, χ̄}. Here n = 2, A = 1, d = 3/2, so any c₀ > 13/2.
  - S₁^±(q) = {χ_d : d fundamental, 1 < |d| ≤ q, χ_d(−1) = ±1}: n = 1, A = 1, d = 1, any c₁ > 7/2.
- (Gen) For fundamental D the real characters of Cl(D) are the χ_{d₁,d₂}, D = d₁d₂ with d₁, d₂ fundamental
  (d₁ = 1 allowed), and L(s,χ_{d₁,d₂}) = L(s,χ_{d₁})L(s,χ_{d₂}), L(s,χ₁) := ζ(s) (classical; see Ditchen
  1312.1502 §7). Each fundamental d | D with D/d fundamental occurs in exactly one such pair.
- (ZF) Dirichlet L-functions (Davenport, Multiplicative Number Theory, Ch. 14): L(s,χ_d) ≠ 0 for
  σ > 1 − c/log(|d|(|t|+2)) except possibly one simple real zero β_d; effectively 1 − β_d ≥ c/(|d|^{1/2}log²|d|);
  Landau–Page: for Q ≥ 3 at most one real primitive χ of conductor ≤ Q has a real zero > 1 − c_LP/log Q.
- (W) A weight w ∈ C_c^∞((1/2,1)), w ≥ 0, w ≢ 0, of Gevrey class 2 (e.g. a rescaled exp(−1/(1−u²))). Then
  |W(σ+it)| ≤ C exp(−c|t|^{1/2}) for −1 ≤ σ ≤ 2 (integrate by parts k times in log v, |(v d/dv)^k w| ≤ C^{k+1}(k!)²,
  take k ≍ |t|^{1/2}). (F6) holds for this w.

### Statement
**Theorem B.** Let A* = 1.724189… be the root in (1,2) of A³ − 3A² − (13/2)A + 15 = 0. For every A > A*
there is a set V_A of negative fundamental discriminants with #{D ∈ V_A : |D| ≤ X} = o(X) such that for
every negative fundamental D ∉ V_A with |D| ≥ D₀(A), every class C ∈ Cl(D) and every x ≥ |D|^A, the number
of prime ideals 𝔭 ∈ C of degree one with x/2 < N𝔭 ≤ x is ≥ c_w x/(h(D) log x). In particular every
primitive positive definite form of discriminant D represents a prime p ∤ D with p ≤ |D|^A.

### Proof
*Parameters.* Fix A ∈ (A*, 1.8) (larger A follow). The function g(A,c) := (1 − A/2)(c/(A−1) − A) is
decreasing in A, increasing in c, and g(A*, 13/2) = 1. Choose c₀ > 13/2, η₁ > 0 and ε' > 0 so small that
φ(σ) := c₀(1−σ)/(2σ−1) − A(1−σ) satisfies φ(σ₀) ≤ 1 − 3ε' at σ₀ := A/2 − η₁ (≥ 3/4). φ is decreasing on
[σ₀,1), since φ'(σ) = A − c₀(2σ−1)^{−2} < 0. Fix ε ∈ (0, η₁(1 − A/2)) for Theorem 1, c₁ = 3.6, η₂ = ε'/A and
the grid σ_j = σ₀ + jη₂ (0 ≤ j ≤ J, J minimal with σ_J ≥ 1 − ε'/(2A)). Let b₁, b₂ be the T-exponents of
(KM) for S₁^± and S₂. For X a power of 2 put T₁ = (log X)³, Y₁ = log log X, Y₀ = (log X)^κ with
κ = (3b₁ + 1)/0.55, and treat the D with X/2 < |D| ≤ X. Write
N^cusp_D(σ) := #{(χ,ρ) : χ² ≠ 1, L(ρ,χ) = 0, β ≥ σ, |γ| ≤ T₁}.

*Step 1: the exceptional set.* Call D (with X/2 < |D| ≤ X) good if:
- (G1) N^cusp_D(σ_j) ≤ |D|^{A(1−σ_j)−2ε'} for 0 ≤ j ≤ J (for j = J this says N^cusp_D(σ_J) = 0);
- (G2) every fundamental d | D with D/d fundamental and |d| ≥ Y₀ (this includes d = D) has L(s,χ_d) ≠ 0
  for β ≥ 9/10, |γ| ≤ T₁;
- (G3) no such d with Y₁ ≤ |d| < Y₀ lies in E := {d : L(s,χ_d) has a real zero > 1 − c_LP/(2 log|d|)};
- (G4) ω(D) ≤ 10 log log X.
The bad D in (X/2, X] number o(X):
- (G1): since A(1−σ_j) − 2ε' ≥ −2ε', failure at j implies N^cusp_D(σ_j) > 2^{−2ε'}X^{A(1−σ_j)−2ε'}. Each
  θ ∈ S₂(X) comes from at most two χ, so Markov and (KM) give
  #{D : (G1) fails at j} ≤ 2^{1+2ε'}X^{2ε'−A(1−σ_j)}Σ_{π∈S₂(X)}N(π;σ_j,T₁) ≪ T₁^{b₂}X^{φ(σ_j)+2ε'} ≤ (log X)^{3b₂}X^{1−ε'}.
  There are J + 1 = O(1) values of j.
- (G2): by (KM) for S₁^± with α = 9/10, #{d : |d| ≤ Z, L(s,χ_d) has a zero with β ≥ 9/10, |γ| ≤ T₁}
  ≪ T₁^{b₁}Z^{c₁/8} = T₁^{b₁}Z^{0.45}. A given d divides at most X/|d| + 1 discriminants D in (X/2, X], so
  #D failing (G2) ≪ X(log X)^{3b₁}Y₀^{−0.55} + (log X)^{3b₁}X^{0.45} = O(X/log X).
- (G3): if Z ≤ |d| < |d'| ≤ Z² and d, d' ∈ E, both real zeros exceed 1 − c_LP/log(Z²), against Landau–Page
  with Q = Z². So E has at most one element in each [Z, Z²], Σ_{d∈E, |d|≥Y₁}1/|d| ≤ 2/Y₁, and
  #D failing (G3) ≤ 2X/Y₁ + O(log log X) = o(X).
- (G4): #{n ≤ X : ω(n) > 10 log log X} ≤ 2^{−10 log log X}Σ_{n≤X}2^{ω(n)} ≪ X(log X)^{−5}.
Let V_A be the union over X of the bad D in (X/2, X]. Then #{D ∈ V_A : |D| ≤ Z} ≤ Σ_{X≤2Z} o(X) = o(Z).

*Step 2: the zero sum for good D.* Let D be good, X/2 < |D| ≤ X, x = |D|^A y with y ≥ 1. The weight
(W) gives |W(ρ)| ≤ C e^{−c|γ|^{1/2}}. By (G4) Cl(D) has at most 2^{ω(D)} ≤ (log X)^7 real characters.
- All zeros with |γ| > T₁: Theorem 1(b) gives Σ ≤ Σ_{|k|≥T₁−1}Ce^{−c|k|^{1/2}}·x·C|D|^{1/2}log³(3|D|)log(|k|+3)
  ≤ x|D|^{−1}.
The zeros with |γ| ≤ T₁ are of four kinds.
- (a) Zeros of the factors L(s,χ_d), |d| < Y₀ (including ζ), of real characters. For |d| < Y₁, (ZF) gives
  1 − β ≥ c/(Y₁^{1/2}log²(Y₁+2)) or 1 − β ≥ c/log(Y₁(T₁+3)); for Y₁ ≤ |d| < Y₀, (G3) and (ZF) give
  1 − β ≥ c/log(Y₀(T₁+3)). In all cases 1 − β ≥ c_κ/log log X. There are ≤ (log X)^{11} such zeros and each
  has x^{β} ≤ x·exp(−c_κ A log|D|/log log X). Total: o(x).
- (b) Zeros of L(s,χ), χ² ≠ 1, with β ≥ σ₀. By (G1), at most |D|^{A(1−σ_j)−2ε'} of them have β ∈ [σ_j, σ_{j+1}),
  and these have x^{β−1} ≤ |D|^{−A(1−σ_{j+1})}; none has β ≥ σ_J. So
  Σ x^{β−1}|W(ρ)| ≤ CΣ_{j<J}|D|^{Aη₂−2ε'} = CJ|D|^{−ε'}.
- (c) Zeros of L(s,χ_d), |d| ≥ Y₀, with β ≥ σ₀. By (G2) β < 9/10, so x^{β−1} ≤ |D|^{−A/10}, and there are
  ≤ (log X)^{11} of them. Total ≤ x|D|^{−0.17}.
- (d) All remaining zeros have β < σ₀. Use x^β ≤ x^{1/2} + log x∫_{1/2}^{β}x^σdσ and count with Theorem 1:
  Σ ≤ CΣ_{|k|≤T₁+1}[x^{1/2}N_D(0;k) + log x∫_{1/2}^{σ₀}x^σN_D(σ;k)dσ]. For σ ≤ 0.6, Theorem 1(b) gives
  ≤ x|D|^{−0.18}y^{−0.4}log T₁, as in the proof of Theorem 2 of the Combos NOTE (A ≥ 1.7). For
  0.6 ≤ σ ≤ σ₀, Theorem 1(a) gives x^σN_D(σ;k) ≤ C(ε)x·y^{−(1−σ)}|D|^{−(1−σ)(A−2σ)+ε}T₁³. The quadratic
  (1−σ)(A−2σ) is convex with vertex (A+2)/4 > σ₀, so on [0.6, σ₀] it is ≥ (1−σ₀)(A−2σ₀) ≥ 2η₁(1 − A/2) > 2ε.
  Keeping y^{−(1−σ)} ≤ y^{−0.1} against log x, the total is ≤ C(ε)x|D|^{−ε}(log X)^{14}.
So Σ_{χ,ρ} x^β|W(ρ)| = o(x), uniformly in x ≥ |D|^A.

*Step 3: conclusion.* For fundamental D, S_χ = S_{χ*}, and (F6) gives
hS_C(x) = xW(1) − Σ_χ χ̄(C)Σ_ρ x^ρW(ρ) + O(h log x), with h ≤ 3|D|^{1/2}log²(3|D|). By Step 2,
S_C(x) = (xW(1)/h)(1 + o(1)). The terms from ramified primes, prime powers and inert primes are
≤ (xW(1)/h)|D|^{−0.2} for x ≥ |D|^{1.7}, as in the proof of Theorem 2 of the Combos NOTE. Each p lies
under at most two degree-one primes of C, which gives the count. A degree-one prime 𝔭 ∈ C of norm p
gives a representation of p by the forms Q with [𝔞_Q] = C^{−1} ((F1), (F2)). □

### Remarks
1. Simpler variant: requiring N^cusp_D(1−δ) = 0 with c₀δ/(1−2δ) < 1, i.e. δ < 2/17, gives the exponent
   2(1−δ) + ε, so any A > 30/17 = 1.7647.
2. Where Theorem 1 enters: it controls every zero left of A/2 for every single D; the family estimate (KM),
   averaged over D, takes over right of A/2 through Markov's inequality. With Pascadi–Thorner's density
   |D|^{3(1−σ)/(3−2σ)} (fundamental D) in place of Theorem 1 the same argument gives only A > 2.32. With
   Kowalski–Michel's estimates alone (per D, exponent 5.5, valid for σ ≥ 3/4, plus Markov over D) it gives
   A > 9, worse than Ditchen. Neither bound is stated anywhere found.
3. A* depends only on the family exponent c₀ = 13/2 of (KM): c₀ = 5.5, 4, 3 would give 1.685, 1.603, 1.524. A
   sharper log-free density for {θ_χ : |D| ≤ X} (for instance from Ditchen's or Duke's large sieve over D)
   would lower A*. GRH gives 1 + ε for every D.
4. The exceptional set is o(X); its worst parts are (G3) (2X/log log X) and (G2) (X/log X). The implied
   constants of (KM) are not explicit, so D₀(A) is not explicit; this does not matter for a density
   statement.
5. Status: proof complete modulo Theorem 1 (ours, unrefereed, reviewed) and the published (KM), (Gen), (ZF).
   Not refereed, not in Lean.

## 3. Novelty audit (2026-10-20)

- arXiv API (`.scratch/lpqf/ax*`, plus the 105-query sweep of `.scratch/combosv/ax/all.json`):
  - abs:"binary quadratic form" AND abs:"least prime": Thorner–Zaman 1803.02823, 1604.01750; Ditchen 1312.1502.
  - abs:"fundamental discriminants" AND abs:"least prime": Ditchen only.
  - abs:"almost all" AND abs:"binary quadratic" AND abs:prime: Zaman 1710.08914 (upper bounds), Luca–Radziwiłł–Shparlinski.
  - abs:"Bombieri-Vinogradov" AND abs:"quadratic form": Lam 1812.08823 (fixed form), Thorner 1401.6677, Ditchen.
  - abs:"primes represented" AND abs:"almost all": Zaman 1710.08914, Teräväinen 1611.08585.
  - abs:"Deuring-Heilbronn": Benli–Goel–Twiss–Zaman 2410.06082 (Dirichlet, explicit DH), Thorner–Zaman 2108.10878, Kadiri–Wong 2105.14181, Zaman 1508.00287, 1502.05679.
  - abs:Weiss AND abs:Chebotarev: Thorner–Zhang 2412.01802 (read: inexplicit, group-theoretic).
  - abs:"least prime" AND abs:Hecke, abs:Linnik AND abs:Hecke, abs:"split prime" AND abs:"ideal class",
    abs:"least prime" AND abs:"quadratic field", abs:"density one" AND abs:discriminants AND abs:prime,
    abs:"Bombieri-Vinogradov" AND abs:"class group", abs:"ideal class" AND abs:"least prime": nothing relevant.
- GXL evidence API, 2 claims: (i) density-one D with exponent below 20/3: only Ditchen (supports 3+ε for most
  forms, 20/3 for all); (ii) explicit all-D exponent below 455: only Thorner–Zaman (694) and citing papers.
- Read: Zaman's thesis (Ch. 1, 4.5, 5.1–5.3, 7.4), Kowalski–Michel 2002 (§1–4), Ditchen 1312.1502 (§1, §7),
  Zaman 1710.08914 (§1), Sardari 1802.06193, 1803.03218 (abstracts, §1), Pascadi–Thorner 2508.14888 (grep).

## 4. Dead ends (NOGO candidates)
- Theorem 1 + Linnik → explicit A for all D: Theorem 1 is |D|^ε-lossy at λ = O(1); explicit A needs N(0.6) explicitly.
- Our large sieve in Zaman's power-sum method: threshold not binding (Y ≈ 5.5 vs 2); B₁, B₂ are Kolesnik–Straus constants.
- Positivity past λ ≈ 0.57: capped by the conductor term f(0)φL_χ/L; large sieve cannot remove it.
- Graham/Heath-Brown mollifier for the family: needs μ_{Q(√D)} cancellation below √|D|.
- Inexplicit log-free density (Kowalski–Michel type) → explicit A: impossible in principle (§1.1).

## 5. Leads
- Fundamental D: Zaman's positivity bounds use L_χ/L ≤ 4/3 (his (3.3)–(3.4) with Q > 1); for Q = 1 the ratio
  is 1. Recomputing his Table 4.14 and Cor 4.5.3 with 1 could push positivity past 0.569 and lower 455 for
  fundamental D. Pure bookkeeping on Zaman's method, not attempted.
- Theorem B with non-fundamental D (old forms; needs Theorem 1 for orders, which holds) and with an explicit
  exceptional-set bound.

# AllDisc: the A > 1 least-prime theorem for all discriminants (every quadratic order, both signs)

Date 2026-10-22. Brief: `.scratch/orca/alldisc-brief.md`. Check scripts and outputs: [checks/](checks/) (§9); the scratch directory
`.scratch/alldisc/` (arXiv listings, abstracts) was deleted. References: [C] = [Combos NOTE](../Combos/NOTE.md),
[D2] = [Density2 NOTE](../Density2/NOTE.md), [R] = [RCApps NOTE](../RCApps/NOTE.md), [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md),
[R3] = [RCApps3 NOTE](../RCApps3/NOTE.md), [R4] = [RCApps4 NOTE](../RCApps4/NOTE.md), [R5] = [RCApps5 NOTE](../RCApps5/NOTE.md).
Paper proofs, not refereed, no Lean. Unconditional (no 2026 half-plane, no Kerr). The inputs are those of [R5]. New here: Lemma M
(ring class conductors and multiplicities), Proposition 1^o (zeros near 1 with multiplicities), Lemma G (real characters of orders),
and the bookkeeping of §2–§6.

## Verdicts

| Item (brief) | Result | Status |
|---|---|---|
| (1) target | Theorems 3^o, 3⁺^o: for every A > 1, all but O(X^{1−c}) non-square discriminants Δ ≡ 0, 1 (mod 4) with \|Δ\| ≤ X (fundamental or not, either sign) satisfy: every ring class of O_Δ (narrow if Δ > 0) has ≥ ν_Δ H/(16h(Δ) log 2x) degree-one primes with norm in (x, x+H], for all x ≥ \|Δ\|^A and H ∈ [x\|Δ\|^{−κ}, x]. So every primitive form of discriminant Δ represents a prime in each such interval, and a prime ≤ 2\|Δ\|^A | **new**, complete proof. [R5] had fundamental D only, about 61% of all discriminants |
| (2) family sieve | Theorem 1^o: [R5] Theorems 1, 1⁺ over all pairs (Δ, χ ∈ Ĝ_Δ), all orders: ≪ X^{3/2+14δ}N^{1−2σ}. No new analysis is needed: Lemma P already sums over all orders, and the left side only gains the non-fundamental Δ | new statement |
| (2) family density | Theorem 2^o: Σ_{Δ∈𝒟(X)}N^♮_Δ(σ;t) ≪ X^{2−σ+60δ} for ½+3δ ≤ σ ≤ 9/10, both signs | new |
| (2) conductors, induced characters, family count | Lemma M: every χ ∈ Ĝ_Δ is inflated from a unique primitive pair (Δ′, χ′), with Δ = m²Δ′. Each primitive pair lies under ≤ 1 + (X/\|Δ′\|)^{1/2} elements of 𝒟(X); each GL₂ form comes from ≤ 6 primitive pairs; Σ_{Δ∈𝒟(X)}h(Δ) ≤ 5X^{3/2}log²(3X). Proposition 1^o: Σ_Δ N^c_Δ(σ,T) ≪ (XT)^ε T(X^{10(1−σ)/(3−2σ)} + X^{1/2}) | new; the extra X^{1/2} comes from small conductors with large multiplicity and is harmless for Markov |
| (3) zero detector | [R3] Lemmas 3, 4 hold for every χ ∈ Ĝ_Δ∖{1}: λ_χ(p) = 0 = (Δ/p) for p \| f, and L(s,χ) = L(s,χ*)E_χ(s) with E_χ zero-free in Re s > 0. The class I coefficients stay independent of Δ (Lemma D) | routine |
| (4) genus characters, Siegel zeros | Lemma G: the real characters of G_Δ correspond to disjoint pairs {d, d′} with L(s,χ*) = L(s,χ_d)L(s,χ_{d′}) and dd′ \| Δ. E₂ counts the Δ with d \| Δ; ν_Δ is defined as in [R]. Only the χ_d can have Siegel zeros | routine (genus theory of orders, proved for both signs) |
| (5) Δ > 0 | narrow classes (P⁺_{K,Z}(f)), the unit η_Δ and the regulator R⁺_Δ of each order enter only through Lemma P⁺ (per order), h⁺R⁺ = √Δ L_Δ (per order) and h⁺ ≤ 5Δ^{1/2}log²(3Δ) | routine (§6) |
| corollaries | for all but O(N^{1−c}) integers n ≤ N (no squarefree or congruence condition): x² + ny² has a prime in every (x, x+x(4n)^{−κ}], x ≥ (4n)^A. The same for x² − ny² and all but O(N^{1−c}) non-square n | new for the n not covered by [R5] (non-squarefree n, n ≡ 3 mod 4); Baier–Zhao 2007 average over squarefree n only |
| (6) failure points | none. Every use of fundamentality in [R]–[R5] has a replacement (table in §8) | — |
| not covered | sectors for non-maximal orders: Lemma P^∠ already covers orders, but Theorem 2^∠ for all orders and [R] Proposition 1′ with multiplicities are not written out. A < 1: NOGO ([R5] §7) | open / NOGO |

## 0. Notation

- 𝒟⁻ (resp. 𝒟⁺): the negative (resp. positive non-square) discriminants Δ ≡ 0, 1 (mod 4). Write Δ = f²d_K, K = Q(√Δ), O_Δ = Z + fO_K.
  𝒟(X) := {Δ ∈ 𝒟⁻ : X/2 < |Δ| ≤ X} and 𝒟⁺(X) := {Δ ∈ 𝒟⁺ : X/2 < Δ ≤ X}. [R5]'s 𝔇(X) ⊂ 𝒟(X) and 𝔇⁺(X) ⊂ 𝒟⁺(X).
  #𝒟(X) = X/4 + O(1).
- G_Δ := I_K(f)/P_{K,Z}(f) ≅ Cl(O_Δ) if Δ < 0, and G_Δ := I_K(f)/P⁺_{K,Z}(f) ≅ Cl⁺(O_Δ) if Δ > 0 ([C] (F1), [D2] §0, [R4] §1).
  h(Δ) := |G_Δ|. For every order of either sign, h(Δ) ≤ 5|Δ|^{1/2}log²(3|Δ|). This follows from [D2] §0: h = w|Δ|^{1/2}L_Δ/2π,
  h⁺log η_Δ = Δ^{1/2}L_Δ, log η_Δ ≥ 0.962 and L_Δ ≤ 4log²(3|Δ|).
- For χ ∈ Ĝ_Δ: χ(𝔞) := χ([𝔞]) if (𝔞, f) = 1 and 0 otherwise; L(s,χ) := Σ_{(𝔞,f)=1}χ(𝔞)N𝔞^{−s}. χ* is the primitive Hecke character
  inducing χ, with conductor 𝔣_χ | fO_K (times ∞ if Δ > 0) and A(χ*) := |d_K|N𝔣_χ ≤ |Δ|. Then L(s,χ) = L(s,χ*)E_χ(s) with
  E_χ(s) := Π_{𝔭|f, 𝔭∤𝔣_χ}(1 − χ*(𝔭)N𝔭^{−s}). E_χ has no zeros in Re s > 0, and |E_χ| ≤ 4^{ω(f)} on Re s ≥ 0 ([C] §0).
- λ_χ(n) := χ(T_n) with T_n := Σ_{N𝔞=n,(𝔞,f)=1}[𝔞] ∈ Z[G_Δ] ([R3] §0). λ_χ(n) = 0 if (n, f) > 1, because an ideal prime to f has
  norm prime to f.
- Ĝ_{Δ,2} := {χ : χ² = 1} and Ĝ_{Δ,c} := Ĝ_Δ∖Ĝ_{Δ,2}. Zero counts, with multiplicity:
  - N^♮_Δ(σ;t) := Σ_{χ∈Ĝ_Δ∖{1}}#{ρ : L(ρ,χ*) = 0, σ ≤ β ≤ 9/10, |γ − t| ≤ 1};
  - N^c_Δ(σ,T) := Σ_{χ∈Ĝ_{Δ,c}}#{ρ : L(ρ,χ*) = 0, β ≥ σ, |γ| ≤ T}.
  For fundamental Δ these are the counts of [R3] and [R].
- F₂(σ) := 10(1−σ)/(3−2σ). Admissible n, Φ, a_V, b_{j,t′}, ‖b‖² := N^{1−2σ}, δ and the parameters of [R5] Theorem 3 are as in [R3] §0, §2
  and [R5] §5.

## 1. Ring class conductors, induced characters, multiplicities, family count (task 2)

**Lemma M.** Let K be quadratic. Write P(m) for P_{K,Z}(m) if K is imaginary and for P⁺_{K,Z}(m) if K is real.
(a) Let f₁, f₂ | f and g := (f₁, f₂). Then P(g) ∩ I_K(f) ⊂ (P(f₁) ∩ I_K(f))·(P(f₂) ∩ I_K(f)).
(b) For χ ∈ Ĝ_Δ (Δ = f²d_K) and f′ | f, say that χ *factors through f′* if χ is trivial on P(f′) ∩ I_K(f). Then χ factors through f′
  iff χ is inflated from a character of G_{f′²d_K}. The divisors f′ | f through which χ factors are exactly the multiples of one
  divisor f_χ. Hence χ is inflated from a unique χ′ ∈ Ĝ_{Δ_χ}, where Δ_χ := f_χ²d_K, and χ′ is *primitive*: it factors through no
  proper divisor of f_χ. Moreover χ′* = χ*, and χ′² = 1 iff χ² = 1.
(c) Conversely, let (Δ′, χ′) be primitive with Δ′ = f′²d_K, let m ≥ 1, and let χ be the inflation of χ′ to G_{m²Δ′}. Then
  f_χ = f′ and its primitive pair is (Δ′, χ′). So (Δ, χ) ↦ (Δ_χ, χ′, f/f_χ) is a bijection from {(Δ, χ) : Δ ∈ 𝒟⁻, χ ∈ Ĝ_Δ} onto
  {primitive pairs} × Z_{≥1} (likewise for 𝒟⁺). For each primitive pair, #{Δ ∈ 𝒟(X) over (Δ′, χ′)} = #{m : X/2 < m²|Δ′| ≤ X}
  ≤ 1 + (X/|Δ′|)^{1/2}.
(d) For a cuspidal automorphic representation π of GL₂(A_Q), at most 6 primitive pairs (Δ′, χ′), with Δ′ of either sign and χ′
  complex, have π_{χ′*} ≅ π. Here π_{χ*} is the automorphic induction of χ* (the theta series of [R] Proposition 1, or the Maass
  form of [R2] Proposition 2).
(e) Family count: Σ_{Δ∈𝒟(X)}h(Δ) ≤ 5X^{3/2}log²(3X), and the same for 𝒟⁺(X). The number of primitive pairs with |Δ′| ≤ Y is
  ≤ 5Y^{3/2}log²(3Y).

*Proof.* (a) Every element of P(g) ∩ I_K(f) is (α) with α ∈ K^× prime to f and α ≡ a mod^× gO_K, where a ∈ Z and (a, g) = 1;
in the real case also α ≫ 0. For each prime p | f₁f₂ put e_i := v_p(f_i) and e := max(e₁, e₂). By the Chinese remainder theorem in
O_K choose β ∈ O_K with:
- β ≡ α mod^× p^e O_K at the p with e₁ ≤ e₂;
- β ≡ 1 mod p^e O_K at the p with e₁ > e₂;
- β ≡ 1 mod pO_K at the p | f with p ∤ f₁f₂.

In the real case replace β by β + tM, with M a large power of f and t ∈ Z large, so that β ≫ 0. Then (β) ∈ P(f₁) ∩ I_K(f):
- at p with e₁ ≤ e₂ we have v_p(g) = e₁ and β ≡ α ≡ a mod^× p^{e₁};
- at p with e₁ > e₂ we have β ≡ 1 mod p^{e₁}.

So β ≡ a₁ mod^× f₁O_K, where a₁ ∈ Z is prime to f₁ (CRT). Also β₂ := α/β satisfies (β₂) ∈ P(f₂) ∩ I_K(f):
- at p with e₁ ≤ e₂, β₂ ≡ 1 mod^× p^{e₂};
- at p with e₁ > e₂, v_p(g) = e₂ and β₂ ≡ α ≡ a mod^× p^{e₂};
- β₂ ≫ 0 if α, β ≫ 0.

Hence (α) = (β)(β₂).
(b) For f′ | f the natural map I_K(f)/(P(f′) ∩ I_K(f)) → I_K(f′)/P(f′) = G_{f′²d_K} is injective, and it is surjective because every
class of G_{f′²d_K} contains ideals prime to any given integer ([C] (F1); Cox, Cor 7.17; the same proof for P⁺). This gives the first
claim. By (a), if χ factors through f₁ and f₂ then it factors through (f₁, f₂). The set of such f′ contains f, so it has a least
element f_χ, which divides all the others. Every multiple f″ | f of f_χ also works, since P(f″) ⊂ P(f_χ) for f_χ | f″. The other
claims are immediate: χ′ and χ agree on ideals prime to f, so χ′* = χ*, and inflation is injective.
(c) χ factors through f′ by construction, so f_χ | f′. Let 𝔞 ∈ P(f_χ) ∩ I_K(f′). By the surjectivity in (b), its class in
G_{f′²d_K} contains some 𝔟 = 𝔞𝔠 ∈ I_K(f) with 𝔠 ∈ P(f′) ⊂ P(f_χ). Then 𝔟 ∈ P(f_χ) ∩ I_K(f), so χ′(𝔞) = χ(𝔟) = 1. Thus χ′
factors through f_χ, and primitivity gives f_χ = f′. The count: m runs over an interval (√(X/(2|Δ′|)), √(X/|Δ′|)].
(d) If π_{χ′*} ≅ π_{ξ′*}, the Artin representations ρ = Ind_K^Q χ′* and ρ′ = Ind_{K′}^Q ξ′* (irreducible, since the characters are
complex) have the same Frobenius traces at almost all p, so ρ ≅ ρ′ (Chebotarev and Brauer–Nesbitt). Fix ρ. If ρ ≅ Ind_K ψ, then
ρ ⊗ χ_{d_K} ≅ ρ, so χ_{d_K} is a nontrivial 1-dimensional constituent of ρ ⊗ ρ^∨ = 1 ⊕ ad⁰ρ. Since dim ad⁰ρ = 3, there are at most
3 fields K. For fixed K, Res_K ρ = ψ ⊕ ψ^c, so ψ ∈ {χ′*, χ̄′*}: 2 choices. Finally ψ = χ′* determines the primitive pair. Suppose
two primitive pairs (f₁²d_K, χ₁), (f₂²d_K, χ₂) have the same χ*. Then χ*|_{I_K(f₁f₂)} factors through f₁ and through f₂, hence
through (f₁, f₂) by (a), and primitivity forces f₁ = (f₁, f₂) = f₂. Then χ₁ = χ₂, since both equal χ* on I_K(f₁). Total ≤ 3·2 = 6.
(e) There are ≤ X discriminants in the range, each with h ≤ 5X^{1/2}log²(3X) (§0). Each primitive pair is a pair (Δ′, χ′) with
χ′ ∈ Ĝ_{Δ′}. □

*Size of the family.* The fundamental Δ alone give Σ h ≫ X^{3/2} for Δ < 0. The non-fundamental Δ are a fraction 1 − 6/π² ≈ 0.39
of 𝒟(X) (§9: 3884 of 9901 for 200 ≤ |Δ| ≤ 20000, Δ < 0). The non-fundamental part of the family sum of Theorem 1^o is about one third of the
total (§9). The analytic conductors are A(χ*) ≤ |Δ_χ| ≤ |Δ|, with equality only for primitive χ with 𝔣_χ = fO_K. Theorems 1^o and 2^o
use only the upper bound |Δ| (the level of the imprimitive L(s,χ)). The large multiplicity of the small conductors matters only near
σ = 1, in Proposition 1^o.

**Proposition 1^o (zeros near the 1-line, all orders).** For X, T ≥ 3, 0 ≤ σ ≤ 1 and ε > 0:

  Σ_{Δ∈𝒟(X)}N^c_Δ(σ,T) ≤ C(ε)(XT)^ε T·(X^{F₂(σ)} + X^{1/2}),

and the same for 𝒟⁺(X).

*Proof.* By Lemma M(b), (c),

  Σ_{Δ∈𝒟(X)}N^c_Δ(σ,T) = Σ_{(Δ′,χ′) primitive, complex}#{m : m²Δ′ ∈ 𝒟(X)}·N(σ,T;χ′*).

Here N(σ,T;ψ) counts the zeros of L(s,ψ) with β ≥ σ and |γ| ≤ T. Split the primitive pairs by |Δ′| ∈ (Y/2, Y] with Y = X/2^k,
0 ≤ k ≤ log₂X + 1. On such a block the multiplicity is ≤ 1 + 2^{(k+1)/2}, and by Lemma M(d) the block contributes
≤ 6Σ_{π∈𝒮_Y}N_π(σ,T), where 𝒮_Y := {π_{χ′*} : (Δ′, χ′) primitive complex, Y/2 < |Δ′| ≤ Y} is a set of distinct cuspidal
representations. The proof of [R] Proposition 1 applies to 𝒮_Y. Each π ∈ 𝒮_Y is cuspidal and satisfies GRC (unitary Satake
parameters), and C(π) ≤ 12Y. By Lemma M(e), |𝒮_Y| ≤ 5Y^{3/2}log²(3Y). Finally I₂ := Σ_{π∈𝒮_Y}∫_{−T}^{T}|L^{ur}(½+iv,π)|²dv is
≤ C(ε)Y^{3/2+2ε}T^{2+ε}, for two reasons:
- |L^{ur}(½+iv,π_{χ*})| ≤ 4^{ω(A(χ*))}|L(½+iv,χ*)|, since at most two local factors are removed at each p | A(χ*), each with
  |L_p^{−1}| ≤ (1 + p^{−1/2})² ≤ 4;
- Σ_{χ∈Ĝ_{Δ′}∖{1}}|L(½+iv,χ*)|² ≤ C(ε)|Δ′|^{1/2+ε}(|v|+3)^{1+ε} for every order. This is [C] Lemma 4 with mollifier length 1,
  whose steps 2–5 bound Σ_χ|L(½+iu,χ*)M_X|² for D = f²d_K. For Δ′ > 0 it is [D2] (2.1) with 𝒵 = Ĝ⁺.

So Σ_{π∈𝒮_Y}N_π(σ,T) ≤ C(ε)(YT)^{F₂(σ)+ε}. Summing over k, 2^{k/2}(X/2^k)^{F₂} = X^{F₂}2^{k(1/2−F₂)} ≤ max{X^{F₂}, X^{1/2}}
for 0 ≤ k ≤ log₂X + 1, and there are ≤ log₂X + 2 values of k. For Δ′ > 0 use the inputs of [R2] Proposition 2: π_{χ*} is a Maass
form of eigenvalue ¼, tempered with GRC, and C(π) ≤ 12Y. Rename ε. □

## 2. Theorem 1^o: the family large sieve over all pairs (Δ, χ ∈ Ĝ_Δ)

**Theorem 1^o.** Take the hypotheses of [R5] Theorem 1: 0 < δ ≤ 1/100, X ≥ X₀(δ), ½+2δ ≤ σ ≤ 9/10, |t′| ≤ X^δ, 2 ≤ V ≤ X^{1/2},
N = 2^j ∈ [V/4, X^{1−4δ}], and b ∈ {b_{j,t′}, ∂_{t′}b_{j,t′}}. Then

  Σ_{Δ∈𝒟(X)}Σ_{χ∈Ĝ_Δ}|Σ_n b(n)λ_χ(n)|² ≤ C(δ)X^{3/2+14δ}N^{1−2σ},

and the same with 𝒟⁺(X) and the narrow groups Ĝ_Δ.

*Proof.* Δ < 0. For each Δ ∈ 𝒟(X) (|Δ| > 4), Parseval on G_Δ gives Σ_χ|Σ_n b(n)λ_χ(n)|² = h(Δ)Σ_{C∈G_Δ}|S_C|², where
S_C := Σ_{𝔞∈C,(𝔞,f)=1}b(N𝔞) ([R3] Lemma 2(i), valid for every order). Since ideals prime to f have norms prime to f,
S_C = S_C^{(Δ)}(b₁) with b₁ := b·1[(·, f) = 1]. So Σ_C|S_C|² is the block e = 1 of Lemma P ([R5] §2), with weight w₁/2 = 1.
All blocks are ≥ 0, hence Σ_C|S_C|² ≤ Q′_Δ(b). Also Φ(|Δ|/X) = 1 on 𝒟(X), Φ ≥ 0, Q′_Δ(b) ≥ 0 for every Δ ∈ 𝒟⁻, and
h(Δ) ≤ 5X^{1/2}log²(3X). Hence the left side is ≤ 5X^{1/2}log²(3X)·𝒬′(b), with 𝒬′(b) = Σ_{Δ∈𝒟⁻}Φ(|Δ|/X)Q′_Δ(b) the quantity of
[R5] Corollary P. The proof of [R5] Theorem 1 shows |𝒬′(b)| ≤ ½Σ_c|Q_c| ≤ C(δ)X^{1+13δ}‖b‖².

Δ > 0. [R4] Lemma 2⁺(i), (ii) hold for every non-square Δ: the order enters only through G⁺_Δ, η_Δ, R⁺_Δ and Dirichlet's formula
h⁺R⁺ = √Δ·L(1,χ_{d_K})Π_{p|f}(1 − (d_K/p)/p). So Σ_χ|Σ_n b(n)λ_χ(n)|² ≤ ŵ(0)^{−1}h⁺(Δ)R⁺_ΔQ_Δ(b) ≤ CX^{1/2}log²X·Q_Δ(b).
By Lemma P⁺ ([R5] §6) the block e = 1 of Q′_Δ(b) is Q_Δ(b) and the other blocks are ≥ 0, so Q_Δ(b) ≤ Q′_Δ(b). The sum
Σ_ΔΦ(Δ/X)Q′_Δ(b) over all non-square Δ > 0 is the one bounded in the proof of [R5] Theorem 1⁺. □

*Remarks.* (1) The inequality chain is that of [R5]; only the left side changes. [R5] dropped the terms of non-fundamental Δ by
positivity, and here they are kept. (2) For a character χ inflated from (Δ′, χ′), λ_χ(n) = λ_{χ′}(n)·1[(n, f) = 1]. So each primitive
character is counted once for every order above it, and positivity in the dual absorbs this multiplicity. Numerically the
non-fundamental part is about one third of the total, and the total stays at the family size from N = X/16 to N = 1.5X (§9).

## 3. The zero detector over all ring class characters; Theorem 2^o

**Lemma D (detector, all orders).** Let Δ ∈ 𝒟(X) ∪ 𝒟⁺(X) with |Δ| ≥ D₀(δ), and χ ∈ Ĝ_Δ∖{1}. Put
𝒢(s,χ) := Σ_{n adm.}a_V(n)λ_χ(n)n^{−s} and 𝒢^F := F𝒢 with F(s) = 1 − 4^{1−s}. Then [R3] Lemma 3 holds with |D| replaced by |Δ|,
and [R3] Lemma 4 holds with ρ any zero of L(s,χ*) with σ ≤ β ≤ 9/10.

*Proof.* The changes are these.
(1) *Euler factors.* For every prime p put L_p(s,χ) := (1 − λ_χ(p)p^{−s} + (Δ/p)p^{−2s})^{−1}. Then L(s,χ) = Π_pL_p(s,χ) for
Re s > 1, as follows. T_{ab} = T_aT_b for coprime a, b. For p ∤ f, χ((p)) = 1 because (p) ∈ P_{K,Z}(f) (P⁺: p ≫ 0), so the local
factor is that of [R3] Lemma 1 (split, inert, ramified), with (Δ/p) = (d_K/p). For p | f there is no ideal prime to f above p, and
λ_χ(p) = 0 = (Δ/p), so L_p = 1. [R3] Lemma 3(a) uses only these identities and |λ_χ(p)| ≤ 2.
(2) *Zeros.* L(s,χ) = L(s,χ*)E_χ(s) with E_χ entire and zero-free in Re s > 0. So 𝒢 = L(s,χ)·E·H·M̃ vanishes at every zero of
L(s,χ*) with real part > ½. Also χ ≠ 1 implies χ* ≠ 1, so L(s,χ) is entire of polynomial growth in vertical strips, which the
contour shift of [R3] Lemma 4 needs.
(3) *Lemma 3(c).* M_t(χ) is a sum over ideals of squarefree norm prime to f; these are primitive ideals. [C] Lemma 2 holds for
every order with constant (N + |Δ|^{1/2})log²(3|Δ|) (for Δ > 0 use [D2] Lemma 5), and V ≤ |Δ|^{1/2}.
(4) *The class II mean value* ([R3] §4): Σ_{χ∈Ĝ_Δ∖{1}}|L(½+δ+iv,χ)|² ≤ C(δ)|Δ|^{1/2+δ}(|v|+3)^{1+δ}. This is [D2] (2.1)
(Proposition 6 holds for every order D = f²d_K, and for 𝒵 = Ĝ⁺ of real orders), used at ½ + δ as in [R3] §4, together with
|E_χ| ≤ 4^{ω(f)} ≤ C(δ)|Δ|^{δ/2}. □

**Theorem 2^o.** Let 0 < δ ≤ 1/100, ½+3δ ≤ σ ≤ 9/10 and T := |t| + 3 ≤ X^δ/4. Then

  Σ_{Δ∈𝒟(X)}N^♮_Δ(σ;t) ≤ C(δ)X^{2−σ+60δ},

and the same for 𝒟⁺(X).

*Proof.* Use [R5] Theorem 2, i.e. [R3] §4 with Y := X^{1−5δ}, V := X^{1/2}/2 and Z := 2Y log²Y, for each Δ ∈ 𝒟(X):
- V ≤ |Δ|^{1/2};
- each χ ∈ Ĝ_Δ∖{1} has ≤ C₀log(XT²) zeros in the window ([C] (F4) for χ*, A(χ*) ≤ X);
- Lemma D splits the χ with a zero in [σ, 9/10] × [t−1, t+1] into 𝒳_I(Δ) ∪ 𝒳_II(Δ).

Class II: |𝒳_II(Δ)| ≤ C(δ)X^{1/2+2δ}Y^{1/2+δ−σ}T by Lemma D(3), (4) and Cauchy–Schwarz over all χ ∈ Ĝ_Δ∖{1}. Summed over the
≤ X discriminants this is ≤ C(δ)X^{2−σ+9δ}. Class I: the coefficients b_{j,t′} do not depend on Δ. So for |τ| ≤ X^δ/2 the sum over
Δ ∈ 𝒟(X) of Σ_{χ∈Ĝ_Δ}sup|P_j(χ,t′)|² is bounded by Theorem 1^o (Sobolev, b_{j,t′} and ∂_{t′}b_{j,t′}). For |τ| > X^δ/2 use the
per-order bound Σ_ΔΣ_χ|P_j|² ≤ CX^{1+δ}(N + X^{1/2})N^{1−2σ} ([C] Lemmas 2, 3, each order). As in [R5],
Σ_Δ|𝒳_I(Δ)| ≤ C(δ)X^{3/2+17δ}V^{1−2σ} ≤ C(δ)X^{2−σ+18δ}. Multiply by C₀log(XT²). For Δ > 0 apply [R4] Theorem 2⁺ with the same
changes and Theorem 1^o for Δ > 0. □

## 4. Genus characters of orders; Siegel zeros (task 4)

**Lemma G (real characters of an order).** Let Δ = f²d_K ∈ 𝒟⁻ ∪ 𝒟⁺ and let L_Δ be its ring class field (narrow if Δ > 0), the class
field of P(f) ⊂ I_K(f). For χ ∈ Ĝ_{Δ,2}, χ ≠ 1, let L_χ ⊂ L_Δ be the fixed field of ker χ (via the Artin map G_Δ ≅ Gal(L_Δ/K)).
(i) L_χ/Q is Galois with group (Z/2)². So L_χ = Q(√d_K, √d) for a fundamental discriminant d ∉ {1, d_K}. If d′ is the fundamental
  discriminant of Q(√(dd_K)), then L(s,χ*) = L(s,χ_d)L(s,χ_{d′}) and |dd′| = A(χ*).
(ii) dd′ divides Δ. In particular d | Δ and |d|, |d′| ≤ |Δ|.
(iii) Put P(χ) := {d, d′} and P(1) := {1, d_K}. The sets P(χ), χ ∈ Ĝ_{Δ,2}, are pairwise disjoint. So
  𝔉(Δ) := ∪_χP(χ) has 2|Ĝ_{Δ,2}| elements: fundamental discriminants (or 1) dividing Δ, at most 4·2^{ω(Δ)} of them.
(iv) Let w ≥ 0, w ≤ 1, be supported in (x, 2x) (x ≥ 2). Put ψ_w(χ_d) := Σ_nΛ(n)χ_d(n)w(n) and
  S_ξ(w) := Σ_𝔞Λ_K(𝔞)ξ(𝔞)w(N𝔞), for ξ a character of I_K(f) or a Hecke character. Then for every C ∈ G_Δ,

    Σ_{χ∈Ĝ_Δ}χ̄(C)S_χ(w) = Σ_{d∈𝔉(Δ)}ε_d(C)ψ_w(χ_d) + Σ_{χ∈Ĝ_{Δ,c}}χ̄(C)S_{χ*}(w) + R_f(C),  |R_f(C)| ≤ 2h(Δ)ω(f)log 2x,

  with ε_d(C) := χ(C) ∈ {±1} for d ∈ P(χ). In particular ε_1 = ε_{d_K} = 1.

*Proof.* (i) *Galois.* P(f) and I_K(f) are stable under the nontrivial automorphism σ of K. So for every extension σ̃ of σ,
σ̃(L_Δ) is the class field of σ(P(f)) = P(f), i.e. σ̃(L_Δ) = L_Δ, and L_Δ/Q is Galois. For τ ∈ Gal(L_Δ/Q) with τ|_K = σ and
𝔞 ∈ I_K(f), functoriality of the Artin map gives τ·Art(𝔞)·τ^{−1} = Art(𝔞̄) = Art(𝔞)^{−1}. Indeed 𝔞𝔞̄ = (N𝔞), and N𝔞 is a
positive integer prime to f, so (N𝔞) ∈ P(f) (also in P⁺). This is Cox, Lemma 9.3, for Δ < 0. Hence τ(ker χ)τ^{−1} = ker χ, so ker χ
is normal in Gal(L_Δ/Q) and L_χ/Q is Galois of degree 4.
*Not cyclic.* Take a prime p ∤ Δ inert in K. The prime ideal pO_K lies in P(f), so its Artin symbol in Gal(L_Δ/K) is trivial and
pO_K splits completely in L_χ/K. So the decomposition group of p in Gal(L_χ/Q) has order 2 (p is unramified with residue degree 2),
and it is not contained in Gal(L_χ/K) because its generator restricts to σ. A cyclic group of order 4 has only one subgroup of order
2, which would have to be Gal(L_χ/K). Hence Gal(L_χ/Q) ≅ (Z/2)².
*Fields.* Its three quadratic subfields are K, Q(√d), Q(√d′), and the three nontrivial characters satisfy χ_{d_K}χ_dχ_{d′} = 1.
*L-functions.* By class field theory χ* is the Hecke character of the quadratic extension L_χ/K, so L(s,χ*) = ζ_{L_χ}(s)/ζ_K(s),
with all Euler factors. Moreover ζ_{L_χ} = ζ·L(s,χ_{d_K})L(s,χ_d)L(s,χ_{d′}) (abelian field) and ζ_K = ζ·L(s,χ_{d_K}).
*Conductors.* By the conductor–discriminant formula over Q and over K (𝔡_{L_χ/K} = 𝔣_χ for a quadratic extension),
|disc L_χ| = |d_K dd′| = d_K²N𝔣_χ. Hence |dd′| = |d_K|N𝔣_χ = A(χ*).
(ii) 𝔣_χ | fO_K gives N𝔣_χ | f², so |dd′| = |d_K|N𝔣_χ divides |d_K|f² = |Δ|.
(iii) Suppose d ∈ P(χ) ∩ P(χ′) with χ ≠ χ′. If both are nontrivial, then L_χ = K(√d) = L_{χ′}, so ker χ = ker χ′ and χ = χ′ (a
character of order 2 is determined by its kernel): a contradiction. If χ′ = 1, then d ∈ {1, d_K}, which (i) excludes. Also d ≠ d′,
since d_K is not a square. For the size: by (ii) every d ∈ 𝔉(Δ) is ±u·2^a with u odd, squarefree and composed of primes dividing Δ,
and a ∈ {0, 2, 3}; the sign is forced for a ∈ {0, 2}. That gives at most 2^{ω(Δ)}(1 + 1 + 2) values.
(iv) Λ_K is supported on prime powers, so S_{χ*}(w) − S_χ(w) = Σ_{𝔭|f,k≥1}Λ_K(𝔭^k)χ*(𝔭^k)w(N𝔭^k). Each of the ≤ 2ω(f) primes
𝔭 | f has at most one power N𝔭^k in (x, 2x), so |S_χ − S_{χ*}| ≤ 2ω(f)log 2x for every χ ∈ Ĝ_Δ. For real χ ≠ 1, comparing the
coefficients of −L′/L in (i) gives S_{χ*}(w) = ψ_w(χ_d) + ψ_w(χ_{d′}). For χ = 1, χ* = 1_K and ζ_K = ζ·L(s,χ_{d_K}) give
S_{1_K}(w) = ψ_w(1) + ψ_w(χ_{d_K}). Now sum with weights χ̄(C) = χ(C) and use (iii). □

So the classical facts used in [R] (genus theory of Cox §6, with L(s,χ_{d,D/d}) = L(s,χ_d)L(s,χ_{D/d})) hold for every order, with
𝔉(Δ) in place of 𝔉(D). For fundamental D, 𝔉(D) is [R]'s set: the pair of χ_{d,D/d} is {d, D/d}.

*Siegel zeros.* Only real characters can have exceptional zeros, and by Lemma G(i) their L-functions are products of two Dirichlet
L-functions L(s,χ_d), d ∈ 𝔉(Δ). For complex χ the zeros near 1 are counted in N^c_Δ and removed by Markov (E₁^o below); no special
treatment is needed. (I4) and (I5) of [R] concern Dirichlet characters and are unchanged.

*Exceptional set E₂^o.* For X a power of 2 and Q₀(X) := (X/2)^δ, let E₂^o(X) be the set of Δ ∈ 𝒟(X) (or 𝒟⁺(X)) such that some
d ∈ 𝔉(Δ) with |d| > Q₀(X) has L(ρ,χ_d) = 0 with β > σ_g = 0.995 and |γ| ≤ T_X. Count:
- for Z ≥ Q₀(X), at most 10^{91}Z^{0.495} fundamental d with |d| ∈ (Z, 2Z] have such a zero ([R] §2.2, via N(σ_g, 2Z) of (I4),
  since T_X ≤ 2Z);
- by Lemma G(ii), d ∈ 𝔉(Δ) only if d | Δ, which happens for ≤ X/(2|d|) + 1 ≤ X/Z + 1 elements of 𝒟(X); and |d| ≤ X.

Summing over dyadic Z ∈ [Q₀(X), X]:

  #E₂^o(X) ≤ C(X·Q₀(X)^{−0.505} + X^{0.495}log X) ≤ CX^{1−0.505δ}.

*ν_Δ.* Let ν_Δ := min{1, (1 − β₁(|Δ|^δ))·δ log|Δ|} if the exceptional character χ₁ of (I4) with Q = |Δ|^δ exists and equals χ_d for
some d ∈ 𝔉(Δ), and ν_Δ := 1 otherwise. Then |d| ≤ |Δ|^δ and (I5) give ν_Δ ≥ c|Δ|^{−δ/2}log^{−1}|Δ|. For a non-maximal order the field
discriminant d_K ∈ 𝔉(Δ) can be small (|d_K| ≤ |Δ|^δ). Then χ_{d_K} belongs to 𝔉_s, and if χ₁ = χ_{d_K}, the main term
W(1) − ε_{d_K}(C)W(β₁) = W(1) − W(β₁) is still ≥ ν_ΔW(1)/2, exactly as in [R] §2.5.

## 5. Theorem 3^o: least primes in every class, all Δ < 0, every A > 1

**Theorem 3^o.** For every A > 1 there are effective c, κ, δ > 0 and X₀ (for 1 < A ≤ 1.2 one may take δ = (A−1)/8000, κ = δ/4,
c = δ/2) with the following property. For X ≥ X₀, all Δ ∈ 𝒟⁻ with |Δ| ≤ X outside a set of size ≪_A X^{1−c} satisfy, for every
C ∈ G_Δ, every x ≥ |Δ|^A and every H ∈ [x|Δ|^{−κ}, x],

  #{𝔭 ∈ C prime of degree one : x < N𝔭 ≤ x + H} ≥ ν_Δ H/(16h(Δ)log 2x).   (5.1)

In particular every primitive positive definite binary quadratic form of discriminant Δ represents a prime p ∤ Δ in every such
interval, and a prime p ≤ 2|Δ|^A.

*Proof.* The statement for A implies it for every A′ > A (same exceptional set), so let 1 < A ≤ 1.2 and take the parameters of [R5]
step 1′: ε′ = (A−1)/80, δ = ε′/100, κ = δ/4, ε₀ = ε₁ = δ/100, σ₀ = ½+4δ, σ_j, J, T_X, Q₀ and σ_g = 0.995. For X a power of 2 let
E₁^o(X) be the set of Δ ∈ 𝒟(X) for which some j has N^c_Δ(σ_j,T_X) > X^{A(1−σ_j)−2ε′}, and E^o(X) := E₁^o(X) ∪ E₂^o(X).

*Step 2′ (Markov).* For σ ∈ [σ₀, 1), Σ_{Δ∈𝒟(X)}N^c_Δ(σ,T_X) ≤ C(δ)X^{F^o(σ)+61δ}, where:
- F^o(σ) := 2 − σ for σ ≤ 9/10. The zeros with β ≤ 9/10 are counted by Theorem 2^o, summed over ≤ 2T_X unit windows; those with
  β > 9/10 by Proposition 1^o at 9/10, which gives C(XT_X)^{ε}T_X(X^{5/6} + X^{1/2}), and 5/6 < 2 − σ.
- F^o(σ) := max{F₂(σ), ½} for σ > 9/10, by Proposition 1^o.

Then φ(σ) := F^o(σ) − A(1−σ) ≤ 1 − 8ε′ on [σ₀, 1): on [σ₀, 9/10] as in [R5], and on (9/10, 1) because φ ≤ max{F₂(σ), ½} ≤ 5/6.
Markov as in [R] §2.2 gives #E₁^o(X) ≤ C(δ)J X^{1−ε′}; #E₂^o(X) ≤ CX^{1−0.505δ} (§4). Take the union over dyadic X and both
sets: c := δ/2 ≤ min{ε′, 0.505δ}.

*The explicit formula.* Fix Δ ∈ 𝒟(X)∖E^o(X) with |Δ| ≥ D₀, C ∈ G_Δ, x ≥ |Δ|^A and H ∈ [x|Δ|^{−κ}, x]. Take w and S_C as in
[R] §2.3 (supp w ⊂ (x, x+H) ⊂ (x, 2x)), Q₀ := |Δ|^δ and T₁ := 2|Δ|^{κ+ε₁} ≤ T_X. By orthogonality on G_Δ and Lemma G(iv),

  h(Δ)S_C = Σ_{d∈𝔉(Δ)}ε_d(C)ψ_w(χ_d) + Σ_{χ∈Ĝ_{Δ,c}}χ̄(C)S_{χ*}(w) + R_f(C).   (5.2)

This is [R] (2.2), with 𝔉(Δ) for 𝔉(D), χ* for χ, and the extra term R_f. Apply (I6) to every ψ_w(χ_d) (conductor |d| ≤ |Δ|) and
every S_{χ*}(w) (conductor A(χ*) ≤ |Δ|). Split 𝔉(Δ) = 𝔉_s ⊔ 𝔉_l (|d| ≤ Q₀ or not; 1 ∈ 𝔉_s) and define M(C) and E_a, …, E_f as in
[R] §2.3. By Lemma G(iii) every primitive character χ_d occurs exactly once in (5.2), so M(C) := W(1) − ε_{d₀}(C)W(β₁) is well
defined. The bounds:
- (a) Δ ∉ E₁^o; as in [R] §2.4(a).
- (b) [R5] step 3: each unit window holds ≤ C|Δ|^{1/2}log³(3|Δ|)log(|k|+3) zeros of all χ* ([C] Theorem 1(b), which uses only (F3),
  (F4) and holds for every order). So E_b ≤ H|Δ|^{−(A−1)/2+6δ}.
- (c) As in [R]: each window holds ≤ Ch(Δ)log(|Δ|(|k|+3)) zeros of the χ*, plus ≤ |𝔉(Δ)|C log(|Δ|(|k|+3)) zeros of the χ_d.
- (d) Δ ∉ E₂^o. There are ≤ |𝔉(Δ)|·CT₁log|Δ| such zeros, all with β ≤ σ_g, and |𝔉(Δ)| ≤ 4·2^{ω(Δ)} = |Δ|^{o(1)}. So
  E_d ≤ H|Δ|^{−0.005A+κ+ε₁+o(1)}.
- (e) As in [R], with [R3] step 4: the χ_d, d ∈ 𝔉_s, are distinct primitive characters of conductor ≤ Q₀, so their zeros are zeros of
  ℒ(s, Q₀).
- (f) Σ|R| ≤ C(h(Δ) + |𝔉(Δ)|)(H/x + x^{−1/2})log(|Δ|x) + |R_f(C)| ≤ H|Δ|^{1/2−A+κ+o(1)}, since |R_f| ≤ 2h(Δ)ω(f)log 2x and
  H ≥ |Δ|^{A−κ}.

Each of these is ≤ H|Δ|^{−δ/2−η}, with the margins of [R5] Theorem 3. [R] §2.5 then gives h(Δ)S_C ≥ ν_ΔH/8. Every 𝔞 counted in
S_C lies in C, so it is prime to f. Prime powers 𝔭^k (k ≥ 2) and inert primes (N𝔭 = p²) contribute ≤ C(Hx^{−1/2} + log x)log x
≤ ν_ΔH/(16h(Δ)). Ramified primes have N𝔭 = p | d_K with p ≤ |Δ| < x, so none lies in the interval. This gives (5.1).

*Forms.* Let Q be a primitive positive definite form of discriminant Δ, reduced without loss of generality, and C := [𝔞_Q]^{−1}. For a
degree-one 𝔭 ∈ C with N𝔭 = p ∈ (x, x+H]: p ∤ f because 𝔭 is prime to f, and p ∤ d_K because p > |Δ|. By [C] (F1), (F2) (Cox,
Prop 7.20, 7.22, Thm 7.7, valid for every order), 𝔭 ∩ O_Δ is an invertible O_Δ-ideal of norm p in [𝔞_Q]^{−1}. So p = Q(x′, −y′)
for some x′, y′. □

## 6. Δ > 0: narrow classes, units and regulators of orders (task 5)

**Theorem 3⁺^o.** Theorem 3^o holds for Δ ∈ 𝒟⁺, with the same c, κ, δ, with G_Δ the narrow ring class group and h⁺(Δ) in (5.1).
Every primitive indefinite binary quadratic form of discriminant Δ represents a prime p ∤ Δ in every (x, x+H], x ≥ Δ^A,
H ∈ [xΔ^{−κ}, x], hence a prime ≤ 2Δ^A.

*Proof.* Use [R5] Theorem 3⁺ (the dictionary of [R2] Theorem 4 and [R4] Theorem 3⁺) with these replacements:
- Theorem 2^o and Proposition 1^o for 𝒟⁺(X) in step 2′;
- Lemma G for the genus part. It is proved for both signs at once; d, d′ have the sign of d_K·□ and may be negative, e.g. Δ = 12:
  𝔉 = {1, 12, −3, −4};
- E₂^o and ν_Δ as in §4;
- in step 3, per window h⁺(Δ) ≤ 5Δ^{1/2}log²(3Δ) zeros (every order);
- (5.2) and (a)–(f) as in §5.

*Forms.* This is [R2] Corollary 5 for orders. A primitive indefinite Q is properly equivalent to some (a, b, c) with a > 0, since it
represents a positive integer primitively. 𝔞_Q := Za + Z(−b+√Δ)/2 is an invertible O_Δ-ideal of norm a (Cox Lemma 7.5 and
Thm 7.7(i); the proofs use only primitivity, not the sign of Δ), and N(ua + v(−b+√Δ)/2) = a·Q(u, −v). If 𝔭 ∩ O_Δ lies in the
narrow class of 𝔞_Q^{−1}, then (𝔭 ∩ O_Δ)𝔞_Q = αO_Δ with α ∈ 𝔞_Q and N(α) = ap > 0, so Q(u, −v) = p. □

*What changes for Δ > 0.*
- Narrow classes: P⁺_{K,Z}(f) replaces P_{K,Z}(f) everywhere. Lemma M(a) makes β ≫ 0 by adding a large multiple of a power of f.
  Lemma G(i) is proved for both signs; the inert-prime argument excludes cyclic quartic L_χ, which could a priori occur only for
  real K.
- Units of orders: O_Δ^{×,+} = η_Δ^Z, where η_Δ is the totally positive fundamental unit of the order. The index
  [O_K^{×,+} : O_Δ^{×,+}] can be as large as ≍ f. The units enter only in three places:
  - the unfolding of [R4] Lemma 2⁺ (Poisson with period R⁺_Δ, per order);
  - Lemma P⁺, where γ ↦ γ/e preserves total positivity and u_γ, and each block is the [R4] quantity of the order Δ/e², ≥ 0 with its
    own η_{Δ/e²};
  - the class number bound h⁺(Δ) ≤ 5Δ^{1/2}log²(3Δ), via log η_Δ ≥ 0.962.
- Regulators: only the product h⁺R⁺ = √ΔL_Δ ≤ C√Δ log(3Δ)loglog(9Δ) is used, in Theorem 1^o. Large or small regulators change
  nothing.
- Square Δ: excluded, as in [R4] (Lemma 2⁺ needs a non-square discriminant).
- Family size: for Δ > 0, Σ_{Δ≤X}h⁺(Δ) is smaller than X^{3/2} by a power of log X (arXiv 2610.28001 bounds the sum of class numbers of
  indefinite forms by X^{3/2}√(loglog X)/(log X)^{3/2}). Only the upper bound X^{3/2+o(1)} is used.

## 7. Corollaries

**Corollary 4^o.** For every A > 1 there are c, κ > 0 with the following properties.
(a) For all but O(N^{1−c}) integers n ∈ [1, N], every primitive positive definite form of discriminant −4n, in particular x² + ny²,
represents a prime in every (x, x + x(4n)^{−κ}] with x ≥ (4n)^A. So the least prime x² + ny² is ≤ 2(4n)^A.
(b) The same holds for all but O(N^{1−c}) non-square n ≤ N and every primitive indefinite form of discriminant 4n, in particular
x² − ny² (primes p = x² − ny²).
(c) The same holds for x² + xy + ((n+1)/4)y² (discriminant −n) for all but O(N^{1−c}) integers n ≡ 3 (mod 4).

*Proof.* Δ = ∓4n (resp. −n) is a discriminant for every n, and n ↦ Δ is injective. Apply Theorems 3^o and 3⁺^o with X = 4N. □

[R5]'s corollary covered squarefree n ≡ 1, 2 (mod 4), where −4n is fundamental. New are all other n, e.g. n = 3, 4, 7, 8, 9, 11, 12.
For the least prime of the principal form with squarefree n, Baier–Zhao (Math. Ann. 338, 2007; primes m² + k on average over
squarefree k) remain stronger.

## 8. Where [R]–[R5] used fundamental discriminants, and why nothing fails (task 6)

| Place | Use of "D fundamental" | Replacement |
|---|---|---|
| [R] §0, Proposition 1 | Ĝ(D) consists of primitive characters, distinct π_χ, conductor \|D\| | Lemma M: unique primitive pair, inflation multiplicity ≤ 1 + (X/\|Δ′\|)^{1/2}, ≤ 6 pairs per π. Proposition 1^o: extra X^{1/2} |
| [R] (G), (2.2) | Cox §6: Ĝ₂ = {χ_{d,D/d}}, L = L(χ_d)L(χ_{D/d}) | Lemma G(i), (iii), (iv); imprimitive Euler factors give R_f |
| [R] §2.2 E₂ | d \| D, ≤ X/\|d\| + 1 discriminants per d | Lemma G(ii): dd′ \| Δ, ≤ X/(2\|d\|) + 1 |
| [R] ν_D | χ₁ = χ_d, d ∈ 𝔉(D) | ν_Δ with 𝔉(Δ); \|d\| ≤ \|Δ\|^δ, (I5) |
| [R3] Lemma 3 | L_p(s,χ)^{−1} = 1 − λp^{−s} + (D/p)p^{−2s} for all p | the same with (Δ/p): for p \| f, λ_χ(p) = 0 = (Δ/p) (Lemma D) |
| [R3] Lemma 4 | zeros of L(s,χ), χ primitive | zeros of L(s,χ*) are zeros of L(s,χ) = L(s,χ*)E_χ |
| [R3] §4, [D2] (2.1), [C] Lemma 2 | per-D mean values and large sieve | stated for every order D = f²d_K (and Ĝ⁺ of real orders) |
| [R5] Corollary P | only D ∈ 𝔇(X) on the left | all Δ ∈ 𝒟(X): Q′_Δ ≥ Σ_C\|S_C\|² (Lemma P, block e = 1) |
| [R] §2.4(b), (c), [R5] step 3 | h(D) ≪ \|D\|^{1/2+o(1)} zeros per window | h(Δ) ≤ 5\|Δ\|^{1/2}log²(3\|Δ\|) for every order |
| [R] §3, [R2] Corollary 5 | forms ↔ ideal classes | Cox Prop 7.20, 7.22, Thm 7.7 (every order); narrow version for Δ > 0 |

The extension therefore holds; no step fails. The scout's worry (the family count and Theorem 3 written for fundamental D; the zero
detector counting Ĝ(D) for fundamental D; Siegel zeros through ν_D) is answered by Lemma M with Proposition 1^o, Lemma D with
Theorem 2^o, and Lemma G with E₂^o and ν_Δ respectively.

## 9. Numerical checks ([checks/](checks/); every script under a hard alarm, all processes < 40 MB)

**`genus_orders.py` (Lemma G).** Scope: every non-square discriminant with 3 ≤ |Δ| ≤ 1500, both signs: 1462 discriminants, 550 of them
non-fundamental. Method: for each fundamental d with |d| ≤ 2|Δ| (d = 1 included), test whether χ_d(n) = 1 for the first 600 values
n > 0 of the principal form prime to Δ. These are norms of principal ideals; for Δ > 0, of totally positive generators. Predicted:
- the passing d number exactly 2^μ = 2|G_Δ[2]|, by Gauss's genus count;
- for Δ < 0, |G_Δ[2]| also equals the number of ambiguous reduced forms;
- the passing d come in pairs {d, fund(dd_K)};
- for each pair, dd′ divides Δ (so |dd′| ≤ |Δ|), and the odd part of each d divides Δ.

Result: 0 failures. Mutation (wide instead of narrow classes for Δ > 0, i.e. also using the norms −n < 0): 465 of 712 positive Δ
fail, as they should.

**`famq_all.c` (Theorem 1^o, Δ < 0).** This is [R5]'s `famq.c` with the left side of Theorem 1^o added. Weights:
fam_all := Σ_{Δ}Φ(|Δ|/X)h(Δ)Σ_C|S_C|², over all 17499 (X = 20000) or 34999 (X = 40000) orders with X/4 ≤ |Δ| ≤ 2X, of which
10636 and 21274 are fundamental. The sequence b and its parameters are those of [R5] §8 (σ = 5/6, V = √X/2, pole-killed detector
weight). The last column is the first inequality of the proof (≤ 1 required).

| X | N | b | fam_all/(X^{3/2}‖b‖²) | non-fundamental part | fundamental part | 𝒬′/(X‖b‖²) | fam_all/(max h·𝒬′) |
|---|---|---|---|---|---|---|---|
| 20000 | 1250 | detector | 0.2385 | 0.0803 | 0.1582 | 0.534 | 0.194 |
| 20000 | 5000 | detector | 0.2242 | 0.0746 | 0.1495 | 0.515 | 0.189 |
| 20000 | 20000 | detector | 0.2202 | 0.0716 | 0.1486 | 0.507 | 0.189 |
| 20000 | 30000 | detector | 0.2174 | 0.0696 | 0.1479 | 0.502 | 0.189 |
| 40000 | 5000 | detector | 0.2341 | 0.0804 | 0.1537 | 0.534 | 0.187 |
| 40000 | 20000 | detector | 0.2210 | 0.0735 | 0.1475 | 0.504 | 0.187 |
| 40000 | 40000 | detector | 0.2224 | 0.0740 | 0.1484 | 0.514 | 0.185 |
| 40000 | 60000 | detector | 0.2206 | 0.0724 | 0.1483 | 0.512 | 0.184 |
| 20000 | 5000 | random signs | 0.2602 | 0.0905 | 0.1697 | 0.613 | 0.185 |
| 20000 | 20000 | random signs | 0.2645 | 0.0937 | 0.1708 | 0.632 | 0.182 |
| 20000 | 5000 | pole kept | 12.46 | 4.61 | 7.86 | 37.7 | 0.144 |
| 20000 | 20000 | pole kept | 49.08 | 18.34 | 30.75 | 148.5 | 0.144 |
| 40000 | 40000 | pole kept | 69.65 | 26.13 | 43.52 | 210.1 | 0.142 |

The all-orders family sum stays at the family size, ≈ 0.22X^{3/2}‖b‖², from N = X/16 to N = 1.5X; the non-fundamental orders
contribute about one third. Sequences that keep the trivial character's pole grow like N/√X, as in [R5] §8.

**`least_prime_all.py` (sanity check of the target; not part of any proof).** Scope: every non-square Δ with 5 ≤ |Δ| ≤ 20000:
19857 discriminants, 7699 non-fundamental. Classes: reduced forms for Δ < 0, and ρ-cycles of reduced indefinite forms (proper
equivalence = narrow classes) for Δ > 0. A self-test against 26 known class numbers (h and h⁺, including the orders −108, −36, −400,
−1984, 20, 316) passes. For each class the script finds the least represented prime p ∤ Δ: the least p with Q ~ (p, ±b, ·).
P(Δ) is the maximum over the classes and e(Δ) := log P(Δ)/log|Δ|; statistics over |Δ| ≥ 200:

| | count | median e | max e (at) | max P(Δ)/\|Δ\| |
|---|---|---|---|---|
| Δ < 0, fundamental | 6017 | 0.858 | 1.231 (Δ = −479, h = 25, P = 1997) | 4.17 |
| Δ < 0, non-fundamental | 3884 | 0.855 | 1.054 (Δ = −236, h = 9, P = 317) | 1.34 |
| Δ > 0, fundamental | 6021 | 0.376 | 0.892 (Δ = 1345, h⁺ = 12, P = 619) | 0.49 |
| Δ > 0, non-fundamental | 3753 | 0.444 | 0.891 (Δ = 2889, h⁺ = 18, P = 1213) | 0.50 |

At these sizes the non-fundamental discriminants behave like the fundamental ones; P(Δ) ≤ 4.17|Δ| throughout. For Δ < 0, A ≥ 1 is
forced: the principal form takes no prime value below |Δ|/4.

## 10. Novelty (2026-10-22)

- arXiv listings math.NT 2026-10 (625 entries, the brief's URL) and 2026-09 (812): all titles scanned. Keywords: quadratic form, class
  group/number, ring class, ideal class, class field, least prime, Linnik, Chebotarev, zero density, zero-free, density hypothesis,
  large sieve, Hecke, discriminant, imaginary quadratic, dihedral, Siegel, represent, almost all, order, binary, primes of the form,
  prime. Abstracts read:
  - 2610.28001: average of class numbers of indefinite forms over all non-square discriminants; no primes;
  - 2610.20214: genus character L-functions of orders in quadratic function fields, a product of two quadratic L-functions with local
    factors; the function-field analogue of Lemma G(i), with no prime-counting;
  - 2610.29347, 2610.20316, 2610.27493: unrelated.

  Nothing on least primes, or primes in classes, for non-fundamental discriminants or on average over orders.
- Web search, 3 queries (least prime × forms × almost all discriminants × orders; Baier–Zhao). Only known items. A summary of
  Baier–Zhao (Wolf, "Some conjectures on primes of the form m² + 1") confirms that they average over squarefree k.
- Inherited: [R5]'s adversarial audit ([NOVELTY2](../RCApps5/NOVELTY2.md)) found nothing that states or implies any almost-all
  exponent below 1.1307, even for fundamental D. Theorems 3^o, 3⁺^o contain [R5]'s theorems (fundamental D are 61% of all
  discriminants), so they are new as statements. The literature for all discriminants is per-Δ with large exponents (Weiss;
  Thorner–Zaman 694 and Zaman 455 for ray and ring classes) or GRH (Sardari 1802.06193, every form represents a prime
  ≪ |D|log^{2+ε}|D| under GRH).
- Method: no new analytic input. Lemma P already sums over all orders, so the extension is bookkeeping: primitive pairs and
  multiplicities, genus theory of orders, imprimitive Euler factors. The new content is that this bookkeeping closes.

# CM3: least primes in every ideal class for almost all CM fields over F: every A > 1 (the exponent 1 is sharp)

Date 2026-10-21. Brief: `.scratch/orca/cm3-brief.md`. Checks: [checks/](checks/) (§6). Paper proofs, not refereed, no Lean. Unconditional
(no 2026 half-planes, no Kerr). References: [CM2] = [CM2 NOTE](../CM2/NOTE.md) (notation of its §0 throughout: F, n, U, R(m), S₀, P_F,
admissible, 𝔯_A, ν_𝔫, 𝒪, O_𝔣, G_{K,𝔣}, j_𝔣, η_𝔣, T_𝔫, λ_χ, ω(K,𝔣), Φ, 𝒦(X)), [R5] = [RCApps5 NOTE](../RCApps5/NOTE.md) (the case F = Q),
[CMQ] = [CMQuartic NOTE](../CMQuartic/NOTE.md), [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md).

## Verdicts

| Item | Result | Status |
|---|---|---|
| task 1: relative Lemma P | Lemma P_F (§1): for every O_F-order O_𝔣 = O_F + 𝔣O_K of every CM K/F, the pair count over all elements of the lattices 𝔠𝔯O_𝔣 (no coprimality to 𝔣) equals Σ over squarefree 𝔢 \| 𝔣 of ω(K,𝔣/𝔢)·Σ_{C∈G_{K,𝔣/𝔢}}\|S_C(b_𝔢)\|², hence is ≥ 0 | new, complete proof; exact in 141/141 tests (F = Q(√5), Q(√3), Q(√10)) |
| task 2: tool | Theorem 1 (§2): Σ_{K∈𝒦(X)}Σ_{χ∈Ĉl_K}\|Σ_𝔫 b(𝔫)λ_χ(𝔫)\|² ≪ X^{3/2+C₁δ}N^{1−2σ} for the [CM2] detector sequences, V/4 ≤ N ≤ X^{1−4δ} | new; [CM2] Thm 1 had + N³X^{−1} |
| task 2: density | Theorem 2 (§3): Σ_{K∈𝒦(X)}N^♮_K(σ;t) ≪ X^{2−σ+C₂δ}, ½ + 3δ ≤ σ ≤ 9/10: the density hypothesis of each class group family, on average over K | new; [CM2]: X^{(13−8σ)/(7−2σ)} |
| task 2: least primes | Theorem 3 (§4): [CMQ] Theorem B / [CM2] Theorem 3 for **every A > 1**: all but O(X^{1−c}) CM K/F with prime (or δ₀-rough, odd) 𝔡_{K/F}, N𝔡 ≤ X, have degree-one primes in every ideal class in every (x, x + xD^{−κ}], x ≥ D^A; every class has a prime of norm ≤ 2d_K^A | **new record** (was (3+√69)/10 = 1.1307, [CM2]); 1 is sharp for every K (§4, Remark) |
| task 3: F ≠ Q | division by 𝔢 is by a lattice, not an element (𝔢 may be non-principal and may share primes with 𝔣/𝔢); block weights ω(K,𝔣/𝔢) (roots of unity, relative unit index); classes of F; R(m); narrow classes; family counting: nothing breaks | §5 |
| task 4: numerics | Lemma P_F exact in 141 (K, 𝔣) tests, 372 blocks, ω ∈ {2,4,6,8,10,24}, h_F = 2 included; Poisson main term for norms with square factors within 0.1–1.6% (1–6% at the smallest sizes) | §6 |
| below 1 | blocked: the GRH limit of the explicit formula ([R5] §7); for the principal class 1 is also a lower bound for every K | NOGO |

## 1. Lemma P_F: positivity without coprimality to the conductor

Fix A ∈ Cl_F, 𝔯 := 𝔯_A ([CM2] §0: 𝔯 ∈ S₀ ∪ {O_F}, [𝔯] = A^{−1}) and b supported on admissible 𝔫 ∈ A, with ν_𝔫O_F = 𝔫𝔯. For (K,𝔣) ∈ 𝒪 put

  Q′_{K,𝔣}(b) := Σ_{𝔫,𝔫′} b(𝔫)b̄(𝔫′) Σ_{𝔠|(𝔫,𝔫′)} η_𝔣(𝔠) Σ_{u∈R(ν_𝔫ν_𝔫′)} #{α ∈ 𝔠𝔯O_𝔣 : N_{K/F}α = ν_𝔫ν_𝔫′u}.

This is the count in the proof of [CM2] Lemma 4(ii) with the condition "αO_K prime to 𝔣" dropped. For squarefree 𝔢 | 𝔣 put 𝔣_𝔢 := 𝔣𝔢^{−1},
b_𝔢(𝔪) := b(𝔢𝔪)·1[(𝔪,𝔣) = 1] (supported on admissible 𝔪 in the class A[𝔢]^{−1}) and, for C ∈ G_{K,𝔣_𝔢}, S_C^{(𝔣_𝔢)}(β) := Σ β(N_{K/F}𝔞) over
the integral 𝔞 ∈ C (so 𝔞 ∈ I_K(𝔣_𝔢)).

**Lemma P_F.** For every (K,𝔣) ∈ 𝒪 and b as above,

  Q′_{K,𝔣}(b) = Σ_{𝔢|𝔣 squarefree} ω(K,𝔣_𝔢) Σ_{C∈G_{K,𝔣_𝔢}} |S_C^{(𝔣_𝔢)}(b_𝔢)|²  ≥ 0.

The term 𝔢 = O_F is ω(K,𝔣)Σ_C|S_C(b)|², the summand of [CM2]'s 𝒬_A. In particular Q′_{K,O_F}(b) = ω(K,1)Σ_{C∈Cl_K}|S_C(b)|².

*Proof.* (a) *Only equal 𝔣-parts survive.* Let 𝔭 | 𝔣 with v_𝔭(𝔫𝔫′) = 1. Then 𝔭 ∉ P_F (it divides an admissible ideal), so 𝔭 ∤ 𝔯. Terms with
𝔭 | 𝔠 vanish (η_𝔣(𝔭) = 0). In the others every α counted lies in O_𝔣 (𝔠𝔯 is integral) and has v_𝔭(N_{K/F}α) = v_𝔭(𝔫𝔫′𝔯²) = 1, which no element
of O_𝔣 has ([CM2] Lemma 2(c)). As 𝔫, 𝔫′ are squarefree, only pairs with (𝔫,𝔣) = (𝔫′,𝔣) =: 𝔢 contribute; 𝔢 is squarefree and admissible. Write
𝔫 = 𝔢𝔪, 𝔫′ = 𝔢𝔪′, so (𝔪𝔪′,𝔣) = 1; the 𝔠 with η_𝔣(𝔠) ≠ 0 are prime to 𝔣, hence 𝔠 | (𝔪,𝔪′) and η_𝔣(𝔠) = η_{𝔣_𝔢}(𝔠) = η(𝔠).

(b) *Division by 𝔢 is a lattice identity.* For such a pair and 𝔠: {α ∈ 𝔠𝔯O_𝔣 : N(α)O_F = 𝔫𝔫′𝔯²} = {α ∈ 𝔢𝔠𝔯O_{𝔣_𝔢} : N(α)O_F = 𝔫𝔫′𝔯²}.
Compare completions (an O_F-lattice is the intersection of its completions). At 𝔭 ∤ 𝔢 they agree. At 𝔭 | 𝔢 (v_𝔭(𝔠𝔯) = 0, f := v_𝔭(𝔣) ≥ 1) they
are O_{𝔭,f} and πO_{𝔭,f−1} ⊆ O_{𝔭,f}, and by [CM2] Lemma 2(c) the α ∈ O_{𝔭,f} with 𝔭² | Nα (here v_𝔭(Nα) = 2) are exactly those in πO_{𝔭,f−1}.
No generator of 𝔢 is needed, and 𝔢 may share primes with 𝔣_𝔢 (when 𝔭² | 𝔣).

(c) *Elements, units, classes.* Fix the pair and 𝔠; put O := O_{𝔣_𝔢}, 𝔩 := 𝔢𝔠𝔯, 𝔐 := 𝔪𝔪′𝔠^{−2} (integral, prime to 𝔣), m := ν_𝔫ν_𝔫′ and
𝒜 := {α ∈ 𝔩O : N(α)O_F = mO_F} (an O^×-set). Choose an integral 𝔩′ prime to 𝔣 in the class of 𝔩 in Cl_F and λ ∈ F^× with 𝔩 = λ𝔩′.
Since λ ∈ F, α ↦ α/λ is an O^×-equivariant bijection 𝔩O → 𝔩′O, and N(α/λ)O_F = 𝔩′²𝔐. For α′ ∈ 𝔩′O with N(α′)O_F = 𝔩′²𝔐: the norm is prime
to 𝔣_𝔢, so α′O_K is prime to 𝔣_𝔢 ([CM2] Lemma 2(e)) and lies in P_{K,O_F}(𝔣_𝔢) (α′ ∈ O); 𝔟 := α′(𝔩′O_K)^{−1} is integral (α′ ∈ 𝔩′O_K), prime
to 𝔣_𝔢, of relative norm 𝔐 and of class j_{𝔣_𝔢}(𝔩′)^{−1}. Conversely, for an integral 𝔟 ∈ I_K(𝔣_𝔢) with N𝔟 = 𝔐 in that class, 𝔟𝔩′O_K is an
integral element of P_{K,O_F}(𝔣_𝔢), so 𝔟𝔩′O_K = α′O_K with α′ ∈ O unique up to O^× ([CM2] Lemma 2(a)), and α′ ∈ 𝔩′O (at 𝔭 | 𝔣_𝔢,
(𝔩′O)_𝔭 = O_𝔭; elsewhere (𝔩′O)_𝔭 = 𝔩′O_{K,𝔭} ∋ α′). Hence the O^×-orbits of 𝒜 (O^× acts freely) correspond bijectively to these 𝔟.
Units exactly as in [CM2] Lemma 4(ii): N(α) = mu₀ with mu₀ ≫ 0 (relative norms from a CM field are totally positive), u₀ = uε² with u ∈ R(m)
unique and ε ∈ U unique up to sign, so each U-orbit of 𝒜 meets 𝒜′ := {α ∈ 𝒜 : Nα ∈ mR(m)} in one pair ±α, and [O^× : U] = ω(K,𝔣_𝔢)/2
(N induces O^×/μ(O)U ≅ N(O^×)/U², norm-one units being roots of unity). Classes: j_{𝔣_𝔢}(𝔩′) depends only on [𝔩] = [𝔢][𝔠]A^{−1} ([CM2] Lemma 2(b)),
so j_{𝔣_𝔢}(𝔩′)^{−1} = j′·j_{𝔣_𝔢}(𝔠)^{−1} with j′ := j_{𝔣_𝔢}(A[𝔢]^{−1}), the common class j_{𝔣_𝔢}(𝔪) of the 𝔪 in supp b_𝔢. Therefore

  Σ_{u∈R(m)} #{α ∈ 𝔢𝔠𝔯O_{𝔣_𝔢} : Nα = mu} = |𝒜′| = ω(K,𝔣_𝔢)·[j′j_{𝔣_𝔢}(𝔠)^{−1}](T_𝔐)   in Z[G_{K,𝔣_𝔢}].

(d) *The block.* By (a)–(c) the pairs with (𝔫,𝔣) = (𝔫′,𝔣) = 𝔢 contribute ω(K,𝔣_𝔢)Σ_{𝔪,𝔪′}b_𝔢(𝔪)b̄_𝔢(𝔪′)Σ_{𝔠|(𝔪,𝔪′)}η_{𝔣_𝔢}(𝔠)[j′j(𝔠)^{−1}](T_{𝔪𝔪′𝔠^{−2}})
= ω(K,𝔣_𝔢)[j′](S S̃), S := Σ_𝔪 b_𝔢(𝔪)T_𝔪, S̃ := Σ_𝔪 b̄_𝔢(𝔪)T_𝔪, by the Hecke relations [CM2] Lemma 3 in Z[G_{K,𝔣_𝔢}] (all 𝔪 prime to 𝔣_𝔢,
squarefree). By [CM2] Lemma 2(b) ([𝔞]^{−1} = [𝔞ᶜ]j(N𝔞)^{−1}; 𝔞 ↦ 𝔞ᶜ preserves norms and coprimality) S^∨ := Σ_C S̄_C C^{−1} = j′^{−1}S̃, so
[j′](SS̃) = [1](SS^∨) = Σ_C|S_C^{(𝔣_𝔢)}(b_𝔢)|². □

*Remarks.* (1) For F = Q (ω = w, U = {±1}, h_F = 1) Lemma P_F is [R5] Lemma P (there 𝔢 = (e) and division by the element e). (2) Orders with
extra units need no special treatment: their blocks carry ω(K,𝔣_𝔢) > 2 (e.g. Q(ζ₅)/Q(√5): ω = 10 for O_K, 2 for every proper order).
(3) Q′_{K,𝔣} ≥ ω(K,𝔣)Σ_C|S_C|² term by term, so the family sum below dominates [CM2]'s 𝒬_A.

**Corollary P_F (dual identity without Möbius terms).** With 𝒩_𝔠 as in [CM2] Lemma 4(ii),

  𝒬′_A(b) := Σ_{(K,𝔣)∈𝒪} Φ(N(𝔣²𝔡_{K/F})/X) Q′_{K,𝔣}(b) = Σ_{𝔫,𝔫′} b(𝔫)b̄(𝔫′) Σ_{𝔠|(𝔫,𝔫′)} Σ_{u∈R(ν_𝔫ν_𝔫′)} 𝒩_𝔠(ν_𝔫ν_𝔫′u; 𝔠𝔯_A; X),

the 𝔴 = O_F part of [CM2] Lemma 4(ii); and for b supported on admissible ideals, with b_A := b·1_{[𝔫]=A},

  Σ_{K∈𝒦(X)} Σ_{χ∈Ĉl_K} |Σ_𝔫 b(𝔫)λ_χ(𝔫)|² ≤ C_F X^{1/2}log^{2n−1}X · Σ_{A∈Cl_F} 𝒬′_A(b_A).

*Proof.* The identity is the definition of 𝒩_𝔠 (finite sums). The inequality is [CM2] §4.1 (Lemma 4(i), Cauchy–Schwarz over A, h_K ≤
CX^{1/2}log^{2n−1}X) with Lemma P_F: for K ∈ 𝒦(X) the term (K, O_F) has Φ = 1 and equals ω(K,1)Σ_C|S_C(b_A)|² ≥ 2Σ_C|S_C(b_A)|², and every
other term Φ·Q′_{K,𝔣}(b_A) is ≥ 0. □

## 2. Theorem 1: the family large sieve at the family size for all lengths below the conductor

**Theorem 1.** Let 0 < δ ≤ 1/100, X ≥ X₀(F,δ), ½ + 2δ ≤ σ ≤ 9/10, |t′| ≤ X^δ, 2 ≤ V ≤ X^{1/2}, N := 2^j with V/4 ≤ N ≤ X^{1−4δ}, and b = b_{j,t′}
or ∂_{t′}b_{j,t′} ([CM2] §3). Then, with C₁ = C₁(n) (C₁ = 2n + 26 suffices),

  Σ_{K∈𝒦(X)} Σ_{χ∈Ĉl_K} |Σ_𝔫 b(𝔫)λ_χ(𝔫)|² ≤ C(F,δ) X^{3/2+C₁δ} N^{1−2σ}.

*Proof.* By Corollary P_F, 𝒬′_A(b_A) = Σ_𝔠 Q_{𝔠,1}, with Q_{𝔠,1} the [CM2] §4.1 sum at 𝔴 = O_F (X′ = X, 𝔩 = 𝔠𝔯_A, L = 4N/N𝔠, M₀ = N(𝔫₂𝔫₂′)),
and it suffices to show Σ_𝔠|Q_{𝔠,1}| ≤ C X^{1+C₁δ/2}‖b‖². [CM2] §4.4–§4.8 bound Q_{𝔠,𝔴} using only the following, which hold at 𝔴 = O_F for every
N ≤ X^{1−4δ}:
- (i) Lemma 9 (Poisson over the ideal lattice): X′ ≥ C₃N𝔠√M₀X^δ and M₀ ≤ X². Here N𝔠√M₀ = √(N𝔫N𝔫′) ≤ 4N ≤ 4X^{1−4δ}, so
  C₃N𝔠√M₀X^δ ≤ 4C₃X^{1−3δ} ≤ X; and M₀ ≤ 16N² ≤ X². This was the only place where [CM2] needed N𝔴 ≤ G₂ = X^{1−2δ}/(4^{n+2}N).
- (ii) Lemma 10 (local densities): 𝔫₂, 𝔫₂′ squarefree, prime to 2S₀, so v_𝔭(𝔫₂𝔫₂′) ∈ {0,1,2}; the case 2 (pairs with common factors, kept by
  Lemma P_F) is in Lemma 10(a) and in (S).
- (iii) §4.5: E₀ = 2^{n+3}N/(N𝔠√X) ≥ 2v^{−1/2} for all M₀ ≤ L² (v = X/(4ⁿM₀)).
- (iv) §4.7 (classical large sieve, Lemma 12): E₀² ≤ L, i.e. 4^{n+2}N ≤ N𝔠X: true.
- (v) §4.8 (main terms): the pole is killed by (Z); Lemma 14 (Jutila over F) with Q = E₀ ≤ X.
Then, exactly as in [CM2] §4.9: the elements of F (Lemma 8(a), 𝔫₂ = 𝔫₂′) give ≪ X^{1+2ε}‖b‖²; the secondary terms ≪ X^{8δ}X‖b‖²; the main terms
≪ X^{(12+n)δ}X‖b‖² (V ≤ X^{1/2}); Poisson errors O(X^{−5}). Corollary P_F multiplies by C_F X^{1/2}log^{2n−1}X and sums over the h_F classes A. □

[CM2] §4.3 (large 𝔴: the source of N³X^{−1}) is no longer used. The two minor text repairs noted by the [CM2] reviewer (imprimitive Euler factors
cost Q^{O(ε)} in Lemmas 12, 14; conductor bookkeeping in Lemma 14) concern §4.7–§4.8 and are unaffected.

## 3. Theorem 2: the density hypothesis on average over CM fields

For K CM over F: N^♮_K(σ;t) := Σ_{χ∈Ĉl_K∖{1}} #{ρ : L(ρ,χ) = 0, σ ≤ β ≤ 9/10, |γ − t| ≤ 1}.

**Theorem 2.** For 0 < δ ≤ 1/100, ½ + 3δ ≤ σ ≤ 9/10, T := |t| + 3 ≤ X^δ/4 and X ≥ X₀(F,δ):  Σ_{K∈𝒦(X)} N^♮_K(σ;t) ≤ C(F,δ)X^{2−σ+C₂δ}, C₂ = C₂(n).

*Proof.* [CM2] §5 (= [R3] §4 over F) with Y := X^{1−5δ} instead of X^{5/(7−2σ)−5δ}, V := X^{1/2}/2 (≤ d_K^{1/2}) and Z := 2Y log²Y, so the
class I lengths satisfy N ≤ 2Z ≤ X^{1−4δ}. Class II (per-field mean value [CMQ] (2.1′) on Re s = ½ + δ): |𝒳_II(K)| ≤ C Y^{1/2+δ−σ}X^{1/2+2δ}Tⁿ,
summed over |𝒦(X)| ≪ X^{1+δ}: ≤ C X^{3/2+3δ}Y^{1/2+δ−σ}Tⁿ ≤ C X^{2−σ+(6+n)δ} (as (1−5δ)(½+δ−σ) ≤ ½ − σ + 3δ for σ ≤ 9/10). Class I: Theorem 1 for
|τ| ≤ X^δ/2, the per-field sieve [CMQ] Lemma 3 with the decay of Φ̂_β for larger |τ| (as in [CM2] §5): Σ_K|𝒳_I(K)| ≤ X^{O(δ)}Σ_{V/4≤2^j≤2Z}X^{3/2}2^{j(1−2σ)}
≤ X^{3/2+O(δ)}V^{1−2σ} ≤ X^{2−σ+O(δ)}. Zeros per window: c₀log(d_KT^{2n}) ([CMQ] (F4)). □

For each K the family Ĉl_K has h_K ≍ X^{1/2±ε} characters of conductor d_K ≍ X; its density hypothesis is Σ_χN(σ;t,χ) ≪ X^{1−σ+ε}, i.e. X^{2−σ+ε}
summed over 𝒦(X). With [CMQ] Proposition C (better for σ > (√41−3)/4 = 0.8508): Σ_KΣ_{χ∈Ĝ_c(K)}N(σ,T;χ) ≪ (XT)^{O(ε)}X^{min{2−σ, 10(1−σ)/(3−2σ)}}.

## 4. Theorem 3: least primes in every ideal class, every A > 1

**Theorem 3.** Let F be totally real of degree n, δ₀ ∈ (0,1) and A > 1. There are effective c, κ > 0 and X₀ (depending on F, δ₀, A) such that for
X ≥ X₀ all K ∈ ∪_{Y≤X}𝒦_F^{δ₀}(Y) (CM, 2 ∤ 𝔡_{K/F}, N𝔭 > D^{δ₀} for all 𝔭 | 𝔡_{K/F}, [CMQ] §4; δ₀ ≥ ½: prime 𝔡) outside a set of size ≤ X^{1−c}
satisfy, for every C ∈ Cl_K, every x ≥ D^A (D := N𝔡_{K/F}) and every H ∈ [xD^{−κ}, x],

  #{𝔓 ∈ C prime : N𝔓 = p prime, x < p ≤ x + H} ≥ H/(16h_K log 2x).

So every ideal class of K contains a prime of degree one with norm ≤ 2D^A ≤ 2d_K^A, some β ∈ O_K has prime norm ≤ 2d_K^A, and for every
integral 𝔞 the relative norm form of 𝔞 takes a value generating 𝔭N_{K/F}𝔞 with 𝔭 of degree one, N𝔭 ≤ 2d_K^A ([CMQ] Corollary B1).
The family is nontrivial: prime 𝔡 alone give ≫_F X/log X fields ([CMQ] §4).

*Remark (sharpness).* For every CM K/F, every β ∈ O_K∖F has N_{K/Q}β ≥ d_K/(4ⁿd_F²) ([CMQ] Remark 2.3, from the relative determinant lemma with
β′ = 1), and a degree-one prime of the principal class is βO_K with β ∉ F. So the least degree-one prime of the principal class is ≥ d_K/(4ⁿd_F²) for every
K, and by Theorem 3 it is ≤ 2d_K^{1+ε} for almost all K: the exponent 1 is optimal (for F = Q: x² + ny² represents no prime below n).

*Proof.* [CM2] Theorem 3 (= [CMQ] Theorem B with the changes of [R3]) with steps 1, 2 replaced. The statement for A implies it for every A′ > A
(same exceptional set), so let 1 < A ≤ 1.2.
1′. *Parameters.* ε′ := (A−1)/80; δ := ε′/(10C₃) with C₃ := C₂ + 2 (so δ ≤ (A−1)/800 and δ ≤ A/4400); κ := min{δ/4, κ₀, δ₀/(n+2)} (κ₀ = κ₀(F,A)
   of [CMQ] step (e)); ε₀ := ε₁ := κ/10; σ₀ := ½ + 4δ; σ_j := σ₀ + jε′/A (0 ≤ j ≤ J, J minimal with σ_J ≥ 1 − ε′/A); σ_g := 0.995;
   T_X := 2X^{κ+ε₁} (≤ X^δ/4).
2′. *E₁ (Markov).* For σ ∈ [σ₀, 1): Σ_{K∈𝒦^{δ₀}(X)}N^c_K(σ,T_X) ≤ C X^{F*(σ)+C₃δ} with F*(σ) := 2 − σ on [σ₀, 9/10] (Theorem 2 over the
   ≤ 2T_X unit windows for β ≤ 9/10, 𝒦^{δ₀}(X) ⊆ 𝒦(X); the zeros with β > 9/10 by [CMQ] Proposition C, X^{F₂(9/10)+O(κ)} = X^{5/6+O(κ)}) and
   F*(σ) := F₂(σ) = 10(1−σ)/(3−2σ) on (9/10, 1) ([CMQ] Proposition C; T-power X^{O(κ)}). Then φ(σ) := F*(σ) − A(1−σ) = 1 − (A−1)(1−σ) ≤ 1 − 8ε′
   on [σ₀, 9/10] and φ ≤ F₂ ≤ 5/6 ≤ 1 − 8ε′ on (9/10, 1). With E₁(X) := {K : some j has N^c_K(σ_j,T_X) > X^{A(1−σ_j)−2ε′}}:
   #E₁(X) ≤ Σ_j X^{2ε′−A(1−σ_j)}·C X^{F*(σ_j)+C₃δ} ≤ C J X^{1−6ε′+ε′/10} ≤ C J X^{1−5ε′}.
3. *E₂* (genus characters ψ ∈ 𝒢(K)∖Ĉl_F⁺ with a zero β > σ_g) as in [CMQ]: #E₂(X) ≤ X^{1−0.9δ₀}.
4. *Decomposition* ([CMQ] (5.2)), for K ∉ E₁ ∪ E₂, D ≥ D₀, T₁ := 2D^{κ+ε₁}:
   (a) non-genus zeros with β ≥ σ₀: E_a ≤ HJ2^AX^{−ε′} (zeros in [σ_j, σ_{j+1}) weighted by x^{σ_{j+1}−1} ≤ 2^AX^{−A(1−σ_j)+ε′});
   (b) zeros of all χ ∈ Ĉl_K with β < σ₀, |γ| ≤ T₁: ≤ 2T₁c₀h_K log(d_KT₁^{2n}) zeros, each with x^{β−1} ≤ D^{−A(1/2−4δ)}, so
       E_b ≤ HD^{−(A−1)/2+4Aδ+κ+ε₁+o(1)} ≤ HD^{−(A−1)/2+6δ}, and (A−1)/2 − 6δ > (A−1)/3;
   (c) |γ| > T₁: E_c ≤ HD^{−1}; (d) rough genus zeros with σ₀ ≤ β ≤ σ_g: E_d ≤ HD^{−0.005A+κ+ε₁+o(1)}; (e) zeros of ζ_{H_F⁺}: ≤ H/16 (κ ≤ κ₀);
   (f) Σ|R| ≤ C(h_K + |𝒢(K)|)(H/x + x^{−1/2})log(d_Kx) ≤ H(D^{1/2−A} + D^{1/2+κ−3A/2})D^{o(1)}.
   Prime powers and primes of degree ≥ 2 contribute ≤ C(Hx^{−1/2} + log x)log^{2n+1}x ≤ H/(16h_K), as x ≥ D^A ≥ D^{1+2δ} and h_K ≤ D^{1/2+o(1)}.
   So h_KS_C ≥ H/2 − H/16 − 4HD^{−η} ≥ H/4 (η = η(F,δ₀,A) > 0), and the count follows as in [CMQ]. Exceptional set: c := min{ε′, 0.9δ₀}/2. □

*Where A > 1 is used.* Only in (b) and in Markov near σ₀ (φ(σ₀) = 1 − (A−1)(½ − 4δ)): h_K ≍ D^{1/2} zeros per window near Re s = ½, each
weighted by x^{−1/2}. This is exactly what GRH gives by the same method ([R5] §7), and by the Remark it is also forced for the principal class.

## 5. Task 3: every place where F ≠ Q changes something

| Place | F = Q ([R5]) | F ≠ Q (here) | Effect |
|---|---|---|---|
| removing 𝔢 | γ ↦ γ/e, e ∈ Z; Δ/e² ≡ Δ mod 8 | 𝔢 is an ideal, possibly non-principal (h_F > 1) and sharing primes with 𝔣/𝔢 (𝔭² \| 𝔣): (b) is a lattice identity via [CM2] Lemma 2(c); (c) moves 𝔢𝔠𝔯 to an 𝔩′ prime to 𝔣/𝔢 by λ ∈ F^×; orders O_{𝔣/𝔢} need no discriminant congruence, 2 \| 𝔣 allowed | none |
| units | w_e/2, w_e = \|O_{Δ_e}^×\| | per block ω(K,𝔣/𝔢) = \|μ(O_{𝔣/𝔢})\|·q(O_{𝔣/𝔢}), q = [N(O^×):U²] ∈ {1,2} (relative unit index; Hasse index for O_K); both may differ between blocks of one order | weights ≥ 2; positivity intact |
| totally positive generators | sign of M | N_{K/F}α ≫ 0 always; R(m), m = ν_𝔫ν_𝔫′, is the same set in every block; λ ∈ F^× needs no sign (N(α/λ) = N(α)/λ²); if 𝔫𝔫′𝔯² has no totally positive generator both sides are 0 | none |
| narrow classes | D > 0: Cl⁺ | K totally complex: G_{K,𝔣} has no sign conditions, j_𝔣 factors through the wide Cl_F ([CM2] Lemma 2(b) for every λ ∈ F^×); narrow classes of F enter only the analytic part (sgn^a in Lemma 11; genus characters in Ĉl_F⁺, ζ_{H_F⁺}) | none |
| classes of F | none | b split by A; b_𝔢 lives in A[𝔢]^{−1}; the class j(𝔩′)^{−1} = j(A[𝔢]^{−1})j(𝔠)^{−1} matches the Hecke relation | none |
| S₀, dyadic primes | none | 𝔣 prime to S₀ ⇒ 𝔣/𝔢 prime to S₀; 𝔢 admissible (odd, prime to S₀) | none |
| dual count | t² − Δy² = 4M | t ∈ 𝔩 ⊂ Rⁿ with t² − 4μ ≪ 0 gives K_t CM and the orders 𝔣 \| 𝔣(t) ([CM2] Lemma 4(iii)); weights Φ(N(t²−4μ)/(N(𝔩𝔶)²X)) | unchanged (𝔴 = 1 only) |
| Poisson hypothesis | X ≥ 64c√M X^δ | X ≥ C₃N𝔠√M₀X^δ ⟸ N ≤ X^{1−4δ} | none |
| family and counting | fundamental D, h ≤ CX^{1/2}log X | 𝒦(X) by N𝔡_{K/F} ∈ (X/2, X], \|𝒦(X)\| ≪ X^{1+ε}, d_K = d_F²N𝔡, h_K ≤ CX^{1/2}log^{2n−1}X; each K ∈ 𝒦(X) is the order (K, O_F) with Φ = 1 | none |
| roughness | all D (Siegel factor ν_D) | Theorem 3 only for δ₀-rough odd 𝔡 (genus characters of small conductor: [CMQ] §7.1); Theorems 1, 2 for all CM K | inherited from [CMQ], not caused by Lemma P_F |

The transfer does not fail anywhere.

## 6. Numerical checks ([checks/](checks/), pari/gp 2.19, hard alarm on every run)

- `lemmaP_rel.gp` → `lemmaP_rel.out` (19 s CPU, 19 MB). For K = F(√−β) and conductors 𝔣: Q′_{K,𝔣}(b) by direct enumeration of α (t in a box,
  √(t²−4μ) ∈ K, membership in the Z-lattice 𝔠𝔯 + 𝔠𝔯𝔣O_K = 𝔠𝔯O_𝔣; no class groups, no Lemma 1), against Σ_𝔢ω(K,𝔣/𝔢)Σ_C|S_C(b_𝔢)|² with
  G_{K,𝔣′} = Cl_{𝔣′O_K}(K)/⟨images of (O_F/𝔣′)^×⟩ (bnrinit) and ω from the roots of unity and the fundamental unit lying in O_{𝔣′}; b random
  Gaussian integers, so equality is exact.
  - F = Q(√5): 21 squarefree odd 𝔫 with N𝔫 ≤ 200 from the primes above 3, 5, 11, 19, 29, 31 (the algebra needs only squarefree, odd, prime
    to S₀; small primes make 𝔢 ≠ 1 frequent), 6 fields (Q(√5,i), Q(√5,√−3), Q(ζ₅), β = 7, 5+2w, 11+3w; w = (1+√5)/2),
    9 conductors (𝔭₁₁, 𝔭₁₁², 𝔭₁₁𝔭₁₉, (3), (3)², 𝔭₁₁·(2), (√5), (11), 𝔭₁₁²(√5)).
  - F = Q(√3) (U⁺ ≠ U², so |R(m)| = 2): 5 fields incl. Q(ζ₁₂) (ω(O_K) = 24, q = 2) and β = 2 (Hasse index 2); 7 conductors incl. 𝔭₂𝔭₁₁.
  - F = Q(√10) (h_F = 2): class A ≠ 1 with reference prime 𝔯 | 3, and A = 1; non-principal 𝔢 above 3, 5, 13; 𝔣 = 𝔭₃′², 𝔭₂𝔭₃′, …; 4 fields.
  - Result: 141/141 (K, 𝔣) tests exact; 372 blocks (360 nonzero, 226 nonzero with 𝔢 ≠ 1); block weights ω ∈ {2, 4, 6, 8, 10, 24} (q = 2 in 16
    blocks); all 12073 pairs with different 𝔣-parts give 0 (step (a)); 22085 per-pair, per-𝔠 counts of step (b) agree; the ring class number
    formula |G_{K,𝔣}|·|(O_F/𝔣)^×|·[O_K^×:O_𝔣^×] = h_K|(O_K/𝔣)^×| holds for all 372 groups (|G| up to 1760, h_K ∈ {1, 2, 4, 24}).
- `poisson_e2.gp` → `poisson_e2.out` (3 min CPU, 13 MB): [CM2] Lemma 9 with Lemma 10 for μ with square factors (e = 2, the pairs Lemma P_F
  keeps), F = Q(√5), 𝔩 = O_F, 𝔠 = 1: brute force 2Σ_tΣ_{𝔶|𝔣(t)}Φ(N(t²−4μ)/(N𝔶²X′)) (𝔣(t) from rnfdisc) against
  (2^{1−n}/√d_F)(X′/√M₀)Σ_𝔶A_μ(𝔶)K₂(N𝔶²v) (A_μ by brute force from the criterion C_k, the dyadic prime included; K₂ by quadrature, Mellin check
  12.06947 vs 12.06955). Ratios brute/main at X′/√M₀ = 20, 40: μ = π₁₁²π₁₉π₂₉: 1.016, 0.989 (80, 160: 0.9965, 1.0008); μ = π₁₁²π₁₉²π₂₉: 0.996,
  1.0015; squarefree control π₁₁π₁₉π₂₉π₃₁: 0.992, 1.0035; μ = π₁₁²·19: 0.942, 0.971, then 1.014, 1.012 at 80, 160 (smallest M₀ = 43681).

## 7. Novelty (2026-10-21)

- arXiv monthly listings math.NT 2026-09 (812 entries) and 2026-10 (625, through 2610.29876): titles scanned for CM field, totally real, Hecke,
  class group/number, ideal/ring class, least prime, Linnik, Chebotarev, zero density/zero-free, density hypothesis, large sieve, almost all,
  families of number fields, Siegel, dihedral, Hermitian, norm form, low-lying (118 hits). Abstracts read for the closest: L. Yang 2610.22728
  (non-split Motohashi formula for toric periods over number fields: first and second moments of class group L-functions of E/F; no zero
  density, no primes), Gao 2610.27493 (moments of angular Hecke L-functions, class number one), Diao 2610.28001 (average of class numbers of
  indefinite forms), Garcia–Grenié–Lee–Molteni 2610.14257 (residues under GRH); known from [R5]: Frei–Sofos–Widmer 2610.18491, Liu 2610.12234
  and 2610.26835, Cossaboom 2609.27715, de Faveri 2610.04045, Kerr 2610.21099. None treats least primes in ideal classes or class-group zero
  density on average over CM fields.
- arXiv search (8 queries): only known items (Yang 1801.05562, Thorner–Zaman 1604.01750, Zaman 1508.00287, 2105.14181, 1902.08640,
  Pierce–Turnage-Butterbaugh–Wood). GXL evidence API: HTTP 429 (daily quota used up).
- Prior bounds: [CM2] A > 1.1307 and [CMQ] A > 3/2 (this project); per field, least prime ideal bounds with large exponents (Thorner–Zaman,
  Zaman), GRH d_K^{1+ε}. Nothing found states or implies an almost-all exponent below 1.1307 for CM fields over any F ≠ Q, nor Theorems 1–2 for F ≠ Q.
- Method: [R5]'s Lemma P mechanism, transferred; new here: the lattice form of the division step (no generator of 𝔢), the per-block unit
  weights ω(K,𝔣/𝔢), and the class bookkeeping j(𝔩′)^{−1} = j(A[𝔢]^{−1})j(𝔠)^{−1}.

## 8. Barrier and open directions

- *A ≤ 1*: NOGO for this method (as [R5] §7: h_K x^{−1/2} ≤ D^{−η} needs x ≥ D^{1+2η}); for the principal class A ≥ 1 is necessary for every K.
- *All CM fields over F* (small genus conductors): still blocked as [CMQ] §7.1; Lemma P_F does not touch it.
- *F varying* (all quartic CM fields): constants depend on F through d_F, h_F, S₀, the units and (HT), (TZ) ([CMQ] §7).
- *Sectors / infinite-type Hecke characters of K; non-CM K/F*: not attempted (angular weights; relative units of positive rank).

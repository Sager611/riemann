# T4 plan

Page: discoveries/prefix-reversals.md (Claim); proofs: ../PrefixRev/NOTE.md; reused Lean: ../PrefixRev/Lean (copied, namespace Riemann.Deep.LeanPrefix).
Build: bash Riemann/Deep/LeanPrefix/tools/mk.sh <Mod>... (one lean at a time, leanrun.sh, .olean into .lake/build/lib/lean).

## Main statement (Main.lean, prefix_reversals)
Hypotheses: (hOff : LZZW2026_Thm3_1) (hJ : Jones2014_Thm1_2_k0). Conclusion, with g = gcd(m+1, k+1), α = k+1, q = n/α:
1. Theorem A: 2 ≤ k < m < n, m + k < n ⇒ (G = ⊤ ↔ n odd ∧ m + k = n - 1 ∧ gcd = 1 ∧ ¬(k+1) ∣ m ∧ ¬(n % 4 = 1 ∧ 4 ∣ k)).
2. B.1 (line, any parity): y ∈ orbit x ↔ y+1 ≡ ±(x+1) (mod g); card (orbitRel.Quotient) = g/2 + 1; IsPretransitive ↔ g = 1.
3. B.2 (n even, g = 1): IsPretransitive ∧ ¬IsPreprimitive ∧ a block B with 2·|B| = n.
4. B.3 (n odd, α ∣ n): ¬IsPreprimitive, interval blocks, G ≃* (ZMod q → DihedralGroup α) ⋊ DihedralGroup q (q-gon action), |G| = (2α)^q·2q.
5. B.4 (n odd, g = 1, α ∤ n): IsPreprimitive ∧ alternatingGroup ≤ G ∧ (n % 4 = 1 ∧ 4 ∣ k → G = alternatingGroup) ∧ (¬… → G = ⊤).
Reading of the Claim: "when α ∣ n … otherwise primitive ⊇ A_n" carries NOTE Thm B's hypotheses (n odd; g = 1): false without them,
e.g. (12,8,3) (α ∣ n, |G| < 8^3·6), (10,6,3) (n even, g = 1, α ∤ n: two blocks), (11,8,2) (g = 3: intransitive).

## Lemma DAG
| Module | Content | Status |
|---|---|---|
| Defs … PrimB2, Small | copied: Lemmas 1-3, (★_d), parity, necessity, interval blocks, kernel G = ⊤ for (7,4,2), (11,6,4), (13,10,2) | compiled |
| Prim0 | IsBlockOf, IsPrimitive, fix_rk_card, rk_ne_one, even_of_coprime | compiled |
| Bridge | IsBlockOf ↔ MulAction.IsBlock; IsPrimitive ↔ MulAction.IsPreprimitive | todo |
| OrbitsB | P^t reaches every y with y+1 ≡ x+1 (mod g) (first return); orbit = O_c; orbit count; transitivity | todo |
| Prop5A | pure arithmetic: Star' (runs of y_s = sΔ mod N, Δ = αd mod N, constant on I = [1, α-1] for runs i ≤ d-2) ⇒ α ∣ n | todo |
| Prop5B | StarD ⇒ Star'; Prop5Core; primitivity (α ∤ n) | todo |
| Lemma7 | involution form: g² = 1, g ≠ 1 induced by a σ-semilinear bijection of F^d (d ≥ 2) on ℙ ⇒ 2·nfix ≤ |ℙ| unless (d,q) = (2,4) | todo |
| Jones | literature input Jones2014_Thm1_2_k0 (verbatim TeX in docstring); A_n ≤ G on the odd primitive line | todo |
| Small23 | kernel G = ⊤ for (23,18,4), (23,16,6), (23,12,10) (Jones case (c)) | todo |
| Wreath | dihedralPerm (q-gon action), Wreath α q, faithful imprimitive action, image = G (Lemma 6), iso, order | todo |
| Main, Check | assembly; #print axioms | todo |

## Prop 5 core (Prop5A), proof used
N = α + β, n = N - 1 = dL (d, L ≥ 3 odd), Δ = αd mod N, Y(s) = sΔ mod N (s ∈ [1, n], a bijection onto [1, n]), LΔ ≡ -α.
Star': for i ≤ d-2 and s, s' ∈ (iL, (i+1)L]: Y(s) ∈ I ↔ Y(s') ∈ I (I = [1, α-1]).
- Case 1 (no run i ≤ d-2 meets I): τ with Y(N-τ) = 1, τ ≤ L-1; induction xτ ≤ L-1 (x < α); ατ = L (Y injective) ⇒ α ∣ L.
- Case 2 (run i ⊆ I): consecutive differences all equal D ≡ Δ, |D| ≤ (α-2)/(L-1); LD + α ∈ {0, N}.
  LD = β: Δ = β/L, run 0 = {Δ, …, β} meets and misses I: contradiction. LD = -α: λ = α/L, Δ = N - λ.
  For 1 ≤ t ≤ λ-1, ρ = tβ mod α, i = t-1 + ⌊tβ/α⌋ ≤ d-2: Y((i+1)L) = ρ ∈ I, Y(iL+1) = ρ + α - λ ⇒ ρ ≤ λ-1.
  λ ≥ 2 ⇒ r = β mod α ∈ [1, λ-1]; r ≥ 2: t = ⌈λ/r⌉ gives ρ = tr ∈ [λ, 2λ-2]: contradiction. So λ = 1 or r = 1 ⇒ α ∣ n.

## Lemma 7 (involution form), proof used
f̄² = id ⇒ f² = μ·id, σ² = id (d ≥ 2). σ = id: fixed points ⊆ ℙ(V_λ1) ∪ ℙ(V_λ2) (λ² = μ), V_λ1 ∩ V_λ2 = 0 ⇒ |V_λ1||V_λ2| ≤ q^d,
each ≤ q^{d-1}: (q-1)|Fix| ≤ q^{d-1}+q-2 (q ≥ 3), ≤ 2^{d-1}-1 (q = 2). σ ≠ id: F₀ = Fix σ, |F| = |F₀|² (x = a + bω);
W = {v : f v = λ₀v}: Hilbert 90 (c = y + ν⁻¹σ y) puts every fixed point in W; (a, b) ↦ a + ωb injective ⇒ |W| ≤ q₀^d;
(q₀-1)|Fix| ≤ |W| - 1. Both ≤ n/2 = (q^d-1)/(2(q-1)) unless (d, q) = (2, 4).

## A_n ≤ G (B.4)
Primitive (Prop 5) + n-cycle P + ¬(A_n ≤ G) ⇒ Jones (a) n prime, G affine: r_k ≠ 1 fixes ≥ 2 points, impossible;
(b) G ≤ PΓL_d(q): r_k involution ≠ 1 fixes m+1 > n/2 points, contradicts Lemma 7 ((d,q) = (2,4) ⇒ n = 5); (c) n ∈ {11, 23}: kernel G = ⊤.

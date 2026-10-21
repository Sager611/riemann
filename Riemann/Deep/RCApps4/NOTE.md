# RCApps4: the family large sieve for indefinite forms (D > 0) and for sectors; both thresholds drop from 3/2 to (3+√69)/10

Date 2026-10-21. Brief: `.scratch/orca/rcapps4-brief.md`. Scratch: `.scratch/rcapps4/` (`qf.py` helpers; `units.gp` → `units.out`
(trace t₀ of the totally positive fundamental unit of every real quadratic order of discriminant ≤ 4000, PARI); `dualpos.py`,
`dualpos_mut.py` (Lemma 2⁺, §1); `smoothpos.py` (the unfolded identity, §1); `poispos.py` (Lemma 6⁺); `dualang.py` (Lemma 2^∠);
`gxl4.py` (novelty); outputs `*.out`). References:
[C] = [Combos NOTE](../Combos/NOTE.md), [D2] = [Density2 NOTE](../Density2/NOTE.md), [R] = [RCApps NOTE](../RCApps/NOTE.md),
[R2] = [RCApps2 NOTE](../RCApps2/NOTE.md), [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md), [RS] = [RSFam NOTE](../RSFam/NOTE.md).
Paper proofs, not refereed, no Lean, no 2026 half-planes. Notation of [R3] §0 throughout (admissible n, Φ, a_V, P, δ, 𝒟 …).

## Verdicts

| Item | Result | Status |
|---|---|---|
| goal 1, tool | Theorem 1⁺: [R3] Theorem 1 for narrow (ring) class characters of real quadratic orders, Σ_{D∈𝔇⁺(X)}Σ_{χ∈Ĝ⁺(D)}\|Σ_n b(n)λ_χ(n)\|² ≪ X^{O(δ)}(X^{3/2} + N³X^{−1})N^{1−2σ}, N ≤ X^{1−4δ} | new; the obstruction named in [R3] §6 (units) is removed by unfolding the unit group with Hecke Grössencharacters (Lemma 2⁺) |
| goal 1, density | Theorem 2⁺: Σ_{D∈𝔇⁺(X)}N^♮_D(σ;t) ≪ X^{(13−8σ)/(7−2σ)+ε}, ½+ε ≤ σ ≤ 9/10 | new |
| goal 1, primes | Theorem 3⁺: [R2] Theorem 4 and Corollary 5 hold for every A > A₀ = (3+√69)/10 = 1.130662…: for all but O(X^{1−c}) fundamental D ∈ (0, X], every narrow class has prime ideals of norm in (x, x+xD^{−κ}] for every x ≥ D^A, and every primitive indefinite form of discriminant D represents a prime ≤ 2D^A; a² − nb² for almost all n | **new record** (was 3/2, [R2]) |
| goal 2, sectors | Theorems 1^∠–3^∠: [R] Theorem 3, [R2] Theorem 3′ and Corollary 4′ (primes p = a² + nb² in (x, x+x(4n)^{−κ}] with arg(a+ib√n) in any arc of length ≥ (4n)^{−κ}) hold for every A > (3+√69)/10 | **new record** (was 3/2, [R2]) |
| goal 3, a_E signs | no change (A > 3, one-sided A > 9/4): at the binding points the detector needs Dirichlet polynomials longer than X or second moments of degree-4/6 L-functions; the family sieve sees neither (§3) | NOGO with these tools |
| goal 4, sharper sieve | no change: the N³X^{−1} term is the family sieve beyond the conductor at scale X/w²; both its Poisson main term and its Poisson error are of size N³X^{−3/2} there, and the error needs cancellation in sums over y of Salié-type sums (§4) | open |

## 1. Positive discriminants: unfolding the units

Let Δ > 0 be a non-square discriminant, O_Δ = Z + fO_K (Δ = f²d_K), G⁺_Δ := I_K(f)/P⁺_{K,Z}(f) (P⁺_{K,Z}(f): the (α) with
α ≫ 0 and α ≡ a mod^× fO_K for some a ∈ Z prime to f), h⁺(Δ) := |G⁺_Δ|. For an integral 𝔞 ∈ P⁺_{K,Z}(f) the generators
γ ≫ 0 with γ ∈ O_Δ form one orbit of O_Δ^{×,+} = η_Δ^Z (η_Δ > 1 the totally positive fundamental unit; η_Δ = (t₀ + u₀√Δ)/2 with
(t₀, u₀) the least positive solution of t² − Δu² = 4); R⁺_Δ := log η_Δ. For γ ≫ 0 put u_γ := ½log(γ/γ′); then
u_{ηγ} = u_γ + R⁺_Δ, and for γ = (T + y√Δ)/2 with Nγ = M: T = 2√M cosh u_γ, y√Δ = 2√M sinh u_γ.
Fix an even ν ∈ C_c^∞([−½, ½]), ν ≥ 0, ν ≢ 0, and put w := ν∗ν: w ≥ 0 is even, supp w ⊂ [−1, 1], ŵ = ν̂² ≥ 0 and ŵ(0) = (∫ν)² > 0.

**Lemma 1⁺ (Grössencharacters).** For k ∈ Z there is a character λ_k of I_K(f) with λ_k((γ)) = e(k u_γ/R⁺_Δ) for every
(γ) ∈ P⁺_{K,Z}(f) (γ ∈ O_Δ, γ ≫ 0). It satisfies λ_k((m)) = 1 for m ∈ Z_{>0} prime to f and λ_k(𝔞̄) = λ̄_k(𝔞).

*Proof.* (γ) ↦ e(k u_γ/R⁺) is well defined (u changes by multiples of R⁺ along the orbit) and multiplicative on P⁺_{K,Z}(f) (the good
generator of a product is the product of good generators), so it is a homomorphism to S¹; S¹ is divisible, so it extends to the free
abelian group I_K(f). The generator of (m) is m, with u = 0. Finally λ_k(𝔞)λ_k(𝔞̄) = λ_k((N𝔞)) = 1 and |λ_k| = 1. □

**Lemma 2⁺ (unfolded dual identity).** Let b be supported on admissible n, and for M ≥ 1 put
𝒲_Δ(M) := 1[(M,f)=1]·Σ_{γ∈O_Δ, γ≫0, Nγ=M} w(u_γ) and Q_Δ(b) := Σ_{n,n′}b(n)b̄(n′)Σ_{c|(n,n′)}(Δ/c)𝒲_Δ(nn′/c²). Then
(i) Σ_{k∈Z}ŵ(k/R⁺_Δ)·Σ_{C∈G⁺_Δ}|Σ_{𝔞∈C,(𝔞,f)=1}b(N𝔞)λ_k(𝔞)|² = R⁺_Δ·Q_Δ(b); in particular Q_Δ(b) ≥ 0 and
  Σ_{χ∈Ĝ⁺_Δ}|Σ_n b(n)λ_χ(n)|² ≤ ŵ(0)^{−1}·h⁺(Δ)R⁺_Δ·Q_Δ(b);
(ii) h⁺(Δ)R⁺_Δ = √Δ·L(1,χ_{d_K})·Π_{p|f}(1 − (d_K/p)p^{−1}) ≤ C√Δ·log(3Δ)·log log(9Δ);
(iii) 𝒬⁺(b) := Σ_{Δ>0 non-square}Φ(Δ/X)Q_Δ(b) = Σ_{n,n′}b(n)b̄(n′)Σ_{c,w: cw|(n,n′), (c,w)=1}μ(w)𝒩⁺_c(nn′/(cw)²; X/w²), where
  𝒩⁺_c(M; X′) := Σ Φ(Δ/X′)(Δ/c)·w(arccosh(T/(2√M))) over T ∈ Z_{>0}, y ∈ Z and non-square Δ > 0 with T² − Δy² = 4M, T ≡ yΔ (2).

*Proof.* (i) Parseval on the coset λ_kĜ⁺_Δ gives Σ_{χ}|Σ_n b(n)λ_{χλ_k}(n)|² = h⁺Σ_C|S_C^{(k)}|², S_C^{(k)} the inner sum. Put
T_n^{(k)} := Σ_{N𝔞=n,(𝔞,f)=1}λ_k(𝔞)[𝔞] ∈ C[G⁺_Δ]. The proof of [R3] Lemma 1 applies verbatim: for split p = 𝔭𝔭̄,
(T_p^{(k)})² = λ_k(𝔭)²[𝔭²] + λ_k(𝔭̄)²[𝔭̄²] + 2λ_k((p))[(p)] and T_{p²}^{(k)} = λ_k(𝔭²)[𝔭²] + λ_k(𝔭̄²)[𝔭̄²] + λ_k((p))[(p)], with
λ_k((p)) = 1 and [(p)] = 1 in G⁺_Δ (p ≫ 0), so (T_p^{(k)})² = T_{p²}^{(k)} + 1; inert (T_p^{(k)} = 0, T_{p²}^{(k)} = 1) and ramified
((T_p^{(k)})² = 1 = T_{p²}^{(k)}) alike. So T_n^{(k)}T_{n′}^{(k)} = Σ_{c|g}(Δ/c)T^{(k)}_{mm′(g/c)²}. By Lemma 1⁺, (S^{(k)})* = Σ_n b̄(n)T_n^{(k)}, so
Σ_C|S_C^{(k)}|² = Σ_{n,n′}b(n)b̄(n′)Σ_{c|g}(Δ/c)[1](T^{(k)}_{M_c}), M_c := nn′/c², and [1](T^{(k)}_M) = 1[(M,f)=1]Σ_{γ mod η^Z}
e(k u_γ/R⁺) (the identity class of G⁺_Δ is P⁺_{K,Z}(f); an integral 𝔞 in it is prime to f iff (M, f) = 1, as in [R3] Lemma 2).
Poisson summation: Σ_k ŵ(k/R⁺)e(ku/R⁺) = R⁺Σ_{j∈Z}w(u − jR⁺); summed over the orbit representatives this is R⁺ times the sum of
w(u_γ) over all γ ≫ 0 with Nγ = M. This is (i). Since ŵ ≥ 0, the k = 0 term is ≤ R⁺Q_Δ, and h⁺Σ_C|S_C^{(0)}|² is the left side
of the inequality. (ii) Dirichlet's formula h(O)R(O) = h(O_K)R(O_K)·fΠ_{p|f}(1 − (d_K/p)/p), h(O_K)R(O_K) = ½√d_K L(1,χ_{d_K}),
and h⁺R⁺ = 2hR (if O_Δ has a unit ε of norm −1 then h⁺ = h, η = ε²; otherwise h⁺ = 2h, η = ε); L(1,χ) ≤ 2 + log d_K.
(iii) As [R3] Lemma 2(ii): if p ∥ M_c and p | f there is no γ; 1[(r,f)=1] = Σ_{w|r}μ(w)1[w²|Δ] for r = g/c; the map
(T, y, Δ) ↦ (T/w, y, Δ/w²) is a bijection onto the solutions for M_c/w², preserves (Δ/c), Φ(Δ/X) = Φ((Δ/w²)/(X/w²)) and
T/(2√M) (hence the weight: w(u_γ) = w(arccosh(T/(2√M))) since w is even); γ ≫ 0 ⟺ T > 0, as γγ′ = M > 0. □

*Checks.* `dualpos.py`: for all 534 real quadratic orders with Δ ≤ 4000 and t₀ < 3000 and 16 pairs (n, n′) (incl. n = n′ and
gcd 15, 105), the number of pairs of ideals of norms n, n′ in one narrow class (proper equivalence of indefinite forms, cycles of
reduced forms) equals Σ_{c|g}(Δ/c)·#{(T,y): T > 0, T² − Δy² = 4nn′/c², T² < (nn′/c²)(t₀+2)} (the orbits of η^Z; boundary ½):
8544 of 8544 exact. `dualpos_mut.py`: dropping the Jacobi factor, counting both signs of T, or the domain T² < M(t₀−2) fails
(5615, 5526, 6906 of 6942).

**Theorem 1⁺.** [R3] Theorem 1 holds with 𝔇⁺(X), Ĝ⁺(D) in place of 𝔇(X), Ĝ(D): for 0 < δ ≤ 1/100, X ≥ X₀(δ), ½+2δ ≤ σ ≤ 9/10,
|t′| ≤ X^δ, 2 ≤ V ≤ X^{1/2}, N = 2^j ∈ [V/4, X^{1−4δ}] and b ∈ {b_{j,t′}, ∂_{t′}b_{j,t′}} of [R3] §2,
  Σ_{D∈𝔇⁺(X)}Σ_{χ∈Ĝ⁺(D)}|Σ_n b(n)λ_χ(n)|² ≤ C(δ)X^{40δ}(X^{3/2} + N³X^{−1})N^{1−2σ}.

*Proof.* Every D ∈ 𝔇⁺(X) is a non-square Δ with f = 1 and Φ(D/X) = 1. By Lemma 2⁺(i) (applied to every Δ, all terms ≥ 0) and
(ii), LHS ≤ ŵ(0)^{−1}Σ_ΔΦ(Δ/X)h⁺R⁺Q_Δ(b) ≤ CX^{1/2}log²X·𝒬⁺(b). By (iii), 𝒬⁺(b) = Σ_{c,w}μ(w)Q⁺_{c,w}, Q⁺_{c,w} :=
Σ_{n₂,n₂′}β(n₂)β̄(n₂′)𝒩⁺_c(n₂n₂′; X/w²), exactly as in [R3] §3.1. It remains to show Σ_{c,w}|Q⁺_{c,w}| ≪ X^{O(δ)}(X + N³X^{−3/2})‖b‖².
[R3] §3.2–§3.8 go through with the following changes; everything else (§3.3 large w, the local densities A_M(y) and Lemma 7, §3.6
secondary term, §3.7 with Lemma 9, §3.8) is used verbatim, because it depends only on b, on the congruence conditions
y² | T² − 4M, (T² − 4M)/y² ≡ 0, 1 (4), on (Δ/c) mod c, and on Lemmas 5–6, none of which sees the sign of Δ.
- *Lemma 5⁺.* y = 0 forces T = 2√M: ≤ 4X′ as before. y ≠ 0: on the support of the weight, |u_γ| ≤ 1, so 2√M < T ≤ 2√M cosh 1, and
  Φ(Δ/X′) ≠ 0 gives y² ≤ 4(T² − 4M)/X′ ≤ 23M/X′. For fixed T, y² | T² − 4M ≠ 0; for fixed y, T lies in ≤ 2^{ω(y)+2}(y²,4M)^{1/2}
  classes mod y² and in an interval of length ≤ y²X′/(2√M). The bound C(ε)M^ε min{√M + 1, √X′ + √(M/X′)} follows as in [R3].
- *Squares.* (iii) excludes square Δ, the Poisson count below includes them. For y ≠ 0 and Δ = a², (T − ay)(T + ay) = 4M leaves
  ≪ M^ε triples; their total contribution is ≪ X^{δ}N‖b‖² ≤ X^{1+δ}‖b‖². (y = 0 is in the first item.)
- *Lemma 6⁺.* Let c be odd squarefree, (c, M) = 1, X′ ≥ 128c√M X^δ, M ≤ X². Then the y ≠ 0 part of 𝒩⁺_c(M; X′) (all Δ ≡ 0, 1 (4))
  equals (μ(c)/c)(X′/(2√M))Σ_{y≥1,(y,c)=1}A_M(y)K⁺(y²X′/(4M)) + O_δ(X^{−10}), with
  K⁺(v) := ∫_0^∞Φ(u)W₁(vu)(1 + vu)^{−1/2}du, W₁(s) := w(arcsinh√s).
  *Proof.* As [R3] Lemma 6 with f(t) := Φ((t² − 4M)/(y²X′))W(t), W(t) := w(arccosh(t/(2√M))) for t > 2√M, f(t) := 0 for t ≤ 2√M.
  W is smooth: w(x) = w₀(x²) with w₀ smooth (w even), and arccosh(τ)² is real-analytic on [1, ∞). On supp f, t ≤ 2√M cosh 1 and
  y²X′ ≤ 23M, so with λ₁ := 8√M/(y²X′): |u′| ≤ λ₁, |u″| ≤ λ₁², u‴ = 0 for u(t) := (t² − 4M)/(y²X′), and |W^{(i)}| ≪_i (2√M)^{−i}
  ≤ (2λ₁)^i. Hence ‖f^{(i)}‖₁ ≪_i √M(3λ₁)^i and |f̂(k/q)| ≪_i √M(96c√M/(|k|X′))^i for q = 4y²c, which is ≪ √M(X^{−δ}/|k|)^i.
  Main term: t = 2√M(1+s)^{1/2} gives ∫f = √M∫_0^∞Φ(s/v)W₁(s)(1+s)^{−1/2}ds = √M·v·K⁺(v), v = y²X′/(4M); the local sums are those
  of [R3] Lemma 6 (they do not involve the sign of Δ); y and −y give the same T. □
  On the ranges of [R3] (w ≤ G₂ := X^{1−2δ}/(512N)) the hypothesis holds; the factor 512 instead of 64 only changes constants in §3.3.
- *§3.5.* g⁺(x) := K⁺(x²) is even, smooth (W₁ is smooth on [0, ∞) since arcsinh(√s)² is analytic at s = 0), supported in
  |x| ≤ 2 sinh 1; Š⁺(z) := Σ_{k odd}ǧ⁺… = (1 − 2^{−z})ζ(z)ǧ⁺(z) has the same poles and decay as Š. The constant 2c₁ becomes ǧ⁺(1).
Collecting: Σ_{c,w}|Q⁺_{c,w}| ≪ X^{12δ}(X + N³X^{−3/2})‖b‖², and the theorem follows. □

**Theorem 2⁺.** For 0 < δ ≤ 1/100, ½ + 3δ ≤ σ ≤ 9/10, |t| + 3 ≤ X^δ/4:
Σ_{D∈𝔇⁺(X)}N^♮_D(σ;t) ≤ C(δ)X^{(13−8σ)/(7−2σ)+60δ}, N^♮_D counting zeros with σ ≤ β ≤ 9/10, |γ − t| ≤ 1 of all χ ∈ Ĝ⁺(D)∖{1}.

*Proof.* [R3] §2 and §4 with χ ∈ Ĝ⁺(D)∖{1}. Lemma 3: for p ∤ D, χ(𝔭)χ(𝔭̄) = χ((p)) = 1 (p ≫ 0) and χ((p)) = 1 at inert p, so
L_p(s,χ)^{−1} = 1 − λ_χ(p)p^{−s} + (D/p)p^{−2s} as before; (c) uses [D2] Lemma 5 (narrow ring class large sieve, (N + D^{1/2})
log²(3D)) in place of [C] Lemma 2. Lemma 4: L(s,χ) is entire of polynomial growth for χ ≠ 1. Class II: the mean value
Σ_{χ∈Ĝ⁺∖{1}}|L(½+δ+iv,χ)|² ≤ C(δ)D^{1/2+δ}(|v|+3)^{1+δ} is [D2] (2.1) for 𝒵 = Ĝ⁺, mollifier length 1 (as in [R2] Prop 2;
the two archimedean types are treated separately). Class I: Theorem 1⁺ for |τ| ≤ X^δ/2, [D2] Lemma 5 otherwise. □

**Theorem 3⁺.** Let A₀ := (3+√69)/10. [R2] Theorem 4 holds for every A > A₀ (same c, κ, δ as in [R3] Theorem 3), and so does
[R2] Corollary 5: for all but O(X^{1−c}) fundamental D ∈ (0, X], every primitive indefinite binary quadratic form of discriminant D
represents ≫ ν_D H/(h⁺ log x) primes in every (x, x+H], x ≥ D^A, H ≥ xD^{−κ}, in particular a prime ≤ 2D^A; and for all but
O(N^{1−c}) squarefree n ≤ N with n ≡ 2, 3 (mod 4), the least prime a² − nb² is ≤ 2(4n)^A.

*Proof.* [R3] Theorem 3 (its four changes to [R] §2) under the dictionary of [R2] Theorem 4 (G⁺, h⁺, Ĝ⁺_c, Ĝ⁺₂, narrow genus theory,
(I6) for Hecke characters of K): step 2 uses Theorem 2⁺ and [R2] Proposition 2 (F₂(σ) = 10(1−σ)/(3−2σ) for D > 0 as well), so
the functions F₁, F₂, F*, φ and the crossing σ_c = 5 − √69/2 are literally those of [R3]; step 3 uses
h⁺ ≤ 5D^{1/2}log²(3D) per window; step 4 and [R] §2.5 as in [R2] Theorem 4 (Σ|R| ≤ C(h⁺ + 2^{ω(D)})(H/x + x^{−1/2})log(Dx)
≤ HD^{1/2−A+o(1)}; prime powers, inert and ramified primes need x ≥ D^{1+2δ}). The corollary is [R2] Corollary 5. □

## 2. Sectors (D < 0)

Let Δ < 0, |Δ| > 4 (so O_Δ^× = {±1}) and ℓ ∈ 2Z. The map (γ) ↦ (γ/|γ|)^ℓ (γ ∈ O_Δ, the generator of (γ) ∈ P_{K,Z}(f), unique
up to sign) is a homomorphism P_{K,Z}(f) → S¹; extend it to a character ψ_ℓ of I_K(f) (S¹ divisible). For f = 1 the coset
ψ_ℓĜ(D) is the set of Hecke characters of ∞-type ℓ used in [R] Theorem 3 (𝒜(D, M) = ∪_{|ℓ|≤M}ψ_ℓĜ(D)). As in Lemma 1⁺,
ψ_ℓ((m)) = 1 for m ∈ Z_{>0} and ψ_ℓ(𝔞̄) = ψ̄_ℓ(𝔞). T_ℓ is the Chebyshev polynomial, T_ℓ(cos θ) = cos ℓθ.

**Lemma 2^∠.** For b supported on admissible n and X ≥ 32,
(i) Σ_{ψ∈ψ_ℓĜ_Δ}|Σ_n b(n)λ_ψ(n)|² = h(Δ)Σ_C|S_C^∠|², S_C^∠ := Σ_{𝔞∈C,(𝔞,f)=1}b(N𝔞)ψ_ℓ(𝔞);
(ii) Σ_{Δ∈𝒟,|Δ|>4}Φ(|Δ|/X)Σ_C|S_C^∠|² = ½Σ_{n,n′}b(n)b̄(n′)Σ_{c,w}μ(w)𝒩^∠_c(nn′/(cw)²; X/w²), with
  𝒩^∠_c(M; X′) := Σ_{(T,y,Δ): T²−Δy²=4M, T≡yΔ (2)}Φ(|Δ|/X′)(Δ/c)T_ℓ(T/(2√M)).

*Proof.* [R3] Lemma 2 with T_n^∠ := Σ_{N𝔞=n,(𝔞,f)=1}ψ_ℓ(𝔞)[𝔞]: [R3] Lemma 1 holds for T^∠ since ψ_ℓ((p)) = 1;
(S^∠)* = Σ b̄(n)T_n^∠; [1](T^∠_M) = ½1[(M,f)=1]Σ_{γ∈O_Δ, Nγ=M}(γ/|γ|)^ℓ, and γ = (T + iy√|Δ|)/2 = √M e^{iθ} with
cos θ = T/(2√M); the solution set is symmetric under y ↦ −y, so (γ/|γ|)^ℓ may be replaced by cos ℓθ = T_ℓ(T/(2√M)). The
bijection (T, y, Δ) ↦ (T/w, y, Δ/w²) maps γ to γ/w and keeps the angle. □

**Theorem 1^∠.** [R3] Theorem 1 holds for Σ_{D∈𝔇(X)}Σ_{ψ∈ψ_ℓĜ(D)}|Σ_n b(n)λ_ψ(n)|², uniformly for even |ℓ| ≤ X^δ, with
X^{40δ} replaced by X^{60δ}.

*Proof.* Positivity and Lemma 2^∠ as in [R3] §3.1. |T_ℓ| ≤ 1 on [−1, 1], so Lemma 5 is unchanged. In Lemma 6 the function is
f_ℓ(t) := Φ((4M − t²)/(y²X′))T_ℓ(t/(2√M)); by Markov's inequality |T_ℓ^{(i)}| ≤ ℓ^{2i} on [−1, 1], and (2√M)^{−1} ≤ 2λ₁
(λ₁ := 4√M/(y²X′), y²X′ ≤ 16M), so ‖f_ℓ^{(i)}‖₁ ≪_i √M((1 + 2ℓ²)λ₁)^i; the hypothesis becomes X′ ≥ 64c(1 + 2ℓ²)√M X^δ, which
holds for w ≤ G₂ := X^{1−4δ}/(768N) (ℓ² ≤ X^{2δ}); §3.3 with this G₂ costs X^{O(δ)}. The main term has
K_ℓ(v) := ∫Φ(u)(1 − vu)^{−1/2}T_ℓ((1 − vu)^{1/2})1[vu < 1]du in place of K; ℓ even makes T_ℓ(√(1−s)) =: P_ℓ(s) a polynomial of
degree ℓ/2 with |P_ℓ^{(i)}| ≤ (ℓ²/2)^i on [0,1] (Markov on [0,1]), so g_ℓ(x) := K_ℓ(x²) is smooth, even, supported in |x| ≤ 2,
with |g_ℓ^{(i)}| ≪_i (1+ℓ²)^i. Two integrations by parts give |ǧ_ℓ(z)| ≪ (1+ℓ²)²(1+|z|)^{−2} on −δ ≤ Re z ≤ 2, |z| ≥ δ; so
Š_ℓ = (1 − 2^{−z})ζ(z)ǧ_ℓ(z) has the same pole structure, ∫_{(−δ)}|Š_ℓ(z)||dz| ≪ C(δ)(1+ℓ²)² and |ǧ_ℓ(1)| ≪ 1. These
factors are ≤ X^{4δ}. §3.6–§3.8 are unchanged (they see only b and the local densities). □

**Theorem 2^∠.** For even 0 < |ℓ| ≤ X^δ and the ranges of [R3] Theorem 2:
Σ_{D∈𝔇(X)}Σ_{ψ∈ψ_ℓĜ(D)}N^♮(σ;t,ψ) ≤ C(δ)X^{(13−8σ)/(7−2σ)+80δ}.

*Proof.* [R3] §2, §4 with ψ for χ. L(s, ψ) is entire (ℓ ≠ 0) and L_p(s,ψ)^{−1} = 1 − λ_ψ(p)p^{−s} + (D/p)p^{−2s} (ψ((p)) = 1);
Lemma 3(c) is [C] Lemma 2 applied to the coefficients c_𝔞ψ_ℓ(𝔞) (Parseval over the coset); class II uses
Σ_{ψ∈ψ_ℓĜ}|L(½+δ+iv,ψ)|² ≤ C(δ)(|D|(|ℓ|+3)²)^{1/2+δ}(|v|+3)^{1+δ} ([D2] (2.1) for 𝒜(D, |ℓ|), as in [R] Proposition 1′), costing
(|ℓ|+3)^{1/2+δ} ≤ X^{δ}; class I uses Theorem 1^∠ (|τ| ≤ X^δ/2) and [C] Lemma 2 otherwise. □

**Theorem 3^∠.** [R] Theorem 3, [R2] Theorem 3′ and [R2] Corollary 4′ hold for every A > A₀ = (3+√69)/10: for all but
O(X^{1−c}) D ∈ 𝔇 with |D| ≤ X, every class C, every arc I ⊂ R/(2π/w)Z of length θ ≥ |D|^{−κ}, every x ≥ |D|^A and
H ∈ [x|D|^{−κ}, x] have ≥ wθν_D H/(128πh log 2x) degree-one 𝔭 ∈ C with N𝔭 ∈ (x, x+H] and arg β_𝔭 ∈ I; for all but O(N^{1−c})
squarefree n ≤ N, n ≡ 1, 2 (4), every (x, x + x(4n)^{−κ}], x ≥ (4n)^A, contains a prime a² + nb² with arg(a + ib√n) in any
prescribed arc of length ≥ (4n)^{−κ}.

*Proof.* [R] Theorem 3 (i)–(v) with the changes of [R3] Theorem 3 (parameters, Markov with F*, β < σ₀, inherited errors) in
place of those of [R2]. In (iv), for ψ ∈ 𝒜(D, M)∖Ĝ, M = 6|D|^{κ+ε₂}: zeros with β ≤ 9/10 are counted by Theorem 2^∠ summed over
the ≤ 2M_X + 1 values of ℓ (M_X := 6X^{κ+ε₂} ≤ X^{δ}), zeros with β > 9/10 by [R] Proposition 1′ (exponent
F₂(σ) + O(κ + ε₂)); so the family exponent is F*(σ) + O(δ) and the Markov step of [R3] Theorem 3 step 2 holds with ε′ as there
(δ ≤ ε′/100). Zeros with β < σ₀ = ½ + 4δ: each unit window holds ≤ C(2M + 1)h log(|D|M(|k|+3)) zeros of all ψ ∈ 𝒜(D, M), so
their contribution relative to ĝ(0)H/h is ≤ |D|^{1/2 − A(1/2−4δ) + 2κ + ε₁ + ε₂ + o(1)} ≤ |D|^{−(A−1)/2 + 6δ}. The rest of (iv)–(v)
is unchanged. The corollary is [R2] Corollary 4′ (D = −4n). □

## 3. The a_E-sign theorem: why the family sieve does not lower 3 or 9/4

[RS] Theorem D needs, for m = 1 (one-sided, A > 9/4) and m = 2 (sign, A > 3), the family counts of zeros of
L_m(s,χ) = L(s, Sym^m f ⊗ θ_χ). Their binding points are σ = 5/6 for both.
- *A twisted family sieve exists for short polynomials.* With the squarefree detector of [R3] §2 applied to L_m (for squarefree n,
  the coefficient of L_m is λ_{Sym^m f}(n)λ_χ(n), and Σ_{n sqfree}λ_{Sym^m f}(n)λ_χ(n)n^{−s} = L_m(s,χ)·H_m with H_m absolutely
  convergent on Re s > ½ for admissible primes p > 16(m+1)²), the class I coefficients are b(n)λ_{Sym^m f}(n), independent of D, and
  the proof of [R3] Theorem 1 should go through with L(s, Sym^m f ⊗ ψ̃) in Lemma 9 (no pole). We do not write this out, because:
- *m = 1.* Class II per D is Y^{1/2−σ}(Σ_χ|L_1|²Σ_χ|M̃|²)^{1/2} ≤ Y^{1/2−σ}|D|^{3/4} ([RS] Lemma B), so X^{7/4}Y^{1/2−σ} over the
  family. With Y ≤ X (the range of the family sieve) this is ≥ X^{9/4−σ} = X^{1.4167} at σ = 5/6, above the bar 1 + (9/4)/6 = 1.375.
  With Y > X, class I for N ∈ [X, Y] only has the per-D sieve, X·Y^{2−2σ}; balancing with class II gives Y = X^{3/(2(3−2σ))},
  i.e. X times [RS] Theorem C exactly. So the best detector bound equals the old one at σ = 5/6 for every choice of Y.
- *m = 2.* Class II is X²Y^{1/2−σ}, ≥ X^{5/3} at σ = 5/6 for Y ≤ X, against the bar X^{3/2}; even the trivial count is X^{3/2} there
  ([RS] §7), so only a nontrivial family bound for degree-6 L-functions at σ = 5/6 would help.
- *What would help.* Either the family sieve for N up to X^{9/8} (m = 1) resp. X^{3/2} (m = 2), i.e. beyond the conductor X, or a
  family second moment Σ_DΣ_χ|L_m(½+δ+iv,χ)|² below X^{(m+3)/2} (AFE length X^{(m+1)/2} ≥ X, same obstruction).
Verdict: no improvement of A > 3, A > 9/4 with these tools.

## 4. Sharpening the family sieve: where N³X^{−1} comes from, and why it stays

In [R3] §3.3 the pairs n = cwn₂, n′ = cwn₂′ with w > G₂ ≍ X/N are family sums at scale X′ = X/w² with lengths L = N/(cw) ≥ X′/64,
beyond their conductor. For these the Poisson count in T splits into the k = 0 term (main term, ≍ √X′ per M, total
N√X/w² ·‖β‖²) and the k ≠ 0 terms (≍ #{y}·#{roots} ≍ √(M/X′) = N/(c√X) per M, total N²/(w√X)‖β‖²); at w ≍ X/N both are
N³X^{−3/2}. Replacing the lattice count by the per-Δ large sieve gives (L√X′ + X′)‖β‖² = (N√X/(cw²) + X/w²)‖β‖², the same
order. A gain needs cancellation, over y ≤ 4N/(c√X) and M = n₂n₂′, in Σ_{T² ≡ 4M (y²)}e(kT/(4y²c)) (Salié-type sums; half-integral
weight Kuznetsov, as in Duke's equidistribution of CM points), uniformly with the coefficients β. Not attempted further. If it worked,
A₀ would drop to 1 ([R3] §6); the numerics of [R3] §7 support the conjecture.

## 5. Numerics (`.scratch/rcapps4/`, every script under a hard alarm)
- Lemma 2⁺, orbit form (`dualpos.py` → `dualpos.out`): 534 real quadratic orders (Δ ≤ 4000, t₀ < 3000, incl. non-maximal), 16 pairs
  (n, n′) each: pairs of ideals in one narrow class (forms, cycles of reduced forms) = Σ_{c|g}(Δ/c)·#orbits, 8544/8544 exact. Mutations
  (`dualpos_mut.py` → `dualpos_mut.out`) fail as they should. Most tested orders have a unit of norm +1, so h⁺ = 2h is exercised.
- Lemma 2⁺(i), the unfolded identity (`smoothpos.py` → `smoothpos.out`; random complex b on 13 values of n, Gaussian w (a Schwartz
  w is allowed in (i)), Grössencharacter values from generators of 𝔞𝔠̄ found by lattice membership, all k with ŵ(kL/R) > e^{−36π},
  general width L ∈ {0.4, 1}: Σ_kŵ(kL/R)Σ_C|S_C^{(k)}|² = (R/L)Q_Δ): 30 orders, 60 tests, max relative deviation 8.9·10^{−16}.
  (A first run with L = 2.5 had a γ-range |u| ≤ 17.5, i.e. y up to ~10^7√M, and was killed after 4.7 h; L ≤ 1 is the case used.)
- Lemma 6⁺ (`poispos.py` → `poispos.out`, bump w, all c odd squarefree): brute force against the main term:
  (M, X′, c) = (392863, 2·10⁶, 1): 44.4781 vs 44.4781; c = 3: −14.826 vs −14.826; c = 35: 1.2687 vs 1.2708; (3065857, 3·10⁶, 15):
  19.3945 vs 19.3948; (10023013, 10⁵, 1): 77.4254 vs 77.447. Where X′/(c√M) is small (c = 2431, 145, 29: ratio ≤ 1.3) the two differ,
  as they should (outside the hypothesis of Lemma 6⁺).
- Lemma 2^∠ (`dualang.py` → `dualang.out`): 399 imaginary quadratic orders (101 ≤ |Δ| ≤ 899, incl. non-maximal), ℓ ∈ {2, 4, 10},
  random complex b on 13 values of n (incl. common factors 3, 5, 15, 105): 1197 tests, max relative deviation 6.7·10^{−15}.
- Thresholds: F₁, F₂ and A₀ are those of [R3] (`.scratch/rcapps3/thresh.py`); D > 0 and sectors change only O(δ)-margins.

## 6. Novelty (2026-10-21, short)
- Prior audits of this project found no average result for indefinite forms at all and nothing at 3/2 ([R2] NOVELTY, adversarial),
  and no family large sieve at the family size ([R3] NOVELTY). The thresholds here are below those.
- GXL evidence API, 3 claims (`gxl4.py` → `gxl4.out`, `gxl4_*.json`): indefinite forms for almost all D (only Ditchen 2013 (definite,
  20/3), Zaman 2018, Gaudet 2024 (imaginary)); a family large sieve over narrow class characters averaged over D (only class number
  distribution papers, Ditchen's sieve, and a 2026 Bombieri–Vinogradov theorem for sectors in one fixed real quadratic field, level 1/5);
  sectors a² + nb² for almost all n (Gaussian primes in narrow sectors, one field). Nothing states or implies Theorems 3⁺, 3^∠.
  Keyless web search returned nothing relevant.
- Method. Unfolding units with Grössencharacters is classical (Hecke; [D2] Lemma 5 also counts in a unit-fundamental domain via
  h⁺log η). New here is its use as a positivity device that turns the dual of the family sieve over D > 0 into a regulator-free
  lattice count of the same shape as for D < 0 (hyperbolic cosh in place of cos), so that [R3]'s Poisson/Mellin analysis applies.

## 7. Dead ends (NOGO candidates)
- a_E signs below 3 (and one-sided below 9/4) via the family sieve: class II (degree-4/6 second moments) or lengths beyond the
  conductor bind (§3).
- Sharpening [R3] Theorem 1 near N ≍ X: the N³X^{−1} term is the sieve beyond the conductor at scale X/w²; per-Δ large sieve gives the
  same order (§4); needs Salié-sum cancellation.
- D > 0 dual with a sharp fundamental domain (T² < M(t₀+2)): depends on t₀(Δ), erratic in Δ, kills Poisson in T; the smooth unfolding
  (Lemma 2⁺) avoids it.

# RCApps5: the almost-all least-prime exponent reaches the GRH limit: every A > 1

Date 2026-10-21. Brief: `.scratch/orca/rcapps5-brief.md`. Check scripts and their outputs: [checks/](checks/) (§8); the scratch
directory `.scratch/rcapps5/` (arXiv listings, search pages) was deleted. References: [C] = [Combos NOTE](../Combos/NOTE.md), [D2] = [Density2 NOTE](../Density2/NOTE.md),
[R] = [RCApps NOTE](../RCApps/NOTE.md), [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md), [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md),
[R4] = [RCApps4 NOTE](../RCApps4/NOTE.md). Paper proofs, not refereed, no Lean. Unconditional: no 2026 half-plane, no Kerr.
Notation of [R3] §0 throughout (𝒟, 𝔇(X), G_Δ, T_M, λ_χ, admissible n, Φ, a_V, δ, ‖b‖² := N^{1−2σ}).

## Verdicts

| Item | Result | Status |
|---|---|---|
| task 1: what fixed (3+√69)/10 | the Markov bar 1 + A(1−σ) against F₁ = (13−8σ)/(7−2σ) and F₂ = 10(1−σ)/(3−2σ), crossing at σ_c = 5 − √69/2 = 0.8467; F₁ is set by the N³X^{−1} term of [R3] Theorem 1, whose only source is the conductor Möbius inversion ([R3] §3.3, common factors w > X/N); no zero-free input is involved | §1 |
| new lemma | Lemma P: for every imaginary quadratic order, the pair count of all elements (no coprimality to the conductor f) is a sum over squarefree e \| f of class-group sums of the order of discriminant Δ/e², hence ≥ 0; it dominates the sum [R3] used | new; exact on 1830 (Δ, b) tests incl. 1449 non-maximal orders (§8) |
| tool | Theorem 1: Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}\|Σ_n b(n)λ_χ(n)\|² ≪ X^{3/2+O(δ)}N^{1−2σ} for V/4 ≤ N ≤ X^{1−4δ} ([R3] detector sequences) | new; [R3] had + N³X^{−1}; family size for every length below the conductor |
| density | Theorem 2: Σ_{D∈𝔇(X)}N^♮_D(σ;t) ≪ X^{2−σ+ε} for ½+ε ≤ σ ≤ 9/10, i.e. the density hypothesis of each D (\|D\|^{1−σ}), on average over D | new; [R3]: X^{(13−8σ)/(7−2σ)} (19/16 against 7/6 at σ = 5/6) |
| least primes, D < 0 | Theorem 3: [R] Theorem 2 (primes in every ring class in (x, x+x\|D\|^{−κ}], x ≥ \|D\|^A, all but O(X^{1−c}) D) holds for every A > 1; every positive definite primitive form of discriminant D represents a prime ≤ 2\|D\|^A | **new record** (was (3+√69)/10 = 1.1307, [R3]); 1 is what GRH gives |
| D > 0, sectors | Theorems 3⁺, 3^∠: [R2] Theorem 4, Corollary 5 (indefinite forms, a² − nb²) and [R] Theorem 3, [R2] Corollary 4′ (arcs, a² + nb²) for every A > 1 | **new record** (was 1.1307, [R4]) |
| other routes (task 2) | DFI/amplification, Heath-Brown refinements, fourth moments, per-D log-free/Huxley/hybrid inputs, Liu's half-plane: not needed, and none goes below 1; Kerr: removes the Siegel factor ν_D only | §7; Kerr part is a COROLLARY |
| A < 1 | blocked: zeros with β near ½ enter with \|x^{ρ−1}\| (h·x^{−1/2} ≤ 1 forces x ≥ \|D\|), the GRH limit of the explicit formula | NOGO |

## 1. What fixed A₀ = (3+√69)/10 (task 1)

The exceptional set E₁ of [R3] Theorem 3 (step 2) is controlled by Markov's inequality over D, which needs, for all σ ∈ [σ₀, 1),

  (M)  Σ_{D∈𝔇(X)} N^c_D(σ, T_X) ≤ X^{1 + A(1−σ) − η}.

Two family bounds were available:
- F₁(σ) = (13−8σ)/(7−2σ), σ ≤ 9/10 ([R3] Theorem 2). (M) holds iff A > 6/(7−2σ): increasing in σ.
- F₂(σ) = 10(1−σ)/(3−2σ), all σ ([R] Proposition 1: Humphries–Thorner detector, GL₂ large sieve threshold Q²|𝒮| ≍ X^{7/2}, second
  moment I₂ ≍ X^{3/2}). (M) holds iff A > 10/(3−2σ) − 1/(1−σ): decreasing on [0.8, 1).

The two thresholds are equal at σ_c = 5 − √69/2 = 0.846688 (the smaller root of 4σ² − 40σ + 31 = 0), with value (3+√69)/10 = 1.130662.
So A₀ is pinned by the σ-range around σ_c and by F₁ just below it.

Inside F₁ ([R3] §4): class II costs X^{3/2}Y^{1/2−σ} (per-D mean squares of L(½+δ+it, χ) and of the mollifier, both ≍ h, sharp);
class I costs X^{3/2}V^{1−2σ} + Y^{4−2σ}X^{−1}. The last term is the N³X^{−1} error of [R3] Theorem 1 at N = Y. Balancing gives
Y = X^{5/(7−2σ)} = X^{0.9422} at σ_c. The N³X^{−1} has exactly one source, [R3] §3.3: the dual terms μ(w)𝒩_c(nn′/(cw)²; X/w²) with
w > G₂ = X^{1−2δ}/(64N), where the T-intervals are shorter than the modulus and Poisson is impossible. These w-terms come from the
inversion 1[(r,f)=1] = Σ_{w|r}μ(w)1[w²|Δ], which [R3] Lemma 2 needs because it adds, for every non-maximal order, the class-group sum over
ideals prime to the conductor f.

Zero-free input: none. The cutoff β ≤ 9/10 (needed by the pole killer F(s) = 1 − 4^{1−s}) is not binding (σ_c < 9/10, F₂(9/10) = 5/6 < 1).
Siegel zeros enter only the factor ν_D and the set E₂ of [R], never the exponent.

So F₁ was limited by the choice of positive terms added in [R3] §3.1. §2 adds different ones.

## 2. Lemma P: positivity without coprimality to the conductor

For Δ ∈ 𝒟 (a negative discriminant, any order) put r_Δ(M) := ½#{γ ∈ O_Δ : Nγ = M} and, for b supported on admissible n,

  Q′_Δ(b) := Σ_{n,n′} b(n) b̄(n′) Σ_{c|(n,n′)} (Δ/c) r_Δ(nn′/c²).

Compared with [R3] Lemma 2, the coprimality factor 1[(nn′/c², f) = 1] is simply dropped.

**Lemma P.** Let Δ ∈ 𝒟 have conductor f, and let b be supported on admissible n. For squarefree e | f put Δ_e := Δ/e² (a
discriminant of conductor f_e := f/e), w_e := |O_{Δ_e}^×| and b_e(m) := b(em)·1[(m, f) = 1]. Then

  Q′_Δ(b) = Σ_{e | f squarefree} (w_e/2)·Σ_{C ∈ G_{Δ_e}} |S_C^{(Δ_e)}(b_e)|²  ≥ 0,   S_C^{(Δ_e)}(β) := Σ_{𝔞∈C, (𝔞, f_e)=1} β(N𝔞).

The term e = 1 is Σ_{C∈G_Δ}|S_C|², the quantity [R3] Lemma 2 sums. For fundamental D with |D| > 4, Q′_D(b) = Σ_C|S_C|².

*Proof.* (a) *Only equal f-parts survive.* Let p | f with p | nn′ for some n, n′ in the support; p is odd. Then p² | Δ, so (Δ/p) = 0,
and for γ = (T + y√Δ)/2 ∈ O_Δ with Nγ = M we have T² = 4M + Δy² ≡ 4M (mod p²), which is impossible if p ∥ M (v_p(T²) is even). In
Q′_Δ(b) only c with (c, f) = 1 contribute, and for them v_p(nn′/c²) = v_p(n) + v_p(n′) ∈ {0, 1, 2} for p | f. Hence a pair (n, n′)
contributes only if (n, f) = (n′, f) =: e. Write n = em, n′ = em′; then (mm′, f) = 1, c | (m, m′), and (Δ/c) = (Δ_e/c).
(b) *Division by e.* Δ_e ≡ Δ (mod 8) since e is odd, so Δ_e is a discriminant, of conductor f/e. If γ = (T + y√Δ)/2 ∈ O_Δ and
Nγ = e²M′, then for p | e both 4e²M′ and Δy² are divisible by p², so p | T; thus e | T and γ/e = (T/e + y√Δ_e)/2 with
T/e ≡ T ≡ yΔ ≡ yΔ_e (mod 2), i.e. γ/e ∈ O_{Δ_e}. Conversely eγ′ ∈ O_Δ for γ′ ∈ O_{Δ_e}. So γ ↦ γ/e is a bijection onto
{γ′ ∈ O_{Δ_e} : Nγ′ = M′}, and r_Δ(e²M′) = r_{Δ_e}(M′) (counts of elements; units play no role).
(c) *The blocks.* By (a), (b), Q′_Δ(b) = Σ_e B_e with B_e := Σ_{m,m′} b_e(m) b̄_e(m′) Σ_{c|(m,m′)} (Δ_e/c) r_{Δ_e}(mm′/c²). For M prime to
f_e, an ideal of P_{K,Z}(f_e) prime to f_e of norm M is (γ) with γ ∈ O_{Δ_e}, Nγ = M, and its generators in O_{Δ_e} form one orbit of
O_K^× ∩ O_{Δ_e} = O_{Δ_e}^× (if γ, uγ ∈ O_{Δ_e} with (Nγ, f_e) = 1 then u ≡ a′/a mod f_eO_K with a, a′ ∈ Z prime to f_e, so u ∈ O_{Δ_e});
conversely every such γ generates an ideal of P_{K,Z}(f_e). So r_{Δ_e}(M) = (w_e/2)[1](T_M) in Z[G_{Δ_e}]. [R3] Lemma 1 (whose
proof uses only [(p)] = 1 in G_{Δ_e}, valid for every order) gives Σ_{c|(m,m′)}(Δ_e/c)[1](T_{mm′/c²}) = [1](T_mT_{m′}) for m, m′
admissible and prime to f_e. Hence B_e = (w_e/2)[1](SS*) with S := Σ_m b_e(m)T_m, S* := Σ_m b̄_e(m)T_m, and [1](SS*) = Σ_C|S_C|²
because 𝔞 ↦ 𝔞̄ preserves norms and coprimality to f_e and [𝔞̄] = [𝔞]^{−1}. □

**Corollary P (dual identity without Möbius terms).** For X ≥ 32 and b supported on admissible n,

  𝒬′(b) := Σ_{Δ∈𝒟} Φ(|Δ|/X) Q′_Δ(b) = ½ Σ_{n,n′} b(n) b̄(n′) Σ_{c|(n,n′)} 𝒩_c(nn′/c²; X),

with 𝒩_c(M; X′) := Σ Φ(|Δ|/X′)(Δ/c) over T, y ∈ Z, Δ ∈ 𝒟 with T² − Δy² = 4M, T ≡ yΔ (mod 2), as in [R3] Lemma 2; and

  Σ_{D∈𝔇(X)} Σ_{χ∈Ĝ(D)} |Σ_n b(n)λ_χ(n)|² ≤ C X^{1/2} log X · 𝒬′(b).

*Proof.* The identity is the definition of r_Δ summed over Δ (γ = (T + y√Δ)/2 lies in O_Δ iff T ≡ yΔ (mod 2)). For the inequality:
[R3] Lemma 2(i) gives Σ_χ|Σ b(n)λ_χ(n)|² = h(D)Σ_C|S_C|² = h(D)Q′_D(b) for D ∈ 𝔇(X); h(D) ≤ CX^{1/2}log X; Φ(|D|/X) = 1 on 𝔇(X); and
every other term Φ(|Δ|/X)Q′_Δ(b) is ≥ 0 by Lemma P. □

*Remarks.* (1) The right side is exactly the w = 1 part of [R3] Lemma 2(ii). Since Q′_Δ ≥ Σ_C|S_C|² term by term, 𝒬′ ≥ 𝒬 ≥ the fundamental
part; the Möbius terms of [R3] are not needed at all. (2) Regrouping Δ = e²Δ′: 𝒬′(b) = Σ_e Σ_{Δ′}Φ(e²|Δ′|/X)(w_{Δ′}/2)Σ_C|S_C^{(Δ′)}(b(e·))|²,
the [R3] sums at the scales X/e² (including |Δ′| ≤ 4 with weight w_{Δ′}/2); [R3]'s alternating w-terms telescope against these
(Σ_{ew=k}μ(w) = 1[k = 1]). Numerically 𝒬′ exceeds 𝒬 by about 1% (§8). (3) For general b (no pole killing) Corollary P with [R3] §3
gives only X(N + X^{1/2})‖b‖², the sum of the per-D sieves: the trivial-character term √X|Σ_n b(n)|² survives in the main term.

## 3. Theorem 1: the family large sieve at the family size, for all lengths below the conductor

**Theorem 1.** Let 0 < δ ≤ 1/100, X ≥ X₀(δ), ½ + 2δ ≤ σ ≤ 9/10, |t′| ≤ X^δ, 2 ≤ V ≤ X^{1/2}, and N = 2^j with V/4 ≤ N ≤ X^{1−4δ}.
Let b be b_{j,t′} or ∂_{t′}b_{j,t′} ([R3] §2). Then

  Σ_{D∈𝔇(X)} Σ_{χ∈Ĝ(D)} |Σ_n b(n)λ_χ(n)|² ≤ C(δ) X^{3/2+14δ} N^{1−2σ}.

([R3] Theorem 1 had (X^{3/2} + N³X^{−1})N^{1−2σ}.)

*Proof.* By Corollary P, 𝒬′(b) = ½Σ_c Q_c with Q_c := Σ_{n₂,n₂′} β_c(n₂) β̄_c(n₂′) 𝒩_c(n₂n₂′; X), the sum over odd squarefree c ≤ 4N,
β_c := b(c·), supported on admissible n₂ prime to c with n₂ ∈ [L/64, L], L := 4N/c. This is [R3]'s Q_{c,w} at w = 1, X′ = X, and it
suffices to show Σ_c|Q_c| ≤ C(δ)X^{1+13δ}‖b‖². [R3] §3.4–§3.8 bound Q_{c,1} with no reference to w other than through the following
hypotheses, which we check for every N ≤ X^{1−4δ}:
- (i) Lemma 6 (Poisson in T) needs (c, M) = 1, M ≤ X² and X′ ≥ 64c√M X^δ, for M = n₂n₂′. The first holds since n = cn₂ is
  squarefree; M ≤ 16N² ≤ X²; and c√M = √(nn′) ≤ 4N ≤ 4X^{1−4δ}, so 64c√M X^δ ≤ 256X^{1−3δ} ≤ X for X ≥ X₀(δ). This is the only place
  where [R3] needed w ≤ G₂ = X^{1−2δ}/(64N); here w = 1 ≤ G₂.
- (ii) Lemma 7 (local densities) needs M odd and cube-free: M = nn′/c² with n, n′ squarefree.
- (iii) §3.5: with v = X/(4M) ≥ c²X/(64N²), S(2^ate√v) = 0 unless 2^ate ≤ E₀ := 23N/(c√X).
- (iv) §3.6 (secondary term, Heath-Brown's quadratic large sieve over q ≤ E₀ with inner length L/κ) needs E₀ + L/κ ≤ 2L/κ, which
  follows from κ ≤ E₀ and E₀² ≤ L, i.e. N ≤ 4cX/529: true.
With these, [R3] gives: the y = 0 terms are ≤ 4XΣ|β_c|² (Lemma 5(a): y = 0 forces n₂ = n₂′); the secondary part is
≪ X^{7δ}c^{−1}X·Σ_{c|n}d(n)⁴|b(n)|² (§3.6); the main part is ≪ X^{11δ}c^{−2}√X·V·N^{1−2σ} (§3.7, Lemma 9: the trivial character's
pole is killed by (Z), the twisted sums are bounded by L(½+δ+iu, ψ(·)(·/q)) and Jutila's mean value over q ≤ E₀). Summing over c,
Σ_c|Q_c| ≪ X^{12δ}(X + √X V)‖b‖² ≪ X^{12δ}X‖b‖² (V ≤ X^{1/2}, Σ_n d(n)⁴|b(n)|² ≪ X^δ‖b‖²). Corollary P multiplies by CX^{1/2}log X. □

(The large-w analysis of [R3] §3.3, which produced N³X^{−1}, is no longer used anywhere.)

## 4. Theorem 2: the density hypothesis on average over D

**Theorem 2.** For 0 < δ ≤ 1/100, ½ + 3δ ≤ σ ≤ 9/10 and T := |t| + 3 ≤ X^δ/4,

  Σ_{D∈𝔇(X)} N^♮_D(σ; t) ≤ C(δ) X^{2−σ+60δ},

N^♮_D as in [R3] §4 (zeros with σ ≤ β ≤ 9/10, |γ − t| ≤ 1, of all χ ∈ Ĝ(D)∖{1}).

*Proof.* [R3] §4 verbatim with Y := X^{1−5δ} instead of X^{5/(7−2σ)−5δ}; V := X^{1/2}/2 and Z := 2Y log²Y, so 4Z ≤ X^{1−4δ} for
X ≥ X₀(δ) and Y ≥ V. Class II ([R3] Lemma 4(II), per-D mean values): Σ_D|𝒳_II(D)| ≤ C(δ)X^{3/2+2δ}Y^{1/2+δ−σ}T, and
3/2 + 2δ + (1−5δ)(½+δ−σ) + δ ≤ 2 − σ + 9δ. Class I: Theorem 1 (for |τ| ≤ X^δ/2; the per-D bound and the decay of Φ̂_β for larger |τ|,
exactly as in [R3]) gives Σ_D|𝒳_I(D)| ≤ C(δ)X^{16δ}Σ_{V/4≤2^j≤2Z}X^{3/2}2^{j(1−2σ)} ≤ C(δ)X^{3/2+17δ}V^{1−2σ} ≤ C(δ)X^{2−σ+18δ}.
Multiply by C₀log(XT²). □

For each D the class group family has h ≍ |D|^{1/2±ε} members of conductor |D|; its density hypothesis is Σ_χN(σ;t,χ) ≪ |D|^{1−σ+ε}
(the trivial count h at σ = ½, and 1 at σ = 1), and summed over D ∈ 𝔇(X) this is X^{2−σ+ε}. [R2] Proposition 1 gives each D only
|D|^{2(1−σ)/(3−2σ)+ε}. Together with [R] Proposition 1 (better for σ > (√41 − 3)/4 = 0.8508):
Σ_D N^c_D(σ,T) ≪ (XT)^ε X^{min{2−σ, 10(1−σ)/(3−2σ)}} for ½ + ε ≤ σ ≤ 1 (T ≤ X^ε).

## 5. Theorem 3: least primes in every class, D < 0, every A > 1

**Theorem 3.** [R] Theorem 2 holds for every A > 1: there are effective c, κ, δ > 0 (for 1 < A ≤ 1.2: δ = (A−1)/8000, κ = δ/4,
c = δ/2) and X₀ such that for X ≥ X₀ all D ∈ 𝔇 with |D| ≤ X outside a set of size ≪_A X^{1−c} satisfy, for every class C ∈ G, every
x ≥ |D|^A and every H ∈ [x|D|^{−κ}, x],

  #{𝔭 ∈ C prime of degree one : x < N𝔭 ≤ x+H} ≥ ν_D H/(16 h(D) log 2x)   ([R] (2.1); ν_D ≥ c|D|^{−δ/2}/log|D|).

In particular every primitive positive definite binary quadratic form of discriminant D represents a prime ≤ 2|D|^A. Corollary: for
all but O(N^{1−c}) squarefree n ≤ N with n ≡ 1, 2 (mod 4), every (x, x + x(4n)^{−κ}] with x ≥ (4n)^A contains a prime a² + nb², and the
least such prime is ≤ 2(4n)^A.

*Proof.* The statement for A implies it for every A′ > A (same exceptional set), so let 1 < A ≤ 1.2. Follow [R3] Theorem 3 (its four
changes to [R] §2) with steps 1 and 2 replaced:
1′. *Parameters.* ε′ := (A−1)/80, δ := ε′/100 = (A−1)/8000 (so δ ≤ A/4400 and δ ≤ (A−1)/80, the other two constraints of [R3] step 1),
   κ := δ/4, ε₀ := ε₁ := δ/100, σ₀ := ½ + 4δ, σ_j := σ₀ + jε′/A (0 ≤ j ≤ J, J minimal with σ_J ≥ 1 − ε′/A); σ_g, T_X, Q₀ as in [R] §2.1.
2′. *E₁ (Markov).* For σ ∈ [σ₀, 1): Σ_{D∈𝔇(X)}N^c_D(σ, T_X) ≤ C(δ)X^{F*(σ)+61δ}, where F*(σ) := 2 − σ for σ ≤ 9/10 (Theorem 2 summed
   over the ≤ 2T_X unit windows for the zeros with β ≤ 9/10; the zeros with β > 9/10 by [R] Proposition 1, F₂(9/10) = 5/6 < 2 − σ),
   and F*(σ) := F₂(σ) = 10(1−σ)/(3−2σ) for σ > 9/10. Then φ(σ) := F*(σ) − A(1−σ) satisfies φ(σ) = 1 − (A−1)(1−σ) ≤ 1 − (A−1)/10
   = 1 − 8ε′ on [σ₀, 9/10], and φ(σ) ≤ F₂(σ) ≤ 5/6 ≤ 1 − 8ε′ on (9/10, 1). Markov as in [R] §2.2:
   #E₁(X) ≤ Σ_j X^{2ε′−A(1−σ_j)}·C(δ)X^{F*(σ_j)+61δ+2(κ+ε₁)} ≤ C(δ)J X^{1−6ε′+62δ} ≤ C(δ)J X^{1−ε′}.
3. *Zeros with β < σ₀* ([R3] step 3, unchanged): E_b ≤ H|D|^{1/2−A(1/2−4δ)+κ+ε₁+o(1)} ≤ H|D|^{−(A−1)/2+6δ}, and (A−1)/2 − 6δ > δ/2.
4. *Inherited errors* ([R3] step 4, unchanged): H|D|^{3δ−A/40+o(1)} and H|D|^{1/2−A+o(1)}, both ≤ H|D|^{−δ/2−η}.
[R] §2.5 then needs ε′ > δ/2 and 0.005A − κ − ε₁ > δ/2 (true) and x ≥ |D|^{1+2δ} for prime powers, inert and ramified primes (true,
A > 1 + 2δ). E₂ and §2.3, §2.4(a), (c), (d), (e) of [R] need only A ≥ 1, δ ≤ A/4400 and the margins above. c := δ/2 ≤ min{ε′, 0.505δ}. The corollary is
[R2] Corollary 4′ (D = −4n, principal class). □

*Where A > 1 is used.* Only in step 3 (the zeros with β < ½ + 4δ, counted trivially: h zeros per window, each weighted by
x^{β−1} ≤ x^{−1/2+4δ}) and in Markov near σ₀ (φ(σ₀) = 1 − (A−1)(½ − 4δ)). Both are exactly the explicit-formula limit that GRH gives.

## 6. D > 0 (narrow classes, indefinite forms) and sectors

**Lemma P⁺.** Let Δ > 0 be a non-square discriminant with conductor f, w = ν∗ν as in [R4] §1 (even, ŵ ≥ 0, supp w ⊂ [−1, 1]),
u_γ := ½log(γ/γ′), and for b supported on admissible n put 𝒲′_Δ(M) := Σ_{γ∈O_Δ, γ≫0, Nγ=M} w(u_γ) and
Q′_Δ(b) := Σ_{n,n′}b(n)b̄(n′)Σ_{c|(n,n′)}(Δ/c)𝒲′_Δ(nn′/c²). Then Q′_Δ(b) = Σ_{e|f squarefree} Q_{Δ_e}(b_e) ≥ 0, where Q_{Δ′}(β) is
[R4]'s Q (with the factor 1[(M, f′) = 1], automatically 1 here) for the order of discriminant Δ_e = Δ/e²; and Q′_Δ ≥ Q_Δ.

*Proof.* Steps (a), (b) of Lemma P verbatim (T² = 4M + Δy² ≡ 4M mod p²; γ ↦ γ/e); γ/e is totally positive iff γ is, and
u_{γ/e} = u_γ. Each block Q_{Δ_e}(b_e) is ≥ 0 by [R4] Lemma 2⁺(i) for the order of discriminant Δ_e (non-square because Δ is):
R⁺_{Δ_e}Q_{Δ_e}(b_e) = Σ_k ŵ(k/R⁺_{Δ_e})Σ_C|S_C^{(k)}|². The block e = 1 is Q_Δ(b). □

**Theorem 1⁺.** [R4] Theorem 1⁺ holds without the N³X^{−1} term: for the parameters and sequences of Theorem 1,
Σ_{D∈𝔇⁺(X)}Σ_{χ∈Ĝ⁺(D)}|Σ_n b(n)λ_χ(n)|² ≤ C(δ)X^{3/2+O(δ)}N^{1−2σ}, V/4 ≤ N ≤ X^{1−4δ}.

*Proof.* By [R4] Lemma 2⁺(i), (ii) and Lemma P⁺: LHS ≤ ŵ(0)^{−1}Σ_{D∈𝔇⁺(X)}h⁺R⁺Q_D(b) ≤ CX^{1/2}log²X·Σ_ΔΦ(Δ/X)Q′_Δ(b), the sum
over all non-square Δ > 0 (Q′_D = Q_D for fundamental D; all other terms ≥ 0). As in [R4] Lemma 2⁺(iii) without the w-inversion,
Σ_ΔΦ(Δ/X)Q′_Δ(b) = Σ_c Q⁺_c with Q⁺_c := Σ_{n₂,n₂′}β_c(n₂)β̄_c(n₂′)𝒩⁺_c(n₂n₂′; X). Then [R4]'s proof of Theorem 1⁺ at w = 1: Lemma 5⁺
(unchanged); square Δ (included by Poisson, excluded here) contribute ≪ X^δN‖b‖² ≤ X^{1+δ}‖b‖²; Lemma 6⁺ needs X ≥ 128c√M X^δ, true as
c√M ≤ 4N ≤ 4X^{1−4δ}; §3.5 with g⁺; §3.6–§3.8 of [R3] as in Theorem 1 (E₀² ≤ L as N ≤ X^{1−4δ}). □

**Lemma P^∠ (sectors, D < 0, ℓ even).** For Δ ∈ 𝒟 put Q′^∠_Δ(b) := Σ_{n,n′}b(n)b̄(n′)Σ_{c|(n,n′)}(Δ/c)·½Σ_{γ∈O_Δ, Nγ=nn′/c²}(γ/|γ|)^ℓ.
Then Q′^∠_Δ(b) = Σ_{e|f squarefree} B^∠_e ≥ 0, where B^∠_e = Σ_{C∈G_{Δ_e}}|S_C^∠(b_e)|² ([R4] Lemma 2^∠(i), with ψ_ℓ for the order
O_{Δ_e}) if |Δ_e| > 4, and for Δ_e ∈ {−3, −4} (h = 1, w_e = 6, 4): B^∠_e = (w_e/2)|Σ_m b_e(m)λ_{ψ_ℓ}(m)|² if w_e | ℓ and B^∠_e = 0 otherwise.

*Proof.* Steps (a), (b) of Lemma P; (γ/e)/|γ/e| = γ/|γ|. For |Δ_e| > 4 the block is [1](S^∠S^∠*) = Σ_C|S_C^∠|², by [R4] Lemma 2^∠ for
the order O_{Δ_e} (all m prime to f_e). For Δ_e ∈ {−3, −4} every ideal is principal; the w_e unit multiples uγ′ contribute
Σ_u u^ℓ(γ′/|γ′|)^ℓ = w_e1[w_e | ℓ](γ′/|γ′|)^ℓ; if w_e | ℓ then ψ_ℓ((γ′)) := (γ′/|γ′|)^ℓ is a Hecke character with ψ_ℓ((p)) = 1, so
[R3] Lemma 1 holds for λ_{ψ_ℓ} and the block is (w_e/2)|Σ_m b_e(m)λ_{ψ_ℓ}(m)|². □

**Theorem 1^∠.** [R4] Theorem 1^∠ holds without the N³X^{−1} term, uniformly for even |ℓ| ≤ X^δ, N ≤ X^{1−4δ} (constant X^{O(δ)}).
*Proof.* Lemma P^∠ in [R4]'s positivity step; the dual is the w = 1 part of [R4] Lemma 2^∠(ii); the Poisson hypothesis
X ≥ 64c(1 + 2ℓ²)√M X^δ holds as c√M ≤ 4N and ℓ² ≤ X^{2δ}; the rest is [R4] Theorem 1^∠ at w = 1. □

**Theorems 2⁺, 2^∠.** Σ_{D∈𝔇⁺(X)}N^♮_D(σ;t) and Σ_{D∈𝔇(X)}Σ_{ψ∈ψ_ℓĜ(D)}N^♮(σ;t,ψ) (even 0 < |ℓ| ≤ X^δ) are ≤ C(δ)X^{2−σ+80δ} for
½ + 3δ ≤ σ ≤ 9/10. *Proof.* [R4] Theorems 2⁺, 2^∠ with Theorems 1⁺, 1^∠ and Y := X^{1−5δ}, as in Theorem 2. □

**Theorem 3⁺ (D > 0).** [R2] Theorem 4 and Corollary 5 hold for every A > 1 (same c, κ, δ as in Theorem 3): for all but O(X^{1−c})
fundamental D ∈ (0, X], every narrow class has ≥ ν_D H/(16h⁺ log 2x) prime ideals of degree one with norm in (x, x+H] for all x ≥ D^A,
H ∈ [xD^{−κ}, x]; every primitive indefinite binary quadratic form of discriminant D represents a prime ≤ 2D^A; and for all but
O(N^{1−c}) squarefree n ≤ N with n ≡ 2, 3 (mod 4) the least prime a² − nb² is ≤ 2(4n)^A.
*Proof.* [R4] Theorem 3⁺ (the dictionary of [R2] Theorem 4) with steps 1′, 2′ of Theorem 3: F*(σ) = 2 − σ from Theorem 2⁺ for
σ ≤ 9/10 and F₂ from [R2] Proposition 2 above; per window h⁺ ≤ 5D^{1/2}log²(3D) zeros in step 3. □

**Theorem 3^∠ (sectors).** [R] Theorem 3, [R2] Theorem 3′ and Corollary 4′ hold for every A > 1: for all but O(X^{1−c}) D ∈ 𝔇 with
|D| ≤ X, every class C, every arc I ⊂ R/(2π/w)Z of length θ ≥ |D|^{−κ}, every x ≥ |D|^A and H ∈ [x|D|^{−κ}, x], there are
≥ wθν_D H/(128πh log 2x) degree-one 𝔭 ∈ C with N𝔭 ∈ (x, x+H] and arg β_𝔭 ∈ I; and for all but O(N^{1−c}) squarefree n ≤ N,
n ≡ 1, 2 (4), every (x, x + x(4n)^{−κ}], x ≥ (4n)^A, contains a prime a² + nb² with arg(a + ib√n) in any prescribed arc of length ≥ (4n)^{−κ}.
*Proof.* [R4] Theorem 3^∠ with steps 1′, 2′: in [R] Theorem 3(iv) the zeros with β ≤ 9/10 of ψ ∈ 𝒜(D, M)∖Ĝ, M = 6|D|^{κ+ε₂}, are counted
by Theorem 2^∠ summed over ≤ 2M_X + 1 values of ℓ (M_X ≤ X^δ), those with β > 9/10 by [R] Proposition 1′ (exponent F₂ + O(κ + ε₂)), so
the family exponent is F*(σ) + O(δ); the zeros with β < σ₀ give ≤ |D|^{1/2−A(1/2−4δ)+2κ+ε₁+ε₂+o(1)} ≤ |D|^{−(A−1)/2+6δ} relative to
ĝ(0)H/h. □

## 7. Task 2: the other routes, and the barrier at A = 1

| Route (brief) | What it would change | Verdict |
|---|---|---|
| Hecke relations / dual lattice identity pushed further | Lemma P: positivity of the all-elements count removes the conductor Möbius terms, hence N³X^{−1}; A₀ = 1 | done (§2–§6) |
| Duke–Friedlander–Iwaniec bilinear forms / amplification | pointwise subconvexity for L(½, χ); class II is a mean square that is already sharp (≍ h per D), and class I is now at the family size | not needed; cannot go below 1 |
| Heath-Brown quadratic large sieve refinements | it enters only the secondary term (§3 (iv)), already ≪ X^{1+ε}‖b‖², i.e. the family size X^{3/2} after the factor h ≍ X^{1/2} | not needed |
| fourth-moment / mollified high-moment family bound | would push the density below X^{2−σ} (toward X^{3(1−σ)}); but at σ → ½ even X^{3(1−σ)} meets the bar 1 + A(1−σ) only for A ≥ 1 | not needed; cannot go below 1 |
| per-D log-free density, Huxley large values, hybrid t-aspect | σ_c no longer exists; per-D input is used only for β < ½ + 4δ, where the trivial count suffices for A > 1 | not needed |
| family sieve beyond the conductor (N > X) | Poisson in T fails already for w = 1, c = 1 when √(nn′) ≫ X; not needed for A; it is what ChebST's 11/3 and the a_E-sign thresholds ([R4] §3) need | open |
| no Landau–Siegel zeros (Kerr 2610.21099, audit [Siegel](../Siegel/REPORT.md)) | (1−β)log q ≥ c_K for all real zeros gives ν_D ≥ min{1, c_K} in [R] (2.1): every class gets ≥ min{1,c_K}H/(16h log 2x) primes, uniformly; the exponent and E₂ ([R] (I4), zeros near 1 of all χ_d) are unchanged | COROLLARY (rests on an unrefereed 2026 proof) |
| quasi-RH half-plane Re s > 0.8749 for all Hecke L (Liu 2610.26835, [LiuHecke](../LiuHecke/REPORT.md)) | removes the zeros with β > 0.875 ([R] Proposition 1, the 9/10 cutoff); the binding point is now σ → ½, so A stays 1 | no effect |

*The barrier at 1.* For A < 1 two things fail. (i) Step 3 of Theorem 3: each unit window holds ≍ h ≍ |D|^{1/2} zeros with β near ½,
and Σ x^{β−1} over them is ≍ h·x^{−1/2} unless the terms χ̄(C)x^{ρ_χ} cancel; h·x^{−1/2} ≤ |D|^{−η} needs x ≥ |D|^{1+2η}. (ii) Markov: at
σ = ½ + ε the bar is X^{1+A/2−Aε}, while even the family density hypothesis X^{3(1−σ)} is X^{3/2−3ε}; so A < 1 fails as ε → 0. Under
GRH (i) is the whole error, and the explicit formula with absolute values gives exactly A = 1 for every D. Going below 1 needs
cancellation among the zeros themselves (pair-correlation-type information over χ ∈ Ĝ(D)). NOGO for this method.

## 8. Numerical checks ([checks/](checks/), every script under a hard alarm, all processes < 10 MB)

- `lemmaP_check.py` (Lemma P, P⁺, P^∠). b random complex on the 27 odd squarefree n ≤ 255 built from 3, 5, 7, 11, 13, 17 (small primes,
  so that they divide many conductors; the algebra needs only n odd). For each Δ: Q′_Δ(b) from a direct count of (T, y), against the sum of
  the blocks (each block computed the same way for Δ/e² with b_e).
  - D < 0, all Δ ∈ [−1600, −7] with f > 1 and 15% of the fundamental ones: 388 tests, 312 non-maximal, 144 with a nonzero block e > 1. Max
    relative deviation 1.8·10^{−15}. In addition every block was recomputed independently as (w_e/2)Σ_C|S_C|² from the reduced forms of
    Δ/e² (S_C from representation numbers of the forms): agreement to 10^{−15}. No negative block.
  - Sectors, ℓ = 2, 4, 6: 1185 tests (936 non-maximal, 408 with a nonzero block e > 1), deviation 6·10^{−14}, smallest block
    −6·10^{−14}|Q′| (rounding).
  - D > 0, non-square Δ ≤ 1200, Gaussian w (width 0.35, ŵ > 0): 257 tests (201 non-maximal, 84 with a nonzero block e > 1), deviation
    8·10^{−16}, no negative block.
- `poisson_check.py` ([R3] Lemma 6 for M with p² ∥ M, the pairs kept by Lemma P; C^∞ Φ): brute-force 𝒩_c(M; X) against the main term,
  e.g. M = 19²·23·29, X = 4·10⁵, c = 1: 1806.5046 vs 1806.5046; M = 19²·23²·29, X = 2·10⁶, c = 3: −1446.2725 vs −1446.2725;
  M = 23²·29·31·37, X = 3·10⁶, c = 7 (X/(c√M) = 102): −243.1299 vs −243.1227.
- `famq.c` (the family sums by brute force over all orders): Φ a sin²-bump on [¼, 2] (= 1 on [½, 1]), n odd squarefree prime to 3,
  σ = 5/6, t′ = 0, V = √X/2; b = a_V(n)n^{−σ}(φ₂(n/N) − 4^{1−σ}φ₂(4n/N)) (the pole-killed detector weight, (Z) holds), φ₂ a smooth bump on
  (¼, 4). 𝒬′ = Σ_ΔΦQ′_Δ (Lemma P), 𝒬 = [R3]'s sum (ideals prime to f), and the fundamental family sum
  F := Σ_{D fund}Φ(|D|/X)h(D)Σ_C|S_C|² (the left side of Theorem 1, weighted). "rand": random signs times n^{−σ}φ₂(n/N);
  "pole": n^{−σ}φ₂(n/N) (trivial-character main term present). ‖b‖² := Σ|b(n)|².

| X | N | 𝒬′/(X‖b‖²) | 𝒬/(X‖b‖²) | F/(X^{3/2}‖b‖²) | rand 𝒬′/(X‖b‖²) | pole 𝒬′/(X‖b‖²) |
|---|---|---|---|---|---|---|
| 20000 | 1250 | 0.534 | 0.528 | 0.158 | 0.579 | 9.85 |
| 20000 | 2500 | 0.530 | 0.525 | 0.154 | 0.589 | 19.2 |
| 20000 | 5000 | 0.515 | 0.509 | 0.150 | 0.613 | 37.7 |
| 20000 | 10000 | 0.500 | 0.493 | 0.146 | 0.606 | 74.6 |
| 20000 | 20000 | 0.507 | 0.500 | 0.149 | 0.632 | 148.5 |
| 20000 | 30000 | 0.502 | 0.496 | 0.148 | 0.642 | 222.6 |
| 40000 | 5000 | 0.534 | 0.530 | 0.154 | – | – |
| 40000 | 20000 | 0.504 | 0.498 | 0.147 | – | – |
| 40000 | 40000 | 0.514 | 0.508 | 0.148 | – | 210.1 |
| 40000 | 60000 | 0.512 | 0.507 | 0.148 | – | – |

  The pole-killed sums stay at the family size from N = X/16 to N = 1.5X (Theorem 1 proves N ≤ X^{1−4δ}); 𝒬′ exceeds 𝒬 by 0.8–1.6%;
  random signs behave alike; the pole-carrying sequence grows like N/√X (148.5 at X = N = 20000, 210.1 at X = N = 40000: ratio √2), the
  term √X|Σb|² of Remark (3). Runtime: X = 40000, N = 60000 in about one minute, 2 MB.

## 9. Novelty (2026-10-21)

- arXiv monthly listings, math.NT, 2026-09 (812 entries) and 2026-10 (590 entries, through 2610.26835): all titles scanned for quadratic
  form, class group/number, ring class, ideal class, class field, least prime, Linnik, Chebotarev, zero density, zero-free, density
  hypothesis, large sieve, Hecke, discriminant, imaginary quadratic, dihedral, theta, Siegel, represent, almost all. Abstracts read:
  Cossaboom 2609.27715 (zero density in the T-aspect for one GL₂ form), Frei–Sofos–Widmer 2610.18491 (effective Chebotarev under
  quasi-RH for Dedekind zeta; least prime splitting; ℓ-torsion), de Faveri 2610.04045 (large sieve for fixed-order characters),
  Liu 2610.12234 (half-plane 0.87496), Trudgian–Zuniga-Alterman 2610.16010, Felber 2609.07133. None treats class-group families on
  average over D or least primes in ideal classes for almost all D.
- arXiv search (HTML interface; the API answered "Rate exceeded"): 16 queries (least prime × forms × almost all discriminants; large sieve ×
  class group characters; zero density × class group L-functions; primes × ideal class × almost all imaginary quadratic; ring class ×
  zero density; primes represented × average over discriminants; density hypothesis × average × discriminant; indefinite forms × least
  prime; Linnik problem × forms × average; least prime ideal class × average; large sieve × theta series weight one; Bombieri–Vinogradov ×
  binary forms; primes represented by every form; dihedral large sieve; class group L-functions × moments × average; short intervals ×
  forms × almost all). Relevant hits only: Ditchen 1312.1502, Zaman 1710.08914, Gaudet 2408.01610, Chen An 2004.14410,
  Pierce–Turnage-Butterbaugh–Wood 1709.09637, Lam 1812.08823 (Bombieri–Vinogradov for one fixed form); all known from the earlier audits.
- GXL evidence API: three claims sent, all refused with HTTP 429 (daily quota used up).
- Earlier audits ([R] NOVELTY, [R2] NOVELTY, [R3] NOVELTY, adversarial; [R4] §6): the literature's best almost-all exponent for every form
  is Ditchen's 20/3 (definite forms, density-zero exceptions); 1.6536, 3/2 and 1.1307 are this project's.
- Nothing found states or implies any almost-all exponent below 1.1307, nor the family sieve of Theorem 1, nor Theorem 2.
- Method: the pair count behind Corollary P (two ideals in one class ⟺ 4nn′ = T² − Δy²) is Soundararajan's device, summed over d by
  Heath-Brown–Pierce; new here is that the all-elements count of a non-maximal order splits into nonnegative class-group sums of the
  orders Δ/e² (Lemma P), which removes the conductor inclusion–exclusion and with it the only loss in [R3].
- Caveat: for a² + nb² alone, Baier–Zhao (Math. Ann. 338, 2007) give a prime m² + n < 4n for all but O(N(log N)^{−C}) squarefree n ≤ N,
  a better exponent for the least prime of the principal form; the corollary here adds power-saving exceptional sets, short intervals,
  arcs, and every class of discriminant −4n.

## 10. Corrections, dead ends (NOGO candidates), pointers

- *Correction to [R3] §6 and [R4] §4.* The N³X^{−1} term is not intrinsic and no Salié-sum cancellation is needed: it comes from summing
  non-maximal orders over ideals prime to the conductor. Any NOGO entry based on that diagnosis ("family sieve near N ≍ X needs
  half-integral weight Kuznetsov") should be replaced by this note.
- *A < 1 by the explicit formula with absolute values*: impossible (§7, barrier). NOGO.
- *Liu's half-plane for the almost-all exponent*: no effect once A = 1 (binding point σ → ½). NOGO.
- *Family sieve beyond the conductor* (N > X): still open; Poisson in T fails for w = 1 itself. Needed by ChebST (lengths X^{3/2}–X^{7/2})
  and by the a_E-sign theorem, not by A.
- *Pointer: CM2.* [CM2 NOTE](../CM2/NOTE.md) Theorem 1 (CM fields over a fixed totally real F) has the same 𝔴-inversion (its Lemma 4(ii),
  §4.3) and states that removing N³X^{−1} would give A₀ = 1. The proof of Lemma P (vanishing at 𝔭 ∥ Nα for 𝔭 | 𝔣, division by 𝔴 as in
  its Lemma 2(c)) looks transferable; not checked here (outside this brief).

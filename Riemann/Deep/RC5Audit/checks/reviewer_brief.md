You are refereeing one decisive step of a paper proof (analytic number theory). No tools needed: answer from the text below only, within 90 seconds of thinking. Be adversarial: look for a false inequality, a hidden hypothesis, or a step that silently needs coprimality to the conductor or w > 1.

NOTATION ([R3] §0): 𝒟 = negative discriminants Δ (≡ 0,1 mod 4, any order), Δ = f²D_K, O_Δ = Z + fO_K, G_Δ = Cl(O_Δ) ≅ I_K(f)/P_{K,Z}(f); T_M := Σ_{N𝔞=M, (𝔞,f)=1}[𝔞] ∈ Z[G_Δ] (O_K-ideals), λ_χ(n) = χ(T_n); [1](·) = coefficient of identity; "admissible" n = odd squarefree with all prime factors ≥ 19; Φ smooth, supp ⊂ [1/4,2], 0 ≤ Φ ≤ 1, Φ = 1 on [1/2,1]; 𝔇(X) = fundamental D<0 with X/2<|D|≤X.
[R3] Lemma 1: for |Δ|>4 and admissible n,n', g=(n,n'), n=gm, n'=gm': T_nT_{n'} = Σ_{c|g}(Δ/c)T_{mm'(g/c)²} (proof: T_p² = T_{p²} + (Δ/p) for p∤f, using [(p)] = 1; T_p = T_{p²} = 0 and (Δ/p)=0 for p | f).

=== RCApps5 NOTE lines 50-122 (verbatim) ===
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


=== [R3] Lemma 1-2 region (verbatim excerpt) ===
- P := Π_{p≤17}p. "n admissible" means n odd squarefree with (n, P) = 1 (all prime factors ≥ 19).
- δ ∈ (0, 1/100] is a small parameter; C, C(δ) are absolute resp. δ-dependent effective constants; ε-losses are X^{O(δ)}.

## 1. Dual identities

**Lemma 1 (Hecke relations in the group ring).** Let Δ ∈ 𝒟, |Δ| > 4, and n, n' admissible, g := (n,n'), n = gm, n' = gm'.
Then T_nT_{n'} = Σ_{c|g}(Δ/c)T_{mm'(g/c)²} in Z[G_Δ].

*Proof.* T_{ab} = T_aT_b for coprime a, b (unique factorisation of ideals prime to f). For a prime p ∤ f: if p splits,
p = 𝔭𝔭̄, then T_p² = [𝔭]² + [𝔭̄]² + 2[(p)] and T_{p²} = [𝔭²] + [𝔭̄²] + [(p)], [(p)] = 1, so T_p² = T_{p²} + 1; if p is inert,
T_p = 0 and T_{p²} = [(p)] = 1, so T_p² = T_{p²} − 1; if p ramifies, T_p = [𝔭], T_p² = [(p)] = 1 = T_{p²}. In all cases
T_p² = T_{p²} + (Δ/p). If p | f then T_p = T_{p²} = 0 and (Δ/p) = 0. Hence T_g² = Π_{p|g}(T_{p²} + (Δ/p)) = Σ_{c|g}(Δ/c)T_{(g/c)²},
and T_nT_{n'} = T_{mm'}T_g². □



=== [R3] Lemma 6 (verbatim) ===

**Lemma 6.** Let c be odd squarefree, (c, M) = 1, X' ≥ 64c√M X^{δ}, M ≤ X². Then the y ≠ 0 part of 𝒩_c(M; X') equals

  (μ(c)/c)(X'/√M)Σ_{y≥1, (y,c)=1}A_M(y)K(y²X'/(4M)) + O_δ(X^{−10}),
  A_M(y) := ¼#{T₀ mod 4y² : y² | T₀² − 4M, (T₀² − 4M)/y² ≡ 0, 1 (mod 4)}.

*Proof.* For fixed y ≥ 1 (and −y), Δ = (T² − 4M)/y² ∈ 𝒟 iff T lies in the classes counted by A_M(y) mod 4y²; then
T ≡ yΔ (2) is automatic (T² ≡ Δy² (4)), and (Δ/c) depends on T mod 4y²c. Poisson summation over each class mod q := 4y²c of
f(t) := Φ((4M − t²)/(y²X')) gives main terms q^{−1}∫f and errors bounded by Σ_{k≠0}|f̂(k/q)|. Put u(t) := (4M − t²)/(y²X');
on supp f, u ∈ [1/4, 2], so y²X' ≤ 16M, |u'| ≤ 4√M/(y²X'), |u''| = 2/(y²X') ≤ (4√M/(y²X'))², u''' = 0; hence


=== [R3] §3.6 (verbatim) ===

  Σ_{n₂,n₂'}β(n₂)β̄(n₂')𝒩_c^{y≠0}(n₂n₂';X') = (μ(c)/c)[2c₁√X'·𝔐_{c,w} + X'(2πi)^{−1}∫_{(−δ)}Š(z)(4/X')^{z/2}𝔖_{c,w}(z)dz] + O(X^{−8}),
  𝔐_{c,w} := Σ_{(a,t,e)∈𝔈}μ(t)(2^ate)^{−1}Σ_{n₂,n₂'}β(n₂)β̄(n₂')A^{(2)}_{n₂n₂'}(2^a)h_{n₂n₂'}(e),
  𝔖_{c,w}(z) := Σ_{(a,t,e)∈𝔈}μ(t)(2^ate)^{−z}Σ_{n₂,n₂'}β(n₂)β̄(n₂')(n₂n₂')^{(z−1)/2}A^{(2)}_{n₂n₂'}(2^a)h_{n₂n₂'}(e).

By (S), each inner double sum equals Σ_ψ â_a(ψ)Σ ϑ·B_{κ,q}(z)·B'_{κ',q}(z), summed over the triples (κ, κ', q) with
T₁⋯T₅² = e, where (L := 4N/(cw) bounds the support of β)

  B_{κ,q}(z) := Σ_{n: κ|n} β(n)n^{(z−1)/2}ψ(n)F_{κ,T₁,T₅}(n),  B'_{κ',q}(z) := Σ_{n': κ'|n'} β̄(n')n'^{(z−1)/2}ψ(n')F_{κ',T₁,T₅}(n').

Put A_κ := Σ_{κ|n}d(n)|β(n)|².

### 3.6 The secondary term (Re z = −δ)
On Re z = −δ, |(2^ate)^{−z}| ≤ E₀^{δ}, |(4/X')^{z/2}| ≤ X'^{δ/2}, and |n^{(z−1)/2}| ≤ 64L^{−(1+δ)/2} on supp β ⊂ [L/64, L]. Fix (a, t, ψ, T₂, T₃, T₄, T₅) and sum over T₁ (it only enters through q): by Cauchy–Schwarz and Heath-Brown's
quadratic large sieve [HB95, Thm 1],

  Σ_{q≤Q odd sqfree}|Σ_{m≤L', m odd sqfree}γ_m(m/q)|² ≤ C(ε)(QL')^ε(Q + L')Σ_m|γ_m|²,

applied with Q = E₀, L' = L/κ (the coefficients γ_m of m = n/κ are fixed: β(κm)(κm)^{(z−1)/2}ψ(κm)·sign), we get
Σ_{T₁}|B_{κ,q}B'_{κ',q}| ≪ X^{δ}(E₀ + L/κ)^{1/2}(E₀ + L/κ')^{1/2}L^{−1−δ}A_κ^{1/2}A_{κ'}^{1/2}. Since κ ≤ e ≤ E₀ and
E₀²/L ≤ 23²Nw/(cX) ≤ 1 (w ≤ G₂), E₀ + L/κ ≤ 2L/κ, so this is ≪ X^{δ}L^{−δ}(κκ')^{−1/2}A_κ^{1/2}A_{κ'}^{1/2}. With
|ϑ| ≤ 2^{ω(T₄)}T₄T₅ and T₄T₅(κκ')^{−1/2} = (T₂T₃)^{−1/2}, and dropping the constraint (a,t,e) ∈ 𝔈 in the remaining sums of


QUESTIONS (answer each: CORRECT or GAP, and SIGNIFICANT or ROUTINE, with a one-line reason; quote the exact phrase if you claim a gap):
Q1. Lemma P (identity and positivity), steps (a), (b), (c).
Q2. Corollary P: the identity and the inequality Σ_D Σ_χ |Σ b λ_χ|² ≤ C X^{1/2} log X · 𝒬′(b).
Q3. Theorem 1's reduction: 𝒬′(b) = ½Σ_c Q_c with Q_c = [R3]'s Q_{c,1} (w = 1, X' = X), and [R3] §3.4–§3.8 applied at w = 1 under (i)–(iv) for all N ≤ X^{1−4δ}. Is any hypothesis of Lemma 6 or §3.6 violated, or any use of w > 1 / coprimality to the conductor hidden?
Q4. Overall: does Theorem 1 (for the detector sequences only) follow? Note: for arbitrary b it is false (trivial and genus characters give ≍ X·N·‖b‖²); the theorem claims only the pole-killed detector sequences b_{j,t'}, ∂_{t'}b_{j,t'}.
Context (my own numerics, X = 2000-32000, N up to X): detector sequences give LHS/(Σ_D h(D)·Σ|b|²) = 0.8-1.1, flat in N and X; without pole killing it grows like N/√X; Lemma P verified exactly on 231 non-maximal orders with my own code.
Final line: OVERALL: CORRECT or GAP, SIGNIFICANT or ROUTINE.
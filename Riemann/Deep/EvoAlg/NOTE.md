# When is A(G) ≅ A_RW(G)? A complete answer in characteristic 0, and what changes in characteristic p

Orca deep dive, 2026-10-22 (brief: .scratch/orca/evoalg-brief.md). Paper proofs, not refereed, no Lean.
Computer checks: [verify.py](verify.py) (default run about 25 s; --full 6–8 min). Summary: [REPORT.md](REPORT.md).

## 0. Problem and results

**Problem.** arXiv 2610.22248 (Kaygorodov, Martín González, Viruel, *Evolution algebras: a survey of results and a
list of open problems*), Question "Graphs versus random walks", TeX l. 1351–1354: *Prove that for every finite
connected graph G one has A(G) ≅ A_RW(G) if and only if G is regular or biregular.* The survey fixes an arbitrary
field 𝕜 (l. 107). This is the conjecture of Cadavid, Rodiño Montoya and Rodriguez (CRR: arXiv 1710.10516 = Linear
Multilinear Algebra 69 (2021) 1858–1877, Conjecture at TeX l. 963–964). CRR work over R (l. 139: "In what follows,
we always assume that A is real") and proved the case of a non-singular adjacency matrix (their Theorem 2.3,
l. 525–526). Supporting singular families: Cadavid, Rodiño Montoya, Rodriguez, Vidal (CRRV: arXiv 2405.12341).

A(G) has natural basis (e_v)_{v∈V} with e_v² = Σ_{w∼v} e_w; A_RW(G) has natural basis (e'_v) with
e'_v² = d_v⁻¹ Σ_{w∼v} e'_w, where d_v = deg v must be invertible in 𝕜. *Biregular* means bipartite with constant
degree d_i on side V_i (CRR, l. 159). Throughout, G is finite and connected with n ≥ 2 vertices (for n = 1, A_RW is
undefined).

**Theorem A (characteristic 0).** Let char 𝕜 = 0. Then A(G) ≅ A_RW(G) over 𝕜 if and only if
1. G is regular, or
2. G is bipartite with every vertex of side V_i of degree d_i (i = 1, 2), and d₁/d₂ ∈ (𝕜^×)³
   (equivalently d₁²d₂ ∈ (𝕜^×)³, since d₁²d₂ = d₁³·(d₂/d₁)).

**Corollary A.**
- Over R, over C, and over every field of characteristic 0 that contains a cube root of every prime, the conjecture
  is true for all finite connected graphs, singular or not.
- Over a field 𝕜 of characteristic 0 the conjecture holds for all G if and only if every positive rational number is
  a cube in 𝕜 (K_{m,n} is biregular with d₁/d₂ = n/m).
- Over Q it is false. The regular or biregular graphs with A(G) ≇ A_RW(G) over Q are exactly the biregular G whose
  degree ratio is not a rational cube: P₃ = K_{1,2} (ratio 2), K_{1,3}, K_{2,3}, the subdivided K₄, …; whereas
  K_{1,8} and K_{2,16} (ratio 8) are isomorphic over Q.

**Theorem P (characteristic p > 0; §6).** In characteristic 2, A_RW(G) is defined exactly when all degrees are odd,
and then A_RW(G) = A(G). In odd characteristic p the answer depends on the degrees mod p and, if some twin class has
size divisible by p, on similitude multipliers of the form x₁² + … + x_k². The equivalence "regular or biregular"
fails in every positive characteristic (double star S(p,p)). Twin classes of size divisible by p produce isomorphisms
that mix twins and exist for no diagonal reason: K_{1,1,3} and a 6-vertex tree are isomorphic to their random-walk
algebras over F₉ but not over F₃.

How the proof goes: natural bases of non-degenerate evolution algebras are known up to orthogonal changes inside
classes of proportional squares (§2; Boudi–Cabrera Casado–Siles Molina). For A(G) these classes are the twin classes,
so an isomorphism is a twin-class-block matrix (§3). Inside a twin class the all-ones vector forces the block's
scale to agree with the scale of its neighbours (norm identity, §4). That kills the orthogonal freedom in
characteristic 0 and leaves the scaling system of the twin-free case (§5).

## 1. Conventions

- N(v) is the neighbourhood of v. For U ⊆ V put 1_U = Σ_{u∈U} e_u ∈ A(G) and 1'_U = Σ_{u∈U} e'_u ∈ A_RW(G). The
  vector (1, …, 1)ᵀ is written 𝟙.
- *Twins:* v ≡ w iff N(v) = N(w). Twins are non-adjacent, since v ∉ N(v). For a twin class S write k_S = |S|,
  N(S) for the common neighbourhood and d_S = |N(S)| for the common degree.
- *Twin quotient* H: its vertices are the twin classes, and S ∼ T iff some (equivalently every) vertex of S is
  adjacent to some (every) vertex of T. N(S) is a union of classes, namely of the H-neighbours of S, so
  d_S = Σ_{R∼S} k_R. H is a simple connected graph. G is bipartite iff H is: twins lie on the same side, and a
  2-colouring of H lifts to G.
- Squares: e_v² = 1_{N(v)} and e'_v² = d_v⁻¹ 1'_{N(v)} are non-zero, so A(G) and A_RW(G) are *non-degenerate*.
  Two squares are proportional iff the neighbourhoods coincide (non-zero 0/1 vectors, resp. positive multiples of
  them). Hence, for both algebras, **the proportionality classes of squares are the twin classes**; "twin-free" is
  the same as property (2LI) of the brief.

## 2. Natural bases of non-degenerate evolution algebras

### 2.1 The Boudi–Cabrera Casado–Siles Molina results, quoted exactly

BCS22 is N. Boudi, Y. Cabrera Casado, M. Siles Molina, *Natural families in evolution algebras*, Publ. Mat. 66 (2022)
159–181, arXiv 2006.14460 (numbering of the arXiv TeX; the survey's citation "[BCS22, Corollary 2.13]" matches it).
Standing hypothesis: "In this paper, 𝕂 will denote an arbitrary field" (TeX l. 254).

- Definitions 2.2 (l. 380–388): "A family of pairwise orthogonal and linearly independent vectors {u₁, …, u_r} of A
  which can be extended to a natural basis of A will be called an *extending natural family*. If {u} is an extending
  natural family for some u ∈ A, we will say that u is a *natural vector*."
- **Theorem 2.4** (l. 400–406): "Let A be an evolution 𝕂-algebra with natural basis B = {e_i}_{i∈Λ} and let u ∈ A.
  Set Supp(u) = {i₁, …, i_r}. Then: (i) If u² ≠ 0, then u is a natural vector if and only if
  rk({e_{i₁}², …, e_{i_r}²}) = 1. (ii) If u² = 0, then u is a natural vector if and only if
  e_{i₁}² = ⋯ = e_{i_r}² = 0."
- **Corollary 2.7** (l. 449–453): "Let A be a non-degenerate evolution algebra over 𝕂. Then the following assertions
  are equivalent: (i) A has a unique natural basis. (ii) There exists a natural basis B such that for any 2 different
  vectors u and v of B, u² and v² are linearly independent." ("Unique natural basis" means unique up to order and
  non-zero scalars, Definition 2.1, l. 367–371.)
- **Theorem 2.11** (l. 493–498): "Let A be an evolution algebra and let r = dim A². Then
  A = ann(A) ⊕ E₁ ⊕ … ⊕ E_r, where E₁, …, E_r are extending evolution subspaces of A satisfying erk(E_i) = 1 for
  all i and if i ≠ j, E_iE_j = 0, dim(E_i² + E_j²) = 2. Moreover, if A is non-degenerate, the decomposition is
  unique." Corollary 2.13 (l. 532–546) gives the resulting block form of a change of natural basis, and the proof of
  Proposition 2.16 (l. 560–573) uses the form b(Σα_ie_i, Σβ_ie_i) = Σλ_iα_iβ_i on a block, with uv = b(u,v)e₁².

Below we use only the direction "natural vector ⇒ rank 1" of Theorem 2.4(i). Lemma 1 proves it over every field.
The other directions need care: Theorem 2.4(i) "⇐", Proposition 2.16 and Corollary 2.17 fail in characteristic 2;
Corollary 2.7 (i)⇒(ii) needs |𝕜| ≥ 4; in Theorem 2.11, r must be the number of proportionality classes of squares,
not dim A². See §8; each claim is checked in verify.py.

### 2.2 Theorem N

Let E be an evolution algebra over an arbitrary field 𝕜 with natural basis (e_i)_{i∈I}, I finite, and assume E is
non-degenerate (e_i² ≠ 0 for all i). Write i ≈ j iff e_i² ∈ 𝕜e_j². For each class T fix u_T ≠ 0 with
e_i² = λ_iu_T (λ_i ∈ 𝕜^×, i ∈ T), put E_T = span{e_i : i ∈ T} and define on E_T the non-degenerate symmetric form
b_T(x, y) = Σ_{i∈T} λ_ix_iy_i. Then E = ⊕_T E_T, E_TE_{T'} = 0 for T ≠ T', and xy = b_T(x, y)u_T for x, y ∈ E_T.

**Lemma 1** (BCS22, Thm 2.4(i) "⇒"; any field). Let (f_a) be a natural basis of an evolution algebra E with natural
basis (e_i), and let f = f_a = Σ_iα_ie_i. Then dim span{e_i² : α_i ≠ 0} ≤ 1.

*Proof.* L(x) = fx is linear and kills the n − 1 independent vectors f_b (b ≠ a), so rank L ≤ 1. Since
L(e_i) = α_ie_i², the image of L is span{e_i² : α_i ≠ 0}. ∎

**Theorem N.** A subset B' ⊂ E is a natural basis of E iff B' = ⊔_T B'_T, where each B'_T is a b_T-orthogonal basis
of E_T. Consequently:
- (a) the decomposition E = ⊕_T E_T and the lines 𝕜u_T do not depend on the natural basis;
- (b) if all classes are singletons (2LI), the natural basis is unique up to order and scaling, over every field;
- (c) if some class has at least 2 elements and |𝕜| ≥ 4, there is a natural basis not obtained from (e_i) by
  ordering and scaling. For 𝕜 infinite, E has finitely many natural bases up to order and scaling iff 2LI holds.

*Proof.* (⇐) Vectors in different blocks multiply to 0, and inside E_T, fg = b_T(f, g)u_T = 0. (⇒) Let B' be a
natural basis and f ∈ B'. Then f² ≠ 0: otherwise f annihilates every vector of B', so f ∈ ann(E), but
ann(E) = {x : x_ie_i² = 0 ∀i} = 0. By Lemma 1 the squares e_i², i ∈ Supp(f), span a line; they are non-zero, so
Supp(f) lies in one class T and f ∈ E_T. Put B'_T = B' ∩ E_T. Independence gives |B'_T| ≤ |T|, and
Σ_T|B'_T| = |I| = Σ_T|T|, so B'_T is a basis of E_T. For f ≠ g in B'_T, 0 = fg = b_T(f, g)u_T.
(a) Each E_T is spanned by the vectors of any natural basis whose square lies in 𝕜u_T.
(b) Every E_T is a line.
(c) Take i ≠ j in one class and α ∈ 𝕜^× with λ_i + λ_jα² ≠ 0 (at most two α fail). Then e_i + αe_j and
−λ_jαe_i + λ_ie_j are b-orthogonal, and their b-norms λ_i + λ_jα² and λ_iλ_j(λ_i + λ_jα²) are non-zero. Replacing
e_i, e_j by them gives a new natural basis, and different values of α give different lines. ∎

*Answer to the survey's question "Natural bases of non-perfect algebras"* (2610.22248, l. 1245–1255, marked new).
Theorem N describes all natural bases of a non-degenerate algebra: they are the products of the sets of orthogonal
bases of the quadratic spaces (E_T, b_T). Over an infinite field the set is finite up to order and scaling iff (2LI).
It can be infinite while Aut is finite: A(C₄) over Q or R has the natural bases {ae₁ + be₃, −be₁ + ae₃, e₂, e₄}
(a² + b² ≠ 0), but |Aut A(C₄)| = 8 (Remark 3.2). This is Theorem N plus BCS22, not a new idea.

## 3. Reduction to twin classes

Let 𝕜 be any field with d_v ∈ 𝕜^× for all v. Write E = A(G) and E' = A_RW(G). By §1 both have the twin classes as
classes; in E', u'_S := 1'_{N(S)} and b'_S(x, y) = d_S⁻¹ Σ_{w∈S} x_wy_w.

**Theorem I.**
1. Let φ : A(G) → A_RW(G) be an algebra isomorphism. Then there are an automorphism τ of H with k_{τS} = k_S,
   scalars c_S ∈ 𝕜^× and matrices P_S (rows indexed by τS, columns by S) such that
   φ(e_v) = Σ_{w∈τS} (P_S)_{wv} e'_w for v ∈ S, and
   - (a) P_SᵀP_S = c_S d_S I,
   - (b) P_S𝟙 = c_T𝟙 for every T ∼ S.

   Here c_S is defined by φ(e_v)² = c_S 1'_{N(τS)} (v ∈ S).
2. Conversely, any c : V(H) → 𝕜^× and P_S ∈ M_{k_S}(𝕜) satisfying (a) and (b) define an isomorphism (take τ = id).

Hence **A(G) ≅ A_RW(G) iff the system (a), (b) has a solution over 𝕜.**

*Proof.* (1) Step 1. The vectors φ(e_v) form a natural basis of E'. By Theorem N each φ(e_v) lies in some E'_{τ(v)}.
Since φ(e_v)² = φ(e_v²) = φ(1_{N(v)}) ≠ 0 and the squares of E'_{T'} lie on the line 𝕜u'_{T'}, we get
φ(1_{N(v)}) ∈ 𝕜^×u'_{τ(v)}. For twins v, v' the vectors φ(1_{N(v)}) and φ(1_{N(v')}) coincide, hence
τ(v) = τ(v'). For non-twins the vectors 1_{N(v)} and 1_{N(v')} are not proportional, and φ is injective, hence
τ(v) ≠ τ(v'). So τ induces an injective, hence bijective, map of classes, and the
counting in Theorem N gives k_{τS} = k_S.

Step 2. For v ∈ S write φ(e_v)² = c(v)u'_{τS}. The left side equals φ(1_{N(S)}) for every v ∈ S, so c(v) =: c_S is
constant on S and non-zero:

  (∗) Σ_{R∼S} φ(1_R) = φ(1_{N(S)}) = c_S 1'_{N(τS)} = c_S Σ_{R'∼τS} 1'_{R'}.

Step 3. φ(1_R) ≠ 0 lies in E'_{τR}. Compare components in E' = ⊕E'_{R'}. If R ∼ S but τR ≁ τS, the τR-component of
the right side of (∗) vanishes, which is impossible. If R' ∼ τS is not of the form τR with R ∼ S, the R'-component of
the left side vanishes while the right side has c_S1'_{R'} ≠ 0, which is impossible. So τ maps N_H(S) onto N_H(τS),
i.e. τ ∈ Aut(H), and

  (∗∗) φ(1_R) = c_S 1'_{τR} whenever R ∼ S.

Since τ preserves sizes and adjacency, d_{τS} = Σ_{R∼S}k_{τR} = d_S.

Step 4. (a): for v, v' ∈ S, φ(e_v)φ(e_{v'}) = b'_{τS}(φ(e_v), φ(e_{v'}))u'_{τS}, which is c_Sδ_{vv'}u'_{τS}. With
b'_{τS} = d_S⁻¹·(dot product) this reads d_S⁻¹(P_SᵀP_S) = c_S I. (b): (∗∗) applied to the class S with a neighbour T
reads φ(1_S) = c_T1'_{τS}, i.e. P_S𝟙 = c_T𝟙.

(2) Each P_S is invertible since c_Sd_S ≠ 0, so φ is bijective. Basis vectors of different classes go to disjoint
supports, so their products vanish. For v, v' ∈ S,
φ(e_v)φ(e_{v'}) = Σ_w P_{wv}P_{wv'}e'_w² = d_S⁻¹(P_SᵀP_S)_{vv'}1'_{N(S)} = c_Sδ_{vv'}1'_{N(S)}, while
φ(e_v²) = Σ_{R∼S}φ(1_R) = Σ_{R∼S}Σ_{w∈R}(P_R𝟙)_we'_w = Σ_{R∼S}c_S1'_R = c_S1'_{N(S)}, using (b) for R with its
neighbour S. ∎

**Remark 3.1.** Theorem I describes every isomorphism: τ, c and the P_S are read off from φ.

**Remark 3.2 (automorphisms of A(G)).** The same proof with E' = A(G) (λ = 1 instead of d⁻¹) gives
Aut A(G) = {τ ∈ Aut(H) size-preserving, P_SᵀP_S = c_SI, P_S𝟙 = c_T𝟙 (T ∼ S)}. For C₄ the classes are {1,3} and {2,4}
and H = K₂. By Lemma 2, c₁² = c₂ and c₂² = c₁, so c³ = 1, and c = 1 over Q or R. Then the P_S are 2×2 orthogonal
matrices fixing 𝟙, i.e. I or the swap. So |Aut A(C₄)| = 2·2·2 = 8 over Q and R (24 over C).

## 4. The norm identity inside a twin class

**Lemma 2.** Let P ∈ M_k(𝕜) and α, γ ∈ 𝕜 with PᵀP = αI and P𝟙 = γ𝟙. Then k(γ² − α) = 0. If char 𝕜 = 2 and α ≠ 0,
then γ² = α for every k. In particular **γ² = α whenever α ≠ 0 and either k ≠ 0 in 𝕜 or char 𝕜 = 2.**

*Proof.* There is the polynomial identity, valid for all P, α, γ (cross terms cancel because 𝟙ᵀPᵀ𝟙 = 𝟙ᵀP𝟙):

  k(γ² − α) = 𝟙ᵀ(PᵀP − αI)𝟙 − (P𝟙 − γ𝟙)ᵀ(P𝟙 + γ𝟙),

and both terms on the right vanish. In characteristic 2, (PᵀP)_{aa} = Σ_iP_{ia}² = (Σ_iP_{ia})². So every column sum
σ satisfies σ² = α, which determines it uniquely (Frobenius is injective), and 𝟙ᵀP = σ𝟙ᵀ. Then
α𝟙 = PᵀP𝟙 = γPᵀ𝟙 = γσ𝟙, so σ² = α = γσ. Since α ≠ 0, σ ≠ 0, hence σ = γ and α = γ². (For α = 0 the claim can fail:
P = [[a, b], [a, b]] in characteristic 2 has PᵀP = 0 and P𝟙 = (a + b)𝟙. In our application α = c_Sd_S ≠ 0.) ∎

This is the step "the orthogonal freedom never produces the 0/1 pattern": inside a twin class the all-ones vector
is an eigenvector of P with eigenvalue γ = c_T, and P is a similitude with multiplier α = c_Sd_S. When 𝟙 is
anisotropic, comparing norms gives γ² = α. The orthogonal part P/γ then has to fix 𝟙 (for instance a permutation of
the twins) and has no effect on the scalars.

**Corollary 4.1 (scaling system).** Suppose every twin-class size is non-zero in 𝕜, or char 𝕜 = 2. Then
A(G) ≅ A_RW(G) iff there is c : V(H) → 𝕜^× such that
- (i) c_T = c_{T'} whenever T and T' have a common H-neighbour, and
- (ii) c_T² = c_Sd_S whenever S ∼ T.

*Proof.* (⇒) (b) gives (i), and Lemma 2 with α = c_Sd_S, γ = c_T gives (ii). (⇐) Put P_S = γ_SI, where γ_S is the
common value of c on N_H(S). ∎

So the isomorphism from (⇐) is diagonal: φ(e_v) = γ_Se'_v (v ∈ S).

**Lemma 3 (char p odd, p | k).** Let char 𝕜 = p be odd, p | k, and α, γ ∈ 𝕜^×. There is P ∈ M_k(𝕜) with PᵀP = αI
and P𝟙 = γ𝟙 iff α ∈ Sim_k(𝕜), the group of multipliers of similitudes of x₁² + … + x_k². Over a finite field,
Sim_k(F_q) = F_q^× for k even and Sim_k(F_q) = (F_q^×)² for k odd.

*Proof.* "Only if" is clear. Conversely, let RᵀR = αI and v = R𝟙. Then v ≠ 0 and vᵀv = αk = 0, and 𝟙 is isotropic too.
The map 𝕜v → 𝕜𝟙, v ↦ γ𝟙, is an isometry between totally isotropic lines. By Witt's extension theorem (non-degenerate
symmetric form, char ≠ 2) it extends to g ∈ O_k(𝕜). Then P = gR works. For finite fields: if k is odd,
det(R)² = α^k forces α to be a square, and conversely R = √α·I. If k is even, use blocks [[a, −b], [b, a]] with
a² + b² = α (every element of F_q is a sum of two squares). ∎

## 5. Proof of Theorem A

Let char 𝕜 = 0. All k_S ≠ 0 in 𝕜, so Corollary 4.1 applies. H is connected (G is).

*Non-bipartite H.* Any two vertices of H are joined by a walk of even length T₀, T₁, …, T_{2m}. Consecutive even-indexed
vertices have the common neighbour T_{2i+1}, so by (i) c is constant, c ≡ a. Then (ii) gives a² = ad_S, so d_S = a
for every class S. The degrees are integers and 𝕜 ⊇ Q, so G is regular.

*Bipartite H* (sides 𝒳, 𝒴; G is bipartite with sides X = ∪𝒳 and Y = ∪𝒴). Even walks connect any two classes on the
same side, so c = x on 𝒳 and c = y on 𝒴. By (ii), d_S = y²/x for all S ∈ 𝒳 and d_S = x²/y for all S ∈ 𝒴. So the
degree is constant on X (=: d₁) and on Y (=: d₂), and G is biregular. Eliminating y = x²/d₂ from y² = xd₁ gives
x³ = d₁d₂², so d₁/d₂ = (x/d₂)³ ∈ (𝕜^×)³.

*Conversely.* If G is d-regular, c ≡ d satisfies (i) and (ii), and φ = d·id. If G is biregular with d₁/d₂ = t³,
t ∈ 𝕜^×, put x = td₂ and y = t²d₂. Then y² = t⁴d₂² = xd₁ and x² = t²d₂² = yd₂, so (ii) holds. The isomorphism is
φ(e_v) = γ_Se'_v, i.e.

  φ(e_v) = t²d₂·e'_v (v ∈ V₁),  φ(e_v) = td₂·e'_v (v ∈ V₂).

Check: for v ∈ V₁, φ(e_v)² = t⁴d₂²d₁⁻¹1'_{N(v)} = td₂1'_{N(v)} = φ(e_v²); for v ∈ V₂ it is symmetric. With real
cube roots this is CRR's map (1710.10516, l. 994–1002: (d₁²d₂)^{1/3} on V₁ and (d₁d₂²)^{1/3} on V₂). ∎

**All isomorphisms in characteristic 0.** By Theorem I and Lemma 2, P_SᵀP_S = γ_S²I and P_S𝟙 = γ_S𝟙, so
P_S = γ_SQ_S with Q_S orthogonal and Q_S𝟙 = 𝟙. Every isomorphism is therefore e_v ↦ γ_S Σ_{w∈τS}(Q_S)_{wv}e'_w, with
τ a size-preserving automorphism of H, Q_S ∈ O_{k_S}(𝕜) fixing 𝟙, and γ_S = c_T (T ∼ S). Here c ≡ d if G is regular
and non-bipartite, and c = x on the V₁-classes, c = y = x²/d₂ on the V₂-classes with x³ = d₁d₂² if G is bipartite
biregular. That gives 3 choices of x over C and 1 over R. Twins can be rotated, but the scalars are rigid.

**Examples** (brief, items (a) and (c)).
- P₃: V₁ = {centre}, d₁ = 2, d₂ = 1. The isomorphism needs t³ = 2, and the centre's coefficient is t² with
  (t²)³ = 4 (the "e³ = 4" of the full Gröbner basis).
- K_{1,3}: t³ = 3, centre coefficient t², (t²)³ = 9.
- Diamond K₄ − e: the classes {u}, {v}, {x, y} form a triangle H with degrees 3, 3, 2. H is non-bipartite and G is not
  regular, so A(G) ≇ A_RW(G) over every field of characteristic 0. The brief's ansatz computation is the case
  P_{x,y} = 2×2 of Theorem I.

## 6. Positive characteristic

Let char 𝕜 = p > 0 and p ∤ d_v for all v (otherwise A_RW(G) is undefined).

**Theorem P.**
1. *p = 2.* Every d_v is odd, so d_v = 1 in 𝕜 and A_RW(G) = A(G).
2. *p odd.* A(G) ≅ A_RW(G) iff there is c : V(H) → 𝕜^×, constant on each N_H(S) with value γ_S, such that
   γ_S² = c_Sd_S for every class with p ∤ k_S, and c_Sd_S ∈ Sim_{k_S}(𝕜) for every class with p | k_S.
3. If no twin class has size divisible by p, then A(G) ≅ A_RW(G) iff the degrees are all congruent mod p, or G is
   bipartite with degrees congruent mod p on each side (residues δ₁, δ₂) and δ₁/δ₂ ∈ (𝕜^×)³.

*Proof.* (1) is clear. (2): Theorem I, then Lemma 2 for classes with p ∤ k_S, and Lemma 3 for classes with p | k_S;
the constancy of c on N_H(S) is (b). (3): the proof of §5 with the degrees read in F_p ⊂ 𝕜. ∎

**What changes, with examples (all checked by brute force in verify.py).**
- *"Only if" fails in every characteristic p.* The double star S(p, p) (two adjacent centres, p leaves each) has all
  degrees ≡ 1 mod p, so A_RW = A. It is neither regular nor biregular. For p = 2 take S(2, 2), degrees 3,3,1,1,1,1.
- *"If" can fail.* P₃ over F₇ is biregular with ratio 2, which is not a cube mod 7, so the algebras are not
  isomorphic. Over F₅ they are.
- *Twin classes of size divisible by p.* K_{1,1,3} (u ∼ v, both adjacent to x₁, x₂, x₃) has degrees
  4, 4, 2, 2, 2 ≡ 1, 1, −1, −1, −1 (mod 3) and H is a triangle. Theorem P(2) gives c ≡ 1 and requires
  c_Xd_X = −1 ∈ Sim₃(𝕜) = (𝕜^×)², i.e. √−1 ∈ 𝕜. So the algebras are **not isomorphic over F₃ but isomorphic over
  F₉ = F₃(i)**, via

    φ(e_u) = e'_u, φ(e_v) = e'_v, φ(e_{x_a}) = Σ_b P_{ba}e'_{x_b},  P = [[0, −1, −1], [1, i, −i], [1, −i, i]]

  (PᵀP = −I, P𝟙 = 𝟙). No isomorphism keeps the twins on lines: that would need γ_X² = c_Xd_X, i.e. 1 = −1.
  The tree T6 (claw u; x₁, x₂, x₃ plus the path u – z – w) behaves the same way: not isomorphic over F₃, isomorphic over
  F₉ via e_u ↦ −e'_u, e_z ↦ e'_z, e_w ↦ −e'_w and the same P on the x's. Its side {x₁, x₂, x₃, z} has degrees
  1, 1, 1, 2, not constant mod 3.

## 7. Homomorphisms: the "moreover" clause of CRR

CRR's conjecture (l. 963–964) and Theorem 2.3 also assert: *if A_RW(G) ≇ A(G), the only homomorphism between them is
the null map.*

**Proposition H.**
1. If 𝕜 is formally real (−1 is not a sum of squares; e.g. Q or R), every non-zero algebra homomorphism
   A_RW(G) → A(G) or A(G) → A_RW(G) is an isomorphism. With Theorem A, CRR's whole conjecture is therefore true over R.
2. If A is singular over 𝕜 and 𝕜 is algebraically closed (e.g. C), there are non-zero non-injective homomorphisms in
   both directions. So over C the "moreover" clause fails for every singular G that is neither regular nor biregular,
   e.g. the diamond and the double stars C_{m,n} (m, n ≥ 2) of CRRV.

*Proof.* (1) Over a formally real field, x² = Σ_vx_v²1_{N(v)} = 0 in A(G) forces x = 0: the coefficient at w is
Σ_{v∈N(w)}x_v², and every v lies in some N(w). The same holds in A_RW(G), with weights d_v⁻¹ > 0. Let f ≠ 0 be a
homomorphism and Z = {i : f(e_i) = 0}. The non-zero f(e_i) are pairwise orthogonal with non-zero squares, hence
linearly independent: multiply a linear relation by f(e_j). For i ∈ Z, 0 = f(e_i²) is a combination with non-zero
coefficients of the f(e_ℓ), ℓ ∈ N(i), so N(i) ⊆ Z. By connectivity Z = ∅, so f is injective, hence bijective.
(2) Take y ∈ ker A \ {0}, x_v with x_v² = y_v, u = Σx_ve_v. Then u² = Σ_wy_w1_{N(w)} = Σ_w(Ay)_we_w = 0. For
β ∈ ker A \ {0}, f(e_i) = β_iu satisfies f(e_ie_j) = 0 = f(e_i)f(e_j) and
f(e_i²) = (D⁻¹Aβ)_iu = 0 = β_i²u². The other direction uses x_v² = d_vy_v. ∎

Example: for twins x, y, e_x ↦ e_x + ie_y, e_y ↦ −(e_x + ie_y), all other e_v ↦ 0 is a homomorphism
A_RW(G) → A(G) over C (verified for the diamond and C_{2,2}).

## 8. Remarks on the literature (side findings; each checked in verify.py)

1. **BCS22 in characteristic 2.** Take 𝕜 = F₂ and e₁² = e₂² = e₃² = e₄, e₄² = e₁. Then u = e₁ + e₂ + e₃ has
   u² = e₄ ≠ 0 and rank-1 support squares, yet lies in no natural basis: u^⊥ ∩ E_T is an alternating plane. So
   Theorem 2.4(i) "⇐" fails. In the 5-dimensional analogue (e₁² = … = e₄² = e₅, e₅² = e₁), e₁ + e₂ + e₃ and e₄ are
   orthogonal natural vectors in no common natural basis, so Proposition 2.16 and Corollary 2.17 fail. The survey
   restates Corollary 2.17 over every field (2610.22248, l. 332–333). All are fine in characteristic ≠ 2, where
   non-degenerate quadratic spaces have orthogonal bases.
2. **BCS22 Corollary 2.7** (i)⇒(ii) needs |𝕜| ≥ 4; its proof says "since 𝕂 has strictly more than three elements".
   Counterexamples: over F₃, e₁² = e₃, e₂² = −e₃, e₃² = e₁; over F₂, e₁² = e₂² = e₃, e₃² = e₁. Both have a unique
   natural basis without (2LI).
3. **BCS22 Theorem 2.11**: "r = dim A²" should read "r = number of proportionality classes of squares". For A(P₅)
   over Q, dim A² = rank = 4, but there are 5 classes and a unique natural basis.
4. **CRRV (2405.12341)**: Corollary "coro:column0" (l. 319–330) deduces t_{ik} = 0 from Σ_{k∈N(r)}t_{ik}² = 0, which is
   valid over R and false over C; Proposition H(2) gives a non-zero homomorphism with a zero column for C_{2,2}. Over
   R their results stand, and Proposition H(1) makes all their "moreover" statements automatic.
5. **Vidal–Cadavid–Rodriguez (arXiv 2111.07399), Proposition "twin construction"** (l. 863–888) claims a non-injective
   surjective homomorphism A_RW(G) → A(G) when G/Π is regular or biregular. The proposed map sends twins i ≠ j to the
   same vector f(e_{i_k}) with f(e_{i_k})² ≠ 0, so f(e_i)f(e_j) ≠ 0 = f(e_ie_j). Moreover, for finite G a surjective
   linear map A_RW(G) → A(G) is injective, and over R no non-zero non-injective homomorphism exists (Proposition H(1)).
   Example: the diamond (G/Π a triangle).

## 9. Computer verification (verify.py)

Default run (≈ 25 s, < 150 MB; hard alarm 60 s), all passing:
- **[N]+[P], finite fields, no theory used.** All natural bases of A(G) and A_RW(G) are enumerated by brute force:
  connected G with n ≤ 5 over F₃, n ≤ 4 over F₅ and F₇ (30 graph–field pairs). Their numbers equal Theorem N's
  product over twin classes. A(G) ≅ A_RW(G) is decided through all natural bases of A_RW(G) and agrees with
  Theorem P. K_{1,1,3} and T6 have 4 natural bases each and are not isomorphic over F₃. The explicit
  F₉-isomorphisms (twin mixing) are verified, and A_RW = A for S(2,2) in characteristic 2.
- **[G], full GL_n systems over Q** (n² unknowns plus 1/det, sympy Gröbner). P₃: the lex basis contains p₀₀³ − 4.
  K_{1,3}: p₀₀³ − 9 lies in the ideal. The diamond has Gröbner basis [1].
- **[L]** The identity of Lemma 2 for k = 1, …, 6 (every class size occurring for n ≤ 7), and the elimination ideal
  (α − γ²) for k = 2.
- **[A], characteristic 0, all 995 connected graphs with 2–7 vertices** (networkx graph atlas). The twin-class system
  of Corollary 4.1 is solved over C, R and Q (sympy Gröbner plus exact solutions) and equals Theorem A's criterion in
  all cases: 15 regular graphs, 9 biregular non-regular ones (all complete bipartite) and 971 others. Isomorphic over
  C and R: 24; over Q: 15. Every predicted isomorphism is built and checked exactly. Extra graphs: K_{1,8} and
  K_{2,16} (isomorphic over Q), the subdivided K₄ and K₅ (over R, not over Q), and the Petersen graph.
- **[H]** The homomorphisms of Proposition H(2) and the literature remarks of §8.

--full (6–8 min) adds:
- the brute force for n = 5 over F₅ and F₇ (30 graphs each);
- K_{1,1,3} and T6 over F₉ by brute force (120 natural bases each, matching Theorem N; isomorphic);
- the full GL₄ systems of P₄ and the paw ([1]) and of C₄ and K₄ (non-empty).

All pass.

## 10. Novelty and sources

- The question is open in the survey 2610.22248 (15 Oct 2026; Question q:rw, "Known for non-singular graphs
  [CRR20, CRR21]; supporting families in [CRRV]").
- CRR (arXiv 1707.05897, 1710.10516) prove the non-singular case. Their key step (Prop. theo:principal, l. 1237–1262)
  needs an invertible A to force monomial isomorphisms. CRRV (arXiv 2405.12341; J. Algebra Appl. 2027) treat
  caterpillars, tadpoles T_{4,m} and one twin-free graph by direct computation; their l. 237 suggests using the twin
  partition but gives no general argument. Vidal–Cadavid–Rodriguez (2111.07399) treat the Hilbert analogue,
  non-singular case only.
- arXiv searches (web search UI, 22 Oct 2026; the API was rate-limited): all 141 title matches of "evolution algebra",
  plus "evolution algebras associated to graphs", "random walk evolution algebra", "evolution algebra graph
  biregular", "evolution algebra regular graph isomorphic", "Markov evolution algebra graph", "twin partition
  evolution algebra" and "natural families evolution algebras". Abstracts were read for 2505.15604 (automorphisms of
  A_Γ, A_RW(Γ)), 2605.12672, 2607.15154, 2404.08752 and 2311.00102. Nothing settles the singular case, notes the
  cube-root obstruction over Q, or describes isomorphisms through twin classes.
- So, as far as we can tell, Theorem A (with the Q counterexample), Theorem P and Proposition H are new. Theorem N is
  BCS22 (Thm 2.4 + Thm 2.11 + the form of Prop. 2.16) restated with a short proof.

Sources (TeX downloaded via arxiv.org/e-print): 2610.22248, 2006.14460, 1707.05897, 1710.10516, 2405.12341,
2111.07399, 2505.15604.

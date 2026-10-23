# Prefix reversals on the boundary m + k = n − 1: a complete answer below the diagonal

Files: [verify.py](verify.py) (machine check, 31 s, 62 MB), [REPORT.md](REPORT.md).
Sources: Blanco–Golubyatnikov–Konstantinova–Maslova–Nikiforov, *Generating the symmetric group by
three prefix reversals*, arXiv 2511.16959v3 (DMTCS 28:1, 2026), Conjecture 4 (TeX l. 884–886);
Lei–Zhao–Zhu–Wu, arXiv 2610.30892v1, Theorem 3.1, the Lemma at TeX l. 277 and the remark at TeX l. 268.
Known pieces are credited where they are used and listed in §8 (see also [NOVELTY.md](NOVELTY.md)).

## Notation

Positions 1..n. For 2 ≤ i ≤ n, r_i(x) = i+1−x for x ≤ i and r_i(x) = x for x > i. Products act
right to left. G = G(n;m,k) = ⟨r_n, r_m, r_k⟩ with 2 ≤ k < m < n.
On the line n = m + k + 1 put α = k+1, β = m+1, N = n+1 = α+β, g = gcd(α, β) = gcd(m+1, k+1).
[t]_q = (q^t − 1)/(q − 1).

## Main results

**Theorem A (answer below the diagonal).** Let 2 ≤ k < m < n and m + k < n. Then
⟨r_n, r_m, r_k⟩ = S_n if and only if

  n is odd, m + k = n − 1, gcd(m+1, k+1) = 1, (k+1) ∤ m, and not (n ≡ 1 mod 4 and 4 | k).

Equivalently: n is odd, α + β = n + 1 with α, β coprime (both are then odd, because α + β = n + 1 is
even; coprimality alone would not force this), β ≢ 1 (mod α), and not α ≡ β ≡ 1 (mod 4).

**Theorem B (the whole line n = m + k + 1, 2 ≤ k < m, any parity).**
1. The orbits of G are O_c = {x : x ≡ ±c (mod g)}, 0 ≤ c ≤ g/2; there are ⌊g/2⌋ + 1 of them.
   G is transitive iff g = 1. (That g = 1 is necessary for G = S_n is known: 2610.30892, Lemma at
   TeX l. 277; new here is the exact orbit decomposition.)
2. n even, g = 1: G is transitive and imprimitive with two blocks (known: 2610.30892, Thm 3.1, and
   earlier in Lean, see §1).
3. n odd, g = 1, α | n: G is imprimitive, the intervals [iα+1, (i+1)α] are blocks, and
   G ≅ D_{2α} ≀ D_{2q} with q = n/α, |G| = (2α)^q · 2q. The nontrivial block sizes are exactly the
   divisors L of n with L | α or α | L.
4. n odd, g = 1, α ∤ n: G is primitive and contains A_n; G = A_n if n ≡ 1 (mod 4) and 4 | k
   (equivalently α ≡ β ≡ 1 mod 4), otherwise G = S_n.

**Corollary (Conjecture 4).** Read literally, Conjecture 4 of 2511.16959v3 (the published DMTCS
text; v1 and v2 state the parity-split form) says m + k < n ⇒ G < S_n.
It fails exactly at the triples of B.4 with G = S_n: infinitely many, e.g. (n, n−3, 2) for every
n ≡ 1 (mod 6) and (n, n−5, 4) for every n ≥ 11 with n ≡ 3, 7, 11 (mod 20); the smallest are (7,4,2), (11,6,4), (13,10,2), (15,8,6),
(17,10,6), (19,16,2), (19,12,6), (19,10,8); there are 318 with n ≤ 101. The parity-split form
(n even: m+k < n; n odd: m+k < n−1), proved in 2610.30892, is therefore sharp, and Theorem B
decides the one remaining line. Theorem A is consistent with Theorem 1 of 2511.16959 (k = 2,
n odd: (n, n−3, 2) generates iff n ≡ 1 mod 3). Its condition g = 1 is necessary by the Lemma of
2610.30892 at TeX l. 277 (G = S_n ⇒ gcd(n−m, n−k) = 1, and here n−m = α, n−k = β), whose special
case n−m | n+1 is Proposition 2 of 2511.16959.

## 1. Off the line (known)

Known: 2610.30892, Thm 3.1. The same result was proved earlier in Lean in GitHub
the-omega-institute/trureturing, files D5/S0/CayleyGrowth/PrefixReversalTriple{Odd,Even}NonGeneration.lean
(2026-09-13). Neither source treats the odd line. The short argument, for completeness:

Let X be the graph on [1,n] with edges {x, r_q(x)} (q ∈ {n,m,k}, r_q(x) ≠ x); its components are
the G-orbits and |E(X)| = ⌊n/2⌋ + ⌊m/2⌋ + ⌊k/2⌋. If n is odd and m+k ≤ n−2 then |E(X)| ≤ (n−1)/2 +
(n−2)/2 < n − 1, so G is intransitive. If n is even and m+k ≤ n−1 then |E(X)| ≤ n − 1; if X is
connected it is a tree, and colouring each vertex by the parity of the number of r_n-edges on its
path to a fixed root gives two classes that r_m, r_k preserve and the fixed-point-free r_n swaps;
so G < S_n. (Transitive examples exist, e.g. (8,4,2); verify.py meets 139 of them.)
Hence for m + k < n only the line n = m + k + 1 with n odd can give S_n.

## 2. The n-cycle and the orbits (Theorem B.1)

**Lemma 1.** On the line, P := r_k r_n r_m is P(x) = x + α (x ≤ m), P(m+1) = α, P(x) = x − β (x ≥ m+2).
Identify [1,n] with Z/N ∖ {0}. Then P is the first-return map to Z/N ∖ {0} of the rotation
R(x) = x + α of Z/N: P(x) = R(x) if R(x) ≠ 0, and P(β) = R²(β) = α.

*Proof.* x ≤ m: r_m(x) = m+1−x, r_n of it = n−m+x = x+α ≥ k+2, fixed by r_k. x = m+1: r_n(m+1) = α,
fixed by r_k. x ≥ m+2: r_n(x) = n+1−x ≤ k, r_k of it = k+1−n−1+x = x−β. Since x+α ≡ x−β (mod N),
P agrees with R except at x = β, where R(β) = N ≡ 0. ∎

*Proof of B.1.* The R-orbits are the cosets of gZ/N, so the P-orbits are the classes
{x ∈ [1,n] : x ≡ c (mod g)}. r_n(x) = N − x ≡ −x, so ⟨P, r_n⟩ has orbits O_c. Each generator maps
x to x or to C − x with C ∈ {N, β, α} ≡ 0 (mod g), so every O_c is G-invariant. ∎
(For g = 1, P is an n-cycle. Necessity of g = 1 for G = S_n is the Lemma of 2610.30892 at TeX
l. 277; the invariant set of multiples of g is the device of 2511.16959 Prop. 2. New here: the
exact orbits.)

## 3. Cyclic coordinates

Assume g = 1 and let θ = α⁻¹ mod N. Put φ(x) = θx mod N; φ is a bijection [1,n] → [1,n].

**Lemma 2.** In the coordinate j = φ(x): P is j ↦ j+1 (j < n), n ↦ 1; r_n is j ↦ n+1−j;
r_k is j ↦ n+2−j on E := φ([1,k]) and the identity off E. (E = {⌈sN/α⌉ : 1 ≤ s ≤ k}, 1 ∉ E.)

*Proof.* φ(x+α) = φ(x)+1, φ(x−β) = φ(x+α−N) = φ(x)+1, φ(β) = −θα ≡ −1, i.e. n, and φ(α) = 1.
φ(N−x) = −φ(x), i.e. N−j. φ(α−x) = 1−φ(x), i.e. N+1−j, for x ∈ [1,k] (j ≠ 1 since φ(x) = 1 iff x = α). ∎

B.2 in these coordinates: for n even, r_k preserves the parity of j (n+2−j ≡ j), r_n and P flip
it, so the classes of j mod 2 are two blocks.

## 4. Blocks for n odd

From now on n is odd and g = 1 (so α, β are odd, 3 ≤ α < β).

**Lemma 3.** Every nontrivial G-invariant partition is 𝒞_d = {C_c : c mod d}, C_c = {j : j ≡ c mod d},
for a divisor d of n with 1 < d < n (so d, L := n/d are odd and ≥ 3). 𝒞_d is G-invariant iff

 (★_d) every class C_c with c ≢ 1 (mod d) satisfies C_c ⊆ E or C_c ∩ E = ∅.

*Proof.* ⟨P⟩ is regular and cyclic, so the block through a point is an orbit of a subgroup ⟨P^d⟩:
a residue class. P and r_n (j ↦ 1−j mod d) permute the classes. r_k maps C_c ∩ E into C_{2−c} and
fixes C_c ∖ E; as d is odd, 2−c ≡ c iff c ≡ 1. ∎

**Lemma 4 (classes in x-coordinates).** Let Δ = αd mod N, y_s = sΔ mod N and
R_i = {y_s : iL < s ≤ (i+1)L} (0 ≤ i < d). Then φ⁻¹(C_c) = R_{d−c} for 1 ≤ c ≤ d. So
φ⁻¹(C_0) = R_0 = {Δ, 2Δ, …, LΔ}, the exempt class is φ⁻¹(C_1) = R_{d−1} = {−Δ, …, −LΔ}, and
LΔ = αn ≡ −α (mod N). Also gcd(L, N) = gcd(Δ, N) = 1 and y_s ≠ 0 for 1 ≤ s ≤ n.

*Proof.* dL = n ≡ −1 gives d⁻¹ ≡ −L, so j = c + td (0 ≤ t < L) has x = αj ≡ (t − cL)Δ, and
t − cL ≡ (d−c)L + 1 + t with (d−c)L + 1 + t ∈ ((d−c)L, (d−c+1)L]. ∎

Write I = φ⁻¹(E) = [1, α−1]. In x-coordinates (★_d) says: every R_i with i ≤ d−2 lies in I or misses I.

**Proposition 5 (primitivity).** If (★_d) holds for a divisor 1 < d < n, then α | n.
Hence G is primitive whenever α ∤ n (equivalently (k+1) ∤ m).

*Step 1 (I ⊆ R_{d−1} ⇒ α | L).* Since α ∉ I ⊆ {−tΔ : 1 ≤ t ≤ L} and −LΔ ≡ α, every x ∈ I is −τ(x)Δ
with 1 ≤ τ(x) ≤ L−1, and τ(x) ≡ x·τ(1) (mod N). By induction xτ(1) ≤ L−1 as integers for x ≤ α−1
(the next multiple is < 2L < N, so it equals its residue). Then α·τ(1) < 2L < N and α·τ(1) ≡ L,
so ατ(1) = L. ∎

*Step 2 (some R_i ⊆ I, i ≤ d−2 ⇒ L | α or L | β).* Consecutive points of R_i differ by Δ mod N
and lie in [1, α−1], so R_i is an integer progression with step Δ' (the representative of Δ in
(−N/2, N/2)) and (L−1)|Δ'| ≤ α−2. From LΔ ≡ −α: LΔ' + α ≡ 0 (mod N) and
|LΔ' + α| ≤ (3/2)(α−2) + α < 2N, LΔ' + α > −N; so LΔ' = −α or LΔ' = N − α = β. ∎

*Step 3 (L | β ⇒ (★_d) fails).* Then Δ ≡ μ := β/L (because Lμ = β ≡ −α ≡ LΔ and gcd(L,N) = 1),
and R_0 = {μ, 2μ, …, Lμ = β} (no wrap). If μ ≤ α−1, R_0 meets I (at μ) and misses it (at β), and
R_0 is not exempt (d ≥ 3). If μ ≥ α, no R_i fits in I (Step 2 needs |Δ'| = μ ≤ (α−2)/(L−1), note
μ ≤ β/3 < N/2), so (★_d) forces I ⊆ R_{d−1}, hence α | L | β by Step 1: impossible as g = 1, α ≥ 3. ∎

*Step 4 (L | α and (★_d) ⇒ β ≡ 1 mod α).* Put λ = α/L, so Δ ≡ −λ and R_i = {v_i − tλ : 1 ≤ t ≤ L}
with v_i = y_{iL} ≡ −iα (taken in [0, N−1]). A direct check (β > α) gives: R_i ⊆ I iff
v_i ∈ [α+1, α+λ−1]; R_i ∩ I = ∅ iff v_i ∈ [0, λ−1] ∪ [2α, N−1] (v_i = λ would put 0 in R_i);
otherwise R_i is split. The sequence v_i starts v_0 = 0 and moves by −α if v_i ≥ α, by +β if
v_i < α; v_{d−1} = α+λ and v_d = λ (as dα = nλ ≡ −λ). Call a lap a stretch from a value
ρ_{t−1} + β down by steps α to the next value ρ_t < α; then ρ_t = tβ mod α, the lap passes through
ρ_t + α ∈ [α, 2α−1], and its other values are ≥ 2α. Counting steps, the index of ρ_t is
(tβ − ρ_t)/α + t, so v_d = λ = ρ_T with d = (TN − λ)/α, i.e. TN = dα + λ = λN: T = λ. The laps
t ≤ λ−1 end before index d−1, so (★_d) needs ρ_t + α ∈ [α+1, α+λ−1], i.e. tβ mod α ∈ [1, λ−1],
for 1 ≤ t ≤ λ−1. If λ = 1 then L = α and L | n − α = β − 1. If λ ≥ 2, r := β mod α ∈ [1, λ−1];
if r ≥ 2 the least t with tr ≥ λ has t ≤ λ−1 and tr ≤ 2λ−2 < α, so tr mod α = tr ∉ [1, λ−1].
So r = 1. ∎

*Conclusion.* If some R_i (i ≤ d−2) lies in I, Steps 2–4 give L | α and α | β − 1 = m, so α | n.
Otherwise (★_d) makes every such R_i miss I, so I ⊆ R_{d−1} and α | L | n by Step 1. ∎

**Lemma 6 (the imprimitive case, B.3).** Let α | n, q = n/α (odd, ≥ 3), B_i = [iα+1, (i+1)α]. In
coordinates x = iα + u (u ∈ Z/α): r_n(i,u) = (q−1−i, 1−u); r_m(i,u) = (q−2−i, 1−u) for i ≤ q−2 and
r_m fixes B_{q−1}; r_k(0,u) = (0, −u) and r_k fixes the other blocks. So the B_i are blocks and
G ≤ D_{2α} ≀ D_{2q} (dihedral on positions u, dihedral ⟨rev_q, rev_{q−1}⟩ on blocks). Conversely
T = r_n r_m maps (i,u) ↦ (i+1, u) for i ≤ q−2 and (q−1,u) ↦ (0, 1−u), so T^q acts on every block
as u ↦ 1−u; then T^q r_k T^{−q} acts on B_0 as u ↦ 2−u and trivially elsewhere, and with r_k it
generates the full D_{2α} on B_0; conjugating by T^i (B_0 → B_i, u ↦ u) gives D_{2α} on each block.
So the kernel of the block action is D_{2α}^q and G = D_{2α} ≀ D_{2q}. ∎
Block sizes: the proof of Proposition 5 shows that (★_d) forces L | α or α | L. Conversely, if
L | α then β ≡ 1 (mod α) gives ρ_t = t ∈ [1, λ−1] in Step 4 (that step is an equivalence), so (★_d)
holds; if α | L then θ ≡ −q (mod N), so E = {N − qx : 1 ≤ x ≤ k} ⊆ C_1 because d | q. (Machine-checked.)

## 5. The alternating group

**Lemma 7 (fixed points of collineations).** Let d ≥ 2, (d, q) ≠ (2, 4), n = [d]_q. A nonidentity
element of PΓL_d(q) fixes at most n/2 points (and at most n/2 hyperplanes) of PG(d−1, q).

*Proof.* Let f(v) = Av^σ induce f̄ ≠ 1. (i) σ = 1: the fixed points are the projective points of the
eigenspaces V_λ, which are independent and proper; since [a]_q + [b]_q ≤ [a+b−1]_q + 1 for a,b ≥ 1,
|Fix| ≤ [d−1]_q + 1 for q ≥ 3 and ≤ [d−1]_2 for q = 2 (only eigenvalue 1); both are ≤ n/2
because n = q[d−1]_q + 1. (ii) σ ≠ 1: let GF(q_0) = Fix(σ), q = q_0^e, e ≥ 2, and N the norm to
GF(q_0). f^e is the linear map B = AA^σ⋯A^{σ^{e−1}}; if f(v) = λv then Bv = N(λ)v. For κ ∈ GF(q_0)^*
pick λ_κ with N(λ_κ) = κ; as c ↦ c^σ/c maps GF(q)^* onto ker N (Hilbert 90), every fixed point with
N(λ) = κ has a representative in W_κ = {v : f(v) = λ_κ v}, a GF(q_0)-space. GF(q_0)-independent
vectors of W_κ are GF(q)-independent (apply f to a shortest relation and subtract), and two vectors
of W_κ give the same point iff they are GF(q_0)-proportional. So those fixed points number
[w_κ]_{q_0}, w_κ = dim W_κ ≤ dim ker(B − κ), and Σ w_κ ≤ d; hence |Fix| ≤ [d]_{q_0}. Finally
[d]_q ≥ [d]_{q_0²} = [d]_{q_0}(q_0^d + 1)/(q_0 + 1) ≥ 2[d]_{q_0} unless d = 2, q_0 = 2, and then
[2]_2 = 3 ≤ (q+1)/2 unless q = 4. Hyperplanes: apply the inverse-transpose automorphism. ∎
(verify.py confirms the bound by brute force for PG(1,q), q ∈ {4,8,9,16,27}, PG(2,2), PG(2,3),
PG(3,2) and the semilinear elements of PΓL_3(4); PG(1,4) attains 3 > 5/2.)

**Proposition 8.** If n is odd, g = 1 and α ∤ n, then G ⊇ A_n.

*Proof.* G is primitive (Prop. 5) and contains the n-cycle P. If A_n ⊄ G, Jones's classification
(G. A. Jones, Bull. Aust. Math. Soc. 89 (2014) 159–165, Thm 1.2, case k = 0, which rests on
Jones, J. Group Theory 5 (2002) 403–407, and CFSG) leaves: (a) n = p prime, G ≤ AGL_1(p);
(b) PGL_d(q) ≤ G ≤ PΓL_d(q), n = [d]_q, d ≥ 2; (c) L_2(11), M_11, M_23 with n = 11, 11, 23. But r_k ≠ 1
fixes n − k = m + 1 > n/2 points (m ≥ k+2 ≥ 4). Nonidentity elements of AGL_1(p) fix ≤ 1 point;
(b) contradicts Lemma 7 (n ≥ 7, so (d,q) ≠ (2,4)); in (c) involutions fix 3, 3, 7 < m+1 points
(and n = 11, 23 are machine-checked). ∎
For n ≤ 101 verify.py also proves G ⊇ A_n without CFSG: each primitive case contains a p-cycle,
p prime ≤ n − 3 (a power of a product-replacement element), and Jordan's theorem applies.

**Parity.** r_i is a product of ⌊i/2⌋ transpositions. With m, k even and m + k = n − 1, all three
generators are even iff (n−1)/2, m/2, k/2 are all even iff n ≡ 1 (mod 4) and 4 | k. This gives B.4,
and with §1 and B.1–B.3 it gives Theorem A. ∎

## 6. Machine check (verify.py, n ≤ 101)

All 2352 line triples (both parities): Lemma 1, the orbit formula B.1. All 1176 odd-n line triples:
classification = Theorem B (692 intransitive, 56 imprimitive, 110 A_n, 318 S_n). For the 484 with
g = 1: Lemma 2, (★_d) versus the action of r_k for every divisor d, the computed block systems
(Atkinson closure) equal {d : (★_d)}; block sizes L | α or α | L when α | n; interval blocks;
Jordan certificates for all 428 primitive cases (primes 2..89); parity. Off the line: all triples
with m+k ≤ n−2, n ≤ 60 (plus 1500 random ones with n ≤ 101) and all even-n line triples: G < S_n.
sympy (optional): |G| for every line triple with n ≤ 25; |G| = (2α)^q·2q for the 21 imprimitive
cases with n ≤ 51. Lemma 7 brute force as above.

## 7. Table

Line triples (m, k), k even (odd k gives g ≥ 2): S = S_n, A = A_n, W = imprimitive D≀D,
Ir = intransitive with r orbits.

| n | (m,k) with m+k = n−1: type |
|---|---|
| 7 | (4,2) S |
| 9 | (6,2) W |
| 11 | (8,2) I2, (6,4) S |
| 13 | (10,2) S, (8,4) A |
| 15 | (12,2) W, (10,4) W, (8,6) S |
| 17 | (14,2) I2, (12,4) A, (10,6) S |
| 19 | (16,2) S, (14,4) I3, (12,6) S, (10,8) S |
| 21 | (18,2) W, (16,4) A, (14,6) W, (12,8) A |
| 23 | (20,2) I2, (18,4) S, (16,6) S, (14,8) I2, (12,10) S |
| 25 | (22,2) S, (20,4) W, (18,6) S, (16,8) A, (14,10) S |
| 27 | (24,2) W, (22,4) S, (20,6) I4, (18,8) W, (16,10) S, (14,12) S |
| 29 | (26,2) I2, (24,4) I3, (22,6) S, (20,8) I2, (18,10) S, (16,12) A |
| 31 | (28,2) S, (26,4) S, (24,6) S, (22,8) S, (20,10) S, (18,12) S, (16,14) S |
| 33 | (30,2) W, (28,4) A, (26,6) S, (24,8) A, (22,10) W, (20,12) A, (18,14) S |
| 35 | (32,2) I2, (30,4) W, (28,6) W, (26,8) I5, (24,10) S, (22,12) S, (20,14) I2, (18,16) S |
| 37 | (34,2) S, (32,4) A, (30,6) S, (28,8) A, (26,10) S, (24,12) A, (22,14) S, (20,16) A |
| 39 | (36,2) W, (34,4) I3, (32,6) S, (30,8) S, (28,10) S, (26,12) W, (24,14) I3, (22,16) S, (20,18) S |
| 41 | (38,2) I2, (36,4) A, (34,6) I4, (32,8) I2, (30,10) S, (28,12) A, (26,14) I2, (24,16) A, (22,18) S |
| 43 | (40,2) S, (38,4) S, (36,6) S, (34,8) S, (32,10) I6, (30,12) S, (28,14) S, (26,16) S, (24,18) S, (22,20) S |
| 45 | (42,2) W, (40,4) W, (38,6) S, (36,8) W, (34,10) S, (32,12) A, (30,14) W, (28,16) A, (26,18) S, (24,20) A |

Counts per odd n (even k only): #S #A #W #I.

| n | S | A | W | I | n | S | A | W | I |
|---|---|---|---|---|---|---|---|---|---|
| 7 | 1 | 0 | 0 | 0 | 55 | 9 | 0 | 2 | 2 |
| 9 | 0 | 0 | 1 | 0 | 57 | 5 | 6 | 2 | 0 |
| 11 | 1 | 0 | 0 | 1 | 59 | 7 | 0 | 0 | 7 |
| 13 | 1 | 1 | 0 | 0 | 61 | 7 | 7 | 0 | 0 |
| 15 | 1 | 0 | 2 | 0 | 63 | 11 | 0 | 4 | 0 |
| 17 | 1 | 1 | 0 | 1 | 65 | 4 | 3 | 2 | 6 |
| 19 | 3 | 0 | 0 | 1 | 67 | 15 | 0 | 0 | 1 |
| 21 | 0 | 2 | 2 | 0 | 69 | 4 | 5 | 2 | 5 |
| 23 | 3 | 0 | 0 | 2 | 71 | 11 | 0 | 0 | 6 |
| 25 | 3 | 1 | 1 | 0 | 73 | 9 | 8 | 0 | 0 |
| 27 | 3 | 0 | 2 | 1 | 75 | 13 | 0 | 4 | 1 |
| 29 | 2 | 1 | 0 | 3 | 77 | 4 | 5 | 2 | 7 |
| 31 | 7 | 0 | 0 | 0 | 79 | 15 | 0 | 0 | 4 |
| 33 | 2 | 3 | 2 | 0 | 81 | 8 | 8 | 3 | 0 |
| 35 | 3 | 0 | 2 | 3 | 83 | 11 | 0 | 0 | 9 |
| 37 | 4 | 4 | 0 | 0 | 85 | 10 | 8 | 2 | 0 |
| 39 | 5 | 0 | 2 | 2 | 87 | 17 | 0 | 2 | 2 |
| 41 | 2 | 3 | 0 | 4 | 89 | 6 | 5 | 0 | 10 |
| 43 | 9 | 0 | 0 | 1 | 91 | 19 | 0 | 2 | 1 |
| 45 | 3 | 3 | 4 | 0 | 93 | 9 | 11 | 2 | 0 |
| 47 | 7 | 0 | 0 | 4 | 95 | 13 | 0 | 2 | 8 |
| 49 | 4 | 4 | 1 | 2 | 97 | 10 | 10 | 0 | 3 |
| 51 | 9 | 0 | 2 | 1 | 99 | 15 | 0 | 4 | 5 |
| 53 | 4 | 4 | 0 | 4 | 101 | 8 | 7 | 0 | 9 |

Totals (n ≤ 101): S 318, A 110, W 56, I 116 (plus 576 odd-k triples, all intransitive).

## 8. Literature and novelty

Known pieces (credit):
1. Off the line and n even (§1, B.2): 2610.30892 Thm 3.1; proved earlier in Lean in
   the-omega-institute/trureturing, D5/S0/CayleyGrowth/PrefixReversalTriple{Odd,Even}NonGeneration.lean
   (2026-09-13). Neither treats the odd line.
2. Necessity of gcd(m+1, k+1) = 1: 2610.30892, Lemma at TeX l. 277 (G = S_n ⇒ gcd(n−m, n−k) = 1);
   its special case n−m | n+1 is 2511.16959 Prop. 2.
3. The n-cycle: the proof of 2511.16959 Thm 1 (TeX l. 484) uses r_n r_{n−3} r_2, conjugate to our P
   when k = 2. The first-return description of P, the Beatty coordinates, the block analysis and the
   D≀D structure were not found anywhere.
4. k = 2 and k = 3 on the line: 2511.16959 Thms 1 and 3, including the blocks of 3 for (n, n−3, 2)
   when 3 | n.
5. A_n versus S_n: the parity observation of Bass–Sudborough (JPDC 63, 2003, Obs. 2.2(4)).

- Wording: v1 (2025-11-21) and v2 (2026-09-01) of 2511.16959 state the parity-split form; only v3
  (= DMTCS 28:1) says "if m + k < n". So the line counterexamples hit a wording slip in v3, not the
  authors' conjecture; v3's literal sentence is also contradicted by its own Theorem 1, which gives
  (n, n−3, 2), n ≡ 1 (mod 6). 2610.30892 calls the odd-n wording ambiguous and leaves the line open.
- Searched: arXiv search (prefix reversals / pancake / generating triples), the full math.GR and
  math.CO listings for 2026-09 and 2026-10, and the GXL evidence API (only 2511.16959 and 1803.01760
  came back, neither treats this line). No characterization of this line was found.
- Significance (my judgment): the counterexamples alone are minor (a wording slip in v3 only);
  the new content is the complete answer on the boundary line — orbits, the D≀D structure,
  primitivity, and the arithmetic criterion — via the n-cycle P = r_k r_n r_m, a first-return map
  of a rotation, which may help with the general open problem (m + k ≥ n).

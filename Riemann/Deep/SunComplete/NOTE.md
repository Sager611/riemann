# Sun's circular-permutation conjectures for primitive elements: complete proofs

Date 2026-10-22. Directory `Riemann/Deep/SunComplete/`. Check everything with `python3 verify.py` (about 10 s; `--cross` also recompiles and reruns the C cross-checks, about 1–2 min).

Numbering follows the printed paper: Z.-W. Sun, *Some new problems in additive combinatorics*, Nanjing Univ. J. Math. Biquarterly 36 (2019), 134–155 (= arXiv:1309.1679). The two statements are **Conjecture 3.5(i)** (joint with Q.-H. Hou) and **Conjecture 3.8**. In the arXiv v9 TeX their labels are Conj3.7 and Conj3.10.

Notation: q is a prime power, P = P_q is the set of primitive elements (generators of F_q^*), φ is Euler's function, ω(m) is the number of distinct prime factors of m, W(m) = 2^{ω(m)}, and θ(m) = φ(m)/m.

## 0. Results

**Theorem A (Sun–Hou, Conjecture 3.5(i)).** For every prime power q > 7, the elements of F_q can be numbered a_1, …, a_q so that all q sums a_1+a_2, …, a_{q−1}+a_q, a_q+a_1 are primitive. The bound is sharp: q = 3, 5, 7 fail. The statement also holds for q = 4, via the order (0, ω, 1, ω²).

**Theorem B (Sun, Conjecture 3.8, corrected).** Let q > 7 be a prime power and a₀ ∈ F_q with (q, a₀) ≠ (11, 2). Then there is a circular permutation (a_1, …, a_{q−1}) of F_q^* such that all q−1 elements a₀ + a_i a_{i+1} (indices mod q−1) are primitive. For (q, a₀) = (11, 2) no such permutation exists (Riemann/Deep/SunPrim/NOTE.md). The bound q > 7 is sharp: for q = 7, every a₀ ≠ 0 fails.

**Earlier results.**
- Theorem A was known for q > 18,888,871: Y.-F. She, arXiv:2609.19114, Thm 1.2(1), using Chvátal–Erdős.
- Theorem A was also known for primes p > 61. This follows from She's Remark together with Cohen–Oliveira e Silva–Trudgian (three consecutive primitive roots exist for every prime p > 61).
- Theorem B was known for q > 18,888,871, and for a₀ = 0 when q > 4: X.-Q. Luo and Y.-F. She, arXiv:2610.06289.

**Proof outline.**
1. Sections 2–3 give a sharper Chvátal–Erdős argument combined with Cohen's prime sieve. Proposition 4 is an explicit per-q criterion; Lemma 5 is a clean version of it.
2. Section 4: the criterion holds uniformly for all q ≥ 30,140 (Lemma 6). Checking each q individually (Proposition 7), it also holds for all q < 30,140 except 79 values (case A) and 83 values (case B). The largest exception is 6,091.
   - So the analytic method alone now covers every q > 6,091, down from 18,888,871.
3. Section 5: the exceptional fields are settled by explicit certificates.
   - She's construction.
   - Two new constructions for B: alternating products (Lemma 9) and geometric products (Lemma 10).
   - Explicit Hamilton cycles.
   - Exhaustive search for (11, 2).
4. Section 6: Luo–She's construction handles a₀ = 0.
5. Section 8: as an independent cross-check, both theorems are verified directly for every q < 30,140.

## 1. The graphs

- **H_q** has vertex set F_q, with x ~ y iff x ≠ y and x + y ∈ P. Theorem A is equivalent to: H_q is Hamiltonian.
- **G_{q,a₀}** has vertex set F_q^*, with x ~ y iff x ≠ y and a₀ + xy ∈ P. Theorem B is equivalent to: G_{q,a₀} is Hamiltonian.

Vertex degrees:
- In H_q, deg x = φ(q−1) − [2x ∈ P].
- In G_{q,a₀}, deg x = φ(q−1) − [a₀ ∈ P] − [a₀ + x² ∈ P]. This holds because, as y runs over F_q^*, a₀ + xy runs over F_q \ {a₀}.

Case parameters:

| case | V | n = \|V\| | β | ε | δ₀ (≤ minimum degree) |
|---|---|---|---|---|---|
| (A) | F_q | q | 0 | 1 if q odd, 0 if q even | φ(q−1) − ε |
| (B), a₀ ≠ 0 | F_q^* | q−1 | 1 | 1 | φ(q−1) − 2 |

In (A) with q even, ε = 0 because 2x = 0 is never primitive.

## 2. Character sums with Cohen's sieve

**e-free elements.** For e | q−1, an element z ∈ F_q is *e-free* if z ≠ 0 and z is not an ℓ-th power for any prime ℓ | e. So P is exactly the set of (q−1)-free elements. The indicator is

1_e(z) = θ(e) Σ_{d|e} μ(d)/φ(d) Σ_{ord χ = d} χ(z),

where every multiplicative character, the trivial one included, is extended by χ(0) = 0. There are exactly θ(e)(q−1) e-free elements.

**Sieve data.** Split the prime divisors of q−1 into a *core* (with product k) and *sieving primes* p_1, …, p_s, so that δ := 1 − Σ 1/p_i > 0. Put
- τ := θ(k)δ;
- Δ := (s−1)/δ + 2 if s ≥ 1, and Δ := 1 if s = 0 (then k = Rad(q−1) and δ = 1);
- **E := (W(k)Δ − 1)√q + 1.**

**Lemma 1 (minorant).** Let G := Σ_i 1_{kp_i} − (s−1)·1_k when s ≥ 1, and G := 1_k when s = 0. Then G(z) ≤ 1_P(z) for every z ∈ F_q.

*Proof.*
- If s = 0, then k = Rad(q−1) and G = 1_P.
- If z is not k-free (for instance z = 0), then G(z) = 0.
- Otherwise let J = {i : z is not p_i-free}. Then G(z) = (s − |J|) − (s − 1) = 1 − |J|, and z ∈ P exactly when J = ∅. ∎

**Lemma 2 (transforms).**
- (A) For every non-trivial additive character ψ, |Σ_{z∈F_q} G(z)ψ(z)| ≤ τE. Also Σ_z G(z) = τ(q−1).
- (B) Let a₀ ≠ 0 and g(z) := G(a₀ + z) for z ∈ F_q^*. For every non-trivial multiplicative character χ, |Σ_{z∈F_q^*} g(z)χ(z)| ≤ τE. Also Σ_{z∈F_q^*} g(z) ≥ τ(q−1) − 1.

*Proof.*

Step 1 (decomposition). For p ∤ k we have θ(kp) = θ(k)θ(p), and the terms with p ∤ ord η cancel, so

1_{kp} − θ(p)1_k = θ(k)θ(p) Σ_{d|k} μ(dp)/φ(dp) Σ_{ord η = dp} η.  (∗)

Since Σ_i θ(p_i) − (s−1) = δ, this gives G = δ·1_k + Σ_i (1_{kp_i} − θ(p_i)1_k).

Step 2 (single sums). Each term is a sum of characters η. In case (A), with ψ ≠ 1:
- Σ_z η(z)ψ(z) = −1 if η is trivial;
- otherwise it is a Gauss sum, of absolute value √q.

In case (B), with χ ≠ 1:
- Σ_{z≠0} η(a₀+z)χ(z) = −χ(−a₀) if η is trivial;
- otherwise it equals η(a₀)χ(−a₀)·J(χ, η), where J(χ, η) = Σ_w χ(w)η(1−w) is a Jacobi sum with |J| ∈ {1, √q}.

So each term has absolute value at most 1 when η is trivial, and at most √q otherwise.

Step 3 (adding up). The order-d characters carry weight 1/φ(d), and there are φ(d) of them. Hence:
- the transform of δ·1_k is at most δθ(k)(1 + (W(k)−1)√q);
- by (∗), each transform of 1_{kp_i} − θ(p_i)1_k is at most θ(k)θ(p_i)W(k)√q.

Using Σ θ(p_i) = s − 1 + δ, the total is

θ(k)[δ + δ(W(k)−1)√q + W(k)√q(s−1+δ)] = τ[W(k)Δ√q − √q + 1] = τE.

Step 4 (main terms).
- (A): Σ_z G(z) = (q−1)θ(k)[Σ θ(p_i) − (s−1)] = τ(q−1).
- (B): Σ_{z≠0} G(a₀+z) = τ(q−1) − G(a₀) ≥ τ(q−1) − 1, by Lemma 1. ∎

As a development check, a brute-force DFT over 16 fields (including q = 125, 343, 625, 729, 1024, 1331, 2401), with every partition and random a₀, gives max |Ĝ|/(τE) = 0.84.

**Lemma 3 (edge counts).** For X, Y ⊆ V, let e(X,Y) count the ordered pairs (x,y) ∈ X×Y with x + y ∈ P in case (A), or with a₀ + xy ∈ P in case (B); x = y is allowed. Put Ñ := q − 1 − β/τ. Then

**e(X,Y) ≥ (τ/n)·[ Ñ|X||Y| − E·√( |X|(n−|X|)·|Y|(n−|Y|) ) ].**

*Proof.*
- By Lemma 1, e(X,Y) ≥ Σ_{x∈X, y∈Y} G(x+y) in case (A), or Σ g(xy) in case (B).
- Expand by Fourier inversion on V (the additive group of F_q, or the cyclic group F_q^*, of order n): this sum equals (1/n) Σ_χ ĝ(χ) S_X(χ) S_Y(χ), where S_X(χ) = Σ_{x∈X} χ(x).
- The trivial character contributes ĝ(1)|X||Y|/n ≥ τÑ|X||Y|/n, by Lemma 2.
- The remaining characters contribute at most (τE/n)·Σ_{χ≠1}|S_X||S_Y| ≤ (τE/n)·√(Σ_{χ≠1}|S_X|²·Σ_{χ≠1}|S_Y|²), by Cauchy–Schwarz.
- By Parseval, Σ_{χ≠1}|S_X(χ)|² = |X|(n − |X|). ∎

## 3. The Hamiltonicity criterion

**Proposition 4.** Work in case (A), or case (B) with a₀ ≠ 0. Fix sieve data with Ñ > E, and put

α* := n(ε/τ + E)/(Ñ + E),  C := ⌊α*⌋ − 1,  a₁ := δ₀ − C + 1,  F(a) := Ñ²·a·(n − a − C) − E²·(a + C)(n − a).

Suppose a₁ ≥ 1, and either a₁ > n/2, or both F(a₁) > 0 and F(n/2) > 0. Then the graph is Hamiltonian.

*Proof.*

(i) α ≤ α*.
- Let I be independent, with s₀ = |I|. Pairs (x, y) ∈ I² with x ≠ y are non-edges. So e(I,I) only counts diagonal pairs x = y with 2x ∈ P (case A) or a₀ + x² ∈ P (case B), and e(I,I) ≤ ε s₀.
- Lemma 3 gives (τ/n)[Ñ s₀² − E s₀(n − s₀)] ≤ ε s₀, that is, s₀(Ñ + E) ≤ n(ε/τ + E).
- This replaces She's bound α ≤ 2W√q·q/(q−1) + 1, which split I into two halves.

(ii) κ ≥ C + 1.
- Suppose S ⊂ V with c := |S| ≤ C and Γ − S disconnected. (S = ∅ is allowed; that is the case where Γ itself is disconnected.) Let A be a smallest component of Γ − S, and put B := V ∖ (A ∪ S). Write a := |A| ≤ b := |B| = n − a − c.
- A vertex v ∈ A has all its neighbours in (A ∖ {v}) ∪ S. So a ≥ δ₀ − c + 1 ≥ a₁, and a ≤ (n−c)/2 ≤ n/2.
- There are no edges between A and B, so e(A,B) = 0. Lemma 3, with n − a = b + c and n − b = a + c, gives Φ(a,c) := Ñ²a(n−a−c) − E²(a+c)(n−a) ≤ 0.
- But ∂Φ/∂c = −Ñ²a − E²(n−a) < 0, so Φ(a,c) ≥ Φ(a,C) = F(a).
- F is a concave quadratic, since its leading coefficient E² − Ñ² is negative. So on [a₁, n/2] it is at least min(F(a₁), F(n/2)) > 0. If a₁ > n/2, no admissible a exists at all.
- Either way we get a contradiction. So κ ≥ C + 1 = ⌊α*⌋ ≥ α.
- This replaces She's bound κ ≥ φ(q−1) − W√q·q/(q−1). Keeping the factors (1−|A|/n)(1−|B|/n) shows that a small component cannot be large.

(iii) By the Chvátal–Erdős theorem, κ ≥ α together with n ≥ 3 implies the graph is Hamiltonian. ∎

**Lemma 5 (clean form (CS)).** Let D := δ₀ + 1, and let A' ≥ α*. One may take A' = ε/τ + E: this always works in case (A), where α* = q(ε/τ + E)/(q − 1 + E) and E ≥ 1, and works in case (B) whenever E ≥ 1/τ. Suppose D > A', 2A' < n and

**(CS)  (D − A')·Ñ²·(1 − 2A'/n) > E²·D.**

Then the graph is Hamiltonian.

*Proof.*
- As in (ii), take c ≤ C < A' and a ∈ [D − c, n/2].
- Since n − a − c ≥ (n−a)(1 − 2c/n), we get Φ(a,c) ≥ (n−a)·[a(Ñ²(1−2c/n) − E²) − E²c].
- (CS) implies Ñ²(1 − 2A'/n) > E², so the bracket increases with a.
- At a = D − c the bracket equals (D−c)Ñ²(1−2c/n) − E²D. Since c < A', this is at least (D−A')Ñ²(1−2A'/n) − E²D, which is positive by (CS). ∎

Asymptotically, Proposition 4 reads φ(q−1) ≳ qE/(q+E), with E ≈ W(k)Δ√q. She's condition was φ(q−1) ≥ 3W(q−1)√q·q/(q−1) + 1; Luo–She's was φ(q−1) ≥ 3W(q−1)√q + 2. The gain has three sources:
- the constant 3 drops to 1, from (i) and (ii);
- the Hoffman-type factor q/(q+E);
- the sieve replaces W(q−1) by W(k)Δ.

## 4. Analytic coverage

**Lemma 6 (uniform range).** For every prime power q ≥ 30,140, (CS) holds in case (A), and in case (B) for every a₀ ≠ 0. Hence H_q and G_{q,a₀} are Hamiltonian.

*Proof.*

Setup and worst-case bounds.
- Let r = ω(q−1), let ℓ_1 < … < ℓ_r be the prime divisors of q−1, and let p_j be the j-th prime, so ℓ_j ≥ p_j.
- Take the t smallest ℓ_j as the core, with t = t(r) from the table below. Then
  - W(k) = 2^t;
  - θ(k) ≥ α_t := ∏_{j≤t}(1 − 1/p_j);
  - δ ≥ δ_{t,r} := 1 − Σ_{t<j≤r} 1/p_j, so Δ ≤ Δ_{t,r} := (r−t−1)/δ_{t,r} + 2;
  - φ(q−1) ≥ α_r(q−1).
- In (CS), replace τ by τ_lb = α_t δ_{t,r}, E by E_r(q) = (2^tΔ_{t,r} − 1)√q + 1, A' by 1/τ_lb + E_r(q), D by α_r(q−1) − 1, and Ñ by q − 1 − β/τ_lb. Each replacement is legitimate:
  - Lemma 3 stays true with any larger E, and E_r(q) ≥ E.
  - 1/τ_lb + E_r(q) ≥ ε/τ + E ≥ α*, so it is an admissible A'.
  - The quantity LHS − RHS of (CS) increases with Ñ, decreases with E and A', and increases with D. Its D-coefficient is Ñ²(1−2A'/n) − E², which (CS) itself makes positive.
- In case (B) we also need E ≥ 1/τ. Since E ≥ √q, the condition √q ≥ 1/τ_lb suffices.
- The worst-case inequality obtained this way is called (CS_{r,t})(q).

Monotonicity in q. Divide (CS_{r,t})(q) by (q−1)²D to get

(1 − A'/D)·(Ñ/(q−1))²·(1 − 2A'/n) > (E/(q−1))².
- A'/D, A'/n and E/(q−1) all have the form f(q) = (a√q + b)/(cq − d) with a, b, c, d ≥ 0 and cq > d:
  - A'/D has a = 2^tΔ_{t,r} − 1, b = 1/τ_lb + 1, c = α_r, d = α_r + 1;
  - A'/n has c = 1 and d ∈ {0, 1};
  - E/(q−1) has b = 1 and c = d = 1.
- The numerator of f'(q) is (a/(2√q))(cq − d) − c(a√q + b) = −a(cq + d)/(2√q) − bc < 0, so all three are decreasing in q.
- Ñ/(q−1) is increasing in q.
- So once (CS_{r,t})(q₀) holds, it holds for every q ≥ q₀.
- In the exact check, √q₀ is replaced by a rational upper bound. That is allowed, because a larger √q only makes (CS) harder.

Table. Checked in exact rational arithmetic (`verify.py [1]`; regenerated by `gen_table.py`):

| r | t | q₀ (A) | q₀ (B) | P_r + 1 |
|---|---|---|---|---|
| 1 | 0 | 23 | 25 | 3 |
| 2 | 1 | 146 | 147 | 7 |
| 3 | 2 | 908 | 909 | 31 |
| 4 | 2 | 3,904 | 3,905 | 211 |
| 5 | 2 | 11,825 | 11,827 | 2,311 |
| 6 | 2 | 30,138 | 30,139 | 30,031 |
| 7–12 | 2 | 65,775 … 1,177,203 | 65,776 … 1,177,204 | ≥ 510,511 (always ≥ q₀) |
| 13–34 | 3 | P_r + 1 | P_r + 1 | |
| 35–60 | 4 | P_r + 1 | P_r + 1 | |

Here P_r = p_1⋯p_r ≤ q − 1. Every q with ω(q−1) = r ≤ 6 and q ≥ 30,140 is therefore covered, and for r ≥ 7 every q is covered.

Tail r ≥ 61. Use no sieve (t = r).
- 1/τ = 1/θ(q−1) ≤ ∏_{j≤r} p_j/(p_j − 1) ≤ ∏_{j≤r} (j+1)/j = r + 1, using p_j ≥ j + 1.
- q > P_r ≥ (r+1)!.
- E ≤ 2^r√q, A' ≤ (r+1) + 2^r√q, D ≥ (q−1)/(r+1) − 1, and Ñ ≥ q − 1 − (r+1).
- Put ρ_r := 2^r(r+1)²/√((r+1)!), so 2^r(r+1)²/√q < ρ_r. Using r + 1 ≤ 2^r√q and D ≥ q/(2(r+1)):
  - A'/D ≤ 4·2^r(r+1)/√q ≤ ρ_r;
  - 2A'/n ≤ 8·2^r/√q ≤ ρ_r;
  - E/(q−1) ≤ 2·2^r/√q ≤ ρ_r;
  - Ñ/(q−1) ≥ 1 − ρ_r.
- So (CS) holds as soon as (1 − ρ_r)⁴ > ρ_r², in particular whenever ρ_r < 1/3.
- ρ_{r+1}/ρ_r = 2((r+2)/(r+1))²/√(r+2) < 1 for r ≥ 6, and ρ_61 < 10⁻²⁰ (checked by `verify.py`). Also √q ≥ r + 1 ≥ 1/τ. ∎

**Proposition 7 (per q, below 30,140).** For every prime power 3 ≤ q < 30,140, Proposition 4 applies with a suitable partition (all 2^{ω(q−1)} partitions are tried, in exact rational arithmetic), except in the following cases:

- E_A = {3, 4, 5, 7, 9, 11, 13, 16, 19, 23, 25, 27, 29, 31, 37, 41, 43, 47, 49, 61, 67, 71, 73, 79, 103, 121, 127, 131, 139, 151, 157, 169, 181, 191, 199, 211, 223, 229, 241, 271, 277, 283, 307, 313, 331, 337, 361, 379, 421, 463, 541, 547, 571, 631, 661, 691, 841, 859, 967, 991, 1021, 1051, 1093, 1171, 1231, 1291, 1321, 1381, 1429, 1471, 1531, 1681, 1951, 2311, 2731, 3571, 4621, 4831, 6091}. This has 79 values; it is the failure set for case (A).
- E_B = E_A ∪ {8, 53, 1123, 1741}. This has 83 values; it is the failure set for case (B), uniformly in a₀ ≠ 0.

`verify.py [2]` re-derives both sets.

Two examples:

| q | q−1 | case | best partition | τ | E | α* | δ₀ | a₁ | result |
|---|---|---|---|---|---|---|---|---|---|
| 6091 | 2·3·5·7·29 | (A) | core {2,3} | 0.2076 | 1550.0 | 1239.6 | 1343 | 106 | F(a₁) < 0: fails |
| 9241 | 2³·3·5·7·11 (smallest possible θ for ω = 5) | (A) | core {2,3}, sieve {5,7,11} | 0.1887 | 2032.1 | 1670.3 | 1919 | 251 | F(a₁) > 0, F(n/2) > 0: passes |

By comparison, She's inequality fails for 1,305 prime powers below 30,140 (the largest is 30,091), and it holds only beyond 18,888,871.

**Corollary (analytic threshold).** By Lemma 6 and Proposition 7, with no computer search:
- H_q is Hamiltonian for every prime power q ∉ E_A;
- G_{q,a₀} is Hamiltonian for every q ∉ E_B and every a₀ ≠ 0.

In particular, both hold for every q > 6,091.

## 5. Certificates for the exceptional fields

**Lemma 8 (She, Remark after Thm 1.2).** Let p be an odd prime with x, x+1, x+2 ∈ P_p, and put m = (p−1)/2. Define
- a_{2k+1} = x/2 + 2k + 1 for 0 ≤ k ≤ m;
- a_{2k} = x/2 − 2k + 1 for 1 ≤ k ≤ m.

Then (a_1, …, a_p) is a Hamilton cycle of H_p.

*Proof.*
- The adjacent sums are x (for a_{2k+1} + a_{2k+2}), x + 2 (for a_{2k} + a_{2k+1}), and x + 1 (for the closing pair a_p + a_1).
- Entries in odd positions are distinct, and so are entries in even positions.
- An odd entry equals an even entry only if 2(i + j) = 0 with 1 ≤ i + j ≤ p − 1, which is impossible. ∎

**Lemma 9 (alternating products, new).** Let q be odd, m = (q−1)/2, and Q the subgroup of squares. Suppose:
- u ∉ Q (u ≠ 0);
- r generates Q;
- a₀ + u ∈ P and a₀ + ur ∈ P.

Put a_{2j+1} = r^j for 0 ≤ j < m, and a_{2j} = u r^{1−j} for 1 ≤ j ≤ m. Then (a_1, …, a_{q−1}) is a Hamilton cycle of G_{q,a₀}.

*Proof.*
- Products: a_{2j+1}a_{2j+2} = u, a_{2j}a_{2j+1} = ur, and the closing product is a_{2m}a_1 = ur^{1−m} = ur.
- The odd-position entries are exactly Q, and the even-position entries are exactly uQ = F_q^* ∖ Q.
- Each entry appears exactly once, because ord r = m. ∎

**Lemma 10 (geometric products, new).** Let q be even, n = q−1 = 2m + 1, and let s be primitive with a₀ + u, a₀ + us, a₀ + us² ∈ P. Put r = s² and t = s·u^{q/2}, so that t² = ur. Define
- a_{2j+1} = t r^j for 0 ≤ j ≤ m;
- a_{2j} = (u/t) r^{1−j} = t r^{−j} for 1 ≤ j ≤ m.

This is a Hamilton cycle of G_{q,a₀}.

*Proof.*
- Products: a_{2j+1}a_{2j+2} = u, a_{2j}a_{2j+1} = ur = us², and the closing product is a_{2m+1}a_1 = t²r^m = u s^{2m+2} = us.
- The entries are t·r^e for e ∈ {0, …, m} ∪ {−1, …, −m}, which is all residues mod n. Since r generates F_q^* (n is odd), every element appears exactly once. ∎

**Computation (`gen_certs.py`, about 1 min; files `certA.json`, `certB.json.gz`).** Fields are realised as F_p[x]/(f). `verify.py` rebuilds the power table of x by honest multiplication by x and requires x to have order q−1, which proves that f is irreducible. Primitivity is tested independently via z^{(q−1)/ℓ} ≠ 1.

- **Case (A)**, 75 fields of E_A with q > 7:
  - Lemma 8 certificate x for 62 primes;
  - explicit Hamilton cycles for q = 9, 13, 16, 25, 27, 29, 49, 61, 121, 169, 361, 841, 1681. These were found by exhaustive search for q ≤ 40 and by Pósa rotation–extension otherwise.
- **Case (B)**, 79 fields of E_B with q > 7, covering all 59,012 pairs (q, a₀) with a₀ ≠ 0:
  - Lemma 9 certificates (u, r) for 58,962 pairs;
  - Lemma 10 certificates (u, s) for 22 pairs (q = 8, 16);
  - explicit cycles for 27 pairs (q = 9: 4 values of a₀; q = 11: a₀ = 3, 4, 8; q = 13: 10; q = 25: 8; q = 31: a₀ = 3; q = 43: a₀ = 5);
  - for (q, a₀) = (11, 2), exhaustive search shows no Hamilton cycle exists.
- `verify.py [3]` checks every certificate condition, rebuilds and checks all (A)-cycles and 541 (B)-cycles (all of those with q < 30, plus a random sample), and redoes the (11, 2) search.

## 6. The case a₀ = 0 (Luo–She)

Write x = g^{X} for a generator g. Then xy ∈ P iff gcd(X + Y, n) = 1, where n = q − 1. So a Hamilton cycle of G_{q,0} is the same thing as a cyclic order of Z/n in which every adjacent sum is a unit.

- **n even.** Use (0, n−1, 2, n−3, …, n−2, 1). The adjacent sums are −1 and +1.
- **n = 7.** Use (6, 0, 5, 1, 3, 2, 4).
- **n = 2h + 1 with h ≥ 4.** Use (n−1, 0, n−2, 1, h, h−1, h+1, h−2, h+2, …, 2, n−3). The adjacent sums lie in {n−1, n−2, n−4, h+1}, which are units because n is odd and gcd(h+1, 2h+1) = 1.

The case n = 5 (q = 6) does not occur. This proves Theorem B for a₀ = 0 and every q > 4. `verify.py [4]` re-checks the sequences for every prime power 4 < q < 3000, and realises them inside F_q for 7 fields.

## 7. Proofs of the theorems

**Theorem A.**
- For q ≥ 30,140: Lemma 6.
- For 7 < q < 30,140 with q ∉ E_A: Propositions 4 and 7.
- For q ∈ E_A with q > 7: the certificates of §5.
- Small cases: q = 3, 5, 7 fail because H_q has a vertex of degree ≤ 1 (for q = 7, vertex 5 has only neighbour 0); q = 4 works with (0, ω, 1, ω²). Both are checked by exhaustive search.

**Theorem B.**
- For a₀ = 0: §6.
- For a₀ ≠ 0: the same three steps with E_B. In those steps the only pair without a Hamilton cycle is (11, 2).
- For q = 7, every a₀ ≠ 0 fails (exhaustive search), so q > 7 is sharp. ∎

## 8. Independent cross-check for all q < 30,140

This section does not rely on §§3–4. It uses `crossB.c` and `crossA.c` with honest F_p[x]/(f) arithmetic; together they take about 90 s and stay under 10 MB. Their outputs, with FNV-64 checksums, are `crossB.out` and `crossA.out`.

- **Case (B).** For all 3,331 prime powers 8 ≤ q < 30,140 and all a₀ ∈ F_q (46,731,497 pairs), a Lemma 9 or Lemma 10 certificate exists, except in two situations:
  - a₀ ≠ 0 with q ∈ {9, 11, 13, 25, 31, 43}: 28 pairs, all in E_B and handled in §5;
  - a₀ = 0 with q = 4^j, 2 ≤ j ≤ 7. Here Lemma 10 cannot apply: 3 | q − 1, so one of log u, log u + log s, log u + 2 log s is divisible by 3. These pairs are covered by §6.
- **Lemma 8.** Its certificate x exists for every prime 7 < p < 30,140 except 13, 29 and 61.
- **Case (A).** Pósa search finds, and fully re-checks, a Hamilton cycle of H_q for all 74 non-prime prime powers 8 ≤ q < 30,140.

So both theorems are verified directly for every q < 30,140. `verify.py [5]` re-derives the crossB checksums in Python for six sample fields, and `--cross` reruns both C programs and compares the outputs.

**Empirical remark.** For q odd and a₀ ≠ 0, Lemma 9 alone works for every q > 43 tested (q < 30,140). Proving that for all large q by character sums would give a construction-only proof of Theorem B, but it is not needed here.

## 9. Files

- `verify.py` — all checks listed above.
- `certA.json` — certificates for Theorem A on E_A.
- `certB.json.gz` — certificates for Theorem B on E_B, a₀ ≠ 0.
- `gen_certs.py` — regenerates the certificate files identically.
- `gen_table.py` — regenerates the Lemma 6 table.
- `crossA.c`, `crossB.c`, `crossA.out`, `crossB.out` — the §8 cross-check.

Total size is under 1 MB.

## References

- V. Chvátal, P. Erdős, A note on Hamiltonian circuits, Discrete Math. 2 (1972) 111–113.
- S. D. Cohen, T. Oliveira e Silva, T. Trudgian, On consecutive primitive elements in a finite field, Bull. LMS 47 (2015) 418–426.
- S. D. Cohen, S. Huczynska, The primitive normal basis theorem — without a computer, J. LMS 67 (2003) 41–56. This is the source of the prime sieve.
- Y.-F. She, Hamiltonicity in graphs defined by primes and primitive elements, arXiv:2609.19114 (v1, 16 Sep 2026).
- X.-Q. Luo, Y.-F. She, A conjecture on circular permutations over finite fields, arXiv:2610.06289 (v1, 5 Oct 2026).
- Z.-W. Sun, Some new problems in additive combinatorics, Nanjing Univ. J. Math. Biquarterly 36 (2019) 134–155; arXiv:1309.1679.

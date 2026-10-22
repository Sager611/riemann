# Audit of arXiv 2610.26686 (Bucić–He–Lee, "Proof of the Erdős–Szemerédi sunflower conjecture"), 2026-10-22

**Verdict: no false step.** Source: arXiv e-print v1 (19 Oct 2026, the only version; no comments field), file `sunflowers.tex`; all line numbers refer to it. Numbering is the compiled one (shared counter): Thm 1 thm:set, Thm 2 thm:vector, Lemma 3 lem:dim-argument, Lemma 4 lem:sunflower-dim, Obs 5, Lemma 6 lem:q_c-interpolation, Lemma 7 lem:degree-counting (the brief's "Lemma 3.7"), Prop 8 prop:transfer.
Claim (l. 254–260): every k-sunflower-free F ⊆ 2^[n] has |F| ≤ 3(2 − 1/(16k³))^n, for all k ≥ 3, n ≥ 1.

## 0. Where the k ≥ 4 barrier is avoided

- The barrier is a statement about testers **on the Boolean cube**: the naive coordinate tester gives Σ_{s ≤ (1−2/k)n} C(n,s) ≥ 2^{n−1} (R138), and Karingula–Lovett (arXiv 2609.27044) prove that for k ≥ 4 every tensor that is nonzero exactly on the cube's k-sunflower relation has full slice rank. The paper never builds a tester on the cube.
- Route: vector problem on [D]^n, D = 16k² (definition l. 272–274) → tester T_c = ∏_r q_c (l. 378–393) → interpolation Lemma 6 (l. 409–532) → common-space Lemma 3 (l. 321–354) → monomial count Lemma 7 (l. 540–595) → Chernoff, Thm 2 (l. 599–605) → random-restriction transfer, Prop 8 (l. 614–649) → constants, Thm 1 (l. 654–683).
- The gain is a power saving in the **alphabet size**: dim P = D^n·Pr[Bin(n,(k−1)/D) ≥ n/(k−1)], exponentially below D^n iff D > (k−1)². With D = 2 nothing is gained, consistent with the barrier. Prop 8 turns a saving C = ρD into (1 + ρ^{1/(D−1)})^N for sets by averaging over random block systems, not by any tester on sets.
- The step that does the work is the **choice of the rainbow values c** in Lemma 6. Exact check (§3a): with c ≡ 1 the reduced polynomial of q has monomials with d_1 = … = d_{k−1} = D−1 (no saving). With the interpolated c, every monomial has d_i ≤ i−1 ≤ k−2 for some i ∈ {1,…,k−1}.

## 1. Line-by-line

**Lemma 3 (l. 321–354), common space with a reference variable.** Suppose f = Σ_{i=1}^{k−1} Σ_j p_j(x_i) G_{ij}(x without x_i) and f is k-diagonal on V. Take any h: V → Q with Σ_V h p_j = 0 for all j. For every a ∈ V, S(a) = Σ_{x_1..x_{k−1}∈V} f(a,x_1,…,x_{k−1}) ∏ h(x_i) vanishes term by term: in term (i,j), the sum over x_i factors out as Σ h p_j = 0. Diagonality gives S(a) = f(a,…,a) h(a)^{k−1}. Hence h ≡ 0, P|_V = Q^V, and |V| ≤ dim P. **Valid.** No symmetry is used, and x_0 never needs a slice. The lemma is sharper than Tao's (d instead of (k−1)d); the standard slice-rank lemma would give the same exponent with a factor k−1 lost.

**Lemma 4 (l. 361–370), tuples with repeated vectors.** Suppose T ≠ 0 on V^k. Then every coordinate is monochromatic or rainbow. If the k vectors are pairwise distinct, they form a sunflower, which is excluded. If two are equal, every coordinate has a repeated entry, so no coordinate is rainbow. Every coordinate is then monochromatic, so all k vectors are equal. **Valid.**

**Observation 5 (l. 397–403).** q_c ≠ 0 exactly on monochromatic and rainbow tuples, and products over Q have no zero divisors. **Valid.**

**Lemma 6 (l. 409–532), the key step.** e(t) = [t = a_0], Δ_h = ∏_{i≥h} e(a_i), R_0 = 1, R_h = R_{h−1}(e(a_h) − P_h(a_h)), with P_h(t) = ∏_{i=1}^{h−1} (t−a_i)/(a_0−a_i) on a rainbow prefix and 0 otherwise.
- Identity: R_{h−1}Δ_h = R_{h−1}P_hΔ_{h+1} + R_hΔ_{h+1}. The first summand has degree ≤ h−1 in a_h and its coefficients do not involve a_h.
- Cases for R_h: prefix not rainbow gives R_{h−1} = 0. a_h = a_0 gives 1 − 1 = 0. a_h = a_i (1 ≤ i < h) gives 0 − 0 = 0. a_h new gives −P_h(a_h) ≠ 0.
- So q_c = Δ_1 − R_{k−1} is 1 on monochromatic tuples, 0 on mixed tuples and nonzero on rainbow tuples. **Valid.**
- Closed form (checked exactly, k = 3…6): on rainbow tuples c = (−1)^k ∏_{1≤i<h≤k−1} (a_h−a_i)/(a_0−a_i).
- k = 3: q = [a_2=a_0] − [a_1≠a_0](a_2−a_1)/(a_0−a_1).
- k = 4: q = e(a_2)e(a_3) + R_1P_2(a_2)e(a_3) + R_2P_3(a_3), with rainbow value c = (a_2−a_1)(a_3−a_1)(a_3−a_2)/((a_0−a_1)²(a_0−a_2)).

**Lemma 7 (l. 540–595).**
- Lagrange in the other variables keeps the exponent j ≤ k−2 ≤ D−1 of a_i, so every monomial of the reduced q_c has some d_i ≤ k−2. The reduced form is unique, so this is a property of q_c itself (§3a).
- For T_c the reduced polynomial is the product of the coordinate polynomials. Pigeonhole over k−1 slots gives a slot i with at least n/(k−1) low coordinates, so the x_i-monomial has at most ⌊(k−2)n/(k−1)⌋ coordinates with exponent ≥ k−1. Grouping gives the form required by Lemma 3.
- Count: Σ_{s≤(1−1/(k−1))n} C(n,s)(D−k+1)^s(k−1)^{n−s}. **Valid.**

**Thm 2 (l. 599–605).** Alon–Spencer A.1.12 with μ = n(k−1)/D and β = D/(k−1)²: Pr ≤ exp(μ(β−1−β ln β)) ≤ (e(k−1)²/D)^{n/(k−1)} ≤ (4k²/D)^{n/(k−1)}. The case D ≤ 4k² is trivial. **Valid.**

**Prop 8 (l. 614–649).**
- t = min(s, ⌊(N−s)/(D−1)⌋) gives (D−1)t + s ≤ N. G consists of D^t distinct s-sets. Blocks are disjoint and each set takes exactly one element per block, so a vector sunflower in G is a set sunflower (and conversely).
- By symmetry, Pr[S ∈ G] = D^t/C(N,s). Taking expectations gives |F_s| ≤ C(N,s)ρ^t.
- ρ^t ≤ ρ^s + ρ^{(N−s)/(D−1)−1}, and the binomial theorem finishes. **Valid.**
- Cosmetic slip: l. 648 credits "the second inequality" to the binomial theorem. The binomial theorem gives the equality; the last inequality uses ρ ≤ ρ^{1/(D−1)}.

**Thm 1 (l. 654–683).** D = 16k² gives ρ = 4^{−1/(k−1)} ∈ [1/2, 1) and 1 + 1/ρ ≤ 3. Also log ρ ≤ −1/k, and e^{−1/x} ≤ 1 − 1/(1+x) with 1 + k(D−1) = 16k³ − k + 1 ≤ 16k³. **Valid.**

## 2. Brief questions (1)

- **Diagonal on F^k, repeated sets, mixed blocks?** The tester is applied only to V^k, where V is the vector image of F_s ∩ G for **one** block system G. Repeated sets are repeated vectors, handled by Lemma 4; checked by brute force on all 4-tuples, repeats included (§3c). Tuples mixing blocks never occur: each realisation of G is bounded separately, and the averaging is linear.
- **Slice-rank lemma for an asymmetric tester with a reference variable?** The paper proves its own Lemma 3, which is correct for any f (see §1). The vanishing mechanism is confirmed numerically (§3d).
- **Interpolation on the alphabet {0,1}^b of a block?** No such alphabet is used. A block is a D-element set, and a transversal picks one of D letters [D] = {1,…,D} ⊂ Q. Lagrange needs D distinct nodes and k−2 ≤ D−1 (D ≥ k−1); both hold for D = 16k².

## 3. Computations (exact arithmetic over Q, except ranks mod 2^61−1, which only certify ranks from below)

(a) **q_c, n = 1** (qc_check.py). Checked for (k,D) = (3,3/5/7), (4,4/5/6/7/9), (5,5/6/7):
- the values (monochromatic 1, rainbow ≠ 0, mixed 0);
- the identity q_c = Σ_h R_{h−1}P_hΔ_{h+1};
- each summand has degree ≤ h−1 in a_h;
- every monomial of the reduced q_c has d_i ≤ i−1 for some i ≥ 1, and max over monomials of min_{i≥1} d_i is k−2.
- Control: the naive c ≡ 1 reaches min_{i≥1} d_i = D−1 (for example 8 at k = 4, D = 9).

(b) **T_c, n = 2** (tester_n2.py). Checked for (k,D) = (3,4) (4096 tuples), (4,4) (65,536 tuples) and (4,5) (390,625 tuples): the tester property holds on every tuple. The reduced polynomial in 2k variables was computed from the values and rebuilt at 300 random points. It has 52,900 = 230² and 254,016 = 504² monomials, every one of which has, in each coordinate, a slot i ≥ 1 with exponent ≤ 2, and its x_i-factor lies in P.

(c) **Dimension count vs actual rank, k = 4** (dimcount.py).

| D | n | dim P (formula) | rank of evaluation matrix on [D]^n | sunflower-free V | rank(P\|_V) |
|---|---|---|---|---|---|
| 4–6 | 1 | 3 | 3 | V_4(D,1) = 3 (tight) | 3 |
| 4 | 2 | 15 | 15 | V_4(4,2) = 9 (exhaustive; = 9 for all D ≥ 3 by König) | 9 |
| 5 | 2 | 21 | 21 | 9 | 9 |
| 6 | 2 | 27 | 27 | — | — |
| 4 | 3 | 63 | 63 | 27 (greedy) | 27 |
| 5 | 3 | 117 | 117 | 25 (greedy) | 25 |
| 6 | 3 | 189 | 189 | — | — |

T_c is diagonal on V^4 (all 4-tuples, repeats included) for every V in the table.

(d) **Contraction test** (the content of Lemma 3 without the monomial bookkeeping). Take V = [D]², which is not sunflower-free, and h ∈ P^⊥ (dim 1 for D = 4, dim 4 for D = 5; 3 random h each). Then S(a) = Σ T_c(a,x_1,x_2,x_3)h(x_1)h(x_2)h(x_3) = 0 for every a. Control: with the naive c ≡ 1, S(a) ≠ 0 for 48/48 and 75/75 pairs (a, h).

(e) **Transfer** (transfer_check.py).
- Set sunflower ⟺ vector sunflower on all k-subsets of G, for (k,D,t,|R|) = (4,4,2,1), (4,5,2,0), (4,4,3,2), (5,5,2,1), (3,4,3,1): 0 mismatches in either direction (32, 650, 912, 130 and 3648 sunflowers).
- Pr[S ∈ G] is uniform and equal to D^t/C(N,s), by enumerating all permutations for (N,D,s) = (7,2,3), (8,3,3), (7,3,2).
- ρ^t ≤ ρ^s + ρ^{(N−s)/(D−1)−1} and the sum identity hold for 5 values of ρ, D ≤ 20 and N ≤ 80.

(f) **Constants.**
- Exact rational checks for k = 3…400 (ln 4 > 1.386, e^{−x} ≤ 1 − x + x²/2): 0 failures.
- Chernoff step as an exact integer inequality S^{k−1} ≤ (4k²)^n D^{(k−2)n}, for k ∈ {3,4,5,6,8}, D ∈ {4k²+1, 8k², 16k²}, n ≤ 150: 0 failures. Per-coordinate slack is 0.45–0.88.
- k = 4: D = 256, C = 161.27, ρ = 0.62996, and 1 + ρ^{1/255} = 1.998189 ≤ 2 − 1/1024 = 1.999023.

## 4. Reproduce

Run `cd Riemann/Deep/Sunflower && PYTHONDONTWRITEBYTECODE=1 python3 <script>`. Runtimes are CPU time; peak memory is maximum resident set size.

| Script | CPU time | Peak memory |
|---|---|---|
| qc_check.py | 11 s | 39 MB |
| tester_n2.py | 83 s | 86 MB |
| dimcount.py | 5 s | 19 MB |
| transfer_check.py | 3 s | 16 MB |

The scripts read only the formulas above and do not need the paper.

## 5. Residual risk

- The general-n step is the tensor-product argument checked by hand in §1. The computations confirm it for n ≤ 2 (decomposition) and n ≤ 3 (rank).
- Only the v1 TeX was read. The claimed Lean 4 proof ("Palomar repository", l. 744) is not public yet.
- Karingula–Lovett's barrier concerns the cube relation and does not touch a tester on [D]^n. No vector analogue of it can hold: already at n = 1, q_c is a sum of k(k−1)/2 slices (Lemma 6), so its slice rank is ≤ k(k−1)/2 < D, not the full D.


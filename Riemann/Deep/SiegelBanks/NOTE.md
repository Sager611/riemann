# Banks, arXiv 2610.31402v1: audit notes (2026-10-23)

Sources read: the TeX e-print (BrauerSiegel-CURRENT.tex, 1846 lines). OpenAI, "Uniform exclusion of Landau–Siegel zeros",
build/paper.tex at github.com/openai/math, preprints/Uniform-exclusion-of-Landau-Siegel-zeros-October-1-2026.
Kerr, arXiv 2610.21099 TeX. Masser, Invent. Math. 63 (1981) 81–95, pp. 81–82 (GDZ PPN356556735_0063/LOG_0011).
Stark, Invent. Math. 23 (1974) 135–152, pp. 137–146 (GDZ PPN356556735_0023/LOG_0016). Scans are not kept (GDZ terms).

## 1. What the "OpenAI argument" is
- OAI Thm 1 (TeX l. 66): "There is an absolute constant c>0 such that every real zero β∈(0,1) of every primitive
  nonprincipal real Dirichlet L-function of conductor q≥3 satisfies (1−β) log q ≥ c."
- How OAI proves it: K = Q(a, b) with a² = d and b² = 2; points θ_n = n1 + n2 a + n3 b + n4 ab, n ∈ {0..N−1}^4; rows
  θ^{α1} σθ^{α2} στθ^{α3}, kept greedily in order of α1 + H(α2+α3); a nonzero interpolation determinant Δ
  (their uniform interpolation lemma); for inert p, θ^p ≡ σθ or στθ (mod p), which gives p^{E_p} | Δ with
  E_p = Σ⌊α1/p⌋; then the size of Δ is compared with its divisibility.
- Kerr (2610.21099) replaces Δ by a pigeonhole (Siegel's lemma) over polynomials of weighted degree ≤ D, and
  replaces the interpolation lemma by Masser's Thm 1. Banks says in his abstract "Kerr's form of the OpenAI argument".
  He generalizes Kerr in three ways: K/k Galois, one variable for each ι ∈ S ∋ 1, and the weights δ_p.
- Inputs actually used in the proofs: Masser 1981; Stark 1974 (Lemma 3, Lemma 4, Theorem 3); Minkowski, the Dirichlet unit
  theorem, Neukirch; Mertens for k (Rosen, Garcia–Lee); Rademacher's convexity bound. All are published. [OAI], [Ko],
  [Bar], [Gr] appear only in the introduction and in rem:almost and rem:compare. cor:siegel for k = Q re-proves OAI.

## 2. The zero estimate (Masser, Thm 1, p. 82, verbatim)
"Write μ = μ(Γ) ... let l ≥ 1 be the rank of Γ, and for 1 ≤ r ≤ n let l_r be the maximum rank of any subgroup Γ' of Γ for
which there exists a subspace V of C^n of dimension r such that v·γ = 0 for all v in V and all γ in Γ' ... μ(Γ) =
min_{1≤r≤n} (l − l_r)/r. Let γ_1, …, γ_l be basis elements of Γ, and for real N ≥ 0 denote by Γ(N) the set of points of
the form n_1γ_1 + … + n_lγ_l for non-negative integers n_1, …, n_l not exceeding N. Theorem 1. For real N ≥ 0, D ≥ 1 let
P be a non-zero polynomial of total degree at most D such that P(γ) = 0 for all γ in Γ(N). Then we have D ≥ (N/n)^μ."
Comparison with Banks's thm:masser (l. 1085). Banks takes the min over proper subspaces V of (l − rk(Γ∩V))/(κ − dim V),
which is Masser's μ with r = κ − dim V. The conclusions are the same.
- Basis: π is injective on Q^{n_K}, because the coordinate at ι = 1 recovers Σ v_t ω_t. So ω^S_1..ω^S_{n_K} is a Z-basis.
- ℘ = f1 − f2 has a nonzero coefficient in O_K ⊂ C; its total degree is ≤ its weighted degree ≤ D.
- Vanishing on {0..N−1}^{n_K} = Γ(N−1), so M = N − 1 ≥ 1.
- The bound depends only on (s, μ). This uniformity is essential, because Γ depends on F.
- prop:mu. Let U be the Q-span of the θ with θ^S ∈ V, so dim U = rk(Γ∩V) = r. The columns of u^S lie in V. The paper shows
  rk u^T ≥ t·r/n_K for every set T of t embeddings: (a) rk u^E = r; (b) rk u^{gT} = rk u^T, since u^{gT} is the entrywise
  image of u^T under g ∈ Gal(K̃/Q); (c) rk u^T ≥ |T∩B| for a set B of r independent rows. Averaging over g gives
  |G̃|·rk u^T ≥ Σ_g |gT∩B| = t r |G̃|/n_K. Hence r ≤ (n_K/s)·dim V, and therefore μ ≥ n_K/s; V = 0 gives equality. Correct.

## 3. The main line, step by step
- lem:frob: θ^{Np} ≡ Frob_P θ (mod P), then CRT. The fibres of P ↦ Frob_P have size |G|/(f|C_p|), so
  N(Q_{p,σ}) = Np^{|G|/|C_p|}, and Σ_{σ∈S∩C_p} log N(Q_{p,σ}) = |G| δ_p log Np. Checked in GP.
- lem:basis: one term of G_α has weighted degree α1 + H|α*| − (ν−H)j, which is ≤ D because ν = Np ≥ H. The case σ = 1 is
  handled separately. The basis change is unitriangular in the z1-degree. Correct.
- lem:vanish, lem:residues: c_α ≡ c'_α (mod Q^{m−b}) implies equal residues mod Q^m. The count is
  Π_α N(Q)^{m−b_α} = N(Q)^{mA−E_p}, CRT over the pairwise coprime Q_{p,σ}; it does not depend on the number of points. Correct.
- lem:pigeon: the number of f with ‖f‖_B ≤ T = P^{1/n_K}/(4A) is ≥ P^A (4 n_K A λ)^{−n_K A} B^{−n_K W}. The values satisfy
  |ρ(f(θ^S))| ≤ A‖f‖_B, using house(ιθ) = house(θ). A nonzero y ∈ 𝒫 has house(y)^{n_K} ≥ N(𝒫). The logarithmic condition
  (1/n_k) Σ δ_p E_p log Np > W log B + A log(4 n_K A λ_K) follows because |G| = n_K/n_k. Correct.
- lem:count: A = (s+1)I/D + O(D^{s−1}), W = (1+(s−1)/H) I + O(D^s), E_p = I/Np + O(D^s) with I = D^{s+1}/((s+1)! H^{s−1}).
  Checked by exact sums (check_counts.py §1–2).
- prop:main: E_p ≥ Σα1/Np − A. The error A·Σ_{Np≤D} log Np ≤ 2 n_k A D ≈ 2 n_k (s+1) I, which is O(1) after dividing
  by I. The terms (log B + log(Dλ))/D are O(1) because log λ_K ≤ N, D ≍ N^{n_K/s} and n_K/s > 1. The bound is
  Σ δ_p log Np/Np ≤ (1+(s−1)/H) n_k (log N + log λ_K) + O(1), and n_k log N = (s/|G|) log D + O(1). Correct.
- thm:precise: N = 2 + ⌈s(2X)^{s/n_K}⌉ gives D ≥ X. X ≥ (log d_K)^{n_K} gives N ≥ log d_K ≥ log λ_K. The case s = n_K
  forces k = Q and S = G, where the claim is trivial. Primes with Np < H give O(1).
- thm:main: H = ⌈4n_K/ε⌉, C0 = ⌈n_K/ε⌉ + 1, λ_K ≤ √d_K (Minkowski's second theorem). Correct, checked in 61,300 cases.
- lem:generator: (τ0) = 𝔞𝔟² with 𝔞 | 𝔡_{F/k}; replace 𝔟 by a small ideal in its class, then balance with squares of units.
  This gives |ιτ| ≤ C_k d_F^{1/n_k}. The basis ω√τ^a ξ^j gives λ_K ≤ C_k d_F^{1/(2n_k)}. GP check: check_generator.gp.
- thm:quadratic: S = {1} ∪ (G∖G_F), s/|G| = 1/2 + 1/|G| ≤ 1/2 + ε/8 because |G| ≥ (ℓ−1)/n_k ≥ 8/ε; (s−1)/H ≤ ε/16;
  δ_p = 1 on inert primes; at most n_k inert primes ramify in K (they lie above ℓ). The bound log d_K ≤ (ℓ−1) log d_F
  + 2n_k(ℓ−2) log ℓ holds. Hermite handles the small fields. Correct.
- lem:bias: Stark (8)–(10) with only the zero β kept, and [F:Q] ≪ log d_F; split primes ≤ (e/2)(−ζ_F'/ζ_F(s0)).
  The constant is absolute. cor:siegel: split ≥ (1 − 5/8 − 1/8) log X − O_k(1) and split ≤ (1/8) log X + c_1(1−β) log²X,
  hence (1−β) log d_F ≥ 1/(16 c_1 C). Small fields: 1−β ≥ L(1,χ)/max|L'|, effective.
- thm:bs: Stark Lemma 3 (p. 139: at most one zero in σ ≥ 1−(4 log|D_M|)^{−1}, |t| ≤ (4 log|D_M|)^{−1}, real and simple),
  Theorem 3 (p. 144: a simple real zero of ζ_B, B/A normal, is a zero of ζ_F with F = A or [F:A] = 2), and Lemma 4
  (p. 140: κ > c_4^{−1}(1−β_1), c_4 effective and absolute). The three-case argument and w_K ≥ 1,
  2^{r1}(2π)^{r2} ≤ (2π)^{n_K} are correct. rem:almost cites Lemma 10 (p. 146, 16 log|D_M|) correctly.

## 4. Effectivity (task 2)
In Siegel's theorem the ineffectivity enters through a case split on a hypothetical second exceptional character
(Deuring–Heilbronn, Landau, Siegel). There is no such split here. thm:precise bounds the log-mass of the inert primes
from above, unconditionally and with constants depending only on (n_K, H). A real zero is used only once, in
lem:bias, which needs only the Euler product and the Hadamard product of ζ_F. For k = Q this gives (1−β) log q ≥ c for
every real primitive χ (F = Q(√d), |d| = q) with c = min(1/(16 c_1 C), min over q < q0 of (1−β_q) log q). The bound q0
comes from the thresholds of thm:precise (ℓ = 67, H = 8448 at ε = 1/8; ℓ = 3 would also do). So c is effective but not
explicit. Explicit values for k = Q exist elsewhere: Banks 2610.15483 gives 3·10⁻³ for q ≥ 10^10, and 2610.21511 gives 5.8·10⁻⁴.

## 5. Numerical tests (all pass)
- check_counts.py: exact lattice sums (s ≤ 5, D ≤ 10^5), Kerr's constants, parameter inequalities in 364,252 cases,
  and the toy (pigeon) inequality in Q(√5, √−3) (lhs/rhs ≤ 0.30). Takes 15 s and 17 MB.
- check_mu.py: rank(u^T) ≥ t r/n_K for about 72k subspaces in 5 fields; the minimum tested ratio is ≥ n_K/s. Takes 9 s and 31 MB.
- check_frob.gp: 2820 congruence checks and the norm formula for S_3 over Q and C_3 over Q(√−3). Takes 64 s.
- check_generator.gp: τ-reduction on 750 quadratic F over Q(√10), Q(√−5) and a cubic field; max R = 0.40.

## 6. Minor points (none affects correctness)
1. Line 891: "ζ_k(1) ≠ 0" should read "(s−1)ζ_k(s) does not vanish at s = 1".
2. The hypothesis X ≥ (log d_K)^{n_K} in thm:precise is stronger than the proof needs (N ≥ log λ_K).
3. The bound "O(n_k D^{s+1}) = O(n_k I)" in prop:main is crude; the real error is about 2 n_k (s+1)·I.
4. rem:compare states a base-field version of Kominers' bound without proof. Kominers [Ko] and Barreto [Bar] are
   unpublished manuscripts. Neither is used in any proof.
5. Novelty is outside this audit. For k = Q the result is known (OAI, Kerr, Green). For k = Q(√−3) it also follows from
   OpenAI's 11/12 half-plane for Hecke L-functions over Q(√−3). It is new for other base fields k, as far as we checked.


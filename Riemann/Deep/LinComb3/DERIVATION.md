# DERIVATION — zeros on the critical line of linear combinations of N ≥ 3 Dirichlet L-functions (lincomb3)

Sources: Lamzouri, arXiv 2610.26841v2 (IsolatedZeros2.tex; TeX line numbers below), proof of Pro:BombieriHejhal l.868–964 and of Thm 1.4 l.1016–1027;
Bombieri–Hejhal, Duke Math. J. 80 (1995) 821–862 ("BH"; Theorem B, Lemmas 7–9, as cited by Lamzouri, not re-read here);
our isolated-zero bound for L(s, χ): Riemann/Deep/IsoZero (DERIVATION §2, §3.7, §8; cert.py hinge_p).
Every numeric constant below is printed by cert.py (logs/cert.log; §6 lists them), except the numerical-test values of §11 (numerics/, logs/analyse.log, logs/argcount.log).

## 0. Setting and statement

q ≥ 3; χ_1, …, χ_N (N ≥ 2) distinct primitive characters mod q with χ_j(−1) = (−1)^a for one a ∈ {0, 1}; b_1, …, b_N nonzero reals;
ε_j = τ(χ_j)/(i^a √q) the root number, θ_j = −(arg ε_j)/2, so e^{2iθ_j} = ε̄_j = i^a τ̄(χ_j)/√q (Lamzouri l.248, l.840–848);
F(s) = Σ_{j≤N} b_j e^{iθ_j} L_j(s), L_j(s) = L(s, χ_j).
H_j(T) = multiset (with multiplicity) of ordinates T < γ ≤ 2T of zeros of L_j with min(γ_{n+1} − γ_n, γ_n − γ_{n−1}) ≤ 2π/((log T)(log log T)^{1/4}) (l.862–865, unchanged).
δ_j = a liminf lower bound for B_j(T, c)/((T/2π) log T), B_j(T, c) the number of c-isolated zeros of L_j with T < γ ≤ 2T, for one fixed c > 0
(a zero is c-isolated if both neighbouring ordinates are more than 2πc/log T away; such a zero is simple and on the line).

Proposition N. For fixed ε > 0 and T ≥ T_0(ε, q, χ, b, N), the number of distinct zeros of F(½ + it) of odd multiplicity with t ∈ [T, 2T] is
 ≥ (1 − ε)(T/2π) log T − Σ_{j≤N} min(|H_j(T)|, (1/N + ε)(T/2π) log T).
Theorem N. liminf_{T→∞} N_0^{odd}(F; T)/N_F(T) ≥ 1 − Σ_j min(1 − δ_j, 1/N) ≥ Σ_j δ_j − (N − 1),
where N_0^{odd}(F; T) counts distinct zeros ½ + iγ of F of odd multiplicity with 0 < γ ≤ T and N_F(T) all nontrivial zeros with 0 < γ ≤ T, with multiplicity.
Corollary 1 (stage 1). δ_j = p(c) ≥ 0.672491 for every j (IsoZero §8: the constants of p(c) do not depend on χ), so:
N = 3: ≥ 3p − 2 ≥ 0.0174737; N = 2: ≥ 2p − 1 ≥ 0.3449825 (the earlier side result); N ≥ 4: nothing (1 − p > 1/4, so the bound is 0).
Corollary 2 (stage 2, the headline; §8–9). Only isolation from the other zeros on the line is needed (Theorem N^R), and for it the simple-zero certificate of arXiv 2610.08965 transfers
(Lemma M): N = 3: ≥ 3s₀ − 2 = 51087/2478195 = 0.0206146; N = 2: ≥ 2s₀ − 1 = 860123/2478195 = 0.3470764; N ≥ 4: 0; s₀ = 1669159/2478195 = 0.6735382.

## 1. Functional equation and the real function G (any N)

1.1 Λ(s, χ) = (q/π)^{(s+a)/2} Γ((s+a)/2) L(s, χ) = ε_χ Λ(1 − s, χ̄) (l.832–837). The gamma factor γ(s) = (q/π)^{(s+a)/2} Γ((s+a)/2) is the same for all j:
same modulus, all primitive (conductor q), same parity a. This is the only place where "same modulus and parity" enters, and it is needed for every N.
1.2 e^{iθ_j} ε_j = e^{−iθ_j} e^{2iθ_j} ε_j = e^{−iθ_j} ε̄_j ε_j = e^{−iθ_j} (|ε_j| = 1). With conj Λ(z, χ) = Λ(z̄, χ̄):
e^{iθ_j} Λ(s, χ_j) = e^{−iθ_j} Λ(1 − s, χ̄_j) = conj(e^{iθ_j} Λ(1 − s̄, χ_j)). Summing with real b_j: Λ_F(s) := γ(s) F(s) = conj Λ_F(1 − s̄) (l.851–855). ✓ for every N.
1.3 ϑ(t) := arg γ(½ + it) (continuous). G(t) := e^{iϑ(t)} F(½ + it) = Σ_j b_j Z_j(t), Z_j(t) := e^{i(ϑ(t) + θ_j)} L_j(½ + it).
Each Z_j is real (1.2 with one term); hence G is real, real-analytic, and its real zeros are the zeros of F on the line (l.946–955). A sign change of G on an interval gives a zero of F of odd multiplicity there.
1.4 i^a τ̄(χ)/√q = ε̄_χ: ε̄_χ = τ̄(χ)/((−i)^a √q) and 1/(−i)^a = i^a for a ∈ {0, 1}. So the brief's/Thm 1.4's phase is θ_j = −(arg ε_j)/2 (mod π; a sign change of e^{iθ_j} is absorbed in b_j).

## 2. Zero counting for F (any N)

F ≢ 0: a_F(n) = Σ_j b_j e^{iθ_j} χ_j(n), and distinct characters are linearly independent on (Z/qZ)^×, so some a_F(n) ≠ 0; the first nonzero term dominates for Re s large, so F ≠ 0 on Re s > A (l.849–850).
N_F(2T) − N_F(T) = (T/2π)(log T + B_F) + O(log T) (Dousselin, arXiv 2311.10285 = Forum Math. 37 (2025), eq. (eq:N), for any N; BH (3.12); Lamzouri l.856–858);
so N_F(T) ~ (T/2π) log T and N(2T, χ_j) − N(T, χ_j) = (T/2π) log T + O(T) (l.859–861).

## 3. Inputs for N functions (changes against Lamzouri's N = 2 marked "N:")

3.1 Joint CLT (Lem:SelbergCLT, l.873–876). N: the vector (log|L_j(½ + it)|/√(½ log log T))_{j≤N}, t uniform in [T, 2T], converges in law to (X_1, …, X_N) i.i.d. standard normal.
This is BH Theorem B (stated there for N functions; Lamzouri uses its N = 2 case). Hypotheses for each pair i ≠ k: Euler products, the common functional equation, Selberg's
density estimate N(σ, T, χ) ≪_χ T^{1−c(σ−½)} log T (l.878–880, fixed χ), and orthogonality Σ_{p≤x} χ_i(p) χ̄_k(p)/p = δ_{ik} log log x + O_q(1):
for i ≠ k the character χ_i χ̄_k mod q is non-principal (χ_i ≠ χ_k), so the sum converges (prime number theorem in progressions). Conjugate pairs χ_k = χ̄_i are allowed:
the covariance of Re P_i, Re P_k (P_j(t) = Σ_{p≤x} χ_j(p) p^{−½−it}) is ½ Re Σ χ_i(p) χ̄_k(p)/p + (terms with p^{∓2it}, mean o(1)) = O(1).
3.2 Separation (Lem:Separation, l.882–891). N: E := {t ∈ [T, 2T] : min_{i<k} |log|L_i(½+it)| − log|L_k(½+it)|| < δ√(log log T)}.
(log|L_i| − log|L_k|)/√(log log T) → (X_i − X_k)/√2 ~ N(0, 1), so meas(E) ≤ Σ_{i<k} (2δ/√(2π) + o(1)) T ≤ (N(N−1)/2)(0.7979δ + o(1)) T ≪_N δT.
3.3 BH Lemma 9 (Lem:BHLem9, l.898–910) is a statement about one L-function at a time ((Eq:StabilityZeros), (Eq:BH6.10) for each j).
N: apply it to each L_j with δ/N in place of δ and intersect: for more than (1 − δ)T(log T)/M blocks ℓ there is A_ℓ = ∩_j A_ℓ^{(j)} with meas(A_ℓ) > (1 − δ)(y_{ℓ+1} − y_ℓ)
on which (Eq:StabilityZeros) and (Eq:BH6.10) hold for every j ≤ N, with error terms O_N(√(log M)/δ²) (constants depend on N, which is fixed).
3.4 Dominance sets. S_j(T) := {t ∈ [T, 2T] : log|L_j(½+it)| > max_{i≠j} log|L_i(½+it)|}. N: meas S_j(T) = (1/N + o(1)) T.
The set {x : x_j > max_{i≠j} x_i} ⊂ ℝ^N is open with boundary of Gaussian measure 0, so by 3.1 (portmanteau) its frequency tends to P(X_j = max_i X_i) = 1/N (i.i.d., symmetric).

## 4. Proof of Proposition N (Lamzouri l.909–964 line by line)

4.1 Blocks y_ℓ = T + ℓM/log T (l.894–896); ℓ is good if (i) 3.3 holds for ℓ and (ii) (y_ℓ, y_ℓ + √δM/log T) ∩ A_ℓ ∩ E^c ≠ ∅ (l.911–915).
For ℓ with (i): meas((y_ℓ, y_ℓ + √δM/log T) ∩ A_ℓ) ≥ √δM/log T − δM/log T ≥ √δM/(2 log T) (δ ≤ 1/4; l.917–923).
If (ii) fails this set lies in E, so #{ℓ with (i), not (ii)} ≤ meas(E)·2 log T/(√δ M) ≪_N √δ T(log T)/M (l.924–927). Bad ℓ: ≪_N √δ T(log T)/M.
4.2 For good ℓ pick t_ℓ in the set of (ii). Markov's inequality on (Eq:BH6.10) for each j (l.929–935): B_ℓ := ∩_j B_ℓ^{(j)} ⊂ [0, M],
meas([0, M] ∖ B_ℓ) ≤ N·O(M√(log M)/(δ²(log log T)^{1/3})), and log|L_j(½ + it_ℓ + iu/log T)| = log|L_j(½ + it_ℓ)| + O((log log T)^{1/3}) for all u ∈ B_ℓ and all j ≤ N.
4.3 Dominance (l.940–949). j_ℓ := argmax_j log|L_j(½ + it_ℓ)| (unique: t_ℓ ∉ E). N: since t_ℓ ∉ E, log|L_{j_ℓ}(½+it_ℓ)| ≥ log|L_i(½+it_ℓ)| + δ√(log log T) for every i ≠ j_ℓ;
by 4.2, for u ∈ B_ℓ and every i ≠ j_ℓ: log|L_{j_ℓ}(½+it_ℓ+iu/log T)| ≥ log|L_i(½+it_ℓ+iu/log T)| + δ√(log log T) − O((log log T)^{1/3}) ≥ … + (δ/2)√(log log T) (T large).
Hence G(t_ℓ + u/log T) = b_{j_ℓ} Z_{j_ℓ}(t_ℓ + u/log T)(1 + η), |η| ≤ Σ_{i≠j_ℓ} |b_i/b_{j_ℓ}| e^{−(δ/2)√(log log T)} ≤ (N − 1)(max|b|/min|b|) e^{−(δ/2)√(log log T)} < 1,
so sign G = sign(b_{j_ℓ} Z_{j_ℓ}) on t_ℓ + B_ℓ/log T (L_{j_ℓ} ≠ 0 there: its log is finite).
4.4 Zeros (l.950–961). I_ℓ = [t_ℓ, t_ℓ + (1 − √δ)M/log T] ⊂ (y_ℓ, y_{ℓ+1}) (disjoint for distinct ℓ). k_ℓ := #zeros of L_{j_ℓ} (on and off the line, with multiplicity) with ordinate in I_ℓ
= (1 − √δ)M/2π + O(√(log M)/δ²) (3.3). Ordinates z_1 ≤ … ≤ z_{k_ℓ}; keep the interior τ-isolated ones (2 ≤ m ≤ k_ℓ − 1, both gaps > 2πτ/log T, τ = (log log T)^{−1/4});
each is a simple zero on the line. Each gap has u-length > 2πτ, larger than meas([0,M] ∖ B_ℓ) = O_N(M√(log M)/(δ²(log log T)^{1/3})) for T large, so it contains a point of B_ℓ;
choose one point per gap (the intervals between chosen points are disjoint). Z_{j_ℓ} is real, nonzero and of constant sign on each open gap (no ordinate of a zero of L_{j_ℓ} inside),
and changes sign at the simple zero z_m; by 4.3, G changes sign between the two chosen points: a zero of F of odd multiplicity, distinct for distinct m.
An interior zero that is not τ-isolated has a neighbouring gap ≤ 2πτ/log T = 2π/((log T)(log log T)^{1/4}) to its true neighbour (I_ℓ is an interval), so it lies in H_{j_ℓ}(T).
Count in I_ℓ ≥ k_ℓ − 2 − #(H_{j_ℓ}(T) ∩ I_ℓ).
4.5 Summation (l.962–964). Σ_{good ℓ} k_ℓ ≥ (1 − O_N(√δ))(T log T/M)((1 − √δ)M/2π − O(√(log M)/δ²)) and the −2 terms total 2T(log T)/M, so with δ = M^{−1/3}, M = M(ε, N) large:
count ≥ (1 − ε/2)(T/2π) log T − Σ_{good ℓ} #(H_{j_ℓ}(T) ∩ I_ℓ) = (1 − ε/2)(T/2π) log T − Σ_j Σ_{good ℓ : j_ℓ = j} #(H_j(T) ∩ I_ℓ).
4.6 Capacity (N: new, harmless for N = 2). For each j, Σ_{good ℓ: j_ℓ=j} #(H_j ∩ I_ℓ) ≤ |H_j(T)| and ≤ Σ_{good ℓ: j_ℓ=j} k_ℓ.
By 4.3, t_ℓ + (B_ℓ ∩ [0, (1 − √δ)M])/log T ⊂ I_ℓ ∩ S_j for good ℓ with j_ℓ = j; these sets are disjoint and each has measure ≥ ((1 − √δ)M − o(1))/log T.
So #{good ℓ : j_ℓ = j} ≤ meas S_j(T)·log T/((1 − √δ)M − o(1)) = (1/N + o(1)) T(log T)/((1 − √δ)M) (3.4), and, with k_ℓ ≤ (1 − √δ)M/2π + C_N√(log M)/δ² (3.3),
Σ_{good ℓ: j_ℓ=j} k_ℓ ≤ (1/N + o(1)) (T/2π) log T · (1 + 2πC_N√(log M)/(δ²(1 − √δ)M)): the factors 1 − √δ cancel. With δ = M^{−1/3} the bracket is 1 + O_N(√(log M) M^{−1/3}) ≤ 1 + ε/2
for M ≥ M_0(ε, N); o(1) is as T → ∞ for fixed M. Order of choices: ε, then M = M(ε, N) and δ = M^{−1/3}, then T ≥ T_0(ε, M). This proves Proposition N.

## 5. Proof of Theorem N

5.1 |H_j(T)| ≤ (1 − δ_j + o(1))(T/2π) log T: as (log log T)^{−1/4} < c for T large, a zero in H_j(T) is not c-isolated, except possibly the first and last zero of the block
(isolation measured inside the block multiset); with 2 (Riemann–von Mangoldt) this is l.1016–1027 / IsoZero DERIVATION 8.1, per j.
5.2 Proposition N: count in [T, 2T] ≥ (1 − Σ_j min(1 − δ_j, 1/N) − 2ε)(T/2π) log T for T ≥ T_0. The same T_0 serves all j (finitely many).
5.3 Dyadic summation: Σ_k over [T/2^{k+1}, T/2^k] with T/2^{k+1} ≥ T_0 gives ≥ (1 − Σ_j min(1 − δ_j, 1/N) − 3ε)(T/2π) log T; distinct blocks give distinct zeros
(a zero at a common endpoint is counted at most twice, O(log T) total); divide by N_F(T) ~ (T/2π) log T (2) and let ε → 0. ∎
5.4 Multiplicity bookkeeping: H_j, k_ℓ, N(T, χ_j), N_F count with multiplicity; c- and τ-isolated zeros are simple; zeros of F are counted once each (odd multiplicity); so the numerator is a
lower bound for the number of zeros of F on the line counted with multiplicity as well.
5.5 Order of quantifiers: q, χ_j, b_j, N, c fixed; δ_j is a liminf statement for that c; ε → 0 last; M, δ depend on ε and N only; T → ∞.
No step uses GRH, RH or the well-spacing hypothesis: BH Theorem B and Lemmas 7, 9 are used only in the form Lamzouri states unconditionally for Dirichlet L-functions
(via Selberg's density theorem), and δ_j comes from unconditional pair correlation (Garunkštis–Paliulionytė; Lamzouri 2609.02882) and Wu's mollified mean value (IsoZero §8.2).

## 6. Constants (cert.py, logs/cert.log)

6.1 p(c) = (2 − λ − C_MT − α_ρ B(ρ))/(1 − λ), the isolated-zero lower bound of IsoZero DERIVATION §3.7 (hinge bound), computed by IsoZero/cert.py hinge_p(c)
(λ = λ_C(c) cluster loss, ρ, α_ρ, B(ρ) as there; C_MT = ½ + cot(1/√2)/√2 = 1.3274992963…). The constants do not depend on χ (IsoZero §8.4(ii)), so δ_j = p(c) for every j.
6.2 Values: see logs/cert.log (p(10⁻³), p(10⁻⁴), p(10⁻⁵); the largest is used; all ≥ 0.672487).
6.3 N-table: threshold (N − 1)/N; naive Np − (N − 1); capacity 1 − N·min(1 − p, 1/N); method ceiling N(2 − C_MT) − (N − 1) (no isolated-zero bound from this Hilbert-space method exceeds 2 − C_MT).
6.4 Headline N = 3: 3p − 2 ≥ 0.0174737 (asserted in cert.py). Margin p − 2/3 ≥ 0.0058 (the N = 3 case needs p > 2/3). Gap to the method's N = 3 ceiling 4 − 3C_MT = 0.0175021: < 3·10⁻⁵.
6.5 With Lamzouri's own inputs (Σ(β − ½) ≤ 0.0806T, α_0): p_A = 2 − C_MT − 0.0806α_0 = 0.6552029…, 3p_A − 2 < −0.034: the N = 3 case is out of reach of 2610.26841's constants and
follows only from the hinge bound of IsoZero.
6.6 Quoted (not certified) comparison: Dousselin (arXiv 2311.10285, Forum Math. 37 (2025), Thm "application"): κ_F ≥ 4.91·10⁻⁹ (N = 3), 3.68·10⁻⁹ (N = 4), 2.94·10⁻⁹ (N = 5),
1.46·10⁻⁹ (N = 10), 2.161·10⁻⁶/(N log N)(1 − o(1)) (N large).

## 7. What the counting step can and cannot give (barrier for the method) and the missing lemma

7.1 What the proof uses about the zeros of L_j: (a) |H_j(T)| ≤ (1 − δ_j + o(1))N_T, N_T := (T/2π) log T (global, per j); (b) per-block zero counts k_ℓ ≈ M/2π (3.3);
(c) the dominance blocks of L_j carry (1/N + o(1)) of the blocks (3.4). Nothing relates the position of H_j(T) to the dominance blocks of L_j.
7.2 Sharpness. Consider any arrangement consistent with (a)–(c) in which, for each j, the non-isolated zeros of L_j fill the blocks where L_j dominates, as far as their capacity (1/N)N_T allows.
For it, the per-block count of 4.4 is k_ℓ − #(H_{j_ℓ} ∩ I_ℓ) exactly (no other sign changes are forced), and the total is (1 − Σ_j min(1 − δ_j, 1/N))N_T: Proposition N is attained.
So no bookkeeping of the dominance mechanism that uses only (a)–(c) gives more than Theorem N. This is a statement about the method and its inputs, not about the zeros:
under GRH plus well-spacing BH obtain 100%, and nothing known excludes independence (7.4), which would give ≥ p for every N.
7.3 Consequences. With δ_j ≤ 2 − C_MT = 0.672501 (ceiling of the Hilbert-space isolated-zero method) Theorem N is positive only for N ≤ 3. For N ≥ 4 positivity needs some δ_j > 1 − 1/N ≥ 3/4:
an isolated (hence simple and on-line) proportion above 3/4 for a Dirichlet L-function, beyond every unconditional simple-zero result (best: 67.35%, arXiv 2610.08965, for ζ and, by transfer, fixed χ).
7.4 Missing lemma (non-concentration on the dominance set). Let Ŝ_j(T) = ∪{I_ℓ : ℓ good, j_ℓ = j} (the blocks in which L_j dominates; equivalently, up to o(T) in measure,
{t : P_j(t) > max_{i≠j} P_i(t)} with Selberg's polynomials P_j(t) = Re Σ_{p≤X} χ_j(p) p^{−½−it}, X = T^{1/(log log T)²}). Needed: for each j,
  limsup_T |H_j(T) ∩ Ŝ_j(T)| / |H_j(T)| ≤ θ_j.
Then 4.5 gives count ≥ (1 − Σ_j min(θ_j(1 − δ_j), 1/N) − ε)N_T (with H^R_j and δ^R_j of §8 in place of H_j, δ_j, the same).
Positivity for N functions with δ_j = s₀ (§9) needs θ < 1/(N(1 − s₀)): θ < 0.7658 (N = 4), 0.6126 (N = 5), 0.3063 (N = 10); with δ_j = p: 0.7633, 0.6107, 0.3053.
Asymptotic independence (θ_j = 1/N, the isolated zeros of L_j equidistributed with respect to which L-function dominates) gives proportion ≥ (1/N)Σ_j δ^R_j ≥ s₀ = 0.6735 for every N ≥ 2
(≥ p = 0.6725 with stage-1 inputs), also improving N = 2 (0.3471 → 0.6735).
7.5 A route to 7.4: Lamzouri's Hilbert-space inequality holds for every finite conjugation-invariant multiset (IsoZero §8.4(i)), in particular for the zeros of L_j with ordinates in Ŝ_j(T)
(a union of blocks of length ≍ M/log T, M → ∞: pairs across block ends cost O(1/M)); the off-line hinge input is global and O(10⁻⁵), so it restricts trivially. What is missing is the
energy of that sub-multiset: a weighted pair-correlation theorem Σ_{γ,γ'} W(γ) K((γ − γ') log T/2π)² ≤ (C_MT + o(1)) Σ_γ W(γ) for W = a smooth function of (P_1(γ), …, P_N(γ))
(short Dirichlet polynomials at the zeros), i.e. Montgomery/BGST/Garunkštis–Paliulionytė pair correlation twisted by functions of Selberg's polynomials. We found no such statement
(arXiv, Zenodo, HAL, GXL searches, logs/novelty.log). With it, δ'_j = p(c) on each Ŝ_j, θ_j = 1/N and 7.4 follows.


## 8. Isolation from on-line zeros suffices (Theorem N^R)

8.1 In 4.4 only the sign pattern of Z_{j_ℓ} on the line is used, and Z_j changes sign exactly at the zeros of L_j on the line of odd multiplicity; zeros off the line do not affect it.
So the count works with a weaker notion. For κ > 0 let Y*_j(T, κ) be the set of simple zeros ½ + iγ of L_j on the line, T < γ ≤ 2T, with no other zero of L_j on the line
(any multiplicity) with ordinate within 2πκ/log T; δ^R_j(κ) := liminf_T |Y*_j(T, κ)| / ((T/2π) log T). Then δ^R_j(κ) ≥ δ_j(κ) (c-isolation with c = κ is isolation from all zeros).
8.2 Theorem N^R. If δ^R_j ≤ δ^R_j(κ) for one fixed κ > 0, then liminf N_0^{odd}(F; T)/N_F(T) ≥ 1 − Σ_j min(1 − δ^R_j, 1/N) ≥ Σ_j δ^R_j − (N − 1).
Proof: 4.1–4.3 unchanged. In 4.4 take, for each y ∈ Y*_{j_ℓ}(T, κ) whose κ/2-neighbourhood lies in I_ℓ, the two open intervals (y − πκ/log T, y) and (y, y + πκ/log T):
they contain no zero of L_{j_ℓ} on the line, so Z_{j_ℓ} has a constant sign on each, opposite signs on the two (y simple); each has u-length πκ, which exceeds the exceptional measure
O_N(M√(log M)/(δ²(log log T)^{1/3})) for T large, so it contains a point of B_ℓ; by 4.3 G changes sign in (y − πκ/log T, y + πκ/log T). Distinct points of Y* are more than 2πκ/log T apart,
so these intervals are disjoint. At most one point of Y* lies within πκ/log T of each end of I_ℓ. Count in I_ℓ ≥ |Y*_{j_ℓ} ∩ I_ℓ| − 2 = k_ℓ − #(H^R_{j_ℓ} ∩ I_ℓ) − 2 with
H^R_j(T) := (zeros of L_j with T < γ ≤ 2T, with multiplicity) ∖ Y*_j(T, κ), |H^R_j(T)| ≤ (1 − δ^R_j + o(1))(T/2π) log T. Then 4.5, 4.6 and 5.2–5.5 verbatim. ∎

## 9. Lemma M: simple-zero certificates give zeros isolated from the other on-line zeros

Setting (arXiv 2610.08965, TeX knausgard839.tex; "the paper" in this section): window f on I = [−½, ½], smoothing f_ε = η_ε² with η_ε ∈ C_c^∞(I) (eq:smooth, l.761), K_ε = f̂_ε (a Schwartz function),
F_z(u) = η_ε(u) e^{−2πiuz}, ⟨F_z, F_w⟩ = K_ε(w − z̄), A = Σ_z ν_z F_z⟨F_{z̄}, ·⟩ over a finite conjugation-invariant multiset Z (l.763–778): A is Hermitian (l.773), tr A = Σ ν_z = |Z|.
For a real point y, P_y := F_y⟨F_y, ·⟩ ⪰ 0 and tr(P_u P_v) = K_ε(u − v)² =: g_ε(u − v).
The paper's proof of its Thm thm:line-general (l.1092–1104) uses about Z only: (eq:key2) (l.1081, pure linear algebra: Lemma lem:threshold at c = 2 and the splitting P/Q, valid for every
conjugation-invariant multiset), the certificate (eq:certificate2) (l.1088, for every finite set of real points; proved for (x, b) = (8722, 4389)·10⁻⁶ by Props prop:gram2, prop:local2, l.1149–1235),
an upper bound tr A² ≤ (C(f_ε) + o(1))|Z| (prop:energy, l.780) and spn(Y) ≤ |Z| + o(|Z|).
9.1 Merging. Fix κ ∈ (0, 1/16] and κ' = √κ. Let Y be the distinct real points of Z with multiplicities m_a. Join two points of Y when their distance is < κ; let the components C be the maximal chains.
Replace every component with at least two points and diameter ≤ κ' by the single real point ȳ_C = Σ_{a∈C} m_a y_a / m_C of multiplicity m_C = Σ_{a∈C} m_a; leave everything else. Call the result Z'.
(i) Z' is conjugation-invariant (only real points move), |Z'| = |Z| = N, tr A' = N.
(ii) The simple real points of Z' are Y* ⊔ L: Y* = the simple real points of Z with no other real point within distance < κ (singleton components of multiplicity 1),
L = the simple points of components of diameter > κ'. Merged points have multiplicity ≥ 2 and enter the paper's Q (rank one, trace m_C), so (eq:key2) holds for Z'.
A merged point ȳ_C lies in [min C, max C], and every real point outside C is at distance ≥ κ from C, so ȳ_C coincides with no other point of Z' (the counting in (ii) is exact).
(iii) Energy. A' − A = Σ_C D_C, D_C = Σ_{a∈C} m_a (P_{ȳ_C} − P_{y_a}), so ‖A' − A‖²_HS = Σ_{C,C'} Σ_{a∈C, b∈C'} m_a m_b M(s_a, t_b; Δ) with s_a = y_a − ȳ_C, t_b = y_b − ȳ_{C'},
Δ = ȳ_C − ȳ_{C'} and M(s, t; Δ) = g(Δ) − g(Δ − t) − g(Δ + s) + g(Δ + s − t) = −∫_0^s∫_0^t g''(Δ + σ − τ) dτ dσ.
Expanding g'' to first order with remainder, the terms g''(Δ)st and g'''(Δ)(s²t − st²)/2 vanish after summation because Σ_a m_a s_a = Σ_b m_b t_b = 0 (centroids);
the remainder gives |Σ_{a,b} m_a m_b M| ≤ 2 ψ_ε(Δ) m_C m_{C'} κ'⁴, ψ_ε(u) := sup_{|v−u|≤2κ'} |g_ε⁗(v)| (|s_a|, |t_b| ≤ κ'). Hence
 ‖A' − A‖²_HS ≤ 2κ'⁴ Σ_{a,b∈Y} m_a m_b ψ̃_ε(y_a − y_b),  ψ̃_ε(u) := sup_{|v−u|≤4κ'} |g_ε⁗(v)|.
(iv) Pair counts (input (P)): for fixed χ and R ≥ 1, the number of ordered pairs of zeros (with multiplicity, on or off the line) in the block with normalised ordinate distance ≤ R is
≤ C_χ (1 + R) N for T ≥ T_0(R). (From Selberg's mean value ∫_T^{2T} (S(t+h, χ) − S(t, χ))² dt ≪_χ T log(2 + h log T), 0 < h ≤ 1 [Selberg 1946; Fujii 1974]: split the block into
intervals of normalised length R; the pairs lie in the same or adjacent intervals, and Σ (counts)² ≪ (N/R)(R² + log(2 + R)) for a suitable offset of the partition
(average the squared sliding-window count over the offset; uniform for 1 ≤ R ≤ log T since Selberg's bound is uniform in 0 < h ≤ 1). Reviewer: sound.) Since g_ε⁗ is a Schwartz function,
Σ_{a,b} m_a m_b ψ̃_ε(y_a − y_b) ≤ Σ_{k≥0} sup_{[k,k+1)} ψ̃_ε · C_χ(k + 2) N + o(N) =: C_χ Ψ_ε N + o(N), Ψ_ε < ∞ (uniform in κ ≤ 1/16); the terms k ≥ log T use only the
classical N(t + 1, χ) − N(t, χ) ≪ log t (at most O(N log T) pairs per unit of distance) against the Schwartz decay of ψ̃_ε, so they are o(N). So ‖A' − A‖_HS ≤ (2C_χΨ_ε)^{1/2} κ √N.
(v) tr A'² = ‖A'‖²_HS ≤ (‖A‖_HS + ‖A' − A‖_HS)² ≤ (√(C(f_ε) + o(1)) + c_ε κ)² N, c_ε := (2C_χΨ_ε)^{1/2} (A, A' Hermitian).
(vi) Long chains: a point of a component of diameter > κ' has at least κ'/(2κ) − 1 ≥ 1/(4√κ) other real points within distance κ'/2 ≤ ½ (consecutive gaps < κ), so by (iv)
|L| ≤ 4√κ · #{ordered pairs at distance ≤ 1} ≤ 4√κ · 2C_χ N = 8C_χ √κ N.
9.2 Conclusion. The paper's 3-line proof of thm:line-general applied to Z' (with spn(Y*⊔L) ≤ spn(Y) ≤ N + o(N); b ≥ 0 by the definition of a certificate, l.1084; ω ≥ 0 without loss) gives
 (1 − x)|Y* ⊔ L| ≥ (2 − (√C(f_ε) + c_εκ)² − b − ω(ε)) N − o(N),
so liminf_T |Y*|/N ≥ (2 − (√C(f_ε) + c_εκ)² − b − ω(ε))/(1 − x) − 8C_χ√κ. Given η > 0 choose ε with (2 − C(f_ε) − b − ω(ε))/(1 − x) ≥ s₀ − η/2 (C(f_ε) → C(f) ≤ C'_0, ω(ε) → 0,
l.761–790, 1088), then κ small. Lemma M: for every η > 0 there is κ > 0 with δ^R(κ) ≥ s₀ − η, s₀ = (2 − C'_0 − b)/(1 − x) = 1669159/2478195 (paper, Thm thm:line, l.172–178; proof l.1229–1235; C'_0 = 6639737/5000000, l.1199).
9.3 Fixed χ (transfer of the paper to L(s, χ)): the only arithmetic input of the paper is prop:energy. Note that it is the energy with complex parameters, tr A² = Σ_{ρ,ρ'} K_ε(z_ρ − z_ρ')², z_ρ − z_ρ' = i(ρ − ρ') log T/2π,
which is not an ordinate-only statement: off-line pairs enter through K_ε at complex arguments. Lem:GaPa is exactly of this form (its left side is Σ r̂(i(ρ − ρ') log T/2π) w(ρ − ρ') over all zeros, on and off the line).
For the dyadic block multiset of zeros of a fixed L(s, χ) it holds with the same C(f_ε)
by Lamzouri's dyadic Lem:GaPa (2610.26841 l.968–975: test functions r even, real, supported in [−1,1], Lipschitz at 0, weight w(u) = 4/(4−u²), all zeros with multiplicity, main term
(T/2π) log(qT)(r(0) + 2∫_0^1 α r(α) dα)) applied to r = R_ε = f_ε * f_ε and R_ε'' exactly as in the proof of prop:energy (l.791–805, the weight removed by R_ε − R_ε''/(4L²)); R_ε is even, real, smooth and supported in [−1, 1]
because η_ε is real and even (l.759), so both are admissible test functions.
Conjugation invariance of the multiset is the functional equation (IsoZero DERIVATION 8.4(i)); spn ≤ (T/2π) log T + O(1) = N + o(N); (P) as in 9.1(iv).
9.4 Corollary. Quantifiers: for every η > 0 Lemma M gives κ_j(η) > 0 with δ^R_j(κ_j) ≥ s₀ − η; Theorem N^R (with κ = min_j κ_j, as δ^R_j(κ) decreases in κ) gives
liminf N_0^{odd}(F;T)/N_F(T) ≥ 1 − Σ_j min(1 − s₀ + η, 1/N); the left side does not depend on η, so letting η → 0: liminf N_0^{odd}(F;T)/N_F(T) ≥ 1 − Σ_j min(1 − s₀, 1/N):
 N = 2: 2s₀ − 1 = 860123/2478195 = 0.3470764…;  N = 3: 3s₀ − 2 = 51087/2478195 = 0.0206146…;  N ≥ 4: 0 (1 − s₀ = 0.3265 > 1/4).
For ζ the same lemma gives: at least s₀ − η of the zeros are simple, on the line, and κ-isolated from the other zeros on the line, for some κ = κ(η) > 0 (not needed here).

## 10. Inputs and their status

| input | where | status |
|---|---|---|
| BH Theorem B (joint CLT), Lemmas 7, 9 | Lamzouri l.873–910 | cited by Lamzouri as unconditional for Dirichlet L-functions (Selberg's density theorem); not re-read (Duke 1995) |
| N_F(T) ~ (T/2π) log T | BH (3.12); Dousselin (eq:N) | classical |
| δ_j = p(c) (stage 1) | IsoZero cert.py, DERIVATION §8 | ours, reviewed (discoveries/linear-combinations-line.md) |
| 2610.08965 certificate (x, b), C'_0 | thm:line-general, prop:gram2, prop:local2 | computer-assisted (535,332,163 boxes, Lean native_decide); our audit Riemann/Deep/ZetaZeros: no flaw |
| dyadic Lem:GaPa | Lamzouri l.968–975 (from Garunkštis–Paliulionytė Lemma 3.1) | as stated by Lamzouri; the same input as the stage-1 side result |
| (P) pair counts | Selberg 1946 (S(t+h,χ) − S(t,χ) mean square); Fujii 1974 | classical, quoted, not re-derived here |
| Theorem N, N^R, Lemma M | §4, §8, §9 | ours (this round) |


## 11. Numerical test (numerics/, logs/zgrid.log, logs/analyse.log, logs/argcount.log)

q = 7, the three odd primitive characters (Conrey labels 3, 5, 6: χ, χ̄ of order 6 and the quadratic character), F with b = (1, 1, 1), (1, −2, ½), (3, 1, −1);
Arb Hardy Z-functions on t = 5 + 0.03k ≤ 2400 (zgrid.py, 3 × ~400 s, one thread); G = Σ b_j Z_j; exact N_F by the argument principle (argcount.py; integers to 10⁻⁴).
11.1 Grid resolution: the sign changes of each Z_j match the Riemann–von Mangoldt main term to ±0.2% on (5, 600], (600, 1200], (1200, 2400] (all zeros of these L_j are found on the line).
11.2 Zeros of F on the line (odd multiplicity) / all zeros, t ∈ (5, 2400]: 2371/2631 = 0.901, 1852/2632 = 0.704, 1839/2631 = 0.699; on (1200, 2400]: 0.899, 0.713, 0.698.
Off-line zeros are 10–30% at these heights; the proven asymptotic lower bounds 0.0206 / 0.0175 are not contradicted (nor is the conjectured 100%, approached like 1/√(log log T)).
11.3 Mechanism: blocks of u-length 10; dominance shares (b = (1,1,1)) 0.329, 0.330, 0.341 (prediction 1/3); the fraction of zeros of F within 0.15 mean spacings of a zero of the
dominant L_j rises with the top-two gap of block medians of log|b_j Z_j| (b = (3, 1, −1): 0.648, 0.871, 0.992 for gap > 0, 1, 2); in all 304 blocks with gap > 2 (three b's) the count of
zeros of F is ≥ (0.35-isolated zeros of the dominant L_j) − 2, the per-block inequality of 4.4.
11.4 Missing-lemma test (7.4): the share of the non-isolated zeros of L_j (κ = 0.35, 0.5) lying in L_j's own dominance blocks is below the share of all its zeros there
(e.g. b = (1,1,1), κ = 0.35: 0.184 vs 0.254, 0.160 vs 0.251, 0.174 vs 0.265): no concentration of the kind 7.2 would need; at this height the evidence favours θ_j ≤ 1/N.


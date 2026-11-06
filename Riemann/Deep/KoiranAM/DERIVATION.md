# Koiran's "Hilbert's Nullstellensatz is in AM" without GRH: derivation

Agent koiran, 2026-11-06. Sources read: Koiran's DIMACS TR 96-27 (the J. Complexity 12 (1996) 273–286 paper with the
same theorem numbers; PostScript from dimacs.rutgers.edu, text extracted locally), Garg–Oliveira–Saxena 2503.20071,
Baldwin–Sivek 2611.06705, Frei–Sofos–Widmer 2610.18491, Rojas–Zhu 1803.04104, Ait El Manssour–Balaji–Nosan–
Shirmohammadi–Worrell 2408.13027, Pak–Robichaux 2412.02064, 2504.03004, 2509.16467, Pak–Robichaux–(Xu) 2508.15715, all in
TeX. "K96" = Koiran 1996, "K97" = Koiran FOCS 1997.

**Status.** The implication "Cureton–Thorner ⇒ Koiran's HN ∈ AM unconditionally" is already stated, as a hedged one-line
remark without details, by Andrews–Garg–Oliveira–Schost, arXiv 2611.05289v1 (3 Nov 2026), sections/introduction.tex l.18,
footnote: CT "claims to prove an effective version of the prime ideal theorem, which would yield an unconditional version
of Koiran's result". Rojas, math/0301111 (2003), had already shown HN ∈ AM under two hypotheses weaker than GRH (DZH,
GIPIT; §9). What follows is the checked derivation, with the analytic constants made numeric, and the descendants.

## 0. Statements

**Problem HN.** Input: a system S: f_1 = ⋯ = f_s = 0 with f_i ∈ ℤ[X_1, …, X_n], each f_i written as its list of nonzero
terms (integer coefficient in binary, exponent vector in binary or unary). |S| is the total bit length. Question: has S
a solution in ℂⁿ (equivalently in ℚ̄ⁿ)? Notation of K96: d = max deg f_i, L = max bit size of a coefficient,
δ = 2 + Σ deg f_i. In every encoding used below (dense, sparse, and the circuits of GOS and Mulmuley),
n, s, log d, log δ, log L = O(|S|), and f_i(z) mod p can be evaluated in time poly(|S|, log p) (repeated squaring for
binary exponents, gate by gate for circuits).

**Theorem 1.** HN ∈ AM, unconditionally (given the Cureton–Thorner theorem, §2). The protocol is Koiran's (K96 Thm 2) with
perfect completeness; only the range of primes changes, to p ≤ x₀ with
x₀ = max{e^{40}, (7.15(log D + 6.43))^{32}, 52.06(eB′)²} (§4: e, D, B′ are K96's degree, discriminant and gap
parameters, all 2^{poly(|S|)}), so log x₀ = poly(|S|). Hence HN ∈ Π₂ᵖ (AM ⊆ Π₂ᵖ, Babai).

**Theorem 2.** DIM (given S and r, is dim Z(S) ≥ r?) ∈ AM, unconditionally (K97 Thm 4.1 minus GRH, §6).

**Corollaries** (§7): HNZ = HNP (parametric HN) ∈ AM; non-vanishing of genus-0 3-point Gromov–Witten invariants of partial
flag varieties ∈ AM; primality of radical or equidimensional Cohen–Macaulay ideals given by circuits ∈ Σ₃ᵖ ∩ Π₃ᵖ (coAM when
the dimension is given); resultant vanishing ∈ AM; Noether normalization and hitting-set certification for VP in PH;
and, through the same prime count used to eliminate complex constants, Bürgisser's characteristic-zero transfer
(VP_ℂ = VNP_ℂ ⇒ #P ⊆ P/poly, i.e. Cook's hypothesis implies Valiant's over ℂ) with its consequences (Koiran–Perifel,
Fournier–Perifel–de Verclos); all without GRH. Schubert vanishing ∈ AM ∩ coAM also follows, but it is already known
unconditionally in the stronger form SV ∈ coRP (Pak–Robichaux 2509.16467).

## 1. Every use of GRH in K96

Dependency chain of K96 (theorem numbers of the TR; "unc." = proved without hypotheses):

| K96 item | content | hypothesis |
|---|---|---|
| Thm 3, Lemma 1, Thm 4 | primitive element: x_i = Q_i(r)/a, r a root of an irreducible R ∈ ℤ[X] of degree ≤ dⁿ, coefficients and log a of size L·d^{nO(1)} | unc. |
| Thm 6, Lemma 2, Thm 7 | quantifier elimination (Fitchas–Galligo–Morgenstern): a satisfiable S has a solution with each x_i a root of a polynomial of degree ≤ 2^{(n log δ)^{c₁}}, coefficient size ≤ L·2^{(n log δ)^{c₂}} | unc. |
| Lemma 3 | if R has a root in 𝔽_p and p ∤ a, then S is satisfiable in 𝔽_pⁿ (g_i := a^{d_i}f_i(Q(X)/a) = R·A_i in ℤ[X], Gauss) | unc. |
| Thm 5 | S unsatisfiable over ℂ ⇒ \|R_S\| ≤ d^{O(n)} s(⌈log s⌉ + L) (arithmetic Nullstellensatz: an integer a ≠ 0 in the ideal; p ∣ a for p ∈ R_S) | unc. |
| **Thm 9** | Adleman–Odlyzko: \|Σ′_{p≤x}(1 − W(p))\| = O(x^{1/2} log(Δ xⁿ)), W(p) = #roots of f mod p, Σ′ over p ∤ Δ | **GRH**, "through an effective version of the Chebotarev density theorem [Lagarias–Odlyzko] (in fact, only the effective prime ideal theorem is used)" |
| Cor 1 | π_f(x) ≥ (1/n)[π(x) − log Δ − c x^{1/2} log(Δ xⁿ)] for f irreducible of degree n, discriminant Δ | from Thm 9 |
| Thm 8 | S satisfiable ⇒ π_S(x) ≥ (π(x) − c₆ n x^{1/2} log x)/(L 2^{(n log δ)^{c₄}}) − L 2^{(n log δ)^{c₅}} x^{1/2} | Thms 4, 7, Lemma 3, Cor 1; minus the ≤ log a primes dividing a |
| Thm 1 | A := d^{c₁n}s(⌈log s⌉ + L), x₀ ≥ L^{c₂}2^{(n log δ)^{c₃}}: unsat ⇒ π_S(x₀) ≤ A; sat ⇒ π_S(x₀) ≥ B := 8A(log A + 3) | Thms 5, 8 |
| Thm 2 | HN ∈ AM (Sipser hashing, §5) | Thm 1 |

R_S is the set of primes p for which S has a solution in 𝔽_pⁿ, π_S(x) = #(R_S ∩ [2, x]). The TR mentions GRH in five places
(searched in the extracted text with line breaks removed): the abstract, the plan of the paper ("[section 5] requires
GRH"), the opening of §4 ("the analysis of satisfiable systems is more involved, and requires [GRH]"), Thm 9, and the
closing remark of §6 (a polynomial bound on |R_S| "would imply under GRH that HN ∈ NP", via Thm 8). So **GRH enters K96 only through Thm 9**, i.e. through a lower bound
for the number of degree-one prime ideals of K = ℚ[X]/(R) of norm ≤ x. The restatements agree: GOS 2503.20071
preliminaries.tex l.455–466 ("this last result is conditional ... and this is also what makes our result conditional"),
Baldwin–Sivek 2611.06705 l.203 ("Koiran finds [the prime] (assuming GRH) by appealing to Lagarias–Odlyzko"). Lemma A below
replaces Thm 9 + Cor 1; the rest of K96 is used verbatim.

## 2. The analytic input

- **(CT)** Cureton–Thorner (Hexagon 2610.00200v2; Lean development github.com/AxiomMath/DedekindZetaZeroFree, theorem
  `dedekindZeta_ne_zero`, no `sorry`, checked in R399): for every number field K, ζ_K(s) ≠ 0 for Re s > 31/32. By the
  functional equation every nontrivial zero satisfies 1/32 ≤ β ≤ 31/32.
- **(FSW)** Frei–Sofos–Widmer 2610.18491 Thm 1.1 (main.tex l.194–206): if ζ_L has no zero with Re s > 1 − δ
  (0 < δ ≤ 1/2) then, with an effectively computable absolute constant c₁, for every Galois L/K, class C and x > 2,
  |π_C(x, L/K) − (#C/#G) Li(x)| ≤ c₁ (#C/#G) x^{1−δ} log(Δ_L^{1/δ} x^{n_L}). With L = K, C = {1}:
  |π_K(x) − Li(x)| ≤ c₁ x^{31/32} log(Δ_K^{32} x^{n_K}) (this is Baldwin–Sivek's Cor. cor:landau, l.379–390).
- **(ChebExp)** For numbers instead of c₁: the smoothed explicit formula of Riemann/Deep/ChebExp/DERIVATION.md §§2–4 with
  L = E = K, H = {1} (no Chebotarev condition; Ram = 0 because L/K is trivial). Its (5.1) reads, for x > e^U,

  Good(x) ≥ xΦ(1) − x^{31/32}M_x(½ log Δ_K + 1/σ₀ + 1/(σ₀ − 1) − n_K g(σ₀)) − (n_K − 1)Φ(0) − n_KΦ(0)τ_x − n_K P_k(x),

  Good(x) := Σ_{𝔭 of degree one} log N𝔭 · k(N𝔭/x), with the kernel k ≥ 0 supported on [e^{−U}, 1], k ≤ ‖k‖.

## 3. Lemma A (replaces K96 Thm 9 and Cor 1)

**Lemma A.** Let f ∈ ℤ[X] be irreducible of degree e ≥ 1 with leading coefficient b, put g(Y) := b^{e−1}f(Y/b) (monic, in
ℤ[Y]) and D := |disc g| = |b|^{(e−1)(e−2)}|disc f|. Let π_f(x) be the number of primes p ≤ x such that f has a root in 𝔽_p.
Then for every x > 2

  (A-FSW)  π_f(x) ≥ (1/e)[Li(x) − c₁x^{31/32}(32 log D + e log x) − (e/2)π(√x)] − ω(b),

and for every x ≥ e^{40}

  (A-num)  π_f(x) ≥ (x/(e log x))·[0.2772 − 0.4953 (log D + 6.43) x^{−1/32} − 0.7233 e x^{−1/2}] − ω(b).

*Proof.* Let θ be a root of f, K = ℚ(θ), so [K:ℚ] = e and bθ ∈ O_K is a root of g.
1. *A degree-one prime gives a root.* If 𝔭 is a prime of K with N𝔭 = p and p ∤ b, then O_K/𝔭 = 𝔽_p and ȳ := bθ mod 𝔭
   satisfies ḡ(ȳ) = 0; since f(X) = b^{−(e−1)}g(bX), x̄ := ȳ/b̄ ∈ 𝔽_p is a root of f mod p. (No condition on p ∣ disc is
   needed: we never use the converse.)
2. *Counting.* Above each p lie at most e primes of degree one, so π_f(x) ≥ (1/e)·#{𝔭: deg 𝔭 = 1, N𝔭 ≤ x, N𝔭 ∤ b}. The
   primes of degree one above the ω(b) primes p ∣ b number ≤ e·ω(b). A prime of degree f ≥ 2 with N𝔭 = p^f ≤ x has
   p ≤ √x, and above each p there are ≤ e/2 of them (Σ e_i f_i = e). Hence
   #{deg 𝔭 = 1, N𝔭 ≤ x} ≥ π_K(x) − (e/2)π(√x).
3. *Discriminant.* g is monic with root bθ ∈ O_K, so disc g = [O_K : ℤ[bθ]]²·disc K and Δ_K ≤ D. The formula
   disc g = b^{(e−1)(e−2)} disc f follows from disc g = Π_{i<j}(bθ_i − bθ_j)² and disc f = b^{2e−2}Π_{i<j}(θ_i − θ_j)².
4. *(A-FSW)*: (FSW) with δ = 1/32, which (CT) allows for this K: π_K(x) ≥ Li(x) − c₁x^{31/32}(32 log Δ_K + e log x).
5. *(A-num)*: (ChebExp) with the kernel (h, s*, σ₀) = (1.5, 0, 1.4), U = 4.5. Arb enclosures (count_cert.py,
   count_cert_output.txt): Φ(1) ≥ 0.46886, Φ(0) ≤ 3.375, ‖k‖ ≤ 1.69088, A₁ ≤ 1.67494, A₂ ≤ 6.09547, g(1.4) ≥ 0.9496,
   1/σ₀ + 1/(σ₀ − 1) ≤ 3.21429, and the prime-power constants S₂, S₃, Φ(½), Φ(⅓). For x ≥ e^{40}: M_x = A₁ (because
   x^{15/16} > A₂/A₁), τ_x ≤ ε_x ≤ 1, and P_k(x) ≤ 1.01624[x^{1/2}(S₂ + Φ(½)/2) + x^{1/3}(S₃ + Φ(⅓)/3) +
   (J − 3)x^{1/4}(‖k‖ + Φ(0)/4)] (ChebExp (4.3), θ(y) < 1.01624y), so that e(Φ(0)(1 + τ_x) + P_k(x)) ≤ 1.22241·e·x^{1/2}.
   Dropping −n_K g(σ₀) ≤ 0 from B(σ₀) and using Good(x) ≤ ‖k‖ log x · #{deg 𝔭 = 1, xe^{−4.5} ≤ N𝔭 ≤ x}:
   #{deg 𝔭 = 1, N𝔭 ≤ x} ≥ (x/log x)[Φ(1)/‖k‖ − (A₁/(2‖k‖))(log Δ_K + 2·3.21429)x^{−1/32} − (1.22241/‖k‖) e x^{−1/2}].
   Here ‖k‖ may be replaced by its certified upper bound ‖k‖⁺ throughout: if G⁻ denotes the lower bound for Good(x),
   then N₁ := #{deg 𝔭 = 1, N𝔭 ≤ x} ≥ G⁻/(‖k‖⁺ log x) both when G⁻ ≥ 0 (‖k‖ ≤ ‖k‖⁺) and when G⁻ < 0 (N₁ ≥ 0). The
   script checks in Arb that Φ(1)⁻/‖k‖⁺ ≥ 0.2772, A₁⁺/(2‖k‖⁺) ≤ 0.4953, 2·3.21429 ≤ 6.43, 1.22241/‖k‖⁺ ≤ 0.7233.
   Primes of degree ≥ 2 are inside Pow, so N₁ counts degree-one primes only, and π_f(x) ≥ N₁/e − ω(b) by step 2. ∎

*Squarefree f.* If f is squarefree, apply Lemma A to any irreducible factor f₁ ∈ ℤ[X]: e₁ ≤ e, D₁ ≤ |disc of g₁|, and
π_f ≥ π_{f₁}. (disc(uv) = disc(u) disc(v) Res(u, v)² for u, v ∈ ℤ[X] shows |disc f₁| ≤ |disc f|.) This is the form
"Lemma A" of the brief, with Koiran's −log disc term replaced by −ω(b).

*Comparison.* K96 Cor 1 under GRH: error c x^{1/2} log(Δxⁿ). Lemma A: error x^{31/32}(32 log D + e log x) up to c₁, or
0.4953(log D + 6.43)x^{31/32}/log x with numbers.

## 4. Theorem 8′ and Theorem 1′ (the gap, with explicit x₀)

Let S be satisfiable. Take R (irreducible, primitive, degree e, leading coefficient b) and a ∈ ℕ from K96 Thms 4, 7 and
Lemma 3. By Lemma 3, every prime p ∤ ab at which R has a root lies in R_S. Lemma A (A-num) for f = R gives, for x ≥ e^{40},

  π_S(x) ≥ (x/(e log x))[0.2772 − 0.4953(log D + 6.43)x^{−1/32} − 0.7233 e x^{−1/2}] − log₂|ab|.   (Thm 8′)

**Theorem 1′.** Let A be K96's bound (Thm 5), B₀ := 8A(⌈log₂A⌉ + 3), B′ := B₀ + log₂|ab|, and

  x₀ ≥ X(S) := max{ e^{40}, (7.15(log D + 6.43))^{32}, 52.06(eB′)² }.

Then S unsatisfiable ⇒ π_S(x₀) ≤ A, and S satisfiable ⇒ π_S(x₀) ≥ B₀.

*Proof.* Unsatisfiable: K96 Thm 5 (every x). Satisfiable: x₀^{1/32} ≥ 7.15(log D + 6.43) makes the second term of the
bracket ≤ 0.4953/7.15 < 0.0693 = 0.2772/4; x₀ ≥ 52.06e²B′² ≥ 109e² (B′ ≥ B₀ ≥ 24) makes the third ≤ 0.7233/10.44 <
0.0693; so the bracket is ≥ 0.1386. For x₀ ≥ e^{40}, log x₀ ≤ √x₀, so π_S(x₀) ≥ 0.1386√x₀/e − log₂|ab|, and
√x₀ ≥ √52.06·eB′ = 7.2153eB′ gives π_S(x₀) ≥ 1.00004B′ − log₂|ab| ≥ B₀. The subtracted log₂|ab| bounds ω(a) + ω(b):
primes p ∣ a are excluded by Lemma 3, primes p ∣ b by Lemma A. (count_cert.py prints the thresholds 7.14719, 10.4372,
7.21501 = 2/0.2772 and 52.0563 = 7.21501².) ∎

With (A-FSW) instead: the same holds for x₀ := (2^{20}Q)^{64}, Q := max{256c₁ log D, 8c₁e, 10.08e, 8e(B₀ + log₂|ab|)}:
each of the three error terms of (A-FSW) is then ≤ x₀/(8 log x₀) (using x₀^{1/32} = (2^{20}Q)² ≥ Q(log x₀)² for Q ≥ 1,
Li(x) ≥ x/(2 log x) for x ≥ 4, π(y) < 1.25506y/log y), so π_S(x₀) ≥ x₀/(8e log x₀) − log₂|ab| ≥ B₀.

**Arthur must compute x₀ from S alone.** Replace e, log D, log₂|ab| by K96's a-priori bounds (Thms 4, 7): with
E := 2^{n(n log δ)^{c₁}} ≥ e, H := L·2^{(n log δ)^{c}} ≥ log ‖R‖, log |a| ≤ H, Hadamard on the Sylvester matrix gives
|disc R| ≤ ‖R‖₂^{e−1}‖R′‖₂^{e}, so log D ≤ (e−1)(e−2)H + (2e−1)H + 2e log(e+1) ≤ 3E²(H + log(E+1)), and log₂|ab| ≤ 3H.
All of log E, log H, log A are poly(|S|), so log x₀ = O(log H + log E + log A) = poly(|S|): x₀ = L^{O(1)}2^{O((n log δ)^{c'})},
the same shape as K96 Thm 1 (x₀ ≥ L^{c₂}2^{(n log δ)^{c₃}}) with new absolute constants. The constants c, c₁, c' of K96's
elimination bounds are unspecified there; nothing number-theoretic remains unspecified with (A-num).

**Size against GRH.** Under GRH the analytic condition is x^{1/2} ≳ (log D + e log x) log x and x/(e log x) ≳ B′, so
x₀ ≈ max((log D)², eB′)·polylog; here x₀ ≈ max((7.15 log D)^{32}, 52(eB′)²). The bit length of Merlin's prime grows by a
factor ≤ 16 + o(1) (1/2 ↦ 1/32); every step stays polynomial.

## 5. The AM protocol (K96 Thm 2, unchanged) and its gap at x₀

E := R_S ∩ [2, x₀], coded as k-bit strings, k := ⌈log₂x₀⌉ + 1 = poly(|S|). Membership in E is in NP: the witness is
z ∈ 𝔽_pⁿ with f_i(z) ≡ 0 (mod p); primality is in P (AKS) or NP (Pratt). Let m be the integer with A ≤ 2^{m−2} < 2A, so
m < log₂A + 3 and m2^m < 8A(log₂A + 3) ≤ B₀.

Arthur: random binary m×k matrices M_1, …, M_m (h_i(x) := M_i x, a linear map 𝔽₂^k → 𝔽₂^m). Merlin: x, x_1, …, x_m ∈ E with witnesses and
h_i(x) = h_i(x_i), x ≠ x_i for all i. Arthur: checks primality, x, x_i ≤ x₀, the witnesses, and the m collisions.
- *Completeness (perfect).* S satisfiable ⇒ |E| ≥ B₀ > m2^m. For any h_1, …, h_m, an element alone in its h_i-fibre is
  one of ≤ 2^m per i, so ≤ m2^m elements are alone for some i; some x ∈ E is alone for no i. Merlin always succeeds.
- *Soundness.* S unsatisfiable ⇒ |E| ≤ A ≤ 2^{m−2}. For x ≠ x′, h_i(x − x′) is uniform on 𝔽₂^m, so
  Pr[∃ valid answer] ≤ Σ_{x∈E} Π_i (|E|2^{−m}) ≤ 2^{m−2}·4^{−m} = 2^{−m−2} ≤ 1/2 (Sipser's lemma as used in K96).
- *Gap.* yes ≥ B₀ = 8A(⌈log₂A⌉ + 3), no ≤ A: ratio ≥ 8(⌈log₂A⌉ + 3) ≥ 24 (Koiran's own B; the Goldwasser–Sipser form
  of Arora–Barak §8.4.1, |S| ≥ K vs |S| ≤ K/2, also applies with K := 2A ≤ B₀). In fact Thm 8′ gives
  π_S(x₀) ≥ 0.1386√x₀/e − log₂|ab|, far above B₀ when the (log D)^{32} term sets x₀.
- *Cost.* All messages and checks are poly(k, |S|). Randomness precedes nondeterminism: this is an AM protocol, with
  one-sided error, which descendants use (2408.13027 l.176–207; 2508.15715 l.1207).

## 6. DIM ∈ AM (K97 Thm 4.1)

K97 decides dim Z(S) ≥ r by a randomised polynomial-time reduction to HN: a random linear change of coordinates, then
random integer values for the first r coordinates (2408.13027 geometric.tex l.5–9 describes it; their own algebraic
version is DIM → HNZ → HN, each step a randomised polynomial-time reduction with errors from the Schwartz–Zippel lemma,
l.27–33). No number theory occurs in the reduction. With Theorem 1 and closure of AM under randomised reductions with
bounded error (BP·AM = AM; 2408.13027 l.217 cites Arora–Barak §8.4.2, using the perfect completeness of §5), DIM ∈ AM.
Hence the exact dimension is computable in P^{AM} ⊆ P^{Π₂ᵖ} (binary search), which is what GOS use.

## 7. Descendants: where GRH enters, and what holds now

| result | GRH enters | without GRH |
|---|---|---|
| K96 §6 remark: a polynomial bound on \|R_S\| would put HN in NP "under GRH" | Thm 8 | the implication is unconditional (moot: K96 Thm 10 shows \|R_S\| can be exponential) |
| K97 Thm 4.1: DIM ∈ AM | only via HN ∈ AM (§6) | **unconditional** |
| 2408.13027 Thm 1: HNZ ∈ AM (parametric HN over the closure of ℚ(x)) | only via Koiran's protocol (specialisations.tex l.174–217: random integer specialisation, then "Koiran's algorithm"); the specialisation bounds (Props sat/unsat-spec-bound, app. A from K96 Thm 6–7) are algebraic | **unconditional** |
| 2412.02064 Thm 1.1: Schubert vanishing SV ∈ coAM (types A–D) | only via HNP ∈ AM (= 2408.13027): Main Lemma ¬SV ≤ HNP is algebraic; other GRH mentions (l.783, 1206, 1710–1798, 3835) are corollaries or discussion | unconditional, **but superseded**: SV ∈ coRP unconditionally (Pak–Robichaux 2509.16467, Sep 2025, l.1141–1149), and coRP ⊆ BPP ⊆ AM ∩ coAM |
| 2504.03004: SV ∈ AM ∩ coAM | via HNP ∈ AM (SV ≤ HNP and ¬SV ≤ HNP, l.1580–1594) | same: superseded by 2509.16467 |
| 2508.15715 Main Thm: non-vanishing of genus-0 3-point GW invariants of partial flag varieties ∈ AM; Thm t:main-HNPE: HNPE ∈ AM | only via HNP ∈ AM (proof of HNPE l.1203–1228 needs perfect completeness: kept) | **unconditional**; hence GW vanishing ∈ coAM, not NP-hard unless PH = Σ₂ᵖ (Boppana–Håstad–Zachos), in coNP under Impagliazzo–Wigderson's assumption alone (l.1418–1432). 2508.15715 says the AM ∩ coAM tools of the Schubert case do not extend to GW; 2509.16467 does not treat GW |
| 2503.20071 Thm primesPH: primality of radical or equidimensional Cohen–Macaulay ideals (circuits of size s) in Σ₃ᵖ ∩ Π₃ᵖ; coAM if dim I is given | three places: Koiran's HN ∈ AM (alg. radical-non-prime, l.160), K97 DIM ∈ AM (l.215, l.264), and **K96 Cor 1 directly** in Lemma zero-dimensional-reducible (l.19–48) and Lemma two-top-components (l.50–87) | **unconditional**: in both lemmas replace Cor 1 by Lemma A; their error terms c x^{1/2}(h 2^{(n log σ)^{c}} + log x) become c x^{31/32}(...), so x₀ = h^{c₉}2^{(n log σ)^{c₁₀}} works with constants ×16 (x₀^{1/32} against the same poly-exponential factor); Lang–Weil (Cafure–Matera) and the height bounds are unconditional |
| Rojas–Zhu 1803.04104: MDZH ⇒ MRH ⇒ DIM ∈ P^{NP^{NP}} | their hypothesis MDZH | superseded: DIM ∈ AM ⊆ Π₂ᵖ. Side remark: MDZH as printed (main.tex l.269–280) fails for K = ℚ for every C > 4 ((log log 3)^C < 10⁻⁴ puts β = 1/2, γ = 14.13 inside the claimed zero-free region), so their Thm (1) "GRH ⇒ MDZH" is false as stated; a statement slip (it holds for d log 3Δ large), harmless for their argument |
| Kuperberg 2014, Zentner, Baldwin–Sivek (knots and 3-manifolds in coNP) | existence of one small prime with a root (Kuperberg Thm 3.2) | done by Baldwin–Sivek 2611.06705 (stated modulo QRH(ε)), not ours |
| Grenet–Koiran–Portier 1210.1451 Cor corr:in-am: resultant vanishing (homogeneous HN) ∈ AM | only via HN ∈ AM ("follows from a simple reduction to the Hilbert Nullstellensatz", l.171–173, 258–263) | **unconditional** |
| Mulmuley GCT V 1209.5993 Thm tnoetherph: constructing a homogeneous system of parameters (Noether normalization) for a variety given by circuits for its defining equations ∈ PH (EXPH variants l.787, 1221, 1817) | via Thm tkoiran (= K96 for circuits, homogeneous variant; l.1744–1750, used l.1975–1980); other GRH mentions (l.257–277, 1732–1739, 1888–1917) are these statements | **unconditional** as far as read (only the cited lines checked) |
| Kopparty–Moshkovitz–Zuiddam 2002.09472 l.666 | quotes K97 (DIM ∈ AM under GRH) as context; no theorem of theirs depends on it | nothing to transfer |
| Qiao et al. 1904.03950 Prop (l.982–986, 4074–4082): maximum isotropic space and isotropic 3-decomposition over ℂ ∈ PH | via K96 HN ∈ AM | **unconditional** (PH part). Their Thm thm:iso_dim_2 (l.794–800: randomized reduction from quadratic residuosity, over ℚ, under GRH) uses GRH elsewhere: not checked |
| **Bürgisser 2000** (book, Ch. 4; Thm 4.4 p.64): characteristic-zero transfer: VP_ℂ = VNP_ℂ ⇒ #P ⊆ P/poly (2406.06217 l.1885, 2606.25121 l.295), hence P/poly = NP/poly = PH/poly = CH/poly; equivalently #P ⊄ P/poly ⇒ VP_ℂ ≠ VNP_ℂ ("Cook's hypothesis implies Valiant's" over ℂ) | only through the elimination of complex constants, i.e. a small prime p at which a satisfiable integer system stays solvable (Bürgisser's own description, 2406.06217 l.2030–2034: "used for the elimination of potentially helpful complex constants, similarly as it was done by Koiran"); the lemma as quoted by Fournier–Perifel–de Verclos 1304.5910 l.355–359 is π_S(a) ≥ π(a)/d^{O(n)} − √a log(wa) (K96 Thm 8 type) | **unconditional** (book not read; mechanism from the quotes): (A-FSW) applied to Koiran's representation (e ≤ E) gives π_S(a) ≥ Li(a)/E − c₁a^{31/32}(32 log D + E log a) − √a − log₂\|ab\| with log E, log log D, log log\|ab\| = poly, so the prime of Thm 4.4 exists below 2^{poly(n)} (KP's p ≤ 2^{n²q(n)}, from Bürgisser's d^{O(n)} representation, becomes 2^{poly(n)} with Koiran's; Lemma B, §10). Still called conditional on 20 Oct 2026 (2606.25121 l.294–296) and on 6 Nov 2026 (2611.07972 sections/introduction.tex l.735–738: "In characteristic zero, the same implication is known only under the Generalized Riemann Hypothesis") |
| Koiran–Perifel cs/0610009, 0706.1477 Prop prop_equivalence: VPSPACE = VP_nb ⇔ [P/poly = PSPACE/poly and VP = VNP] (reals, ℂ) | the constants elimination via Bürgisser Thm 4.4 (cs/0610009 l.567–625); the consequence NC/poly = PSPACE/poly (l.149–153) via Bürgisser 2000 | **unconditional** |
| Fournier–Perifel–de Verclos 1304.5910: polynomials with coefficients in MA without O(n^k) circuits over ℂ; Lemma collapse (VP = VNP over ℂ ⇒ CH = MA) | Bürgisser's lemma (l.355–359), K96 (HN ∈ PH, l.386–388), Bürgisser Ch. 4 (l.363–367); "GRH is used to eliminate the complex constants" (l.161) | **unconditional** |
| Bürgisser 2406.06217 Thm (2) (l.1857–1888): VP^F = VNP^F ⇔ VP_nb^F = VNP_nb^F in characteristic 0 | via Bürgisser 2000 (#P ⊆ P/poly) | **unconditional** |
| 2309.07612 (succinct hitting sets) Thm l.208, l.289, Prop l.878–887 | via Bürgisser 2009 ("if VP = VNP, then CH collapses to P/poly", l.886) and Koiran–Perifel Prop 3 | probably unconditional (Bürgisser 2009 not read; same constants elimination expected) |
| Hitting-set certification for VP in PH (2610.00644 l.205–215, from Mulmuley and Forbes–Shpilka) | via K96 | **unconditional** (PH; 2610.00644 notes only the unconditional CH bound) |
| Bürgisser–Cucker 2006 Cor 8.7 (cited in 2412.02064 l.3835) | the arXiv version cs/0312007 mentions GRH only when quoting K96 (l.664) | not checked |

## 8. What is not claimed

- Not the observation itself: 2611.05289 states it (Status above). Not DIM either: the same footnote covers it implicitly
  (2611.05289 l.27 says K97 needs the same "effective prime ideal theorem").
- No new analytic input: the theorem is K96 with Thm 9 replaced by Lemma A. Baldwin–Sivek did the same substitution for
  the existence of one prime (Kuperberg's lemma); we need the count, which their (eq:lower-bound-unconditional) and
  Prop. goal-general also give (up to the 1/d and the root bound W(p) ≤ d).
- (A-FSW) rests on FSW's unrefereed but detailed proof; (A-num) rests on the ChebExp derivation (reviewed, correct) and
  its Arb kernel class; both rest on (CT), Lean-checked but unrefereed.
- The elimination-theory constants of K96 (c₁, …, c₆) are not made numeric; they are GRH-free.
- K96's encoding: the TR speaks of systems of size |S| with the parameters d, L, δ; GOS apply it to circuits and
  Pak–Robichaux to sparse systems; the bounds above need only log d, log δ, log L ≤ poly(|S|) and fast evaluation mod p.

## 9. Prior hypotheses weaker than GRH (Rojas 2003, Rojas–Zhu 2018) and the strip

- Rojas math/0301111 Thm thm:big (l.341–349): HN ∈ AM if DZH or GIPIT holds. GIPIT (l.596–602): with ℓ := 1 + log(n_K log d_K)
  and C := 1.01 + ℓ^κ, π_K((x+1)^C) − π_K(x^C) ≥ 1 + (n_K/2)(x+1)^{C/2} for all x ≥ ℓ^κ. (CT)+(FSW) give GIPIT for every
  K ≠ ℚ once κ is large: the interval has length ≥ Cx^{C−1}, the main term is ≥ x^{C−1}/log(x+1) and the error is
  ≤ 2c₁(x+1)^{31C/32}(32 log d_K + n_K C log(x+1)); since every K ≠ ℚ has ℓ ≥ 1 + log(2 log 3) = 1.787, a large κ makes
  C/32 − 1 ≥ ℓ^κ/64 and x^{C/32−1} ≥ 2^{ℓ^κ/64} beats e^{ℓ}·poly(x, C) ≥ the error/main ratio (log d_K ≤ e^{ℓ−1}).
  (GIPIT is undefined for K = ℚ, where log(n_K log d_K) = log 0; there the prime number theorem in short intervals is
  unconditional.) DZH (l.555–595) asks for |S_K(x, T)| ≤ x^{1−1/(1.01+ℓ^{κ})}T^{1.99/(1.01+ℓ^κ)} for all x, T ≥ ℓ^κ; the
  strip does not give this near x ≈ ℓ^κ (it saves x^{1/32} against a factor log d_K), but Rojas only needs GIPIT.
  So Rojas 2003 + CT is a second route to Theorem 1.
- Rojas–Zhu 1803.04104 Thm thm:main: GRH ⇒ MDZH ⇒ MRH ⇒ DIM ∈ P^{NP^{NP}} (Stockmeyer counting, not AM). Superseded.
  Their MDZH (l.269–280) is false for K = ℚ as printed (§7 table), hence their implication (1) is false as stated; the
  intended statement (for n_K log Δ_K above a constant) is fine and is implied by (CT).

## 10. Bürgisser's constants elimination without GRH

**Lemma B** (replaces Bürgisser 2000, p.64, as quoted in 1304.5910 l.355–359: under GRH,
π_S(a) ≥ π(a)/d^{O(n)} − √a log(wa)). Let S be a system over ℤ in n unknowns, solvable over ℂ, and let R (degree e),
a, b, D be as in §4. Then for every a₀ ≥ e^{40}

  π_S(a₀) ≥ (a₀/(e log a₀))[0.2772 − 0.4953(log D + 6.43)a₀^{−1/32} − 0.7233 e a₀^{−1/2}] − log₂|ab|,

and with c₁: π_S(a₀) ≥ Li(a₀)/e − c₁a₀^{31/32}(32 log D + e log a₀) − √a₀ − log₂|ab|. This is Thm 8′ (§4); only the
bookkeeping is new. With Koiran's bounds (log e, log log D, log log|ab| polynomial in n, log δ, log L), the number of
primes p ≤ a₀ at which S stays solvable is at least any prescribed N ≥ 2 once (proof of Theorem 1′ with B′ := N + log₂|ab|)
a₀ ≥ max{e^{40}, (7.15(log D + 6.43))^{32}, 52.06e²(N + log₂|ab|)²}, which has poly(n, log δ, log L, log N) bits.

*Use in the transfer.* In the constants-elimination step (cs/0610009 l.588–625, "we proceed as in [Bürgisser 2000]"),
S expresses that a circuit with unknown constants Y computes the target polynomial on all Boolean inputs: 2ⁿ equations,
degree ≤ 2^{q(n)}, weight ≤ 2^{2^{q(n)}}, so log δ ≤ n + q(n) + 1 and L ≤ 2^{q(n)}, and a₀ has poly(n) bits. A prime
p ≤ a₀ with S solvable over 𝔽_p gives a Boolean circuit of size poly(n) computing the target modulo p; taking N := π(2^{t}) + 1
with t ≥ the bit length of the target values (e.g. t = n log n for the permanent of a 0/1 matrix) forces p > 2^t, so the
values are recovered exactly. Hence every statement whose only GRH input is this step holds without GRH. Not read:
Bürgisser's book itself; the claim that its GRH use is this step rests on 2406.06217 l.2030–2034 (Bürgisser), cs/0610009
l.567–625 and 1304.5910 l.159–162, 355–367.

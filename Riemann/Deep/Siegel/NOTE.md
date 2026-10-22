# Kerr, arXiv 2610.21099v1: step-by-step audit (2026-10-21)

Paper: B. Kerr, "A simpler proof that there are no Landau–Siegel zeros" (v1, 15 Oct 2026, the only version; no
arXiv comment). Read from the TeX e-print (`siegel_zero_proof.tex`, 419 lines). Equation numbers below are the
compiled ones: (1) bias, (3)/(4) Frobenius, (9) = `eq:PS`, (11) = `eq:count-condition`.

Other sources read: OpenAI, "Uniform exclusion of Landau–Siegel zeros" (1 Oct 2026, TeX at
github.com/openai/math/…/build/paper.tex); Banks, arXiv 2610.15483 (TeX, structure and abstract);
D. W. Masser, Invent. Math. 63 (1981) 81–95, pp. 81–83 (scan at GDZ, PPN356556735_0063 / LOG_0011; not kept,
GDZ terms forbid redistribution); Tao's blog post of 19 Oct 2026 and its comments; arXiv API listings.
Arithmetic checks: `check.py` (pure Python, < 2 s, < 50 MB) → `check_output.txt`.

## 0. Shape of the proof
Siegel zero ⇒ χ(p) = −1 on almost all log-mass of primes ≤ Q (Lemma 1). Pigeonhole over polynomials
f ∈ L_D (O-coefficients, monomials x1^i x2^j x3^k with i + 8(j+k) ≤ D), using Frobenius θ^p ≡ σθ or στθ (mod p)
at those primes, gives F ≠ 0 with F(θ_n, σθ_n, στθ_n) ≡ 0 (mod P) and |values| ≤ P/2, hence exact vanishing at
all N⁴ points. Masser: deg F ≥ ((N−1)/3)^{4/3} > D ≥ deg F. Contradiction.

## 1. Lemma 1 (analytic input): correct
- s = 1 + 1/log X ∈ (1, 2) for X ≥ 3. Davenport ch. 14: −Re L'/L(s,χ) ≤ c log q − Σ_ρ Re 1/(s−ρ) for 1 < s ≤ 2,
  c absolute, every term ≥ 0 for real s > 1, keep only ρ = β. With −ζ'/ζ(s) = 1/(s−1) + O(1):
  Σ Λ(n)(1+χ(n)) n^{−s} ≤ C log q + 1/(s−1) − 1/(s−β); left side ≥ (2/e) Σ_{p≤X, χ(p)=1} log p/p.
- 1/(s−1) − 1/(s−β) = (1−β)/((s−1)(s−β)) ≤ (1−β) log²X since s−β ≥ s−1. Right error term (sharper:
  min(log X, (1−β) log²X)). Implied constant absolute.
- "Davenport, eq. (2), p. 89": content standard; exact page/number not checkable (book not openly available).
  OpenAI §2 proves the identical bound in full, so nothing hangs on the citation.

## 2. Number field and Frobenius: correct for q ≠ 8
- K = Q(√q, √2) has degree 4 and {1, √q, √2, √2q} is a Z-basis of O = Z[√q, √2] iff q ∉ {1, 8} (up to squares
  q ≠ 2·□). Fundamental discriminants: only q = 8 fails (K = Q(√2), O = Z[√2]; also Lemma 7's hypothesis fails:
  v = (4, √2, −2√2, −1) is Q-dependent). q = −8, −4, 4m are fine (√q = 2√m). Harmless: the proof is only used for
  |q| ≥ q_0. OpenAI uses squarefree d and excludes d = 2 explicitly.
- For q < 0 one must read √(2q) as √2·√q (sign convention). Harmless.
- (√q)^p = q^{(p−1)/2} √q ≡ (q/p) √q, (√2)^p ≡ (2/p) √2 (Euler's criterion; χ(p) = (q/p) for odd p ∤ q),
  p-th power is additive mod p and fixes Z mod p ⇒ (3), (4). The Milne citation is unnecessary.
- check.py §2: 1386 random θ, 12 discriminants (−3, −4, 5, −7, −8, 12, 13, −15, −23, −163, 21, −1155),
  inert p < 400: 0 failures.

## 3. Lemma 2 (basis): correct
Term a of x1^v (x1^p − x_ℓ)^e x2^j x3^k has weighted degree (i + 8(j+k)) − (p−8)a ≤ D (needs p > 8). Ordered by
x1-degree the change of basis is unitriangular over Z, so {G_ijk} is an O-basis of L_D.

## 4. Lemma 3 (count): correct
‖c‖ ≤ U(1+√2)(1+√|q|) ≤ 8√|q|·U for |q| ≥ 1 (check.py §4); (U+1)^4 > x^4 also when x < 1, so the strict
product bound holds.

## 5. Lemma 4 (residues): correct
G_ijk(points) ∈ p^{⌊i/p⌋}O; k_p ≥ ⌊D/p⌋ ≥ ⌊i/p⌋; |O/p^m O| = p^{4m}; product over (i,j,k) = p^{4k_p r − 4E_p};
CRT over p ∈ S. The bound is independent of the number N⁴ of points, as it must be.

## 6. Lemma 5 (pigeonhole) and (9) ⇒ (11): correct
- #{‖f‖_B ≤ H} > (P/(32r√|q|))^{4r} B^{−4W} with H = P/(4r); beats P^{4r} Π p^{−4E_p} iff (11). ✓
- ‖θ_n‖ ≤ (N−1)(1+√2)(1+√|q|) < B; ‖F(pt)‖ ≤ r·2H = P/2; F(pt) = Pα, α ∈ O; α ≠ 0 ⇒ |Nm| ≥ P⁴ > (P/2)⁴.
  Uses [K:Q] = 4 and O ⊂ O_K. ✓
- Leading terms, exact sums (check.py §1): at D = 10⁵, r/(D³/384) = 1.00027, W/(5D⁴/6144) = 1.00031,
  I/(D⁴/1536) = 1.00034 (I = Σ i); (W − 5I/4)/D³ → −0.0023; |E_p − D⁴/(1536p)| ≤ 0.0011·D³ on sampled p.
- Main terms: Σ E_p log p ≈ I Σ_{p∈S} log p/p, W log B ≈ (5/4) I log B; so (11) needs Σ log p/p > (5/4) log B
  + o(log B). (9) has 13/10 > 5/4: slack (1/20) log B absorbs O(D³Q) = O(D⁴/log²N), O(D³ log B), r log(r√|q|).
  Rigorous numeric slack using E_p ≥ I/p − r and θ(Q) ≤ 1.04 Q (check.py §5): lhs/rhs = 1.030 (N = 10³),
  1.037 (N = 10⁴), 1.039 (N = 10⁶). "N sufficiently large" is mild.
- Q = D/log²N only keeps the O(D³Q) error small; p ≤ Q < D gives k_p = ⌊D/p⌋ ≥ 1.
- Why 4/3 matters: Σ_{p∈S} log p/p ≤ log Q ≈ (4/3) log N, cost (5/4)(log N + ½ log|q|). Need 5/4 < 13/10 < 4/3.

## 7. Lemma 6 versus Masser's published Theorem 1
Masser (1981), p. 82: l = rank Γ; for 1 ≤ r ≤ n, l_r = max rank of a subgroup of Γ lying in a subspace of
dimension n − r; μ(Γ) = min_r (l − l_r)/r. Let γ_1, …, γ_l be basis elements of Γ; Γ(N) = {Σ n_i γ_i : integers
0 ≤ n_i ≤ N}. "Theorem 1. For real N ≥ 0, D ≥ 1 let P be a non-zero polynomial of total degree at most D such
that P(γ) = 0 for all γ in Γ(N). Then we have D ≥ (N/n)^μ."
- Kerr's μ (min over proper L of rank(Γ/(Γ∩L))/(n − dim L)) is the same number (r = codim L). ✓
- Same box (0 ≤ m_j ≤ M), same bound (M/n)^μ. ✓ deg F ≥ 1 is automatic (0 ∈ Γ(M)).
- Difference: Kerr allows any generators γ_1, …, γ_s; Masser assumes a Z-basis. Repair: a maximal Z-independent
  subfamily spans a finite-index Γ' with μ(Γ') = μ(Γ) and Γ'(M) ⊂ the vanishing set. In the application the
  γ_i are a basis (π is injective on Z⁴). Harmless.
- Uniformity: the bound depends only on n and μ, not on Γ. This is essential (Γ depends on q) and it is what
  the published statement gives. Masser remarks the exponent is best possible (c(γ)·N^μ upper bound).
- No erratum found (zbMATH title search). Tao (blog, 19 Oct 2026) endorses this use of Masser; a commenter
  (D. Rui, 21 Oct) gives an independent sumset-to-subspace proof of the needed bound.

## 8. Lemma 7 (μ ≥ 4/3): correct, μ = 4/3 exactly
π injective on Z⁴ (else ker π is a rational line). For proper L ⊂ C³: W = span_C(π^{−1}(L) ∩ Z⁴) is rational and
proper (dim π^{−1}(L) ≤ 3), so ker π ⊄ W, π|_W injective, rank(Γ ∩ L) = dim W ≤ dim L; ratio ≥ (4−h)/(3−h) ≥ 4/3
with equality at L = 0. Application: π(v) = 0 and det(γ_1, γ_2, γ_3) = 4√(2q) verified exactly in Q(a, b)
(check.py §3); Q-independence of (√2q, √2, −√q, −1) ⇔ [K:Q] = 4 ⇒ ker π lies in no proper rational subspace.
Implicit: ((N−1)/3)^μ ≥ ((N−1)/3)^{4/3} needs (N−1)/3 ≥ 1 (true for N large).

## 9. Final contradiction: correct, one unstated finite step
- Σ_{p∈S} log p/p ≥ log Q − C(log|q| + (1−β) log²Q) − O(1); log Q = (4γ/3) log|q| − 2 log log N + O(1);
  log B = (γ + ½) log|q| + O(1). With δ = (1−β) log|q| < ε: need γ/30 > C + 13/20 + (16/9) C γ² ε + o(1).
  Choose γ > 30(C+1), then ε ≲ 1/(C γ²), then |q| ≥ q_0(γ, ε). γ, ε depend only on absolute constants, so q_0
  is absolute.
- Not in the paper: the finite range |q| < q_0. It follows from finiteness: finitely many χ; each L(s, χ) has
  finitely many zeros in [½, 1) and L(1, χ) ≠ 0; zeros β < ½ give (1−β) log|q| > ½ log 3. OpenAI §5 does this
  via a sequence argument. So "absolute c" is consistent; effective in principle, but Kerr's c is not explicit
  (explicit versions of the OpenAI route: Banks δ₀ = 3·10⁻³ for q ≥ 10¹⁰; arXiv 2610.21511: 5.8·10⁻⁴, Lean).

## 10. Davenport–Heilbronn sanity test: passes
Arithmetic input beyond "a Dirichlet series with a real zero":
(a) Lemma 1: Euler product, so −ζ'/ζ − L'/L has coefficients Λ(n)(1+χ(n)) ≥ 0 (χ completely multiplicative,
    values 0, ±1), plus the Hadamard product and functional equation.
(b) Lemma 4: χ(p) is the Legendre symbol (q/p), so Euler's criterion gives the Frobenius action on Z[√q, √2].
(c) Lemma 5 / Lemma 7: integrality and the norm in a quartic field; Q-independence of 1, √q, √2, √(2q).
The real zero enters only through (a), as "χ(p) = −1 on ≥ ~97.5 % of the log-mass up to Q". A DH-type linear
combination of L-functions has neither (a) nor (b); the argument does not apply to it.

## 11. Comparison with OpenAI (1 Oct 2026)
| | OpenAI | Kerr |
|---|---|---|
| auxiliary object | interpolation determinant Δ ≠ 0 (greedy rows, weights 1:H:H) | Siegel-lemma polynomial F (weights 1:8:8) |
| non-vanishing | self-proved uniform interpolation lemma (nearest-point + Lagrange) | Masser 1981, Theorem 1 |
| Frobenius gain | row replacement, p^{E_p} ∣ Δ | basis G_ijk, residue count |
| exponent gap | log U vs (3/4)(1 + S₂/S₁) log U | (4/3) log N vs (13/10) log B, needs 13/10 > 5/4 |
| small q | explicit sequence argument | not stated |
Same field, same points, same analytic lemma. Kerr trades the self-contained interpolation lemma for a
45-year-old published zero estimate; the OpenAI interpolation lemma would also give a Masser-type bound
(with a worse constant), which Kerr's proof could use instead (only log D ≈ (4/3) log N matters).

## 12. Public record (fetched 2026-10-21)
arXiv 2610.21099: v1 only, no comments. Tao, "A sketch of the proof of non-existence of Siegel zeroes"
(19 Oct 2026): lists Banks, Kerr, Green simplifications; praises Kerr's use of Masser; no error raised against
Kerr in the fetched comments. arXiv 2610.21511 builds on a Lean formalization of the OpenAI argument.

## 13. Minor slips (none affects correctness)
1. "K has degree four", Z-basis: false for q = 8 (one modulus; covered by the finite step).
2. Finite range |q| < q_0 not addressed (trivial repair, §9).
3. Lemma 6 stated for generators; Masser assumes a basis (trivial repair, §7; application uses a basis).
4. √(2q) for q < 0 should be defined as √2·√q.
5. Implicit (N−1)/3 ≥ 1 when lowering μ to 4/3.
6. Lemma 1 citation page/number unverified; Milne citation unnecessary.
7. Typos: "there exists distinct f1, f2", "For an integer N ≥ 4" (unused restriction).

Verdict: NO REFUTATION. Every step is correct; the only external input (Masser 1981) is quoted faithfully up to
the generator/basis point.

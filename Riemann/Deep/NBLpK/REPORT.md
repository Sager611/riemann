RESULT: For every number field K and 1 < p < 2: ζ_K(s) ≠ 0 on Re s > 1/p ⟺ κ_K(1/x) ∈ L^p(0,1) and 1_(0,1) lies in the L^p(0,1)-closure of span{λ_𝔟 : 𝔟 ≠ O_K}, λ_𝔟(x) = κ_K(1/(N𝔟·x)) − κ_K(1/x)/N𝔟, κ_K(y) = c_K·y − #{𝔞 : N𝔞 ≤ y} (equivalently with all k ∈ ℕ; ⟺ that span is dense in the step-function space S_p); the Riesz–Möbius means 1 + Σ_{N𝔟≤N} μ_K(𝔟)(1 − N𝔟/N)λ_𝔟 → 0 at rate N^{−η}, η < 1/p − θ, when ζ_K ≠ 0 on Re s > θ; hence unconditionally (Cureton–Thorner's Lean-checked 31/32) for every K and 1 < p ≤ 32/31, and for p ≤ 1/B_2 = 1.14497 conditional on Liu (new)

New to our search; a routine extension of our NBLp theorem with one genuinely new step (Lemma MEM_K, §2).

Node N5 (Nyman–Beurling / Báez-Duarte), number-field analogue. Agent `nflp`, 2026-11-05. Paper proof, no Lean.
Full proof: [PROOF.md](PROOF.md). It adapts [NBLp/PROOF.md](../NBLp/PROOF.md) (K = ℚ) lemma by lemma.

## 1. Statement and why this formulation (PROOF.md §1, §9.1)

* **Functions.** κ_K(y) = c_K y − A_K(y), A_K(y) = #{ideals of norm ≤ y}, c_K = Res_{s=1} ζ_K. This is de Roton's "complementary function"
  Ψ_F for F = ζ_K (JNT 129 (2009), eq. (1)); its Mellin transform is −ζ_K(s)/s. The dilations λ_𝔟 vanish on (1,∞) and are step functions:
  on (1/(m+1), 1/m], λ_𝔟 = A_K(m)/N𝔟 − #{𝔠 ⊆ 𝔟 : N𝔠 ≤ m}, the expected minus the actual number of multiples of 𝔟 (for ℚ: {m/k}).
* **Ideals or norms.** λ_𝔟 depends only on N𝔟, and μ_K summed over ideals of norm k gives b_K(k) (zero off the norm set). So the ideal
  family (smallest; it carries the explicit approximants) and the full natural-number family (de Roton's B_F^nat) give the same criterion;
  the full family also gives the density statement in S_p.
* **Membership (M_p): κ_K(1/x) ∈ L^p(0,1).** For n ≥ 2, κ_K is unbounded (for ℚ(i) it is ¼ of the circle-problem error), so the
  functions may fail to lie in L^p at all; for large n and p near 2 this is an open mean-value question. It holds unconditionally for
  p < (n+1)/(n−1) (Landau's ideal theorem), for all p ≤ 2 when n ≤ 3 (de Roton's mean-square result, quoted), and whenever ζ_K ≠ 0 on
  Re s > 1/p (new Lemma MEM_K). De Roton's p = 2 theorems carry the same condition (Ψ^{(1)}_F ∈ L²(0,∞)).
* **Necessity** (zero ⟹ no approximation) holds for every 1 < p < ∞ under (M_p): x^{w₀−1} kills every λ_k but not 1_(0,1). The Mellin
  identity is first valid only on α_K < Re w < 1 (α_K = 1 − 2/(n+1)) and is continued to Re w > 1/p using (M_p). For p > 2 the
  criterion fails for every K.

Two corrections to the brief: (1) ζ_K is the Hecke L-function of the trivial character of K, so Liu's half-plane applies to it
directly; Brauer induction is not needed. (2) The criterion needs the membership condition (M_p) as part of the statement for n ≥ 2
(it is implied by the zero-free half-plane, so the corollaries are unaffected).

## 2. Proof summary (what changes against NBLp)

| NBLp step | For ζ_K |
|---|---|
| Structure, D_m-invariance, density in S_p | values A_K(m)/k − A_K(⌊m/k⌋); D_m preserves the ideal span only for m in the norm set, so density in S_p is stated for the full family |
| Lemma T (Mellin of the dilations) | only on α_K < σ < 1 (Landau) |
| **Lemma MEM_K (new)** | under a zero-free half-plane, the Mellin identity moves to every line a > θ (resp. a ≥ 1/p − ε/2) by Paley–Wiener line shifting, using Lindelöf in the zero-free region and the Hadamard midline bound; gives (M_p) and puts all lines of the argument in reach, for every degree |
| Lemma Γ (gamma ratio) | degree n: \|γ_K(s+ε)/γ_K(s)\| ≤ e^{c₁ε}(2+\|s\|)^{nε/2}, c₁ = 3 + ½log\|d_K\| + 2n; forces ε ≤ 1/(2n) |
| Lemma H (Hadamard midline \|ζ_K(s)\| ≤ \|ζ_K(s+ε)\|·\|γ_K ratio\|) | Hadamard product of ξ_K (Re B = −Σ Re 1/ρ) |
| Lemma A (Lindelöf in a zero-free half-plane), Lemma P (Riesz–Perron) | convexity input of degree n/2 (Rademacher); μ → μ_K |
| A4 (key estimate ‖Φ_ε‖² ≤ Kε via the exact residue integral) | pole term c_K/ζ_K(1+ε) ≤ C(K)ε (ζ_K(1+ε) ≥ c_K/ε can fail when the Euler–Kronecker constant is negative); R(2σ') ≥ 1 still holds (ζ_K decreasing) |
| A5, Theorem B | correction ‖κ_K(1/x)‖_p·\|P_N(1)\| needs (M_p) |
| Order of limits | unchanged: N → ∞ at fixed ε, then ε → 0; nothing is used on the line Re s = 1/p |

Degree and discriminant enter only through ε₀ = min{½, (2/p−1)/2, 1−1/p, 1/(2n)} and the constants; they never restrict p.

## 3. Self-review

* **Pass 1 (own, targeted):** necessity (boundedness of the functional, continuation of the Mellin identity, removable point w = 1);
  the A4 residue bookkeeping with the pole c_K; every ε₀ constraint; Lemma PW hypotheses in MEM_K; the A1 identities.
* **Pass 2 (4 parallel adversarial reviewers, one per block):** Lemma T_K/MEM_K/Theorem B_K: SOUND (they confirmed that no σ > α_K
  requirement survives); Lemmas Γ_K, H_K, A_K, P_K: SOUND (K = ℚ constants e^{5ε}, e^{8ε} recovered); A1–A7 and corollaries: SOUND
  (orientation, both residues, decay on the horizontal sides, order of limits); structure and necessity: SOUND with one FIXABLE GAP,
  fixed: §0 allowed real k but Prop. C_K(a) needs integer k (counterexample k = 3/2), now restricted to integers (real dilations only in
  Remark 9.1). Also fixed: Θ_K ≤ 1 only (≤ 31/32 by Cureton–Thorner), not < 1 from the classical theory.
* **Residue identity, numerically** (jcheckK.py, mpmath; Simpson to |t| = 100, step 0.0005 on |t| ≤ 2 and 0.025 beyond, plus the
  tail 2/T):

| K | p | ε | J numeric | J formula | difference | sign of the s = 2σ' residue flipped | pole term dropped |
|---|---|---|---|---|---|---|---|
| ℚ(i) | 8/7 | 0.1 | 3.25092993 | 3.2509308 | −8.8·10⁻⁷ | −5.035 | 4.143 |
| ℚ(√5) | 8/7 | 0.1 | 3.19141856 | 3.19142129 | −2.7·10⁻⁶ | −4.881 | 4.036 |
| ℚ(i) | 32/31 | 0.03 | 3.14592323 | 3.14592362 | −3.9·10⁻⁷ | −3.551 | 3.349 |
| ℚ(√5) | 32/31 | 0.03 | 3.12812606 | 3.12812653 | −4.7·10⁻⁷ | −3.527 | 3.327 |

  A first run with the uniform step 0.025 gave +2.0·10⁻⁴ in both p = 32/31 rows; extending to |t| = 400 left it unchanged
  (2.01·10⁻⁴); the cause was the quadrature near t = 0, where the pole of ζ_K at s = 1 lies only 1 − σ' = 0.046 from the line. The p = 8/7
  rows used the uniform step (pole at distance 0.175); the fine step reproduces ℚ(i) (−8.74·10⁻⁷).
  Numerically B(ε)/ε ≈ 0.39 (ℚ(i)), 0.42 (ℚ(√5)) at p = 8/7 and 0.349, 0.353 at p = 32/31, far below the proven K_{p,K}
  (C_{Γ,K} = 5374 for ℚ(i), 5760 for ℚ(√5), 5993 for ℚ): the constants are qualitative, as in NBLp.

## 4. Numerics (nbk.c, exact 64-bit integer sieve, M = 10⁸ intervals, one process, 2 min; log run1e8.txt; tails nbk_tail.py)

‖1 + Σ_{N𝔟≤N} μ_K(𝔟)(1 − N𝔟/N) λ_𝔟‖_{L^p(0,1)} = (Σ_m |v(m)|^p/(m(m+1)))^{1/p}, v(m) = 1 + P_N(1)A_K(m) − Σ_{n≤m} Σ_{k|n,k≤N} b_K(k)(1 − k/N)a_K(n/k),
a_K = 1 * χ_D, b_K multiplicative (b(p) = −1 − χ(p), b(p²) = χ(p)). Brackets: rigorous [partial sum to 10⁸, partial sum + tail bound],
the tail bound from |A_K(u) − c_K u| ≤ C₁√u + C₀ (ℚ(i): C₁ = π√2/4 by lattice squares; ℚ(√5): C₁ = 5 by the hyperbola method).
Unbracketed: partial sum plus a geometric continuation of the last decade increment (decade ratios 0.18–0.19 at p ≤ 1.1, matching
|v(m)| ≈ m^{1/4}); the rigorous bound is too weak there.

| K | p | N = 10 | N = 100 | N = 1000 |
|---|---|---|---|---|
| ℚ(i) | 32/31 | [0.52709, 0.52762] | [0.11739, 0.12010] | 0.031555 (≤ 0.0406) |
| ℚ(i) | 1.1 | [0.53057, 0.53170] | [0.12336, 0.13047] | 0.037303 (≤ 0.0644) |
| ℚ(i) | 1.5 | 0.55362 | 0.17798 | 0.10857 |
| ℚ(i) | 2 | 0.58992 | 0.30552 | 0.31659 |
| ℚ(√5) | 32/31 | [0.83214, 0.83333] | [0.27276, 0.28077] | 0.063437 (≤ 0.0920) |
| ℚ(√5) | 1.1 | [0.83385, 0.83646] | [0.27879, 0.30016] | 0.068471 (≤ 0.1558) |
| ℚ(√5) | 1.5 | 0.84438 | 0.32455 | 0.12429 |
| ℚ(√5) | 2 | 0.85874 | 0.41438 | 0.29440 |
| ℚ (NBLp, M = 10⁹) | 1.1 | 0.34595 | 0.05923 | 0.01505 |

* Regression: with K = ℚ the code reproduces NBLp's partial sums digit for digit (M = 10⁶, 10⁷).
* At p = 32/31 the decrease over each decade is rigorous (upper bound at N = 10ᵏ⁺¹ below the lower bound at 10ᵏ). Local exponents at
  p = 1.1: 0.63, 0.52 (ℚ(i)), 0.48, 0.61 (ℚ(√5)), consistent with the GRH-conditional rate N^{−η}, η < 1/p − 1/2 = 0.41 (Theorem B_K;
  both fields satisfy GRH numerically far beyond the heights that matter). At p = 1.5 the decay slows; at p = 2 the norms stall
  (ℚ(i) even rises from N = 100 to 1000), as for ζ and Dirichlet L-functions (logarithmic rate at p = 2).
* Max |A_K(m) − c_K m|/√m for m ≤ 10⁸: 0.48 (ℚ(i)), 0.57 (ℚ(√5)); max |v(m)| grows like m^{0.3}: the functions are unbounded,
  unlike K = ℚ, which is why (M_p) is part of the statement.

## 5. Novelty pass

Question: is an L^p (p ≠ 2) Nyman–Beurling or Báez-Duarte criterion known for Dedekind zeta functions (or the Selberg class)? Not found.
* **de Roton** (read in full text from HAL): Trans. AMS 359 (2007) Thm III (Beurling–Nyman, L²(0,1), Selberg class, with the
  condition Ψ^{(1)}_F ∈ L²(0,∞), known for degree < 4) and JNT 129 (2009) Thm I/Cor. 1 (natural-number Báez-Duarte form, L²(0,∞),
  rate (log log N)^{−1/3}): p = 2 only. Her introduction cites Báez-Duarte and Balazard–Saias for the L^r (r < p) ζ-results only.
* **Dimitrov–Oliveira** (1608.07887, TeX read): their L^p theorem (Thm GNB, continuous dilations, all p > 1) assumes κ bounded (class
  𝒟, TeX l. 317–323), which excludes ζ_K for n ≥ 2; natural numbers only at p = 2. So even the continuous (Beurling-type) L^p criterion
  for ζ_K with n ≥ 2 was not covered; Theorem A_K gives it too (Remark 9.1).
* **Delaunay–Fricain–Mosaki–Robert** (1101.1199, 1112.0166): Hilbert-space (weighted L²) criteria and zero-free discs for general
  Dirichlet series.
* Citing works (Semantic Scholar): de Roton 2009: 6, de Roton 2007: 11; none treats L^p or number fields beyond p = 2.
* arXiv listings math.NT 2026-11 (140) and 2026-10 (854), math.CA and math.FA for both months: no Nyman–Beurling or Báez-Duarte title
  (keyword scan of all titles). arXiv API (abs: Nyman ∧ Dedekind, Beurling ∧ "Dedekind zeta", Nyman-Beurling ∧ L-functions,
  Baez-Duarte, Beurling-Nyman): nothing on number fields or L^p. OpenAlex (3 queries): only the works above.
* Zenodo (5 queries, incl. most recent Nyman–Beurling records): 6982337, 17957671, 22647730, 22650745, 23129278, 23370435 are ζ-only;
  none mentions Dedekind, number fields or L^p.
* Our files: DISCOVERIES (nyman-beurling-lp: ζ and Dirichlet L), COROLLARIES l. 10, 13 (continuous L^p for ζ), R350 (this question,
  listed as open). Evidence API: quota exhausted (HTTP 429), not used.

Assessment: new to our search. It is a routine extension of our NBLp theorem (like NBLp2), with one genuinely new step: for n ≥ 2 the
functions are unbounded and, for large degree, not known to lie in the weighted L² spaces on the lines the argument needs; Lemma MEM_K
shows that the zero-free hypothesis itself supplies this, so the endpoint theorem and the unconditional p ≤ 32/31 statement hold for
every number field. The unconditional corollary inherits the status of Cureton–Thorner's 31/32 (Lean-checked, unrefereed). Best home:
a generalization line on the nyman-beurling-lp discovery page or COROLLARIES; the orchestrator decides.

## 6. Limits

The theorem moves no RH node: for p < 2 the criterion is equivalent to the zero-free half-plane Re s > 1/p, so it is no stronger than
the known half-planes and no weaker than GRH for ζ_K as p → 2. Constants depend on K (no uniformity). p = 2 is not covered (de Roton).
Membership (M_p) for p near 2 and large degree remains a mean-value question that the theorem does not settle (it follows from the
zero-free half-plane).

## 7. Files

* [PROOF.md](PROOF.md): full proof, lemma-by-lemma status against NBLp, K-dependence table (§9.2).
* [nbk.c](nbk.c) (cc -O2 nbk.c -lm; ./nbk M type N…, type 0 = ℚ, 1 = ℚ(i), 2 = ℚ(√5)); [run1e8.txt](run1e8.txt); [nbk_tail.py](nbk_tail.py) (tails).
* [jcheckK.py](jcheckK.py) (uv run --no-project --with mpmath python3 jcheckK.py); [jcheckK.txt](jcheckK.txt).

RH relevance: N5 (Nyman-Beurling for number fields)

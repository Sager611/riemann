DISCOVERY: Over every field 𝕜 of characteristic 0, for every finite connected graph G (≥ 2 vertices), A(G) ≅ A_RW(G) iff G is regular, or G is bipartite biregular with part degrees d₁, d₂ such that d₁/d₂ is a cube in 𝕜. So the open problem of arXiv 2610.22248 holds over R and C (in characteristic 0 exactly over the fields with cube roots of all primes) and fails over Q (P₃, K_{1,3}, K_{2,3}, …). The 2021 CRR conjecture holds over R; over C its isomorphism part holds but its "moreover" clause fails for every singular G that is neither regular nor biregular (e.g. the diamond). Pending novelty check.

# EvoAlg report (2026-10-22)
Proofs: [NOTE.md](NOTE.md). Checks: [verify.py](verify.py) (25 s; `--full` 6–8 min; all pass, < 150 MB). Brief: `.scratch/orca/evoalg-brief.md`. Paper proofs, not refereed, no Lean.

## Results (complete proofs in NOTE §2–7)
- **Theorem N (task 1).** In a non-degenerate evolution algebra over any field, the natural bases are exactly the unions of b_T-orthogonal bases of the blocks E_T (classes of proportional squares, which are the twin classes for A(G) and A_RW(G)). Key input: BCS22 Theorem 2.4(i) "⇒" (Publ. Mat. 2022, arXiv 2006.14460), quoted exactly and reproved in one line.
- **Theorem I.** The isomorphisms A(G) → A_RW(G) are exactly a size-preserving automorphism τ of the twin quotient plus blocks P_S with P_SᵀP_S = c_S d_S I and P_S𝟙 = c_T𝟙 (T ∼ S). Holds over any field.
- **Lemma 2 (the new step).** k(γ² − α) = 𝟙ᵀ(PᵀP − αI)𝟙 − (P𝟙 − γ𝟙)ᵀ(P𝟙 + γ𝟙), so the all-ones vector forces γ² = α: rotating twins never changes the scalars. What remains is the twin-free scaling system, which gives Theorem A (NOTE §5; the isomorphism is CRR's diagonal map with 𝕜-cube roots).
- **Characteristic p (task 2).** In characteristic 2, A_RW(G) = A(G) whenever it is defined. For odd p the criterion uses degrees mod p, plus similitude multipliers on twin classes whose size is divisible by p. "Iff regular or biregular" fails in every characteristic p (double star S(p,p)). K_{1,1,3} and a 6-vertex tree are isomorphic over F₉ only by twin mixing, P = [[0,−1,−1],[1,i,−i],[1,−i,i]], and are not isomorphic over F₃.
- **CRR's "moreover" clause** (non-isomorphic ⇒ only the zero homomorphism). It holds automatically over formally real fields, so CRR's whole conjecture holds over R. It fails over C for every singular G that is neither regular nor biregular, e.g. the diamond.
- **By-product.** For the survey's question on natural bases, Theorem N gives the description. Over infinite fields there are finitely many natural bases iff (2LI). A(C₄) has infinitely many natural bases but |Aut| = 8 over Q and R.

## Checks (task 3)
- All 995 connected graphs with 2–7 vertices: the twin-class system solved over C, R and Q agrees with the criterion (24 isomorphic over C and R, 15 over Q). Every predicted isomorphism was checked exactly.
- Full GL_n Gröbner bases over Q: P₃ gives p₀₀³ = 4, K_{1,3} gives p₀₀³ = 9, the diamond gives [1] (and in --full, P₄ and the paw give [1]).
- Brute force with no theory: over F₃, F₅ and F₇ (n ≤ 5), natural-basis counts equal Theorem N and isomorphism decisions equal Theorem P; K_{1,1,3} and T6 checked over F₃ and F₉.

## Novelty (task 4)
- The question is open in 2610.22248 (15 Oct 2026). CRR (1707.05897, 1710.10516) settle only the non-singular case over R; CRRV (2405.12341) gives examples. No solution was found in 141 "evolution algebra" titles, 8 further arXiv queries, or the math.RA and math.CO listings for Sep–Oct 2026, and none mentions the obstruction over Q. No authors were contacted.

## Side findings (NOTE §8, verified)
- BCS22: Theorem 2.4(i) "⇐", Proposition 2.16 and Corollary 2.17 fail in characteristic 2, and the survey repeats Corollary 2.17 over every field. Corollary 2.7 needs |𝕜| ≥ 4. In Theorem 2.11, "r = dim A²" fails for A(P₅).
- CRRV: Corollary coro:column0 holds over R only.
- 2111.07399: the proposition "twin construction" is false.

## Not done
- Infinite graphs and the Hilbert-space version.

## Review (2026-10-22)
> INCOMPLETE / GAP — source verification only; the self-contained mathematical proofs are CORRECT.
> Theorem N’s rank-one multiplication argument works over every field and independently proves the needed BCS22 implication.
> Theorem I correctly derives the twin-quotient permutation, preserved sizes/degrees, and block equations.
> Lemma 2 and the connectedness argument correctly yield Theorem A’s precise cube-root condition.
> The characteristic-p criterion and the complex counterexamples to CRR’s “moreover” clause are valid.
> `verify.py`: exit 0; all checks passed.
> Significance: establishes the stated characteristic-zero classification, including singular graphs; independent novelty remains unverified.
> `REPORT.md:1`: over C, only the isomorphism part of CRR’s conjecture holds—not its “moreover” clause; “DISCOVERY” asserts unverified priority.
> TODO @scout: authenticate the BCS22 quotation against the original; the permitted NOTE/REPORT-only review cannot establish quotation fidelity.
- Fixed: line 1 (CRR's "moreover" clause over C; priority pending the novelty check). The BCS22 quotation is checked separately; Theorem N's own argument proves the implication used.

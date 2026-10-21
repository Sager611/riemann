NO DISCOVERY: cubic (S₃) Chebotarev from our ring class results is weaker than published almost-all bounds (Ellenberg–Pierce–Wood 2017, Pierce–Turnage-Butterbaugh–Wood 2020, Lemke Oliver–Thorner–Zaman 2024); the strongest version our tools give is routine.

Date 2026-10-21. Details and proofs: [NOTE.md](./NOTE.md). Paper proofs, not refereed, no Lean, no 2026 half-planes.

## Short answer
- Cubic fields need only 3 of the h(D) class group characters (trivial, ψ, ψ̄) plus the quadratic character χ_D.
  Entry 2's family large sieve (all h(D) classes at once) plays no role; a zero-free box near Re s = 1 suffices.
- Entry 2 does give a corollary (Theorem B): for A > 1.1307, all but O(X^{1−c}) fundamental D, every cubic field of
  discriminant D has primes of each S₃ class in (x, x+x|D|^{−κ}], x ≥ |D|^A. But c ≈ 10^{−4} < 1/3, so the cubic fields over the
  exceptional D are only bounded by X^{1−c}·X^{1/3+ε} (Ellenberg–Venkatesh): nothing about almost all *fields*.
- Strongest version with our tools (Theorem A): any A > 0, all but O(X^ε) fundamental D, asymptotic count of each S₃ class in
  (x, x+x|D|^{−κ}]; all but O(X^{1/3+ε}) cubic fields with fundamental discriminant (73% of all cubic fields). Inputs: our family
  density near Re s = 1 ([R] Prop 1, [R2] Prop 2) for the cubic characters, Thorner–Zaman's log-free density for χ_D, Ingham/VK for ζ.
- Non-fundamental d_K (27%): not covered by entry 2 or [R] Prop 1; a sketch via Humphries–Thorner gives O(X^{1/2+ε}) fields, not claimed.

## Why not new (sources read, arXiv API 8 queries, GXL 5 claims)
- Ellenberg–Pierce–Wood, ANT 11 (2017), Thm 2.1 + Rem 2.1: all cubic fields, all but O(X^{1−δ+ε}) (δ ≤ 2/25) have ≫ X^δ/log X
  primes ≤ X^δ of each splitting type. Exponent 0.08 against our 1.13, and counted in fields.
- Pierce–Turnage-Butterbaugh–Wood, Invent. 219 (2020), Thm 3.3(2): cubic fields with squarefree tame discriminant part (contains
  all fundamental ones), all but O(X^{1/3+ε}): asymptotic Chebotarev for x ≥ exp(κ(log log D)^{5/3+o(1)}). Same 1/3 as Theorem A.
- Lemke Oliver–Thorner–Zaman, Invent. 235 (2024), Thm 1.1, Cor 3.8: all cubic fields, all but O(X^ε), the 2-dim Artin L-function
  is zero-free in a strip of fixed width up to height D^{O(1)}; with classical facts for χ_D this gives Theorem A, short intervals included.
- Lemke Oliver–Smith (arXiv 2405.08383): every G, but fields over a bad quadratic resolvent are not controlled (their Rem 1.11).
- Every cubic field: 3-cycles and transpositions ≪ |d_K|^{3/4+ε} (Cho–Lemke Oliver–Zaman, arXiv 2512.24963, Dec 2025), identity
  ≪ |d_K|^{3/2+ε} (Ge–Milinovich–Pollack 2019). So for 3-cycles our almost-all 1.13 is worse than the every-field bound.
- Conditional on the unrefereed 2026 half-planes: every field, least prime ≤ (C log D)^{7.99} (Frei–Sofos–Widmer 2610.18491; ArtinHalfPlane).
- Short intervals for almost all cubic fields: not found stated anywhere, but routine from the LOTZ strip. At most a COROLLARIES entry.

## Checks
- Dictionary (`.scratch/cubiccheb/dict.gp`): 14 fundamental D (3-rank 1 and 2, both signs), 26 cubic fields, primes 3 ≤ p ≤ 3000, p ∤ D·disc:
  splitting type in K equals the prediction from (D/p) and ψ(𝔭); 11,135 tests, 0 mismatches; nfdisc(K) = D always.
- Least primes (`least.gp`): 2716 D with 3 | h(D), |D| ≤ 20000, 5564 fields; worst log P/log|D|: identity 1.30 (D = −23; median
  0.31), 3-cycles 0.39, transpositions 0.52. Consistent with polylog least primes.

## Dead ends (NOGO candidates)
- Cubic fields from the every-class theorem: exceptional D may carry large 3-torsion; needs c > 1/3. Use the 3 characters instead.
- Fewer than X^{1/3} exceptional fields at small A: χ_D is shared by all cubic fields of discriminant D; needs moments of |Cl(D)[3]|.

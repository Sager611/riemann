NEW: no source found proves or states, for singular graphs, any criterion for A(G) ≅ A_RW(G). That includes Theorem A (char 0: regular, or biregular with d₁/d₂ ∈ (𝕜^×)³), the failure over Q (P₃, K_{1,3}, …), Theorem P (char p) and Proposition H. The CRR conjecture (over R, 1710.10516 l. 963) is still called open by CRRV (J. Algebra Appl., online 9 May 2026) and by the survey 2610.22248 (15 Oct 2026). No source treats the question over a field other than R or C.

# EvoAlg novelty check (adversarial, 2026-10-22)
Scratch `.scratch/evonov/` was deleted after the run. Claims are as in [REPORT](REPORT.md) and [NOTE](NOTE.md).

## Verdicts, main results
- **Theorem A / Corollary A: NEW.** Closest work: CRR21, Theorem 2.3 (non-singular graphs, over R); CRR20 (regular graphs and K_{m,n}, with real cube roots); CRRV (caterpillars, tadpoles T_{4,m}, the bull graph). Cadavid–Rodriguez, Proc. CNMAC 7 (2020) 010423, and CRRV l. 237 only *suggest* using the twin partition.
- **Theorem I and Lemma 2: NEW as stated.** The block form of isomorphisms and base changes is known: BCS22 Cor. 2.13, and Calderón–Fernández Ouaridi–Kaygorodov (JPAA 2024, arXiv 2303.00862), Lemma "princi" and Remark "baseshape".
- **Theorem P: NEW.** No source treats A_RW(G) in positive characteristic. Reis–Cadavid (2006.11934) cover only derivations of A(G).
- **Proposition H: NEW.** Its positivity argument is the one in CRRV coro:column0. 2111.07399 proves injectivity only for non-singular graphs, and only in the Hilbert setting.
- **Theorem N: PARTLY KNOWN.** Known: the blocks, their uniqueness and block-diagonal base changes (BCS22 Thm 2.11 + Cor. 2.13; CFK24 Remark baseshape; ACV 2607.15154, proof of the balance theorem). Known in char ≠ 2: "a block's natural bases are exactly its b_T-orthogonal bases" (the proof of BCS22 Prop. 2.16; for dim A² = 1, Brache–Martín Barquero–Martín González–Sánchez-Ortega, 2103.01625, Mediterr. J. Math. 2025; the survey's Prop. prop:On is the one-block case). New: only the all-characteristics statement and the A(C₄) example (infinitely many natural bases, finite Aut).

## Verdicts, side findings
- **BCS22 Thm 2.4(i) "⇐" fails in char 2: KNOWN.** Cabrera Casado, Cardoso Gonçalves, Gonçalves, Martín Barquero, Martín González, Ruiz Campos, "Solenoids in automorphism groups of evolution algebras", arXiv 2505.14584 (May 2025), §5, Remark after the natural-vector proposition. Their example: K³ in char 2 with e₁² = e₂² = e₃² ≠ 0, where u = e₁+e₂+e₃ is not natural. They conclude that "char(K) ≠ 2 must be added". Their example also refutes Prop. 2.16 with |C| = 1, but they do not say so.
- **BCS22 Cor. 2.17 fails in char 2: NEW.** It is not noted in 2505.14584 or in any other source found. The survey (l. 332–333) restates it for every field, although it cites 2505.14584 (l. 590).
- **BCS22 Cor. 2.7 needs |𝕜| ≥ 4: KNOWN.** 2505.14584 §5 says the cardinality hypothesis "is not present in the statement … however, it is used in its proof". Our F₂ and F₃ counterexamples, which show the hypothesis is necessary, are not given there.
- **BCS22 Thm 2.11 "r = dim A²": KNOWN, already fixed in print.** The phrase is only in arXiv v1 (TeX l. 493). The published statement (Publ. Mat. 66 (2022), p. 167; UAB DDD PDF) has no "r = dim A²". NOTE §8.3 and REPORT should say "arXiv v1".
- **CRRV coro:column0 holds over R only: NEW.** CRRV work over "a field 𝕂" (l. 94, 224), but the proof needs a formally real field. Crossref lists no erratum or update, and no later work (2505.15604, the survey) mentions it.
- **2111.07399, "twin construction" proposition false: NEW.** arXiv has only v1. Crossref lists no erratum, and no citing work mentions it. CRRV (same group) restate the "only the null map" conjecture, which the proposition would contradict.
- By-product: Cadavid–Rodriguez, CNMAC 2020, Example 2.2 states "only the null map" for all caterpillars C_{a₂…a_{n−1}} with a_i ≥ 0. For n = 3 these are the stars K_{1,m}, which are biregular, so A ≅ A_RW over R (CRR20). Presumably n ≥ 4 was meant.

## What was searched
- **Citing works, with abstracts.** Semantic Scholar: 12 seeds (CRR20, CRR21, BCS22, CRRV, 2111.07399 and their DOIs, EL15, EL21, CCRR21, CFNT22, the survey), 104 citing works. OpenAlex: 11 seeds, 85 works (the CFNT22 and survey lists returned 429). 112 distinct works in total, keyword-scanned.
- **arXiv math.RA, 2024-01 to 2026-10 (up to 21 Oct).** Listing pages 2024-01…08 gave 857 ids; "Rate exceeded" came twice at 2024-09, so the listings were stopped there. API month queries (cat:math.RA) then gave 4,238 records with abstracts; 6 listing-only ids came from one id_list call. There are 20 "evolution" titles; all abstracts were scanned. All-time API queries gave 103 "evolution algebra(s)" records, plus 148 more from random walk/graph queries.
- **Full text read.** The survey 2610.22248 and the survey it cites, Ceballos–Falcón–Núñez-Valdés–Tenorio (Expo. Math. 2022 = arXiv 2604.13126, also Zenodo 17748510); Kaygorodov 2306.00425 (nothing relevant). BCS22 (arXiv TeX and published PDF); CRR20, CRR21, CRRV, 2111.07399; 2505.14584, 2505.15604, 2303.00862, 2607.15154, 2012.06533, 2103.01625, 2401.05924, 2609.32784, 2311.08387, 2605.12672, 2404.08752, 1312.4685, 1902.02191, 2012.05124, 2104.02470, 2503.07187; the CNMAC 2020 PDF.
- **Other sources.**
  - zbMATH Open: 156 "evolution algebra(s)" documents plus 65 with MSC 17D92 (2023–26); 32 cite CRR20, CRR21, BCS22 or 2111.07399 (Sib. Math. J. version); reviews scanned.
  - Crossref: no "updates" (errata) for BCS22, CRR20, CRR21, CRRV or 2111.07399 (Sib. Math. J. version). Posted content 2023–26: 3 unrelated Preprints.org items.
  - Zenodo: 14 queries, nothing relevant. HAL: 3 records. Hexagon: all 381 works. openai/math CONTENTS.md. Nothing relevant in any of them.
  - UdeA repository: Henao Acevedo's 2020 thesis (abstract only).

## Residual risk
- Not seen: the published versions of 2111.07399 (Sib. Math. J. 63 (2022) 995–1011, paywalled; the zbMATH review does not mention the proposition) and of CRRV (J. Algebra Appl.).
- Not seen: Henao Acevedo's unpublished 2022 master's thesis (cited by CRRV).
- The source of 2311.00102 was unavailable, so only its abstract was read.
- Semantic Scholar keyword search was rate-limited (3 of 6 queries ran). Google Scholar was not used.

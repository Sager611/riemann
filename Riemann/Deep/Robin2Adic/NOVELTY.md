# NOVELTY — Robin2Adic (v_2(n) ≤ K unconditionally), scan of 12 Nov 2026, 18:00–19:12 UTC

Verdict: no unconditional v_2(n) result beyond Axler's v_2(n) ≤ 20 (published) or teal-sea's ν_2(n) ≤ 24 (unreviewed) found. Our K = 72 (KLN) and the conditional K = 75 (Chourasiya–Simonič) are new as statements. The method is Gracia Llorente's (Zenodo 24044336); our part is the inputs (our Robin range, the Lean half-plane, KLN density with RH(3e12)) and an Arb certificate. A reviewer may call it a routine combination.

## Prior art (v_2 and t-free Robin)
| result | status | source |
|---|---|---|
| v_2(n) ≤ 16 | published | Wu (arXiv 1812.00987) |
| v_2(n) ≤ 19 | published | Hertlein, Integers 18 (2018) |
| v_2(n) ≤ 20; every 21-free n | published | Axler, Ramanujan J. 61 (2023) 909–919 (arXiv 2110.13478, Thm 1.2/1.4) |
| ν_2(n) ≤ 24, 25-free, valuation list (Cor. 2) | unreviewed, "awaiting outside review" | zeta.teal-sea.com/read/hunts/robin_tfree/RESULTS.html (10 Oct 2026) |
| v_2(n) ≤ 73 | conditional (Re s > 7/8 zero-free; CS Table 1 with b_4 = 0.54) | Gracia Llorente, Zenodo 24044336 (7 Nov 2026), authenticated in our R427 |

## Searches (scripts .scratch/robin2adic/nov.py, oa.py; deleted after use; outputs summarized here)
- Zenodo 24044336: one version only (concept 24044335; zver.py).
- arXiv monthly listings math.NT Oct–Nov 2026, title /Robin/ (lst.py): 2611.13606 (threshold-crossing condition for the least CA counterexample), 2610.12816 (reciprocal LCM powers and Robin's criterion), 2611.02058 (unrelated). None states a valuation range.
- arXiv API since 1 Sep 2026: abs:Robin AND abs:free (13, all off-topic except 2609.26787, a finite arithmetic equivalent of Robin); abs:Robin AND abs:valuation (2611.13606 only); abs:"Robin inequality" / "Robin's inequality" (2611.13606, 2610.12816, 2609.33794 higher-order CA numbers, 2609.26787).
- Zenodo: "74-free" → only 24044336; "Robin inequality" AND valuation → 0; "Robin's inequality" → 31 records, newest 24122732 (5130 circuit, CA margins), 23765626 (expository), 22014299, 21725356 (valuation bound for superabundant numbers only), 20931450 (Kaneko–Lagarias), 20733714 (least counterexample); "Robin" AND "free" → 608 records, none on Robin's inequality in Oct–Nov; Gracia Llorente's 7 records: only 24044336 is on Robin.
- OpenAlex: Axler's "On Robin's inequality" (W3209842564, doi 10.1007/s11139-022-00683-0) has 4 citing works (2023–2024: arithmetic-function inequalities, three univalent-function papers); none extends v_2.
- Our files: TREE.md N4 and R427 (GL authenticated; robin2adic open); no v_2 result of ours. REFUTATIONS.md: CS's b_4 = 0.54.

## What is new here
1. Unconditional (no RH, no CS): Robin for every n > 5040 with v_2(n) ≤ 72 (δ from KLN Table 1; 71 with δ = 1), every 73-free n in particular.
2. Conditional on CS 2507.15184 (b_4 corrected to 0.5501): v_2(n) ≤ 75.
3. Slot form: a counterexample has p_j^{a_j+1} > 1.14·10^22 in every slot (KLN), e.g. 2^73 | n, and a_j ≥ 2 for every prime-factor slot with p_j ≤ 1.07·10^11 (CS: 2^76 | n).

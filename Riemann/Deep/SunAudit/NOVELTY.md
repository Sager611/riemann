NEW: (1) 1208.2683 Conj. 2.1; (2) 1208.2683 Conj. 3.1 (Lucas part); (3) 1211.1588 Conj. 4.30(iii); (4) 1809.07766 Conj. 6.8 (6.15) at n = 7 (typo level); (5) 1402.6641 Conj. 2.25 (wrong formula in Remark 2.25); (6) the further 1811.10503 Conj. 4.9 instances F_7, F_11, F_13 only ; KNOWN: (7) 2607.07638 Domb conjecture at p = 3 (K. Calderon, arXiv 2610.04772v1, 3 Oct 2026, Remark "rem:p3"); (6) base case n = 2, 3-term progressions (the-omega-institute/trureturing issue #16084, 13-14 Oct 2026, Lean refutation; earlier than our page of 22 Oct) ; CORRECTED-IN-PRINT: none

# Adversarial novelty check of the sunaudit counterexamples (sunnov4, 2026-11-05)

Verification: `python3 Riemann/Deep/SunAudit/verify.py` gives ALL PASS (1.1 s, 17 MB). `python3 Riemann/Deep/SunCx/verify_2607.py` passes (3.3 s, 16 MB). `verify.py` has no row for the F_7 instance of item (6); I checked it by hand: in F_7 the sums a+2b over distinct a, b in {0,1,2,4} give {1,...,6}, so 6 < min{7, 8} = 7.

"Print" below means Sun's own homepage copy of the journal version (maths.nju.edu.cn/~zwsun, PDF with the journal header). Every printed statement is identical to the final arXiv TeX, so nothing was corrected in print.

## Action needed (outside this file)
- discoveries/sun-counterexamples.md, row "arXiv 1811.10503, Conj. 4.9": the case n = 2 is not new. the-omega-institute/trureturing preregistered the refutation n = 2, A = {0,1,2} over Q and in every allowed field (all 3-term progressions, with a proof in prose) on 2026-10-13 (issue #16084). It merged a Lean proof that the claim is false (PR #16566, file `Problems/sun-2021-conjecture-4-9-weighted-distinct-sums.md`) on 2026-10-14. Our first record is R151 on 2026-10-22. The same file also notes that the group version in Remark 4.9 fails for Z/p.
- Same page, entry "excluded as not new" for 1309.1679 Conj. 3.8: the earliest public refutation at (11, 2) is trureturing issues #14536 and #14581 (2026-10-08), four days before the Zenodo preprint the page cites.

## Sources searched for every item (2026-11-05)
- Semantic Scholar citations with citing sentences. Citer counts: 1208.2683: 49, 1211.1588: 11, 1809.07766: 43, 1402.6641: 9, 1811.10503: 2, 2607.07638: 1.
- OpenAlex works citing each paper (34, 6, 26, 4, 5, 0) and zbMATH `ci:` (Zbl 1308.11031: 6 citers, Zbl 1458.11004: 3, others 0).
- The GXL Evidence API (10 claims: no supporting or refuting passage outside Sun's own papers, except Calderon for (7)).
- OEIS: JSON searches, the entries named below, and revision histories of A261583 and A002895.
- MathOverflow through the Stack Exchange API: Sun's questions by topic, a global search, and questions 497120, 368072, 316142, 316143, 315702 and 337896 with their answers and comments.
- The Number Theory List archive (NMBRTHRY): subject lists for Aug-Dec 2012, Feb-Mar 2014, Aug-Sep 2015, Sep-Dec 2018, Jul 2025 and Jul 2026.
- Zenodo API (about 23 queries).
- GitHub code and issue search (about 50 queries). This covers trureturing, openai/math (`CONTENTS.md` and `overview.tex`: no Sun), atlas-lean, formal-conjectures, and the other AI-math repositories that the search returned.
- The arXiv API, plus the math.NT and math.CO listings for Oct and Nov 2026.
- Sun's homepage: recent publications and 11 talk slide sets (Sequence, ObservPrime, IntroLucasSeq, Representation, PrimeProblems, QuadraticResidue, DetJacobiSymbol, DetPer-Problem, Permutation, CNT-talk6, ConjDet).
- Not reachable: Google Scholar; the web_search arXiv engine (timed out); the full text of Kalita and Chaliha, Finite Fields Appl. 2026, 10.1016/j.ffa.2026.102905 (no abstract, no arXiv copy); the content of Sun's 2021 book.

## (1) 1208.2683 Conj. 2.1 (fails at n = 9 for 0 < alpha < 0.0023364): NEW
- Versions read: TeX v1-v9.
  - v1-v4 (13 Aug-11 Sep 2012): absent.
  - v5 (13 Sep 2012): "there is n_0(alpha)". This form is not refuted, since it holds for large n.
  - v6 (19 Sep 2012) to v9 (31 Oct 2013): "For any alpha > 0 ... for all n >= 2ceil(alpha)^2+ceil(alpha)+6". Print (143p.pdf, World Sci. 2013, pp. 244-258) is identical.
  - The same statement is Conj. 1.1 of 1209.3729v2 (18 Sep 2012), "alpha any positive real number". It was removed in v3 (19 Sep 2012).
  - Sun, Sci. Sin. Math. 54 (2024) 1357-1364 (243p.pdf, received 2023-11-30), introduction: it restates [6, Conj. 2.1] for any positive real alpha with the same threshold, as open. So Sun still knew no exception in 2023-24.
- Sun's own remarks: Remark 2.1 says "verified for alpha = 2, 3, ..., 700 and n <= 10^6" and gives sharper thresholds for alpha = 2..10. It covers alpha = 1 through Mandl's inequality. It lists no exceptions and checks no alpha < 1 or non-integer alpha.
- Item-specific sources (TeX read):
  - Hou-Sun-Wen 1208.3903, Wang-Zhu 1303.5595, Chen-Guo-Wang 1304.5160, 1309.6025, 1208.5213, Mansour-Shattuck 1410.6943, Sun 1409.5685, Sun-Zhao 2004.01080, Sun 1207.7059.
  - Axler 1409.1777, 1606.06874 and 1609.07934; Visser 2505.04951 and 2507.14410; Hassani math/0606765 (abstract only).
  - Luca-Stanica, J. Comb. Number Theory 4 (2012) (zbMATH review 1339.11038: Bernoulli, Euler and Apery-type numbers only).
  - MO 337896 (Mandl analogues; does not mention it). No OEIS entry. trureturing's screen of 1208.2683 (file `SCREENED-OUT-R71.md`) covers Remarks 2.6, 2.8, 2.10 and Conj. 2.12, 2.13 only.
- Verdict: no prior report. NEW.

## (2) 1208.2683 Conj. 3.1, Lucas part ((A, B) = (3, -26), n = 5): NEW
- Versions read:
  - v1 (Conjecture 2) and v2 (Conjecture 3): "B in Z".
  - v3-v9 (Conjecture 3.1): "B != 0" and "with limit 1". All versions have "A > 1, A^2 > 4B and (A > 2 or B >= -9) ... (u_{n+1}^{1/(n+1)}/u_n^{1/n})_{n>=4} strictly decreasing".
  - Print (143p.pdf) is identical to v9.
- Sun's own remarks: Remark 3.1 says only that (u_n) is increasing. v3-v4 add "See [HSW] for recent progress". The only stated exclusion is A = 2 with B <= -10. Nothing is excluded for A >= 3.
- Item-specific sources:
  - Hou-Sun-Wen (Publ. Math. Debrecen 2014), Thm 1.1: "there exists N" for general (A, B); N = 4 only for Fibonacci.
  - Mansour-Shattuck 1410.6943: "n >= N(seq)".
  - Wang-Zhu, Chen-Guo-Wang, Luca-Stanica: no Lucas part.
  - OEIS A000045, A000129, A000032: no Sun comment on this.
  - Sun's talk slides "Sequence.pdf" (2012) cite only the Fibonacci part.
  - MO, Zenodo, GitHub: nothing.
- Verdict: NEW.

## (3) 1211.1588 Conj. 4.30(iii) (m = 1, n = 2: 11^2 * 3^3 = 3267 > 5^5): NEW
- Versions read: TeX v1-v29.
  - v1-v27 (7 Nov 2012-23 Jul 2015): absent.
  - v28 (26 Jul 2017) and v29 (1 Dec 2017): Conj. 4.30, dated 2015-08-25, with (iii) "for any m, n in Z^+". The v29 changes are bibliography only.
  - Print (176r.pdf, CANT II, Springer 2017, pp. 279-310) is identical.
- Sun's own remarks: Remark 4.30 calls (iii) "an analog of Firoozbakht's conjecture" and gives no exceptions.
  - OEIS A261583 (Sun, 25 Aug 2015; all 14 revisions read) states only the weaker part, "p(m,n)^{1/n} strictly decreasing". This is the right-hand inequality, which the counterexample does not touch.
  - A006450 (Sun, 17 Aug 2015): only p_{p_n}^{1/n} decreasing.
- Item-specific sources: 2609.07371 (2026; cites Conj. 4.4 only); A038580; MO, Zenodo, GitHub: nothing.
- Verdict: NEW. The left-hand inequality appears only in the 2017 paper.

## (4) 1809.07766 Conj. 6.8 (6.15) (n = 7: det = -2 != 0): NEW (misprint level)
- Versions read:
  - v1-v2: absent.
  - v3 (2 Oct 2018), Conj. 4.2(i): "when p = 1 mod 4, det[floor(i^2j^2/p)] = 0".
  - v4-v8 (to 8 Nov 2018): same claim, for primes p = 1 mod 4 only. This is consistent with n = 7.
  - v9 (27 Jun 2019) and v10 (9 Jul 2019, final): Conj. 6.8, "!= 0 iff n is either 9 or a prime greater than 7 and congruent to 3 modulo 4".
  - Print (194s.pdf, FFA 59 (2019) 246-283) is identical.
  - So the failure exists only in v9, v10 and print.
- Item-specific sources (TeX read; none treats (6.15)):
  - 2512.24012, 2512.18581, 2507.18589, 2407.04556, 2404.11547, 2401.14301, 2210.14741, 2401.05853.
  - 2407.07085 (= Chaliha-Kalita, Indian J. Pure Appl. Math. 2026), 2606.22548; 2605.19486 (abstract only).
  - Sun's ConjDet.pdf (2013).
  - MO 315702 is about the other determinant, det[floor((i^2+j^2)/p)]. K. Conrad's comment there gives n = 7 and 25 as exceptions for that matrix. It is a different matrix and does not mention (6.15).
  - OEIS, Zenodo, GitHub: nothing.
- Residual risk: Kalita-Chaliha, FFA 2026 (title "Concerning some determinants involving quadratic residues") is unread.
- Verdict: NEW, minor ("prime > 7" is probably a slip for "prime > 3").

## (5) 1402.6641 Conj. 2.25 (fails at n = 5, 6, 7, 13, 17, 21): NEW (minor)
- Versions read: TeX v1-v12 (25 Feb 2014-25 Feb 2016) are identical; v7+ add "(Cf. [IR, p. 120])". Print (160p.pdf, World Sci. 2015, pp. 169-187) is identical, including the formula pi(sqrt x) + #{sqrt x < p <= x : p != 3 mod 4}. That formula counts each split prime once, not twice.
- Item-specific sources:
  - Sun-Zhao 2004.01080 (pi(kn) only) and Sun 1409.5685 (TeX).
  - OEIS A237578 (rational analogue; no Gaussian version), A055025.
  - Sun's talk slides PrimeProblems.pdf: no Gaussian conjecture.
  - The Evidence API finds only Sun's own formula.
  - MO, Zenodo, GitHub: nothing.
- Verdict: NEW. The statement is false as written, and the Remark's formula is wrong.

## (6) 1811.10503 Conj. 4.9, further instances: NEW; base case KNOWN
- Versions read: TeX v1 (absent), v2 (Conj. 4.9), v3-v5 (Conj. 4.10), v6 (Conj. 4.9): same statement throughout. Print (214p.pdf, J. Algebraic Combin. 54 (2021) 893-912) is identical.
- MO 316142, re-read 2026-11-05: 0 answers, one 2018 comment (Alekseyev: A126972). Last activity 2018. MO 316143: none. OEIS A126972 and Sun's S18 entries (A321597 ... A322099): no restatement.
- trureturing #16084 (2026-10-13/14) covers n = 2 over Q and all allowed fields (every 3-term progression).
  - Their scans: progressions with n <= 6 and |A| <= n+4; all subsets of {0..9} for n = 2, 3, over Q only; and F_p for p <= 29 with A = {0,1,2} only.
  - They list as open a candidate amendment: add (n, |A|) = (2, 3) as a second exception.
- Our instances are not in their scan, and they refute that amendment:
  - F_7, n = 2, A = {0,1,2,4}: 6 < 7.
  - F_11, n = 3 and 4, A = {0,1,2,5}: 9 < 11.
  - F_13, n = 4, A = {0,1,3,4}: 10 < 11.
- Other sources: citers 2206.05021, 2108.07723 (TeX); LAMA 2023, ERA 2020, AAM 2025 (abstracts): they do not discuss Conj. 4.9. Zenodo: nothing.
- Verdict: the base case is KNOWN (trureturing, priority over our page). The F_7, F_11 and F_13 instances are NEW, as minor extensions.

## (7) 2607.07638 Section 4, Conj. (i) (p = 3: 312 == 6 mod 9): KNOWN
- Versions read: TeX v1 (8 Jul 2026), v2 (16 Jul), v3 (21 Jul): identical, "Let p be an odd prime". The identical statement is MO 497120 (Sun, 3 Jul 2025): 0 answers, 0 comments.
- Source: K. Calderon, "Supercongruences for Domb Numbers", arXiv 2610.04772v1, 3 Oct 2026.
  - Abstract: "Sun's original formulation has the exceptional counterexample p=3".
  - Remark rem:p3: "U_3 = 4 + 28*11 = 312 == 6 (mod 9). Thus p=3 must be excluded."
  - The paper proves the corrected statement for p > 3.
  - OEIS A002895 links the paper (revision #291, Michel Marcus, 2026-10-06).
- Verdict: KNOWN. Calderon's report predates our R264 (2026-11-05).

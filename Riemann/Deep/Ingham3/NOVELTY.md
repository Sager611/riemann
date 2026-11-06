# NOVELTY — Ingham3 (6 Nov 2026, 08:00–08:30 UTC)

**Verdict: new.** (a) A prime between n^7 and (n+1)^7 for every n ≥ 1 is the first all-n result below exponent 8. (b) n₀(3) = e^{627.01} is the first single-exponential explicit threshold for Ingham's cubes theorem (before: doubly exponential). (c) k = 4, 5, 6 had no explicit threshold at all.

## Previous records found
| statement | source | status |
|---|---|---|
| prime in (n³,(n+1)³) for n ≥ exp(e^{33.217}) (arXiv v1), exp(e^{33.3}) (published) | Dudek, arXiv 1401.4233, Funct. Approx. 55 (2016) | refereed |
| same for n ≥ exp(e^{32.892}) (x ≥ exp(e^{33.990})) | Cully-Hugill, arXiv 2107.14468 (J. Number Theory 2023) | refereed |
| same for n ≥ exp(e^{32.537}) | Cully-Hugill, PhD thesis (cited in 2601.15564 l. 156) | thesis |
| every n with n³ ≤ 1.649·10^{40}; P₂ between all cubes | arXiv 2601.15564 (Jan 2026; "open problem" l. 142, explicit route "futile" l. 156) | preprint |
| prime between n^k and (n+1)^k for all n: k ≥ 155, 140, 90 | Cully-Hugill 2107.14468; Cully-Hugill–Johnston 2111.10001, 2402.04272 (Thm powerthm: k = 90) | refereed / preprint |
| same for k ≥ 86 | Lee, arXiv 2602.14340 (Feb 2026) | preprint |
| same for k ≥ 8 (from the 7/8 half-plane) | Visser, arXiv 2610.20058 (14 Oct 2026, v1); our ShortPrimes corollary (R323) | preprint / ours |
| "k = 7 would need a half-plane below 6/7" | our R381 (true for the half-plane without density) | ours |

## Searches
- **arXiv monthly listings** 2026-10 and 2026-11 (python3 .scratch/lst.py, regex 'cubes?\b|consecutive (powers|cubes)|Ingham|Legendre|short intervals|between (consecutive|powers)|primes between|power.?s and primes', math.NT; a second regex pass over math.CA, math.GM, math.HO, math.CV): the only relevant titles are 2610.20058 (k = 8) and 2610.27830 (RH-conditional, audited in Riemann/Deep/Visser). Nothing on k = 7 or on explicit Ingham thresholds.
- **arXiv API** (.scratch/axq.py): abs:"consecutive cubes" (2601.15564, 1611.07251 = Dudek's thesis, 0810.2113, 1010.4883; rest unrelated), abs:"consecutive powers" AND abs:prime (2402.04272, 2111.10001, 2107.14468), abs:"between consecutive" AND abs:powers AND abs:prime (adds 2602.22502, RH-conditional, and 2602.14340), author Cully-Hugill (8 papers, newest 2402.04272), Johnston (nothing relevant). Current versions: 2610.20058 v1, 2602.14340 (Mar 2026 source), 2402.04272 (Feb 2026 source), 2601.15564 (Mar 2026 source).
- **web_search (arXiv mode):** "primes between consecutive seventh powers every n zero density half-plane": no relevant result.
- **Zenodo** (API, phrase queries): "consecutive cubes" (1 hit, unrelated), "primes between consecutive" (10: Legendre/squares notes), "consecutive powers" (10, unrelated), "between consecutive cubes" (0), "primes between cubes" (0), "seventh powers" AND primes (0), "7th powers" (1, unrelated).
- **HAL** (API): "consecutive cubes" → hal-02680213 (M. Diouf 2020, "Primes between two consecutive cubes for all positive integers", an elementary g_n = O(p_n^{2/3}) claim; also arXiv 2004.10528). Not accepted: 2601.15564 (2026) treats the problem as open. Not audited here (it would supersede (a) and (b) if correct; it claims an elementary proof of a result far beyond current methods).
- **GXL evidence API** (3 claims: all-n k = 7; single-exponential cube threshold; density + half-plane for k-th powers): only the known results above, a Mills-type "m ≥ 1 438 989" (2004.01285) and a non-rigorous Collatz preprint (1110.3465).
- **Our files:** rg over Riemann/Deep, TREE.md, NOGO.md, DISCOVERIES.md, COROLLARIES.md, REFUTATIONS.md, discoveries/: no earlier consecutive-powers threshold; R381 and Harvest3 (candidate 2 = this task); ShortPrimes (k ≥ 8).

## What is not new
- The method (smoothed explicit formula + zero-free region + explicit zero density) is Ingham's, made explicit by Dudek and Cully-Hugill; the half-plane and KLN are inputs. The novelty is the numbers: a half-plane removes the zeros near Re s = 1 that forced the double exponential.
- KLN's density formula is evaluated by us at every σ (log C₁ is affine in σ for fixed parameters); this re-evaluation is routine.

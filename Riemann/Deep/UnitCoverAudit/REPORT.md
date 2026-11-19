AUDIT PASSED: every step of DERIVATION.md and cert.py checks out against its sources, and independent code (own lattice enumeration, 256-bit interval union bounds, the ExplProg auditor's quadrature for Lemma D, own exact closure) reproduces every number of the claim; each item holds as written, conditional on the Lean-checked half-plane (H) and the cited refereed inputs.

# UnitCoverAudit — referee report on discovery unit-group-cover (Riemann/Deep/UnitCover), 19 Nov 2026

Agent unitaudit (Orca). Read-only web (two arXiv e-prints); nobody contacted; no Lean; UnitCover/ and discoveries/
untouched. Chain with line numbers: FULLPROOF.md. Every process < 300 MB (memcap peaks 15–136 MB; one Lemma D run that
did all 1328 cells in one process grew to 295 MB and was killed after its first half; rerun in two halves, 134/136 MB).

## Verdict per item (claim of discoveries/unit-group-cover.md, verbatim order)

| Item | Holds? | Decisive evidence |
| --- | --- | --- |
| Half-plane: ζ and every L(s, χ) ≠ 0 for Re s > B = 0.8733819666132076, from Bebel's Re s > B₂ + 2.2·10⁻²⁵⁴ | yes | HPAudit (B): every q ≥ 1 (so ζ) and χ, 2.16·10⁻²⁵⁴; B − B₂ = 8.0·10⁻¹⁷ |
| ≤ 2.1067 log Q + 58 primes, all in (e⁻⁷y, y], y = (820 L)^{7.8978}, cover every 2 ≤ q ≤ Q | yes | sup U ≤ 0.546406 (iv, 256 bits); a₁ + b₁ ≤ 2.1066815; δ ≤ 10⁻³ on 664 cells + tail (0.00098366, margin 0.344339); exponent 7.8978 ≥ E only lowers δ |
| = (1/log 2 + 0.664) log Q + 58 | yes | 1/log 2 + 0.664 = 2.106695 ≥ 2.1066815 |
| c' = 820000: 1.9492 log Q + 98 | yes | a₁ + b₁ ≤ 1.949149, a₀ = 96, sup U ≤ 0.477007; δ ≤ 10⁻⁶ (worst cell 9.825817·10⁻⁷, tail 0.345212) |
| Least certified counts 27, 44, 77 (≤ 10^33.33), 164, 498 (≤ 10^39.91), 4858 | yes | each row U < 1 at its (k, r, s, α), locally least; Lemma D at the row's c' and largest N; log₁₀ y 33.325, 39.9028 |
| ≥ log₂ φ(q) primes needed | yes | 2^{|P|} subset products |
| First 17 primes (≤ 59) cover every q ≤ 2000 | yes | own closure; 16 fail (1547, …, 1995) |
| Existence only; random draw succeeds with probability ≥ 0.45 | yes | 1 − 0.5464 = 0.4536 |

## What was checked
- **A. Derivation, line by line** (FULLPROOF §0–§6): construction and no-ops (Lemma 1.1); Lemma D (conditional law,
  ν(A_j), prime powers via RS62, explicit formula with trivial zeros incl. s = 0 for even χ, pairing ρ ↔ 1 − ρ̄ and the
  two envelopes, ζ denominator, Σ1/γ² ≤ 2(1 + B²/196)Σ Re 1/ρ); monotonicity in y and the cell/tail logic; Lemma C
  (Fourier/Parseval, Cauchy–Schwarz, the k-th moment expansion with weights 2^{−k−|e|}, tower property, Markov, closed
  relation sets, kernel count Hom(ℤ^k/span R, Ĝ), BFS completeness, torsion and divisor sums, hole contraction, union);
  the linear theorem's sup over Q (maximiser 1 + x = f/ε, generic class exactly 1); table; small Q.
- **B. Sources.** Half-plane: statement and range (all moduli, primitive or not, ζ as q = 1) from HPAudit l.1, l.20.
  ExplProg inputs (kernel, Lemma Z from BMOR, K1, Zζ): ExplProgAudit PASSED; usage here matches E:229–237. Hughes
  2610.04678 TeX: covering lemma l.132–146, collision step l.149–158, hole contraction l.160–170 = D §3.1/§3.5 with no bad
  characters (B = 0, ρ = (1 + δ)/2, β = δ + (1 − δ)α); Thm main l.71–79 and Thm structured l.87–96 match the page's
  "Before". KST 2609.25393 TeX: Thm l.109–123 ((log Q)^{1+ε} primes ≤ (log Q)^{C_*(A+1)}, all but O(Q(log Q)^{−A})),
  exceptional sets: moduli divisible by chosen primes (l.126), the zero-density set E₀ (Lemma l.338–346: prime sums
  ≪ P^{1−σ}, σ ≤ 1/7, conductors > Y), ineffective Linnik for conductors ≤ Y (l.158, l.252).
- **C. Certificate code** (FULLPROOF §7): every proof step a certified Arb sign; floats only choose parameters, each
  re-certified; outward rounding of printed constants; coverage from log 9 (cells) and log 3 (sup over Q) with no gaps;
  undecided Arb comparisons fall back safely. Rerun of `cert.py all`: identical log, 51.1 s, 23 MB (logs/cert_rerun.txt).
- **D. Independent recomputation** (code in this folder, logs/): indep_classes.py (own extended-gcd lattice insertion,
  canonical HNF, own Smith form: 1, 4, 48, 2340 closed sets, 57 classes for k = 4, tables and sets identical to cert.py;
  brute force on 10 groups), indep_union.py (mpmath.iv 256-bit: every a₁, b₁, sup U and table row), indep_delta.py (all
  664 cells for both c', 20 interior/exponent checks, both tails, 9 table c' values; ExplProgAudit/indep.py quadrature,
  not cert.py), indep_small.py (all §6 claims). Every value equals the certificate's.
- **E. Falsification** (falsify.py): the §3.4 sums against exact sums for all 24 (k − ρ, d-tuple) types up to Q = 4000
  (max ratio 0.64); the Lemma D-tight law ν = (1 − δ)·uniform + δ·1 on 8 small groups, k = 2–4: exact E[C_r^k] ≤ (3.1) and
  exact failure probability (subset DP) ≤ the Lemma C bound. No counterexample; smallest moduli (H = 2) and Q < 9 covered.
- **F. Reviewer** (verdict-only, Lemma C chain pasted inline): passes 1–2 GAP only for hypotheses my paste omitted (N = r + s
  with common draws; 0 < α < 1; both explicit at D:125, D:132); pass 3 with them stated: every step CORRECT.

## Issues found (none load-bearing; nothing to retract)
1. Caveat, prime sizes: "the primes are huge (10^26 to 10^48) and exceed Q for Q below 10^38.4" describes y, not the
   window: P lies in (e⁻⁷y, y], whose lower end is 10^22.67 at Q = 3, and the whole window lies above Q only for
   Q < 10^35.03 (y > Q for Q < 10^38.39). Suggested: "primes in (10^22.7, 10^47.6]; all exceed Q for Q below 10^35".
2. Lemma Zζ silently uses that ζ has no zero with |γ| < 14 (to get 1 + B²/196); classical, but DERIVATION §0 should cite it.
3. Two a₀ are not least for U < 1 (k = 2: a₀ − 1 gives 0.927; k = 4, δ = 10⁻⁶: 0.954), because the search threshold is 0.9
   (C:313); a prime could be shaved; the stated constants stay valid.
4. cert.py prints "N <= 4 log Q + 200" while the proof uses N − 1 ≤ 4L + 200, L = max(log Q, log 9); the code checks the
   latter (C:101, C:90), so only the print is loose.
5. "Empties the exceptional sets of KST and Hughes" (page) is descriptive: the proof follows Hughes's lemma (no bad
   characters by Lemma D; no-ops replace his and KST's exclusion of moduli divisible by chosen primes; the union bound over
   all q replaces his third set). KST's zero-density set E₀ (σ ≤ 1/7) is the kind the half-plane (b = 0.1266 < 1/7) removes,
   but KST's argument is not run.

## Observation: what the result rests on
KST l.129 say the exceptional set is "unavoidable with present technology": a cover of every prime modulus p ≤ Q by
primes ≤ (log Q)^C forces a quadratic nonresidue ≤ (log p)^C (otherwise all subset products are residues), far beyond
Burgess. The discovery has exactly this strength (nonresidue ≤ (820 log p)^{7.8978}), consistent only because (H) is a
quasi-GRH half-plane (ExplProg item 3 gives the same nonresidue bound). So the whole content is (H) plus the
combinatorics audited here; any doubt about (H) transfers to this result in full. The page's caveat states the reliance.

## What survives, dependents
Everything as written. Dependents (rg 'unit-group-cover|UnitCover'): DISCOVERIES.md l.11 (N7 consequences list),
TREE.md l.15 (N11 record), Riemann/Tree/R404.lean l.5, l.17; none in COROLLARIES.md, RECORDS.md, REFUTATIONS.md.
Outside the claim: mc.py's Monte Carlo rerun independently (mc_indep.py, own seed): medians 10, 12, 15 for Q = 100,
300, 1000 (page 10, 12, 14; 40/30/12 draws, within noise). NOVELTY.md was not in this brief's scope.

## Files
FULLPROOF.md; indep_small.py, indep_classes.py, indep_union.py, indep_delta.py, falsify.py, mc_indep.py; logs/ (outputs,
classes_k*.json, cert_rerun.txt, reviewer_passes.txt). Run each with `uv run --no-project --with python-flint --with mpmath [--with numpy]
python3 <file>` (indep_delta.py modes: time | cells [cp i0 i1] | select).

RH relevance: N11

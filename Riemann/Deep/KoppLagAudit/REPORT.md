AUDIT PASSED: every item of unit-generated-orders holds as written; the chain (Lemma A, Corollary B, the bound D₀ from the Littlewood bound, the sieve, the exact class groups, the tables) was refereed line by line and recomputed with independent code, which reproduces D₀ exactly, the 82 survivors and all 90 table entries; no false step, only one cosmetic line-number slip.

# KoppLagAudit: line-by-line referee of unit-generated-orders (Kopp–Lagarias Problem 1.7), 19 Nov 2026

Agent ugoaudit (Orca), N7 → N11. Read-only web access (arXiv e-print 2512.11311 v3, the latest version). Nobody contacted. No Lean, no commits.
Riemann/Deep/KoppLag/ and discoveries/ untouched (cert.py was rerun on a copy in .scratch). Chain with source line numbers: FULLPROOF.md.

## Verdict per item

| Item (page wording) | Holds? | Evidence |
| --- | --- | --- |
| Unit-generated real orders = Δ = n² − 4 (n ≥ 3), n² + 4 (n ≥ 1) | yes | KL Prop. 2.3, TeX l.543–556 |
| Problem 1.7 quoted as KL state it | yes | KL l.329–335 verbatim (page cites l.329–336; \end{prob} is l.335) |
| Cl = Cl[2] exactly for the 53 + 37 table discriminants | yes | own class groups n ≤ 6000 (11998 orders) + own sieve to n = 10⁹ + 1 + own D₀ |
| 89 distinct (Δ = 5 twice), 60 maximal, 29 non-maximal | yes | own TeX parse: 2A 55 entries (53 > 0, plus −4, −3), 2B 37; f column |
| largest Δ = 68640 = 262² − 4, f = 2, Cl (ℤ/2)⁴, Cl⁺ (ℤ/2)⁵ | yes | own: Δ₀ = 17160, h = 16, h⁺ = 32, μ = 6 |
| one class per genus for exactly 84; exceptions 221, 725, 1517 (ℤ/4), 3965, 7565 (ℤ/2 × ℤ/4) | yes | own h⁺ vs 2^{μ−1}; TeX 2A l.1647–1648, 1655 |
| Δ < 10¹⁸: elementary + exact computation, no GRH, no zero-free region | yes | Lemma A, Cor. B, sieve and forms use none (forms ↔ Cl⁺ and the reduced-cycle theorem are classical, cited) |
| Δ ≥ 10¹⁸ from \|L(1, χ)\| ≥ 0.11/max(3, log log q), q ≥ 3 primitive | yes | D₀ = 20485358601780401 < 10¹⁸ (own mpmath intervals, exact match); (L) = LittlewoodAudit item (a), PASSED |
| any c ≥ 0.037 suffices | yes | own: D₀(0.037) = 9.6959·10¹⁷ ≤ 10¹⁸ (base k = 16); D₀(0.036) = 1.0279·10¹⁸ |
| Caveat (F2, genus theory cited; PARI GRH-conditional) | yes | F2 rechecked with own L(1, χ) on 889 orders (\|Δh\| ≤ 1.7·10⁻¹⁰); μ = KL l.1527–1537 = Cox; genus count asserted on 11998 orders |
| Side remark: KL l.1609 "86" vs 89 distinct | yes | 90 positive entries, 89 distinct |

## Refereed steps (details in FULLPROOF §1–§3)

- Lemma A (DERIVATION l.46–68): normalisation, (A1)–(A3), parity, concavity, both endpoint constants (1.5639, 3.1440 at Δ = 100) and their
  monotonicity rederived. Reviewer subagent (verdict-only, Lemma A + Corollary B pasted inline): **CORRECT**, recomputing F(10) = 3.8 − √5 =
  1.5639…, G(10) = 8.6 − 2.44√5 = 3.1440… and both derivatives positive.
- Corollary B (l.71–79): splitting, primitivity of 𝔭², norm p² ✓.
- D₀ (l.81–101): R(O) ≤ log ε_n < ½ log(Δ + 4) for both signs, (L) at q = Δ₀ ≥ 5 with Δ₀ ≤ Δ, ∏(1 − χ(p)/p) ≥ 1/P_{ω(Δ)}, g increasing on [16, ∞),
  h ≤ |Cl⁺/Cl⁺²| = 2^{μ−1} ≤ 2^{ω}, the k-split and the induction for k ≥ 15 (each inequality and its direction) ✓.

## Code audit (FULLPROOF §4)

- cert.py / forms.py read line by line: exact integers and arb balls only; no coverage gap (n ≤ 1000 exact: 998 + 1000 = 1998 orders;
  sieve n ∈ [1001, 10⁹ + 1] for N = +1 and [1001, 10⁹] for N = −1, contiguous chunks, chunk primes valid for every member, every chunk survivor
  re-tested with all primes 4p⁴ < Δ; lesson 277 ✓); primitive reduced forms by exact integer tests; ρ-cycles asserted closed and disjoint;
  composition test and composition-free inverse test asserted to agree; the KL table transcription in cert.py equals our TeX parse.
- sievecheck.c ✓ (same range, 128-bit). checks.py and the PARI scripts read; cross-checks only, not load-bearing.
- cert.py rerun on a copy (.scratch/ugoaudit/certcopy): logs/cert_output.txt identical line by line after masking timings (130 lines,
  logs/cert_rerun_output.txt), survivors.txt identical; 137 s, memcap peak 149 MB.

## Independent recomputation (own code; logs/, footprints in logs/footprints.txt, all < 40 MB except the cert.py rerun)

- d0check.py (mpmath interval arithmetic, not flint): D₀ = 20485358601780401 exactly, worst k = 14, base k = 15, step 3.02206 > 2.03846;
  c = 0.037: D₀ = 9.6959·10¹⁷ (base 16); c = 0.036: 1.0279·10¹⁸ > 10¹⁸.
- ugo.c sieve (Euler-criterion tables, early exit): 58 + 24 survivors on n ∈ [1001, 10⁹ + 1], identical to KoppLag's.
- ugo.c class 6000: all 11998 orders n ≤ 6000, h⁺ by ρ-cycles, 2-torsion decided by genus theory instead of composition (h⁺ = 2^{μ−1};
  h = 2^{μ−1} or 2^{μ−2} as J ∈ Cl⁺² or not, via the genus characters at −1): the 2-torsion orders are exactly the 90 table entries, with
  Cl, Cl⁺ and f as in the TeX; the 82 survivors have the same h and h⁺ as survivors.txt (21 ≤ h ≤ 192) and none is 2-torsion.
- Class number formula with our own L(1, χ) (log-sine sum): 889 orders (all n ≤ 300, the 82 survivors, every n > 1000 with 47 | n),
  max \|h_an − h\| = 1.7·10⁻¹⁰ (double precision; a check, not a proof).

## Falsification

- Smallest cases: Δ = 5 (both families), 8, 12, 13, 20, 21: in the tables, recomputed. Δ near 100 (Lemma A's threshold): the proof never uses Lemma A
  below Δ = 10⁶; Lagrange's form of Lemma A holds on all 11998 orders (Δ ≥ 100) and on 4008 principal cycles (n ∈ [10⁶, 10⁹], n = 10⁹ − 1, 10⁹, 1000, 1001).
- Thresholds: n = 1000/1001 and the sieve end are covered on both sides; the last survivor (n = 5968) lies inside our exact range n ≤ 6000;
  D₀ is 49 times below 10¹⁸; the constant threshold 0.037 is tight at 3 decimals for this method.
- Tight families: μ − 1 = ω (the equality case of the genus step in §3.4) holds exactly when 32 | Δ: 18 listed orders (32, 96, …, 27552, 68640;
  at 68640: ω = 5, μ = 6, h⁺ = 2⁵, h = 2⁴), all recomputed; the five narrow exceptions are exactly the listed orders with J ≠ 1 but J ∈ Cl⁺² (class6000.log: N = +1, JinSq = 1, h⁺ = 2h; Δ odd, all primes ≡ 1 mod 4: 13·17, 5²·29, 37·41,
  5·13·61, 5·17·89), where Cl⁺ ≅ (ℤ/2)^{μ−2} × ℤ/4. No counterexample to any lemma or to the list.

## Slips (cosmetic)

- S1 The page (Before) cites Problem 1.7 at TeX l.329–336; it is l.329–335 (\begin{prob} … \end{prob}). No content change.

## Dependents (all stand)

DISCOVERIES.md l.11 (half-plane consequences line), TREE.md l.15 (N11 entry), TREE.md l.450 and Riemann/Tree/R414.lean (round 14 of N11).
Nothing in COROLLARIES.md, RECORDS.md, NOGO.md or other discoveries/ pages builds on it. Had (L) failed, the complete list for Δ < 10¹⁸ would
survive unconditionally, and any constant c ≥ 0.037 restores the full statement.

## Not re-proved (cited classical inputs, numerically confirmed)

Forms ↔ Cl⁺(O) for orders (Halter-Koch Thm 6.4.2), one ρ-cycle of reduced forms per proper class (Cohen §5.6), genus theory for orders
(Halter-Koch Thm 6.5.2, KL Prop. 3.4), the class number formula for orders (Neukirch I.12.12); (L) via LittlewoodAudit and HPAudit.
Lemma A's sharpness remark (checks.py C3b) and the novelty search were not redone (outside the brief).

## Files

FULLPROOF.md (chain with TeX and DERIVATION line numbers); ugo.c (modes `sieve`, `class 6000 300 47`, `pcyc 2000`; `cc -O2 -o ugo ugo.c -lm`);
d0check.py (`uv run --no-project --with mpmath python3 d0check.py`); parse_tables2.py, compare.py (TeX parse, comparison);
logs/: class6000.log, sieve.log, d0check.log, pcyc.log, compare.log, cert_rerun_output.txt, footprints.txt.

RH relevance: N11


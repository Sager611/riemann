AUDIT PASSED: every load-bearing step holds (theorem and hypotheses, kernel, worst-case tail, zero count, code), an independent 28 MB recomputation and mpmath falsification reproduce both bounds, and the condition narrows to Weisz's 14,400 ordinates being right to 8·10⁻⁹⁶³⁹ (his c values match ours to 7·10⁻²⁵ for all zeros, enough to drop them from the condition).

Audit of the conditional record mertens-limsup (discoveries/mertens-limsup.md, Riemann/Deep/MLimsup/), node N10.
Agent mlimaudit (Orca), 2026-11-19/20 UTC. No Lean, no commits. Chain with source lines: FULLPROOF.md.

## 1. Method, line by line (FULLPROOF S1–S3)
| Item | Holds as written? | Evidence |
|---|---|---|
| Theorem (OtR form): K ∈ C², K ≥ 0, even, O((1+y²)⁻¹); k(0) = 1, k = 0 on \|t\| ≥ T; zeros \|γ\| < T on the line, simple ⇒ lim inf ≤ h_K(y₀) ≤ lim sup for every y₀ | yes | OtR p. 144, read from the scan; (2.10) p. 143 gives k(t) = ∫K e^{−ity} |
| "If RH fails or a zero is multiple, both limits are infinite" | yes | OtR p. 141: RH + simplicity follow from either one-sided bound; BT TeX l. 103–107 |
| Page's short form omits K ∈ C², K even, K = O((1+y²)⁻¹) | wording gap only | REPORT §1 of MLimsup states them; the kernel satisfies them (S2) |
| Zero hypothesis for \|γ\| < 18292.5 | yes, unconditionally | vdL–tR–W 1986 (1.5·10⁹ zeros); not even needed given p. 141 |
| k = (g∗g̃)/‖g‖², g step, 1350 cells, h = 13.55: K = \|ĝ\|²/(2π‖g‖²) ≥ 0, k(0) = 1, supp [−18292.5, 18292.5] | yes | S2; b_j > 0; R(d) rechecked by an fmpz_poly product |
| Rayleigh quotient / top eigenvector | algebra correct; not load-bearing | S2: it only chooses b; any integer b is admissible |
| "the best positive-definite kernel" | overstated | best only within 1350-cell step kernels of width 13.55, rounded to integers; harmless |
| Zeros 14401–20311 charged −2\|k(γ)\|\|c\| | yes | S3; valid for every phase; \|c\| upper bounds rigorous (Arb), mpmath check §4 |
| N(18292.5) = 20311; zeros isolated by acb.zeta_zeros | yes | Arb zeta_nzeros exact; zero 20312 > T₁ and 20311 < T₁ asserted (L); mpmath bracket §4 |
| Role of Weisz's balls | as stated | zeros 1–14400: γ (9860 digits, radii ≤ 8.6·10⁻⁹⁸⁴⁰) and c (60 digits); y± are exact decimals, need no trust |

## 2. Code audit of certify.py (FULLPROOF S4)
- Input balls follow Weisz's rule (mid ± r ± one unit in the last digit); all arrays are strings (no float underflow).
- Phases: p_hi = 32,818 bits (page: 32,819, off by one). |γy| < 2^31994, so ~248 digits survive after the point; the data
  radius gives |y|·rad γ ≤ 4·10⁻²¹³. Reduction θ − kk·2π is rigorous for any integer kk. Enough bits: yes.
- Integer kernels exact; node-straddle and support asserts; sums in 320-bit balls; final tests are certain arb comparisons
  of exactly the two claimed inequalities. Extension: exact T₁ and N(T₁), |c| upper bounds in Arb. No bug found.
- Cosmetic: certify.log's ± 4·10⁻⁴⁰ / ± 3.2·10⁻²⁴ are print rounding; true radii 9.9·10⁻⁶⁰ (known part) and
  2.5·10⁻³⁴ (final, 128 bits). Margins over the claimed constants: 5.98·10⁻¹¹ (lim sup), 1.09·10⁻¹¹ (lim inf).
- 0.9 GB = numpy <U9860 gamma array (568 MB) + zip bytes; nothing else needs memory.

## 3. Recomputation under the cap (lean_certify.py, lean_certify.log)
- Independently written: rows streamed from the npz (npzstream.py), θ/(2π) at 33,200 bits and cos_pi/sin_pi of the
  fractional part (C: θ − kk·2π then cos), R(d) by an fmpz_poly product (C: double sum), extra zeros at 128 bits (C: 96).
- Result (370 s, one core): S₊ = 2.0766736799154281625607508194758677440 ± 9.9·10⁻⁶⁰, S₋ = −2.0789843657637535331593121186068539176
  ± 9.9·10⁻⁶⁰, P₊ = 0.00277353815562681239, P₋ = 0.00256651915283701111 (N(T₁) = 20311, 5,911 zeros):
  h_K(y₊) ≥ 2.073900141759801350167752 > 2.0739001417 and h_K(y₋) ≤ −2.076417846610916522047037 < −2.0764178466.
  Matches certify.log to all printed digits.
- Weisz's JP values: overlap with his certificates; |Δmid| = 3.6·10⁻⁶⁰ (plus), 4.6·10⁻⁵⁹ (minus), inside his radii.
- Footprint: memcap peak 28 MB; /usr/bin/time -l: 19.6 MB peak footprint, 49.6 MB max RSS (cap 290 MB).

## 4. Falsification (mpmath, independent of Arb; mp_falsify.json; FULLPROOF S5)
- 50 extra zeros (n = 14401 + 120k, k < 49, and 20311): mpmath γ within 4.7·10⁻²¹ of Arb's; every Arb upper bound for
  |c| = 1/|ρζ′(ρ)| lies above mpmath's value, at most 5.7·10⁻¹⁵ relative (file print rounding). No violation.
- Count: γ₂₀₃₁₁ = 18291.652 < T₁ = 18292.5 < γ₂₀₃₁₂ = 18292.886, so N(T₁) = 20311.
- Phases θ/π mod 2 at 33,300 bits for 9 zeros (1, 2, 3, 72, 1000, 5000, 10000, 14399, 14400) at both points: Arb agrees
  to 4.9·10⁻⁵⁰ (the precision of the comparison file).
- Kernel values by a direct breakpoint sweep of ∫g(x)g(x−γ)dx at the same zeros, both kernels: agree to 3.7·10⁻³¹; the
  terms 2k Re(c e^{iθ}) agree to 3.4·10⁻³¹.
- c = 1/(ρζ′(ρ)) at 75 digits from Weisz's γ, zeros 1, 2, 72, 14400: his c strings agree to 4.8·10⁻⁶² (convention confirmed).
- Full sums in mpmath (mp_fullsum.py: θ/(2π) mod 1 at 33,300 bits for all 14,400 zeros, kernel from R, 80-digit sums;
  25 MB): S₊, S₋ and both JP sums agree with Arb to ≤ 4.4·10⁻⁴⁵ (the print precision of the Arb values).
- Fresh full-precision data check (fullprec_newton.py, own code, 33,100 bits, random zeros 7338, 8578, 9656, seed 20261120):
  interval Newton N(B) ⊂ B on twice Weisz's box, Im N(B) in his γ ball, c in his c balls: 3/3 (memcap peak 189 MB).
- Weisz's own certified JP values: reproduced by both certify.py (to 10⁻⁴⁰, its test tolerance) and lean_certify.py
  (inside his radii, §3).

## 5. Reviewer (most delicate step: theorem application with the 1350-cell kernel, the worst-case split and the case split)
- Pasted inline with exact formulas (FULLPROOF S1–S3 condensed). Verdict: **CORRECT**. Checked: hat-function
  interpolation of g∗g̃ and normalisation hR(0); K = ĝ(y)ĝ(−y)/(2π‖g‖²) = |ĝ|²/(2π‖g‖²) ≥ 0 under k(t) = ∫K e^{−ity};
  C² and O((1+y²)⁻¹); conjugate pairing; |E| ≤ P uniformly in y; completeness via N(T₁) = 20311 counted with
  multiplicity; case split applied separately to lim sup and lim inf. It validates the deduction, not the numbers
  (those are §3–§4). Session: reviewer run 2b176bc8.

## 6. Dependents (rg -n -i 'mlimsup|mertens-limsup|2.0739' discoveries/ TREE.md COROLLARIES.md RECORDS.md DISCOVERIES.md FORMALIZATIONS.md)
- RECORDS.md l. 19 (N10 conditional record), TREE.md l. 14 (N10 node text: "limsup M(x)/√x > 2.07390 and liminf < −2.07641,
  conditional on Weisz's published …"), discoveries/mertens-limsup.md. Not in DISCOVERIES.md, COROLLARIES.md,
  FORMALIZATIONS.md. Riemann/Tree/R389.lean and Riemann/Deep/MertensX/REPORT.md cite Weisz only, not this record.
- Nothing downstream uses the constants; no dependent breaks.

## 7. The condition, narrowed (lowprec_all.py, lowprec_all.json; FULLPROOF S6)
- All 14,400 Weisz zeros rechecked against Arb acb.zeta_zeros at 96 bits (326 s, 24 MB): every index n matches
  (|γ^W − γ^Arb| ≤ 5.2·10⁻²⁵), every c midpoint is within 6.8·10⁻²⁵ of our rigorous c, no ball disagrees.
- So the c half of the condition is discharged: replacing Weisz's c by the true c moves h by ≤ 7.9·10⁻²³ (margins:
  5.98·10⁻¹¹ lim sup, 1.09·10⁻¹¹ lim inf).
- What remains: the published ordinates must be right to ε < margin/(|y|·2Σ|k||c|): 1.7·10⁻⁹⁶³⁶ (y₊), 8.0·10⁻⁹⁶³⁹ (y₋),
  i.e. ≈ 9,643 significant digits of each γ₁ … γ₁₄₄₀₀ (Weisz states ≤ 8.6·10⁻⁹⁸⁴⁰). About 55 distinct zeros are checked at
  full precision (spotcheck2.log 50, spotcheck.log 3, §4 3), all pass; a full check costs ~6 core-days. The record stays conditional.

## 8. Corrections to the page (none affects the claim)
- "32,819 bits" → 32,818 (certify.py's p_hi; Weisz's certificate also says 32818).
- "± 3.2·10⁻²⁴", "± 4.5·10⁻²⁴", "± 4·10⁻⁴⁰": print rounding of str(30)/str(40) (valid but loose enclosures); the lean
  rerun's true radii: 9.9·10⁻⁶⁰ (known part, the same for certify.py's 320-bit sums) and 2.5·10⁻³⁴ (final, 128 bits).
- "reproduces Weisz's certified values to 10⁻⁴⁰": that is the test tolerance; the agreement is inside his ±1.3·10⁻⁵⁹ / ±5.6·10⁻⁵⁹.
- Method sentence omits K ∈ C², K even, K = O((1+y²)⁻¹) (OtR hypotheses; satisfied). "Best positive-definite kernel":
  best within 1350-cell step kernels only. "Earlier, zeros 7213 and 14400": spotcheck.log also has zero 3.
- "|y| ≈ 10^9626": y₊ ≈ −1.27·10⁹⁶²⁵, y₋ ≈ +4.96·10⁹⁶²⁶ (sign irrelevant). "lattice-found": Weisz does not say how.
- Condition can be weakened to the γ digits alone (§7).

## 9. Not done
- Full-precision recheck of the other ~14,345 ordinates: the remaining condition (§7). Novelty not re-audited.
- kernel.py (the optimiser) not rerun: not load-bearing (§1). The Krein / Boas–Kac optimality remark not audited.

## 10. Files (Riemann/Deep/MLimAudit/)
- lean_certify.py + npzstream.py (streamed recomputation; lean_certify.log, lean_dump.json, ext_zeros_128.txt);
  mp_falsify.py (mp_falsify.log, mp_falsify_cde.log, mp_falsify.json); mp_fullsum.py (mp_fullsum.log, .json);
  lowprec_all.py (lowprec_all.log, .json); fullprec_newton.py (fullprec_newton.log); FULLPROOF.md.
- Rerun (data: Zenodo 23651459 zip, unzip the npz next to the extracted package in DATADIR; each < 290 MB, one core):
  `uv run --no-project --with numpy --with python-flint==0.9.0 python3 Riemann/Deep/MLimAudit/lean_certify.py DATADIR` (6 min);
  mp_falsify.py and mp_fullsum.py need `--with mpmath` (5 min, 3 min); lowprec_all.py (6 min); fullprec_newton.py n… (65 s/zero).

RH relevance: N10

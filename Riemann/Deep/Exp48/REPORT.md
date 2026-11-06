NO RECORD: KNOWN — Jaswal (Zenodo 23251444, 9 Oct 2026; v5 23547020, Thm 5.2) already proved that, given no real zero of L(s, χ) in (7/8, 1) (implied by the Lean-checked half-plane), the lists of Elsenhans–Klüners–Nicolae of imaginary quadratic fields with class group of exponent 2, 4, 8 (56, 203, 778 fields) are complete; the extension to all discriminants (Tables 4 and 8 of arXiv 2605.26020: 485 and 2329) and its boundary-CM theorem follow in one line (routine corollary, not stated in any source found).

# Class groups of exponent 4 and 8 without a Siegel-zero hypothesis (node N11)

Agent exp48, 2026-11-06, 07:41–08:40 CET. Brief: `.scratch/orca/exp48-brief.md`. Files at the end.

## Claim

A class group has exponent E if every class raised to the power E is trivial and E is the least such number.
Elsenhans–Klüners–Nicolae (EKN, arXiv 1803.02056, Acta Arith. 193 (2020)) listed the imaginary quadratic fields of small
exponent but proved completeness for exponents 2, 4, 8 only up to "at most one missing field", which would carry a Siegel
zero. Given the half-plane (TREE.md N7), through our Littlewood78 bound or directly:

- **Theorem A** (= Jaswal's Thm 5.2 with his hypothesis (H) supplied). A fundamental discriminant D < 0 has E(D) = 4 iff it is
  one of EKN's 203 (largest |D| = 435,435), and E(D) = 8 iff it is one of EKN's 778 (largest |D| = 430,950,520); 1046
  fields have E(D) | 8 (9 + 56 + 203 + 778).
- **Theorem B** (routine corollary; not found stated). The tables of Aiken–Ross–Shvydkoy–Xue (ARSX, arXiv 2605.26020) over
  all negative discriminants are complete for E = 2, 4, 8: 88, 485, 2329 discriminants (largest |D| 7392, 887,040,
  1,723,802,080). Their Thm all-CMpts-on-bdry holds unconditionally: the 37 discriminants of its table are all D whose CM
  points in the fundamental domain all lie on its boundary.

Which published statements become unconditional, and why:

| source | statement | status now | reason |
|---|---|---|---|
| EKN Thm 2, 2nd sentence | "assuming no Siegel zeros … all fields with exponents 2, 4, 8" | unconditional | the class-number bound EKN's search assumes holds for every field (Littlewood78; or (H) as in Jaswal) |
| EKN Thm 2, 3rd sentence | "at most one missing field with exponent 2, 4 or 8" | no missing field | same |
| EKN Thm 2, 1st sentence | ERH: all fields of exponent ≤ 5 and 8 | E = 1, 2, 4, 8 unconditional; E = 3, 5 still ERH | odd exponents have no 2-rank bound; they need a small split prime (ERH via Bach–Sorenson) |
| ARSX Thm exp-tables (2), (3) | E = 2, 4, 8 tables complete under no Siegel zeros; else one more fundamental D | unconditional | ARSX l.666–668: EKN's fundamental lists are their only conditional input |
| ARSX Thm exp-tables (4)–(6) | E = 3, 5, 6, 7 under ERH | unchanged | as above |
| ARSX Thm all-CMpts-on-bdry | table conditional on no Siegel zeros | unconditional | needs only the E | 2 tables (exponent 2: OpenAI's 5 Oct note / Weinberger) |

## Before

- EKN (2018): complete under "no Siegel zeros", otherwise at most one missing field; Thm 1 (unconditional sieve) covers
  |D| ≤ 3.1·10^20 for every exponent ≤ 8.
- ARSX (May 2026): the same two alternatives for all discriminants.
- OpenAI, "The Quasi-Riemann Hypothesis", 5 Oct 2026 (§1): EKN Thm 2 plus the half-plane gives exponent ≤ 2, hence Euler's
  65 idoneal numbers; exponents 4 and 8 not mentioned.
- Jaswal, Zenodo 23251444 (v1, 9 Oct 2026) to 23547020 (v5, 21 Oct): Prop. 5.1, unconditional, exactly 1046 fundamental D
  with |D| < 2·10^9 and E(D) | 8; Thm 5.2 under (H): EKN's exponent-2, -4, -8 lists are complete. ARSX's all-discriminant
  tables are mentioned (l.1033–1035 of v5) but not claimed complete.

## Method

1. Read EKN's TeX (DERIVATION.md §1): a Siegel-type input enters in exactly two places, the bound N_c on the number ω(D) of
   prime factors (Thm upp_bound; Lemma Siegel or Siegel_one) and the per-ω, per-Rédei-rank caps B_s of the §6 search (Lemma
   Chen's inequality (b_ch), l.695–696). EKN's possible missing field is a field violating (b_ch) (l.705–710). Everything
   else (Lemma 1, genus theory, Rédei ranks, the §7 direct search) has no analytic input.
2. Replace both by Littlewood78, h(D) ≥ 0.0350·√|D|/max(3, log log|D|) for every fundamental D < −4 (cert.py):
   it exceeds Chen's √|D|/(7.732π log|D|) for |D| ≥ 35 (ratio ≥ 3.9 for |D| ≥ 10^6) and EKN's 0.655/(πe)·√|D|/log|D|
   for |D| ≥ e^{6.6}; with genus theory (h ≤ c^{ω−1} if E | c) and |D| ≥ p_1⋯p_ω it gives ω(D) ≤ 52, 20, 9 for E | 8, 4, 2,
   inside EKN's searched ranges ω ≤ 58, 24, 11. So EKN's search misses no field. Jaswal's route instead uses (H) to verify
   the hypothesis of EKN's Lemma Siegel and Tatuzawa's lemma to exclude the exception; ours avoids Tatuzawa's zero
   localization. Both rely on EKN's §6 computation for large |D|.
3. Independent check without L-functions or GRH (scan.c): every discriminant D < 0, fundamental or not, with |D| ≤ 10^10;
   a split prime q with 4q^8 < |D| or a prime form with f_q^8 non-principal rejects D; survivors get all reduced forms
   counted (exact h) and every 8th power tested.
4. Theorem B (conductor.py): D = D_0 f^2 with E(D) | E forces E(D_0) | E and, by ARSX's Lemma property-of-f
   (Fan–Pollack Prop. 2.3), lcm_{p^k‖f} p^{k−1}(p − (D_0/p)) | 12·|O^×|·E. Every such f with |D_0 f^2| > 10^10 (964
   candidates after cuts) is excluded by an exact test on a divisor.

## Certificate

One command (python-flint arb, 256 bits, about 1 s): `cd Riemann/Deep/Exp48 && uv run --no-project --with python-flint --with mpmath python3 cert.py`
→ `logs/cert_output.txt`: ALL PASS — (1) Littlewood78 > Chen for d ≥ 35; (2) > EKN's Siegel bound for d ≥ e^{6.6}, and
1/(4(1 − B)) = 1.9744; (3) K_2, K_4, K_8 = 10, 21, 53 (exact to k = 399, induction beyond); (4) EKN's N = 11, 24, 58 and
d_N ≤ 2.01·10^11, 2.38·10^34, 3.17·10^110 recomputed; (5) (1) and (3) again with OpenAI's 7/8 fallback (constant 0.0334:
Chen exceeded for d ≥ 41, same K_c).
Scan: `cc -O2 -o scan scan.c -lm; ./scan 10000000000 > list.txt` (or in pieces `./scan X Xlo`; 9.5 min single
thread for 10^10, 2 MB), then `python3 conductor.py list.txt 10000000000` (5 s, 18 MB) and `python3 lemma_check.py`.

## Checks

- Scan to 10^10 equals ARSX's Tables 1, 2, 4, 8 discriminant by discriminant (13, 88, 485, 2329) and EKN's counts and
  maxima for fundamental D (9, 56, 203, 778; 163, 5460, 435,435, 430,950,520); nothing with E | 8 in (2·10^9, 10^10]
  (`logs/cmp_output.txt`, `logs/scan_stats.txt`). EKN's nine fields with |D| > 10^7 all have exponent 8.
- PARI (`spot.gp`, `logs/spot_output.txt`): the largest fields of exponent 2 and 4 and the nine of exponent 8 with |D| > 10^7,
  structures certified by bnfcertify, e.g.
  −430950520 = −2^3·5·7·11·13·47·229: (Z/8)^3 × (Z/2)^3, h = 4096; −435435: (Z/4)^2 × (Z/2)^3; −1723802080 (conductor 2 over
  −430950520): (Z/8)^3 × (Z/2)^4; all agree with the scan's h and 2-, 4-torsion counts.
- conductor.py self-test: its exclusion test never fires on the 3602 known pairs (D, E) with E(D) | E and finds an element
  of order 8 in all 2329 exponent-8 discriminants; without cuts D_0 = −3 has 786,432 candidates for E = 8 (ARSX README:
  0.8M). lemma_check.py: L(f) | 12uE(D) holds on all 1869 non-fundamental D with E(D) | 8 up to 10^10.
- Boundary-CM table: the 37 odd-or-(−4) discriminants with E | 2 found by the scan equal ARSX's table.
- Self-review (DERIVATION.md §6): every Siegel-type use in both TeX sources listed; D = −3, −4 (units) have h = 1 and
  Littlewood78 is used only for D < −4; D ≡ 0 mod 4 and mod 8 are covered (χ_D primitive; EKN's four 2-parts; the scan runs
  every D ≡ 0, 1 mod 4); d < 10^6 checked by the scan; no step uses GRH (PARI's GRH-conditional quadclassunit appears only in
  spot.gp, certified there). `constcheck.py`: of 169 numeric tokens in DERIVATION.md, the 17 not found in cert.py, the code or
  the logs are years, identifiers, facts quoted from the TeX (EKN's maxima for exponents 3, 5, 6, 7, Chen's 1.5·10^6, the
  sieve count 1,002,279, 3.14·10^20), the sum 4025 = 2 + 4023, and roundings of printed balls (0.0411679, 0.0767003,
  9.34·10^−4).
- Observation (not a flaw that matters): EKN l.602–604 say Lemma Chen gives ω bounds 11, 24, 58; with Chen's constant the
  three conditions of Thm upp_bound give 12, 24, 59 (the stated numbers use Lemma Siegel's constant), and the unconditional
  half of Thm upp_bound has its own Tatuzawa exception; EKN do not show that it coincides with Lemma Chen's. Both
  exceptions vanish here.

## Novelty

KNOWN: Jaswal's Zenodo preprint (concept 23251443), §5.1, in every version from v1 (9 Oct 2026). Searched: the TeX of
2610.21099, 2610.30476, 2610.15483, 2610.31402, 2610.21511, 2610.18491, 2610.26835 (no mention); OpenAI QRH 30 Sep, 5 Oct
(exponent 2 only) and Landau–Siegel 1 Oct; Zenodo 23607492 (E(D) → ∞ under a half-plane, no lists); arXiv listings
math.NT and math.GM for 2026-10 and 2026-11 (titles: exponent, idoneal, imaginary quadratic, Siegel zero) and Zenodo
'"class group" AND exponent' (`.scratch/lst.py`); arXiv API full archive (3 queries); HAL (no hit); our files (Jaswal
recorded only as a class-number source). Theorem B and the boundary-CM corollary were not found stated anywhere, but
they are one line from Jaswal's Thm 5.2 and ARSX l.666–668. Details: NOVELTY.md.

## Caveat

- The main statement is Jaswal's; our contributions are a second proof route (Littlewood78, no Tatuzawa localization), the
  unconditional scan of all discriminants to 10^10, and the conductor check that closes ARSX's tables.
- Beyond 10^10, fundamental fields rest on EKN's unpublished §6 computation (refereed paper; Jaswal relies on it beyond
  2·10^9); EKN's unconditional sieve (Thm 1) covers |D| ≤ 3.1·10^20.
- Littlewood78 is a reviewed paper proof, not formalized; Theorem B uses Fan–Pollack's Prop. 2.3 through ARSX's lemma
  (tested above, not re-proved). No Lean. No reviewer verdict on this report yet.
- Suggested bookkeeping (orchestrator): COROLLARIES.md row for Theorem B and the boundary-CM corollary, citing Jaswal Thm
  5.2 + ARSX l.666–668, verified by exp48's scan; Harvest3 candidate 1 → KNOWN (Jaswal).

RH relevance: N11

## Files

- `REPORT.md` (this file), `DERIVATION.md` (every step and constant), `NOVELTY.md`.
- `cert.py` → `logs/cert_output.txt`: arb certificate.
- `scan.c` → `logs/scan_E_divides_8_upto_1e10.txt` (2915 lines: D, fundamental, h, E, #2-torsion, #4-torsion),
  `logs/scan_stats.txt`; comparison `cmp.py` → `logs/cmp_output.txt` (fetches ARSX's data, prints its sha256).
- `constcheck.py` (self-review of DERIVATION constants); `logs/novelty_lst.txt`, `logs/novelty_websearch.txt`.
- `conductor.py` → `logs/conductor_output.txt` (X1 = 10^10) and `logs/conductor_output_X1_2e9.txt`; `lemma_check.py` → `logs/lemma_check_output.txt`;
  `spot.gp` → `logs/spot_output.txt`.

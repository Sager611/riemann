AUDIT PASSED: every step of PROOF.md checks against Washington–Yang's TeX and every computation reproduces (reruns under the 290 MB cap, two independent recounts of the gap, falsification hunts empty); no error, six harmless remarks (F1–F6); inherited inputs listed below.

# Referee audit of the odd Robin criterion (discoveries/odd-robin.md, Riemann/Deep/OddRobin) — oddraudit, 14 Nov 2026

Scope: PROOF.md §§1–6 line by line ([FULLPROOF.md](FULLPROOF.md), with source line numbers), all certificates rerun, the gap recounted by a third method, falsification hunts. Sources: WY arXiv 2008.04787 — arXiv has only v1 (11 Aug 2020); the journal version (IJNT 17 (2021) 843–870) is paywalled and was NOT compared.

## Verdict by part of the claim

| part | status | basis |
|---|---|---|
| Gap: every odd n, N < n < c0, satisfies (OR) | holds | §1 refereed; gap.py, dfs_gap.py rerun; enum_s.py (new method) finds the same 1,098 n, 0 violators; values.py recomputes every constant |
| N violates (sharpness) | holds | δ0 = 0.00121349015854 > 0 (mpmath interval, 80 digits) |
| RH ⟺ (OR) for every odd n > N | holds | = WY Thm 1.2 (TeX l. 87; = Thm 3.6, l. 440) + the gap; both directions checked against l. 440–456 |
| (OR) for N < n ≤ e^(5.0·10^8), arithmetic only | holds | oddca rerun to x_1 = 5·10^8 (full range, under the cap), identical output |
| (OR) for N < n ≤ 10^(10^46.44) | holds given RobinCH | transfer rerun (1221 boxes, min 1.07·10^−5); Lemmas O1, O2 and the use of RobinCH's margins checked; RobinCH's own lemmas inherited |
| (OR) for N < n ≤ 10^(10^46.79) | holds given RobinCH + half-plane hypothesis | transfer rerun (1209/1235 boxes); conditional, as PROOF.md says |
| Under RH, f_{2,1}(n) < e^γ/2 for every n > N; N the largest exception | holds | (b) + even-n reduction checked; f_{2,1}(N) = 0.8916175219 > e^γ/2 = 0.8905362091 |

Inherited, not re-derived here: RobinCH's analytic lemmas (9–13, 11_k; rerun by STCAudit and HP141, no line-by-line referee audit in the repo — the natural next audit) and inputs (RH to 3·10^12, Chirre–Helfgott ψ bound arXiv 2512.15709 — checked by CHAudit, Büthe θ < 10^19); WY's analytic part of Prop. 3.5 (x_1 ≥ 120409, under RH; constants spot-checked: wyconst.log); the half-plane inputs of the 46.79 rows.

## Reruns (copies in .scratch/oddraudit/run, one process at a time, memcap 290 MB + /usr/bin/time -l; rerun.log)

| run | result (= check_all.log?) | max RSS |
|---|---|---|
| gap.py | ALL PASS, 10 vectors, 3 in the gap (=) | 40 MB |
| dfs_gap.py | 1,098 n, largest log-margin −3.7812489·10^−4 at 3N (=) | 42 MB |
| oddca.py 1e6 / 5e7 | PASS (9,882 + 68,906; 10,972 + 2,991,542) | 43 / 110 MB |
| oddca_seg.py 5e8 (oddca.py with a preallocated prime buffer, no other change) | PASS, 13,236 exact + 26,346,275 in 79,866 blocks, log N_last = 500016054.5165 (=) | 168 MB |
| transfer.py 5e8 1.43e47 0.87338196661320752 | PASS, 1209 boxes, min 1.07·10^−5, 10^(10^46.7931) (=) | 43 MB |
| transfer.py 5e8 1.42e47 7/8 | PASS, 1235 boxes, 10^(10^46.7900) (=) | 44 MB |
| transfer.py 5e8 6.46e46 1 | PASS, 1221 boxes, 10^(10^46.4480) (=) | 43 MB |

Largest x_1 rechecked under the cap: 5·10^8 (the whole §3 range). The original 0.42 GB peak came from oddca.py's np.concatenate copies; oddca_seg.py gives identical counts and margins.

## Independent recomputations (no arb; mpmath interval arithmetic, γ hard-coded to 60 digits; values.log, indep.log)

- σ(N) = 1015178899037727378898944000 (multiplicativity; N's factorisation re-derived by trial division), σ(N)/N = 3.66366513980014; δ0 = 0.00121349015854 > 0 (N violates); log-margin(c0) = −9.23561414276·10^−5 (c0 satisfies); ε* = 0.00352346464443 < 1/(t_c0 log t_c0) = 0.00367908685049; c_73(1) = 0.00151163009349 > δ0. All as in PROOF.md.
- ε-cost enumeration redone (exponents 0..59): options within δ0 only at p = 3 {4, 5}, 5 {2, 3}, 11 {1, 2}, 67 {0, 1}, 71 {0, 1}; of the 32 combinations exactly 10 have cost ≤ δ0; in the gap 3N (−3.7812489·10^−4), 11N (−5.9604394·10^−4), 67N/5 (−5.0950437·10^−4), all certified < 0. Identical to PROOF.md l. 32–41.
- Gap by a third method (enum_s.py: best-first search over COMPLETE odd exponent vectors, exponent 0 a branch, exact integer bound, exact rational T0 ≤ (e^γ/2) log log N; no ε-cost, no arb; 5,539 nodes, 0.7 s, 40 MB): S = {odd m < c0 : σ(m)/m ≥ T0} has 1,099 elements: N itself and 1,098 in (N, c0) — the same set as dfs_gap.py's list, element by element, same σ. Every one certified (interval) to satisfy (OR) except N. Byproduct: N is the only odd m ≤ N with σ(m)/m ≥ (e^γ/2) log log N.

## Falsification (falsify.log, interp.log) — nothing found

- Neighbours of N: N p^k/q^j and N p r/q, odd p, q, r ≤ 200, k ≤ 3, j ≤ 2: 19,333 n > N; largest margin −9.2356·10^−5 (c0 itself).
- Neighbours (p, q, r ≤ 100) of all 309 odd CA numbers from c0 to x_1 = 2000: 2,820,131 n > N; largest margin −9.2356·10^−5 (again c0, as (71·c0)/71).
- Odd superabundant records in (N, c0): 7 (2.2N, 3N, 6.6N, 11N, 13.4N, 33N, 40.2N), largest margin −3.78·10^−4 (3N).
- §2 numerically: odd CA sequence rebuilt in mpmath (332 events up to x_1 = 2000, min relative gap 1.0·10^−4, no ties; N = event (61,1), next event (67,1) = c0); g(N_k) = g(N_{k+1}) at each tie to 10^−50 (331 pairs); the line bound log(σ(n)/n) ≤ ℓ(log n) holds for every odd n ≤ 2·10^6 and every n in S; all 273 odd violators ≤ 2·10^6 and N satisfy Lemma 3.4's conclusion. Every odd CA number below c0 violates (OR) (c0 is the first that satisfies), so "both endpoints satisfy" first occurs at [c0, next].
- Inherited RobinCH certificate against the truth at its weakest end (rtest.log; oddca_seg.py extended to x_1 = 5.06·10^8: PASS, 13,260 + 26,645,237 numbers, log N_last = 506006330.79, 179 MB): at five CA numbers M = 2^32 N_ε with x_1 ∈ [5.0095·10^8, 5.0500·10^8], true (√x log x) R(M) = 0.7554–0.7560, certified m = 0.50 ± 0.006 (b0 = 1 and 7/8): m < truth, consistent.
- Lemma O2 against the exact α_2 at x = 41 … 1.43·10^47 (corr_check.log): ⌊log_2 x⌋ ≤ α_2 ≤ A/log 2 at every sample; normalised odd correction 0.0029 at 5·10^8, ≤ 10^−8 beyond 10^20.

## Findings (none affects the claim)

- F1. WY numbering (TeX l. 31–40, one counter per section): the main theorem MainThm1 (l. 87) is Theorem 1.2 in arXiv v1; LEAN.md calls it "Theorem 1.1". Lemma 3.4, Prop. 3.5, Thm 3.6 as in PROOF.md ✓.
- F2. WY l. 153 says ties take "the smaller power", but its formula (l. 144) takes the larger (at ε = F(p,a) the floor argument is exactly a+1). PROOF.md l. 9 follows the formula and covers both maximisers; harmless.
- F3. WY l. 399 "satisfied for all odd CA numbers with 61 < x_1 < 120409" is false for real x_1 ∈ (61, 67) (there N_ε = N violates); read "largest prime > 61". Our §3 recheck replaces WY's computation anyway.
- F4. PROOF.md l. 70 "Robin margins ≥ 0.49 there" depends on the box width: RobinCH's m is 0.468 on log x ∈ [20.03, 20.04], 0.50 on log x ± 0.001 at x ≈ 5.0·10^8 (the true value is 0.756, rtest.log); against a correction of 0.0029 immaterial.
- F5. PROOF.md l. 79 proves (⇒) via Prop. 3.5 + §§1–3; shorter and equivalent: (b) = WY Thm 1.2 + §1 (the gap) alone. "Large margin" (l. 88) means the sign of −0.7702 + 7.1476/log x, negative for x ≥ 10724 (wyconst.log: 2 + c − √2 − 0.99154√2 = −0.77027, α(120409) = 2.66160 ≤ 2.6616).
- F6. check_all.py's 0.42 GB peak: oddca.py's np.concatenate + boolean-mask copies of the 26 M primes; a preallocated buffer (oddca_seg.py, 8 changed lines) halves it with identical output.

Suggested edits (orchestrator's call, OddRobin not touched): LEAN.md "Theorem 1.1" → "Theorem 1.2 (arXiv v1)"; check_all.py: oddca.py's sieve → oddca_seg.py's preallocated buffer (fits the 300 MB cap); PROOF.md l. 70 "≥ 0.49" → "≈ 0.47–0.50 (box-dependent; true 0.756)".

## Reviewer

reviewer (verdict-only, §4 Lemmas O1–O2 + the use of RobinCH's margins, and the §2 interpolation step, exact values pasted): CORRECT ("the strict Arb test m − corr > 0 proves R(M_ε) > C_odd throughout the box"; interpolation by concavity of h − ℓ).

## What survives

All of discoveries/odd-robin.md as written. Weakest links are inherited: RobinCH (above 5·10^8) and, for the 46.79 rows, the half-plane hypothesis; the gap, the RH equivalence with its sharp threshold N, and the range e^(5.0·10^8) rest on arithmetic, arb/mpmath intervals and WY's Theorem 1.2 only.

Files: FULLPROOF.md; rerun.log, indep.log, values.log, interp.log, falsify.log, corr_check.log, wyconst.log; rtest.log; code/ (enum_s.py, values.py, interp.py, falsify.py, corr_check.py, wyconst.py, oddca_seg.py, rtest.py, dfs_dump.py, chain.sh; run them from a copy of OddRobin/ + RobinCH/ certify*.py in one directory, transfer.py's sys.path line pointed there).

RH relevance: N4


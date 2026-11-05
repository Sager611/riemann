RECORD: every odd n with 277093800961598968947975375 < n ≤ 10^(10^46.44) satisfies σ(n) < (e^γ/2) n log log n (finite gap of arXiv 2008.04787 closed; inputs RH to 3·10^12 + Chirre–Helfgott; 10^(10^46.79) with the Lean-checked half-plane), and RH ⟺ this holds for every odd n > 277093800961598968947975375 (sharp: that n fails)

# Odd Robin inequality: the gap N < n < c0 of Washington–Yang closed (node N4, side result)

Agent oddrobin, 2026-11-05 (21:32–22:20 UTC). Files in this folder:
- PROOF.md: the full argument (gap lemma, interpolation, direct check, transfer lemma O1/O2, theorem, residual risks).
- gap.py: the gap by an ε-cost enumeration (10 exponent vectors, 3 in the gap), arb.
- dfs_gap.py: the gap again, independently, by exhaustive branch-and-bound over all odd n < c0 (exact integers), arb.
- oddca.py: every odd colossally abundant number from c0 up to largest prime 5·10^8, arb + exact integer products.
- transfer.py: odd colossally abundant numbers with largest prime in [5·10^8, 1.43·10^47] from our Robin certificate (Riemann/Deep/RobinCH) plus an explicit odd correction.
- check_all.py: reruns all of the above; check_all.log: its output (ALL PASS).

## 1. Short answer

1. **The gap is closed, and it needs no RH.** Washington–Yang (arXiv 2008.04787; Int. J. Number Theory 17 (2021) 843–870) proved, assuming RH, the odd analogue σ(n)/n < (e^γ/2) log log n for odd n ≥ c0 = 3^4·5^3·7^2·11···67 = 18565284664427130919514350125, and left open (TeX l. 400) the odd n with N < n < c0, N = c0/67 = 277093800961598968947975375 (N itself fails). Those n are finitely many but up to 1.9·10^28, so brute force is out. A cost identity at the ε where N and c0 tie shows any failure in the gap must deviate from N by a total "cost" below N's own excess δ0 = 0.0012135; only 10 exponent patterns qualify, 3 lie in the gap (3N, 11N, 67N/5), and all 3 satisfy the inequality (closest: 3N, log-margin −3.78·10^−4). An independent exhaustive search (1,098 odd n in the gap with σ(n)/n above N's threshold) agrees.
2. **Sharp RH criterion.** With Washington–Yang's RH argument for large numbers: **RH ⟺ σ(n) < (e^γ/2) n log log n for every odd n > 277093800961598968947975375**, and the bound is best possible (n = N fails). This is the exact odd counterpart of Robin's "n > 5040".
3. **Range without RH.** Every odd n with N < n ≤ 10^(10^46.44) satisfies it (RobinCH's inputs: RH to height 3·10^12, Chirre–Helfgott, Büthe; 10^(10^46.79) if one adds the half-plane Re ρ ≤ 0.87338… or 7/8 for zeros above 3·10^12): odd colossally abundant numbers are checked directly up to largest prime 5·10^8 (26.4 million numbers), and beyond that by our Robin certificate (RobinCH) plus a correction of size ≈ 2/x for removing the power of 2 (Lemma O1/O2). With no analytic input at all (arithmetic only): up to e^(5.0·10^8).
4. **Their question about f_{2,1} (TeX l. 85).** Under RH, f_{2,1}(n) = σ(odd part of n)/(n log log n) stays below its limsup e^γ/2 for every n > N, and N is the last n where it does not (corollary in PROOF.md §5(c)).

## 2. Statement and before/after

| | before | now |
|---|---|---|
| odd n in (N, c0) | open ("we do not know", WY 2020) | all satisfy (OR), unconditionally (pure arithmetic) |
| RH-equivalent form | RH ⟺ (OR) for odd n ≥ c0 | RH ⟺ (OR) for odd n > N; sharp (N fails) |
| range without RH | odd colossally abundant numbers with largest prime ≤ 10^6 checked by computer by WY (precision not stated), i.e. c0 ≤ n ≤ e^(≈10^6) | N < n ≤ e^(5.0·10^8) by arithmetic alone; N < n ≤ 10^(10^46.44) with RH to 3·10^12 + Chirre–Helfgott; 10^(10^46.79) adding the half-plane (certified, arb) |
| last odd colossally abundant failure | N (WY) | N (arb, rechecked) |

(OR): σ(n)/n < (e^γ/2) log log n.

## 3. Method

- **Gap (gap.py, PROOF.md §1).** ε* = log(68/67)/log 67, where N and c0 both maximise σ(m)/m^{1+ε} over odd m. For odd n: log(σ(n)/n) = log(σ(N)/N) + ε*(log n − log N) − L(n), L(n) = Σ_p c_p(a_p) ≥ 0 a sum of per-prime costs (concavity in the exponent). n fails iff L(n) ≤ D(log n), and D decreases on the gap from δ0 = 0.0012135. Every prime ≥ 73 costs ≥ 0.0015116 > δ0; for p ≤ 71 only 3^5, 5^2, 11^2, 67, 71 are cheap enough. 10 patterns, 3 in the gap, all pass.
- **Independent gap check (dfs_gap.py).** Branch and bound over odd n < c0 with exact integers and the bound σ(k)/k < ∏_{p|k} p/(p−1) over the fewest admissible primes: 35,792 nodes, 1,098 odd n in the gap with σ(n)/n ≥ (e^γ/2) log log N, all satisfy (OR) at their own n.
- **Odd colossally abundant numbers to largest prime 5·10^8 (oddca.py).** Events (p, a) sorted by ε = log(1 + 1/(p + … + p^a))/log p; merge positions and near-ties certified in arb; below p = 10^5 each number checked exactly, above in blocks (S = log σ/n and T = log N rise along a run of new primes, so S_after vs log log T_before covers the block; exact integer products). 26,359,511 numbers, last failure N, all from c0 on pass.
- **Transfer (transfer.py, Lemma O1/O2).** N_ε = M_ε/2^a with M_ε colossally abundant: (OR) at N_ε ⟺ Robin's log-margin R(M_ε) exceeds C_odd ≤ 2/x + log(((x+1) log x/log 2 + 2)/2)/(L_N log L_N), L_N = x(1 − 1/log x) − log 2. RobinCH's Level 3 boxes give (√x log x) R(M_ε) ≥ m; all boxes keep m above the normalised correction (≤ 0.0029 at x = 5·10^8, against m ≥ 0.49 there; smallest remaining normalised margin over all boxes 1.07·10^−5, the same in all three runs).
- **Interpolation (WY Lemma 3.4).** Between consecutive odd colossally abundant numbers that both satisfy (OR), every odd n does.

## 4. Certificates and rerun

One command, 5.7 minutes, 0.42 GB, one core (final run 22:12–22:18 UTC, ALL PASS):

    cd Riemann/Deep/OddRobin && uv run --no-project --with python-flint --with mpmath --with numpy python3 check_all.py

Output (check_all.log):

| part | result |
|---|---|
| gap.py | ALL PASS (δ0 = 0.00121349015854, c0 margin −9.2356·10^−5; 3 gap candidates, margins −3.78·10^−4, −5.10·10^−4, −5.96·10^−4) |
| dfs_gap.py | PASS (1,098 candidates, closest 3N) |
| oddca.py 5e8 | PASS (13,236 one by one + 26,346,275 in 79,866 blocks; last failure N) |
| transfer.py, no half-plane | PASS, 1221 boxes, x_1 ≤ 6.46·10^46 ⇒ all odd N < n ≤ 10^(10^46.4480) |
| transfer.py, + Re ρ ≤ 7/8 above 3·10^12 (conditional) | PASS, 1235 boxes, x_1 ≤ 1.42·10^47 ⇒ 10^(10^46.7900) |
| transfer.py, + Re ρ ≤ 0.87338196661320752 above 3·10^12 (conditional) | PASS, 1209 boxes, x_1 ≤ 1.43·10^47 ⇒ 10^(10^46.7931) |

## 5. Inputs and tradeoffs

| part | inputs | status |
|---|---|---|
| gap (N, c0) | none (arithmetic + arb) | two independent certificates |
| c0 ≤ n ≤ e^(5.0·10^8) | none (arithmetic + arb) | certificate |
| up to 10^(10^46.44) | RobinCH: RH to 3·10^12 (Platt–Trudgian), Chirre–Helfgott pointwise ψ bound (preprint 2512.15709), Büthe below 10^19 | certificate |
| up to 10^(10^46.79) | the above + half-plane Re ρ ≤ 0.87338… (Lean-checked per QuasiRiemannTracker, axioms unchecked) or 7/8 (OpenAI 2026, Lean port) | certificate, conditional on the half-plane |
| all odd n > N | RH (WY Prop. 3.5 for largest prime > 1.43·10^47) | conditional, as Robin's criterion |

## 6. Novelty pass (2026-11-05)

- arXiv listings math.NT 2026-10 and 2026-11 (python3 .scratch/lst.py 'Robin|odd|sigma' math.NT --zenodo 'odd Robin inequality'): only 2610.12816 (Robin's criterion for reciprocal LCM powers; unrelated). Zenodo: nothing relevant.
- GXL evidence API (3 claim queries, 240+ passages each): the only papers on the odd inequality are WY 2008.04787 itself and Koya, arXiv 2202.08390 (2022), who proves the unconditional second-type bound σ(n)/n ≤ (e^γ/2) log log n + 0.7398…/log log n for all odd n ≥ 3 and quotes WY's threshold c0 unchanged; Fan–Kobayashi–Molnar arXiv 2511.02106 (2025) cite WY (IJNT version) for σ^{[k]} analogues, not the gap. No source states the gap closed or N as the sharp threshold.
- arXiv API and Semantic Scholar citation queries were rate-limited (HTTP 429) during this run; the published IJNT version was not read (paywalled), but Koya (2022) and the arXiv v1 (only version) both state the threshold c0.
- Our NOVELTY/discoveries files: robin-nicolas-ranges.md covers Robin and Nicolas only.

## 7. Review

Reviewer (verdict-only, 2026-11-05, on PROOF.md + gap.py): "FLAWED as written; the finite-gap argument itself is correct." §§1–2 (cost identity, monotonicity, candidate reduction, interpolation): sound. Lemma O1–O2 and the transfer sign: correct. §5(b) (sharp RH equivalence) and the even-n reduction of §5(c): correct. Two wording flaws, both fixed: (1) 10^(10^46.79) was called unconditional although its rows assume the half-plane, so the headline now uses 10^(10^46.44) and labels 46.79 as conditional; (2) the tie convention (larger maximising exponent, WY's floor formula) is now stated where p | N_ε ⟺ p ≤ x_1 is used.

## 8. Residual risks

- The range above 5·10^8 inherits RobinCH's certificate (reviewed, rerun ALL PASS) and the unrefereed Chirre–Helfgott bound; the gap and everything up to e^(5.0·10^8) use none of it.
- Float sorting of odd colossally abundant events: every near-tie and every merge position is certified in arb; F(p, a) is evaluated in log1p form (no cancellation).

RH relevance: N4 (Robin-type inequalities; side result).

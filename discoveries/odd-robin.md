# The odd Robin criterion made sharp (node N4)

**Claim.** Every odd n with 277093800961598968947975375 < n < c0 = 3^4·5^3·7^2·11···67 = 18565284664427130919514350125 satisfies σ(n) < (e^γ/2) n log log n (unconditionally, by arithmetic and arb). Hence RH ⟺ σ(n) < (e^γ/2) n log log n for every odd n > N = 277093800961598968947975375, and the threshold is sharp (N fails). Unconditionally the inequality holds for every odd n with N < n ≤ 10^(10^46.44) (RH to height 3·10^12 + Chirre–Helfgott through our Robin certificate RobinCH; 10^(10^46.79) with the Lean-checked half-plane), and up to e^(5.0·10^8) with arithmetic alone. Under RH, f_{2,1}(n) stays below its limsup e^γ/2 for every n > N (Washington–Yang's question at TeX l. 85).

**Before.** Washington–Yang (arXiv 2008.04787; Int. J. Number Theory 17 (2021) 843–870) proved RH ⟺ the inequality for odd n ≥ c0 and left the odd n in (N, c0) open ("we do not know", TeX l. 400). Koya (arXiv 2202.08390) quotes c0 unchanged.

**Method.** At the ε where N and c0 tie as odd colossally abundant numbers, every odd n carries a nonnegative cost (a sum of per-prime costs of deviating from the CA exponents); a violation inside the gap forces cost(n) < δ0 = 0.0012135, N's own excess. Only 10 exponent patterns qualify, 3 in the gap (3N, 11N, 67N/5), and all 3 satisfy the inequality (closest: 3N, log-margin −3.78·10^−4). An independent exhaustive branch-and-bound over odd n < c0 (1,098 odd n in the gap above N's σ(n)/n threshold) agrees. The range: a failure forces one at an odd CA number; those are checked directly up to largest prime 5·10^8 (26.4 million numbers, arb) and beyond by RobinCH plus an explicit correction for removing the power of 2 (Lemmas O1, O2).

**Certificate.** Riemann/Deep/OddRobin/check_all.py (gap.py, dfs_gap.py, oddca.py, transfer.py; python-flint arb and exact integers): ALL PASS in 5.7 minutes, 0.42 GB, one core; rerun by the orchestrator on a copy: identical.

**Review.** Agent-side reviewer: argument correct after two wording fixes (the 46.79 range is conditional on the half-plane; the tie convention stated). Orchestrator's reviewer (verdict-only, PROOF.md): CORRECT as a deduction from the certified bounds; SIGNIFICANT (closes the stated published gap, sharp cutoff); the arb computations and the RobinCH certificate were not audited by that reviewer.

**Novelty.** arXiv math.NT listings 2026-10 and 2026-11 (lst.py): only 2610.12816 (Robin's criterion for reciprocal LCM powers, unrelated); Zenodo "odd Robin": nothing; HAL: nothing; GXL evidence API (agent): only Washington–Yang and Koya treat the odd inequality, Fan–Kobayashi–Molnar (arXiv 2511.02106) cite it for σ^[k] analogues, not the gap; our files: robin-nicolas-ranges.md covers Robin and Nicolas only. Citation databases were rate-limited; the IJNT version was not read (paywalled).

**Caveat.** Unrefereed; the range above largest prime 5·10^8 inherits RobinCH's certificate and the unrefereed Chirre–Helfgott bound (the gap and everything up to e^(5.0·10^8) use neither). Lean (oddrlean, 12 Nov 2026): the gap theorem, the sharpness at N and the equivalence are machine-checked from Washington–Yang’s theorem (named) and Mathlib’s log/γ facts used through an abstract interface; all numerics by norm_num, sorry-free, standard axioms, reviewer FAITHFUL; the final instantiation with Mathlib’s Real.log is not compiled because its imports exceed the 300 MB cap ([LEAN.md](../Riemann/Deep/OddRobin/LEAN.md)).

**RH relevance.** N4 (Robin-type RH criteria): an RH-equivalent with a sharp explicit threshold, the odd counterpart of Robin's n > 5040.

**Files.** [OddRobin](../Riemann/Deep/OddRobin/REPORT.md) (REPORT.md, PROOF.md, check_all.py and its parts, check_all.log).

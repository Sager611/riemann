AUDIT PASSED: Theorem A (a prime in (x, x + x^{0.856}] for every real x ≥ 1) holds as written on its stated inputs — every step that is new relative to the audited Ingham3 pipeline checks line by line, the θ = 0.856 certificate reruns identically, and its worst cell recomputes independently at 0.9509 ≤ 0.9711.

# Ingham4Audit — line-by-line referee report on Theorem A (ingham4audit, 20 Nov 2026)

Agent ingham4audit (Orca, node N10). Read-only web (arXiv e-prints via curl -o: KLN 2101.12263, Cully-Hugill–Johnston 2111.10001; deleted afterwards); nobody contacted; no Lean, no commits. Audited: discoveries/primes-short-intervals.md (Theorem A), Riemann/Deep/Ingham4/ (DERIVATION §1 (I4′), §2, §3′, §5, §6; cert.py part main), against Ingham3Audit (whose checks of the common pipeline transfer: FULLPROOF §0) and R501. FULLPROOF.md has the chain with DERIVATION sections and KLN TeX lines.

## Findings
- F1 (KLN, harmless, unstated). KLN's argument bound (TeX l. 1439–1448) turns "(…)log X + (…)log T + b₉" into "log T·[(…) + b₉/log H₀]" using H₀ ≤ T; this needs b₉ ≥ 0, which KLN do not state. True: b₉ = 48.17 at both heights (kln_indep.py).
- F2 (prose gap in (I4′), not in the code). DERIVATION's reading "H₀ enters only through RH(H₀), T ≥ H₀ ≥ H, X = kT ≥ kH₀ ≥ 10^9, log H₀, log log H₀ and δ < log H₀(log log H₀ − 1)/2" is correct for the proof (all 94 H₀-lines classified, FULLPROOF §1) but omits that the row conditions U(α,k,H₀) > 1, J-coefficients > 0 (U decreasing) and V(α,k,δ,∞) > 1 must be re-checked at the new H₀. cert.py does re-check them (KP.ok), which is why row 0.60 is dropped.
- F3 (wording). KLN's Lemma 4.14 fixes H₀ = 3.0610046·10^{10} in its statement (l. 1550). Theorem A uses KLN's proof re-run with H₀ = 3·10^{12} (as Cully-Hugill–Johnston, refereed, did: 2111.10001 main.tex l. 734), not the lemma as printed. Suggested input wording: "KLN's proof of Lemma 4.14 with H₀ = 3·10^{12}".
- F4 (cosmetic). TrapKernel's far branch (λ₀ ≥ max(2π/a, 63.93)) is valid but 30% loose on [62.8, 92] (C₂/λ² is attained only beyond 2/(as₁) = 92.0); it is not used near the binding cells (λ₀ = 3.02 there).
- No false step found: trapezoid estimates, G_a monotone and ≤ C₂/λ², moment tables and the three lookup branches, fup, KLN at RH(3·10^{12}), real-x weight/explicit formula/criterion, prime-power bound with m = 2 ≥ 10/9, δ-cells (i)–(iv), coverage of [42, 300], the tail lemma and its monotonicity, small x.

## Verdict per task
| task | verdict | evidence |
|---|---|---|
| (1) §3′ kernel, TrapKernel | HOLDS | both estimates re-derived; 10 stored moments at 256 bits by independent quadrature: stored/definition 1.0006–1.0018 (all ≥ 1); 10 stored G(λ_i) ≥ G_a; mom() branches vs definition 1.0015–1.30; 9000 random z: max |Φ̃_a|/(ab G_a) = 0.99999999938 (kernel_check.py) |
| (2) (I4′) KLN at RH(3·10^{12}) | HOLDS (F1–F3) | 94 H₀-lines read; own 256-bit implementation from the TeX: Table 1 at KLN's H₀ within 4.4·10^{−4}; at 3·10^{12}: row 0.60 dropped (U = 0.976 < 1; 1.029 at KLN's H₀); row 0.85: A = 8.51220384, B = 3.377370681; cert.py agrees to 10^{−37} on 19 rows × 4 σ (kln_indep.py); reviewer CORRECT |
| (3) §2 for real x | HOLDS | re-derived (FULLPROOF §4); R501: PP ratio ≤ 0.2053 |
| (4) §5 cells, tail | HOLDS | (i)–(iv) and every tail step re-derived; 4932 contiguous cells, 42.0 → 300.0 exact; R_an(300) = 0.014177 independently (cert 0.014196); all monotonicity conditions true; R_an non-increasing to L = 10^4 (tail_small.py) |
| (5) small x | HOLDS | 1476^{1000} ≤ 9999991^{856}; 2e^{42} = 3.48·10^{18} < 4·10^{18} (Bertrand); log 4/log 7 = 0.712414; R501 sieve |
| (6) certificate, worst cell | HOLDS | cert.py part main on a byte-identical copy: results identical (111 s, memcap 89 MB, max RSS 109 MB); worst cell [191.82, 191.84] by direct σ-integration of KLN's F(σ) with the true envelope, best row at every σ, 256 bits: 0.950941 ≤ 0.971108 (worstcell.py) |

## What survives if an input fails (from DERIVATION §8; not rerun here)
- KLN's proof not accepted at H₀ = 3·10^{12} (only Lemma 4.14 as printed): pipeline A0 (triangle, KLN at 3.0610046·10^{10}) certifies x^{0.857} for every real x ≥ 1.
- KLN Theorem 1.1's "δ ≥ 1" taken literally (I3A F3: misprint): x^{0.861} (A1 with δ = 1 in every row).
- Bebel's half-plane (I1) false: nothing here survives (HPAudit: AUDIT PASSED, unrefereed).

## LiDens (Li-range record 110.4757 T², 20 Nov)
- Task (2)'s verdict covers LiDens's KLN-at-RH(3·10^{12}) input: LiDens REPORT l. 47–48 takes A(σ), B(σ) from Ingham4/cert.py's KP and build_kln_H0('3e12'), the 19 rows audited here; LiDens's own σ-ranges per row (each needs σ > 1/2 + d/log H₀ ∈ [0.51152, 0.51183]) were not rechecked.

## Optional: θ₁₄₁ (TREE N7)
- θ₁₄₁ = 0.8733458561… in place of B only trims the σ-integral by 3.6·10^{−5} at the top (estimated < 0.01 of Z near log x = 190.6, from the band breakdown of DERIVATION §9): θ = 0.855's peak 1.152 stays > 1.14, so the certified exponent stays 0.856.

## Recomputations (logs/)
- cert.py part main on a byte-identical copy (sha1 1a25901e…): main, main_minus, main_minus_fine == Riemann/Deep/Ingham4/logs/cert_result.json; 111 s, memcap 89 MB, max RSS 109 MB (rerun_main.log, rerun_part_main.json, rerun_main_footprint.txt).
- kln_indep.py: KLN Lemma 4.14 from the TeX at 256 bits, I(A,n) by quadrature; Table 1 (KLN's H₀) within 4.4·10^{−4}; at 3·10^{12}: 19 admissible rows, row 0.60 fails U > 1 (0.976213); cert.py equal to 1.1·10^{−37}; 68 s, memcap 18 MB, max RSS 44 MB (kln_indep.log).
- kernel_check.py: envelope, monotonicity, C₂/λ², 10 moment entries + 10 mom() calls at 256 bits; 23 s, memcap 52 MB, max RSS 98 MB (kernel_check.log).
- worstcell.py (+ kln_core.py): the worst cell by direct σ-integration at 256 bits: 0.950940998468 (λ-tails two ways, identical to 12 digits; 480 s and 166 s; memcap ≤ 23 MB) (worstcell_test.log, worstcell.log); with 6 instead of 3 σ-nodes per cell: 0.950941031702 (Δ = 3.5·10^{−8}; 261 s; worstcell_ns2.log).
- tail_small.py: A_max, B_max, M*, M₁ independently; R_an(300) = 0.014177 (cert 0.014196); monotonicity conditions; R_an(L) for L up to 10^4; small-x integers; 116 s, memcap 19 MB (tail_small.log).

## Reviewer
- One reviewer subagent on the most delicate step (KLN's proof at H₀ = 3·10^{12}: the H₀ inventory, the formulas of C₁, C₂, M, b₁₂ and every monotonicity step pasted inline): CORRECT ("U(T)^{e(T)} ≤ U(H₀)^{e(T)} ≤ U(H₀)^{e(H₀)} … does not require U(T) > 1"; "drop the σ₀ = 0.60 row"; "the rectangle, Jensen and Gabriel steps introduce no additional height restriction").

## Files
REPORT.md, FULLPROOF.md; kln_indep.py, kln_core.py, kernel_check.py, worstcell.py, tail_small.py (kln_indep.py and kernel_check.py import Riemann/Deep/Ingham4/cert.py read-only with sys.dont_write_bytecode = True; the runs logged here used a byte-identical copy). Run each as: bash .scratch/memcap.sh 290 900 uv run --no-project --with python-flint --with mpmath --with sympy python3 Riemann/Deep/Ingham4Audit/<script> (worstcell.py takes NG NS, e.g. 3 1); logs/ (rerun_main.log, rerun_part_main.json, rerun_main_footprint.txt, kln_indep.log, kernel_check.log, worstcell_test.log, worstcell.log, worstcell_ns2.log, tail_small.log, *_footprint.txt).

RH relevance: N10

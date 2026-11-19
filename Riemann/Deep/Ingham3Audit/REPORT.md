AUDIT PASSED: (a) and (b) hold as written on the four stated inputs — no false step in the mathematics, and the one certificate gap (k = 5: the scan started at the double 288.73000000000001819, leaving log x ∈ [288.73, 288.73 + 1.8·10⁻¹⁴) unchecked) is closed here by one extra cell (R = 0.99628 < 1 at 128 and 256 bits).

# Ingham3Audit — referee report on "a prime between consecutive seventh powers for every n; explicit thresholds for k = 3–6" (pbpaudit, 14 Nov 2026)

Agent pbpaudit (Orca). Read-only web (arXiv e-prints via curl -o); nobody contacted; no Lean, no commits. Audited: discoveries/primes-between-powers.md (claims (a), (b)), Riemann/Deep/Ingham3/ (DERIVATION.md, cert.py, logs). Not re-audited: Bebel's half-plane (Riemann/Deep/HPAudit, AUDIT PASSED), KLN's proof beyond the statement level, Ramaré's b₁–b₄. FULLPROOF.md has the chain with TeX line numbers.

## Findings
- F1 (certificate gap, k = 5, closed here). find_L0 returns the double nearest 288.73, which is 288.73000000000001819 > 288.73; the scan starts there. The claim "n ≥ e^{57.746}, i.e. log n⁵ ≥ 288.73" therefore had the sliver log x ∈ [288.73, 288.73 + 1.8·10⁻¹⁴) unchecked (≈ 4·10¹⁰ integers n). The other starts are below their decimals (217.449999…, 468.839999…, 1881.029999…). Closed by one extra cell (Recomputations). Same pattern in HP141's refresh at k = 4 (double 468.68 > 468.68 while the statement is n ≥ e^{117.17}, log x = 468.68): closed by the same one-cell check (section HP141 theta141 slivers).
- F2 (KLN slip, harmless). KLN Lemma 3.2 justifies a₂ = 2.851 by "a₁t^{1/6}log t + a₂ ≥ 1.461"; at t = e⁻⁶ the left side is 1.460416 < 1.461. The printed statement (3.5) is still true: checks2.py proves |ζ(1/2+it)| ≤ min_{T≥t}(0.63T^{1/6}log T + 2.851) on [0, 3] (Arb, margin 6.6·10⁻⁷; |ζ(1/2)| = 1.4603545).
- F3 (KLN misprint, hedged). Theorem 1.1 prints δ ≥ 1; everything else in KLN (l. 453, 626, 1099, 1549; Table 1 with δ = 0.303–0.314) has δ > 0, and none of the δ-conditions in the TeX (rg of all δ-inequalities) is δ ≥ 1. Lemma 4.14 is the statement used. If δ ≥ 1 were really needed: (a) still holds; (b) becomes log x₀ = 1947.12, 491.86, 305.60, 231.99.
- F4 (cosmetic). Ingham3's cell_bound does not assert L_a(1 − B) > 1 (needed for the x-monotonicity); true for every cell (L_a ≥ 32.236). cert.log says "every n ≥ 100 certified" for the cells; in fact n = 100 itself (log x = 32.2361913019166396 < double L₁ = 32.2361913019166408) is covered only by the direct search, which includes it.
- No false step found in the mathematics: weight, Mellin transform, explicit formula, trivial zeros, prime powers, kernel majorant, moments, pairing and σ₀ split, KLN at every σ (log C₁ affine in σ: KLN (4.72) l. 1569–1574), case majorants, σ- and x-cells, tail lemma, small n (FULLPROOF §2).

## Verdict per item
| item (discoveries/primes-between-powers.md) | as written | evidence |
|---|---|---|
| (a) k = 7, every n ≥ 1 | HOLDS | n ≤ 100: cert.py direct search; n ≤ 10⁴: FLINT-proven primes (checks1.py, largest p − n⁷ = 549); 3592 cells on [L₁ = 32.2362, 240], worst [205.52, 205.54]: 0.154624 (128 = 256 bits; direct 0.1520); tail R_an(240) = 0.2437 |
| (b) k = 3, n ≥ e^{627.01} | HOLDS | first cell 0.999891 (direct 0.9862); 2076 cells; tail 0.01171 |
| (b) k = 4, n ≥ e^{117.21} | HOLDS | 0.999672 (0.9879); 1712 cells; tail 0.01695 |
| (b) k = 5, n ≥ e^{57.746} | HOLDS with this audit's extra cell (F1) | sliver [nextafter(288.73, 0), 288.73…]: 0.996281 (0.9819); first cell 0.999567; tail 0.03683 |
| (b) k = 6, n ≥ e^{36.242} | HOLDS | 0.999679 (0.9817); 1651 cells; tail 0.08765 |
| inputs as quoted | HOLD at statement level | FULLPROOF §1: KLN Lemma 4.14 (δ > 0, σ > 1/2 + d/log H₀, T ≥ H₀), Table 1 transcription, a₁ via HPY 0.618, a₂ (F2), δ (F3); PT; HSW (its C₃ misprint is irrelevant) |
| "1/(1 − B) = 7.898 > 7" | HOLDS | 1/(1 − 0.8733820) = 7.8975 |

## Recomputations (logs/summary_cells.json; cells rebuilt exactly by mkcells.py, counts = cert_result.json)
- cert.py on a copy: cert_result.json byte-identical (logs/rerun_cert_result.json), cert.log identical but timings; 369 s, memcap 186 MB.
- 20 cells by cert.py's own code at 128 and 256 bits (cellcheck.py: module re-executed with ctx.prec changed): identical to 12 digits. Independent float quadrature (indep.py): certified ≥ direct in all 20; certified/direct = 1.011–1.022 for every cell with log x ≥ 200 (all four thresholds included), 1.04–1.42 for the k = 7 cells with log x ≤ 180, where R ≤ 0.0044.
- Tails at 256 bits: 0.24370 (k = 7, L* = 240), 0.08765 (6, 300), 0.03683 (5, 400), 0.01695 (4, 640), 0.01171 (3, 2400) = cert.log.
- KLN rows (rows.py, 256 bits): row 0.60 on [0.515, 0.6275] with s − s_min ≥ 8.77·10⁻⁴; … row 0.87 on [0.8675, B]; c₁ > 0, q₁ > 0, B_P(B) > 0 in all 20 rows, both δ runs.

## Falsification attempts (all failed to break it)
- Smallest n past the direct search: n = 101–103 lie in the first cell [32.2362, 32.4862]: R = 2.97·10⁻⁴; and directly (n ≤ 10⁴).
- Cell joins: [179.986, 180] 0.004383, [180, 180.02] 0.004400, [239.98, 240] 0.02068, then the tail; scans are contiguous double lists (no gaps, no overlaps needed).
- σ = 0.515 (row 0.60, s_min = 0.514123) and σ near B (row 0.87, margin 0.354): admissible; Theorem 1.1's δ ≥ 1 taken literally: (a) still holds (cert.py's δ = 1 run, rerun identical).
- Zeros just above 3·10¹²: N_σ(T₀) = 0 needs RH on (0, T₀] (PT: 0 < γ ≤ 3·10¹²); above T₀ only KLN's count enters. One zero at β = B, γ = T₀⁺ would add ≈ 7·10⁻¹² at the k = 7 worst point; KLN allows ~10¹⁰ such zeros and the bound already pays for them.
- Precision: 128 vs 256 bits identical; sup |sin y|/y check, moment grid and Mom lookup rounding re-read (FULLPROOF §2, Ingham4 CODE-AUDIT).

## What survives if an input fails
- KLN Theorem 1.1's δ ≥ 1 binding: (a) unchanged; (b) with log x₀ = 1947.12, 491.86, 305.60, 231.99.
- (I1) false: nothing here survives (Ford's region gives back exp(e^{32.90}) for cubes); (I1) is the half-plane audited in Riemann/Deep/HPAudit (unrefereed).
- KLN's internal inputs: b₃, b₄ (Ramaré, X ≥ 10⁹) not rechecked; a₂ holds by a margin of 6.6·10⁻⁷ (F2).

## Reviewer
One reviewer subagent on the core chain (pairing and σ₀ split, S via N(t) ≤ t log t, F(σ) via KLN at every σ with Stieltjes integration from T₀, both case majorants, σ-cell freezing, criterion), formulas pasted inline: CORRECT (its notes: the σ₀ split, vanishing boundary terms, row validity per cell, exponents < 2, and the trivial-zero bound ∫w/(u³−u) ≤ M₀/(x³−x) all check).

## Files
REPORT.md, FULLPROOF.md; checks1.py (Table 1 vs TeX, |ζ(1/2+it)| on [0,3] coarse, proven primes n ≤ 10⁴), checks2.py (KLN (3.5) on [0,3]), mkcells.py, cellcheck.py, indep.py, rows.py; logs/ (checks1.log, checks2.log, rerun_cert.log, rerun_cert_result.json, rerun_footprint.txt, cells*.json, cellcheck128/256.*, tail256.*, indep_part1.log, indep_part2.log, indep_rest.json, rows.log, summary_cells.json, footprints.txt).
Fix to adopt in Ingham3 (not applied: read-only): run_k should start the scan at a double ≤ the stated decimal (e.g. math.nextafter(L0, 0)) or the statement should quote the double; add assert L_a(1 − B) > 1 in cell_bound; HP141's slivers are closed below.

## HP141 theta141 slivers (follow-up, 14 Nov 2026)
HP141's θ₁₄₁ cert.py (B = 0.87334585610671858, C₃ = 9.4926), run as a byte-identical copy (.scratch/pbpaudit/hp141/cert.py) through cellcheck.py (4th argument = cert path), one zero-width cell [nextafter(L₀, 0), L₀] per k, L₀ = the first double of HP141's scan (mkcells_hp141.py, logs/cells_hp141.json). KLN Table 1 reproduced 20/20 in both runs.

| k | stated x₀ | first double of the scan | sliver | R, 128 bits | R, 256 bits |
|---|---|---|---|---|---|
| 6 | e^{217.37} | 217.37000000000000455 | 4.5·10⁻¹⁵ | 0.9960125143 | 0.9960125143 |
| 5 | e^{288.65} | 288.64999999999997726 | none (double below) | 0.9953716982 | 0.9953716982 |
| 4 | e^{468.68} | 468.68000000000000682 | 6.8·10⁻¹⁵ | 0.9903416409 | 0.9903416409 |
| 3 | e^{1880.45} | 1880.4500000000000455 | 4.5·10⁻¹⁴ | 0.9706120671 | 0.9706120671 |

All R < 1, so every sliver is closed: HP141's x₀ and its rounded-up statements n ≥ e^{626.82}, e^{117.17}, e^{57.73}, e^{36.23} hold from the stated decimals. Only k = 4's n-statement had a real sliver (4·117.17 = 468.68 = L₀); k = 3 and 6 round up by 0.01 in log x, and k = 5's double lies below its decimal. Peaks: 128 bits memcap 14 MB, max RSS 97.0 MB, 12.3 s; 256 bits memcap 83 MB, max RSS 120.4 MB, 17.1 s. Logs: logs/hp141_cellcheck{128,256}.{log,json}, logs/hp141_footprint{128,256}.txt.

RH relevance: N10

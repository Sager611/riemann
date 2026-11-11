NO COUNTEREXAMPLE: confined-pair law liminf (2π/log(g_n/2π))·μ(g_n) = 2, lower-scale law s·μ ≳ 2, trunk-slope law, and μ > 4 at confined pairs all survive finite tests on t ∈ [10^8, 5·10^11]. The trend runs against the liminf law: the per-window confined minima of s·μ are 1.70–2.06, below 2 in 56 of 64 windows (Arb-certified 1.6978 at n = 2,767,213,045), and their frequency below 1.948 grows 100-fold with t. The two Davenport–Heilbronn numerics are not reproducible: no tested section length N ≤ 59 gives both r_* = 0.72346 and μ = 1.45.

# JerbyConf: the finite laws of Jerby's 2610.21952 tested

Orca agent jerbyconf, 11 Nov 2026 (21:44–23:01 UTC). Every computation ran single-threaded and one at a time, except one 3 s awk extraction that overlapped a scan. Peak memory was 180 MB for jerby_data.py and at most 76 MB for everything else (/usr/bin/time -l). Statements are quoted exactly in [STATEMENTS.md](STATEMENTS.md), from the TeX (`Clouds.tex`, e-print of 15 Oct 2026) and the Zenodo data 23390338. Notation: s_n = 2π/log(g_n/2π), the paper's Δ_n (l. 2016–2021). **Confined** means the block (g_{n−1}, g_{n+1}] holds exactly two zeros of Z. This is the paper's own operational reading (l. 630–636: Rosser's rule holding is taken as confinement). The dynamical, path-dependent definition (l. 503–511) is checked separately in §4.

## Summary

| Statement | Test | Result |
|---|---|---|
| F1: μ > 4 at isolated bad points (l. 323–326) | jerbyfals | **false**: three Arb-certified zero-free blocks |
| F2: liminf s·μ = 2 over confined isolated bad points (l. 651–663) | minima and frequencies in 64 windows and 4 resonance scans; the author's own data; CUE model | not refutable by finite data. The trend moves **away** from 2: window minima 1.70–2.06 at every decade above 10^8 (2.88 and falling in the author's range); the frequency below 1.948 rises ×100; the CUE model (N = 64, 128) puts 1.7% of isolated bad points at confined s·μ < 2, with a tail down to 0.29 and no barrier at 2 (heuristic) |
| F4: pair term \|1/a_1 + 1/a_2\| ≥ 2, hence s·μ ≳ 2 (l. 2010–2051) | zeros located around the 16 lowest confined points | the pair term is 2.07–2.36 (≥ 2 holds). The heuristic omits one term: a zero just beyond the empty Gram interval, with the other near zeros, subtracts 0.12–0.59, giving s·μ = 1.70–1.81 |
| F3: \|Z′(g_n)\| ≍ A_act/(g_{n+1} − g_n) (l. 1826–1851) | ratio ρ on 6,516 low-viscosity points | in the author's range ρ ≈ 1 (0.86–4.7 over the 40 lowest). Above 10^8, ρ = 0.83–1129: no uniform constants |
| F5: DH collision parameter r_* ≈ 0.72346 (l. 778–783) | diagonal deformation of the section family, N = 2–59, exact and author's θ | not reproduced by any tested N: closest is 0.71305 (N = 9, 10) |
| F6: μ_DH(g_44^DH) ≈ 1.45 (l. 840–846) | same family and the exact D(s) | not reproduced: exact 2.178 (normalised Z_DH), 2.960 (ξ_DH); Zenodo table 2.7435. The section gives ≈ 1.45 only for N ≥ 21 with the Γ drift, where r_* = 0.40–0.47 |
| (iii) μ > 4 at confined pairs (l. 651–655) | all samples | holds: minimum confined μ = 4.8099 (n = 273,364,058, Arb) |
| Collision-avoidance reading of the bound (l. 594–600, 803–813) | argument principle on Jerby's family Z_N(t; r) | **fails at the extremes**: all 16 lowest confined points (μ = 4.8–7.0 > 4) lose both real zeros on an interval of r (a collision); 11 have two real zeros again at r = 1. Controls: 0 of 24 typical points at 10^9 collide, 5 of 200 in the author's range |

## 1. F2, the confined-pair law: data and trend

**The author's range** (Zenodo table, n ≤ 250,000; every isolated bad point there is confined, since the first Rosser violation is at n = 13,999,525). Minimum s·μ per range of n: 6.256 ([1, 10^3)), 3.705 ([10^3, 10^4), n = 8,569), 3.186, 3.003, 3.018, 2.890, then **2.881** ([200,000, 250,000), n = 237,330). It is falling, and no value is below 2.2 (logs/jerby_data.txt).

**64 log-uniform windows**: jerbyfals's scanner and windows, 28,990,518 Gram points, 3,076,101 isolated bad. Jerbyfals's CAND lines are complete for μ/θ′ < 0.62, i.e. s·μ < 1.948. The 21 windows without a confined point below 1.948 were re-run complete to s·μ < 2.513 (run_M.log, 874 s, 11 MB). Per-window table: logs/analysis_conf.txt.

| decade of t | windows | isolated bad | windows with confined min s·μ < 2 | range of window minima | lowest (window) | lowest confined μ |
|---|---|---|---|---|---|---|
| 10^8 (incl. L00 at 9.98·10^7) | 18 | 1,856,463 | 13 | 1.897–2.063 | 1.8974 (L14) | 5.295 |
| 10^9 | 17 | 825,673 | 17 | 1.714–1.994 | 1.7143 (L31) | 5.435 |
| 10^10 | 17 | 266,322 | 15 | 1.805–2.020 | 1.8051 (L45) | 6.346 |
| 10^11 | 12 | 127,643 | 11 | 1.716–2.025 | 1.7160 (L57) | 6.634 |

The window minima follow s·μ_min = 2.120 − 0.021 log₁₀ t (r = −0.30). They do not rise toward 2, even though samples shrink tenfold with height (10^6 down to 10^5 Gram points). The sample-size-free statistic is the frequency of confined s·μ < 1.948 per isolated bad point:

| log₁₀ t bin | 8.0 | 8.5 | 9.0 | 9.5 | 10.0 | 10.5 | 11.0 | 11.5 |
|---|---|---|---|---|---|---|---|---|
| isolated bad | 846,403 | 904,391 | 504,903 | 320,770 | 159,818 | 106,504 | 85,115 | 42,528 |
| confined s·μ < 1.948 | 1 | 12 | 18 | 14 | 5 | 8 | 11 | 7 |
| per 10^5 | 0.12 | 1.33 | 3.57 | 4.36 | 3.13 | 7.51 | 12.9 | 16.5 |

**Resonance scans** (jerbyfals R1–R3, t ∈ [1.02·10^8, 1.20·10^9], ≈ 3·10^6 isolated bad in checked regions): the confined minima of s·μ are 1.781 (R1), 1.748 (R2a), **1.698** (R2b, n = 2,767,213,045) and 1.741 (R3). Every confined CAND line in these scans has s·μ < 1.948: 135, 275, 347 and 234 points.

**Random-matrix (CUE) model** (cue_visc.py; tool transfer from random-matrix theory). Eigenangles of Haar unitaries play the zeros and the lattice 2πm/N plays the Gram points. Good means S_m = #{θ_j ≤ x_m} − m is even (ζ: (−1)^n Z(g_n) > 0 ⟺ S(g_n) even). The viscosity is (2π/N)|Σ ½cot((x_m − θ_j)/2)|. With N = 64 (24,000 matrices, 120,698 isolated bad) or N = 128 (1,500 matrices), 1.7–1.9% of isolated bad points are confined with s·μ < 2 and 0.06–0.09% with s·μ < 1. The quantiles 10^{-2}, 10^{-3}, 10^{-4} are 1.79, 1.10, 0.64, and the minimum over 116,368 confined points is 0.285. In this random-matrix model the confined lower tail shows no barrier at 2 and reaches 0.285. That is a heuristic, not a proof that the liminf is 0, which would need the recurrence of arbitrarily small values. ζ approaches it from below: its frequency below 1.948 is 1.2·10^{-6} at 10^8 and 1.6·10^{-4} at 10^11.5, against 1.5·10^{-2} in CUE. This is consistent with the slow convergence of ζ's local statistics to GUE. (The "random matrix Gram points" of 2610.39856 are the same lattice up to a configuration-dependent rotation, which is locally uniform.)

**Verdict F2.** A liminf cannot be refuted by finitely many values. But every decade from 10^8 to 5·10^11 has confined values below 2 (lowest 1.698, 15% below). The excursions grow more frequent with t, and the CUE model has no barrier at 2. The law's empirical basis (minima at 2.9 and falling for n ≤ 250,000) was extrapolated past the point where the data cross 2.

## 2. F4, the lower-scale heuristic: what it omits

gram_win.c -zeros locates every zero in [g_{n−6}, g_{n+6}] (64 grid points per Gram interval, bisection to 10^{-11}). decomp.py splits s·Z′/Z(g_n) into the pair term P = s·Σ_pair 1/(g_n − γ) (Jerby's |1/a_1 + 1/a_2|), the other zeros within six Gram intervals, the drift s·π/4 (Z′/Z = π/4 + Σ 1/(t − γ) + O(1/t)), and the remainder. Signs are taken relative to Z′/Z, so the columns add up to s·μ (logs/decomp_low.txt):

| n | t | s·μ | side | a_1, a_2 | P | other zeros (±6 intervals) | drift + far zeros | nearest zero beyond the near wall (Gram units) |
|---|---|---|---|---|---|---|---|---|
| 2,767,213,045 | 9.74·10^8 | 1.6978 | left | 0.971, 0.865 | 2.186 | −0.523 | +0.035 | 0.001 |
| 2,415,933,307 | 8.56·10^8 | 1.7056 | left | 0.994, 0.942 | 2.068 | −0.386 | +0.024 | 0.093 |
| 20,798,828,979 | 6.61·10^9 | 1.7143 | left | 0.968, 0.914 | 2.127 | −0.445 | +0.032 | 0.058 |
| 823,483,757,101 | 2.22·10^11 | 1.7160 | left | 0.961, 0.867 | 2.194 | −0.506 | +0.028 | 0.001 |
| 2,813,656,271 | 9.89·10^8 | 1.7468 | right | 0.924, 0.987 | 2.096 | −0.370 | +0.021 | 0.087 |
| 1,917,816,785,061 | 5.00·10^11 | 1.7600 | left | 0.987, 0.744 | 2.358 | −0.460 | −0.137 | 0.037 |
| 273,364,058 | 1.10·10^8 | 1.8125 | left | 0.994, 0.919 | 2.094 | −0.125 | −0.157 | 0.129 |
| 1,565,734,083 (zero-free) | 5.68·10^8 | 1.3610 | – | – | 0 | +1.327 | +0.034 | 0.106 / 0.062 |

All 16 lowest confined points (the full set is in logs/decomp_low.txt) share one configuration. One zero of the pair sits at the far wall (a_1 = 0.94–0.999), so P ≥ 2 holds, close to equality. The next zero beyond the empty interval sits at the near wall (0.001–0.13 Gram units outside), and together with the other near zeros it subtracts 0.12–0.59. The drift π/4 and the far zeros nearly cancel: their sum lies in [−0.16, +0.04]. So the lower scale is set by the pair *and* the opposing zero; Jerby's l. 2044–2051 caveat names this omission, and it is quantitatively decisive.

## 3. (iii) μ > 4 at confined pairs

Not violated in any sample. The minimum confined μ is **4.8099**: n = 273,364,058, t = 109,580,510.618, R1, Arb μ ∈ [4.80992529117390 ± 4.9·10⁻¹⁶]. Next come 5.0625 (R1), 5.0845 (R2b) and 5.2953 (windows). A confined μ < 4 needs s·μ < 8π/log(t/2π), which is 1.52 at 10^8, 1.33 at 10^9 and 1.07 at 10^11. Observed confined minima are ≥ 1.698 everywhere, and below 10^8 no isolated bad point at all has μ < 4 (ConjFals exhaustive to n = 2.54·10^8; Jerby to 2.5·10^7). In the CUE model the event behind a confined μ < 4 has frequency ≈ 4·10^{-3} per isolated bad point at t = 10^8 (s·μ < 1.52) and ≈ 8·10^{-4} at 10^11 (s·μ < 1.07). ζ is still 10–100× below CUE there, so none is expected in reachable ranges, but the bound is not structural.

## 4. Dynamical confinement and collisions along the diagonal

gram_win.c -dyn implements Jerby's family Z_N(t; r) = cos θ(t) + r Σ_{k=2}^{M} k^{−1/2} cos(θ(t) − t log k), M = ⌊√(g_n/2π)⌋ (l. 448–467, diagonal path a = r1). For r = 0, 1/NR, …, 1 it counts zeros in the rectangle g_{n−1} < Re t < g_{n+1}, |Im t| < 2.5 by the argument principle, with adaptive arg steps < 0.25 rad and θ(x + iy) = θ(x) + iyθ′(x) + O(y²/x). |Im t| ≥ 2.5 holds no zero, since the core dominates: Σ_{k≥2} k^{−1/2−|y|} ≤ ζ(3) − 1 < ½. It also counts real sign changes on 401 points. Validation: Jerby's g_126 (Fig. 2) is confined, with two real zeros for all r, as his figure shows.

- The 16 lowest confined points: 12 are confined along the diagonal (two zeros in the strip for every r), including the lowest (1.6978, n = 2,767,213,045) and L57/L63 at 2.2·10^11 and 5·10^11. Four (n = 273,364,058, 2,415,933,307, 2,293,437,789, 2,520,246,131) are not: at r = 1 the section's far-wall zero lies just outside the strip. Z's zero there sits 0.001–0.019 Gram units inside the wall, within reach of the Riemann–Siegel remainder by which Z_N(t; 1) differs from Z/2. So the dynamical definition is unstable exactly at the extremes the liminf law is about.
- **All 16 collide.** For n = 2,767,213,045 (μ = 5.096, NR = 100) the pair is real for r ≤ 0.28 and non-real for r ∈ [0.29, 0.98]. There min |Z_N| on the real segment reaches 0.33 and the strip still holds two zeros. It is real again at r ≥ 0.99. All 16 pass through an interval of r with no real zero in the strip. At r = 1, 11 have two real zeros and 4 have one (the other has left through the wall). For n = 20,798,828,979 the pair is still non-real at r = 1: the section misses Z's shallow dip, Z(midpoint) = −0.0031, against a remainder of ≈ 0.003. Logs: dyn_low.log, dyn_high.log, dyn_fine_2767213045.log.
- Controls: of 24 isolated bad points sampled at t ≈ 9.7·10^8 (μ = 9.9–488), none collides and 23 are confined along the diagonal (logs/dyn_typ_sample.log). Of 200 isolated bad points sampled from the author's range, 5 collide and 195 are confined (logs/dyn_gbg_sample_summary.txt).

So the interpretive reading "whenever the first-order variation drives a confined pair toward collision, the accompanying transverse motion is forced to be sufficiently strong to resist that tendency" (l. 597–600) fails numerically for confined Z-points with μ ≈ 5–7 > 4. That reading is not a stated theorem. Caveats: the r-grid is 1/50–1/100, the real grid has 401 points, and the computation is double precision, not interval arithmetic. The contrast "DH collides, Z does not" (l. 803–813) is not a dichotomy, since Z's extremal confined pairs collide too and in 11 of 16 come back as two real zeros by r = 1.

## 5. F3, the trunk-slope law

trunk.py computes ρ = |Z′(g_n)|(g_{n+1} − g_n)/A_act, with spacing π/θ′, on all 6,516 isolated bad points with μ/θ′ < 0.62 (windows and resonance scans; logs/trunk.txt). In the author's range: his five lowest-viscosity points have ρ = 0.91, 0.99, 1.44, 1.13, 1.07, and his 40 lowest span 0.86–4.7 (median 1.10). So there the estimate holds with constant ≈ 1. Over all 24,269 of his isolated bad points ρ = 0.13–154 (median 1.19). Above 10^8 the 40 lowest per half-decade have median 1.6–8.7 and maximum 17–358, and over all points ρ = 0.83–1129. The extremes at large t are humps with *both* neighbours near zero (|Z(g_n)|/A_act = 352 at n = 2,767,213,045), not his coherent one-small-one-large pattern. Jerbyfals's three counterexamples have ρ = 6.1, 4.3, 4.1. The law carries no constants and is declared qualitative (l. 1858–1863), so it is not refutable. As a quantitative law it holds only at the heights where it was observed.

## 6. F5, F6: Davenport–Heilbronn

The section is not defined in the paper (l. 727–744 gives no N and no terms). We use the analogue of the Z case, S_N(t; r) = cos θ_DH(t) + r Σ_{k=2}^{N} c_k k^{−1/2} cos(θ_DH(t) − t log k), with c = (1, κ, −κ, −1, 0) mod 5, θ_DH(t) = Im log Γ(3/4 + it/2) − (t/2) log(π/5) and D(s) = Σ c_k k^{−s} (dh.py at 30 digits, dh2.py in double; logs/dh.log, dh2.log):

- **Gram point**: exact g_44^DH = 85.58474348258131. The Zenodo table gives 85.58485882085794. The table's DH Gram points solve θ_DH − 1/(48t) = nπ (Stirling without the 1/t term): θ_DH(g_2^table)/π = 2.00063, and g_44 is shifted by 1.153·10^{-4} = 1/(48 t θ′).
- **Exact viscosity at g_44**: normalised Z_DH = e^{iθ_DH} D(½ + it) gives Z = −0.411255, Z′ = 0.895576, μ = **2.1777**. The completed ξ_DH(½ + it) = |G| Z_DH gives μ = **2.9601** (d/dt log|G| = −0.78248). The Zenodo table's ZDH(g_n) is proportional to ξ_DH (ratios at n = 43, 44, 45 agree to 4 digits), but its derivative is 7.3% too small: mu_exact = **2.7435**, not 2.9601. The paper's ≈ 1.45 is "for the section used in Figure 3".
- **Section family** (two θ variants agree to 10^{-3}): μ_sec = |S′/S|(g_44; 1) and r_* (first zero of Δ_44(r) = S(g(r); r), g(r) the continued extremum):

| N | 6 | 7 | 8 | 9–10 | 11 | 12 | 13 | 14–15 | 21 | 28 | 42 | 59 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r_* | 0.9706 | 0.8766 | 0.9136 | **0.7131** | 0.6425 | 0.6639 | 0.6323 | 0.5433 | 0.4544 | 0.4670 | 0.4086 | 0.4555 |
| μ_sec | 7.91 | 2.75 | 3.78 | 1.044 | 0.886 | 1.012 | 0.868 | 0.558 | 0.638 | 0.704 | 0.660 | 1.014 |
| μ_sec + d log\|G\| | 8.69 | 3.54 | 4.56 | 1.826 | 1.669 | 1.795 | 1.650 | 1.340 | 1.420 | 1.487 | 1.443 | 1.796 |

N ≤ 5 has no collision. No N gives r_* = 0.72346. The closest is N = 9, 10 (0.71305, 1.4% off), where μ = 1.04 normalised or 1.83 with the Γ drift. μ ≈ 1.45 occurs only with the Γ drift and N ≥ 21, where r_* = 0.40–0.47. **Discrepancy**: for every tested N the two printed DH numbers are mutually inconsistent. Without the author's section they cannot be checked. Whatever the section, the qualitative claim survives: the pair collides at μ_DH < 4.

## 7. Side numerics (Zenodo data)

- 220,373 good and 24,269 isolated bad for n ≤ 250,000 (l. 342–343): **reproduced** from the table's Z values. The README's GBG = 24,270 counts n = 250,000, which cannot be classified without g_250,001.
- Smallest isolated-bad μ "≈ 4.26178" (l. 353): the table and mpmath give **4.26161122355** (n = 8,569). The printed value is off by 1.7·10^{-4}.
- μ_exact = |Z′/Z| holds internally in the ζ table (0 mismatches). In the DH table it does not (§6).

## 8. Certified points (Arb, python-flint 200 bits; certify_conf.py, logs/certify_conf.log)

Each entry is certified isolated bad. μ and s·μ are enclosed. Z at the midpoint of the screened pair has the sign opposite to Z(g_n), which proves ≥ 2 zeros in that Gram interval. "Exactly 2" rests on screening (64 points per interval). mpmath siegelz agrees to all printed digits.

| n | g_n | Z(g_{n−1}), Z(g_n), Z(g_{n+1}) | μ | s·μ |
|---|---|---|---|---|
| 2,767,213,045 | 973,586,544.73709663760849 | 0.01797, 6.51865, 0.01851 | 5.09598420420750 | **1.69784494919686** |
| 2,415,933,307 | 856,157,979.08783297077634 | 0.00183, 8.33911, 2.79367 | 5.08445085093938 | 1.70562712731325 |
| 20,798,828,979 | 6,608,889,900.79222123876440 | 0.01048, 6.63584, 1.31782 | 5.66778785350560 | 1.71426285543296 |
| 823,483,757,101 | 222,170,877,817.55668166534 | 0.01990, 5.88266, 0.02823 | 6.63363570657363 | 1.71602966613692 |
| 1,917,816,785,061 | 499,999,992,481.22706577579 | 0.02081, 6.28618, 0.82952 | 7.03088850356622 | 1.76001517446895 |
| 273,364,058 | 109,580,510.61801599035710 | −0.00290, −10.00992, −4.90534 | **4.80992529117390** | 1.81246976219157 |

## Reviewer

One `reviewer` subagent saw the three central claims with the statements quoted inline. Its verdicts:

- **A, GAP.** The F2 trend is fairly stated as "no finite refutation, unfavourable trend". Finite-N CUE samples do not establish liminf 0, and sign changes give "≥ 2" zeros, not exactly 2.
- **B, GAP.** The collisions are strong numerical evidence, but discrete r- and t-grids are not a proof, and l. 594–600 is an interpretation, not a theorem.
- **C, GAP (narrow).** Say "no tested N", not "no N".

The text above has been worded accordingly. None of the three is reported as a counterexample or as a false law.

## Novelty

- arXiv monthly listings Oct–Nov 2026, 7 math categories (`.scratch/lst.py '\bGram\b|Jerby|Rosser|confined|viscosity' --zenodo '"Gram point"'`): the Gram-titled papers are 2610.17349 (derivatives of the Gram function), 2610.21952 (Jerby), 2610.22107 (second moment along Gram points) and 2610.39856 (moments at "random matrix Gram points", the CUE lattice used here, for moments not viscosity). None tests a liminf of s·μ, a dynamical confinement or collisions in Jerby's family. Zenodo "Gram point": 4 unrelated records.
- Our files: jerbyfals (μ > 4 false; confined min 4.81 and refined conjecture open), ConjFals, R441 (statement harvest). Not done before here: confined s·μ minima per height, the frequency trend, the decomposition, the dynamical counts and the DH reproduction. No random-matrix model of this statistic was found.

RH relevance: N8 (Gram/Rosser counting), N13 (Jerby's mechanism). N8: at Rosser-respecting length-2 blocks, the lower edge of the viscosity is set by one configuration: a pair at the far wall plus a zero at the near wall. With height it drifts toward the CUE tail, with no barrier at 2. N13: Jerby's collision-avoidance mechanism is contradicted numerically at the extremal confined points, where his diagonal deformation collides and in most cases recovers. RH itself is untouched.

## Files

STATEMENTS.md (exact quotes); gram_win.c (jerbyfals's scanner plus -zeros and -dyn modes); analyze_conf.py, decomp.py, trunk.py, jerby_data.py, cue_visc.py, dh.py, dh2.py, certify_conf.py, cert_args.py (pair midpoints for certify_conf.py from zeros_low.log), dyn_sample.py (summary of the 200-point sample); copies of reso.c, certify_cx.py, certify_gram.py, certify_zerofree.py, make_windows.py, validate.py. logs/: run_M.log (+ .time), windows_M.txt, analysis_conf.txt, zeros_low.log, decomp_low.txt, dyn_low.log, dyn_high.log, dyn_fine_2767213045.log, dyn_jerby_examples.log, dyn_typ_sample.log, dyn_gbg_sample.log, dyn_gbg_sample_summary.txt, trunk.txt, jerby_data.txt, cue_64.log, cue_128.log, cue_64_test.log, dh.log, dh_out.json, dh2.log, certify_conf.log (*.time: max RSS ≤ 76 MB, except jerby_data.py at 180 MB).

## Suggested TREE line (N8)

- R### node N8 (Jerby 2610.21952 finite laws): μ > 4 already false. The confined-pair liminf s·μ = 2 is contradicted by its trend: window minima 1.70–2.06 at 10^8–5·10^11 (Arb 1.6978, n = 2,767,213,045), frequency below 1.948 up ×100, CUE-model tail to 0.29. The lower scale omits the zero at the near wall (−0.12 to −0.59). The trunk slope holds in the author's range only (ρ ≈ 1 below 1.7·10^5; ρ up to 1129 above 10^8). The DH numerics r_* = 0.72346 and μ = 1.45 are not reproducible (no tested N ≤ 59 gives both). Extremal confined pairs collide along Jerby's diagonal (numerically) [JerbyConf](Riemann/Deep/JerbyConf/REPORT.md)

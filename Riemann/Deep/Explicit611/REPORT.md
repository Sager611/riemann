RECORD: |ζ(1/2+it)| ≤ 0.568 t^{1/6} log t for all t ≥ 3, and |ζ(1/2+it)| ≤ 0.4821 t^{1/6} log t for all t ≥ 8.97·10^17 (arXiv 2602.05614 Thm 1: 0.611 and 0.566), by the same method and lemmas as 2602.05614, with its wasteful steps removed; ball-arithmetic certificate.

# Explicit611 — lowering the explicit constant in |ζ(1/2+it)| ≤ C t^{1/6} log t

## Theorem (new)
For real t:
- (a) t ≥ 3: |ζ(1/2+it)| ≤ 0.568 t^{1/6} log t. (Before: 0.611, Revers 2602.05614 Thm 1(a); 0.618 Hiary–Patel–Yang 2207.02366.)
- (b) t ≥ 8.97·10^17: |ζ(1/2+it)| ≤ 0.40913 t^{1/6} log t + 3.01756 t^{1/6} − 2.920 ≤ 0.4821 t^{1/6} log t. (Before: 0.566, Thm 1(b), from 0.469028 t^{1/6} log t + 3.99399 t^{1/6} − 21.4623.)
- (c) t ≥ 10^10: |ζ(1/2+it)| ≤ 0.45002 t^{1/6} log t + 2.58840 t^{1/6} − 2.905 ≤ 0.5598 t^{1/6} log t.
- Leading coefficient of the method: 0.40913 (r1 = 80) and → min_c √(1/c + (32/(15√π))√c)/(3π^{1/6}) = 0.40278 (c = 1.4029) as t → ∞ (paper: 0.46903).
- Literal-hypotheses variant (van der Corput only where the paper's Lemma 2.7 states t ≥ 5·10^6): (a) holds with 0.5782 (`certify.py 0.5782 ABCEF 892 39894 1e10 6000`, log `run_literal.log`).

Dependencies (same as the paper): Lemmas 2.3–2.9 of 2602.05614 (2.3, 2.4 are Hiary–Patel–Yang's Lemmas 2.3, 2.4, refereed; 2.5–2.9 unrefereed), Gabcke's Riemann–Siegel remainder in the paper's form (`MainSumEstimate`, TeX l. 750-760), and arb's ζ.

## Where 0.611 binds (feasibility, done first)
- 0.611 binds at the junction of the trivial Riemann–Siegel bound and the van der Corput bound, t ≈ 3.3–3.8·10^7: 4t^{1/4}/(2π)^{1/4} − 2.0895 = 0.61087 t^{1/6} log t at 3.3·10^7 (TeX l. 194-205); region 1 (TeX l. 963-976) gives 195.52 vs 195.53 at 3.8·10^7. Every one of the 10 regions (l. 958-1081) also hits exactly 0.61100 at its right end t1, and the last region (one φ = 0.33711 for [3.92·10^8, 8.97·10^17]) hits 0.611 at 8.97·10^17. None of this is a directly computed range. The paper's l. 96-101 shows that no constant below 0.5075 is possible for t ≥ 3 (the ratio exceeds 0.5075 near t = 45.62).
- How 0.611 and 0.566 arise: [3, 200] is Hiary's interval-arithmetic computation, quoted (l. 180-192). [200, 3.8·10^7] uses the trivial Riemann–Siegel bound (l. 194-270). [3.8·10^7, 8.97·10^17] uses 10 regions, each with one hand-tuned (φ, c) and uniform constants M1–M7. Above 8.97·10^17 the paper uses φ = 1/3, c = 1.34 and the T4/T5 split, with ratio 0.56512 at 8.97·10^17. The free parameters are φ, c, the region ends and r0 = 4. No code is published (the TeX has no links; the source holds only ZetaRevers.tex).
- At t = 3.8·10^7: K = 350, n1 = 2459, R = 7. The last block [KR, n1] has 10 terms but is bounded as a full block r = 7 (≈ 9.3 of the total 195.5). The trivial part 2Σ_{n<4K} n^{-1/2} ≈ 147 dominates.
- Hence the two refinements of round R366 are not marginal. At t = 3.8·10^7 with the paper's own K = 350, c = 1.82 and μ0 (`consistency.log`), the per-block formula gives ratio 0.59615, against 0.61096 for the paper's uniform bound. A short last block alone (T3, l. 2045-2049) lowers it to 0.56875; adding per-block μ_r gives 0.56719. The T5 refinement (l. 1230-1234), together with removing the √ split it sits in, is worth ≈ 0.04 at 8.97·10^17.

## Changed inequalities (line numbers: TeX source of 2602.05614v1, `ZetaRevers.tex`)
1. l. 589, 597-608 (`Intermediate3`). Old: μ ≤ μ0 = ½(1+1/r0)(2+1/r0+1/(t0^φ−1)) = 1.40625 for every block. New: each block r has its own μ_r = ½(1+1/r)(2+1/r+1/(K−1)) (Lemma 2.5 with K0 = K, which the lemma allows), and μ_r → 1. Lemmas 2.7 and 2.9 are statements about one block, so per-block constants are allowed.
2. l. 610-747 (M1–M4). Old: uniform bounds through t0, t1, R0, R1 and W^{-1/3} ≤ 1/(5π^{1/3}t0^{φ−1/3}). New: for t ≤ 10^10, each block is evaluated with its own r and K over the exact interval t ∈ [2πn1², 2π(n1+1)²) in ball arithmetic. For the tail: W^{-1/3} ≤ 1/(π^{1/3}(r+1)) and K/t^{1/3} ∈ [1, 1+t0^{-1/3}] per r.
3. l. 798-816 (T3, the paper's own remark l. 2045-2049). Old: max_{1≤L≤n1−KR+1} is replaced by max_{1≤L≤K}, so the last block is bounded as K terms. New: with L_R = n1 − KR + 1, 2|Σ_{n=KR}^{n1}| ≤ min(2L_R/√(KR), (2/√(KR))·√((L_R+M−1)(L_R/M + Y))). Here Y is the paper's bound (l. 486-497) for (2/M)Σ(1−m/M)|s'_m|. Lemma 2.4 applied with L = L_R is increasing in L, and Lemma 2.7 is uniform in L ≤ K.
4. l. 570 and l. 958-1081 (10 regions, one (φ, c) each, K = ⌈t^φ⌉). New: every n1 (39,394 intervals) gets its own integer K ∈ [⌈t^{1/3}⌉, ⌈t^{0.35}⌉] (= ⌈t^φ⌉ for some φ ∈ [1/3, 0.35], checked as K³ > t_hi and K−1 < t_lo^{0.35}). K is chosen so that the last block is short, and each block gets its own c_r ∈ (0, 3). The paper's concluding remark (l. 2035-2040) calls hundreds of regions "a tremendous amount of additional computational work"; the whole check here runs in about 6 minutes.
5. l. 434-453 and 486-520 (Lemma 2.8 (a),(c),(d) and M ≤ cW^{1/3}+1, M^{-1/2} ≤ (cW^{1/3})^{-1/2}). New: for n1 ≤ 6000 (t < 2.3·10^8) and for the tail blocks r < r1, the sums are computed exactly: Σ_{m≤min(M,L−1)} (1−m/M)·X_m, with X_m the Lemma 2.7 bound. M is any integer with M − 1 < 3W^{1/3} (equivalently M = ⌈cW^{1/3}⌉ for some c ∈ (0,3), so M ≤ K−1 by Lemma 2.7(b)), and s'_m(L) = 0 for m ≥ L. Each block takes the smaller of this bound and the relaxed one.
6. l. 856-880 (`SumIntegralEstimate`) is no longer needed for t ≤ 10^10, where the block sum is exact for each n1. In the tail, the per-r terms are exact for r < r1, followed by a decreasing-summand integral (item 7).
7. l. 1133-1148 (`T2T3T4T5`: √(a+b) ≤ √a + √b) and l. 1163-1188 (T5: 1+1/r ≤ 5/4, Σ_r 1 = R−3). New: no split. For r ≥ r1, s(r) = √(A/(r(r+1)) + B(r+1)/(r t^{1/3})) is decreasing, so Σ_{r1}^{R} s ≤ s(r1) + ∫_{r1}^{R} √(a²/x² + b²) dx = s(r1) + [√(a²+b²x²) − a·asinh(a/(bx))]_{r1}^{R}, using R ≤ t^{1/6}/√(2π) and asinh y ≤ log 2y + 1/(4y²). Gain at 8.97·10^17: the paper's remark l. 1230-1234 (sharper Σ√(1+1/r)) alone gives ≈ 0.006; removing the split and using per-block μ_r gives the rest of 0.566 → 0.482.
8. l. 1864-1885 (Lemma 2.7 proof, part (a)). The hypothesis t ≥ 5·10^6 is used only to show R ≥ 4. Parts (b)-(e) need only t ≥ 2012, r ≥ 4 and K ≥ t^{1/3}, and Lemma 2.9 adds nothing else. We use the block bounds for 1.58·10^6 ≤ t < 5·10^6 only where R ≥ 4 (asserted). Without this step, the binding point is the trivial Riemann–Siegel bound at n1 = 892 (ratio 0.57813), giving 0.5782.
9. l. 180-192 (Lemma 2.1, [3, 200] quoted from Hiary). New: our own ball evaluation of ζ(1/2+it) on [3, 300] (5,972 balls). For 300 ≤ t < 1.58·10^6 we use the exact sum 2Σ_{n≤n1} n^{-1/2} + 1.463 t^{-1/4} + 0.127 t^{-3/4} (the first line of the proof of Lemma 2.2, l. 216-222), with no integral comparison.

## Where 0.568 binds
t ≈ 2.30·10^6 (n1 = 605, K = 133, R = 4: the trivial part 2Σ_{n<532} n^{-1/2} plus a single 74-term block r = 4): ratio 0.56772. There the van der Corput bound and the trivial Riemann–Siegel bound cross. Next binding points: 0.56329 (trivial bound, n1 = 7, t ≈ 308) and 0.5598 for t ≥ 10^10. The remaining gap to the floor 0.5075 belongs to the method: the trivial part 2Σ_{n<4K} n^{-1/2} with K ≥ t^{1/3} and r ≥ 4 is forced by the lemma hypotheses.

## Certificate (one command, 6.5 min, single thread, max RSS 49 MB, uv + python-flint + numpy)
```
uv run --no-project --with python-flint --with numpy python3 Riemann/Deep/Explicit611/certify.py 0.568 ABCEF 500 39894 1e10 6000
```
Expected last line: `ALL PARTS OK` (log: `run_final.log`). Parts:
- A: ball evaluation of ζ on [3, 300].
- B: exact-sum trivial Riemann–Siegel bound for n1 ∈ [6, 500] (t < 1.577·10^6).
- C: one interval per n1 ∈ [501, 39894] (t < 1.00004·10^10), with T1 = 2Σ_{n<4K} n^{-1/2} exact plus blocks r = 4..R, the last one partial.
- E: uniform tail bound for all t ≥ 10^10.
- F: the same bound for t ≥ 8.97·10^17 (part (b)).

Rigor:
- python-flint 0.9.0 arb balls at 96 bits; every final comparison is arb "certainly <", with the bound's upper end against the lower end of C·t^{1/6} log t at the left end of the interval (the right side increases in t).
- Decimal constants (0.568, 1.463, 0.127, 0.35) are balls containing them.
- K, c_r, M_r and r1 are chosen in floating point, and every inequality is then rechecked in balls. The hypotheses are asserted on every interval: K³ > t_hi, K−1 < t_lo^{0.35}, R ≥ 4, t ≥ 2012, 0 < c < 3, M − 1 < 3W^{1/3}, M ≤ K−1.
- Partial sums Σ n^{-1/2} are exact ball sums, or Hurwitz ζ(1/2, N+1) above 2·10^5.
- Tail monotonicity is checked: b2 > 0 and b2 > |b3|(log t0/6 + 1)/t0^{1/6} make the ratio decreasing for t ≥ t0.

Sanity and cross-checks:
- `sanity_authors.py` reproduces all 40 coefficients a1..a4 of the paper's 10 regions and b1..b3 at t0 = 8.97·10^17 and 7·10^11 to every printed digit. Each region's bound equals 0.61100 at its t1 (log `sanity.log`).
- `crosscheck_tail.py` (float) compares direct per-block sums at 11 values of t ∈ [10^10, 10^20] (K = ⌈t^{1/3}⌉, every block the better of the relaxed and exact-Weyl bounds) with the closed-form tail bounds of Parts E/F. The direct sums are always below: max direct/formula 0.9942; at 8.97·10^17 the ratios are 0.47925 and 0.48205 (log `crosscheck_tail.log`).
- `consistency.py` (float): at t = 3.8·10^7 the per-block formula with the paper's K, c, μ0 and a full last block gives 190.78 ≤ 195.52, the paper's region-1 bound. The per-block formula is therefore the paper's bound before its uniform relaxations.
- First run, without exact Weyl sums: 0.5725 for t ≥ 3 and 0.4841 for t ≥ 8.97·10^17 (log `run_05725.log`).

## Novelty pass (2026-11-05/06)
- `python3 .scratch/lst.py 'zeta.*(explicit|estimate)|explicit.*zeta|subconvex|van der Corput|exponent pair' math.NT math.CA --zenodo '"zeta" "explicit" "critical line"'`: Oct–Nov 2026 listings give only 2611.05267 and 2610.07765 (van der Corput property of cubes and squares, combinatorial). Zenodo (most recent first): nothing relevant.
- HAL (`zeta explicit`, 20 newest): nothing on ζ(1/2+it) bounds.
- arXiv API:
  - `all:Revers AND all:zeta`: only 2602.05614, still v1.
  - `all:"sub-Weyl" AND all:zeta`: only Patel–Yang 2302.13444, still v1.
  - `ti:explicit AND ti:zeta` and `abs:"van der Corput" AND abs:zeta AND abs:explicit`: newest relevant is Dhiman–Kadiri–Quesada-Herrera 2609.00537 (explicit B-process and approximate functional equations), which has no t^{1/6} log t constant.
- web_search (arxiv): nothing better than 0.611.
- Repo: no earlier NOVELTY or DISCOVERIES entry on this constant (TREE N9, R343, R366 only track it).

## Caveats
- The certificate is only as sound as 2602.05614's Lemmas 2.5–2.9 (unrefereed, v1). An error there would also void 0.611. The new steps use those lemmas per block, exactly as stated. The only exception is item 8, which has the literal variant 0.5782.
- Part A trusts arb's ζ (standard, rigorous ball arithmetic).
- Needs the reviewer verdict before a DISCOVERIES entry (AGENTS.md).

RH relevance: N9 (explicit subconvexity on the critical line; quantitative record)

## Files
- `certify.py`: the certificate.
- `sanity_authors.py`: re-derivation of the paper's constants.
- `crosscheck_tail.py`: float cross-check of the tail formula.
- `consistency.py`: float comparison with the paper at t = 3.8·10^7.
- `novelty.log`: listing, Zenodo and HAL scan output.
- `run_final.log`, `run_0568_v2.log`, `run_05725.log`, `sanity.log`, `crosscheck_tail.log`, `consistency.log`: run logs. `run_final.log` is the certificate run with the final code.

## Suggested DISCOVERIES line (after reviewer verdict)
|ζ(1/2+it)| ≤ 0.568 t^{1/6} log t for t ≥ 3 and ≤ 0.4821 t^{1/6} log t for t ≥ 8.97·10^17 (before 0.611, 0.566); ball-arithmetic certificate.

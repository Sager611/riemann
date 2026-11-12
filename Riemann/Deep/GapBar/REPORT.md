NO BARRIER: GapCeil's 3.30279 is a real-amplifier artefact; BM's mollifier moved to the midpoint, M(s+iκπ/2L) (equivalently complex P), reaches 3.32658 (r = 2, float); certified Λ > 3.324; failure of every amplifier certified only at grid points (r ≤ 3 at κ = 3.335, 3.35, 3.5, 4, 5; r = 2 also 3.34-3.45; r = 1 also 3.30)

# gapbar: the ceiling of Bui–Milinovich's criterion was set too low (node N8, falsification)

## Answer
- **The barrier asked for is false.** GapCeil's ceiling 3.30279 (r = 2) holds only for real amplifiers in BM's asymmetric
  set-up F = e^{ivtL}ζ(s)ζ(s+iκπ/L)M(s). BM's method allows two more things:
  - BM's own Lemma main covers mollifier shifts, so the mollifier can sit at the midpoint, M(s + iκπ/(2L)), with real P;
  - BBLR's theorem is for complex coefficients, so P may be complex.

  Both lift the r = 2 ceiling to **3.32658** (float). No barrier can hold at κ ≤ 3.324, where success is certified.
- **Mechanism (exact identity).** The midpoint shift multiplies BM's kernel by e^{ic(r2−r1)/2} and makes it real and of
  product form: **U_L·conj(U_R)·e^{ic(r2−r1)/2} = 4 sin(c(X1−X2+σ)/2)·sin(c(X3−X4+σ)/2)** (c = κπ, σ = τ1 − τ2,
  r1 = X1+X2, r2 = X3+X4).
  - Consequence: the sup over all complex amplifiers equals the sup over midpoint real amplifiers.
  - GapCeil kept Re of the unshifted kernel, which is exact for real P. Extended to complex or shifted amplifiers, it
    violates an exact symmetry of the true mean values; the Hermitian kernel satisfies it (§(ii)).
- **New record (Arb): Λ > 3.324.**
  - Witness: the midpoint mollifier with a real degree-80 P, r = 2, ϑ = 2499/10000, v = 1.367258.
  - Certificate: h1 = 1.00056020750 ± 4.5·10⁻¹⁵. It uses BM's Lemma main as printed, with mollifier shifts α3 = −β3 = iκπ/(2L).
  - Also certified: Λ > 3.323 (midpoint, degree 64), and with complex P Λ > 3.322 (degree 64) and Λ > 3.315 (degree 32).
  - Old record: 3.30 (GapCeil). The new value is 0.0026 below the float ceiling.
- **Corrected barrier points (Arb).** At each point, every complex amplifier (hence every real one, midpoint or not)
  fails for every v and every ϑ ≤ 1/4:
  - r = 1 at κ ∈ {3.30, 3.335, 3.35, 3.5, 4, 5};
  - r = 2 at κ ∈ {3.335, 3.34, 3.35, 3.36, 3.38, 3.42, 3.45, 3.5, 4, 5};
  - r = 3 at κ ∈ {3.335, 3.35, 3.5, 4, 5};
  - r = 4, 5 at κ = 3.335: not certified. The high-mode bound broke at N = 40, and larger N does not fit under the
    300 MB cap.

  So at κ ∈ {3.335, 3.35, 3.5, 4, 5} the criterion fails for every amplifier with r ≤ 3, and the r = 2 method lies
  between 3.324 (certified success) and 3.335 (certified failure).
- **Still open for a BARRIER:**
  - κ between grid points: a naive Arb ball in κ still blows up, with entry radii of 10⁵⁷-10⁶⁴ on [3.4, 3.45];
  - κ beyond the grid (5 for r = 1, 2, 3);
  - r ≥ 4 (float only: ceilings 3.147 and 3.074 for r = 4, 5), and any bound uniform in r.

## (i) Why the barrier fails: two admissible extensions (DERIVATION.md §1, §3)
1. **Mollifier shifts (real P, as printed).**
   - BM's Lemma main (l.370-384) gives I_{α,β}(M) with M(½+α3+it)M(½+β3−it), uniformly for shifts ≪ 1/L.
   - In c(α, β) the mollifier shifts enter only through y^{−(α3+β3)x−α3(x1+x2)−β3(x3+x4)}.
   - With α3 = −β3 = iκπ/(2L) this factor is e^{ic(r2−r1)/2}.
   - F(t) = e^{ivtL}ζ(s)ζ(s+iκπ/L)M(s+iκπ/(2L)) still vanishes at t_n and t_n − κπ/L, which is all the reduction
     (l.164-202) needs. F' is handled by Q(L⁻¹(∂_{α1}+∂_{α2}+∂_{α3})) as in l.402.
2. **Complex P.**
   - BBLR's theorem as quoted (l.500-517) is sesquilinear in complex a_h. BM's evaluation is linear in b_ib_j
     (l.531, l.537), so complex P gives b_i·conj(b_j) with the same kernel.
   - Wirtinger's inequality holds for complex-valued f (l.140-155).
   - The kernel is Hermitian, not real: Γ = Γ^† with a large imaginary part (r = 1, κ = 3.3: a[0,1] = −2.156 − 10.061i).
3. **Equivalence.** With f = e^{icU/2}g (U = log h/L, so e^{icU/2} = h^{iκπ/(2L)}) the identity above turns the complex
   form into the real form of the product kernel. Hence complex sup = midpoint real sup.
   - Float check: r = 2, κ = 3.335, 14 modes: 0.9932086 = 0.9932086.
   - Arb check: the imaginary parts of the twisted arrays contain 0 in every run.

## (ii) Consistency check: an exact symmetry the unshifted real kernel violates (DERIVATION.md §4)
- For real P, conj G_β(−t−α) = G_{α−β}(t) with G_β = ζ(s)ζ(s+iα)M(s+iβ). So the true main terms satisfy
  c_j(v; β) = c_j(v; α−β).
- Float test, r = 2, κ = 3.31, degree-10 g, amplifier g·e^{iμU} (`mutest.py`):

| μ/c | 0 | 1 | 0.2 | 0.8 | 0.5 |
|---|---|---|---|---|---|
| Hermitian kernel K | 0.9610710 | 0.9610661 | 0.9698746 | 0.9698755 | **0.9778261** |
| Re K (GapCeil) | 0.9610710 | 0.9380853 | 0.9589095 | 0.9425451 | 0.9508246 |

- K is symmetric about μ = c/2 and peaks there. Re K is symmetric about μ = 0 instead, so it cannot be the complex
  extension of the true main term. For real P (μ = 0) the two kernels agree, so GapCeil's real-P numbers stand.

## (iii) Certificates
**Records** (`certmid.py` = GapCeil certr.py on the twisted arrays; `certc.py` = certr.py with the Hermitian form):

| amplifier | r | deg P | ϑ | κ | v | h1 (Arb) | time, peak |
|---|---|---|---|---|---|---|---|
| midpoint, real P | 2 | 80 | 2499/10000 | 831/250 = **3.324** | 1.367258 | **1.00056020750 ± 4.5·10⁻¹⁵** | 59 s, 216 MB |
| midpoint, real P | 2 | 64 | 2499/10000 | 3323/1000 | 1.366889 | 1.00071833023 ± 1.1·10⁻¹² | 59 s, 196 MB |
| complex P | 2 | 64 | 2499/10000 | 1661/500 | 1.366893 | 1.00153083832 ± 4.8·10⁻¹² | 31 s, 168 MB |
| complex P | 2 | 32 | 2499/10000 | 663/200 | 1.364289 | 1.00261659841 ± 3.7·10⁻¹² | 24 s, 172 MB |

- Exact coefficients: `logs/certmid_r2_d80_k3.324.txt` and `_d64_k3.323.txt` (A_n, shifted Legendre basis), and
  `logs/certc_*.txt` (A_n + iB_n).
- Peak memory is maximum resident set size from /usr/bin/time; the scripts' own "peak" lines print MiB.
- For complex P the orientation (which side is conjugated) does not matter: the other orientation is conj P.
- Checks at degree 32: Re P alone gives h1 = 0.9595, conj P gives 0.9277.

**Barrier grid** (`gtw.py` + `runcerttw.py` = GapCeil gcert_r.py/runcert_r.py on the twisted kernel; DERIVATION.md §7):
- Inequality certified: κ²M̃1(v) − M̃0 ⪰ 0 on L²(X^{r²−1}dX) for all v (M̃0 ⪰ 0, MS̃ − 0.8M̃0 ⪰ 0, v_hi M̃0 − MS̃ ⪰ 0
  and v-cells).
- A real symmetric kernel gives the same inequality for complex functions, so it covers every complex f = e^{icU/2}g.

| r | κ | float sup_v λ/κ² (complex) | N | D | prec | operators (v-cells) | max C/δ_D | certified | peak MiB | s |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 3.30 | 0.9978166 | 72 | 50 | 289 | 53 (25) | 0.15 | yes | 192 | 153 |
| 1 | 3.335 | 0.9679570 | 72 | 50 | 289 | 23 (10) | 0.15 | yes | 118 | 115 |
| 1 | 3.35 | 0.9555881 | 72 | 50 | 292 | 15 (6) | 0.14 | yes | 136 | 110 |
| 1 | 3.5 | 0.8448332 | 72 | 50 | 292 | 13 (5) | 0.14 | yes | 170 | 111 |
| 1 | 4 | 0.6026901 | 72 | 52 | 308 | 9 (3) | 0.13 | yes | 200 | 95 |
| 1 | 5 | 0.3728260 | 72 | 55 | 332 | 9 (3) | 0.11 | yes | 221 | 130 |
| 2 | 3.335 | 0.9932086 | 72 | 50 | 289 | 35 (16) | 0.29 | yes (cells shrink 0.4) | 170 | 147 |
| 2 | 3.34 | 0.9892151 | 72 | 50 | 289 | 25 (11) | 0.30 | yes | 169 | 139 |
| 2 | 3.35 | 0.9813094 | 72 | 50 | 292 | 25 (11) | 0.29 | yes | 157 | 148 |
| 2 | 3.36 | 0.9735110 | 72 | 50 | 292 | 19 (8) | 0.29 | yes | 138 | 163 |
| 2 | 3.38 | 0.9582295 | 72 | 50 | 292 | 15 (6) | 0.29 | yes | 170 | 130 |
| 2 | 3.42 | 0.9288859 | 72 | 50 | 292 | 15 (6) | 0.29 | yes | 114 | 143 |
| 2 | 3.45 | 0.9078990 | 72 | 50 | 292 | 13 (5) | 0.28 | yes | 178 | 105 |
| 2 | 3.5 | 0.8747533 | 72 | 50 | 292 | 13 (5) | 0.27 | yes | 165 | 147 |
| 2 | 4 | 0.6383383 | 48 | 76 | 377 | 9 (3) | 0.40 | yes | 245 | 144 |
| 2 | 5 | 0.4019926 | 40 | 76 | 397 | 9 (3) | 0.38 | yes | 209 | 122 |
| 3 | 3.335 | 0.9267478 | 48 | 66 | 330 | 15 (6) | 0.74 | yes | 210 | 92 |
| 3 | 3.35 | 0.9162626 | 48 | 66 | 334 | 21 (9) | 0.72 | yes (cells shrink 0.4) | 178 | 91 |
| 3 | 3.5 | 0.8225136 | 48 | 66 | 334 | 13 (5) | 0.70 | yes | 166 | 105 |
| 3 | 4 | 0.6127076 | 48 | 76 | 377 | 11 (4) | 0.59 | yes | 162 | 139 |
| 3 | 5 | 0.3876087 | 40 | 84 | 422 | 9 (3) | 0.57 | yes | 208 | 136 |
| 4 | 3.335 | 0.8644294 | 40 | 80 | 366 | 13 (5) | 1.28 | **no** (δ_W < 0 in 4 operators) | 227 | 134 |
| 5 | 3.335 | 0.8174263 | 40 | 88 | 387 | 13 (5) | 1.77 | **no** (7 operators); **316 MB, above the cap** | 302 | 232 |

Columns: N φ-modes, Taylor depth D, Arb precision in bits; logs in `logs/certTW_r<r>_<κ>.txt|json`, table from `summarize.py`
(`logs/summaryTW.md`). The shrink-0.7 runs at r = 2, κ = 3.335 and r = 3, κ = 3.35 each failed one cell (Cholesky).
They are kept as `*_shrink0.7_fail.*`; narrower cells passed.

- A direct Hermitian certificate (`gcertc.py`, real 2N embedding) failed at r = 2, κ = 3.335
  (`logs/certC_r2_3.335.txt`): the Im Γ_01 terms are linear in t_1, and C exceeded δ_D in 10 of 16 cell operators.
  The twist makes them vanish.

## (iv) Float landscape (`cland.py`, `ceil45.py`, `cscan.py`, `copt.py`, `optmid.py`)
- **Closure ceilings, complex = midpoint (real, GapCeil):**
  - r = 1: 3.297505 (3.272097);
  - **r = 2: 3.326575** (3.302788);
  - r = 3: 3.238189 (3.226093);
  - r = 4: 3.147354 (3.141892);
  - r = 5: 3.073662 (3.071102).
- **sup_v λ/κ²**, complex (real):

| κ | 3.30 | 3.335 | 3.35 | 3.4 | 3.5 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| r = 1 | 0.99782 (0.97614) | 0.96796 (0.94747) | 0.95559 (0.93560) | 0.91613 (0.89774) | 0.84483 (0.82937) | 0.60269 (0.59664) | 0.37283 (0.37093) | 0.23387 (0.23348) | 0.12113 (0.12109) |
| r = 2 | 1.02194 (1.00229) | 0.99321 (0.97419) | 0.98131 (0.96257) | 0.94336 (0.92550) | 0.87475 (0.85857) | 0.63834 (0.62848) | 0.40199 (0.39567) | 0.26102 (0.25713) | 0.13906 (0.13769) |

  - Monotone decreasing in κ.
  - The complex gain shrinks with κ: 0.019 at 3.335, 0.0014 at 8 (r = 2).
- **κ_max by degree, r = 2, ϑ = 0.2499:**
  - complex P: 3.216396 (4), 3.268823 (8), 3.291084 (12), 3.302541 (16), 3.309199 (20), 3.313404 (24), 3.318217 (32),
    3.322267 (48), 3.323885 (64);
  - midpoint real P: the same to 6 digits at 32, 48 and 64, and 3.324691 at degree 80. The two classes differ (degree 8,
    κ = 3.25: 1.01560266 against 1.01558240) but approach the same closure optimum.
  - Degree 16 already beats GapCeil's real degree-60 record (3.300010).

## (v) The three open pieces
1. **κ between grid points.**
   - Treating κ as an Arb ball [3.4, 3.45] in the twisted certificate (`runcerttwball.py`) fails: entry radii reach 10⁵⁷-10⁶⁴
     through the cancelling Taylor convolutions (`logs/certTWB_r2_3.4_3.45_FAIL.txt`).
   - The high-mode bound itself survives (C ≤ 1 ± 0.5).
   - What is needed: a Taylor model in κ (exact ∂_κ arrays at κ0 plus a majorant remainder), or an analytic monotonicity.
     Neither was done. M1 ⪰ 0 gives κ²M1(κ) − M0(κ) ⪰ κ0²M1(κ) − M0(κ) for κ ≥ κ0, so only the κ-dependence of the operators matters.
2. **Large κ.** No analytic tail. The float ratio falls like ≈ κ^{−2} (0.139 at κ = 8 for r = 2), and the complex/real
   difference vanishes.
3. **r ≥ 4 / uniform in r.**
   - Float ceilings fall after r = 2 (3.238, 3.147, 3.074 for r = 3-5).
   - Certified at grid points for r ≤ 3 (common grid {3.335, 3.35, 3.5, 4, 5}).
   - The r = 4 and r = 5 certificates at κ = 3.335 failed. With N = 40 φ-modes the high-mode constant C exceeds δ_D
     in 4 (r = 4) and 7 (r = 5) operators: C/δ_D reaches 1.3, because the r-weights inflate Γ.
     N ≈ 60 would need about 280 MB or more.
   - The r = 5 run peaked at 316 MB (302 MiB), above the 300 MB cap: the 2 s watchdog missed the spike.
     No further heavy r ≥ 4 runs were made.
   - No bound uniform in r.

## Novelty
- `python3 .scratch/lst.py 'large gap|gaps between (consecutive )?zeros|Wirtinger|mollif|zero spacing|consecutive zeros|twisted fourth moment' --zenodo '"large gaps" zeta'`
  (2026-10/11, 7 math categories; Zenodo) finds the same hits as GapCeil:
  - 2611.14486 and 2610.26038 (small/half-sized gaps);
  - 2610.22927 and 2610.31744 (moderate gaps, biases);
  - 2610.28527 (Dirichlet L, DPC);
  - unrelated Zenodo records.
- **TeX checked** (`.scratch/gapbar/lit/`):
  - BM 1410.3635: M(s) unshifted, P real in all examples (l.333, l.341). No ceiling is stated.
  - Bredberg 1101.3197: M(½+it) unshifted (l.61-71), real coefficients.
  - BHSJ 2511.14415 (l.492-500) proposes a different change, Z(t)Z(t+2πκ/log T)·A(½+it)·A(½+it+2πiκ/log T): two
    amplifiers of length T^{ϑ}, ϑ < 1/8. It has no midpoint mollifier and no numbers.
  - IKT 2510.14309 cites BM's 3.18 only (l.423, l.494).
- GXL evidence API (1 query on complex/midpoint amplifiers in the Wirtinger method): 15 neutral passages
  (BM, BMN 2010, Bredberg, BHSJ, Feng–Wu 1003.0752), none supporting or refuting.
- No midpoint/complex amplifier, and no 3.32-3.33 ceiling or record, was found.
- **New:** the twist identity, the corrected ceilings, Λ > 3.324, and the complex/midpoint barrier certificates.
- **Known:** BM's criterion and Lemma main (with shifts). The idea of shifting factors in Hall-type arguments is classical
  (Hall's Z(t)Z(t+2πκ/log T)), but not for the mollifier.

## Review
- Reviewer (decisive argument pasted inline, async run f06df76c): overall **SIGNIFICANT**.
  - (i) complex-P admissibility: **CORRECT**. Note: conj P keeps the value only when the conjugation convention is
    reversed too. That is what the orientation argument uses.
  - (ii) twist reduction: **GAP**, in claim (c) only. The twist does not keep polynomials, so "the optimum is the midpoint
    mollifier" held only at the L²/sup level. (a) and (b) were judged correct.
    - Repair: the midpoint class with real polynomial P is admissible directly via BM's Lemma main with α3 = −β3 = iκπ/(2L).
      Its kernel is exactly the twisted one (§(i).1), and the record Λ > 3.324 is certified in that class.
      Finite-degree optimisers of the two classes differ (§(iv)).
  - (iii) conclusions: **CORRECT**, taking the certificates as premises.
    - Float crossings (3.30279, 3.32658) are estimates.
    - Failure at one κ does not exclude success at a larger κ: that is open piece 1.
- One correction to the reviewer brief: it quoted "real pencil of Re K = 0.9813" at r = 2, κ = 3.335. The computed value
  is 0.9741944 (0.9813094 is the complex value at κ = 3.35). This does not affect the verdicts.
- Both reviewers saw the degree-64 record (3.323). The degree-80 certificate (3.324) uses the same code.
- **Reviewer 2** (the repair, pasted inline with BM l.360-420, async run bcf576f8): **CORRECT, ROUTINE.**
  - (a) Lemma main allows all six shifts uniformly; α3, β3 enter only through the y-power.
  - (b) The sign of the factor is as stated, and the product form follows.
  - (c) The derivative step keeps (v − S)², because shift differences cancel under Σ∂_α and Σ∂_β.
  - Reviewer 2 calls it "a corollary of the printed lemma and reduction"; the improvement itself is a record.
  - Neither reviewer checked the polynomial or the Arb run.

## Caveats
- **Inputs.** BM's Lemma main and Theorem smoothed as printed, including the mollifier shifts; BBLR's twisted fourth
  moment (complex a_h); M1(v) ⪰ 0 (main term of a mean square), used only in the remark of (v).1.
- **Grid only.** Barrier statements hold at the listed κ, for the listed r.
- **Machinery reused unchanged.** The certificate is GapCeil's: weighted Jacobi basis, high-mode bound, Schur step and
  v-covering, reviewed there. Only the Γ arrays changed.
  - For the new twisted arrays: all phases satisfy |α| ≤ c, so the tail majorant still holds.
  - r = 3 runs used N = 48, D = 66; the r = 3 Jacobi basis was not re-verified with verify_q.py.

RH relevance: N8 (zero spacing, Λ = ∞ direction).
- This corrects a method ceiling: the BM/Hall Wirtinger route stops near 3.327 (r = 2, ϑ < 1/4), not 3.303.
- Large gaps do not decide RH.

## Files (`Riemann/Deep/GapBar/`)
- `DERIVATION.md`: admissibility, Hermitian kernel, twist identity, symmetry check, ceilings, certificates, gaps.
- **Twisted kernel (core):**
  - `gtw.py`: twisted Arb Γ arrays plus IMCHECK;
  - `runcerttw.py`: barrier driver, e.g. `uv run --no-project --with python-flint --with numpy --with scipy python3 runcerttw.py 3.335 72 2 0.8 1.7 0.4`;
  - `batchTW.sh`, `batchTW2.sh`, `batchTW3.sh`: run the grids with an RSS watchdog (290 MB, then 280 MB at 1 s);
  - `scan_deg.sh` (the degree scans);
  - logs: `logs/certTW_r<r>_<κ>.txt|json`, `logs/batchTW.log`, `summarize.py` (writes the table).
- **Records:**
  - `optmid.py` (float, midpoint real P) and `certmid.py` (Arb), e.g.
    `uv run --no-project --with python-flint python3 certmid.py logs/optmid_r2_d80_k3.324.json 3.324 80`;
  - `copt.py` (float, complex P) and `certc.py` (Arb, complex P);
  - logs: `logs/certmid_*.txt`, `logs/certc_*.txt`, `logs/optmid_*.json`, `logs/copt_*.json`.
- **Failed routes, kept:**
  - `gcertc.py` + `runcertc.py` (direct Hermitian certificate), `logs/certC_r2_3.335.txt`;
  - `runcerttwball.py` (ball κ), `logs/certTWB_r2_3.4_3.45_FAIL.txt`.
- **Float:**
  - `ctest.py` (complex matrices), `cscan.py` (closure ceilings), `cland.py` (landscape), `ceil45.py` (r = 4, 5);
  - `mutest.py` (phase family, symmetry test), `twcheck.py` (twist check), `herm.py` (Γ Hermitian check);
  - run, e.g., `uv run --no-project --with numpy --with scipy python3 twcheck.py 2 3.335 14`;
  - logs: `logs/cscan_r123_N14.txt`, `logs/cdeg_scan.txt`, `logs/cland.txt`, `logs/ceil45.txt`.
- **Copies of GapCeil code** (unchanged): gcore.py, gcore_r.py, gcert.py, gcert_r.py, runcert.py, runcert_r.py, certr.py, verify_q.py.

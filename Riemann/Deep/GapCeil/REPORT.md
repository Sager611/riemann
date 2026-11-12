PARTIAL: float ceiling 3.30279 (r = 2); certified in Arb: no amplifier works for r = 1 at 16 κ in [3.2725, 8] or for r = 2 at {3.305, 3.31, 3.35, 3.5, 4, 5}; certified record Λ > 3.30 (r = 2); κ between grid points, large κ and r ≥ 3 remain

# gapceil: the ceiling of Bui–Milinovich's large-gap criterion (node N8)

## Answer
- **The method cannot give much more than 3.2648, but a little more than we had.** Over all amplifiers (every degree and the
  L² closure), all v and ϑ ≤ 1/4, the main-term criterion h1 = c0/(κ²c1) > 1 holds only for κ < κ*_r:
  κ*_1 = 3.27210 (r = 1), **κ*_2 = 3.30279 (r = 2, best of r = 1-4)**, κ*_3 = 3.22609, κ*_4 = 3.14189 (float, converged to 6
  digits; κ*_r falls after r = 2).
- **Certified barrier points (r = 1).** At each of 16 κ (3.2725, 3.273, 3.275, 3.28, 3.29, 3.30, 3.35, 3.40, 3.50, 3.60,
  3.80, 4, 4.5, 5, 6, 8) an Arb certificate proves κ²c1(v) − c0 ≥ 0 for every φ ∈ L²[0,1/4] and every v ∈ ℝ. So h1 ≤ 1 for
  every admissible r = 1 amplifier. The point 3.2725 lies 0.0004 above the float ceiling.
- **Certified record (side result, r = 2).** Λ > 3.30 unconditionally: degree-60 P, r = 2, ϑ = 2499/10000, v = 1.353909,
  h1 = 1.00000835431 ± 2.8·10⁻¹² (`certr.py`). Also certified: Λ > 3.296 (degree 36) and Λ > 3.29 (degree 24).
  The record was 3.2648 (r = 1, degree 14). GapLambda's "r = 2 gave no gain" was a degree-8 artefact: r = 2 overtakes r = 1
  near degree 12. 3.30 is 0.0028 below κ*_2.
- **Certified barrier points (r = 2, follow-up).** No r = 2 amplifier works at κ ∈ {3.305, 3.31, 3.35, 3.5, 4, 5} (any degree, any v, any ϑ < 1/4).
  These run in the weighted Jacobi basis of DERIVATION.md §7. With the record, the r = 2 method is pinned between 3.30
  (certified success) and 3.31 (certified failure), at those points.
- **What remains** (each is needed for BARRIER):
  1. κ between grid points: no certified monotonicity or Lipschitz bound in κ, for r = 1 or r = 2.
  2. κ beyond the grid (8 for r = 1, 5 for r = 2): no analytic large-κ tail.
  3. r ≥ 3: float only (κ*_3 = 3.22609, κ*_4 = 3.14189, κ*_5 = 3.07110); no uniform-in-r bound or r-monotonicity.

## (i) Formulation (DERIVATION.md §1-3)
1. In X_i = ϑx_i and τ1 = at1, τ2 = bt2 (ab is the Jacobian), the 4-term integrand of GapLambda factors as
   U·conj(U'), U = e^{icX1} − e^{icX2}e^{−ic(τ1−τ2)}, c = κπ. Then c_j ∝ Q_j[f] = ∫ f(X+X1+X2)f(X+X3+X4)ρ_j with an explicit
   entire ρ_j, f(U) = P(1−U/ϑ). ϑ appears only in supp f, so ϑ = 1/4 with f ∈ L²[0,1/4] dominates every ϑ < 1/4.
2. With φ = I²f (I = ∫_X^{1/4}), Q_0 = ⟨φ,M0φ⟩ and Q_1(v) = ⟨φ,(v²M0 − 2vMS + MSS)φ⟩. Each operator is multiplication by
   ρ_j(X,0) = ∫∫(v−X−τ1−τ2)^{2j}4sin²(c(τ1−τ2)/2) plus a compact Volterra series Σ Γ_ij (I^i)^T X^p I^j.
3. **Ceiling = top of a pencil.** κ*_r = sup{κ : sup_v λ_max(M0, M1(v)) > κ²}. Polynomials P are dense in the φ-closure,
   so the sup over all degrees equals the L² value. The L² optimiser violates φ(1/4) = φ'(1/4) = 0, i.e. it puts point
   masses δ, δ' at h = y (the top of the mollifier). Polynomials approach it only algebraically: this explains GapLambda's
   ≈ 0.001 gain per degree.

## (ii) Dual certificate (DERIVATION.md §4; `gcert.py`, `runcert.py`, `logs/cert_*.txt|json`)
- **Inequality certified:** κ²M1(v) − M0 ⪰ 0 on L²[0,1/4] for all real v (r = 1). This is an operator inequality c0 ≤ κ²c1
  on the whole admissible class, so BM's criterion fails at that κ for every P, v and ϑ < 1/4.
- **Exact data:** Γ arrays (Taylor data of ρ_j) as 2-D convolutions of derivative arrays, in Arb.
- **Low modes:** Legendre algebra exact; V = degree < 24.
- **High modes:** integration by parts ((Iw)(0) = ⟨w,1⟩ = 0) gives ⟨w,Kw⟩ ≥ −C|w|². Over all runs C ≤ 0.221·δ_D and δ_W ≥ 0.35.
- **Schur step:** Arb Cholesky of a 24×24 matrix.
- **Covering of v:** M0 ⪰ 0, MS − 0.8M0 ⪰ 0, 1.6M0 − MS ⪰ 0 and 2-13 v-cells; 7-29 operator certificates per κ.
- **Runs:** each κ takes 20-51 s and < 200 MB (peak per run in `logs/summary.txt`).

| κ | float sup_v λ/κ² | certified (N = 24 Legendre modes of φ) |
|---|---|---|
| 3.2725 | 0.9996492 | yes: 29 operators (13 v-cells), D = 50, 289 bits, 130 MB, 26 s |
| 3.273 | 0.9992138 | yes: 25 operators (11 v-cells), D = 50, 289 bits, 120 MB, 26 s |
| 3.275 | 0.9974749 | yes: 23 operators (10 v-cells), D = 50, 289 bits, 129 MB, 25 s |
| 3.28 | 0.9931487 | yes: 19 operators (8 v-cells), D = 50, 289 bits, 129 MB, 24 s |
| 3.29 | 0.9845848 | yes: 19 operators (8 v-cells), D = 50, 289 bits, 129 MB, 25 s |
| 3.30 | 0.9761375 | yes: 15 operators (6 v-cells), D = 50, 289 bits, 123 MB, 24 s |
| 3.35 | 0.9355951 | yes: 15 operators (6 v-cells), D = 50, 292 bits, 129 MB, 24 s |
| 3.40 | 0.8977374 | yes: 13 operators (5 v-cells), D = 50, 292 bits, 123 MB, 24 s |
| 3.50 | 0.8293683 | yes: 11 operators (4 v-cells), D = 50, 292 bits, 133 MB, 25 s |
| 3.60 | 0.7697185 | yes: 9 operators (3 v-cells), D = 50, 295 bits, 123 MB, 23 s |
| 3.80 | 0.6720298 | yes: 9 operators (3 v-cells), D = 50, 297 bits, 115 MB, 20 s |
| 4.00 | 0.5966355 | yes: 9 operators (3 v-cells), D = 52, 308 bits, 146 MB, 27 s |
| 4.50 | 0.4655106 | yes: 9 operators (3 v-cells), D = 54, 322 bits, 154 MB, 29 s |
| 5.00 | 0.3709250 | yes: 9 operators (3 v-cells), D = 55, 332 bits, 161 MB, 30 s |
| 6.00 | 0.2334784 | yes: 7 operators (2 v-cells), D = 59, 358 bits, 152 MB, 37 s |
| 8.00 | 0.1210863 | yes: 11 operators (4 v-cells; cell shrink 0.3), D = 67, 411 bits, 193 MB, 51 s |

## (iii) Float ceiling (r = 1 unless stated)
| deg P (ϑ = 1/4) | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | L² closure |
|---|---|---|---|---|---|---|---|---|---|---|---|
| κ_max | 3.183185 | 3.228562 | 3.246264 | 3.255011 | 3.259958 | 3.263028 | 3.265064 | 3.266483 | 3.267512 | 3.268282 | **3.272097** |

- **Convergence.** The closure value is stable from N = 6 to N = 40 φ-modes; 2 modes already give 3.268304.
- **Record gap.** The old record 3.2648 (r = 1, degree 14, ϑ = 0.2499) sits 0.0073 below the r = 1 ceiling and 0.038 below the r = 2
  ceiling. The new record 3.30 sits 0.0028 below the r = 2 ceiling.
- **r = 2, ϑ = 0.2499, κ_max by degree:** 3.198259 (4), 3.248343 (8), 3.269497 (12), 3.280311 (16), 3.286566 (20), 3.290507 (24),
  3.293148 (28), 3.296359 (36), 3.298163 (44), 3.299276 (52), 3.300010 (60).
- **Closures:** r = 2: 3.302788 (ϑ = 1/4); r = 3: 3.226093; r = 4: 3.14189.
- **Ceiling vs mollifier length** (r = 1, same main-term formula; BBLR proves it only for ϑ < 1/4):

| ϑ | 1/8 | 3/16 | 1/4 | 5/16 | 3/8 | 7/16 |
|---|---|---|---|---|---|---|
| κ* | 2.973414 | 3.133791 | 3.272097 | 3.378425 | 3.446444 | 3.478148 |

  Even a twisted fourth moment up to ϑ = 7/16 would lift this route only to ≈ 3.48 (r = 1).
- **Landscape.** sup_v λ(κ)/κ² decreases monotonically: 1.2875 (3.0), 0.99999 (3.2721), 0.8977 (3.4), 0.5966 (4), 0.3709 (5),
  0.2335 (6), 0.1211 (8). λ itself falls from 11.59 to 7.75; the multiplication-part ratio is 6.1-6.6.
  The κ = 5, 6 values agree to 7 digits between the FFT pipeline and the Arb-exact matrices.

## (iv) RH version
Under RH, BM's criterion gives λ > κ with the same pencil and constants, so the same ceiling applies (3.27-3.30). The RH
record 3.63846723 (Durkan–Pearce–Crump 2610.28527) uses a different method: symmetric reduced-fraction coefficients,
a weighted zero-count mean and finite quadratic forms (scout check of their TeX, l.76, l.144, l.792, l.941-968).
It is not this pencil, so it was not analysed (brief).

## Novelty
- `python3 .scratch/lst.py 'large gap|gaps between (consecutive )?zeros|Wirtinger|zero spacing|consecutive zeros' --zenodo '"large gaps" zeta'`
  (Oct-Nov 2026, 7 math categories, Zenodo) finds no work on the limits of BM's method. The hits are 2611.14486 (small
  gaps), 2610.22927, 2610.26038 and 2610.31744 (other gap questions), 2610.28527 (DPC), plus unrelated Zenodo records.
- No ceiling is stated in BM 1410.3635 (l.70 states 3.18; l.237 says only that general amplifiers improve earlier work),
  in BHSJ 2511.14415 (l.126: "no attempt to optimize"; their own L⁴ theorem), in IKT 2510.14309 (only Montgomery–Odlyzko
  limits, l.488-495) or in DPC 2610.28527.
- Our GapLambda NOVELTY.md and REPORT.md estimated the ceiling at "near 3.266-3.27" from degree extrapolation. That
  estimate missed both the closure value and r = 2.
- **New:** the pencil formulation in φ = I²f, the closure ceilings (r = 1-4), the ϑ-dependence and the operator-inequality
  certificate.
- **New record:** Λ > 3.30. No unconditional bound above 3.2648 appears in GapLambda's search (OpenAlex citers of BM,
  arXiv to Nov 2026, Zenodo, HAL) or in this month's listings.
- **Known:** the criterion itself (BM).

## Evidence for the open pieces (float; `evid.py`, `evid_arb.py`, `logs/evid/`)
- **No corner mode.** The r = 2 closure ceiling is 3.3027883 for N = 8, 20, 40 and 60 φ-modes.
- **r = 2 landscape.** sup_v λ/κ² falls monotonically: 0.99412 (3.31), 0.96257 (3.35), 0.92550 (3.4), 0.80012 (3.6),
  0.62848 (4), 0.39567 (5), 0.25713 (6), 0.13769 (8).
- **r = 1 at κ = 10.** The Arb-exact matrices give 0.07634 (cap-limited: κ = 14 hit 252 MB and was stopped).
- **r = 5.** κ*_5 = 3.07110, continuing the fall after r = 2.

## Review
- **Reviewer 1 (barrier argument, pasted inline): GAP, SIGNIFICANT if repaired.**
  - The gap: the degree cutoff for P_WBp was stated as N + 2D instead of N + 2D + 2, and parameters were missing inline.
  - The code already used a larger window: N_w = N + 2D + 4 exceeds the maximal degree N + 2D − 1 (0 ≤ i, j ≤ D − 1, p ≤ 2).
  - The parameters (D, γ, precision, cells, δ_D, C, δ_W, ε) are now in DERIVATION.md §4 and `logs/cert_*.json`.
  - Steps 2-4 were judged structurally sound.
- **Reviewer 3 (weighted Jacobi high-mode bound, r = 2, pasted inline): GAP in step (d) as written, ROUTINE.**
  - Steps (a)-(c) work.
  - The gap: the exact Parseval window holds only for the truncated series B_T. The tail R needs ‖X^{−b}R‖ ≤ ε in every block.
  - The code already does this: ε is an entrywise radius on A_VV, enters the V–W term as (1+1/η)ε², and is summed into C.
    DERIVATION.md §7.5 now says so.
  - t_k uses the shifted bound only for k ≤ N/2, then the global ‖I^k‖ ≤ ϑ^k/k!.
- **Reviewer 2 (r = 2 record at 3.296, pasted inline): GAP, ROUTINE.**
  - The gap: the degree-36 witness P and the Γ formulas were not in the inline text.
  - They are in `logs/record_P_r2_d36.txt` (exact coefficients), `certr.py` and DERIVATION.md §5. The degree-60 witness for
    3.30 is in `logs/record_P_r2_d60.txt`.
  - The reviewer found "the arithmetic implication sound" and classed the result as a constant improvement for RECORDS.md.

## Caveats
- **Grid only.** The barrier is exact at each listed κ. Between grid points and for κ > 8 the evidence is the smooth,
  monotone float pencil. A κ-Lipschitz bound was not achieved: ball-κ and Cauchy bounds lose 10⁵-10¹⁰ through the
  cancelling Taylor series (DERIVATION.md §6).
- **r = 1 and r = 2 at grid points only; r ≥ 3 float only.** For r ≥ 2 the multiplication part X^{r²−1}d vanishes at X = 0.
  Working in L²(X^{r²−1}dX) with the Jacobi basis, I stays tridiagonal (DERIVATION.md §7), so the r = 1 certificate carries over.
  - Correction to the first report: there is no non-compact "Hardy part" at the corner. The kernels are homogeneous of
    degree 0, so the operators are compact like I.
  - The r = 2 ceiling 3.30279 itself is float only.
- **Not analysed:** BM's second criterion h_2 = 4c1/(κ²c2) (Λ* > 4.05 for multiple zeros, l.83, l.329) is the sibling pencil
  (M1, M2); the same machinery applies.
- **Inputs as printed:** BM's Theorem smoothed / Lemma main and BBLR's twisted fourth moment, as in GapLambda.
- **Validation of the transformation:**
  - BM's P gives h1 = 1.002177814 here versus 1.00217781433 in the GapLambda Arb certificate.
  - At ϑ = 0.2499, the best r = 1 amplifiers of degree 0-4 (κ = 3.2645) and r = 2 of degree 0-6 (κ = 3.25) reproduce
    GapLambda's direct 7-D quadrature of BM's printed integrand to all 6 digits (`logs/float/check_r2.txt`).
  - Arb (`certr.py`) and float (`gcore_r.py`) agree to 7-8 digits on the record amplifiers. They are independent codes
    (derivative-array convolutions versus FFT Taylor data).

RH relevance: N8 (zero spacing, Λ = ∞ direction). The ceiling closes BM's Wirtinger route at ≈ 3.27-3.30. Progress
on large unconditional gaps needs a new F (other shifts or several ζ factors), ϑ > 1/4 (worth at most ≈ 0.2), or another method.

## Files
- `DERIVATION.md`: reduction, φ-representation, pencil, certificates, gaps.
- **Barrier certificate (r = 1):** `gcert.py` (Arb operators and the PSD test) and `runcert.py` (v-covering driver), e.g.
  `uv run --no-project --with python-flint --with numpy --with scipy python3 runcert.py 3.2725 24 auto 0.8 1.6`
  (κ = 8 needs the extra argument `0.3`). Logs: `logs/cert_<κ>.txt` and `.json` (each operator's δ_D, C, δ_W, ε and
  Cholesky result), `logs/summary.txt`.
- **Record certificate (r = 2):** `r2opt.py` (float optimiser, writes `logs/r2opt_*.json`) and `certr.py` (Arb), e.g.
  `uv run --no-project --with python-flint python3 certr.py logs/r2opt_r2_d60_k3.30.json 3.30 80`.
  Logs: `logs/certr_r2_d60_k3.30.txt`, `logs/certr_r2_d36_k3.296.txt` and `logs/certr_r2_d24_k3.29.txt`;
  witnesses: `logs/record_P_r2_d60.txt` and `logs/record_P_r2_d36.txt`.
- **r = 2 barrier certificate:** `gcert_r.py` (weighted Jacobi basis) and `runcert_r.py`, e.g.
  `uv run --no-project --with python-flint --with numpy --with scipy python3 runcert_r.py 3.31 72 2 0.8 1.7`.
  κ = 3.305 needs the extra argument `0.4`; κ = 4 and 5 use `48 2 0.8 1.7 0.7 76`.
  `verify_q.py` checks the basis (`logs/verify_q.txt`); `diag_r2.py` explains the D = 52-55 failures.
  Logs: `logs/certR_r2_<κ>.txt|json`, `logs/summaryR.txt`.
- **Float pencil:** `gcore.py` (r = 1, FFT Taylor data), `gcore_r.py` (general r), `scan.py` (degree and closure),
  `scan_r.py`, `scan_r_deg.py`, `scan_th.py`, `landscape.py`, `check_r2.py` (against GapLambda's scans) and
  `test_bm.py` (BM's P), `evid.py` and `evid_arb.py` (open-piece evidence). Logs: `logs/float/`, `logs/evid/`.

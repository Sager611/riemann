BARRIER: the one-piece Levinson–Conrey method for GL(2) cannot give κ_f > 0.07005633 at ν < 25/128 (record 0.0700562)

# gl2ceil — certified ceiling of the one-piece Levinson–Conrey method for GL(2) L-functions (nodes N8, N12)

## Claim
1. Computer-assisted (Arb) barrier for the functional of Bernard's Thm th:mm2 (l.770–783) / AT21 l.205–218: for every R > 0, every real polynomial P
   with P(0) = 0, P(1) = 1, every Q in Bernard's class (any degree) and every ν < 25/128 (Andersen–Thorner's range at θ_KS = 7/64),
   1 − (1/R) log c(P, Q, R, ν/2) ≤ κ*(ν/2) ≤ κ*(25/256) ≤ 0.07005633 (ξ = ν/2; κ* is non-decreasing in ξ, §1.7 of DERIVATION).
   With our record (degree2mol, rigorous, ξ = 25/256 − 10⁻⁷): 0.0700562390 ≤ κ*(25/256) ≤ 0.07005633. Float value of the ceiling
   (exact Euler–Lagrange optimum, H¹ closure): 0.0700563167 at R = 9.867773, P = sinh(λx)/sinh λ, λ = 0.770589, Q(1) = −2.375·10⁻⁵.
   The record is 9.1·10⁻⁸ below the certified ceiling, 7.8·10⁻⁸ below the float ceiling, 5.9·10⁻⁹ below the float optimum at its own ξ.
2. Under the generalised Ramanujan conjecture (ν < 1/4): 0.0897101789 (record) ≤ κ*(1/8) ≤ 0.08971026 (certified); float 0.0897102510.
3. Riemann ζ at Conrey's θ = 4/7 (same functional: KRZ 1609.03738 l.405–407, "precisely the term recovered by Conrey"): the one-piece
   ceiling is ≤ 0.40896 (certified; float 0.4089593260 at R = 1.291586, λ = 1.003074). Conrey 1989 had 0.4088: within 1.6·10⁻⁴ of the
   method's limit. Historical; later ζ records (BCY 0.4105, Feng 0.4128, PRZZ 0.4173) use two-piece or Feng mollifiers.
4. Consequence for N8/N12: no (P, Q, R) of any degree adds more than 9.1·10⁻⁸ to the GL(2) record at ν < 25/128 (8.1·10⁻⁸ under
   Ramanujan). More needs a longer admissible ν (Andersen–Thorner's 1/4 − θ/2 comes from shifted-convolution bounds; AT21 l.223 says so
   heuristically, this proves the quantitative form), a smaller θ, or another mollifier shape.

## Formulation (DERIVATION.md §0–1)
κ*(ξ) = sup over R > 0, P ∈ R[x] (P(0)=0, P(1)=1), Q in Bernard's class of 1 − (1/R) log c. Proved reductions: complex Q never helps
with real P (the |·|² splits into real and imaginary parts); the real class is Q(0) = 1, Q(y) + Q(1−y) ≡ β; c = (1 + w(1)²)/2 +
(C_P/ξ)∫w² + ξB_P∫w′² (w = e^{Ry}Q); c is H¹-continuous and polynomials are H¹-dense in both classes, so "any degree" = the H¹
closure; the inner minimum over P is the sinh family (B_P, C_P enter linearly); C_λ ↑, B_λ ↓ (exchange argument). So
κ*(ξ) = sup_{R>0, λ≥0} 1 − log M(R,λ)/R, M(R,λ) = inf_Q of a positive quadratic functional on an affine H¹ space (its minimiser solves
a linear Euler–Lagrange BVP, 2508.11108's ODE in t = R(2y−1)). κ* is non-decreasing in ξ, so ξ = 25/256 bounds the open range.

## Dual and certificate (DERIVATION.md §2–4; ceil.py)
Lemma D (weak duality): for every continuous piecewise-polynomial σ and every symmetric μ (μ(y) = μ(1−y)), M = ∫μ,
  M(R,λ) ≥ 1/2 + M − 2σ(0) − (2σ(1) + Me^{−R})²/2 − (1/α)∫(σ′ + μe^{−Ry})² − (1/γ)∫σ²,   α = C_λ/ξ, γ = ξB_λ
(pointwise Fenchel + integration by parts + the symmetric pairing 2∫μe^{−Ry}w = βM + minimisation over β). At the Euler–Lagrange
solution σ = γw′, μ = e^{Ry}(αw − γw″) close the gap (floats: gap/J ≤ 10⁻¹³). Cover of (R, λ) ∈ [1/2, 28] × [0, 3] by adaptive cells;
per cell σ, μ (32 pieces, degree 12, double coefficients, σ made exactly continuous) come from the float optimum at the cell centre;
the dual is evaluated in Arb (256 bits; ∫xʲe^{ax} = ₁F₁(j+1; j+2; a)/(j+1)); λ by monotone subcells (1/α ≤ ξ/C_{λ_i},
1/γ ≤ 1/(ξB_{λ_{i+1}}), dyadic grid 2⁻²⁰); R by a third-order Taylor bound with explicit Cauchy–Schwarz remainders (concave lower bound minus
convex target: the two endpoints suffice). Tails by hand: R ≤ 0.5367 (closed form at R = 0: M(0,λ) ≥ 2 − 4ξ/(3ξ+1) = 1.69789, and
d log c/dR ≥ −ξ/√3); λ ≥ 3 and R ≥ 28 (c ≥ 1/2 + e^{2R}K/(e^{2R} + 2K), K = (C/ξ)e^R gd(R)/(4R), increasing in R and C).
Run: `cd Riemann/Deep/GL2Ceil && HLBITS=20 uv run --no-project --with numpy --with scipy --with python-flint python3 ceil.py 25 256 0.07005633`
(Ramanujan: `HLBITS=21 ... ceil.py 1 8 0.08971026 22 3`; ζ: `R0=0.125 ... ceil.py 4 7 0.40896 6 6`)

| ξ (setting) | κ_up | result | cells (dual + crude) | time, peak | log |
|---|---|---|---|---|---|
| 25/256 (θ_KS = 7/64) | 0.070057 | CERTIFIED | 1141 + 19 | 185 s, 83 MB | logs/run_25_256_0.070057.log |
| 25/256 | 0.0700565 | CERTIFIED (HLBITS=17) | 1184 + 19 | 199 s, 85 MB | logs/run_25_256_0.0700565.log |
| 25/256 | 0.0700564 | CERTIFIED (HLBITS=18) | 1211 + 19 | 235 s, 83 MB | logs/run_25_256_0.0700564.log |
| 25/256 | **0.07005633** | CERTIFIED (HLBITS=20) | 1257 + 19 | 198 s, 83 MB | logs/run_25_256_0.07005633.log |
| 1/8 (Ramanujan), RMAX 22 | 0.089711 | CERTIFIED | 905 + 3 | 114 s, 83 MB | logs/run_1_8_0.089711.log |
| 1/8 | 0.0897105 | CERTIFIED (HLBITS=17) | 930 + 3 | 139 s, 85 MB | logs/run_1_8_0.0897105.log |
| 1/8 | 0.0897103 | CERTIFIED (HLBITS=18) | 978 + 3 | 152 s, 86 MB | logs/run_1_8_0.0897103.log |
| 1/8 | **0.08971026** | CERTIFIED (HLBITS=21) | 1014 + 3 | 157 s, 83 MB | logs/run_1_8_0.08971026.log |
| 4/7 (ζ), R0 1/8, RMAX 6, Λ 6 | 0.40896 | CERTIFIED | 429 + 97 | 73 s, 84 MB | logs/run_4_7_0.40896.log |
| 4/7 | 0.4089595 | not certified: margin 2·10⁻⁷ below the λ-grid loss at 2⁻¹⁷ | — | 35 s | logs/run_4_7_0.4089595.log |
Cell lists: logs/ceil_<ξ>_<κ_up>.json. One thread, < 90 MB everywhere.

## Float ceiling vs degree of Q (max over R and λ; Ritz with odd Chebyshev terms; floatceil.py → logs/floatceil.log)
| deg Q | 1 | 3 | 5 | 7 | 11 | 15 | 19 | 27 | 39 | 59 | ∞ (EL ODE) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ξ = 25/256 | −0.13478 | −0.04236 | 0.01146 | 0.03657 | 0.05824 | 0.06582 | 0.06860 | 0.069915 | 0.0700530 | 0.07005631 | 0.07005632 |
| ξ = 1/8 | −0.07260 | 0.00572 | 0.04989 | 0.06939 | 0.08429 | 0.08836 | 0.08941 | 0.0896975 | 0.0897102 | 0.08971025 | 0.08971025 |
AT21's published degree-27 parameters give 0.0696050 / 0.0896557 (their Q coefficients not optimal); degree 1 (simple zeros) is negative.

## Checks
- Duality gap in floats at the Euler–Lagrange optimum: gap/J ≤ 10⁻¹³ at (R, λ) = (9.87, 0.77), (6, 0.77), (14, 1.2), (20, 2), (2, 0.3)
  (gapfloat.py → logs/gapfloat.log); duals generated from a global-polynomial Ritz Q fail at R ≳ 10 (μ = e^{2Ry}L needs exponentially
  small relative errors near y = 1), hence piecewise duals from the Euler–Lagrange ODE.
- Fixed-dual reach (measured): relative gap 0.94·ΔR² in R and 0.0055·Δλ² in λ (reach.py → logs/reach.log); cells are sized accordingly.
- Cross-check against independent code (xcheck.py, imports Degree2Mol/cert.py read-only; logs/xcheck.log): at the record parameters the exact
  record c (cert.py) and our Arb dual bound D at the same (ξ, R, λ) satisfy c ≥ D with (c − D)/c = 5.8·10⁻⁸ (ν < 25/128; κ_D = 0.0700562449 =
  the float optimum at that ξ) and 2.3·10⁻¹⁰ (Ramanujan).
- Coverage (covcheck.py, exact rationals): every certified cell list tiles its core box ([1/2, 28] or [1/2, 22] or [1/8, 6] × [0, Λ]) with
  matching area and no overlapping pair: COVER OK for all certified runs (the failed ζ run: COVER FAIL, as it should).
- Taylor and integration test (taylorcheck.py, logs/taylorcheck.log): at three cells (R = 9.87, 20, 6.02) the direct Arb dual at R_c ± d
  is above the Taylor lower bound used by the certificate for d = 2⁻¹², 2⁻¹⁰, 2⁻⁸ (gap ~ d³), and an mpmath quadrature of the dual
  integrals agrees with the Arb moment evaluation to 10⁻¹⁶.
- Reproducible: ceil.py is deterministic; the float stage only proposes (σ, μ), every inequality is decided in Arb. Peak 83–87 MB, one thread.
- Independent replay (replay.py; logs/replay_*.log): a re-implementation of the whole check (exact-rational integrals by the finite
  formula ∫₀¹e^{ax}G = e^aX − Y instead of ₁F₁ moments, bounds, Taylor assembly, λ-subcells, crude cells, tails, strip) replays every cell
  of the cell lists whose tiling covcheck.py verifies: REPLAY OK for ξ = 25/256 at 0.07005633 (1257 dual + 19 crude cells, 0 failures, tails and strip OK, 219 s, 90 MB) and for ξ = 1/8 at
  0.08971026 (1014 + 3 cells, 0 failures, 180 s). It shares only the untrusted float proposer build().
- Reviewer 1 (decisive argument pasted inline with exact values, .scratch/gl2ceil/reviewer_brief.txt): "VERDICT: CORRECT ... no
  mathematical gap found ... SIGNIFICANCE: SIGNIFICANT — a useful method-specific barrier, placing the supplied record within 7.61·10⁻⁷
  of its ceiling. Scope: mathematical review only; unseen implementation and certificate output were not independently verified."
- Reviewer 2 (final numbers + covcheck/xcheck/taylorcheck, .scratch/gl2ceil/reviewer_brief2.txt): "GAP — implementation caveat
  substantially reduced, not independently closed; no new mathematical objection. SIGNIFICANT barrier; the last-digit tightening alone is
  ROUTINE ... Remaining todo: independent verifier audit and certificate replay, including tails and correspondence between certified
  cells and coverage"; wording: say computer-assisted, state ξ = ν/2 with monotonicity, do not imply other one-piece variants.
  Response: replay.py (above) and the wording of Claim 1. Reviewer 3 (replay evidence): "STILL OPEN — ... does not
  replace the requested soundness audit of replay.py's inequalities and rigorous rounding". Reviewer 4 (code audit, replay.py pasted inline,
  .scratch/gl2ceil/reviewer_brief4.txt): "SOUND under the supplied lemmas/domain assumptions and exact Q conversion ... KP factors,
  reflection, iexp, exponential indices, and all first/second derivative signs are correct ... Floats select witnesses or Rs, not trusted
  bounds ... No false-pass defect found". Net: mathematics CORRECT (reviewer 1), implementation independently replayed and audited SOUND.

## Novelty
- Scouts: python3 .scratch/lst.py 'Levinson|mollifier|critical line|zeros on the critical|proportion of zeros' --zenodo '"critical line"
  "Levinson"' (math.NT, CA, CV, SP, PR, FA, CO; 2026-10 and 2026-11): no GL(2) critical-line paper, no ceiling (hits 2610.08965, 2610.23467
  concern ζ, not this method's limit). TeX searched (optimal|best possible|limit|cannot|upper bound|Euler–Lagrange|variational):
  KRZ 1609.03738, PRZZ 1802.10521 (l.343–345: 100% would follow from θ = ∞; l.363: only an interim optimisation), BCY 1002.4127,
  2508.11108 (l.952–956 leaves the joint (P, S) optimisation open), 2607.00282 (l.268: no optimisation attempted), AT21 (l.223: "any
  significant theoretical improvement that uses Levinson's method comes from extending the permissible range" — a heuristic; we prove
  the quantitative form), Bernard. OpenAlex citers of AT21 (4) and of Bernard (8): no ceiling, no later GL(2) record.
- Not checked in full text: Conrey 1989 (Crelle; not on arXiv) and Farmer 1993 "Long mollifiers" (not on arXiv); Farmer's statement
  concerns θ → ∞, not a fixed-θ ceiling; Conrey §4 optimises P for fixed Q (the sinh step, used here), not Q.
- Our repo: no earlier ceiling for this method (rg in discoveries/, NOGO.md, RECORDS.md, Riemann/Deep/*/REPORT.md).
- Verdict: no statement of the exact optimum or of a certified ceiling of the one-piece Levinson–Conrey functional found; the method
  (Lagrangian dual of the inner quadratic problem + Arb cell cover) is standard, the result is new as far as checked.

## Caveat
- Scope is the proven method: real P (Bernard l.161, AT21 l.194), Q in Bernard's class, one mollifier piece, the mean value c(P,Q,R,ξ)
  at lengths ν < 25/128 (or 1/4 under Ramanujan). Complex P is not covered by these mean value theorems; AT21 l.996's inf_P formula with
  Im(B)² is the infimum over complex P and is ≤ the real-P infimum (harmless for AT21, who restrict to real Q). Degree2Mol's cplx.log found
  numerically that complex (P, Q) do not beat the real optimum; not certified here.
- The ceiling is a statement about this functional, not about κ_f: other mollifiers (two-piece, Feng-type, Hecke-eigenvalue twisted),
  longer ν, or other zero-counting devices (Levinson with the derivative of L, Atkinson-type) are outside it.
- The R → 0 and λ → ∞ tails use crude analytic bounds; they are far from tight but cover their regions with room.

RH relevance: N8, N12 — this caps the critical-line proportion obtainable for every GL(2) L-function by the one-piece Levinson–Conrey
method at the current mollifier length: the record is already at the method's limit, so the node needs a longer admissible ν (shifted
convolution sums beyond Blomer–Harcos) or a new mollifier, not more optimisation. Suggested NOGO line: "gl2ceil: one-piece
Levinson–Conrey for GL(2) is capped at κ_f ≤ 0.07005633 (ν < 25/128), ≤ 0.08971026 (Ramanujan); record within 1·10⁻⁷."

## Files
- ceil.py — certificate (args: ξ numerator, denominator, κ_up [, RMAX, LAM]; env R0, HLBITS, TMAX); writes logs/run_*.log and logs/ceil_*.json
  (cell list). Runs: table above.
- replay.py — independent replay of a cell list (logs/replay_*.log); covcheck.py — exact coverage check (logs/covcheck.log);
  xcheck.py — comparison with Degree2Mol/cert.py's exact record values (logs/xcheck.log); taylorcheck.py — Taylor/integration test.
- floatceil.py — float ceilings and the degree table (logs/floatceil.log); gapfloat.py, reach.py — float duality gap and fixed-dual reach.
- DERIVATION.md — reductions, Lemma D, cell bounds, tails, with TeX line numbers.
- Scratch (not needed): .scratch/gl2ceil/inner.py (κ(R, λ) landscape), reviewer_brief*.txt (the four review briefs), STATE.md.

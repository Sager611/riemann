RECORD: certified h(0.51) > 0.0237842773 (> 0.016352), with an optimised Levinson–Conrey detector inside E1's proven class (ν = 132499/10⁶ < 0.1325, R = 7260151/10⁶, deg Q = 49, deg P = 23; par_051.json, cert_051.log), plus θ = 0.52, 0.55, 0.60 and the uniform constant h > 0.0092120 on (1/2, 1).

# E1Arb — Arb certificate for the optimised-detector h(θ) of E1 (N8; e1arb, 14 Nov 2026)

## 1. The functional and E1's admissible class (reconstruction)
Functional (E1 §4 = e1_constants.py = SantiAudit FULLPROOF l. 24 = E1Rate kappa_opt.py, same formula):
  c(P,Q,R,ν) = 1/2 + (C_P/ν)·∫₀¹e^{2Rv}Q(v)²dv + ν·B_P·∫₀¹e^{2Rv}(RQ+Q')²dv,  C_P = ∫₀¹P'², B_P = ∫₀¹P²  (cross term ∫PP'·∫(e^{2Rv}Q²)'/2 = −1/2);
  κ = 1 − R⁻¹ log c;  h(θ;κ) = (3 + κ − √((1−κ)(9 − κ − 8c(θ))))/4,  c(θ) = 2 − θ/2 − cot(θ/√2)/√2 (Wang); h ↑ in κ and θ.
Constraints E1's proof uses (E1.md §1–§3, SantiAudit FULLPROOF, E1Audit FULLPROOF):
| # | constraint | source | optimised detector (kappa_opt.py) |
|---|---|---|---|
| a | 0 < ν < (3θ−1)/4 strict (Steps 1–5), ν < 1/2 (K-branch, Young's residues) | E1-M | float value taken AT ν = (3θ−1)/4 = 0.1325: open endpoint, NOT admissible (supremum only) |
| b | ∃ η > 0: θ−η > 1/2, ν < (3(θ−η)−1)/4 (compact windows); J > (1+ν−θ+δ)/(2θ−1) fixed | E1 §2 end | automatic for any fixed ν < (3θ−1)/4 |
| c | Q real polynomial, Q(0) = 1, Q(y)+Q(1−y) ≡ 1 (⇒ Q(1) = 0), any fixed degree (E1 used deg 201) | E1 §3, (Emean): Thm 1.2 at Q = 1, 1+x^j | Q = (1−y) + Σx_j(T_{2j+1}(1−2y) − (1−2y)): exact for every x (T_{2j+1} odd) ✓ |
| d | P real polynomial, P(0) = 0, P(1) = 1; R > 0 | E1 §3, SL Thm 1.2 | P = sinh(λx)/sinh λ is NOT a polynomial: value = inf over polynomials, not attained |
| e | κ < 1, c(θ) < 1; Wang 0 < λ < θ < 1 (λ ↑ θ) | Thm 5.1 | κ ≈ 0.095, c(0.51) = −0.130 ✓ |
Verdict on (1): the optimised detector lies in E1's class except at two limits, both suprema: ν at the open endpoint and P = sinh.
Certificate repairs both: ν = (3θ−1)/4 − 10⁻⁶ (rational) and P = S_N(rx)/S_N(r), S_N the odd Taylor polynomial of sinh of degree 2N−1,
r rational (EXP-010's own P has this form, N = 7). No other constraint is touched; κ(ν) is continuous, so the loss is O(10⁻⁶).

## 2. Certificate at θ = 51/100 (cert.py, cert_051.log, par_051.json, par_051.cert.json)
Optimiser (opt.py): for fixed (R, P), c is a quadratic form in Q's free coefficients (one linear solve); for fixed (R, Q) the best P is
sinh(λx)/sinh λ (closed form); alternate, then maximise κ in R (61-point grid + bounded Brent). Reproduces kappa_opt.py exactly at deg Q = 9, 17
(0.0861612746, 0.0946362626 at ν = 0.1325). E1Rate's 0.094636 is NOT converged in deg Q (kappa_opt.py stops at 8 Chebyshev terms):
| deg Q | 17 | 25 | 33 | 41 | 49 | 65 | 81 |
|---|---|---|---|---|---|---|---|
| κ (float, ν = 0.132499) | 0.0946355197 | 0.0950629072 | 0.0950780422 | 0.0950785356 | 0.0950785515 | 0.0950785520 | 0.0950785520 |
Exact parameters (par_051.json): θ = 51/100, ν = 132499/10⁶ < 53/400 = (3θ−1)/4, R = 7260151/10⁶, deg Q = 49 (24 Chebyshev coefficients,
denominators 10¹²), P = S₁₂(rx)/S₁₂(r), r = 769390197/10⁹, deg P = 23. All constraints of §1 checked exactly in rationals (cert.py asserts).
Arb (python-flint 0.9.0, 2048 bits): C_P, B_P exact rationals; ∫₀¹e^{2Rv}F = e^{2R}G₁ + G₀ by finite integration by parts (exact rationals);
balls only for e^{2R}, log, sqrt, cot. No quadrature, no float input to the certified value.
  c(P,Q,R,ν) = 713.274514253319486067643369964 ± 5e−28;  κ = 0.0950785514937877112672558915643 ± 2e−32
  c(51/100) = −0.1300381157880715855894510;  h(51/100; κ) = 0.0237842773156573537489249864635 ± 4e−32
  **CERTIFIED: κ > 0.0950785514, h(51/100; κ) > 0.0237842773** (Arb comparisons of the balls against the rationals; peak footprint 15 MB).
Control (same code, par_control.json = EXP-010 detector, ν = 3/25, R = 177/20, deg Q = 201, 4000 bits): c = 3302.28692320442756563901243594,
κ1 = 0.0844779079697484032032447429603, h(0.51; κ1) = 0.0163521461332344899947500 — E1 §4 / e1_constants.log reproduced to all printed digits.
Gain at θ = 0.51: 0.0237843 vs E1 §4 0.0163521 (+45 %), vs E1's ν = 1249/10000 row 0.0182315 (e1_constants.log), vs E1Rate float 0.0234738.
Lower bound max{0, c, (c+2κ)/3, h} = h here ((c+2κ)/3 = 0.0200397, c < 0).

## 3. θ = 0.52, 0.55, 0.60 (cert_052/055/060.log, par_*.json; same code, 2048 bits, deg Q = 49, ν = (3θ−1)/4 − 10⁻⁶)
| θ | ν | R | κ > (certified) | h(θ;κ) > (certified) | E1 §4 h(θ;κ1) | E1 ν = 1249/10⁴ row | Wang c(θ) | (c+2κ)/3 | LB = max{0,c,(c+2κ)/3,h} |
|---|---|---|---|---|---|---|---|---|---|
| 0.51 | 132499/10⁶ | 7.260151 | 0.0950785514 | **0.0237842773** | 0.0163521 | 0.0182315 | −0.1300381 | 0.0200397 | 0.0237843 (was 0.0182315) |
| 0.52 | 139999/10⁶ | 6.864774 | 0.1004766731 | **0.0380040012** | 0.0268628 | 0.0287282 | −0.0956189 | 0.0351115 | 0.0380040 (was 0.0287282) |
| 0.55 | 162499/10⁶ | 5.888850 | 0.1167173932 | **0.0787094809** | 0.0566885 | 0.0585144 | −0.0005773 | 0.0776192 | 0.0787095 (was 0.0585144) |
| 0.60 | 199999/10⁶ | 4.723602 | 0.1440186685 | **0.1408715866** | 0.1013780 | 0.1031456 | 0.1345543 | 0.1408639 | 0.1408716 (was c = 0.1345543, Wang) |
"was" = best value in our files (E1 §4 / e1_constants.log, Wang); E1Rate's float table (deg Q = 17, at the inadmissible endpoint) had 0.023474 and
0.140866. Converged in deg Q: κ changes by < 3·10⁻¹⁰ from deg 49 to 65 at every θ (opt_sweep.log, opt_051.log).
At θ = 0.60 the margin of h over (c+2κ)/3 is only 7.7·10⁻⁶ and over Wang's c(θ) 0.0063; both are certified separately (cert_060.log balls).

## 4. Falsification (falsify.py, falsify.log; peak footprint 16 MB)
- F1 independent quadrature (mpmath tanh-sinh, 90 digits, from the exact coefficients) of C_P, B_P, W2, W2': agrees with the Arb closed form to
  ≤ 1.3e−35 relative. (At 40 digits it disagreed at 1e−13: Horner on the monomial Q, coefficients ~5.6e29 — a quadrature artefact, not Arb.)
- F2 exact identities: ∫e^{2Rv}Q(RQ+Q') = 0·e^{2R} − 1/2 and ∫PP' = 1/2 as rationals; Conrey's form |P(1)Q(0)|² + (1/ν)∫∫e^{2Rv}|QP'+ν(RQ+Q')P|²
  expanded = certified c to all 25 printed digits.
- F3 24 admissible perturbations (R, r by ±10⁻³ relative, all x_j by ±10⁻⁴ relative): max κ − κ_cert = −5.1e−12 < 0 (local maximum);
  same detector at ν = 0.1324, 0.13249, 0.132499: κ = 0.0950073, 0.0950721, 0.0950786 (κ increasing in ν, loss of the 10⁻⁶ margin ≈ 7e−7).
- F4 precision: G₀ ≈ −2.2e87 against W2 ≈ 79 (≈ 86 digits cancel): at 128 and 256 bits Arb returns a non-finite ball (no false value), from 384
  bits on the balls agree; the certificate uses 2048 (control 4000).
- F5 inadmissible inputs are refused by the exact checks: ν = 53/400 (endpoint), Q + y²/10⁶ (breaks Q(1) = 0 and the symmetry), P·(1+10⁻⁹).
- F6 float optimiser vs Arb: |Δκ| ≤ 3.4e−15 at all four θ (≫ 6 digits). F7 dh/dκ = 0.702, 0.698, 0.686, 0.668 > 0 (h increasing in κ, as used).
- F8 the optimiser's inf over P (√(ab)·coth λ) is 9.0e−16 below the certified polynomial P: the sinh step costs nothing measurable.

## 5. Uniform constant on (1/2, 1) (cert_uniform.log, par_uniform.json, opt_uniform.log)
ν = 124999/10⁶ < 1/8 ≤ (3θ−1)/4 for every θ > 1/2, so this one detector is admissible at every θ ∈ (1/2, 1) (E1 §1's argument with ν = 3/25);
h increases in θ (c'(θ) = ½cot²(θ/√2) > 0), hence liminf S/N ≥ h(θ;κ) ≥ h(1/2;κ) on (1/2, 1). R = 7701341/10⁶, deg Q = 49, deg P = 23:
  κ = 0.0896858962562198102871552399633 ± 5e−32,  **h(1/2; κ) > 0.0092120308** (was 0.0055353 with κ1, 0.0074291 with E1's ν = 1249/10⁴ detector).

## 6. Review
Reviewer subagent (14 Nov 2026, constraints (a)–(e), the exact parameters and checks of §1–§2 pasted inline): **CORRECT.**
(1) fixed degrees 49/23, R > 0, ν < 1/2, exact endpoint and symmetry conditions are admissible; no degree or R bound in the argument; the auxiliary
moments (Q = 1, 1+x^j) need no symmetry; κ < 1, c(θ) < 1 hold. (2) any fixed ν < (3θ−1)/4 is legitimate: e.g. η = 10⁻⁶ gives θ−η = 0.509999,
window threshold 0.13249925 > ν; J fixed large; equality at the endpoint is not. (3) ∂h/∂κ = [1 + (5−κ−4c)/√((1−κ)(9−κ−8c))]/4 > 0 for κ, c < 1;
evaluate at the certified κ ball with outward rounding (done: h is computed as an Arb ball from the κ ball; the decimals are Arb comparisons).

## 7. Prior art for the compared values (scan 14 Nov 2026, 21:01 UTC)
- Ours: E1 §4 / e1_constants.log (Arb): h(θ;κ1) and the ν = 1249/10⁴ row (the "was" values of §3); E1Rate float 0.023474, 0.140866 (deg Q = 17,
  at the inadmissible endpoint; now superseded and certified higher).
- Wang 2609.07918: c(θ) (c(0.6) = 0.1345543; negative below θ0 = 0.5502). SL Zenodo 23134787 v0.03: max{0, c, (c+2κθ)/3, h(θ;κθ)}, κθ ≈ 0.0358,
  θ ≥ 0.527 — below Wang at 0.6, below E1 elsewhere. SL Zenodo 22860012 v0.07 (PDF read: Thm 1.1 N_odd/N ≥ (θ−1/2)/(4eC[q]), C3 = 0.6568;
  (14) max{0, c, (c+2κ)/3}; Table at θ = 3/4 gains ≤ 1.3e−5 over Wang): at θ = 0.6 its κ ≤ 0.0140, so ≤ Wang. Steuding 2002: unspecified constant.
- lst.py, math.NT listings 2026-10 and 2026-11 ('short interval|simple zero|Levinson|critical line|mollif|proportion of zeros|zeros of the Riemann
  zeta'): 19 title hits, none on proportions in short windows (2611.12438: survey of ζ and ζ' zeros, abstract read); Zenodo '"simple zeros"
  "short intervals"': only SL's records 23134787, 23135354 relevant; newest (24227401, 14 Nov) is conditional and off-topic.

## 8. Draft RECORDS.md row (orchestrator; RECORDS.md not edited)
- (N8) Record: explicit proportion of simple critical zeros in every window (T, T + T^θ]: liminf S/N > 0.0092120 uniformly for θ ∈ (1/2, 1)
  (was 0.0074291, E1), and > 0.0237843, 0.0380040, 0.0787095, 0.1408716 at θ = 0.51, 0.52, 0.55, 0.60 (was 0.0182315, 0.0287282, 0.0585144, E1;
  and Wang's c(0.6) = 0.1345543); E1 with an Arb-certified optimised Levinson–Conrey detector (deg Q = 49, ν = (3θ−1)/4 − 10⁻⁶).
  [details](Riemann/Deep/E1Arb/REPORT.md)
Caveat: inherits E1's inputs (SL Zenodo 23134787 v0.03 audited, Wang 2609.07918 unrefereed); a constant record, no new theorem: no discovery, no Lean.

## 9. Files (Riemann/Deep/E1Arb/)
opt.py (float optimiser; opt_0.1325.json = reproduction of kappa_opt.py at the inadmissible endpoint; opt_<ν>.json; opt_051.log, opt_sweep.log,
opt_uniform.log); mkpar.py (float → exact rationals); cert.py (certificate; cert_<θ>.log, par_<θ>.json exact parameters, par_<θ>.cert.json balls;
cert_control.log = E1 §4 reproduced); falsify.py/.log (F1–F9). Peaks: ≤ 28 MB footprint (/usr/bin/time -l; memcap ≤ 23 MB). Runtime: seconds each.

RH relevance: N8

# FULLPROOF — the Mertens power-saving chain, step by step, with source line numbers

Agent mertensaudit, 2026-11-19/20. M78 = Riemann/Deep/Mertens78/PROOF.md, MCH = Riemann/Deep/MertensCH/DERIVATION.md (line numbers l.N).
Verdict per step: HOLDS / HOLDS (note) / FAILS. Notation as in M78: B the half-plane abscissa, H = {Re s > B}, Z(s) = (s−1)ζ(s)/s, u = log|Z|.

## Part I. External inputs (statement checked against the source)

| input | where used | source statement (fetched) | verdict |
|---|---|---|---|
| Rademacher convexity | M78 l.9 (I2), MCH l.9 | Kadiri–Lumley–Ng arXiv 2101.12263 TeX l.396–403 (Lemma bndLschi, citing Rademacher Thm 4): for −1/2 ≤ −η ≤ σ ≤ 1+η ≤ 3/2, ζ(s) ≤ 3·(1+s)/(1−s)·((1+s)/2π)^((1−σ+η)/2)·ζ(1+η) in absolute values, all t | HOLDS: (I2) is exactly this (0 < η ≤ 1/2 ⟺ the outer inequalities; t unrestricted). Original Math. Z. paper not fetched (paywalled). Numerical check (rademacher.py, logs/rademacher.log): max ratio ∣ζ∣/RHS over t ∈ [0.013, 10^4] (step 0.5) ≤ 0.113 for η ∈ {0.05, 0.1} and σ ∈ {−η, 0.5, 0.85, 0.92, 1+η}, ≤ 0.29 for η = 1/2; 0.9994 at s = −1/2 + 0.013i, η = 1/2 — the bound is an equality at s = −1/2 (ζ(−1/2) = −ζ(3/2)/4π), a sign the quotation is exact |
| Herglotz on a half-plane | M78 l.12 (I4) | v ≥ 0 harmonic on {x > 0} ⟹ v(x+iy) = cx + ∫P_x(y−w)dν(w), c ≥ 0, ν ≥ 0, ∫dν/(1+w²) < ∞ (Axler–Bourdon–Ramey, positive harmonic functions on the upper half-space, rotated) | HOLDS (standard; rotation x ↔ y harmless) |
| 7/8 half-plane | M78 l.8 (I1) | Lean, Fayed.QRH.zeta_ne_zero_of_seven_eighths_lt_re (refereed at statement level by hpaudit, TREE N7) | HOLDS (input) |
| θ₁₄₁ half-plane | HP141 l.5–6, l.67–69 | R466 l.10–12: QRT PR 141, L(s,χ) ≠ 0 for Re s > θ₁₄₁ = p/q (exact, 0.8733458561067185768011…), all χ incl. ζ; admitted with standard axioms, three kernels; our static scan clean | HOLDS (input); B′ = 0.87334585610671858 ≥ θ₁₄₁ re-verified with the exact p/q (Part IV) |
| Platt–Trudgian | MCH l.11 (I4) | arXiv 2004.09765 TeX l.61–62 Thm 1: RH true up to height 3 000 175 332 800 (abstract l.46: all zeros with 0 < γ ≤ 3·10^12 have β = 1/2) | HOLDS. Note: the abstract's 3·10^12 alone would NOT cover the Harnack discs (they reach Im s = T′ + 0.52 > 3·10^12); Thm 1's exact height does (T′ + 0.52 < 3 000 175 332 800, asserted in cert_ch.py l.26) |
| Chirre–Helfgott | MCH l.12 (I5) | arXiv 2511.14736v1 Cor. mertensimple, ∣M(x)∣ ≤ 3x/(π·10^10) + 11.39√x for all x ≥ 1; refereed by CHAudit2 (AUDIT PASSED, line 1), open there: Platt's private residue sums above 1.9·10^6 | HOLDS at statement level (not redone; cited) |
| Arb / python-flint | M78 l.14, MCH l.14 | ball arithmetic, 120 bits | HOLDS (two observations on ball construction, Part III) |

## Part II. M78 (half-plane certificate, all x) — line by line

- **Lemma A1 (M78 l.20–29) HOLDS.** u harmonic on H: (s−1)ζ(s) entire, = 1 at s = 1, zero-free on H by (I1); s ≠ 0 on H. For B < σ ≤ 1+η, (I2) applies (σ > B > −η) and gives u ≤ log 3ζ(1+η) + log|1+s|/|s| + α(ℓ − log 2π), α = (1+η−σ)/2 ∈ [0, κ). The affine-in-α comparison with κ(ℓ − log(1+B)) is ≤ 0 at α = 0 (ℓ ≥ log(1+σ) > log(1+B)) and at α = κ (log(1+B) < log 2π); so u ≤ h. For σ > 1+η: u ≤ log ζ(σ) + log|(s−1)/s| ≤ log ζ(1+η) < h (log 3 > 0, |1+s| > |s|, |1+s| > 1+B). So v = h − u ≥ 0 harmonic on ALL of H (zeros of ζ excluded by (I1), not by the bound). The 8/15 = 1/(1+B) at B = 7/8; with B′ the same proof gives A = log 3ζ(1+η) + κ log(1/(1+B′)) (MCH l.67), as coded (bprime.py l.28).
- **Lemma A2 (M78 l.31–35) HOLDS.** With w = t + w′ and evenness of P, Σλ_jP_{x_j}(t+y_j−w) − P_{x₁}(t−w) = LHS(C1) − RHS(C1) at w′ ≥ 0, integrated against ν ≥ 0; the linear term c(Σλ_jx_j − x₁) ≥ 0 by (C2), c ≥ 0. Negative λ_j are harmless here: only the combined kernel's sign enters. Each P_x(·−w) is O(1/(1+w²)), so all integrals converge separately.
- **Harnack case (M78 l.37) HOLDS** ((x₀² − x₁²)w² ≥ 0; ρx₀ = x₀²/x₁ ≥ x₁).
- **Theorem A3 (M78 l.39–44) HOLDS.** v(s₁) ≤ Σλ_jv(s_j) ⟹ −u(s₁) ≤ Σλ_j(h(s_j) − u(s_j)) − h(s₁). λ_j > 0 needs an upper bound for −u(s_j): −u = log|1/ζ| + log|s/(s−1)| ≤ log(ζ(σ)/ζ(2σ)) + log|s/(s−1)| = L₀ (I3). λ_j < 0 needs an upper bound for +u(s_j): u ≤ log ζ(σ) + log|(s−1)/s| = L₁. Directions correct; σ_j > 1 needed for both (all lines ≥ 1.005, after the B′ shift ≥ 1.0034; Part III). log|1/ζ(s₁)| = −u(s₁) + log|(s₁−1)/s₁|.
- **Corollary A4 (M78 l.46–57) HOLDS.** Re-derived each 1/t term: λ_j > 0: (1+κ)(c_j+1+σ_j) + 2c_j from h, + 2 from log|s_j/(s_j−1)| ≤ 1/(t−c_j) ≤ 2/t; λ_j < 0: (1+κ)·2c_j + (c_j+σ_j) = (3+2κ)c_j + σ_j, and log|(s_j−1)/s_j| ≤ 0 (σ_j > 1/2); s₁: σ₁. log(1−u) ≥ −2u on [0, 1/2] justifies log(t−c) ≥ log t − 2c/t for t ≥ 2c. κ log t collects Σλ_j − 1. R ≥ 0, so R/t ≤ R/T₁. Λ ≥ 1 by integrating (C1).
- **Lemma B1 (M78 l.61–62) HOLDS.** (1/2πi)∫_(c) z^s ds/(s(s+1)) = (1 − 1/z)·1[z ≥ 1] (residues 1 and −1/z; for z < 1 shift right; z = 1 gives 0 both ways). Termwise integration: Σ|μ(n)|n^{−c}·∫|ds|/|s(s+1)| < ∞.
- **Lemma B2 (M78 l.64–65) HOLDS.** 1/ζ is analytic on H (1/ζ(1) = 0: no residue at s = 1; s = 0, −1 ∉ H). Horizontal segments: Harnack from σ₀ = 1.02 at height σ − B gives |1/ζ(σ+iT)| ≤ C T^{a′} uniformly for σ ∈ [σ₁, σ₀) (ratio ≤ (σ₀−B)/(σ₁−B), R bounded), a′ = ((σ₀−B)/(σ₁−B) − 1)(1+η−B)/2 at η = 0.05; integrals O(T^{a′−2}) → 0. Re-computed: a′ ≤ 0.759 (B = 7/8, σ₁ ≥ 0.89: 8.6667·0.0875 = 0.7583); a′ ≤ 0.776 (B′, σ₁′ − B′ ≥ 0.015: 8.7769·0.088327 = 0.7752). Note: cert_ch.py --bp builds the B′ certificates with bprime.make() and does not execute bprime.run()'s "a′ < 2" assertion; the analytic bound 0.776 makes this harmless. The aside for 7/8 < σ₁ < 0.89 (Lean qualitative bound) is never used (σ₁ ≥ 0.89, B′: ≥ 0.888, asserted).
- **Lemma B3 (M78 l.67–68) HOLDS.** #{n ∈ (x, u]} ≤ |u − x| + 1; mean of |u − x| over the window is h/4; M(u) = 0 for 0 < u < 1.
- **Lemma B4 (M78 l.70–72) HOLDS.** |∫(s+1)w^s dw| ≤ |s+1|hX^{σ₁} needs σ₁ ≥ 0 (true); the second branch 2X^{σ₁+1}; division by h gives X^{σ₁}·min(1/|s|, τ/(|s||s+1|)), τ = 2X/h.
- **Lemma B5 (M78 l.74–76) HOLDS.** ∫_{T₁}^∞ t^{a−1}min(1, τ/t)dt = (τ^a − T₁^a)/a + τ^a/(1−a) (τ ≥ T₁), τT₁^{a−1}/(1−a) ≤ τ^a/(1−a) (τ < T₁); factor 2 for t < 0; needs 0 < a < 1 (asserted, batch.py l.16).
- **Theorem B6 (M78 l.78–84) HOLDS.** X^{σ₁}/x^θ = K(x), X^{σ₁}τ^a/x^θ = Q(x) with τ = 2x^{1−θ}/c_h + 1 (re-derived). Monotone: σ₁ − θ < 0 because θ ≥ θ* = σ₁ + a(1−σ₁)/(1+a) > σ₁; σ₁ + a(1−θ) − θ ≤ 0 asserted in Arb (cert.py l.141–142); (1 + (c_h/2)x^{θ−1}) and (2/c_h + x^{θ−1}) decrease. h ≤ x for x ≥ x_a asserted (l.144). All-x form: |M(x)| ≤ x ≤ x_a^{1−θ}x^θ on [1, x_a] (trivial bound joins at x_a with no gap); run.py l.25 takes max(log₁₀C, (1−θ)log₁₀x_a) and minimises over the trial x_a.

## Part III. MCH (hybrid contour, the current record) — line by line

- **Lemma H0 (MCH l.31–33) HOLDS.** R(y) = −yΣ_{d>√y}μ(d)/d² − Σ_{d≤√y}μ(d){y/d²}; Σ_{d≥m}d^{−2} < 1/m² + 1/m ≤ 2/m with m > √y; so |R| ≤ 3√y (also for y < 1). |M(u) − M(x)| ≤ Q(u) − Q(x) (u > x) ≤ (6/π²)|u − x| + 6√X; mean over the window gives βh, β = 3/(2π²).
- **Lemma H1 (MCH l.35–39) HOLDS.** Majorant with B = 1/2 re-derived (ℓ > log(3/2); A′ = log 3ζ(1.1) + 0.3·log(2/3) = 3.336 > 0). Harnack in the disc |s − (1.02 + it)| < 0.52: Z has no zero there (Re s > 1/2 and t − 0.52 < Im s < t + 0.52 ≤ T′ + 0.52 < 3 000 175 332 800, Platt–Trudgian Thm 1; Re s ≥ 1 classical); Λ_L = (ρ+r)/(ρ−r) = (2σ₀′ − σ − 1/2)/(σ − 1/2) ✓. Then M78 A3/A4 with one positive mass, c = 0 ✓. Recomputed independently (evalind.py): Λ_L = 1.9714286, a_L = 0.29142857, log C_L = 10.048373 (Arb: 1.9714, 0.29143, 10.048).
- **Lemma H2 (MCH l.42–43) HOLDS.** Γ stays in Re s ≥ 0.85 for |t| ≤ T′ (zero-free: PT for 0 < |t| ≤ T′, ζ(σ) ≠ 0 on (0, 1) for t = 0) and in Re s ≥ σ_H > B for |t| ≥ T′; tails at height U → ∞ as in B2.
- **Lemma H3 (MCH l.45–52) HOLDS.** Three pieces re-derived: (σ_L line) J_L0 + 2C_L∫_{T₁}^{T′}t^{a_L−1}dt ≤ J_L0 + 2C_LT′^{a_L}/a_L; (horizontal) |s|, |s+1| ≥ T′, X^σ ≤ X^{σ_H}, C_L(σ)T′^{a_L(σ)} ≤ C_L*T′^{a_L(σ_L)} because Λ_L decreases, A′ > 0, log(ζ(1.02)/ζ(2.04)) > 0 and σ ≤ σ_H in R (all asserted, cert_ch.py l.30–31); 2 segments × (σ_H − σ_L) × 1/(2π) = (σ_H − σ_L)/π ✓; (σ_H line) G(τ, T′) re-integrated in both cases ✓, G ≤ τ^a(1/a + 1/(1−a)) ✓. The kernel bound needs σ ≥ 0 on Γ (σ ≥ 0.85) ✓.
- **Lemma H4 (MCH l.54–57) HOLDS.** Each factor at its worst endpoint: r, x^{1/2−θ}, x^{σ_L−θ} at x_k; x^{σ_H−θ} at x_k or x_{k+1} by the sign of σ_H − θ; τ (increasing in x) at x_{k+1}, and G and min(1/T′, τ/T′²) increase in τ. Needs 1/2 < σ_L ≤ θ < 1 (asserted l.220).
- **Lemma H5 (MCH l.59–60) HOLDS** (M78 B6 argument; σ_H ≤ θ and σ_H + a(1−θ) − θ ≤ 0 asserted, l.227–228).
- **Theorem (MCH l.62–64) HOLDS.** CH(x)/x^θ = (3/(π10^10))e^{(1−θ)L} + 11.39e^{(1/2−θ)L}, L = log x: convex, maximum on [1, X₀] at an endpoint ✓. Junctions: X₀ = 10^{n₀/10^6} is the same exact Arb number in the CH piece and the first range; ranges contiguous in exact integers n_k; X₁ = last range's end (asserted l.221–224). No gap anywhere on [1, ∞).
- **§4 B′ translation (MCH l.66–71) HOLDS.** σ_j′ − B′ = σ_j − 7/8 exactly (Decimal shift, bprime.py l.16), so (C1)/(C2) are the same inequalities; all σ_j′ ≥ 1.0034 > 1; κ, A, a, x₁ recomputed with B′ in make() (l.25–31); every other Cert attribute that depends on B (x₁, F, the Harnack mass) reads the module global B, which bprime sets to B′ before building the B′ certificates; the 7/8 certificates are built before the import and their tail() does not read B (checked).


## Part IV. Code audit (Mertens78/cert.py, run.py, batch.py; HP141/MertensCH/cert_ch.py, bprime.py — the record's scripts)

| item | code | check | verdict |
|---|---|---|---|
| A1–A4 formulas | cert.py l.14–27 (κ, A, Λ, a), l.28–33 (h, L₀, L₁), l.90–105 (tail = log C_∞) | term by term against M78 l.20–57 and, for B′, MCH l.67 | HOLDS |
| (C2) | cert.py l.52–53 | Σλ_j(σ_j − B) − x₁ > 0 in Arb | HOLDS (margin 0.392 / 0.394 for ML092g8 / ML093g5, my mpmath.iv) |
| (C1) tail | cert.py l.58–69 | bracket ≤ w²π[LHS−RHS] and monotone for \|w\| > max c_j (re-proved; reviewer CORRECT) | HOLDS |
| (C1) on [−Y, Y] | cert.py l.70–89 | 4000 uniform pieces of [−Y₀, Y₀] + geometric pieces to ±Y; pieces share float endpoints; last uniform end = Y₀ exactly (checked) | HOLDS |
| ball construction | cert.py l.80, l.121: arb((lo+hi)/2, (hi−lo)/2) from floats | could in principle miss an ulp of [lo, hi]; counted: 0 of 28 634 (ML092g8), 0 of 18 678 (ML093g5), 0 of 95 000 J_L0 balls miss an endpoint (the radius is rounded up to a 30-bit mag, which absorbs the midpoint rounding) | HOLDS (observation only) |
| J₀ / J_L0 | cert.py l.115–128: t from 0 to exactly T₁ = 1000, widths 0.02 (t < 100), 0.01; sup over the ball of 1/\|ζ\| (abs_lower > 0 asserted) ÷ \|σ+it_left\| | upper Riemann sum; covers [0, 1000] | HOLDS; J_L0 ≤ 19.5963 vs my trapezoid 17.483 (consistent) |
| θ, c_h, x_a | run.py l.9–24, final() l.130–158 | floats only choose θ (then σ₁ + a(1−θ) − θ ≤ 0 asserted in Arb), c_h (exact decimal; any value valid) and x_a = 10^k; h ≤ x_a asserted | HOLDS |
| all-x constant (M78) | run.py l.25–29 | max(log₁₀C, (1−θ)log₁₀x_a) in float; every printed C, C′ of M78 §5 is ≥ the computed value (msanity.log: e.g. 7.908e5 → 7.91e5, 1.838e6 → 1.84e6, 1.007e382 → 1.01e382) | HOLDS |
| RH-strip constants | cert_ch.py l.17–36 | T′ + 0.52 < T_PT, A′ > 0, log(ζ(1.02)/ζ(2.04)) > 0 asserted; Λ_L, a_L, C_L, C_L* as MCH l.35–39, l.50 | HOLDS |
| hypotheses | cert_ch.py l.56–57, l.202–203, l.220–228 | 0 < a < 1, σ_H > σ_L, σ_H ≥ 0.888 (B′), h ≤ x on every range and at X₁, x ≥ 2·10^8, σ_L ≤ θ < 1, s ≤ θ and θ ≥ θ* at X₁ | HOLDS |
| tiling | cert_ch.py l.221–224 | X₀ = first range start, contiguous integer exponents n/10^6, X₁ = last end; independently re-asserted in evalind.py | HOLDS (no gap) |
| H4/H5 terms | cert_ch.py l.200–233 | Ta βc_h, Tb 6x^{1/2−θ}√r, Td, Te with min(1/T′, τ₁/T′²) (falls back to 1/T′ when τ₁ ≤ T′ is not certain), Tc with Garb (both branches, else the τ^a bound) | HOLDS |
| rounding of c | cert_ch.py l.236–244 | 3 digits rounded up, re-compared in Arb against every piece | HOLDS |
| B′ | bprime.py l.14–31 | B′ = 0.87334585610671858 ≥ θ₁₄₁ (exact p/q: excess 3.2·10^−18); D = 82707194664071/5·10^16 exact; module-global B set before the B′ Cert objects are built | HOLDS |
| thresholds (lesson 277) | — | every scan starts at or below its threshold: t from 0 (J), w from −Y (C1), x from 1 (CH, trivial), ranges from X₀ exactly, final from X₁ exactly | HOLDS |

## Part V. Falsification and independent recomputation (own code; logs/)

- **Key lemma |1/ζ(σ₁ + it)| ≤ C|t|^a** (falsify.py, flint acb ζ at 64 bits on t-grids, top-12 local maxima refined by golden section; mpcheck.py repeats the peaks in mpmath, identical to 4 digits). Largest ratios:

| bound tested | line σ | t-range | grid | max ratio (at t) |
|---|---|---|---|---|
| Φ(t) of Thm A3, ML092g8_e05_Bp (the record's final certificate) | 0.918345856 | [0, 10^4] | 0.05 | 1.251·10^−9 (14.1217, next to γ₁ = 14.1347) |
| C_∞(1000)t^a, same certificate (a = 0.106441, log C_∞ = 21.4908) | 0.918345856 | [1000, 10^4] | 0.05 | 8.789·10^−10 (1083.2225) |
| Φ(t), ML093g5_e10_Bp (146/95 ranges of the record) | 0.928345856 | [0, 10^4] | 0.1 | 7.483·10^−7 (14.1208) |
| C_∞(1000)t^a, same (a = 0.132675, log C_∞ = 15.0037) | 0.928345856 | [1000, 10^4] | 0.1 | 4.701·10^−7 (1083.2247) |
| C_L t^{a_L}, Lemma H1 (a_L = 0.291429, log C_L = 10.0484) | 0.85 | [1000, 10^4] | 0.05 | 2.651·10^−5 (1083.2099) |
| C_L T′^{a_L} and C_∞(T′)T′^a at t = 3·10^12, 3·10^12 + 0.37 (mpmath Riemann–Siegel) | 0.85; 0.918346 | 2 points each | — | 6.0·10^−9; 1.3·10^−11 |

  Every peak sits at a zero ordinate (γ₁ = 14.13, γ near 1083.2), as expected; no ratio comes near 1. Lemma A2 itself, with true ζ values (not the L₀/L₁ bounds): max over 160 t (first 40 zero ordinates + 120 spread to 10^4) of v(s₁) − Σλ_jv(s_j) = −4.05 < 0 (ML092g8_e05_Bp). (C1) independently (c1check.py, mpmath.iv, exact float endpoints, own tail at Y = 16): HOLDS for ML092g8 (27 148 intervals) and ML093g5 (16 694); true minimum of π(LHS − RHS) is 1.54·10^−4 (w = 3.10) and 1.62·10^−4 (w = 2.93): the LP certificates are tight but strictly positive.
- **End-to-end recomputation of both record pairs** (evalind.py: mpmath, 40 digits, no flint and no import of the audited code; formulas typed from MCH l.35–60 and M78 l.46–50 with B′; only the float plan (X₀, X₁, ranges, certificate and c_h per range) is taken from cert_ch.py's planner, dumpplan.py; any plan gives a valid bound):

| θ | CH piece at X₀ | max over ranges | final [X₁, ∞) | c (mine) | c (Arb, HP141 log) | claimed |
|---|---|---|---|---|---|---|
| 0.9409 (X₀ = 10^240.25, X₁ = 10^327.10, 1737 ranges) | 15 091.9856 | 15 192.6885 (range 1114) | 15 177.3065 | 15 192.689 | 15 192.6885 | 1.52·10^4 ✓ |
| 0.9382 (X₀ = 10^242.80, X₁ = 10^333.65, 1817 ranges) | 96 607.6198 | 96 825.4855 (range 1146) | 96 776.9360 | 96 825.485 | 96 825.4855 | 9.69·10^4 ✓ |

  Agreement to 10 digits. Independent certificate constants: Λ = 2.205079136 / 2.170725378, a = 0.10644111 / 0.13267488, log C_∞(T′) = 21.46999417 / 14.9866546 (ML092g8_e05_Bp / ML093g5_e10_Bp); RH strip Λ_L = 1.9714286, a_L = 0.29142857, log C_L = 10.048373, log(J_L0 + 2C_LT′^{a_L}/a_L) = 20.347116 (Arb: 20.347).
- **Exact M(x)** (msanity.py): own sieve to 10^7 reproduces OEIS A084237 M(10^k), k ≤ 7; max |M(x)|/(9.69·10^4x^0.9382) = 1.03·10^−5 (x = 1), 4.6·10^−6 on [2, 10^7] (x = 5); for 1.52·10^4x^0.9409: 6.6·10^−5, 2.9·10^−5. OEIS M(10^k) for k ≤ 23 (M(10^23) = 62 467 771 689) and Hurst's extremes to 10^16 (arXiv 1610.08551 TeX l.479–512: M = −35 629 003, 40 371 499; q = 0.571 at 7 766 842 813, −0.525 at 71 578 936 427 177): largest ratio 6.2·10^−10 and 3.8·10^−9. The bound holds (vacuous, as expected: it exceeds √x by a factor ≥ 10^4 everywhere).


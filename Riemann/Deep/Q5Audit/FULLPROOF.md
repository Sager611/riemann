# FULLPROOF — Question 5 of arXiv 2508.11108 (κ(θ) > 2θ/3 on (0, 1/2]), refereed chain

Source: arXiv 2508.11108v1 (only version, API checked 2026-11-19), e-print file `1_Latest_9_Aug.tex` (1178 lines); "l." = its lines.
Every step below was re-derived by hand by the auditor and compared with the TeX line cited.

## A. Statement (l.964-966, item 5 of the enumerate opened at l.943; items at l.945, 954, 956, 962, 964)
"Prove that κ(θ) > 2θ/3 for any θ ∈ (0,1/2], using P(x) = x and Q(y) = Q_R(y) = S_R(2R(y−1/2)) with S_R defined in (eqn: SRT)."
l.966: "It suffices to establish a good lower bound for C₁(R) for R ≥ 2√(3/5)." R is tied to θ by (simplipara) l.372: θ = R⁻¹√(3/5);
§sec: proof restates it (l.724-728): κ(θ) := 1 − log c/R, θ = R⁻¹√(3/5), P = x, Q = Q_R. So θ ∈ (0,1/2] ⇔ R ∈ [2√(3/5), ∞), 2√(3/5) = 1.549193338.
κ here is the value of the formula (a ζ-function statement; θ is the ζ-mollifier length, y = T^θ, l.139). Q5 is a pure inequality about it.

## B. Reduction (each step re-derived; ✓ = TeX agrees)
B1 (kappaprop l.270) κ = 1 − log c(P,Q,R)/R. (eqn cpqroptim l.274, cpqroptimwy l.278) c = 1 + θ⁻¹∫₀¹∫₀¹(w(y)P'(x) + θw'(y)P(x))² dx dy, w = e^{Ry}Q.
    Check: d/dx[e^{Rθx}P(x+u)Q(y+θx)]|₀ = e^{−Ry}(P'(u)w(y) + θP(u)w'(y)), so the e^{2Ry}-weighted Conrey integrand is (P'w + θPw')². ✓
B2 class (l.288-292, Q-funcopt l.289): Q ∈ C¹[0,1] real, Q(0) = 1, Q(y) + Q(1−y) = 1 (hence Q(1) = 0, Q(1/2) = 1/2); P = x (linearP l.285).
B3 (eqn:cpr l.322, eqn:J l.326): expand, ∫PP' = 1/2, ∫ww' = (e^{2R}Q(1)² − 1)/2 = −1/2: c = 1/2 + J, J = ∫₀¹(w²/θ + (θ/3)w'²) dy
    (P = x: ∫P² = 1/3, ∫P'² = 1). ✓
B4 (eqn:S l.331) S(t) = Q(1/2 + t/(2R)), t ∈ [−R,R]: e^{2Ry} = e^{R+t}, dy = dt/(2R), Q' = 2R S', w' = e^{R+t}R(S + 2S'):
    J = (e^R/(2R))∫_{−R}^{R}(e^tS²/θ + (θ/3)R²e^t(S+2S')²) dt (l.334). ✓  (S-funcopt l.339): S(−R) = 1, S(t) + S(−t) = 1, so S(R) = 0.
B5 l.343: ∫e^tSS' = [e^tS²/2] − ½∫e^tS² = −e^{−R}/2 − ½∫e^tS². ✓  Hence (newJ l.346) J = −θR/3 + (e^R/(2R))K(S),
    K = ∫_{−R}^{R}e^t(c₀S² + c₁S'²), c₀ = 1/θ − θR²/3, c₁ = 4θR²/3 (c0c1 l.352). ✓ (θR²/3·(−2e^{−R})·e^R/(2R) = −θR/3.)
B6 fold (eqn:K l.356): t → −t on [−R,0], S(−t) = 1 − S(t), S'(−t) = S'(t): K = ∫₀ᴿ e^t(c₀S² + c₁S'²) + e^{−t}(c₀(1−S)² + c₁S'²). ✓
B7 c₀ = c₁ ⇔ 1/θ = 5θR²/3 ⇔ θR = √(3/5) (simplipara l.372) ✓; then c₀ = c₁ = 4/(5θ) ✓ and K = (4/(5θ))I_R(S) (K1 l.376),
    I_R(S) := ∫₀ᴿ e^t(S² + S'²) + e^{−t}((1−S)² + S'²) dt, boundary data S(0) = 1/2, S(R) = 0 (init l.404).
B8 Constants: θR/3 = √(3/5)/3 = 1/√15; (e^R/(2R))·4/(5θ) = 2e^R/(5θR) = 2e^R/(5√(3/5)) = (2/√15)e^R (5√(3/5) = √15). So
    (Q5)  c(S) = 1/2 − 1/√15 + (2/√15) e^R I_R(S),   1/√15 = √(3/5)/3 = 0.2581989, 2/√15 = 0.5163978.
    Agrees with (simplcPRQ) l.420-424 and (crRdef) l.728 after the EL identity (finalK l.415, re-derived: I_R(S_R) = 1 − e^{−R} − S_R'(0) − ∫₀ᴿe^{−t}S_R). ✓
B9 κ > 2θ/3 ⇔ log c < R(1 − 2θ/3) = R − (2/3)√(3/5) ⇔ c < e^R·ρ, ρ := e^{−(2/3)√(3/5)} = 0.59666600 (c > 0 always). ✓
    Dividing (Q5) by e^R:  κ > 2θ/3 ⇔ (1/2 − 1/√15)e^{−R} + (2/√15)I_R(S_R) < ρ.

## C. The paper's Q_R is the minimiser over H¹ (needed: c(Q_R) ≤ c(S) for the non-C¹ zero extension)
C1 EL (diffeq l.397): (cosh t S')' = cosh t S − e^{−t}/2 ⇔ (sdiffeq l.401) S'' + tanh t S' − S + 1/(1+e^{2t}) = 0, S(0) = 1/2, S(R) = 0
    (mainODEbdry l.433). Re-derived from δI = 0. ✓ The paper derives only the necessary condition (l.390-397: "Suppose K* attains a minimum").
C2 Sufficiency (auditor's proof; also DERIVATION §4): let S_R ∈ C²[0,R] solve C1, S ∈ H¹(0,R) with the same end values, φ = S − S_R ∈ H¹₀.
    I(S) = I(S_R) + 2∫[2cosh t S_R'φ' + (2cosh t S_R − e^{−t})φ] + ∫[2cosh t (φ² + φ'²)]; the middle term is
    2[2cosh t S_R'φ]₀ᴿ + 2∫[−(2cosh t S_R')' + 2cosh t S_R − e^{−t}]φ = 0 (φ abs. continuous, S_R' ∈ C¹). So I(S) ≥ I(S_R), equality iff S = S_R.
    No density argument is needed; DERIVATION's density remark is correct but superfluous.
C3 Existence/uniqueness: a homogeneous solution h with h(0) = h(R) = 0 has ∫2cosh t(h² + h'²) = 0 by the same identity, so h ≡ 0; hence
    f(R) ≠ 0 in (const: C1R) l.468 and the BVP has exactly one solution for every R > 0.
C4 (SRT) l.440 is that solution: it is variation of parameters, S_R = C₁f + g₀ + g₁w₁ + g₂w₂ with v₁ = g₂/(W(1+e^{2u})), v₂ = −g₁/(W(1+e^{2u})),
    f(0) = 0, g₀(0) = 1/2, w_i(0) = 0, so S_R(0) = 1/2 for every C₁. Solving C₁f(R) + g₀(R) + g₁w₁ + g₂w₂ = 0 and multiplying by
    2F⁺(0)e^{R/φ} (φ + 1/φ = √5) gives exactly (const: C1R) l.468 (auditor's algebra ✓). Requires: g₁, g₂ homogeneous solutions and
    W = g₁g₂' − g₁'g₂ — checked numerically by the auditor (indep.py, §F). The odd/even fold is automatic: U(t) := 1 − S(−t) solves the same ODE
    (1 − 1/(1+e^{−2s}) = 1/(1+e^{2s})) with U(0) = 1/2, U'(0) = S'(0), so U ≡ S; Q_R(y) = S_R(2R(y − 1/2)) satisfies (Q-funcopt) on all of [0,1].
C5 Q_R ∈ C^∞ is admissible in the paper's C¹ class (l.288-292); the C¹ → polynomial step (Weierstrass, l.292) is the paper's, not needed for Q5.

## D. Zero extension and one node's R-range (DERIVATION §4 l.~64-67; q5cert.py docstring)
D1 S_a = Q_a(1/2 + t/(2R_a)) on [0,R_a], Q_a polynomial with Q_a(0) = 1, Q_a(1) = 0, Q_a'(y) ≡ Q_a'(1−y) ⇒ S_a(0) = 1/2, S_a(R_a) = 0, and
    I_{R_a}(S_a) uses (1 − S_a)² correctly (fold B6). Extend S = 0 on [R_a, R]: continuous, piecewise C^∞ ⇒ S ∈ H¹(0,R), S(0) = 1/2, S(R) = 0.
D2 On [R_a,R] the integrand is e^t(0+0) + e^{−t}(1+0), so I_R(S) = I_a + e^{−R_a} − e^{−R} (I_a := I_{R_a}(S_a)). ✓
D3 By C2 and B9: κ(θ) > 2θ/3 ⇐ (1/2 − 1/√15)e^{−R} + (2/√15)(I_a + e^{−R_a} − e^{−R}) < ρ
    ⇔ L(R) := (2/√15)(I_a + e^{−R_a}) + b·e^{−R} < ρ, b = 1/2 − 3/√15 = −0.27459667 < 0. ✓
D4 Monotonicity: b < 0 and e^{−R} decreasing ⇒ L increasing in R. So L(R_b) < ρ covers [R_a, R_b]; for R_b = ∞, sup L = (2/√15)(I_a + e^{−R_a})
    (not attained), so (2/√15)(I_a + e^{−R_a}) < ρ covers [R_a, ∞). The 1/R in κ and the e^R in c are absorbed exactly (B9 is an equivalence).
D5 Coverage: θ ∈ (0,1/2] ⇔ R ≥ 2√(3/5) = 1.5491933; nodes [1.5491, 2.0386] ∪ [2.0386, 3.3906] ∪ [3.3906, ∞) ⊇ [1.5491, ∞), contiguous
    (shared endpoints), first node 1.5491 < 1.5491933 (Arb-checked). θ → 0+ (R → ∞) is covered analytically by the last node (D4), not by a scan.


## E. Transfer (not needed for Q5 as stated; needed for the uses)
E1 Q5 is about the value κ(θ) of the construction. A genuine zero proportion needs a polynomial Q: c is continuous in Q for the C¹ norm on the
    affine class {Q(0)=1, Q(y)+Q(1−y)=1} and polynomials are dense there (approximate Q' by polynomials, symmetrise p(y) ↦ (p(y)+p(1−y))/2,
    integrate from Q(0) = 1); the strict inequality c(Q_R) < e^Rρ therefore holds for some polynomial Q (the paper's own remark, l.288-292).
E2 2607.00282v1 (main.tex): class P, Q real polynomials, P(0)=0, P(1)=1, Q(0)=1, Q'(y)=Q'(1−y) (l.545-547); κ(θ) := sup of (kappaprop l.558)
    over P, Q, R (l.565); Cor. fullresultpgl3per (l.587-603) with ϰ₂(ϖ) = κ((1−ϖ)/(2(1+ϖ))), ϖ ∈ (0,1), ϰ₃(ϖ) = κ((1−3ϖ)/(6(1−ϖ))), ϖ ∈ (0,1/3);
    l.604-611 (ourproportionlev): "In [ConShortMollif = 2508.11108] we showed κ(θ) > 0 ... More quantitatively ϰ₂ > (1−ϖ)/(3(1+ϖ)),
    ϰ₃ > (1−3ϖ)/(9(1−ϖ))", i.e. κ(θ) > 2θ/3 at every θ₂(ϖ) ∈ (0,1/2) and θ₃(ϖ) ∈ (0,1/6). 2508.11108 proves this only for θ < min{θ₀, θ₁}
    (Thm precimain l.303-305; l.774 and l.786-790 define θ₀, θ₁ by "there exists", l.307 calls θ₀ "explicitly computable" but gives no value).
    By E1 and sup ≥ construction, Q5 (C-D) supplies exactly the missing step at all lengths of 2607.00282. 2607.00282 uses the INEQUALITY;
    it never names the construction.

## F. Independent numerics (auditor's own code; logs in logs/)
F1 indep_nodes.py (mpmath tanh-sinh quadrature, not the exact-sum formula of q5cert.py) from q5nodes.json: I_a = 1.01189685074579979,
    1.04289774882065615, 1.1097967743420091 = q5cert's Arb balls to all printed digits; L(R_b) = 0.596488884, 0.596539952, 0.590493235;
    ρ − L = 1.771e−4, 1.261e−4, 6.173e−3 (same as q5cert). Q(0) = 1, Q(1) = 0, Q(1/2) = 1/2, Q'(.3) = Q'(.7) at 40 digits.
F2 indep_bvp.py (scipy solve_bvp of (mainODEbdry), c via (crRdef) and via the direct functional; they agree to ≤ 2·10⁻¹²): κ(θ) − 2θ/3 > 0 at all
    27 grid points in (0,1/2] (22 by direct BVP: θ = 1/2, 0.49999, near both node boundaries 0.37997 and 0.22845, ..., 0.02; 5 via the R → ∞ limit, θ = 0.01 ... 10⁻⁴:
    I_R + e^{−R} → I* = 1.1414593260695, constant to 1.6·10⁻¹⁵ over R = 20, 30, 40). Margins: θ = 1/2: 2.510·10⁻² (7.5%); θ = 0.38: 2.240·10⁻²;
    θ = 0.2285: 7.73·10⁻³; θ → 0: relative margin → 2.357% (κ/θ → 0.682381 vs 2/3; paper l.771: lim(1 − S'(0) − ∫) ≤ 1.1421, ours 1.1414593; l.743, l.758: ≥ 0.248, ≥ −0.39006, ours 0.2485939, −0.3900533).
    Tolerance 10⁻⁷ vs 10⁻⁹ changes κ by ≤ 2·10⁻¹¹.
F3 Paper's table (l.828-842) vs F2: κ(1/4) = 0.17638, κ(1/6) = 0.11467, κ(1/8) = 0.085450, κ(5/54) = 0.063196, κ(1/100) = 0.0068238,
    κ(1/500) = 0.0013648: all agree with the printed digits; κ(1/2) = 0.358437 and κ(2/3) = 0.446367 disagree with the printed 0.334…, 0.364…
    (as DERIVATION §4 says). Rigorously: node 1 zero-extended to R = 2√(3/5) (D2 with q5cert's I_a) gives κ(1/2) ≥ 0.358428 (float evaluation of the certified bound), so the printed 0.334… is wrong.
F4 indep_closed.py (mpmath hyp2f1, 30 digits): g₁, g₂ of (SRT) solve the homogeneous ODE (residual ≤ 10⁻²⁹), W(u) of (spec: Wronsk) equals
    g₁g₂' − g₁'g₂ (and W cosh u = √5/2, Abel), the paper's f'(0), g₀'(0), g₁(0), g₂(0), v₁(0), v₂(0) reproduced, C₁(R) of (const: C1R) gives
    S_R(R) = 0 (≤ 2·10⁻²⁹) and S_R'(0) equal to the BVP value to 10 digits at θ = 1/2, 1/4, 1/10. So (SRT) is the BVP solution (C4).
F5 q5cert.py rerun on a copy (.scratch): log byte-identical to Degree2Mol/logs/q5cert.log; 1.1 s, 35 MB max RSS. flintsem.py: arb(fmpq) is
    an enclosing ball, exp encloses, strict < is True only when certain, fmpq_poly(poly) is composition (used by Q_cheb).


AUDIT PASSED — Question 5 of arXiv 2508.11108 holds as claimed: every step (reduction, minimiser, zero extension, R-coverage, certificate) was re-derived or recomputed independently and none fails; only wording slips remain.

# q5audit — referee of our claim that Question 5 of arXiv 2508.11108 holds (discoveries/gl2-critical-line.md; Degree2Mol/q5cert.py)

Sources: 2508.11108v1 (only version, arXiv API 2026-11-19), `1_Latest_9_Aug.tex`; 2607.00282v1 `main.tex`. Full chain with TeX lines: [FULLPROOF.md](FULLPROOF.md).
Code: indep_nodes.py, indep_bvp.py, indep_closed.py, flintsem.py (own code, none from Degree2Mol); logs/.

## 1. Item-by-item (claim and method sentences of the page)

| # | Page statement | Verdict | Evidence |
|---|---|---|---|
| 1 | Question 5 is "κ(θ) > 2θ/3 for every θ ∈ (0, 1/2]" for P = x, R = √(3/5)/θ | HOLDS | Q5 verbatim l.964 (item 5 of l.943-966): also fixes Q = Q_R = S_R(2R(y−1/2)), S_R from (SRT) l.440; R from (simplipara) l.372, l.724-728 |
| 2 | Question 5 is true | HOLDS | FULLPROOF B–D; reviewer CORRECT; q5cert rerun identical; independent node check F1 |
| 3 | 2508.11108 proves it only for θ < min{θ₀, θ₁}, defined, not evaluated (l.307, 774-790) | HOLDS | Thm precimain l.303-305; θ₀ l.774, θ₁ l.786, conclusion l.790; l.307 "explicitly computable", no value |
| 4 | 2607.00282 uses it at all its lengths | HOLDS (wording) | main.tex l.604-611 asserts κ(θ) > 2θ/3 at θ₂(ϖ) ∈ (0,1/2), θ₃(ϖ) ∈ (0,1/6) citing 2508.11108, with κ := sup over P, Q, R (l.565). It uses the inequality, not the construction; sup ≥ construction (E1) |
| 5 | c(S) = 1/2 − √(3/5)/3 + (2/√15)e^R I_R(S) | HOLDS | B1–B8 re-derived; = (eqn:cpr), (newJ), (K), (K1), (simplcPRQ), (crRdef) l.322-424, 728 |
| 6 | near-optimal S on [0,R_a] extended by 0 is admissible in H¹ | HOLDS | D1; it is NOT in the paper's C¹ class (S' jumps at R_a) nor polynomial, and need not be: it is only an H¹ competitor in the convexity inequality C2 (I(S_R) ≤ I(S), no density argument needed); the mean value theorem is applied to Q_R via polynomial approximants (E1) |
| 7 | ... "covers every R ≥ R_a" | WORDING | admissible at every R ≥ R_a, but proves the inequality only on [R_a, R_b]; nodes 1, 2 alone fail for large R (L(∞) = 0.6322, 0.6058 > ρ = 0.5967). q5cert.py does use R_b correctly |
| 8 | three nodes 1.5491, 2.0386, 3.3906 cover θ ∈ (0, 1/2] | HOLDS | D4-D5: L increasing in R (coefficient 1/2 − 3/√15 < 0); 1.5491 < 2√(3/5) = 1.5491933; contiguous; last node to ∞ covers θ → 0+ analytically |
| 9 | Certificate q5cert.py (1 s) | HOLDS | rerun on a copy: byte-identical log, 1.1 s, 35 MB (logs/q5cert_rerun.*) |
| 10 | Checks: Q5 at l.964; table l.823-830 | HOLDS | rows 2/3 and 1/2 are l.828, l.830 |
| 11 | Side: table entries κ(2/3), κ(1/2) wrong (0.44636, 0.35843 vs 0.364…, 0.334…) | HOLDS | F3: own BVP 0.446367, 0.358437; rigorous κ(1/2) ≥ 0.358428 from node 1 |

## 2. Derivation referee (details: FULLPROOF A–E)
- Conrey's integrand (l.274-278) → (P'w + θPw')², expansion, ∫PP' = 1/2, ∫ww' = −1/2 (Q(1) = 0 forced by Q(0) = 1, Q(y) + Q(1−y) = 1): c = 1/2 + J. ✓
- Substitution y = 1/2 + t/(2R), the integration by parts of ∫e^tSS', c₀ = 1/θ − θR²/3, c₁ = 4θR²/3, the fold with S(−t) = 1 − S(t), c₀ = c₁ ⇔ θR = √(3/5),
  and the constants θR/3 = √(3/5)/3 = 1/√15 and (e^R/(2R))·4/(5θ) = (2/√15)e^R: all ✓ (B3–B8).
- κ > 2θ/3 ⇔ (1/2 − 1/√15)e^{−R} + (2/√15)I_R(S_R) < ρ = e^{−(2/3)√(3/5)}: an equivalence, so e^R and 1/R are absorbed exactly (B9). ✓
- The paper only derives the Euler–Lagrange equation as a necessary condition (l.390-397). Sufficiency over H¹ holds by the identity
  I(S) − I(S_R) = ∫2cosh t(φ² + φ'²) ≥ 0 (C2); the BVP has a unique solution for every R (C3); (SRT) with (const: C1R) is that solution (C4, F4).
- Zero extension: I_R = I_a + e^{−R_a} − e^{−R} (integrand e^{−t} on [R_a, R]) ✓; monotone increasing bound, worst case R_b ✓; θ → 0+ by the last
  node for all R ≥ 3.3906, no finite scan involved ✓. Q_R(y) = S_R(2R(y−1/2)) needs S_R on [−R,0]: the analytic solution satisfies
  S(t) + S(−t) = 1 automatically (C4). ✓

## 3. Code audit (Degree2Mol/q5cert.py, q5nodes.json; read line by line)
- Inputs exact: R_a, R_b, h_n are decimal strings → Fraction → fmpq (exact). h_last is reset to 1/2 − Σ others (shift ≤ 9·10⁻¹⁴), so Q(1) = 0
  exactly; asserts Q(0) = 1, Q(1) = 0, Q'(y) ≡ Q'(1−y) as exact polynomial identities (fmpq_poly call = composition, flintsem.log). ⇒ S_a(0) = 1/2.
- I_exact: t = 2R(y − 1/2), dt = 2R dy, e^{±t} = e^{∓R}e^{±2Ry}, S'² = Q'²/(4R²) ✓; ∫_{1/2}^1 e^{ay}G = Σ_k (−1)^k (e^aG^{(k)}(1) − e^{a/2}G^{(k)}(1/2))/a^{k+1}
  (finite, a = ±2R, sign of a handled by ap *= a) ✓; the two sums are exact rationals; only exp is rounded, in Arb at 4000 bits, outward.
  Independent tanh-sinh quadrature reproduces all three I_a to all 15 printed digits (F1).
- Inequality: lhs = (2/√15)(I_a + e^{−R_a}) [+ (1/2 − 1/√15 − 2/√15)e^{−R_b}] < rhs = e^{−(2/3)√(3/5)}, tested with Arb strict < (True only if
  certain) ✓; exactly D3–D4 at the worst point R_b of each node ✓. Chain contiguity asserted (R_b of node k = R_a of node k+1) ✓; "PROVED" needs the
  last node open-ended ✓. Margins ρ − L: 1.771·10⁻⁴, 1.261·10⁻⁴, 6.173·10⁻³ (small but ≫ the true ball radii ≤ 7·10⁻¹⁰²⁷, q5cert_radii.py, logs/q5cert_radii.log).
- Lesson 277: the first node 1.5491 ≤ 2√(3/5) = 1.5491933 (Arb: True). Weakness: this check is printed but not part of the "PROVED" gate
  (ok is not and-ed with it); harmless here because it prints True. Docstring cites Question 5 at "l.958-960"; it is l.964.
- The info column "kappa(theta_a)" uses c(S_a) at R_a only; it is not part of the proof.

## 4. Falsification (own code; FULLPROOF F)
- No counterexample: κ(θ) − 2θ/3 > 0 at 27 θ in (0, 1/2] incl. θ = 1/2, both node boundaries (0.37997, 0.22845) and θ = 0.1, 0.01, 0.001, 10⁻⁴
  (indep_bvp.log). Smallest relative margin 2.357% (θ → 0, κ/θ → 0.682381); at θ = 1/2 the margin is 0.0251 (7.5%).
- The paper's own limits (l.743, 758, 771: ≥ 0.248, ≥ −0.39006, ≤ 1.1421) and its constants f'(0), g₀'(0), g₁(0), g₂(0), v₁(0), v₂(0) agree with ours;
  its table agrees at 1/4, 1/6, 1/8, 5/54, 1/100, 1/500 and is wrong at 2/3 and 1/2 (as our page says).
- The closed form (SRT) is the BVP solution (F4): homogeneous residual ≤ 10⁻²⁹, Wronskian exact, S_R(R) = 0, S_R'(0) = BVP value to 10 digits.
- 2607.00282 does not define the construction; its κ is the sup (l.565), so Q5 + density (E1) gives its (ourproportionlev) at every ϖ.

## 5. Reviewer (one, verdict-only, chain B9 + C2 + D1–D5 + the certified numbers pasted inline)
CORRECT: algebra of (1) from c_R; minimality over H¹ (EL cancels the cross term; a derivative jump does not prevent H¹); I_R(S) = I_a + e^{−R_a} − e^{−R};
L' = −d e^{−R} > 0 so the right endpoint covers each finite interval; last interval L < 0.5904932 < ρ for all finite R ≥ 3.3906; θ → 0+ covered.

## 6. Dependents (rg 'question 5|q5cert|2θ/3|2508.11108|2607.00282'); none is affected since the audit passed
- discoveries/gl2-critical-line.md (Claim; Method; Certificate; Review) — suggested wording: "covers every R ∈ [R_a, R_b] (all R ≥ R_a for the
  last node)" instead of "covers every R ≥ R_a"; "which 2607.00282 uses" → "an inequality 2607.00282 (l.604-611) uses".
- RECORDS.md:30 ("arXiv 2508.11108's Question 5 settled"); Degree2Mol/REPORT.md item 5 (l.16-17), NOVELTY.md l.4, l.24; TREE.md:454 and
  Riemann/Tree/R418.lean l.11-13 (Question 5 l.964 verbatim, θ₀/θ₁); COROLLARIES.md:47 (table erratum: confirmed; its twist-family shares use cert.py
  at fixed lengths and do not depend on Question 5).

## 7. Not checked / out of scope
- The mean value theorems (Conrey's Theorem 2, the ζ asymptotic) are inputs of the paper's κ, not of Question 5 (an inequality about the formula).
- The paper's Figure diff_2third_act (PDF data not extracted). Novelty was not re-audited (not in the brief).

## 8. Process notes
- Memory: a first indep_bvp.py run (solve_bvp, tol 10⁻¹⁰, 300 000 nodes, a nonlinear functional component) reached 1467 MB before memcap killed it, and a
  large-R test (tiny2) reached 646 MB max RSS in 4 s, under memcap's sampling interval. Fixed: linear BVP, ≤ 5000 nodes, RSS guard 230 MB,
  R ≤ 40 directly and the R → ∞ limit beyond; all later runs ≤ 93 MB. Lesson: test the largest case of a solver alone first, with a node cap.
- Timeline 22:27–22:49 UTC (well inside 1.5 h); downloads in .scratch/q5audit deleted at the end.

## 9. Files
REPORT.md, FULLPROOF.md (chain with TeX lines); indep_nodes.py, indep_bvp.py, indep_closed.py, flintsem.py, q5cert_radii.py; logs/ (q5cert_rerun.*, indep_*.log/.err, flintsem.log, q5cert_radii.log).

RH relevance: N8, N12

RECORD: at least 7.005% of the zeros of the L-function of every holomorphic or Hecke–Maass newform on Γ₀(N) lie on the critical line (was 6.96%, Andersen–Thorner 2021); 8.971% under the generalised Ramanujan conjecture (was 8.96%).

# degree2mol — joint (P, Q) optimisation of Levinson–Conrey for GL(2) L-functions (nodes N8, N12)

## Claim
1. κ_f > 0.07005 (certified κ ≥ 0.0700562390188) for every holomorphic or Hecke–Maass newform f on Γ₀(N), any N, unconditionally.
   Input: Andersen–Thorner's mean value theorem (arXiv 2004.03581 Thm thm:nu: ν < 1/4 − 7/128 = 25/128) with Bernard's
   zero-counting bound; we use ν = 2ξ, ξ = 25/256 − 10⁻⁷.
2. κ_f > 0.08971 (certified ≥ 0.0897101789336) under the generalised Ramanujan conjecture for the Maass spectrum (ν < 1/4).
3. At Bernard's unconditional length ν < 5/27 (his Thm th:mm2; holomorphic f, squarefree level): κ_f > 0.06642 (2508.11108: 6.32%
   with P = x; Bernard: 2.97%).
4. Families of arXiv 2607.00282 (July 2026; its Cor. fullresultpgl3per and simplpergl2 are stated with κ(θ) = sup over
   (P, Q, R), main.tex l.565, l.587-603): twists L(s, Π×χ), |γ| ≤ Q^ε, on average: PGL(3) self-dual Π: 0.1199 − O(ε) on the line
   (stated 1/9); PGL(2): 0.3658 − O(ε) on the line (stated 1/3) and 0.3562 − O(ε) simple and on the line (stated 1/3). Certified
   κ(1/6⁻) ≥ 0.1199187386, κ(1/2⁻) ≥ 0.3658158115, and with linear Q κ(1/2⁻) ≥ 0.3562685057.
5. Question 5 of 2508.11108 (l.964) is true: κ(θ) > 2θ/3 for every θ ∈ (0, 1/2] for their construction (P = x, Q = Q_R,
   R = √(3/5)/θ); three certified inequalities (q5cert.py). 2607.00282 l.605-611 uses exactly this inequality for all ϖ
   ((1−ϖ)/(3(1+ϖ)) and (1−3ϖ)/(9(1−ϖ)) are 2θ/3 at its two lengths) while 2508.11108 proved it only for θ below a non-explicit θ₀;
   our proof closes that gap.
6. Questions 1-3 of 2508.11108: the joint optimum is P = sinh(λx)/sinh λ (λ ≈ 0.77 at the GL(2) lengths) with S solving their ODE for
   general c and free β (Q(1) ≈ −2·10⁻⁵, not 0); at Conrey's θ = 4/7 it gives 0.408959 vs 0.408770 with Conrey's Q (Question 3: +0.0002).
7. Correction: 2508.11108's table (l.828, l.830) prints κ(2/3) = 0.364…, κ(1/2) = 0.334…; their construction gives 0.44636 and
   0.35843 (their BVP (mainODEbdry) solved by shooting; their (simplcPRQ) and the direct functional agree to 10 digits; our Ritz agrees).
   Rows θ ≤ 1/4 are right. Their closed form (eqn: SRT) is right: evaluated with mpmath it matches the shooting solution (S'(0) =
   −0.4415648493 at θ = 1/2) and gives 0.3584369085 (closedform.py). So only the two table entries are wrong; cause unknown (Levinson's
   Q = 1 − y at the same R gives 0.3337 and 0.4366, so not that). The paper's "numerical evidence" for Question 5 was thus far stronger than printed.
8. No simple-zero share for a fixed GL(2) L-function at these lengths: linear Q (needed, Heath-Brown) with optimal P, R gives κ < 0
   (−0.135 at ν = 25/128, −0.073 at ν = 1/4), as Bernard l.83 and KRZ19 l.1783 expected.

## Before
| setting (ζ-normalised ξ = ν/2) | best published | by | ours | gain |
|---|---|---|---|---|
| fixed GL(2), uncond., ν < 25/128 | 0.0696 (their parameters: 0.0696049) | Andersen–Thorner 2021 (optimal P, Q of degree 27) | 0.0700562 | +0.00045 |
| fixed GL(2), Ramanujan, ν < 1/4 | 0.0896 (their parameters: 0.0896556) | same | 0.0897101 | +0.00005 |
| fixed GL(2), uncond., ν < 5/27 | 0.0632 (P = x) / 0.0297 | 2508.11108 / Bernard 2015 | 0.0664220 | +0.0032 |
| PGL(3) twist family, ξ → 1/6 | 1/9 (from 0.114663) | 2607.00282 | 0.1199187 | +0.0088 |
| PGL(2) twist family, ξ → 1/2 | 1/3 (all zeros; from 2θ/3) | 2607.00282 | 0.3658158 | +0.032 |
| PGL(2) twist family, simple, ξ → 1/2 | 1/3 (Levinson Q, R = 1.3: 0.3428) | 2607.00282 | 0.3562685 | +0.023 |
2508.11108's "8.54%" is at ν = 1/4, i.e. conditional on Ramanujan (l.849); both its GL(2) figures were below AT21's.
Contributions at ν = 25/128: paper's construction (P = x, β = 1, R = √(3/5)/ξ) 0.066660; best R 0.069441; β free 0.069443;
optimal P 0.070053; joint 0.070056 (logs/opt2.log).

## Method
c(P,Q,R) (Conrey/Bernard formula) = (1 + e^{2R}Q(1)²)/2 + (∫P'²/ξ)∫w² + ξ∫P²∫w'², w = e^{Ry}Q (DERIVATION §2). For fixed P it is a
positive quadratic in Q; for fixed Q the best P is sinh(λx)/sinh λ, λ = ξ√(∫w'²/∫w²) (§3). Ritz in Q (30 odd Chebyshev terms, degree 59,
β free; 1 term for simple zeros), alternated with the P-step, golden-section search in R (opt2.py, untrusted). Converged: 30, 40 and 60
terms agree to 10⁻⁸ (logs/sanity.log). Then R, λ, h_n are rounded to rationals, P = degree-25 Taylor polynomial of sinh(λx) normalised
to P(1) = 1, and c is evaluated exactly (x-integral exact over Q, y-integral by the finite exp-polynomial formula; one exp and one log in
Arb, 4000 bits). Question 5: for P = x, β = 1, ξR = √(3/5), c = 1/2 − √(3/5)/3 + (2/√15)e^R I_R(S); a near-optimal S at R_a extended
by 0 beyond R_a covers [R_a, R_b] with one inequality; nodes R_a = 1.5491, 2.0386, 3.3906 cover [2√(3/5), ∞) (DERIVATION §4).

## Certificate
`cd Riemann/Deep/Degree2Mol && uv run --no-project --with python-flint --with mpmath python3 cert.py && uv run --no-project --with python-flint --with mpmath python3 q5cert.py`
(cert.py 52 s, 54 MB; q5cert.py 1 s; single-threaded). Output: logs/cert.log, logs/q5cert.log; exact R, λ, h_n printed (section C)
and stored in params.json, q5nodes.json.

## Checks
- Same Arb code reproduces: AT21's tables 0.0696049654, 0.0896556716 (claims 0.0696, 0.0896); 2508.11108's construction 0.0631964,
  0.0854501 (table 0.0632, 0.0854); Levinson θ = 1/2: 0.34276; Conrey's Q (2508.11108 l.974) with optimal P at θ = 4/7: 0.40877 (Conrey 0.4088).
- cert.py asserts P(0) = 0, P(1) = 1, Q(0) = 1, Q'(y) ≡ Q'(1−y) exactly, ξ strictly inside each range, deg Q = 1 for the simple-zero
  entry. AT21's closed form inf_P c is ≤ our exact c by < 10⁻¹⁸ relative; 2508.11108's (kappa: general) with K(S) by quadrature agrees
  to 10⁻⁵⁶ (formula re-derived in DERIVATION §4: correct).
- Complex Q in Bernard's class: with real P it can only raise c (DERIVATION §3); with complex P too, random starts return to the real optimum.
- κ(θ) = sup is non-decreasing in θ (inf_P c = … coth(ξ√(W₁/W₀)) decreases in ξ), so certifying at ξ = bound − 10⁻⁷ covers the limits.

## Novelty
NOVELTY.md: arXiv listings (2026-10/11), arXiv API, OpenAlex citers of AT21 (4) and KRZ19 (7), Zenodo, repo; 2508.11108 and
2607.00282 have one version each and no follow-up. No later GL(2) record, no joint optimisation, Question 5 open. HAL not checked.

## Caveat
- The fixed-form gains are small (+0.045 and +0.005 percentage points) and come from a standard optimisation; a referee may call them
  routine. Most of the gain is convergence in the degree of Q (degree 27 with our Ritz already gives 0.069915 vs AT21's 0.069604).
- All numbers rest on the quoted theorems (AT21 Thm thm:nu, Bernard's zero-count proposition, 2607.00282's corollaries with κ = sup);
  none is re-proved here. The family items assume 2607.00282's mean value theorems hold for its general P, Q (stated so at l.545, l.565).
- 7.0056% looks like the ceiling of the one-piece Levinson–Conrey method at ν = 25/128 (numerically converged, not proved): more needs a
  longer mollifier (better shifted-convolution bounds) or several pieces.

RH relevance: N8, N12

## Files
- cert.py, params.json, at21_tables.json: certificate for items 1-4 (logs/cert.log)
- q5cert.py, q5nodes.json: Question 5 (logs/q5cert.log)
- opt2.py, opt1.py, q5.py: untrusted optimisers (logs/opt2.log, logs/family.log); odecheck.py, closedform.py, tablediag.py: table check
  (logs/odecheck.log, logs/closedform.log, logs/tablediag.log); sanity.py (logs/sanity.log); cplx.py (logs/cplx.log)
- DERIVATION.md (formulas and constants with source lines), NOVELTY.md

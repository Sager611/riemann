CM WRONG: Connes–Moscovici PNAS 2022 Thm 5.1, Eq. 33, N(E) = (E/2π)(log(E/2π) − 1) − (log E)/2π + O(1), is false — the −(log E)/2π term is absent: the computed count residual tends to −1/4 with fitted log coefficient b = 2·10⁻⁷ ± 6·10⁻⁸ (needed −0.159155) up to t = 2·10⁵, and the term is exactly what Heywood's formula gives when misapplied at the critical −1/(4y²) endpoint. Bhattacharjee (and earlier Taira–Willems–Wrochna, Carver) are right; CM's arXiv law (Thm 6.1) is right.

# ProlateFals — node N13: Connes–Moscovici prolate counting law vs Bhattacharjee

Orca agent prolatefals, 11 Nov 2026, 15:39–16:45 UTC. Brief: `.scratch/orca/prolatefals-brief.md`.

## Short answer

| Claim | Source | Verdict | Evidence |
|---|---|---|---|
| N(E) = (E/2π)(log(E/2π) − 1) + O(1) | CM21 arXiv 2112.05500 Thm 6.1 (l.904–908) | **true** (CM's proof is a semiclassical area, not a count; proved by TWW26, B26) | b = 0, c = −¼ to 10⁻⁷ |
| … − (log E)/2π + O(1) | **CM22 PNAS Thm 5.1, Eq. 33** | **false** | residual would fall 0.759 from t = 1.7·10³ to 2·10⁵; it falls 0.082, decelerating like (log t)/t |
| … + O(log E) | CM22 Eq. 2; Connes 2402.13082 l.391–393 | true (weak) | — |
| ϑ(t_n) = (n + n₀ − 3/8)π + O(t_n^{−2/3}), n₀ = 0 | B26 Thm B | **true**; bound not sharp: error ≈ −((κ² − ½)/2π)(log t)/t = −25.05(log t)/t | exponent 0.60 → 0.90 → 1 |
| Θ(E) = E log(4E/κ) − E + π/4 + O(E^{−2/3}), κ = 2πλ², all λ, both parities | B26 Thm A | **true**; sharp: Θ − Φ = (κ²/4 − 1/8)(log E)/E + O(1/E) | slopes 39.3529/9.7446/0.4918/157.7892 vs κ²/4 − 1/8 = 39.3534/9.7446/0.4919/157.7887 (λ = √2, 1, ½, 2) |

## 1. Statements

Full side-by-side with definitions and line numbers: [STATEMENTS.md](STATEMENTS.md). Both sides count the same eigenvalues in the same normalisation: N(E) = #{n : 2E_n⁺ ≤ E}, −(E_n⁺)² the negative eigenvalues of the even part of W_sa (λ = √2) on (√2, ∞), regular at √2, ~ sin(2πλx)/x at ∞ (CM22 Prop 4.3(ii) and Thm 5.1 proof; B26 b.377–384). Only t = 2E links them: linear, so no log term can be reconciled away. The two laws contradict each other (difference −(log E)/2π → −∞).

## 2. Computations

### Method (code/prolate.c)
- ODE on (λ, ∞): ((x² − λ²)u′)′ + (4π²λ²x² + E²)u = 0, u analytic at λ (Frobenius start), then Taylor-series steps in x (polynomial coefficients → 5-term recurrence; step ≤ 1.5 rad and ≤ 0.3(x − λ); exact step bookkeeping h = fl(x + h) − x).
- Prüfer angle φ = atan2(k u, √(x² − λ²) u′), k = √(E² + ω²x²), ω = 2πλ, unwrapped from φ(λ⁺) = π/2 (B26's normalisation, eq:prufer). Θ(E) = lim(φ − ωx) evaluated at ωX = max(10⁵, 30E) plus tail ∫(k − κ sinh t) + [a cos 2φ/(2φ̇)] − ∫a²/(2k) (error O(k(X)⁻²)).
- Eigenvalues: even ⇔ Θ ∈ πℤ, odd ⇔ Θ ∈ π/2 + πℤ (B26 Prop 2.5 = CM boundary functionals at ∞), bracketed on a 0.2 grid, Illinois refinement to |ΔΘ| < 10⁻¹².
- Single-threaded; prolate RSS 1.4 MB (ps; top 976 K); ≤ 9 min per run.

### Validation

| Check | Result | Log |
|---|---|---|
| mpmath arbitrary-precision ODE (odefun), raw phase at X = 12, E = 3, 10 | agree to 4.7·10⁻¹⁴ | logs/check_mpmath_out.txt |
| Ramis–Richard-Jung–Thomann certified eigenvalues (C. R. Math. 363 (2025), §4.2.3) μ₋₂ = −39.3832165744668224, μ₋₁₄₈ = −9104.3331714128495040 | ours −39.3832165743689, −9104.3331714132372 (rel. 2.5·10⁻¹², 4.3·10⁻¹⁴) | logs/rrt_check_out.txt |
| CM's own figures (arXiv compareplots1–3.pdf, dots parsed from the PDF streams, axis calibrated on the red zeta dots) | blue dots = our t₁…t₁₀₀, index shift 0, rms 0.17/0.24/0.37 = 0.3 rounding units | logs/cmfigs_out.txt |
| B26's code (theta.c + check_levels.py, t < 50) | t₁ = 12.551192, t₂ = 19.419652, t₃ = 24.703794 vs ours 12.55120975, 19.41966362, 24.70380351 (agree to 2·10⁻⁵, his RK4) | — |
| Independent count: sign changes of raw Wronskians [u,s](X), [u,c](X) (no unwrapping), E-step 0.01 | λ = √2: 200/201 even (last root shifted past E = 200 by finite X), 201/201 odd; λ = 1: 100/100, 101/101 | logs/signchange_*.txt |
| X and step convergence (E = 10⁴, ωX = 10⁵…10⁶, step 1.5/0.75) | Θ stable to 10⁻¹⁰ | — |
| Monotonicity of Θ on the 0.2 grid, interlacing | no violation to E = 1000 (1516 even, 1516 odd eigenvalues), strict even/odd interlacing | logs/eig_*_sqrt2_E1000.txt |

### Counting law (λ = √2, Dirac 2D, t = 2E)
R_n := n − ½ − (t_n/2π)(log(t_n/2π) − 1) (count residual at the midpoint of each step); smoothed R_s(t) = (Θ − Φ)(t/2)/π − ¼.

| t | n | R | CM22 shape −(log t)/2π |
|---|---|---|---|
| 12.55 | 1 | 1.1154 | −0.403 |
| 237.3 | 100 | 0.1322 | −0.870 |
| 1999.5 | 1516 | −0.1779 | −1.210 |
| 2.0·10⁴ | 22491–22495 | −0.23990 | −1.576 |
| 1.0·10⁵ | 138068–138071 | −0.24758 | −1.832 |
| 2.0·10⁵ | 298200–298202 | −0.24870 | −1.943 |

Local slope b = ΔR/Δlog t: −0.0296 (1.7·10³→2·10⁴), −0.0048 (→10⁵), −0.0016 (→2·10⁵); CM22 needs −0.159 forever. Fit R_s = c + b log t + (α log t + β)/t on E ∈ [10³, 10⁵] (logs/analysis3_out.txt):

| λ | b (free) | c | rms, b = 0 | rms, b = −1/2π forced |
|---|---|---|---|---|
| √2 | 2.2·10⁻⁷ ± 6·10⁻⁸ | −0.250003 | 8·10⁻⁸ | 4.7·10⁻² |
| 1 | 1.7·10⁻⁸ ± 5·10⁻⁹ | −0.250000 | 6·10⁻⁹ | 4.7·10⁻² |
| ½ | 5·10⁻¹⁰ ± 1·10⁻⁹ | −0.250000 | 9·10⁻¹⁰ | 4.7·10⁻² |
| 2 | 3·10⁻⁶ ± 8·10⁻⁷ | −0.250036 | 1·10⁻⁶ | 4.7·10⁻² |

Trap: on t ≤ 2000 alone a plain fit gives b = −0.187 ± 0.002 (n ≤ 100: −0.396), i.e. it "confirms" CM22. The decay of R is (κ²/4 − 1/8)(log E)/(πE) (from I(E) and the Bessel edge), which mimics a log over two decades; only t ≳ 10⁴ separates the laws. Carver (below) reports the same trap.

### Θ − Φ, all λ (logs/theta_loggrid_*.txt)

| E | λ = √2 | λ = 1 | λ = ½ | λ = 2 |
|---|---|---|---|---|
| 10 | 4.2355 | 1.7180 | 0.13763 | 8.4014 |
| 10² | 1.3567 | 0.40242 | 0.025163 | 4.3290 |
| 10³ | 0.22658 | 0.062704 | 0.0036490 | 0.79968 |
| 10⁴ | 0.031720 | 0.0085142 | 4.7815·10⁻⁴ | 0.11631 |
| 10⁵ | 0.0040781 | 0.0010758 | 5.9140·10⁻⁵ | 0.015264 |

Split (λ = √2, logs/analysis2_out.txt): Θ − Φ = [I(E) − (E log(4E/κ) − E)] + [η(∞) + π/4]; the first part ≈ (κ²/4)(log E)/E dominates; the Bessel-edge part η(∞) + π/4 = −1.40·10⁻⁵ at E = 10⁵ ≈ −(log E)/(8E) = −1.44·10⁻⁵, which is TWW26's remainder R(E) = −(1/16π)E_TWW^{−1/2} log E_TWW in our variables. So B26's "η(∞) = −π/4 + O(E^{−2/3})" holds with O((log E)/E).

### Location rule (B26 Thm B)
e_n = ϑ(t_n)/π − (n − 3/8), n₀ = 0 (Θ(0) = 0.746 < π; first even eigenvalue t₁ = 12.5512). e₁ = −1.365, e₁₀ = −1.004, e₁₀₀ = −0.382, e₁₅₁₆ = −0.0721, e at t = 2·10⁴, 10⁵, 2·10⁵: −1.0096·10⁻², −2.4225·10⁻³, −1.2981·10⁻³. Local exponent p (|e| ~ t^{−p}): 0.60 (n ∈ [10,50]), 0.72, 0.78, 0.81 (n ≤ 1516), 0.855, 0.887, 0.900 (to 2·10⁵) → 1; e·t/log t = −20.4, −21.0, −21.3 → −(κ² − ½)/(2π) = −25.05 slowly. B26's own check (t < 50) stays at e ≈ −1.4, outside the asymptotic regime.

### Zeros
t_n − γ_n for n ≤ 100: max 2.17 (n = 1: 12.551 vs 14.135). Spectrum = rigid ϑ-lattice; zeros = same lattice + S(T): N(T) − N_ζ(T) = −S(T) + O(1), unbounded (Selberg; Tsang Ω±((log T/log log T)^{1/3})).

## 3. Proof check

**CM22 (loser), proof of Thm 5.1** (PNAS §5 "UV Behavior of Spectrum of Dirac", proof paragraph after Eq. 33; Europe PMC PMC9295779 XML): "we apply Heywood's formula (1.4) of ref. 19 … A Liouville transformation shows that W′⁺_sa is unitarily equivalent to S on (0,∞), Sφ = ∂²φ − qφ, q(y) = −(4π)²cosh²y − ¼(coth²y − 2) [Eq. 34]. A delicate computation which is a refinement of the semiclassical estimate of Proposition 3.2 gives the additional logarithmic contribution −(log E/2π)." **This step is wrong.**
- Heywood (1954) (1.4), as stated by Nursultanov–Rozenblum, Opuscula Math. 37 (2017) pp. 112–113: N(H,(0,λ)) = π⁻¹∫₀^∞[(λ + h)^{1/2} − h^{1/2}]dx + O(1), h = −q, for "q(x) bounded near zero and tending to −∞ at infinity". Here q(y) = −1/(4y²) + O(1) as y → 0⁺: unbounded, critical (LCNO) inverse square. Ref. 20's own singular-endpoint theorem (NR Thm 5.4) needs h″ = O(|h′|^γ), γ < 4/3, which holds for h ~ y^{−α} only if α > 2; α = 2 is excluded. TWW26 (l.233) and B26 (b.67) state the same.
- Evaluated anyway, (1.4) gives J(E)/π, J = ∫₀^∞[√(E² + V) − √V]dy, V = −q. Computed (mpmath quad): J(E) − I(E) + ½ log E = 1.295121, 1.295892, 1.295991 at E = 10³, 10⁴, 10⁵ (converging like (log E)/E to ≈ 1.29600), so J/π = true mean count − (log E)/2π + 0.6625 + o(1), and in t = 2E exactly CM22's −(log t)/2π. The "delicate computation" is right; the formula is not applicable: near y = 0 the regular solution is √y J₀(Ey), whose phase converges (to π/4, Langer: drop the 1/(4y²) from the phase integrand) instead of growing like ½ log E.
- CM21 Thm 6.1's proof (l.909–920) uses the semiclassical area (Prop 4.2) as the count — a heuristic (TWW26 l.233), but the stated law is true.

**B26 (winner)**: read in outline (Lemmas 2.1, 3.1–3.4, Props 2.5–2.6, Thm 3.5); no error found at that level; every conclusion is confirmed numerically above, including continuity, monotonicity, simplicity, the criterion, n₊ = 0, interlacing. Weak points: the O(E^{−2/3}) remainder is far from sharp; the figure argument ("against the paper's own figures", b.67) cannot decide (−(log t)/2π varies by 0.47 over CM's range t ≤ 240, inside O(1)); the attached numerical check stops at t < 50.

## 4. Literature and novelty

| Source | Date | Content | Relation |
|---|---|---|---|
| Carver, Zenodo 21461931 (AI-assisted, author-unverified) | 20 Jul 2026 | certified Prüfer computation, 73 even eigenvalues (checked vs RRT25), N_even(t) = (t/2π)(ln(t/2π) − 1) + ¼ + O(ln t/t) for t ≤ 1265; "the logarithmic term in [33] appears to be erroneous"; Bessel J₀ / Maslov-quarter explanation; ϑ-lattice; notes the log-fit trap | first to flag Eq. 33 (numerically) |
| Taira–Willems–Wrochna, arXiv 2609.32639 | 26 Sep 2026 | rigorous Bohr–Sommerfeld rule (Thm 1.1), Weyl law with O(1) (Cor 1.2), any μ (Rem 1.4); Heywood, NR, Atkinson–Fulton inapplicable (l.233) | rigorous refutation of Eq. 33 (Dirac consequence not drawn) |
| Bhattacharjee, Zenodo 24127214 | 11 Nov 2026 | elementary Prüfer proof, ϑ-lattice (n − 3/8), N − N_ζ = −S + O(1) | correct; main claim not new (cites TWW26, not Carver) |
| Ramis–Richard-Jung–Thomann, C. R. Math. 363 (2025) 1065–1081 | Sep 2025 | analytic spectral determinants, negativity of the non-classical spectrum, certified μ₋₂, μ₋₁₄₈ | no counting law; our benchmark |
| de Andrade Silva, Zenodo 21407212, 20753641 | Jun–Jul 2026 | simplicity + interlacing of the negative spectrum; τ → 0 asymptotics | consistent |
| Connes, arXiv 2402.13082 (l.391–393) | Feb 2024 | cites CM22 Thm 5.1 as + O(log E) | silent weakening, no correction |
| CCM 2310.18423, 2511.22755; 2611.13129 | 2023–2026 | semilocal prolate, zeta spectral triples, uncertainty-operator numerics | no counting-law statement |

Scans: lst.py '(prolate|Connes|Moscovici|Sonin)' over math.NT/CA/CV/SP/PR/FA/CO/AP/math-ph/OA, 2026-10 and 2026-11: no counting-law paper; Zenodo '"prolate"', '"Connes" AND "Moscovici" AND prolate', '"Connes-Moscovici"', '"UV prolate"': the three records above; HAL: 0 hits; repo: no prior entry (R419 lists B26 only).

Novelty: none for us beyond records (largest verified range t ≤ 2·10⁵ vs Carver's 1265; second-order term (κ²/4 − 1/8)(log E)/E, a corollary of TWW26 + the expansion of I). Not a discovery. Suggested lines (orchestrator decides):
- REFUTATIONS.md: Connes–Moscovici, PNAS 119 (2022) e2123174119, Thm 5.1 Eq. 33 (−(log E)/2π): false; proof applies Heywood (1.4) to q ~ −1/(4y²) at y = 0, outside its hypotheses; refuted by TWW26, Carver 21461931, Bhattacharjee 24127214; our numerics b = 2·10⁻⁷ ± 6·10⁻⁸ to t = 2·10⁵.
- NOGO.md (N13): CM prolate UV spectrum is a rigid ϑ-lattice (ϑ(t_n) = (n − 3/8)π + O((log t)/t)); N − N_ζ = −S + O(1) is unbounded, so the operator matches only the smooth zero count, never the zeros.

## Plain words

Connes and Moscovici's published paper says the count of their operator's eigenvalues falls behind the zeta-zero count by a slowly growing amount, (log E)/2π. It does not: we computed the eigenvalues to 12 digits, matched them against certified values, against Connes–Moscovici's own plots and against an independent high-precision solver, and followed them out to the 300,000th eigenvalue. The lag stops growing and settles at exactly a quarter. Their extra term comes from a textbook counting formula used at a singular endpoint where it does not hold. Bhattacharjee is right, but two earlier papers (Taira–Willems–Wrochna in September, Carver in July) had already found this. The corrected picture also shows the limit of the approach: the eigenvalues sit on a rigid lattice and miss the zeros' wobble S(T), so the gap to the actual zeros is unbounded.

RH relevance: N13 (spectral realisation of the zeros; the prolate operator route)

## Files
- [STATEMENTS.md](STATEMENTS.md) — the two laws side by side, definitions, line numbers.
- code/prolate.c — Θ(E), eigenvalues, raw Wronskians (cc -O2 prolate.c -lm). Modes: theta, grid, loggrid, eig, wr.
- code/analysis.py, analysis2.py (uv run --no-project --with mpmath), analysis3.py — fits, windows, split, Heywood integral.
- code/check_mpmath.py, cmfigs.py (needs .scratch PDF parse helper; dots of CM's figures), rrt_check.py, signchange.py.
- Reproduce (in a scratch dir): cc -O2 -o prolate code/prolate.c -lm; ./prolate 1.4142135623730951 eig 0 1000 50000 0.2 > e0_1000.txt (and eig 1 → e1_1000.txt; loggrid 1 100000 26 > lg1.txt; lg_1/lg_0.5/lg_2 likewise; eig 0 10002 -1 0.25 10000 > ew_1e4.txt, 50001 … 50000 > ew_5e4.txt, 100001 … 100000 > ew_1e5.txt; wr 0.01 200 0.01 1500 > wr1.txt); then python3 code/analysis*.py <dir>.
- logs/eig_even_sqrt2_E1000.txt, eig_odd_sqrt2_E1000.txt (j, E_j, Θ − target), eig_even_sqrt2_window_{1e4,5e4,1e5}.txt, eig_{even,odd}_lam1_E100.txt, theta_loggrid_{sqrt2,lam1,lam0.5,lam2}.txt, tn_even_sqrt2_vartheta.txt (n, t_n, e_n), analysis*_out.txt, cmfigs_out.txt, check_mpmath_out.txt, rrt_check_out.txt, signchange_*_out.txt.

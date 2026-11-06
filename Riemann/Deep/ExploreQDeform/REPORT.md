NO IDEA: every q-deformation found moves hyperbolicity (or circle zeros) only away from the RH endpoint — Laguerre's q^{k²} multiplier is the de Bruijn–Newman flow on the log-cylinder (known sector-reducing, Piotrowski–Sorensen–White 2018) and only repackages the N6 record, Rogers–Szegő q-binomials are its inverse, the 2026 q-Laguerre–Pólya classes contain LP, every test passes equally for Davenport–Heilbronn, and Cherednik's DAHA "RH" is collision-freeness near cyclotomic limits that already fails at q = 2…5 for every torus knot T(n, n+1), 3 ≤ n ≤ 7.

# Exploration explore16: q-deformation and quantum groups → RH tree (2026-11-06)

Field: q-orthogonal polynomials (Rogers–Szegő, continuous q-Hermite), q-analogues of the GORZ Hermite limit, q-Laguerre–Pólya classes, Cherednik's RH for DAHA superpolynomials. Item (e) of the brief (Jones zeros, Hall–Littlewood, Hecke traces) was not run: explore11 covered colored Jones polynomials and quantum modularity, and Rogers–Szegő polynomials, which come from Hall–Littlewood functions (Warnaar 2008), are covered by T2.
Read first and not repeated: R46 (q-zeta deformation dead), R34, R56 (multiplier sequences), R338/R341 (Rodriguez-Villegas), R354 (superpolynomials), R378, explore11 (knot invariants), explore15 (screened sector lemma).
Seven scripts, one single-threaded process at a time, nice 19, hard alarm, all under 70 MB and 5 min. Derivations: [DERIVATIONS.md](DERIVATIONS.md). Raw numbers: [outputs.txt](outputs.txt).

## Short answer
1. **The field's q-families all flow the wrong way.** Real-rootedness-preserving q-operations add damping, so they carry hyperbolicity from the RH end (q = 1) down to a trivial end; none carries it up.
   - Laguerre/theta family Σ C(d,j) q^{j²} γ(n+j) X^j (T1): exactly the backward heat flow on the cylinder v = log(−X), i.e. the de Bruijn–Newman flow in the multiplicative variable; its zeros obey v_k' = d + Σ coth((v_k − v_j)/2). Sector shrinking cos Δ' = cos Δ / q is a 2018 theorem (PSW 1802.02641). It gives the "Jensen–Laguerre constant" Λ_JL ∈ [0, 3.098·10⁻²⁶], RH ⟺ Λ_JL = 0 — the N6 record d ≤ 1.6139·10²⁵ in other words (same number to 12 digits).
   - Rogers–Szegő family Σ [d j]_q γ(n+j) X^j (T2): [d j]_q = C(d,j)·e^{+εj²/2}·(geometric), q = e^{−ε}, an anti-damping. It loses real zeros at ε = 2·s*(d,n) (measured ratio 1.005–1.026), and its zeros then sit near a circle. Circle zeros are the log-convex side; Ξ's Taylor coefficients are log-concave.
   - Schur–Szegő family Σ [d j]_q q^{j(j−1)/2} γ(n+j) X^j and Dimitrov–Shapiro's weak q-LP class (2606.17864, June 2026): real-rooted whenever J^{d,n} is (Schur composition; their theorem B_q(LP) ⊆ LP). For q ≤ 1/4 it holds for any positive Turán sequence (Hutchinson), so it is DH-blind.
2. **DAHA "RH" is a collision count, not a mechanism (T3).** For T(n, n+1), n = 3…7, the L-polynomials of Cherednik's superpolynomials keep all zeros on |t| = q^{−1/2} only for q < q_c = 1.185 (= 32/27), 1.058, 1.019, 0.908, 0.850 (his ϖ₀ = 1/q_c < 2 holds), near the two cyclotomic points q → 0 ((1 − u^{2δ+2})/(1 − u²)) and q = 1 (MacMahon q-Catalan). At q_c two circle zeros collide — R46's mechanism. It fails at q = 2…5, where q could be a field size, so it is not a Weil-type statement.
3. **Adversarial tests.**
   - Davenport–Heilbronn: its margins s*(d,n) exceed Ξ's at all 45 computed (d,n) (ratio 1.26–2.78), though J^{d,0}_DH fails near d ≈ 1.93·10⁴. Same Rogers–Szegő law, all q-Jensen and Schur–Szegő checks hyperbolic, Turán ratios > 1 for k ≤ 139.
   - Known wedges: s*/δ(n)² → 1 as n grows (GORZ); Λ_JL's bound is Chasse's sector bound; Dimitrov–Shapiro at q ≤ 1/4 is Hutchinson + Csordas–Norfolk–Varga.
   - One negative for explore15's lead: the Jensen multiplier d!/((d−j)! d^j) is not "Gaussian × multiplier sequence" at the Chasse scale (residual non-real for every d ≥ 3), so the screened sector lemma cannot be had from the exact cylinder flow; it needs a direct finite-degree proof.

## Survey (arXiv API, abstract and e-print pages; listings math.CA, CV, NT, QA for 2026-10 and 2026-11, 1897 entries; one Evidence API call)
- Rogers–Szegő / q-Hermite: RS_n(z) = Σ [n k]_q z^k, H_n(cos θ | q) = e^{inθ} RS_n(e^{−2iθ}); unit-circle zeros for 0 < q < 1 from q-Hermite orthogonality (Szegő 1926; Mazel–Geronimo–Hayes 1990 as cited in the brief, not on arXiv, not checked); up to scaling, RS are the circle orthogonal polynomials of the Jacobi theta weight; real negative zeros for q > 1.
- Laguerre multiplier q^{k²}: Piotrowski–Sorensen–White 1802.02641 (2018), Thm 2.5: e^{−α²k²/2} maps the sector |arg z| ≤ θ into arccos(min{1, e^{α²/2} cos θ}), sharp; corollary: exp(−α²D²/2) is periodic-strip decreasing on exponential polynomials. This is T1's lemma (my derivation duplicated it; Evidence API top passage). Deformed exponential and partial theta function (1709.04357, 2412.02462; Kostov 2605.29991, 2605.29903): the same damping applied to e^x and 1/(1 − x).
- q-analogues of LP: Dimitrov–Shapiro 2606.17864: strong (q-separated zeros) ⊊ LP ⊊ weak (B_q-preimage of LP); q-Jensen criterion; no mention of ζ. Their open question "which q-special functions belong to the classes" has no ζ-relevant answer: Ξ lies in the weak class for every q < 1 iff RH, and for q ≤ 1/4 trivially (Hutchinson + Turán).
- GORZ side: GORZ PNAS 2019; 2108.06775 (log-polynomial sequences), 2301.02492 (Wright circle method), 2511.02628 (Hermite–Jensen limits of q-multinomial coefficient sequences) — Jensen polynomials of q-series coefficients, never a q-deformation of Ξ's.
- Cherednik: 1709.07589 (strong RH conjecture: for ω = 1/q > ϖ'_i the t-zeros have |ξ| = √ω; ϖ₀ < 2 for uncolored algebraic knots; q → 0 limits are roots of unity; numerics only), 2304.02200 (similarity with Ono–Rolen–Sprung zeta-polynomials; "RH similar to the Lee–Yang circle theorem"), 2410.14703 (motivic weak RH justified; Galkin–Stöhr L-functions).
- Circle zeros from L-values: period polynomials (1601.03114, 1607.04699, 2305.03951, 2501.18024, 2610.18065) and zeta-polynomials (1602.00752): dominance plus functional equation, covered by R7.5, R338, R341.
- 2611.03907 (Nov 2026): the Askey–Wilson function is the Laplace transform of a positive measure — positivity gives log-convexity, not zero locations. No paper ties a q-deformation to ζ's zeros.

## Transfers, tests, verdicts
### T1 Laguerre/theta family → N6, N2 (q1_margins.py, q2_decay.py, q4_qlp.py)
- Margin s*(d,n) = sup{s: Σ C(d,j) e^{sj²} γ(n+j) X^j hyperbolic}, binary search with certified sign-change counts (Arb, 1200–1600 bits). RH ⟺ s* ≥ 0 for all (d,n).
- Ξ: s* = δ² exactly at d = 2 (Turán). s*(d,0) = 1.73·10⁻², 5.25·10⁻³, 1.73·10⁻³, 6.76·10⁻⁴ at d = 3, 8, 24, 64, so 2d·s* ≈ 0.083–0.087. s*/δ(n + d/2)² at d = 3 rises 0.57, 0.66, 0.71, 0.76, 0.81 for n = 0, 4, 10, 20, 40 (GORZ limit 1).
- Exact facts: q^{j²} = Gaussian average over rotations of X (the theta kernel); in v = log(−X) it is e^{−τ∂²} on a cylinder; zero ODE above; real neighbours pull complex zeros in (Im coth((x+iβ)/2) < 0, explore15's screening term).
- Λ_JL := −inf s*(d,n): ≥ 0 (anti-damped Ξ has a divergent series); ≤ −log cos Δ₀ = 3.09809·10⁻²⁶ by PSW, with Δ₀ = 2.48921·10⁻¹³ from RH to 3·10¹² and the Lean-checked half-plane. 1/(−2 log cos Δ₀) = 1/sin² Δ₀ = 1.61389709728·10²⁵ (TREE's 1.6141·10²⁵ uses Platt–Trudgian's exact height).
- Verdict: known lemma, RH-equivalent only at q = 1, no new region. Dictionary only: N6 sector lemmas = N2 strip lemmas on the cylinder.
### T2 q-binomial families → N6 (q1, q2, q4)
- Rogers–Szegő: ε_RS/(2s*) = 1.0118, 1.0080, 1.0073, 1.0049, 1.0050 (Ξ: (d,n) = (3,0), (8,0), (24,0), (8,20), (24,40)) and 1.0063–1.0256 for DH; no recovery up to ε = 6 s*; afterwards zeros near a circle (Ξ, d = 12: log|X| ∈ [4.035, 4.059] while |arg(−X)| ∈ [0.24, 2.56] at ε = 4).
- Schur–Szegő: 190/190 certified hyperbolic (q = 0.05…0.95, 5 pairs (d,n), Ξ and DH). Dimitrov–Shapiro q-Jensen J_{n,q}, n ≤ 40, q = 0.25…0.999: 39/39 for each q, Ξ and DH.
- Verdict: dead. RS is anti-damping; the others propagate downward and are strictly weaker than RH for each q < 1.
### T3 DAHA "RH" → N12, N13 (q3_daha.py, q3b_intervals.py)
- At a = 0 the superpolynomial of T(n, n+1) is the q,t-Catalan polynomial; normalisation checked against Cherednik's trefoil zeta and his q → 0 limit S⁰ (DERIVATIONS §4). Functional equation holds for n ≤ 7.
- q_c as in the short answer (T(3,4) exact: 32/27). RH set for T(6,7): (0, 0.905] ∪ [0.975, 1.005]; T(7,8): (0, 0.85] and q = 1. At q_c two circle zeros collide (minimum zero distance 7·10⁻⁶–3·10⁻⁵) and four leave together.
- Sanity test against function fields: Weil's RH holds for smooth curves at every field size; Cherednik's a = 0 polynomials fail at q = 2, 3, 4, 5 for every n ≥ 3. (The Galkin–Stöhr zeta of the singularity z^n = w^{n+1} is a different specialisation, L(Γ, q, t) = Ĥ(q, t, a ↦ −t/q)|_{q↦qt}, which needs the full a-dependence; not computed. Trefoil: 1 − t + q t², on the circle for q ≥ 1/4.)
- Verdict: dead. The mechanism is R46's collision argument near a degenerate cyclotomic point; ζ has no such point on a functional-equation-preserving path.
### T3' Jensen truncation versus Gaussian damping → N6, explore15's lead (q5_jensen_vs_gauss.py)
- J^d(F)(X/d) = F ⊙ Π_{i<j}(1 − i/d): leading Gaussian e^{−j(j−1)/(2d)}, matching Chasse's bound to leading order (Δ⁻² + 1/3 versus Δ⁻² − 1/6).
- Residual Π(1 − i/d)·e^{j(j−1)/(2cd)} tested on (1 + X)^d: non-real zeros for c = 1 at every d ≥ 3 (half of them), c = 2 from d = 6, c = 4 from d = 32; c = 8 real up to d = 40.
- Verdict: no exact "PSW step + multiplier sequence" factorisation at the Chasse scale; a rigorous screened lemma must work with Jensen polynomials directly.
### Side identity (q6_laguerre_mixture.py; q = 1)
J^{d,n}(X) = 8 d!/(4ⁿ (1/2)_n (n+1/2)_d) ∫ Φ(u) u^{2n} L_d^{(n−1/2)}(−u²X/4) du, a positive mixture of dilated Laguerre polynomials (checked to 10⁻⁴⁰). Dilations over (0, ∞) have no common interlacer, so nothing follows; probably classical, not claimed.

## Adversarial tests
- **Davenport–Heilbronn** (same γ code as explore15, functional equation, off-line zero 0.808517 + 85.699348i): margins 3.65·10⁻², 8.37·10⁻³, 2.35·10⁻³, 8.47·10⁻⁴ at d = 3, 8, 24, 64, above Ξ's at all 45 computed (d,n); Λ_JL(DH) > 0 (F_DH is not LP), set by its off-line zeros (the first has angle 0.0072, −log cos = 2.6·10⁻⁵). Every q-family test above passes for DH at computable sizes: the transfers see only zero angles, which is what Chasse's lemma already uses.
- **Known equivalents and wedges**: T1 reproduces the N6 record (Chasse), GORZ (s*/δ² → 1) and de Bruijn's strip lemma (PSW); T2's q ≤ 1/4 result is Hutchinson + Turán; T3 is R46's collision mechanism. Holland's and Savva's wedges (n³ log²(n+2) ≥ 7000 d⁵, n² log(n+2) ≥ 24 d³) are untouched: no q-family gives hyperbolicity of J^{d,n} itself in any new (d,n) region.

## Ranking (RH relevance × chance × gain)
| rank | transfer | node | relevance | chance | gain | result |
|---|---|---|---|---|---|---|
| 1 | T1 Laguerre/theta = cylinder de Bruijn–Newman flow; Λ_JL | N6, N2 | medium | 0 for RH: damping only, RH only at q = 1 | none (= N6 record; lemma known) | dictionary |
| 2 | T3' screened lemma via the exact cylinder flow | N6 (record) | low–medium | 0 by factorisation (q5) | explore15's ×1.6 only by a direct proof | route closed |
| 3 | T2 Rogers–Szegő / Schur–Szegő / weak q-LP | N6 | medium | 0 | — | dead |
| 4 | T3 DAHA "RH" | N12, N13 | low | 0 | — | dead: collision count, fails on [2, 5] |

## Suggested NOGO line
- explore16 (Orca, exploration: q-deformation, Rogers–Szegő/q-Hermite, q-Laguerre–Pólya classes, DAHA superpolynomials → N6, N2, N12, N13) NO IDEA: q-families carry hyperbolicity only away from q = 1 — Laguerre's q^{k²} is the de Bruijn–Newman flow on the log-cylinder (sector-reducing, cos Δ' = cos Δ/q, known: 1802.02641), giving Λ_JL ∈ [0, 3.1·10⁻²⁶] = the N6 record 1.614·10²⁵; Rogers–Szegő q-binomials are its inverse (real zeros lost at ε ≈ 2·s*, ratio 1.005–1.026); Schur–Szegő and the weak q-LP class (2606.17864) contain LP; all DH-blind (DH margins exceed Ξ's to d = 64); the Jensen multiplier is no Gaussian × multiplier product (residual non-real for d ≥ 3); Cherednik's DAHA RH is collision-freeness near cyclotomic q → 0, q = 1 (T(n,n+1): q_c = 0.85–1.19) and already fails at q = 2…5 [ExploreQDeform](Riemann/Deep/ExploreQDeform/REPORT.md)

## Files
- `q1_margins.py` (40 s, 45 MB): γ for Ξ and DH (Arb), Laguerre margins s*(d,n), Rogers–Szegő thresholds, Schur–Szegő checks, near-circle zeros.
- `q2_decay.py` (61 s, 64 MB): s*(d,0) to d = 64, fine ε_RS/(2s*) ratios, DH circle check.
- `q3_daha.py` (265 s, 39 MB) and `q3b_intervals.py` (109 s, 44 MB): q,t-Catalan L-polynomials of T(n, n+1), n ≤ 7, RH thresholds, collisions, q = 2…5.
- `q4_qlp.py` (30 s, 64 MB): Turán ratios to k = 139, Dimitrov–Shapiro q-Jensen polynomials, Λ_JL versus the N6 record.
- `q5_jensen_vs_gauss.py` (1 s): Jensen multiplier versus Gaussian factorisation.
- `q6_laguerre_mixture.py` (13 s, mpmath): numeric check of the Laguerre-mixture identity.
- `DERIVATIONS.md`, `outputs.txt`.
- Commands: `uv run --no-project --with python-flint python3 -u <script>` (q1, q2, q5; q4 adds `--with mpmath`), `--with numpy` for q3, q3b, `--with mpmath` for q6; run inside this directory (q2, q4, q5 import q1, q3b imports q3).

RH relevance: N6 (dictionary only; closes the cylinder-flow route to explore15's screened lemma), N2 (multiplicative de Bruijn–Newman flow, no gain); N12, N13 dead ends.

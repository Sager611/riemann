VERDICT: MECHANISM — PR 137 lowers only the class count R_d at the contact, by proving the marked inverse and two-factor ordinary moments with slot capacity (1−r)/2 + 0.0552429 − 0.0722785·r instead of r + 2S < 1, through a recursive second-Poisson energy induction that is blind to row phases (worth σ = 3.27·10⁻⁴ or λ = 3.53·10⁻⁵ in BeyondB2 units); PR 141 is the same lever at its continuum limit; BARRIER LINES (against PR 137 and PR 141 alike): R310 FALSE, beyondb2 FALSE, R363 FALSE as worded (its dual reading CONFIRMED), TREE N7 sentence FALSE; CEILING θ_∞ = 0.873345856106718576801… = PR 143's degree-9 root (the fifth power of the continuum gain's closed form at the return's KKT contact, derived exactly in poly143.py), binding constraint gate 1 (k ≥ thr, the HIGH child's long-branch margin) in the continuum limit, together with the zero-b-coefficient contact of the return. The anchor induction converges to a linear ODE, every admissible chain stays below its solution, and PR 141 (merged 13 Nov) already sits at that limit, so the mechanism is exhausted 3.611·10⁻⁵ below B_2.

# B2Below: how QuasiRiemannTracker PR 137 goes below B_2, and how far its mechanism reaches (node N7, b2below, 12–14 Nov 2026)

Sources: PR 137 head 9d81a4c0 (521 files, `.scratch/b2below/src`; paths below relative to `Submissions/marked_anchor_induction/Qrh/`); audited scalar model [Hexagon2](../Hexagon2/REPORT.md); levers [BeyondB2](../BeyondB2/REPORT.md). Details: [MECHANISM.md](MECHANISM.md). Session paused 12 Nov 22:36 → 14 Nov 17:23 UTC (machine/API); the tracker state below is from 14 Nov 17:25 UTC.

## Short answer
1. **Mechanism.** One inequality changes: the slot capacity of the two moments that the count window of the PR103 return consumes. The count R_d at the contact falls from 0.6265692817 to 0.6262976997, and nothing else in the scalar system moves.
2. **Barrier lines.** All four claimed that going below B_2 needs cross-row phase information. PR 137 does it with a mean square over rows, so all four are false as worded. The quantitative content of BeyondB2 (B(σ), B_2 − λ) and R310's sensitivity stay correct and predict PR 137's gain.
3. **Reach.** The chain is an explicit Euler scheme for r z′ = z − (1/thr − 2)/6. As steps and pads go to 0 it converges, monotonically, to a linear ODE. The resulting contact θ_∞ is the ceiling of every freedom of PR 137's own lemmas. PR 141 (Bebel, merged 13 Nov 08:28 UTC) proves a rational just above θ_∞, and PR 143 (open) proves θ_∞ itself. Reaching 0.873 by this lever would need 10.7 times the limiting gain.

## 1. Mechanism (MECHANISM.md §1–3)
| | old (PR 103/129) | new (PR 137) |
|---|---|---|
| marked inverse moment, free width r, slots S | r + 2S < 1, 2r + 8S < 3 (Interface/Native.lean:51, :69) | S ≤ markCap r − ν, markCap w = (1−w)/2 + g0 − g1·w, r ∈ [1467/2000, 311/400] (Interface/FiberImproved.lean:55, :62–71) |
| ordinary two-factor moment at doubled width w = 2m | same gates | S ≤ markCap w − ν (FiberImproved.lean:73–82) |
| count window at (δ, x) | R_d = 1 − δ(x + 2(1−x)m_d), m_d = (5/6 − δx)/(5/2 − δ − 2δx) | N/J: faces 1 − δ(r + 2x·markCap r) at the crossing r = 2t/3, long piece 1 − δ + (5/6 − δ)(t − 1) (Return/ImprovedWindowQ.lean:34–44; ImprovedScheme.lean:56–67) |
- (g0, g1) = (552429/10⁷, 7227847/10⁸) is a line under the chain's gains. It touches them at r = 0.753, which is the contact crossing 0.75241. Our exact-rational rerun of the 176 stages (chain.py) finds a minimum margin of 1.6·10⁻⁷ there.
- The rest of the return is PR 103's system. Data.lean gives θ = 1 − c, ℓ = 4c − 1/3, y = h − 2ℓ, h = 81447/100000, which is the Hexagon2 model with b = 2h − 1 − 3ℓ.
  - In that model: native 0.873381967 (= B_2), PR 137 line 0.8733466805 (b = 0.11039).
  - At PR 137's own h, θ = 0.8733468 is feasible with active margin 1.55·10⁻⁷.
  - Every other entry of the full system stays ≥ 2.7·10⁻⁴.
- **Phases.**
  - The proved statement is Σ_{u∈rows}|poly(u)·slots(u)|² ≤ C U^{1+ε}, invariant under row phases.
  - The proof (Inverse/Poisson/*) expands a Gauss-row energy Σ_k|Σ_U c_U χ_U(k)|² and applies a second, exact Poisson summation in k (Grouping.lean).
  - It keeps the column Gauss phases through Phase.lean (χ_A(d)·conj χ_A(h)·χ_A(g)⁴·χ_A(f)⁴ = conj χ_A(y)) and GaussRoot.lean.
  - It then drops all phases: by the triangle inequality over labels and dual frequencies (FDecomp.lean), and by Cauchy–Schwarz on each physical row y of radius 1 − 6u (InBall.lean). The child is again a mean square, paid by the previous anchor.
  - **Kept:** the rows k of the moment's own dual. **Dropped:** the family rows u and the physical rows y.
  - It is a Heath-Brown-type recursive Poisson large-sieve bound for a selected sextic energy. It is not a signed class sum (BeyondB2 S1/S2), and it does not touch the direct pairing S7 (ℓ = 11/3 − 4B is unchanged).
- **Lever map** (levers.py; exact closed form of B(σ)):
  - Dual lever: σ = 3.270·10⁻⁴, against the 0.003518 that 0.873 needs.
  - Direct lever: λ = 3.529·10⁻⁵.
  - R310's F_d multiplier 1 + 7.27·10⁻⁴, times its sensitivity −0.0482, gives −3.5·10⁻⁵, which checks.
  - At the contact, ΔR = −0.31591·G(c) and dB/dG = −0.0408.

## 2. Our barrier lines
| line | verdict | evidence | corrected wording (≤ 2 lines) |
|---|---|---|---|
| R310 (NOGO l.433): "below the B_2 contact, stronger fourth/sixth moments, amplification or detector changes do not help; only cross-row cancellation could move toward 13/15" | **FALSE** | PR 137 goes below B_2 with stronger moments alone, a marked second moment and a two-factor fourth moment, using absolute values over rows throughout; the size matches R310's own F_d sensitivity | R310 N7: stronger marked moments lower B_2 only through the class count R_d (dB/dF_d-multiplier −0.0482); the anchor induction's limit is θ_∞ = 0.8733458561 (PR 141/143, B2Below); 13/15 is out of reach of this lever (window shape keeps R_d ≥ 0.6217) |
| beyondb2 (NOGO l.441): "below Liu's contact B_2 only by a phase-sensitive estimate … every input is phase-blind" | **FALSE** | PR 137's new input is row-phase-blind and monotone under subclasses, yet it moves below B_2; BeyondB2 §1 filed moments as "positive by nature" and held their capacities fixed. B(σ) and B_2 − λ remain correct conditional sensitivities | beyondb2 N7: below B_2 by (i) larger marked-moment capacity → count R_d (row-phase-blind; PR 137 = σ 3.3·10⁻⁴; exhausted at θ_∞, B2Below), (ii) signed class sums B(σ), (iii) a direct-pairing saving B_2 − λ; (ii)–(iii) need phase input that no current input gives |
| R363 (NOGO l.499): "no phase-sensitive large sieve below the contact B_2 (cubic large sieve non-orthogonal, 2607.07911) …" | **FALSE as worded**; dual reading **CONFIRMED** | PR 137's inverse half is a recursive second-Poisson energy bound with exact Gauss-root phases, i.e. a large-sieve-type inequality for the selected sextic family, and it is used below B_2. Still no large sieve gives the dual signed class saving (BeyondB2 §2.3: LS1 loses Z^{0.0333}) | R363 N7: no large sieve gives the dual signed class saving σ below B_2 (signed LS loses Z^{0.0333}; cubic LS non-orthogonal, 2607.07911); a selected recursive-Poisson energy bound does lower the count (PR 137) |
| TREE N7: "below the contact B_2 only by a phase-sensitive estimate (… B_2 − λ; BeyondB2; PR 137 is the first such estimate, an exact second Poisson step …)" | **FALSE** | PR 137 changes neither the direct pairing (Data.lean: ℓ, y unchanged) nor the dual class sum; its Poisson step sits inside a mean-square moment | below B_2: larger marked-moment slot capacity lowers the class count (PR 137 0.8733468; PR 141 0.8733458561 = θ_∞, the exact limit of this anchor induction, B2Below); phase levers σ, λ (BeyondB2) unproved |

## 3. Reach: convergence and ceiling (MECHANISM.md §4–5; chain_limit.py, ceiling.py, perturb.py, ode.py)
- **The stage lemma.** From an anchor (R, z), PR 137 proves the gain Y(r) = (r/R)z + (1/k − 2)(1 − r/R)/6. This requires k > thr with thr = (5wb/4 − 5/6)/(3wb/2 − (1+R)/2 − z); gate 1, H·F0 > 3·10⁻⁵, has F0 = 0 exactly at k = thr (ChainData.lean:75, :83, :102; High/Algebra.lean:59, :74).
- **Admissible HIGH lines give wb ≥ R** (Stg.Hyp: l(R) ≤ z_R and 0 ≤ a1 ≤ 2/5; ChainCheck.lean:52). Since d thr/d wb > 0, the best line passes through (R, z). Then thr = (5R/4 − 5/6)/(R − 1/2 − z), so ∂k/∂R > 0 and smaller steps always win.
- **The limit** is r z′ = z − (14 − 18r − 12z)/(30(3r − 2)) with z(7/9) = 0 (the gain is 0 at thr = 1/2). The ODE is linear, with integrating factor r^{−4/5}(3r − 2)^{−1/5}.

| freedom (PR 137's own lemmas) | tested values | effect |
|---|---|---|
| number and size of stages | 176 × 1/4000 (PR 137); 1/1000 … 1/400000; one step | z(0.7524): 8.436, 8.705, 8.786, 8.794·10⁻⁴ (zero pads), rising monotonically to 8.796·10⁻⁴ (ODE); one step gives 1.6·10⁻⁵ |
| pads (k + 1/5000 ↑10⁻⁶, decrement 4·10⁻⁸, pad 10⁻⁷, ↓10⁻¹⁰, gate 1's 3·10⁻⁵) | PR 137's / 0 | with PR 137's pads, 1/40000 steps are worse (8.350·10⁻⁴); zero pads is the supremum |
| threshold k ≥ thr | k = thr | the gain decreases in k |
| base anchor (311/400, 0), start | 0.775, 311/400, 7/9 | B 0.8733462229, …, 0.8733458617 (7/9 is the last point with thr ≤ 1/2) |
| HIGH-child line | chord vs any line with 0 ≤ a1 ≤ 2/5, l(R) ≤ z | wb ≥ R; the chord attains wb = R (reviewer GAP → CLOSED) |
| count line vs exact gain curve | line, PR 137 curve | 0.8733466805 vs 0.8733466253 |
| window δ ∈ [189/512, 15/32], x ∈ [15/32, 1/2] | PR 137's vs all | 1.1·10⁻¹⁰ (its left edge binds in the limit) |
| h (b) | re-optimized each time | the contact stays at x = 1/2 |
| **all together (ODE)** | | **B* = 0.8733458561067; PR 143's root 0.8733458561067185768; difference 1.1·10⁻¹⁶** |

- **Binding constraint.** Gate 1 (k ≥ thr), with wb ≥ R and the start pinned at thr = 1/2 (R = 7/9). Every other gate has slack ≥ 10⁻⁴ at the contact crossing 0.7523. Iterating adds nothing: the HIGH child already uses every earlier anchor, z(7/9) = 0, and θ enters no gate.
- **Lemma constants (not freedoms; sensitivities).**
  - Poisson radius 1 − 5u instead of 1 − 6u: 0.8733384.
  - F0's 5/6 → 5/6 + 0.001: 0.8733417; 5/4 → 5/4 − 0.001: 0.8733427.
  - Twice the limiting gain: 0.8733173 (levers.py μ = 2; its μ ≥ 5 rows are artefacts, because the crossing leaves the curve's domain [0.7328, 0.7775]).
  - 0.873 needs G(c) ≈ 0.0094, which is 10.7 × z_∞(c).
- **Falsification check.** Since 12 Nov, the tracker has PR 141 (Bebel, merged 13 Nov 08:28, "continuum marked anchors", 0.873345856106718576801) and PR 143 (open, "Exact continuum contact", θ_∞). Nothing is below θ_∞. Our float ceiling of 12 Nov 22:26 (0.8733458652) was 9.1·10⁻⁹ high (δ-grid); ode.py removes that error.

## 4. New barrier (method ceiling) and what a Lean check needs
**Statement.** Use PR 137's stage lemma any finite number of times, from the native anchor (R_0, 0):
- Y(r) = (r/R)z + (1/k − 2)(1 − r/R)/6 with H = r/R ∈ [47/50, 1];
- the hypotheses of High/Algebra.lean Stg.Hyp: l(R) ≤ z_R, 0 ≤ a1 ≤ 2/5, k ≥ thr(R, z, a0, a1) (gate 1);
- any pads ≥ 0, and any lines admitted by ChainLine.

Feed any minorant of the resulting gains to PR 137's improved count window of the PR 103 return, with any h. The resulting half-plane is then Re s > θ with θ ≥ θ_∞ = 0.8733458561067185768…, the root of P(t) = 3888·19⁵(15t−13)⁴(1−t)⁵ − 7(7−6t)⁴(9t−7)⁵ in (13/15, 7/8).

**Proof sketch** (ours; two reviewer rounds, see §5). Write T(ρ, z) = (14 − 18ρ − 12z)/(30(3ρ − 2)). This is (1/thr − 2)/6 at wb = R.
1. **Lines.** Stg.Hyp gives wb ≥ R, and thr increases in wb. With k ≥ thr, each stage gives at most Y(r) ≤ (r/R)z + T(R, z)(1 − r/R).
2. **Monotone in z.** ∂Y/∂z = H − (1 − H)·0.4/(3R − 2) ≥ 0.94 − 0.06·2.02 > 0 on R ≥ 0.7328.
3. **Comparison.** The limit z_∞ solves (z/ρ)′ = −T/ρ². Along it, dT/dρ < 0 on all of (2/3, 7/9]. This uses the second reviewer's repair, because a stage from R ≥ 0.7328 with H = 47/50 can end at r = 0.6888. Put w = (33ρ − 14)/48 − z_∞. Then w′ = (3ρ − 8/5)w/(ρ(3ρ − 2)) with w(7/9) = 35/144 > 0, and dT/dρ = −16w/(25ρ(3ρ − 2)²). Both identities hold to 1.5·10⁻¹⁵ (check_identity.py). So T(R, z_∞(R)) ≤ T(ρ, z_∞(ρ)) on [r, R], which gives E(R, z_∞(R); r) ≤ z_∞(r). With step 2, induction over the stages gives every gain ≤ z_∞, and pads only lower it.
4. **Count.** The balanced count decreases in the gain, and the margin g increases as R_d falls. So B* ≥ B*(z_∞) = θ_∞. Numerically ode.py gives 0.8733458561067; PR 143's root is 0.8733458561067185768.

**Status.** The steps are elementary inequalities, checked numerically at their extremes. The equality with PR 143's polynomial root is numerical (10⁻¹⁶).

**For a Lean check** (no submission; we contact no one):
- (a) An explicit rational super-solution U(ρ) ≥ z_∞ on [0.7328, 7/9] with U(7/9) = 0, and the one-step inequality (r/R)U(R) + T(R, U(R))(1 − r/R) ≤ U(r) for H ∈ [47/50, 1], by posBox/bernCheck boxes (the Poly2 machinery of ImprovedWindowQ.lean). This replaces the ODE.
- (b) Monotonicity lemmas for steps 1, 2 and 4.
- (c) A Bernstein certificate that g < 0 near the contact for B < θ_∞ with the count paid by U. PR 143 already contains the exact-contact algebra.
- This is a theorem about PR 137's scalar system, not about L-functions.

## 5. Checks
- **Reviewer 1** on the ceiling's freedom coverage gave GAP: a HIGH line with slope a1 = 1 gives wb = R − 0.045. That line violates Stg.Hyp's 0 ≤ a1 ≤ 2/5 and ChainCheck's 0 < a1 < 2/5, so the re-review said CLOSED. The reviewer asked us to call the chord ceiling a continuum-limit supremum (done) and noted that it did not re-verify the numbers.
- **Reviewer 2** on the comparison lemma (§4, steps 1–4) gave GAP: stage endpoints can fall below 0.7328. Its own repair, the identity for w, extends the monotonicity to (2/3, 7/9]. We verified both identities and inserted the repair. The reviewer judged lines, z-monotonicity, pads and the start above 7/9 sound, and said the bound follows for every allowed finite chain.
- **Exact-rational rerun of the 176 stages** (chain.py): left 175 = 1467/2000, line minorant margin ≥ 1.6·10⁻⁷, gate minima match the Lean gates.
- **Independent agreement with the tracker:**
  - The scalar model reproduces PR 137's θ (feasible at h = 0.81447) from the PR 103 system.
  - ode.py reproduces PR 143's exact root to 1.1·10⁻¹⁶.
- **Caveats:**
  - The grids in ceiling.py and perturb.py are float; their 9·10⁻⁹ bias is removed only in ode.py.
  - The structural constants of the lemmas (6 in u, 5/6, 5/4, 3/2) were read from the Lean statements, not re-derived.
  - The full-system entries come from BeyondB2's float model.
- **Novelty.**
  - θ_∞ is not ours: Bebel's PR 141/143 state it (13 Nov). PR 143 proves the bound at θ_∞; neither PR says that no schedule beats it.
  - Our parts: the claim that the method's ceiling equals θ_∞, steps 1–4 above, and the corrections in §2.
  - No arXiv or Zenodo scan was done here.


## 6. PR 141 and PR 143 against this ceiling (added 14 Nov 17:45 UTC; PR 141 sources read by raw URL at head 03a1e5ad: TwoSlot/AnchorContinuum/{Curve, IdealCurve, ContactBalance, ContactTangent, AlgebraicContact, Convexity}.lean)
**Same mechanism.**
- PR 141's Curve.lean defines q0(r) = (r⁴(r − 2/3))^{1/5} and P(r) = 11r/16 − 7/24.
- It defines flow(r, z) = A(r)z + b(r), with A = (4/5)/r + (1/5)/(r − 2/3) = (15r − 8)/(r(15r − 10)) and b = (3r − 7/3)/(r(15r − 10)).
- IdealCurve.lean defines ideal(r) = P(r) − (P(7/9)/q0(7/9))·q0(r).
- This is exactly our limit ODE r z′ = z − T(r, z) and our closed form z_∞ = (33r − 14)/48 − (35/144)(9r/7)^{4/5}(9r − 6)^{1/5}.
- PR 141 reaches the ideal curve from below with the strict sub-solutions λ·Z(T, ·), λ < 1, T < 7/9. Their reserve D − flow(W) = (λ − 1)b > 0 pays the finite stages ("zero stage decrement"). The count line is the tangent at the contact width on [3/4, 311/400] (ContactTangent.lean).
- The "averaged physical floor" (AnchorExtension/Library/JointDetector/BoundaryPhysical/PhysicalDirect) keeps the return's other entries feasible. It does not enter the contact value: the derivation below uses only the count and the active margin.

**Binding constraints.**
- **Chain:** gate 1 (k = thr) with wb = R, started at 7/9. Curve.lean forcing_neg says b(r) < 0 exactly for 2/3 < r < 7/9, i.e. thr < 1/2.
- **Return:** the active margin G_N − E(R_d) at x = 1/2, at the amplitude where its b-coefficient vanishes. There R_d = 5/9 + δ/6 (ContactBalance.lean count).
- **Count:** the three-way count balance (short_balance, long_balance).
- The balance reads δ = (5w/4 − 5/6)/(w − 1/2 − ideal(w)) = thr(w, ideal(w)): the contact amplitude equals the HIGH threshold at the crossing width.

**Where the degree-9 polynomial comes from** (poly143.py; exact sympy, 8 s, 16 MB):
- **KKT contact:** δ = (45w − 14)/(6(9w + 1)) (PR 141's density). The balance then forces ideal(w) = −(45w² − 2w − 24)/(2(45w − 14)). The zero margin gives θ = (99w − 14)/(12(9w − 1)), a Möbius map with inverse w = 2(6θ − 7)/(9(12θ − 11)). Equivalently θ = (14 + 25δ)/(18 + 24δ) (theta_density_relation).
- **Fifth power:** the closed form satisfies (P(w) − ideal(w))⁵ = (5/16)⁵·7·w⁴(w − 2/3) (PR 141 ideal_gap_pow_five). Here P − ideal = (95/48)(27w² − 12w − 4)/(45w − 14), which is where 19⁵ comes from.
- **Substitution** gives −48600000·(15θ − 13)·Q₉(θ) = 0. Q₉ has exactly the coefficients of PR 143's root([−274959111914881, …, 487371101862528]), and θ = 13/15 is a spurious factor. PR 143's form 3888·19⁵(15t−13)⁴(1−t)⁵ − 7(7−6t)⁴(9t−7)⁵ equals −Q₉.
- **Values:** θ_∞ = 0.873345856106718576801184796065 (mpmath 40 digits, theta_inf.py; it equals the orchestrator's value), w* = 0.752322053814, δ* = 0.425830047445, R_d = 0.626527230130.

**Ceiling verdict.**
- θ_∞ is the ceiling of the marked-anchor mechanism in the sense of §4. Every admissible finite chain has gain ≤ ideal (our comparison lemma, repaired after the second review), and the count and the margin are monotone in the gain.
- PR 141 (0.873345856106718576801184796, the upper Farey neighbour) and PR 143 (θ_∞ itself, open) sit at the ceiling. Neither claims optimality; the comparison argument is ours.
- **Our own ceiling:**
  - ode.py: 0.8733458561067, 1.1·10⁻¹⁶ from the root.
  - The float value 0.8733458652 from ceiling.py (12 Nov 22:26 UTC, before PR 141 was opened at 13 Nov 04:04 UTC) was 9.1·10⁻⁹ high.
- **Equivalent levers at θ_∞:** λ = B_2 − θ_∞ = 3.61105·10⁻⁵; σ = 3.34630·10⁻⁴ (B(σ) closed form). That is still 1/10.5 of what 0.873 needs.
- **What could still move it:** the stage lemma's constants (Poisson radius 1 − 6u, F0's 5/6, 5/4, 3/2), the count window's shape (native top at 311/400), or another entry of the return. Nothing among PR 137's freedoms does.

**Barrier lines against PR 141.** PR 141 adds a continuum supplier, an averaged floor and an exact return to PR 137's lever, and no phase input. So each verdict of §2 holds for PR 141 unchanged:
- R310: FALSE.
- beyondb2: FALSE.
- R363: FALSE as worded; its dual reading CONFIRMED.
- TREE N7 sentence: FALSE.

The corrected wordings in §2 already cite θ_∞ (PR 141/143) as the exhausted limit.


## 7. Uniform obstruction: referee gap A closed (obstruction.py, 23 exact or interval checks, 0 failed; sympy + mpmath.iv, 5 s, 15 MB)
**Theorem.** Take θ ∈ [17/25, θ_∞) and any h > 0. Pay the return's count by any certificate of the class 𝒞 below, with every improved capacity coming from anchors of PR 137's stage lemma. Then the return's top-edge margin fails at the single point (δ*, x) = (0.4258300474448933…, 1/2), whatever h is. So no θ < θ_∞ is certified. Below 17/25 the count-free 13/15 floor of the same system (R109) applies instead.
- **Witness.** w* = 2(6θ_∞ − 7)/(9(12θ_∞ − 11)) = 0.752322053813543, δ* = (45w* − 14)/(6(9w* + 1)), R* = 5/9 + δ*/6 = 0.6265272301297044.
  - All three are enclosed through PR 141's rational isolation of θ_∞ (width 10⁻²⁴; E1).
  - The KKT gain equals the closed form z_∞(w*) = 0.000885513400196 (E2).

**Proof.**
1. **Margin.**
   - Per J, the Lean top-edge margin is θ − reserve − Ebase(d, C, δ, p) − (3/4)d·R, with d = h + 10⁻¹⁰ and C = min(0.075d, (2y − (2 − physB)d)/physL) (Data.lean; Formulas.lean Ebase, rowMargin).
   - With d = h and C = T/2 it equals the audited margin g exactly (S4). The 10⁻¹⁰ shift changes it by 10⁻¹⁰(9/100 − 3R/4) < 0, the min only lowers C below T/2, and the reserve is > 0. So the Lean margin is < g.
   - At x = 1/2 and R = 5/9 + δ/6: ∂g/∂b = 0 (S1) and ∂g/∂θ = 3 + 4δ (S2). Everywhere, ∂g/∂R = −3h/4 (S3).
   - θ_∞ is the zero of g at (δ*, 1/2, R*) (§6). Hence, for every b, i.e. every h: g(θ, b, δ*, 1/2, R) = (3 + 4δ*)(θ − θ_∞) − (3h/4)(R − R*).
2. **Every certified count at the witness is ≥ R*.** 𝒞 is the count envelope of PR 137 (ImprovedScheme.lean: the minimum over the old count, the native windows and the improved windows) and of PR 141 (ContactBalance.lean short_balance and long_balance).
   - **(a) Native windows** (countN/countJ at the witness):
     - C, weakC7 and capP give 0.6268068944, 0.6268990091 and 0.6270237337, all > R*.
     - capC would give 0.6252219 < R*. But its admission guard 9, upperM·DJ − MN = −0.0544 < 0 there, so no admitted capC window contains the witness.
     - The old count stays > R* for κ = 2θ − 1 ∈ [9/25, 2θ_∞ − 1] (400 interval cells; smallest lower end 0.626612).
   - **(b) Improved windows.** A value V with:
     - V = φ(c), the inverse face at the crossing c;
     - V ≥ φ(2(t − c)), the ordinary face (or, when 2(t − c) > Top, the native plain curve at Top/2);
     - V = L(t) = 1 − δ + (5/6 − δ)(t − 1);
     - where φ(w) = 1 − δ(w + 2x·cap(w)).
     - By §4 and the downward closure of anchor rectangles, cap(w) ≤ (1 − w)/2 + z_∞(w), with only the native capacity below the anchor range.
     - F_env(w) = 1/2 + w/2 + z_∞(w) is convex (S5: q0″/q0 = −(16/225)/(r²(r − 2/3)²) < 0). F_env(363/500) = 0.86756 < F* = F_env(w*) = 0.87705 (P1), and the native F is < F* (P2). So F ≤ F* at every width ≤ w*.
     - Suppose V < R*. Since L is increasing (P3), L(t) < L(t*) gives t < t* = 3w*/2. But φ(c) < R* and φ(2(t − c)) < R* force c > w* and 2(t − c) > w*, so t > 3w*/2. Contradiction.
     - No continuity, balance uniqueness or window shape is used. Every competing balance (short, long, both faces), every window (δ, x) range and every boundary is therefore covered, as is a window that misses the witness, which leaves the native count there.
3. **Conclusion.** For h > 0 (the top conductor edge d = h + 10⁻¹⁰), steps 1–2 give margin ≤ (3 + 4δ*)(θ − θ_∞) − reserve < 0. ∎

**Invariants (ii).**
- **Starts.** Any native anchor (R_0, 0) is admissible.
  - z_∞ ≥ 0 on (2/3, 7/9], because T(ρ, z_∞(ρ)) ≥ 0 ⇔ 9ρ/7 ≥ 3(3ρ − 2) ⇔ ρ ≤ 7/9 (I5).
  - Above 7/9 every stage keeps z ≤ 0, since z ≤ 0 gives Y ≤ z(2.2H − 1.2) ≤ 0 for H ≥ 47/50. A stage straddling 7/9 ends with Y ≤ 0 ≤ z_∞.
- **Threshold denominators** are positive along every admissible chain:
  - 3wb/2 − (1+R)/2 − z ≥ R − 1/2 − z ≥ 363/500 − 1/2 − z_∞(363/500) > 0 (I1);
  - thr increases in wb, because (1 − R)/2 − z ≥ 1/9 − z_∞(363/500) > 0 (I2).
- **Stage endpoints.**
  - Gate 2, (1 − r)/2 + Y(r) < 137/1000 at both endpoints, puts every parent at R ≥ 363/500. There order preservation holds with dE/dz ≥ 3583/4450 (I3, exact).
  - Endpoints r ≥ (47/50)(363/500) = 0.6824 > 2/3 stay in the comparison domain (I4), where T decreases along z_∞ (S6, S7: w > 0 there).
- **Closure.** Every anchor width and every downward-closed width lies in (2/3, 7/9]. The count bound needs only anchors in [w_c, 7/9] with w_c = 0.676973817149 (W1). That range contains PR 137's [363/500, 7/9] and PR 141's [3/4, 7/9].

**Still conditional (exactly):**
- 𝒞 — a count face or detector outside it would be a new mechanism;
- θ < 17/25 — R109's floor, not the witness;
- §4's comparison lemma — two reviewer rounds plus the exact identities S6–S7, read against Stg.Hyp, ChainCheck and H ∈ [47/50, 1] from the Lean source;
- not formalized in Lean.

So this lever is exhausted: no admissible choice of chain, h or window certifies any θ below θ_∞.

## Files
| file | role |
|---|---|
| chain.py | exact-rational replay of the 176-stage chain |
| system.py, infimum.py | scalar model with the improved count |
| chain_limit.py, ceiling.py, perturb.py | continuum limit and the freedom table |
| ode.py | exact ceiling and comparison with PR 143 |
| levers.py | σ/λ equivalents and full-system entries |
| poly143.py, theta_inf.py | exact derivation of PR 143's degree-9 polynomial; θ_∞ and levers to 30 digits |
| obstruction.py | §7: uniform obstruction and chain invariants, 23 exact or interval checks |
| check_identity.py | the reviewer's w-identity |
| run_all.sh | reruns everything into run.log, with peak footprints of 8–21 MB (memcap ≤ 51 MB, single thread) |

Run: `bash Riemann/Deep/B2Below/run_all.sh` (about 2 min; uv with mpmath and numpy).

Sources:
- PR 137 head 9d81a4c0 (local copy)
- [PR 141](https://github.com/GTimG/QuasiRiemannTracker/pull/141), [PR 143](https://github.com/GTimG/QuasiRiemannTracker/pull/143), [PR 103](https://github.com/GTimG/QuasiRiemannTracker/pull/103)
- [Hexagon2](../Hexagon2/REPORT.md) · [BeyondB2](../BeyondB2/REPORT.md) · [R310](../../Tree/R310.lean)

RH relevance: N7

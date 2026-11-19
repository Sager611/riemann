# SGFAudit: the chain behind small-gap-floor, with source lines (sgfaudit, 19 Nov 2026)

Notation as the page: f^(u) = ∫ f(α)e^{−2πiαu}dα, J_θ(f) = f(0) + 2∫_0^θ f(α)α dα, Δ(x) = (1 − |x|)_+, 1_θ = 1[|α| ≤ θ].
TeX downloaded from arxiv.org/e-print (2611.14486 v-current SmallGapsPaper.tex, 574 lines; 1810.08843 LPBandZETAV_17.tex; 2610.28527 Gaps_for_Dirichlet.tex).

## 1. Imported statements (read in the TeX)
- I1 Lamzouri 2611.14486 Thm Main2, l.392–397: "Let 1/2 < θ < 4/7 be fixed and suppose that there exists a constant 0 < λ ≤ 1 and a
  nonnegative even compactly supported function f ∈ C_c²(R) such that f^(0) = 1, f^(u) ≤ 0 for all |u| > λ and
  f(0) + 2∫_0^θ f(α)α dα > 1 (Eq:conditionf). Then D(λ) > 0." Matches the page verbatim, plus the extra restriction 0 < λ ≤ 1 (harmless).
- I2 Fourier normalisation, l.307–309: f^(ξ) = ∫ f(u)e^{−2πiξu}du. Same as the page and as Mathlib's 𝓕 (Lean). So the scale of λ in
  "f^(u) ≤ 0 for |u| > λ" is the page's scale; no 2π rescaling is hidden.
- I3 D(λ), l.183: D(λ) = liminf_T N(T)^{−1}#{0 < γ_n ≤ T : γ_{n+1} − γ_n ≤ 2πλ/log T}, zeros with multiplicity (l.177). Matches the page.
  The barrier is a statement about the hypothesis set of I1, so D's normalisation enters only through I1 itself.
- I4 Lamzouri's own use, l.452–465: θ = 4/7 − 10⁻⁵, f^(ξ) = R(ξ/λ) (Selberg minorant): an instance of I1, no other freedom used.
- I5 CGdL 1810.08843 Thm thmpaircorreltation, l.317–325: "Assuming RH and N*(T) ~ N(T) (l.215–216), N(0.6039, T) ≫ N(T)", where
  N(x, T) (l.219–221) counts PAIRS 0 < γ' − γ ≤ 2πx/log T. Proof, Lemma thmsmallgapszerosofxi l.527–562: for g with ĝ ≥ 0,
  ∫ĝF ≥ ĝ(0) + 2∫_0^1 ĝ(x)x dx + o(1) (l.540–543, F ≥ 0 off [−1, 1], Goldston–Montgomery l.474–477), g ≤ 0 for |x| ≥ λ, g(0) = 1.
  Criterion: ĝ(0) + 2∫_0^1 ĝ x dx > g(0), i.e. J_1(ĝ) > 1 with ĝ ≥ 0: Lamzouri's functional at θ = 1 (identify Lamzouri's f = ĝ).
  Class A_LP, l.420–430: even continuous g ∈ L¹, g(0) = ĝ(0) = 1, ĝ ≥ 0, g eventually ≤ 0. Broader regularity than C_c² (see §3.4).
- I6 DPC 2610.28527 Thm thm:positive-proportion, l.61–69: under RH, ∀ε > 0 ∃c_ε > 0: #{n : T ≤ γ_n ≤ γ_{n+1} ≤ 2T,
  (γ_{n+1} − γ_n)log γ_n/2π < 1/2 + ε} ≥ c_ε T log T, multiplicity counted. Implies D(λ) > 0 for every λ > 1/2 in I3's normalisation:
  apply at T/2; log γ_n ≥ log T − log 2, so gap < 2π(1/2 + ε)/(log T − log 2) ≤ 2πλ/log T for ε = (λ − 1/2)/2 and T large;
  count ≥ c_ε(T/2)log(T/2) ≫ N(T). Page claim "D(λ) > 0 for every λ > 1/2 under RH" holds as written.

## 2. The freedoms of I1 and the LP (lesson 191)
Parameters of Thm Main2 as printed (l.392–394), each with how the dual covers it:
| freedom in I1 | covered by |
|---|---|
| θ ∈ (1/2, 4/7), fixed | certificate at θ = 4/7 exactly (logs/verify_t47exact_0775.log prints theta=4/7) + J_θ ≤ J_{4/7} for f ≥ 0 (J_mono) |
| λ ∈ (0, 1] | certificate at λ = 31/40 exactly (k0·h = 7·31/280) + admissible(λ') ⊂ admissible(λ) for λ' ≤ λ (Admissible.mono) |
| f ≥ 0 on all of R | used twice: ∫f(P − |α|1_θ) ≥ 0 and Σ_{m≠0} f(m/h) ≥ 0; this is the hypothesis that makes the LP an inequality LP |
| supp f any compact set (no supp ⊂ [−θ, θ]) | P − |α|1_θ ≥ 0 on all of R by periodicity (§3 step 6) |
| f ∈ C_c² (not Schwartz) | only continuity, compact support and f^ = O(u⁻²) are used; Schwartz: Lean AdmissibleDecay (Decay.lean l.19–27) |
| f even | f^ real and even; folding 2∫_0^θ fα = ∫_{−θ}^θ f|α| |
| f^(0) = 1, "J > 1" | homogeneous; Lemma 1 gives J ≤ f^(0), so J > 1 = f^(0) is impossible (even J = 1 is the most) |
| "f^(u) ≤ 0 for |u| > λ" (strict) | atoms at |u| = λ: f^(λ) ≤ 0 by continuity of f^ |
| extra factors/shifts in Eq:conditionf | none: l.394 is exactly f(0) + 2∫_0^θ f(α)α dα > 1 |
No freedom of I1 escapes the dual. θ = 4/7 is excluded by Lamzouri but certifying it only covers more.

## 3. Lemma 1 (weak duality), line by line — HOLDS
Data: h > 0, k0 ≥ 1, k0h = λ, η_k = −1 (1 ≤ k < k0), η_k ≥ −1 (k ≥ k0), Σ|η_k| < ∞; σ = Σ_{k≠0} h(1 + η_|k|)δ_{kh};
P(α) = (1 − h) + 2hΣ_{k≥1}η_k cos(2παkh); hypothesis P(α) ≥ α·1[α ≤ θ] on [0, 1/(2h)], 0 ≤ θ ≤ 1/(2h).
1. f ∈ C_c²: |f^(u)| ≤ ‖f''‖₁/(4π²u²), so Σ_k |f^(kh)| < ∞. ✓
2. Poisson hΣ_k f^(kh) = Σ_m f(m/h): the periodisation F(α) = Σ_m f(α + m/h) is continuous (locally finite), its Fourier
   coefficients are h f^(kh) (absolutely summable), so its Fourier series converges uniformly to F, everywhere; evaluate at 0. ✓
   (Normalisation check: ∫f^(xh)e^{−2πimx}dx = h⁻¹f(−m/h) = h⁻¹f(m/h).) Numerically confirmed in §6 (F1).
3. ∫f^dσ = hΣ_{k≠0}f^(kh) + 2hΣ_{k≥1}η_k f^(kh) = Σ_m f(m/h) − hf^(0) + ∫f·2hΣη_k cos(2παkh) (Fubini: Σ|η_k| < ∞, f ∈ L¹)
   = f(0) + Σ_{m≠0}f(m/h) + ∫f(P − 1). ✓ (page's display, DERIVATION §1, reproduced term by term)
4. J_θ(f) − f^(0) = f(0) + ∫f(α)(|α|1_θ − 1)dα (evenness). ✓
5. Difference = Σ_{m≠0}f(m/h) + ∫f(P − |α|1_θ). ✓
6. P ≥ |α|1_θ on R: [0, 1/(2h)] hypothesis; [−1/(2h), 0] evenness; |α| > 1/(2h) ≥ θ: target 0 and P(α) = P(α − m/h) ≥ 0 with
   α − m/h in the base period (period 1/h since cos(2π(α − m/h)kh) = cos(2παkh − 2πmk)). At θ = 1/(2h) the endpoint is shared,
   and P(θ) ≥ θ there by the hypothesis at α = 1/(2h). ✓ So with f ≥ 0: ∫f^dσ ≥ J_θ(f) − f^(0).
7. Weights h(1 + η_k) = 0 for k < k0, ≥ 0 for k ≥ k0, = h for k > K; atoms carrying weight sit at |u| ≥ k0h = λ, where
   f^ ≤ 0 (|u| > λ by I1, |u| = λ by continuity). So ∫f^dσ ≤ 0 and J_θ(f) ≤ f^(0) = 1. □
Direction of every inequality checked; no step uses more than I1 gives. Lean: weak_duality (WeakDuality.lean l.307–312) states
exactly this with general summable η (§5).

### 3.1 Theorem A (discrete sine dual) — HOLDS
a = 1/(1+θ), h = a, k0 = 1 (hard-core condition vacuous), λ = a, η_k = −sinc²(ka) ∈ [−1, 0], Σ|η_k| < ∞.
P(α) = 1 − aΣ_{k∈Z}sinc²(ka)e^{2πiαka} = 1 − Σ_mΔ(α − m/a) (Poisson; sinc² ∈ L¹, (sinc²)^ = Δ). Triangle m is supported on
|α − m/a| < 1; for 0 ≤ α ≤ 1/a − 1 = θ only m = 0 meets α, so P = 1 − (1 − α) = α (θ ≤ 1). Centres are 1/a = 1 + θ ∈ (1, 2] apart:
at most two triangles meet any α, their sum on an overlap is 2 − 1/a ≤ 1, so P ≥ 0. θ ≤ 1/(2a) = (1+θ)/2 ⇔ θ ≤ 1 (equality
at θ = 1: the allowed boundary case of step 6). Lemma 1 ⇒ J_θ ≤ 1 for λ ≤ 1/(1+θ): 7/11 at θ = 4/7, 1/2 at θ = 1. Lean
theoremA' (Triangle.lean l.217) has no hypothesis beyond 0 < θ ≤ 1, λ ≤ 1/(1+θ) (TriangleSeries discharged by triangleSeries).
Numerics (F3): series vs triangle form agree to 8·10⁻⁷ (a = 7/11) and 1.0·10⁻⁶ (a = 1/2), inside the tail bound 2/(π²aN).

### 3.2 Theorem B (lattice duals) — HOLDS: see §4 (certificates) and §5 (Lean statement).
(i) CertB: θ = 4/7, h = 31/280, k0 = 7, λ = 31/40 = 0.775, K = 235. (ii) CertB1: θ = 1, h = 119/600, k0 = 3, λ = 119/200 = 0.595,
K = 202. Both satisfy θ ≤ 1/(2h) (4/7 ≤ 140/31, 1 ≤ 300/119), η_k = −1 for k < k0 and η_k ≥ −1 (asserted in both checkers).

### 3.3 Arithmetic of the claim
0.7851 − 0.775 = 0.0101 ✓. Floors certified in [0.775, 0.7851] (θ < 4/7; upper end = SmallGap's certified kernel, not re-audited
here) and [0.595, 0.6039] (θ = 1; upper end = CGdL l.320). "float LPs ... near 0.784": float only (DERIVATION §3 λ*(X) ≈
0.7839 − 0.207/X; logs/dual2_t47k3_X70.log: float-feasible at 0.78, certification FAILS, verify_t47exact_078.log min −0.0021);
labelled float on the page ✓.

### 3.4 θ = 1 and CGdL's class (gap in the page's coverage, closed here)
The page, DERIVATION and Lean cover C_c² (Admissible) and decaying kernels (AdmissibleDecay: f, f^ = O(|x|⁻²), incl. Schwartz).
CGdL's class A_LP (I5) only asks g even, continuous, L¹, ĝ ≥ 0: not inside either class. Closing argument (mollify):
g_ε = g·e^{−πεx²}, f_ε := ĝ_ε = ĝ ∗ G_ε ≥ 0 (G_ε = ε^{−1/2}e^{−πξ²/ε}), f_ε^ = g_ε ≤ 0 for |x| ≥ λ, g_ε(0) = 1, |g| ≤ g(0) = 1
so Σ_k|g_ε(kh)| < ∞; the periodisation of f_ε is ∫ĝ(β)Σ_mG_ε(α + m/h − β)dβ (Tonelli), continuous, with absolutely summable
Fourier coefficients h g_ε(kh): Poisson holds; steps 3–7 of §3 apply to f_ε, so J_1(f_ε) ≤ 1; ε → 0: f_ε(0) → ĝ(0) (ĝ
continuous, bounded), f_ε → ĝ in L¹. Hence J_1(ĝ) ≤ 1 on all of A_LP, i.e. CGdL's p_f(λ) ≤ 0 (l.513) and 𝒫(f) ≥ 0.595 for
every f ∈ A_LP. So "the method of CGdL can never reach λ ≤ 0.595" holds, with this routine step added (reviewer: §7).
Note: CGdL's conclusion concerns pairs N(λ, T) under RH and N* ~ N, not D(λ); the page's wording ("the same linear program")
is accurate; the LP value does not depend on which conclusion is drawn from it.

## 4. Certificate code audit (verify_dual.py) and independent recomputation
verify_dual.py read line by line:
- Exactness: λ = Fr(str(lam)) (0.775 → 31/40, 0.595 → 119/200), h = λ/k0, θ exact; η_k = round(10¹²·e)/10¹² rational, −1 forced for
  k < k0, clamp to ≥ −1, three asserts (no weighted atom inside (−λ, λ), weights ≥ 0, hard core). ✓
- Outward rounding: arb/acb balls throughout, PASS iff min lower endpoint > 0; 360 bits only narrows balls. ✓
- P = (1 − h) + 2h·Re Q(z), Q = Ση_kz^k, z = e^{2πiαh}: Re z^k = cos(2παkh) ✓; P' = 2h·(−2πh)·Im(zQ'(z)) ✓; M2 = 2hΣ|η_k|(2πkh)² ≥ |P''| ✓.
- Coverage: points {iδ ≤ 1/(2h)} ∪ {θ, 1/(2h)}, sorted; consecutive closed cells cover [0, 1/(2h)] with no gap, scan starts at
  α = 0 ✓ (i = 0) and ends at 1/(2h) ✓. θ is a grid point; a cell with right end ≤ θ checks g = P − α (so P(θ) − θ is checked as
  the closed right end of the last left cell), every other cell starts at ≥ θ and checks g = P. ✓
- Cell bound min(g(a0), g(a0) + g'(a0)δ') − M2δ'²/2 ≤ min over the cell of the Taylor minorant ✓ (weaker than the endpoint
  minimum of the concave quadratic, hence valid).
- Thresholds (lesson 277): the claims are at λ = 31/40 and 119/200 exactly (k0h = λ), lower λ by monotonicity; θ at 4/7 and 1
  exactly. One fragility: without argv[3] the script reads θ from the JSON float 0.5714285714285714 < 4/7, which would not cover
  θ ∈ (0.5714285714285714, 4/7); the logged run (verify_t47exact_0775.log: theta=4/7) and REPORT.md l.57's rerun command pass 4/7. ✓
Independent recomputation (indep_verify.py, own code, 256 bits; log indep_verify.log; memcap peak 15 MB, footprint 14.8 MB):
- Object: the integer lists certB/certB1 parsed from Lean/TheoremB.lean (= the JSON etas after verify_dual's rounding, all 235 and
  202 entries equal: cmp_cert.log), i.e. exactly the Lean props CertB/CertB1.
- Method (different from theirs): each cos/sin from an exact rational argument (arb.sin_cos_pi_fmpq), centred third-order Taylor
  cells g(c) − |g'(c)|r + min(0, g''(c))r²/2 − M3r³/6, M3 = 2hΣ|η_k|(2πkh)³, adaptive bisection from width 1/500, θ a forced breakpoint,
  covered length summed exactly (= 1/(2h)).
- CertB: PASS, 3430 cells (4601 evaluations); smallest centre value of G = 5.7575·10⁻⁴ at α ≈ 1.25508 (theirs: worst cell
  [1.25505, 1.2551]). CertB1: PASS, 3154 cells; smallest centre value 6.8644·10⁻⁴ at α ≈ 1.95394 (theirs: [1.953875, 1.9539]).
- Extreme entries: their worst-cell bounds recomputed with our P: 0.00020445060658898454668201948285665997403740 (CertB) and
  0.00023892856886554410666361113655773630445300 (CertB1), equal to their 360-bit logs in all 40 printed digits.
- Random entries: P at 0, θ, 1/(2h) and 25 random grid points per certificate, Horner (acb) vs direct sums: balls overlap, midpoints
  differ ≤ 3.5·10⁻⁷⁴. Edges (F6): CertB P(0) = 0.008102, P(4/7) − 4/7 = 8.954·10⁻⁴, P(140/31) = 0.01307; CertB1 P(0) = 0.006062,
  P(1) − 1 = 1.1556·10⁻³, P(300/119) = 3.372.
Verdict §4: both certificates hold as stated; true margins are ≈ 2.8× the logged lower bounds.

## 5. Lean statements (read, not rebuilt: brief says no Lean) — FAITHFUL
- Admissible λ f (WeakDuality.lean l.28–34): ContDiff ℝ 2 f, HasCompactSupport f, even, 0 ≤ f, 𝓕f(0) = 1, Re 𝓕f(u) ≤ 0 for λ < |u|.
  Mathlib's 𝓕 is ∫e^{−2πi⟨v,w⟩}f(v)dv = I2's normalisation. "Re" only weakens the hypothesis (𝓕f is real for real even f). = I1's f-class.
- J θ f = f 0 + 2∫_{0..θ} f α·α (l.36) = Eq:conditionf's left side.
- Main2Hyp θ λ f (l.415–416) = 1/2 < θ ∧ θ < 4/7 ∧ 0 < λ ∧ λ ≤ 1 ∧ Admissible λ f ∧ 1 < J θ f: I1's hypotheses verbatim.
- weak_duality (l.307–312): hypotheses exactly Lemma 1's (0 < h, k0h = λ, Summable η, hard core, η ≥ −1, 0 ≤ θ ≤ 1/(2h),
  (if α ≤ θ then α else 0) ≤ P h η α on [0, 1/(2h)]) ⇒ J θ f ≤ 1 for Admissible λ f. Inclusive branch α ≤ θ, as certified.
- CertB (TheoremB.lean l.116–117): ∀α ∈ [0, 140/31], (if α ≤ 4/7 then α else 0) ≤ Pfin (31/280) certB α, with Pfin the finite cosine
  sum of etaOf certB (= L[k−1]/10¹², 0 beyond the list; P_etaOf identifies it with P). CertB1 (l.190–191): same at 119/600, 300/119, θ = 1.
  These Props are hypotheses in Lean; their truth is the Arb certificate (§4, re-verified independently).
- theoremB (l.124–125): CertB → 0 ≤ θ ≤ 4/7 → λ ≤ 31/40 → Admissible λ f → J θ f ≤ 1. not_main2Hyp (l.136–137): CertB → λ ≤ 31/40 → ¬Main2Hyp θ λ f.
  theoremB1 (l.196–197): CertB1 → θ ≤ 1 → λ ≤ 119/200 → J θ f ≤ 1. theoremA' (Triangle.lean l.217–219) and not_main2Hyp_of_le_7_11
  (l.223–224): no extra hypotheses. Decay versions (Decay.lean l.214–219, 343–345) for AdmissibleDecay, incl. Schwartz (of_schwartz).
- Axioms: Check.lean prints axioms for all of these; LEAN.md l.152–175 records [propext, Classical.choice, Quot.sound]; rg over
  ExploreDPP/Lean: no sorry, admit, axiom, native_decide, implemented_by, extern (this audit).
Page's Lean sentence ("Lemma 1, monotonicity and Theorem A machine-checked without hypotheses ..., Theorem B from the named Arb
certificates CertB, CertB1, also for Schwartz kernels; sorry-free, standard axioms") matches the files.

## 6. Falsification (falsify.py, falsify.log; float64 + Arb; peak 21 MB)
- F1 (Poisson with the m ≠ 0 term live): f = B₃(α/5) (cubic B-spline, C², supp [−10, 10] ∋ m/h = 9.032 resp. 5.042, f^ = 5sinc⁴(5u)),
  ∫f^dσ vs f(0) + Σ_{m≠0}f(m/h) + ∫f(P − 1) by direct quadrature: equal to 6.6·10⁻¹⁵ (CertB) and 1.2·10⁻¹⁴ (CertB1). Step 3 ✓.
- F2 (constructible kernel just above the floor): Lamzouri's own f = λ_fR^(λ_f·), λ_f = 0.7892 (l.463–469). J at θ = 4/7 − 10⁻⁵ =
  1.000100147643 (his l.469: 1.000100147…, reproduced, quadrature = closed form); at θ = 4/7: 1.000105535. With CertB the identity
  ∫f^dσ = (J − 1) + ∫f(P − |α|1_θ) holds to 10⁻¹⁶; ∫f^dσ = +0.011781 > 0 comes from the atom at u = 7h = 0.775 (weight 0.7573,
  f^(0.775) = R(0.775/0.7892) = 0.0094 > 0: the kernel is not admissible at λ = 0.775), contributing +0.014243. No contradiction;
  the dual's binding atom sits exactly at λ.
- F3 (corner θ = 1/(2h), Theorem A at a = 1/2 and a = 7/11): series = triangle form within tail bound (§3.1).
- F4 (corner: is f ≥ 0 load-bearing?): yes — step 6 integrates f(P − |α|1_θ) with P > |α| strictly somewhere on [0, θ] for CertB/B1
  (e.g. P(θ) − θ = 8.95·10⁻⁴), so a sign-changing f could make the slack negative. I1 and CGdL both impose f ≥ 0 (l.392 "nonnegative";
  CGdL l.423 "ĝ ≥ 0"), so the hypothesis is available; Theorem A alone (P = |α| exactly on [−θ, θ]) would survive dropping it.
- F5 (smallest failing case of Lemma 1's shape): if θ > 1/(2h) the target's support leaves the base period and step 6 breaks;
  both certificates have θ ≤ 1/(2h) with room (4/7 < 4.516, 1 < 2.521). Not triggered.
No counterexample to any step.

## 7. Reviewer (one, on the most delicate steps; verdict pasted)
Inline brief: Lemma 1 with steps (i)–(vii) and the coverage question (every freedom of Thm Main2), and §3.4 (CGdL class).
- Item 1 (Lemma 1 + coverage): CORRECT within θ' ≤ θ, λ' ≤ λ; "no failing line in (i)–(vii)"; alias terms Σ_{m≠0}f(m/h) ≥ 0 cover
  arbitrary compact support; C_c² gives summability; continuity gives f^(±λ) ≤ 0; monotonicity excludes the strict > 1. Its caveat
  (parameters λ ∈ (0.775, 1] untreated) concerns a stronger claim than the page makes; the page claims only λ ≤ 0.775.
- Item 2 (§3.4): CORRECT, given the full dual hypotheses incl. h ≤ 1/2 at θ = 1 (CertB1: h = 119/600).

## 8. Upper end (float sanity only; SmallGap's own exact certificate not re-audited)
smallgap_sanity.py (mpmath 50 digits) on Riemann/Deep/SmallGap/cert.json: f^(0) = 1.00000000000006, J_θ(f) = 1.0000219 (θ = 399993/700000),
f ≥ 0 sampled on [0, 10], f^ ≤ 0 sampled on [0.7851, 9.79] (max −1.2·10⁻⁷⁴), f^(0.775) = +0.0076 (not admissible at the floor, as §4
requires). Consistent with the floor interval [0.775, 0.7851].

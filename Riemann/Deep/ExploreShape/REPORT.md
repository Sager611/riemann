IDEA: per-gap Wirtinger inequalities of order 4 and 5, fed with BM's own mean values of |F''|², …, |F⁽⁵⁾|², push Bui–Milinovich's large-gap criterion past its Wirtinger ceiling (degree-64 amplifier, float mean values: 3.3235 → 3.3419 for order 4 with an Arb-certified per-gap inequality, 3.3584 for order 5; record 3.324; orders ≤ 2 provably give nothing) (node N8)

# explore53: spectral shape optimisation → N8 large gaps (Bui–Milinovich criterion beyond Wirtinger)

## Transfer
- Tool: the extremal (Dirichlet-spectral) problem for a string/beam on one gap: which quadratic forms
  Σ_j a_j‖f^{(j)}‖² (+ terms that telescope across gaps) are ≥ 0 for every f vanishing at the ends of a gap of length ≤ c.
  This is the eigenvalue-optimisation side of Krein's and Keller–Tadjbakhsh's problems, posed as a cone/SDP.
- Node N8 (large gaps): BM's criterion uses only Wirtinger, Σ_gaps ∫|F|² ≤ (c/π)²∫|F'|², c = κπ/L, for
  F = e^{ivtL}ζ(s)ζ(s+iκπ/L)M(s+iκπ/(2L)). GapCeil/GapBar's ceilings (3.30279; 3.32658, certified failure at 3.335–5)
  fix this Wirtinger step. Here the Wirtinger step is replaced by the best J-th order per-gap inequality.

## Statement (claim tested)
Let u = tL, c_j = lim mean of |d^jF/du^j|² (main terms). BM's Lemma main differentiated in the shifts (as BM do for
j = 1, Q(x) = v + x, l.402) gives c_j(v) = ⟨(v − S)^{2j}⟩, the 2j-th moment of BM's main-term density in the
frequency variable S (GapLambda DERIVATION_BM §2–3): F has an asymptotic spectral measure ν (compact support).
- **Criterion J.** Suppose a_0..a_J, φ_pq (p odd, q even, p, q < J) satisfy, for every gap length g ≤ π and every real
  f ∈ H^J(0, g) with f(0) = f(g) = 0,
  Σ_j a_j‖f^{(j)}‖² + Σ φ_pq [f^{(p)} f^{(q)}]_0^g ≥ 0.                                          (V_J)
  If Σ_j a_j κ^{2j} c_j(v) < 0 for some amplifier and twist v, then Λ > κ.
- Proof sketch: if Λ ≤ κ, consecutive zeros of F are ≤ (1+o(1))κπ/L apart (BM l.180-183). Apply (V_J) (rescaled) on
  each gap; the bracket terms telescope (F's jets are continuous at the zeros) up to O(1) end terms; complex F splits
  into real and imaginary parts (all terms are real-bilinear); end pieces are o(T) by BM's Weyl-bound argument
  (l.194), for every fixed j. Summing gives Σ_j a_j κ^{2j} c_j ≥ −o(1).
- J = 1 is BM's criterion (a = (−1, 1)). Not conditional on RH (same inputs as BM: BBLR ϑ < 1/4, Lemma main).

## Results in one table (float; r = 2, ϑ = 0.2499, midpoint mollifier, degree-24 real amplifier re-optimised for each J)
h_J > 1 means "criterion J proves Λ > κ". Forms from an SDP over polynomial bases with cutting planes, validated on a
degree-29 basis and 300 gap lengths (residual smallest relative eigenvalue: J = 4: −4·10⁻⁹, J = 5: −1.7·10⁻⁷,
J = 6: −3·10⁻⁵, so J = 6 is provisional). Script `scanJ.py`, log `logs/scanJcut_d24.log`, data `logs/scanJ_r2_d24.json`.

| κ | h1 (Wirtinger = BM) | h3 | h4 | h5 | h6 |
|---|---|---|---|---|---|
| 3.32 (no cuts) | 0.994669 | 0.996045 | 1.008756 | 1.020759 | 1.036769 |
| 3.33 | 0.986678 | 0.988107 | 1.001117 | 1.012794 | 1.022013 |
| 3.36 | 0.963350 | 0.964796 | 0.979091 | 0.991470 | 1.000699 |
| 3.39 | 0.940961 | 0.942023 | 0.957971 | 0.971038 | 0.980376 |
| crossing h = 1 | 3.3129 | 3.3147 | 3.3315 | 3.3480 | ≈ 3.355–3.361 |

- J = 6 rerun with 12 cutting rounds (residual −2.7·10⁻⁶) and the Wirtinger amplifier: h6 = 0.99022 at κ = 3.37,
  0.97010 at 3.40, so its crossing is ≈ 3.355, not 3.361 (`logs/scanJ67_d24.log`). J = 7 did not converge
  (residual −0.40: float conditioning of 7th derivatives); no J = 7 number is claimed.

- Same amplifier class, same mean values: orders 4 and 5 add about 0.017 each to κ, order 6 about 0.007–0.013
  (increments +0.0018, +0.0168, +0.0165, ≈ +0.01). The limit J → ∞ is not known.
- GapBar's degree-80 Wirtinger float ceiling is 3.32658 (+0.0137 over degree 24); its certified-failure points
  κ = 3.335, 3.34, 3.35 (every amplifier, r ≤ 3) are failures of the Wirtinger step only: at degree 24, J = 5 already
  succeeds at κ = 3.33 with h5 = 1.0128 (J = 6, provisional, near 3.355).

## Lemma 0 (orders ≤ 2 give nothing): a barrier, proved
For f(0) = f(g) = 0 and j ≤ 2, ‖f^{(j)}‖² = (g/2)Σ_k ω_k^{2j} b_k² exactly (b_k sine coefficients, ω_k = kπ/g:
f'' has sine coefficients −ω_k²b_k because the boundary terms of the two integrations by parts contain f(0), f(g)).
At J = 2 there is no telescoping jet form that survives the reflection x → g − x (|f'|² is reflection-even).
So a form a0‖f‖² + a1‖f'‖² + a2‖f''‖² is valid on all gaps ≤ π iff p(x) = a0 + a1x + a2x² ≥ 0 on x ≥ 1
(x = ω², normalised). By Markov–Lukács duality a symmetric J = 2 criterion succeeds iff (m0, m2, m4) lies outside the
closed moment cone of [1, ∞), i.e. iff m2 < m0 (BM's own criterion) or m0m4 < m2². The latter never happens:
∫|F'|² = −Re∫F''F̄ + o(T) and Cauchy–Schwarz. Hence higher-order Wirtinger with derivatives ≤ 2 never beats BM;
the brief's "sine-series symbol test" is exactly this case. Numerically: SDP J = 2 value = m2 − 1 (logs/sdpJ_d24_k3.32.log).
- For j ≥ 3 the sine-symbol test is false: f = x(g − x) has f''' = 0 but Σ ω_k⁶ b_k² = ∞; the correct per-gap
  problem has two free boundary values f''(0), f''(g) (Lemma 1).


## Lemma 1–2 (exact per-gap validity for J = 3, 4; `crit3.py`, `exact4.py`)
Let f ∈ H⁴(0, g), f(0) = f(g) = 0, b_k its sine coefficients, ω_k = kπ/g, A0 = f''(0), A1 = f''(g), and
d_k = (2/g)((−1)^k A1 − A0). Integrating by parts against sin, cos (only f(0) = f(g) = 0 is used):
f'' ↔ −ω²b (sine), f''' ↔ d_k − ω_k³b_k (cosine; constant term (A1 − A0)/g), f'''' ↔ ω_k⁴b_k − ω_k d_k (sine),
f'(0) = Σω_kb_k, f'(g) = Σ(−1)^kω_kb_k. Hence for
Q(f) = Σ_{j≤4} a_j‖f^{(j)}‖² + φ12[f'f'']_0^g + φ32[f'''f'']_0^g:
Q = (g/2)Σ_k [p(ω_k)b_k² + β(ω_k)b_kd_k + γ(ω_k)d_k²] + (a3 + φ32)(A1 − A0)²/g,
p = Σ_j a_jω^{2j}, β = φ12ω − (2a3 + φ32)ω³ − 2a4ω⁵, γ = a3 + φ32 + a4ω².
Minimising over each b_k (A0, A1 fixed) leaves (2/g)[(A0 + A1)²S_odd + (A1 − A0)²(S_even + (a3 + φ32)/2)],
S_odd/even = Σ_{k odd/even} ρ(ω_k), ρ = (4pγ − ω²q²)/(4p) with β = ωq; the x⁵ and x⁴ terms (x = ω²) of the
numerator cancel identically, so ρ = (n0 + n1x + n2x² + n3x³)/(4p(x)) = O(1/x).
**Q ≥ 0 on every gap g ≤ π ⇔ p > 0 on [1, ∞), S_odd(g) ≥ 0 and S_even(g) + (a3 + φ32)/2 ≥ 0 for all g ∈ (0, π].**
(J = 3: a4 = φ32 = 0; equivalently −α‖f‖² + γ‖f'‖² + a3‖f''' + μf'‖², `crit3.py`.)
- Checks: Lemma 1 (`crit3.valid`) and Lemma 2 (`exact4.check4`) agree on the J = 3 optimum; the SDP's J = 3 point
  (degree-17 basis) is invalid by −4.5·10⁻⁴ in S_odd(π) only (γ_min 0.986056 vs 0.986055): the relaxation is tight.
- All J = 4 forms of the tables satisfy Lemma 2 with margin ≥ 7.6·10⁻¹¹ (binding at g = π, `logs/ex4_test.log`).
  A naive float evaluation of ρ = γ − β²/(4p) cancels catastrophically at small g (it showed −1.3·10⁻⁴); the
  cancelled numerator is used instead.

## Degree-64 amplifier (closer to the record's degree 80)
| κ | h1 | h3 | h4 | h5 |
|---|---|---|---|---|
| 3.34 | 0.987074 | 0.988495 | 1.001440 | 1.013200 |
| 3.36 | 0.971424 | 0.972889 | 0.986656 | 0.998872 |
| crossing h = 1 | 3.3235 | 3.3249 | 3.3419 | 3.3584 |

Degree 80 (the record's amplifier degree; `logs/scanJcut_d80.log`, `logs/scanJ_r2_d80.json`, Lemma 2 check `logs/ex4_d80.log`):

| κ | h1 | h4 | h5 |
|---|---|---|---|
| 3.342 | 0.986129 | 1.000540 | 1.012320 |
| 3.36 | 0.972048 | 0.987238 | 0.999440 |
| crossing h = 1 | 3.32658 (GapBar) | ≈ 3.3427 | ≈ 3.3592 |

The degree-80 J = 4 forms satisfy Lemma 2 in float (margin 3.6·10⁻¹⁰ at g = π).

- J = 4 witness at κ = 3.34 (exactly valid by Lemma 2, float mean values): a = (−0.96172372, 1, 0.015590782,
  −0.022892551, 0.0025713852), φ12 = 0.035361146, φ32 = 0.022892551, twist v = 1.362816, h4 = 1.00144
  (`logs/scanJ_r2_d64.json`). So in floats Λ > 3.34 by BM's own mean values, against the record 3.324 and GapBar's
  certified Wirtinger failure at κ = 3.34 (every amplifier, r ≤ 3).
- J = 5 (validated on polynomial bases only, residual −2·10⁻⁷): κ = 3.3584 at degree 64.

## Arb certificate of the per-gap inequality (J = 4 witness, `cert4.py`, `logs/cert4_k3.34.log`, `logs/cert4_k3.34.json`)
- Form (exact rationals): the κ = 3.34 degree-64 witness with a0 raised by 2·10⁻⁵ (dρ/da0 = β²/(4p²) ≥ 0, so this
  only helps), φ12 raised by 2·10⁻⁶ (telescoping, no global effect; makes the small-gap coefficient n3 > 0, which the
  SDP had left at −8·10⁻¹⁰), φ32 := −a3:
  a = (−5899179108889603/6134092019850000, 1, 14789027537/948575124080, −4037122085/176350902122,
  1982999731/771179571628), φ12 = 3691572375143311/104390384155500000, φ32 = 4037122085/176350902122.
- Checked in Arb (128 bits): (1) p > 0 on [1, ∞) (real roots −5.441, 0.96554; a complex pair); (2) n3 = 1.98·10⁻⁸ > 0,
  largest real root of N is X_N = 1.285668…, so ρ > 0 for x > X_N and every gap g < g* = π/√X_N = 2.7707 is valid
  termwise; (3) g ∈ [g*, π]: 204 adaptive Arb cells, k ≤ 61 (the dropped terms have x_k ≥ 62² > X_N, positive):
  min S_odd ≥ 4.0·10⁻⁶ (at g → π), min S_even ≥ 0.01605. **CERTIFIED: the inequality holds on every gap g ≤ π.**
- With BM's float mean values this rational form gives h4 = 1.00141933 at κ = 3.34 (Wirtinger: h1 = 0.98707445),
  `h4_rational.py`, `logs/h4_rational.log`. Missing for a certified Λ > 3.34: Arb values of c_0, …, c_4.

## Independent check of the higher mean values (`crosscheck.py`, `logs/crosscheck_r1.log`)
GapLambda's float code of BM's printed 9-fold integrand (l.306-312), extended from (v − S)^{≤2} to (v − S)^{≤4}, and
`moments.py` (GapBar's kernel + τ-moment recursion, twist off) agree to all printed digits at r = 1, ϑ = 1/4,
κ = 3.18, degree 5, v = 1.26: h1 = 1.04053590 and ⟨(v − S)^k⟩/⟨1⟩ = 1, 0.030934125, 0.095036118, 0.0088906238,
0.020059772 (k = 0..4) in both. This checks the arithmetic of c_2; the derivative structure itself is the audit item.

## Falsification tests of the per-gap claims (all forms, `falsify.py`, `falsify2.py`)
Every function below has all zero gaps ≤ π, so a valid form must give Σ_j a_j m_{2j}(F) ≥ 0.
- F = sin(x)H(x), H any trigonometric polynomial of period 2πM (M = 4, 8; frequencies ≤ 12), complex (twists
  included): minimum over H = p(1) > 0 for every form (0.0063 J3, 0.0335 J4, 0.035 J5), attained at H = 1.
- F = sin x · tanh(sin x/ε) (smooth |sin|: 81% of its spectrum at frequency 0): positive for ε ∈ [0.02, 3].
- F = H(x) sin φ(x), φ' = 1 + β(1 + cos(x/R)) ≥ 1 (chirps; gaps π/φ'), H free (81 modes): minima 0.025–13.8 > 0.
- Wirtinger itself: minimum 0 (± 2·10⁻¹⁴) on sin x, as it must be.

## Why orders ≥ 3 help and orders ≤ 2 do not (mechanism)
- In u = tL units F has an asymptotic spectral measure ν (BM's main-term density pushed to S); its Hankel matrices
  are positive definite up to order 10 (smallest eigenvalue 6.8·10⁻⁸ at degree 24), as they must be.
  Standardised central moments (degree 24, κ = 3.32): skewness −5·10⁻⁴, kurtosis 2.192, m6 6.81, m8 27.4, m10 ≈ 134:
  a unimodal, compactly supported, nearly symmetric bump, far from the two-atom measure (δ_{±π/c}) that makes
  Wirtinger sharp.
- Wirtinger uses only Var ν. A gap only forces its Dirichlet frequencies to be ≥ π/c, so any symbol p ≥ 0 there is
  admissible for the sine modes; a degree-2J symbol can stay just above 0 on [1, 3] (scaled), where ν has its mass,
  and dig below 0 on [0, 1). The price is the curvature of f at the zeros: the free boundary values f''(0), f''(g)
  enter through S_odd, S_even (Lemma 2), and these are the binding constraints (at g = π).
- At order ≤ 2 the curvature does not enter and the moment-cone condition is automatic (Lemma 0), which is why the
  brief's "symbol test" sees nothing.
- Structural limit of the whole local class: per-gap quadratic data cannot see a sign or phase change of F from one
  gap to the next (jets enter the telescoping terms quadratically). So no local inequality of any order separates F
  from functions glued from independent pieces on the gaps; the J → ∞ limit of the method is that extremal problem,
  not computed here.

## Adversarial test (iv)
- The new step (V_J) holds for every function whose zero gaps are ≤ π (tested above on three families); it says
  nothing about ζ by itself. All ζ-specific content enters through ν, i.e. through BM's Lemma main and BBLR's
  twisted fourth moment (Euler product, d_r(h)-mollifier, ϑ < 1/4), exactly as in BM. For the Davenport–Heilbronn
  function those main terms are different (no Euler product; its own mean values would have to be computed), so the
  inequalities used here distinguish ζ from DH only as BM's criterion does. ν is a genuine positive measure
  (Hankel PD), so no positive-definite-distribution loophole is used.

## Status: what is proved, what is float
- Proved here: Lemma 0 (barrier for orders ≤ 2), Lemma 1–2 (exact per-gap criteria for J = 3, 4), the reduction
  (V_J) ⇒ criterion J (same summation and end-piece argument as BM l.164-202), and (Arb) the validity of the
  rational J = 4 witness on every gap ≤ π.
- Float only: the mean values c_j, j ≥ 2 (BM's Lemma main differentiated j times in the shifts: Q(x) = (v + x)^j
  in place of BM's v + x; c_j(v) = ⟨(v − S)^{2j}⟩, checked against GapBar at j ≤ 1 to 10⁻⁸: d = 8, κ = 3.25 gives
  h1 = 1.01560266 = optmid's value); all h_J values; J = 5, 6 validity (polynomial bases).
- Audit item before any certificate: that BM's main term c(α, β) depends on the shifts only through the exponentials
  y^{−…}, T^{−…} (so each L⁻¹∂ contributes one factor v − S); this is what BM use for j = 1 (l.402) and GapLambda
  re-derived (DERIVATION_BM §2).
- Certificate path (not done, ~1 agent-day): (a) Arb τ-moments (τ1 + τ2)^m, m ≤ 8, in GapBar's `certmid.py`
  Γ-arrays, giving c_0..c_4 for the degree-64 or 80 midpoint amplifier; (b) rationalise the J = 4 form with a 10⁻⁶
  margin and verify Lemma 2 in Arb on g ∈ (0, π]: DONE for κ = 3.34 (section above). Remaining: (a) only.
  Expected certified record ≈ Λ > 3.34 (J = 4), ≈ 3.355 if a J = 5 exact criterion (4 boundary values) is added.

## Novelty (NOVELTY.md)
- In-field "generalized Wirtinger" inequalities (Hall 2002 JNT 93 and 97, Beesack, Saker, BHSJ 2511.14415) are
  first-derivative inequalities in higher powers (L⁴, L^{2k}, ∫H(y'/y)y^{2k}); BHSJ (Nov 2025, l.431) call
  Wirtinger-type inequalities in f, f' "the only known method". No use of second or higher derivatives per gap
  found (arXiv API queries, 2026-10/11 listings, Zenodo, our repo). Classical higher-order Wirtinger inequalities
  (Böttcher–Widom 0412269) assume all derivatives vanish at the ends, which F does not satisfy.

## Answer to the brief
- (i) Tool: per-gap Dirichlet extremal problem for higher-order forms (beam/string spectral optimisation, solved as
  an SDP + exact sine–cosine reduction). Node N8, large gaps.
- (ii) Smallest case: the J = 2 sine-symbol test gives nothing, provably (Lemma 0), and the symbol test is false for
  j ≥ 3 (Lemma 1–2 replace it). Bounded scale: BM's mean values from GapBar's own kernel (degree 8–64 amplifiers).
- (iii) New float bound Λ > 3.34 (J = 4) and ≈ 3.358 (J = 5) unconditionally, beyond the record 3.324 and beyond
  GapBar's certified Wirtinger failures at 3.335, 3.34, 3.35. Per-gap step Arb-certified (J = 4, κ = 3.34); mean
  values c_j, j ≥ 2, float only, so no certified Λ yet. Not conditional on RH.
- (iv) The per-gap step is function-generic; ζ enters only through BM's mean values, as in BM.
- (v) No prior use found.
- Next step (one Orca agent): Arb certificate along the path in "Status" (target Λ > 3.34), plus an exact J = 5
  criterion (four boundary values f''(0), f''(g), f''''(0), f''''(g)) and r ∈ {1, 3}, degree 80.

## Suggested NOGO line
- explore53 (Orca, exploration: spectral shape optimization → N8 large gaps) NOT A BARRIER: GapCeil/GapBar's ceilings fix BM's Wirtinger step; per-gap inequalities with derivatives ≤ 2 provably add nothing (Markov–Lukács + Cauchy–Schwarz), but order-4/5 inequalities with BM's mean values of |F''|², …, |F⁽⁵⁾|² reach κ ≈ 3.342/3.358 in floats (degree 64; Wirtinger 3.3235), order-4 per-gap inequality Arb-certified at κ = 3.34, mean values c_j (j ≥ 2) float and to audit [ExploreShape](Riemann/Deep/ExploreShape/REPORT.md)

RH relevance: N8 (large gaps between critical zeros, the Λ → ∞ direction). It reopens the Λ record route that
GapCeil/GapBar had closed at ≈ 3.33 and shows their ceilings cover only the first-derivative step.

## Files
- `moments.py` (BM main-term spectral moments ⟨S^k⟩, k ≤ 12; reuses GapBar gcore/gcore_r read-only), `joint.py`
  (central forms, joint amplifier/form loop), `scanJ.py` (κ scans), `sdpJ.py`, `sdp3.py`, `sdpcut.py` (SDP
  relaxation + cutting planes), `validate.py` (large-basis check), `crit3.py` (Lemma 1), `exact4.py` (Lemma 2),
  `falsify.py`, `falsify2.py`, `cert4.py` (Arb), `h4_rational.py`, `crosscheck.py`, drivers `ex4_test.py`, `ex4_d80.py`, `val_scan.py`, `val_J5.py`, `cut_test.py`.
- Logs: `logs/scanJcut_d24.log`, `logs/scanJcut_d64.log`, `logs/scanJ67_d24.log`, `logs/scanJ_r2_d{24,64}.json`,
  `logs/ex4_test.log`, `logs/falsify*_*.log`, `logs/mom_*.json`, `logs/crit3_r2_d24_k3.32.json`. Peak memory
  ≤ 263 MB (memcap), single-threaded.

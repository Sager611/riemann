RECORD: Λ_f > 1.744 (= 218/125) unconditionally for every primitive L-function of a holomorphic or Maass cusp form on GL(2) over Q, by an Arb-certified order-6 per-gap inequality fed with the printed mixed moments of all derivatives (arXiv 1410.7765, Thm mixed); was √3 = 1.7320508 (Barrett–McDonald–Miller–Ryan–Turnage-Butterbaugh–Winsor 2015). No certified gain elsewhere (reasons per problem below).

## Statement
- Λ_f := lim sup (γ̃_{n+1} − γ̃_n)/(π/log(√q γ̃_n)) over critical-line zero ordinates of L(s, f). Then Λ_f > 218/125.
- Proof: ExploreShape's Criterion J with J = 6. g(t) = e^{it log T}L(½+it, f); 1410.7765 Thm mixed gives
  ∫_T^{2T}|g^{(j)}|² = 2c_f T L^{2j+1}/(2j+1)·(1 + o(1)), L = log T, for every fixed j (the spectral measure is uniform on
  [−1, 1]). An order-6 form valid on every gap ≤ π (exact rationals a_0..a_6, φ_pq; Arb, Lemma 4 / cert5.py) has
  H := −Σ_j a_j Λ^{2j}/(2j+1) = 7.7047·10⁻⁴ > 0 at Λ = 218/125 (exact). If Λ_f ≤ 218/125, fix A > 218/125 with H(A) > 0
  (H is continuous); gaps of g in [T, 2T] are then eventually ≤ Aπ/L; summing the form over the gaps (Re g and Im g,
  telescoping jets, end pieces o(TL) after scaling, by subconvexity) gives −2c_f H(A) T L ≥ −o(TL), a contradiction.
- Float optimum of the method (uniform measure): order 3 no gain, order 4 ≈ 1.740, order 5 ≈ 1.7441, order 6 ≈ 1.7452.
- By-products (float only): Bui–Hall's Z^{(k)} criterion re-computed for k ≤ 10 (higher orders add ≤ 0.002; their
  printed k = 2 amplifier class does not reach 1.606).

## Targets (Wirtinger-type large-gap theorems; "line" = TeX line of the Wirtinger step)
| problem | paper | old record | Wirtinger step | input for j ≥ 2? | float gain (orders 3–6) | certified |
|---|---|---|---|---|---|---|
| GL(2) L-functions | Barrett et al. 1410.7765 | √3 (uncond.) | l.309-325, applied l.343 | yes, printed for all μ, ν (Thm mixed l.239) | 1.7321 → 1.7452 (J = 6) | **1.744** |
| zeros of Z^{(k)}, k ≤ 10 | Bui–Hall 2304.05181 | Λ^{(1)} > 1.9; C(k) table | Thm Wi l.233, l.255-264 | Prop main1 stated j ≤ 1; its proof (l.420-433, R^j) gives every j | k = 1: 1.9009 → 1.9026; k = 2: +0.0004 in h | none (no Arb mean values) |
| ζ, unamplified Z (Hall 2005) / BM / Bredberg | 1101.3197, 1410.3635 | 3.367 (ours) | — | — | excluded (brief: no ζ certificate) | — |
| ζ, Hall's L⁴ route | Bui–Hall–Subirà Jorge 2511.14415 | 2.64 (P = 1) | l.505-507 (Hall's L⁴ form) | quartic form; our tool is quadratic: would need ∫|(Z²A²)^{(j)}|² (amplified 4th moments with shifts) | not computed; ζ (excluded); far below 3.367 | — |
| Dedekind ζ_K, K quadratic | Bui–Heap–Turnage-Butterbaugh 1410.3888 | 2.866 | l.161-195 | shifted twisted moments (l.212), so every j by Cauchy as in BM; not checked line by line | not computed (time) | — |
| some Dirichlet L(s,χ), χ even (q-aspect) | Bredberg 1003.2290 | 3.54 | l.190, l.227 | not checked (time) | not computed | — |
| zeros of Z' via Opial/Sobolev | Saker–Steuding 1109.3855 (RH) | — | Opial-type | PDF only; superseded by Bui–Hall | — | — |
| conditional on conjectured moments | Saker 1001.0494, 1002.2695 | 6.14, 9.64 | Opial/Boyd | conditional on moment conjectures: not a target | — | — |
| not Wirtinger | 2610.28527, 1112.6038, 1305.4068, 1704.05834 | — | Montgomery–Odlyzko / heuristic | not a target | — | — |

## Method
- Read each target's TeX (Wirtinger line, mean-value lemma, which derivatives it gives). Where every ∫|F^{(j)}|² follows,
  compute the spectral moments m_{2j} and run ExploreShape's order-J SDP with cutting planes (sdpJ.py, sdpcut.py,
  validate.py, copied unchanged), then certify the per-gap form with ExploreShape's cert5.py/exact5.py/exact6.py
  (copied unchanged; Lemma 4 audited in GapCertAudit REPORT2) and the criterion H > 0 in exact rationals (gl2crit.py).
- GL(2): moments exact, m_{2j} = Λ^{2j}/(2j+1); form from gl2form.py at Λ = 1.744 with SDP margin 10⁻⁵ (the margin-free
  forms at 1.744, 1.745 fail Arb validity by ~10⁻⁹ at g = π, then at g ≈ 1.11: logs/cert6_gl2_*.log).
- Z^{(k)}: zk.py evaluates Bui–Hall's 9-fold main term (re-derived; working notes below) by 7-D Gauss–Legendre; zjoint.py
  runs the joint amplifier/form loop (linear P = BH's class, twist v, midpoint η = λ/4).

## Evidence (logs in Riemann/Deep/GapFamily/logs/)
- GL(2): uniform_134.log, uniform_56.log (bisection: J=1 [1.73125, 1.73516] ∋ √3; J=3 same; J=4 [1.73906, 1.74297];
  J=5 [1.74414, 1.74465]; J=6 [1.74516, 1.74566]); gl2form.log; gl2_J6_1.744m.json (float form);
  cert6_gl2_1.744m.log + cert5_gl2J6_1.744m.json (exact rational form; (1) p(1) = 0.0146, no real root ≥ 1; (2) det C =
  5.8·10⁻¹¹; (4) X* = 2.2195, g* = 2.1087, 110 cells, min det S_odd ≥ 2.64·10⁻⁷, min det(S_even + C/2) ≥ 6.0·10⁻⁷;
  CERTIFIED); gl2crit_1.744m.log (H = 7.704685931·10⁻⁴ exact rational, CRITERION HOLDS). Peaks: SDP 163 MB, cert 22 MB.
- Z^{(k)}: bhtable_d1.log, zjoint_k1_l1.92_J4.log, zjoint_k1_l1.90_J356.log, zjoint_k2_l1.606_J6.log, bhopt_k2_d4.log,
  zjoint.jsonl. Peaks ≤ 176 MB.

## Checks
- Clean rerun (.scratch/gapfamily/clean, fresh copies of cert5.py, exact5.py, exact6.py, gl2crit.py and the float form):
  CERTIFIED, H identical, output JSON byte-identical (cmp).
- Sanity: the J = 1 bisection brackets √3 (their theorem); order 3 gives nothing, as for ζ (ExploreShape: order-3 gain tiny).
- Moment identity: (−1)^{μ+ν}2^{μ+ν+1}/(μ+ν+1) = ∫_0^2(−s)^{μ+ν}ds, so the binomial expansion of g^{(j)} = e^{itL}(iL + d/dt)^j L
  recombines to ∫_0^2(1−s)^{2j}ds = 2/(2j+1). Λ_f's normaliser log(√q t) = log T + O_f(1) on [T, 2T].
- zk.py (Z^{(k)}): k = 0 reproduces ExploreShape's independent crosscheck to all printed digits (h1 = 1.04053590 and
  four moments); the Q-weights reproduce Hall's ∫Z'⁴/∫Z⁴ ~ L⁴/560 exactly (z4check.py, mpmath).

## Review
- reviewer (20 Nov 2026; inputs, reduction, moment formula and certificate pasted inline): CORRECT, ROUTINE ("improved
  constant record"). Requested clarification adopted in Statement: apply the form at a fixed A > 218/125 with H(A) > 0,
  not at gaps π(1 + o(1)). It verified the analytic implication, not the Arb run (that is the audited cert5.py pipeline).

## Novelty
- python3 .scratch/lst.py (math.NT 2026-10/11, 'Wirtinger|large gap|gaps between (consecutive )?zeros|derivatives of
  Hardy|Z^{(k)}'): 6 titles (2611.14486, 2610.14127, 2610.22927, 2610.26038, 2610.28527, 2610.31744): none on GL(2) or
  higher-order inequalities. Zenodo '"large gaps"': 190 records, none on L-functions.
- arXiv API: 'gaps zeros "cusp form"', '"large gaps" "modular form"', '"Wirtinger" "L-function"', 'gaps zeros GL(2)',
  'gaps zeros automorphic': only 1410.7765 (and LeClair 1704.05834, heuristic). No later GL(2) large-gap bound found.
- Our files: the per-gap tool is ours (discoveries/large-gaps-higher-order.md); its transfer to GL(2) is new here.
- Zenodo '"Wirtinger"' (52 records) and '"Hardy Z-function" zeros': none on zero gaps; listing titles matching
  'GL\(2\)|cusp form' (7, Oct–Nov 2026): none on zero gaps.

## Z^{(k)} (Bui–Hall 2304.05181): no certified gain, and a falsification note
- Input: available by BH's own proof for every j (R^j in l.420-433), although Prop main1 is stated for j ≤ 1.
- Float gain tiny: k = 1, linear P: crossings ≈ 1.9009 (Wirtinger, re-optimised) → ≈ 1.9026 (order 6); k = 2: +0.0004 in h.
  Mechanism: ν_k is nearly two-atomic (kurtosis 1.294 for k = 1, 1.145 for k = 2; ζ: 2.19; two atoms: 1), where Wirtinger is sharp.
- Not certified: c_{k,j} would need an Arb version of zk.py's 7-D integral (≈ GapCeil/GapBar's Γ-array machinery with the
  Q-weights); not attempted (gain ≤ 0.002).
- BH's printed numerics do not reproduce (float, two node sets, code validated as in Checks): k = 1, P = 1−2.5x, v = 0.2,
  η = 0.5 gives h_1 = 0.99315 (printed > 1.02); re-optimised linear P gives 1.00083, so Λ^{(1)} > 1.9 stands. k = 2: max over
  all linear P and all v of h_2(1.606) = 0.99581 (η = 0.5), 0.99588 (midpoint); printed > 1.003. A degree-4 P gives 1.0175,
  so C(2) = 1.606 is true by BH's method but not by its printed amplifier class. k = 3..10: C(k) hold (h = 1.0050–1.0132).
  BH's printed c_{k,j} formula (l.345) also omits the η-dependence of the exponent. Reviewer not run on this (float only).

## Caveat
- The order-6 per-gap inequality is ours and unrefereed (audited once, GapCertAudit); the GL(2) record inherits that.
- Thm mixed of 1410.7765 is used as printed (all μ, ν, uniform error); not re-derived here.
- The gain is a constant record (√3 → 1.744, +0.7%); method ceiling ≈ 1.7452 at order 6; J → ∞ open (gapjinf).
- Dedekind (1410.3888) and q-aspect Dirichlet (1003.2290) targets were identified but not computed.

RH relevance: N8 (large gaps between critical zeros; the per-gap tool transfers to every L-function whose mean values of all
derivatives are known; for uniform spectral measures the gain is √3 → 1.7452 at order 6).

## Files
- Riemann/Deep/GapFamily/: REPORT.md; uniform.py, gl2form.py, gl2crit.py (GL(2)); zk.py, bhcheck.py, bhopt.py, bhtable.py,
  zjoint.py, k0check.py (Z^{(k)}); copied unchanged from ExploreShape: sdpJ.py, sdpcut.py, validate.py, cert5.py, exact5.py,
  exact6.py; logs/ (all runs above). .scratch/gapfamily/: STATE.md, z4check.py, meta.py, arxq*.py.
- Reproduce: `python3 gl2form.py 6 1.744m 1e-5` (uv --with numpy --with scipy --with cvxpy), then
  `JJ=6 SRC=logs/gl2_J6_1.744m.json python3 cert5.py gl2J6_1.744m 0 0 0 61` (uv --with python-flint --with numpy; 10 s,
  22 MB), then `python3 gl2crit.py gl2J6_1.744m 218/125`.

## Ready lines
- RECORDS.md: "GL(2) large gaps: Λ_f > 1.744 = 218/125 unconditionally for every primitive holomorphic or Maass cusp-form
  L-function over Q (gapfamily, 20 Nov 2026: order-6 per-gap inequality + printed mixed moments of all derivatives, 1410.7765
  Thm mixed; Arb + exact H = 7.7·10⁻⁴; reviewer CORRECT, ROUTINE); was √3 (Barrett et al. 2015) [GapFamily](Riemann/Deep/GapFamily/REPORT.md)"
- REFUTATIONS.md candidate (float; two implementations; reviewer CONFIRMED): "Bui–Hall 2304.05181 Thm 2, k = 2: no linear P and no v gives h_2(1.606) > 1
  (max 0.9959, printed > 1.003); a degree-4 P repairs it (1.0175). k = 1 printed parameters give 0.993 (printed > 1.02) [GapFamily]"

## Bui–Hall discrepancy: second implementation (orchestrator follow-up, 20 Nov 2026)
- Outcome: CONFIRMED (reviewer). BH's printed numerical witnesses are not given by their own main-term formula, read either
  way; the theorem statements Λ^{(1)} > 1.9 and C(2) = 1.606 survive with other amplifiers (float). arXiv has only v1
  (o47_ver.py: v1, 2023-04-11); journal version Mathematika 69(3) (2023) 780–794, DOI 10.1112/mtk.12206 (CC BY-NC-ND), full
  text not compared (Wiley returned HTML to curl).
- bh2.py, independent of zk.py: BH's general formula l.418-433 (operators Q⁴R^{2j}, Q(x) = (½+x)^k l.371, R(x) = v+1+x l.396)
  at the point of l.434 (α1 = β1 = 0, α2 = −β2 = πiλ/L, α3 = −β3 = 2πiη/L), integrand typed term by term in BH's 9 variables
  (x, x1..x4, t1..t4; nested simplex map, tensor Gauss–Legendre in all 9, no substitution, no closed forms); h = c_{k,0}/(λ²c_{k,1})
  (l.441-444); P(1−x−x1−x2)P(1−x−x3−x4) with P = 1 − a·(arg), ϑ1 = 0.2499. Mode 'literal' = exponent printed at l.345 (no η).
- Results (logs/bh2.log; nodes (6,8,10) and (7,9,12) agree to 2·10⁻⁷; zk.py agrees to 2·10⁻⁷):

| case | printed | general (η in exponent) | literal l.345 | max over linear P, v: general / literal |
|---|---|---|---|---|
| k=1, λ=1.9, P=1−2.5x, v=0.2, η=0.5 | > 1.02 | 0.9931524 | 0.9864064 | 1.0008115 / 0.9944235 |
| k=2, λ=1.606, P=1−2x, v=0.2, η=0.5 | > 1.003 | 0.9939420 | 0.9931845 | 0.9958050 / 0.9947325 |

- Reading R = v + x (BM's twist) instead: h ≈ 0.205 (k=1), 0.270 (k=2). Scaling P cancels in h.
- So the literal l.345 formula (η dropped) does not reproduce the printed values either: the discrepancy is not that slip.
- k = 2, every η (zk.py, logs/bh_eta_k2.log; max over linear P and v): η = 0, 0.15, 0.3, 0.4015 (= λ/4), 0.55, 0.8:
  0.994732, 0.995412, 0.995800, 0.995878, 0.995713, 0.994749 (symmetric about λ/4, as the β ↔ α−β symmetry requires):
  no linear P reaches h = 1 at 1.606 for any η. Degree-4 P (general, η = λ/4): 1.017538.
- Reviewer (verdict-only, formulas and numbers pasted inline): CONFIRMED; "these findings invalidate the stated numerical
  certificates, not the theorem statements"; recovery is float (no rigorous error bounds for the replacement witnesses).

## Run (from Riemann/Deep/GapFamily; env OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1; each under bash ../../../.scratch/memcap.sh 290 900)
- GL(2) form: uv run --no-project --with numpy --with scipy --with cvxpy python3 gl2form.py 6 1.744m 1e-5   (writes logs/gl2_J6_1.744m.json)
- GL(2) per-gap: JJ=6 SRC=logs/gl2_J6_1.744m.json uv run --no-project --with python-flint --with numpy python3 cert5.py gl2J6_1.744m 0 0 0 61
- GL(2) criterion: uv run --no-project python3 gl2crit.py gl2J6_1.744m 218/125   (H = 7.704685931e-4, CRITERION HOLDS)
- BH, second implementation: uv run --no-project --with numpy --with scipy python3 bh2.py 1 1.9 0.5 general 7,9,12; same with
  "1 1.9 0.5 literal 6,8,10", "2 1.606 0.5 general 7,10,12", "2 1.606 0.5 literal 6,9,11"
- BH, zk.py: uv run --no-project --with numpy python3 bhcheck.py 1 1.9 2.5 0.2 0.5 8,8,6,8,6,12,12;
  uv run --no-project --with numpy --with scipy python3 bhopt.py 2 1.606 <eta> 1 8,8,6,8,6,14,14 (eta = 0 0.15 0.3 0.4015 0.55 0.8);
  same with "2 1.606 mid 4 8,8,6,8,6,14,14"; bhtable.py 1,2,3,4,5,6,7,8,9,10 1

# Appendix: working notes (derivations, as written during the run)

## Working notes (derivation; appended as I go)
### Z^{(k)} (Bui–Hall 2304.05181): the input is BM's Lemma main with Q-weights
- BH l.237-264: f_k = e^{2ivθ(t)} Z^{(k)}(t) Z^{(k)}(t+πλ/L) A(1/2+it+2πiη/L), A = Σ_{n≤y} P(log(y/n)/log y) n^{-s},
  y = T^ϑ, ϑ = 0.2499, r = 1. Λ^{(k)} ≤ λ ⇒ zeros of f_k in [T,2T] are ≤ (1+o(1))πλ/L apart; Wirtinger (Thm Wi, l.233)
  gives h_k = c_{k,0}/(λ² c_{k,1}) ≤ 1 + o(1) (l.262). Records (l.59-130, 466-487): Λ^{(1)} > 1.9 (P = 1−2.5x, v = 0.2,
  η = 0.5, h_1 > 1.02); C(2..10) = 1.606, 1.451, 1.365, 1.306, 1.265, 1.231(table: 1.232), 1.205, 1.184, 1.167.
- Mean values: Prop main1 (l.338-356) is stated for j = 0, 1, but its proof (l.420-433) applies R(·)^j R(·)^j to
  c(α,β) for general j (operator identities l.410-415 exact for every j; Cauchy step l.403 as in BM). f_k^{(j)} =
  (iL)^j QQ R^j(…) + O(L^{j-1}·…) (θ' = L/2 + O(1), θ^{(m)} ≪ T^{1-m}). So c_{k,j}, j ≤ J, follow by the same proof.
- Printed c_{k,j} (l.345) omits the η-dependence (exponent has no η). Re-derived from the general formula (l.425):
  exponent at α1=β1=0, α2=−β2=iπλ/L, α3=−β3=2πiη/L is i[c(X2−X4) − 2πη(X1+X2−X3−X4) − c t3 D3 + c t4 D4],
  c = πλ, X = ϑx, X_i = ϑx_i, τ1 = (1−X1−X3)t1, τ2 = (1−X2−X4)t2, σ = τ1−τ2, D3 = X3−X4+σ, D4 = X1−X2+σ.
- BH's linear factors are exactly D4, D3 and ab dt1dt2 = dτ1dτ2; Q-arguments are p3+t3D3, q3−t3D3, p4+t4D4, q4−t4D4
  with p3 = 1/2−X3−τ1, q3 = 1/2−X4−τ2 (q3−p3 = D3), p4 = 1/2−X1−τ1, q4 = 1/2−X2−τ2 (q4−p4 = D4). Substituting u = t3D3:
  c_{k,j} ∝ ∫ P(w1)P(w2) e^{iE0} G3 G4 (v+1−S)^{2j}, E0 = c(X2−X4) − 2πη(X1+X2−X3−X4), S = X+ΣX_i+τ1+τ2,
  G3 = ∫_0^{D3}(p3+u)^k(q3−u)^k e^{−icu}du, G4 = ∫_0^{D4}(p4+u)^k(q4−u)^k e^{icu}du; domain X+X1+X2 ≤ ϑ,
  X+X3+X4 ≤ ϑ, τ1 ∈ [0,1−X1−X3], τ2 ∈ [0,1−X2−X4], w1 = 1−(X+X1+X2)/ϑ, w2 = 1−(X+X3+X4)/ϑ.
- Check k = 0: e^{iE0}G3G4 (η=0) = c^{-2}[e^{ic(X1−X3)}+e^{ic(X2−X4)}−e^{ic(X1−X4)}e^{icσ}−e^{ic(X2−X3)}e^{−icσ}] = GapCeil's
  ρ_j bracket (GapCeil DERIVATION §1). At X_i = 0: G4 = conj G3, kernel |G3|² ≥ 0.
- Midpoint twist (GapBar §3): 2πη = c/2, i.e. η = λ/4 (BH use η = 0.5 fixed). Symmetry β ↔ α−β holds for Z^{(k)} too
  (Z^{(k)}(−t) = (−1)^k Z^{(k)}(t), real A-coefficients).
- So the order-J criterion (ExploreShape Criterion J) applies verbatim with ν_k = law of S under this kernel and v+1 in
  place of v: Σ_j a_j λ^{2j} c_{k,j} < 0 ⇒ Λ^{(k)} > λ.
### Float code and its checks (zk.py, k0check.py, bhopt.py; logs/)
- zk.py: 7-D tensor Gauss–Legendre over (X, r1, θ1, r2, θ2, τ1, τ2) (X_1 = r1θ1, X_2 = r1(1−θ1), …), inner GL for G3, G4.
- k = 0 check against ExploreShape crosscheck (r = 1, ϑ = 1/4, κ = 3.18, degree 5, v = 1.26, unshifted mollifier):
  h1 = 1.04053590 and |⟨(v−S)^n⟩/⟨1⟩| = 1, 0.03093412489, 0.09503611834, 0.008890623706, 0.02005977184 (n ≤ 4),
  identical to all printed digits (nodes 10,10,7,10,7,14,14; 8-node run differs by 4·10⁻⁵).
- Q-weight check (independent, mpmath, .scratch/gapfamily/z4check.py): ϑ → 0, λ = 0 limit gives
  ∫Z'^4/∫Z^4 ~ L⁴·6∫∫G_1² = L⁴/560 exactly, Hall's ∫Z⁴ ~ TL⁴/(2π²), ∫Z'⁴ ~ TL⁸/(1120π²).
- Bui–Hall's printed parameters do NOT reproduce: k = 1, λ = 1.9, P = 1−2.5x, v = 0.2, η = 0.5 gives
  h_1 = 0.99315 (nodes 6 and 8 agree to 2·10⁻⁷), not > 1.02. Re-optimised linear P (a = 2.376), v = 0.237, η = λ/4:
  h_1 = 1.00083 > 1, so Λ^{(1)} > 1.9 still holds by their method. k = 2, λ = 1.606, η = 0.5: max over all linear P and
  all v is h_2 = 0.99581 < 1 (printed h_2(1.606) > 1.003): C(2) = 1.606 is not reached by the printed amplifier class.
### Bui–Hall Table 2 recomputed (bhtable.py, logs/bhtable_d1.log; max over all linear P and all v, Wirtinger)
- h at λ = C(k), η = 0.5 | η = λ/4: k=1 1.000812 | 1.000832; k=2 0.995805 | 0.995878; k=3 1.005739; k=4 1.004997;
  k=5 1.006592; k=6 1.005355; k=7 1.009117; k=8 1.010967; k=9 1.012480; k=10 1.013147 (midpoint values ≤ 5·10⁻⁵ higher).
  So C(k) holds for k ≠ 2 with re-optimised linear P (BH's printed P, v give smaller h: k=1 0.99315). k = 2 fails in the
  linear class even with the order-6 tool (h6 = 0.996278, logs/zjoint_k2_l1.606_J6.log); a degree-4 P gives Wirtinger
  h = 1.017538 at 1.606 (logs/bhopt_k2_d4.log), so C(2) = 1.606 is true by BH's method, not by its printed amplifier.
### Higher-order gain for Z^{(k)} is tiny (zjoint.py)
- k = 1, linear P, η = λ/4, λ = 1.90: h1 1.000832, h3 1.000839, h5 1.002052, h6 1.002520; λ = 1.92: h1 0.981828,
  h4 0.982346 (dh/dλ ≈ −0.95). Crossings: J=1 ≈ 1.9009, J=6 ≈ 1.9026 (+0.0018). k = 2 at 1.606: +0.0004 in h.
- Mechanism: ν_k is near the two-atom Wirtinger extremal: kurtosis 1.294 (k=1), 1.145 (k=2) against 2.19 for ζ
  (ExploreShape); a two-atom measure has kurtosis 1 and leaves no room for any per-gap form (Wirtinger is sharp).
### GL(2) (Barrett–McDonald–Miller–Ryan–Turnage-Butterbaugh–Winsor, arXiv 1410.7765): input available as printed
- Thm mixed (l.239-246): ∫_T^{2T} L^{(μ)}(½+it,f)L^{(ν)}(½−it,f̄)dt = (−1)^{μ+ν}2^{μ+ν+1}/(μ+ν+1)·c_f T(log T)^{μ+ν+1}
  + O(μ!ν!T(log T)^{μ+ν}) for ALL μ, ν ≥ 0. Wirtinger step l.309-325 (Bredberg's Cor. 1), applied to
  g = e^{iρt log T}L(½+it,f) (l.343-352): Λ_f ≥ √3 (Thm gaps, l.164; ρ = 1, l.384-389).
- Spectral form: (−1)^{μ+ν}2^{μ+ν+1}/(μ+ν+1) = ∫_0^2(−s)^{μ+ν}ds, so ∫|g^{(j)}|² = c_f T L^{2j+1}∫_0^2(ρ−s)^{2j}ds (1+O(1/L)):
  ν = uniform on [−1, 1] at ρ = 1; in gap units (gaps ≤ Λπ/L, mean spacing π/L) m_{2j} = Λ^{2j}/(2j+1).
- End pieces: subconvexity |g^{(j)}| ≪ T^{1/2−δ+ε} (l.336-341) makes them and the end jets o(TL^{2j+1}).
- So ExploreShape's Criterion J applies verbatim: a valid (V_J) with Σ_j a_j Λ^{2j}/(2j+1) < 0 gives Λ_f > Λ.

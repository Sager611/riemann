OPEN: B3 (short-window simple zeros, barrier θ > 1/2) — keeping SL's exp(iφ) exact instead of the Taylor step (11) keeps E1's dual-term bound below 1/2 (|𝒯(n)|/(H X0^{−1/2}) = 0.29–0.59, flat to θ = 0.349, T ≤ 1e8, complex shifts), so with the frozen detector at ν = 0.121 the method gives a positive proportion for every θ > 0.49467 (h(0.4947) = 7.1e−5; one-piece optimum ≈ 0.4940), reviewer GAP (uniform dual estimate unwritten), not certified. ROBUST: B1 (every freedom of Bernard's th:mm2 is inside the certified class; Littlewood + Jensen sharp); B2 (the (7/11)ℤ sine DPP shares all bandwidth-4/7 data; deeper density estimates lose (log T)^c); B4 (TEN in its named class, λ < 1 loses 0.685(1 − λ); AF for bandwidth one; Prop O analytic).

(ceilaudit, Orca, 20 Nov 2026, 02:33–04:25 UTC)

# CeilAudit — falsification of our own N8 method ceilings

Scope: B1 GL(2) one-piece Levinson–Conrey ceiling; B2 small-gap 7/11; B3 short-window θ > 1/2; B4 simple-zero ceilings
(TEN 153/227, AF bandwidth-one 0.6818287, Prop O). Lessons 191/309: a ceiling must cover every freedom of the method's own
lemmas and name the step it fixes. All runs single-threaded under bash .scratch/memcap.sh 290 <s>; peaks in logs/.

## Notes (written as I go)
- B3 step: SL (11) Taylor-expands exp(iφ), φ = −2πξ(e^{2y/H} − 1 − 2y/H) ≍ T/H² (E1Audit FULLPROOF l.9–10); for θ < 1/2, T/H² > 1
  and the expansion diverges. The exact profile needs no expansion: with u = H(e^{2y/H} − 1),
  A(ξ) = ξ^{−β}π^{−1/2}∫ G(u)e^{iwu}du, G(u) = e^{−y²+(1−2β)y/H}e^{−iT(u/H − log(1+u/H))}·dy/du, w = (T − 2πξ)/H (FFT in u).
  Heuristic: the profile is a chirp of amplitude ≍ H/√T on |T − 2πξ| ≲ T/H (L¹ mass ≍ √T, not H), but the Bessel phase
  meets the chirp's curvature (2π)²/T at a stationary point for n ≲ N0, so |𝒯(n)| ≍ H X0^{−1/2} as for θ > 1/2 (tested below: confirmed).
- B3 float test 1 (b3_profile.py; logs/b3_profile_real*.log; peak 134 MB): exact-profile dual term 𝒯(n), shifts α = 0.01, β = 0.02.
  FFT profile vs brute y-quadrature 2·10⁻⁶–1.6·10⁻⁵, vs Taylor J = 8 (τ = T/H² = 0.11) 1.7·10⁻⁵. max_{n≤N0}|𝒯(n)|/(H X0^{−1/2}):
  | T | H | θ | τ | L¹(A)/H | pq = 6 | pq = 30 | decay ≥ 8N0 |
  |---|---|---|---|---|---|---|---|
  | 1e6 | 1000 | 0.500 | 1 | 0.42 | 0.290 | 0.467 | ≤ 5.5e−6 |
  | 1e6 | 500 | 0.450 | 4 | 0.80 | 0.425 | 0.526 | ≤ 4.0e−6 |
  | 1e6 | 250 | 0.400 | 16 | 1.60 | 0.371 | 0.542 | ≤ 5.4e−6 |
  | 1e6 | 125 | 0.349 | 64 | 3.19 | 0.375 | 0.546 | ≤ 4.2e−6 |
  | 1e8 | 10⁴ | 0.500 | 1 | 0.42 | 0.419 | 0.324 | ≤ 3.1e−6 |
  | 1e8 | 2500 | 0.425 | 16 | 1.60 | 0.536 | 0.340 | ≤ 5.6e−6 |
  | 1e8 | 1250 | 0.387 | 64 | 3.19 | 0.542 | 0.372 | ≤ 5.1e−6 |
  The profile's L¹ mass grows like 0.4√τ (chirp), but the per-term size stays ≍ H X0^{−1/2} (ratio 0.29–0.55, flat in τ ≤ 64 and T),
  with E1's decay beyond N0: E1's Step 4–5 bound, hence the dual ≪ T^{1/2+ε}H^{−1/2}y², survives below θ = 1/2.
  Mechanism (proof route): swap the y- and ξ-integrals; the inner ξ-integral is Weber's ∫₀^∞ J₀(b√ξ)e^{−pξ}dξ = p^{−1}e^{−b²/(4p)}
  (p = 2πi(e^{2y/H} − 1), b = 4π√(n/(pq))), a pure chirp in y; stationary point y_s² = πn/(2N0) gives |𝒯(n)| ≍ H(pq/(nT))^{1/4}e^{−O(n/N0)}
  with no condition on T/H² (the cutoff χ(ξ/T) and the Y/K branches need the usual bookkeeping).
- B3 float test 2 (b3_onset.py, logs/b3_onset.log): frozen SL/EXP-010 detector (as nu12.py; check ν = 3/25 → R = 177/20,
  κ = 0.0844779079697, as E1), R ∈ [7.80, 9.00]: onset θ_c(ν) = max{(4ν+1)/3, θ_det(κ(ν))}, θ_det solving k*(θ) = κ.
  ν = 0.1200: κ = 0.0844779 → 0.494993; ν = 0.1205: 0.0847592 → 0.494814; ν = 0.1210 (R = 44/5): κ = 0.0850388, θ_det = 0.494637,
  moment side 0.494667 → θ_c = 0.494667; ν = 0.1215: 0.495333. So with the moment valid for ν < (3θ−1)/4 below 1/2,
  E1's chain gives liminf S/N ≥ h(θ; 0.0850388) > 0 for every θ > 0.49467 (float; a detector re-optimised at ν ≈ 0.121 lowers it further).
- Asymptotics behind test 1 (y-stationary phase after the swap, d = 0): |𝒯(n)| ≈ 2^{−1/2}·H X0^{−1/2}·e^{−πn/(2N0)} for every θ < 1
  (y_s² = πn/(2N0), y-phase curvature ≍ T/H ≫ 1); matches the logs' decay: e^{−π} ≈ 0.04 (≥ 2N0), e^{−2π} ≈ 1.9e−3 (≥ 4N0), e^{−4π} ≈ 3.5e−6 (≥ 8N0).
- B2 test (b2_density.py, logs/b2_density.log): 7/11-DPP form factor = min(|α|, 1) exactly on |α| ≤ 4/7, = 4/7 on [4/7, 1].
  Non-uniform density freedom (Conrey 8/7 near the line, Ingham, Huxley, Guth–Maynard deeper, constants 1): max_a e^{qθa}n(a) for
  θ = 0.58, q = 2.01: 10^{0.33}, 10^{0.97}, 10^{1.78}, 10^{2.88}, 10^{3.85} at L = 10², 10³, 10⁴, 10⁶, 10⁸ (peak depth a ≈ 33 … 387 ≍ log L·const);
  θ = 0.56, q = 2.02 (qθ < 8/7): 10^0 throughout. Grows like a power of L whenever qθ > 8/7: the log losses of the deeper estimates block θ > 4/7.
- B1 test (b1_complex.py, logs/b1_complex.log, 46 MB): independent low-degree Ritz (deg P ≤ 7, deg Q ≤ 9), ξ = 25/256, R ∈ {9.5, 9.87, 10.3}:
  jointly complex (P, Q) BFGS from 12 starts returns the real optimum exactly (gain 0); basis too small (κ = 0.028 < 0.070), so weak;
  gl2audit's high-degree tests (cplx.py second variation, cplxglobal.py) stand.
- B3 float test 3 (b3_profile.py with complex shifts α = 0.05 + 0.1i, β = −0.03 + 0.07i; logs/b3_profile_cplx.log, 115 MB):
  T = 1e6, θ = 0.500/0.450/0.400/0.349, pq = 6, 30: ratios 0.313/0.505, 0.460/0.568, 0.400/0.585, 0.405/0.590; decay ≥ 8N0 ≤ 5.5e−6. Same picture.
- B3 float test 4 (b3_ceilnu.py, logs/b3_ceilnu.log, 81 MB): detector side at its one-piece optimum (gl2ceil's exact Euler–Lagrange
  evaluation imported read-only, Q(1) = β − 1 ≈ −1.5·10⁻⁴ at the optimum, so SL's class Q(1) = 0 is essentially attained):
  κ*(3/25) = 0.0861107 (frozen detector 0.0844779); ν = 0.120: κ* = 0.0861107, detector θ > 0.493955, moment θ > 0.493333 → 0.493955;
  ν = 0.121: 0.0868303 → 0.494667. One-piece detector limit θ_c* ≈ 0.4940 (frozen detector 0.49467).
- h at ν = 121/1000, κ = 0.0850388 (logs/b3_h.log): h(0.49467) = 3.7·10⁻⁵, h(0.4947) = 7.07·10⁻⁵, h(0.495) = 4.05·10⁻⁴, h(1/2) = 5.93·10⁻³.
- Reviewer (reviewer_b3.txt; formulas pasted inline): VERDICT GAP — C1 (exact-profile dual estimate) lacks a uniform, summable proof;
  "its stationary-phase mechanism is plausible"; exact joint saddle y_σ = σH·arsinh√(πn/(2pqT)), ξ_σ = (T/2π)e^{−2y_σ/H};
  Hessian det D²(Φ/K) = −2(e^{4s} + e^{2s}) (u = 2πξ/T, K = T/H, s = y_σ/H): nondegenerate, independent of T/H², "no local
  stationary-phase obstruction occurs at θ = 1/2"; both Hankel branches interfere (upper envelope suffices); needs localisation and
  interchange, uniform Hankel/stationary-phase expansions in p, q, n (small n/N0 too), summable remainders and large-n tails, complex
  shifts and shift derivatives, explicit residues; (3) no intrinsic sub-1/2 failure; N(T+H) − N(T) = (H/2π)log(T/2π) + O(H²/T + log T);
  C2 follows if C1 is proved ((1 + 4ν)/3 = 0.494666…).

## B3 — short windows: "the method's barrier θ = 1/2" (OPEN)
Step fixed: SL Thm 1.2 = E1-M, the Gaussian-weighted mollified second moment in (T, T + T^θ], at the Taylor step SL (11): exp(iφ),
φ = −2πξ(e^{2y/H} − 1 − 2y/H) ≍ T/H², expanded to order J; remainder T^{1+ν+ε}(T/H²)^J, small only for θ > 1/2.
Class: Taylor profiles A_j (factor (T/H²)^j), E1's trivial-Voronoi dual bound, SL's frozen degree-201 detector, Wang's c(θ), SL's parity transfer.
Freedoms the method's own lemmas admit but the barrier excludes: (i) SL's (10) is exact, so exp(iφ) can be kept (exact chirp profile);
(ii) ν re-chosen per θ; (iii) the detector re-optimised at the new ν; (iv) R. The barrier statement ("below θ = 1/2 the method stops on
both sides … the detector needs θ > 0.49384, the Taylor step θ > 1/2"; "1/2 is the method's barrier on both sides") uses the expansion.
Test: (i) keeps E1's per-term size and decay to θ = 0.349, τ = 64, T ≤ 1e8, real and complex shifts (tests 1, 3); with ν < (3θ − 1)/4
below 1/2 the detector side alone binds: θ_c = 0.49467 (frozen detector, ν = 0.121, R = 44/5), ≈ 0.4940 (one-piece optimum).
Verdict OPEN: float gain θ > 0.49467 for a positive proportion (h(0.4947) = 7.1·10⁻⁵) of simple critical zeros in windows T^θ, past 1/2.
Certificate steps: (1) Step 4″, a uniform bound |𝒯(n)| ≪ H(pq/(nT))^{1/4}(1 + n/N0)^{−K}T^ε for the exact profile, θ ∈ (1/3, 1), complex
shifts ≪ 1/L and shift derivatives, by joint (y, ξ) stationary phase at the reviewer's saddle with localisation, both Hankel branches,
remainders and the n-tail; (2) residues = exact Mellin value of the exact profile (main term √πH·c unchanged); (3) compact windows of
width T^{θ−η} with θ − η > 0.4947; (4) recheck SL §3's profile-dependent branches (ε = +1 primal, K-Bessel, off-band) for the exact profile;
(5) Arb: κ(121/1000) at R = 44/5 (nu12.py's construction) and k*(θ), h(θ; κ) at the chosen θ. Gain is small (a range record), but the
stated barrier is false as a method barrier: the moment side now reaches θ > 1/3 + 4ν/3; the limit is the detector (θ_c* ≈ 0.494 one-piece).
Novelty: lst.py 'short interval|simple zero|Levinson|critical line|mollif|stationary phase|zeros of the Riemann' math.NT (2026-10, -11;
20 titles, none on simple zeros in windows below √T); Zenodo '"short intervals" "simple zeros"': only SL's 23134787, 23135354 (known);
E1.md §6's scan (12 Nov; Karatsuba's θ > 27/82 is for zeros of odd order, not simple); rg of our files: the Taylor step is recorded as a
barrier (E1.md §5, E1Audit FULLPROOF l.81, the discovery page), no exact-profile treatment. Logs: b3_profile_real*.log, b3_profile_cplx.log,
b3_onset.log, b3_ceilnu.log, b3_h.log, reviewer_b3.txt; peaks ≤ 134 MB (one run killed by memcap at 432 MB before a grid fix, no result used).

## B1 — GL(2) one-piece Levinson–Conrey ceiling 0.07005633 (0.08971026 under Ramanujan) (ROBUST)
Step fixed: mollified second moment ⇒ zero count, κ_f ≥ 1 − R⁻¹log c(P, Q, R, ν/2): Littlewood on [a, ·] × [T, 2T], a = ½ − R/L, then
Jensen (mean log|VM| ≤ ½ log mean|VM|²), with Bernard 1404.6432 th:mm2 at Andersen–Thorner's ν < 25/128.
Class = Bernard's input verbatim (gl2audit): real P, P(0) = 0, P(1) = 1; Bernard's complex Q of any degree, any β; any R; any ν < 25/128.
Freedoms of the cited lemmas and their fate:
- real mollifier shifts n^{−c/L}: P(x)e^{−cν(1−x)} is real with P(0) = 0, P(1) = 1, so inside the class (the ceiling is a sup over H¹);
- Q of any degree and β, R, ν: inside; derivatives of V: Q's degree; derivative scales other than 1/L break Q(y) + Q(1−y) ≡ β (counting device);
- imaginary shifts = complex P: not admitted (Bernard: real P); sesquilinear extension: complex P with real Q adds ∫∫e^{2rv}B² ≥ 0 (proved);
  jointly with complex Q: no gain in floats (gl2audit cplx.py, cplxglobal.py; my low-degree Ritz b1_complex.log: optimum stays real);
- a different inequality in the same step: none, given only second moments — Jensen is equality for constant |VM|, and Littlewood at
  several a gives min_a U(a)/(½ − a), attained by the configuration with all n₊ zeros of V at β = ½ (meets every bound U(a) at once);
- several pieces (KRZ's GL(2) two-piece theorem, a different input): ≤ 2.1·10⁻¹¹ (gl2two); Feng-type pieces, ν ≥ 25/128: no GL(2) mean value.
Verdict ROBUST: every freedom of Bernard's lemma is inside the certified class; the step's inequalities are sharp for the data used;
extensions need new mean values and add nothing in floats.

## B2 — small gaps: 7/11 = 1/(1 + θ), θ = 4/7 (ROBUST in its class)
Steps fixed: (i) Theorem U′'s length θ < κ/2: Lemma M gives ∫|R_T|^q ≪ T for q < κ/θ, Step 7′ needs q > 2 (hairpin weight e^{2θa} per
depth bin against Conrey's e^{−κa}, κ = 8/7 − δ = 2·(4/7), Conrey's mollifier length); (ii) DPC's sign condition φ(1 + θ) > 1; (iii) the
DPP limit of Lamzouri-type LPs (explore34 Theorem A floor, explore40's realizability cap).
Class: zero data = band-limited statistics of the ordinates of bandwidth ≤ 4/7 (pair, two-line, n-level), plus a density slope ≤ 8/7.
Mechanism: the discrete sine DPP on (7/11)ℤ (unit intensity) shares every correlation of bandwidth ≤ 1/a − 1 = 4/7 with the sine
process (F_DPP = min(|α|, 1) on |α| ≤ 4/7, = 4/7 on [4/7, 1]; b2_density.log) and has no gap below 7/11; DPC's φ(1 + θ) > 1 is this
limit at every θ (θ = 1 under RH: ½). So the three limits are one quantity: Conrey's length 4/7.
Freedoms: (a) non-uniform density (Conrey near the line, Ingham/Huxley/Guth–Maynard deeper): fails, hairpin mass ~ (log T)^c once
qθ > 8/7 (θ = 0.58, q = 2.01: 10^{0.33} → 10^{3.85} for L = 10² → 10⁸; θ = 0.56: O(1)); (b) q = 2 in 7′ (mean value valid for θ < 1): no
uniform integrability of |R_T|² on small sets; large-value estimates (Huxley, Guth–Maynard) give, by a rough count, ∫|R_T|^{2+δ} ≪ T^{1+δθ/4+ε}, a power loss;
(c) BGSTB's complex-γ F ~ |α| on [4/7, 1], the one input the 7/11-DPP violates, reaches real gaps only through off-line weights
T^{2α(β−½)}, controlled only for α < κ/2 = 4/7 (cglunc's Prop O mechanism); (d) the arithmetic Λ-terms of the two-line moments are not
zero statistics the DPP shares, but no route uses them beyond main terms at bandwidth θ (not tested).
Verdict ROBUST (class: bandwidth ≤ 4/7 zero data, slope 8/7): passing 7/11 needs κ > 8/7 near σ = ½ with one log T (a longer
unconditional mollifier) or a uniform-integrability input for |R_T|² at θ > 4/7. Caveat (the page's own): θ < κ/2 is proved sufficient,
its necessity is the isolated-deep-zero heuristic.

## B4 — simple zeros on the line: TEN 153/227, AF 0.6818287, Prop O (ROBUST in their classes)
TEN 153/227 = 0.6740088. Step fixed: the universal configuration-wise affine-defect estimate D∘ ≥ c·r − B·span of the simple on-line
block, on windows cos(√2 s) + Σ_{j≤96}c_j cos(2πjs) (complete for even windows as J → ∞; the frame's real-symmetric Gram needs even φ),
λ = 1; certificate Prop D: one periodic configuration X with e_W(X) ≤ ν − H_W for every window (Arb LDLᵀ, min pivot 3.85·10⁻⁴, stable in J).
Freedoms excluded (cglunc's own table): one-window refinements (Q-side defects Σ(q_j − 2)², cross term 2tr(P₁Q₊), von Neumann slack),
float cap ≤ 0.6777 (record window), no family-wide cap; several windows as constraints on pair statistics (global LP): AF's 0.6818;
higher moments tr G̃^k: dead (AF §7(e)); CGL: dead (Prop O); bandwidth λ < 1 (admitted: X = e^L ≤ T): my test b4_lambda.log —
H loss (A − B)(1 − λ) = 0.685(1 − λ) (MT window; 0.667 for ψ₀; λ* = √(A/B) ≈ 1.77 > 1), against TEN's whole defect gain ≤ 1.5·10⁻³ and
AF's rank cap κ ≤ λ: helps only if the defect gain grows faster than 0.685 per unit of 1 − λ (not seen; not certified).
Verdict ROBUST in its class, correctly named; not a frame ceiling: 0.6740 → 0.6777 (one window) and → 0.6818 (global) stay open, and no
certificate shape beyond TEN exists even in floats (NPAFrame: none for stationary laws).
AF 0.6818287. Step: bandwidth-one data (trace, every window's HS norm, all cross terms); a periodic law (Lean PairCeiling given EnclOK).
Freedom: level-2 (triple) words with Σ|ξ| ≤ 1: float gain 6.0–6.7·10⁻⁴ over the bandwidth-one value in periodic models (NPAFrame,
N = 10, 12), killed unconditionally by K ≍ NT^{−2δ} off-line zeros at depth δ, allowed by slope 8/7 < 2 (ExploreNPA). Verdict ROBUST
for its class (a theorem); the unconditional escape needs off-line control (slope > 2 or exception-free boxes).
Prop O. Step: CGL's diagonal bound S(K) ≥ K(0)Σm² for K̂ ≤ 0 beyond |α| = 1; K(2iy) = ∫K̂cosh(4παy) is dominated by the negative edge
(≍ e^{4πAy}, A > 1), so a density slope > 2A > 2 is needed; cglunc already tested the lemma-191 freedoms (analytic weights, complex shifts,
sign-changing or complex windows, several windows, Gaussian and compact tails, more points). Verdict ROBUST (analytic).

## Suggested NOGO lines
- ceilaudit B1 (N8): gl2ceil's one-piece GL(2) ceiling 0.07005633 (0.08971026 Ramanujan) covers every freedom of Bernard's th:mm2 — real P
  with P(0) = 0, P(1) = 1 (real mollifier shifts included), complex Q of any degree and β, any R, ν < 25/128; Littlewood + Jensen are sharp
  given second moments; only new mean values (complex P with complex Q, Feng-type pieces, ν ≥ 25/128) could pass it.
- ceilaudit B2 (N8): the small-gap limit 7/11 holds for zero data band-limited to 4/7 (pair, two-line, n-level) plus density slope ≤ 8/7:
  the (7/11)ℤ discrete sine DPP shares all of them and has no gap below 7/11; deeper density estimates (Ingham, Huxley, Guth–Maynard)
  lose (log T)^c once qθ > 8/7; passing needs slope > 8/7 near σ = ½ with one log T, or uniform integrability of |R_T|² at θ > 4/7.
- ceilaudit B4-TEN (N8): 153/227 holds for universal configuration-wise affine-defect estimates on even windows at bandwidth λ = 1;
  λ < 1 loses 0.685(1 − λ) in H against a total defect gain ≤ 1.5·10⁻³; one-window refinements (float cap 0.6777) and global
  pair-statistics LPs (AF 0.6818) are outside the class.
- ceilaudit B4-AF (N8): 0.6818287 holds for bandwidth-one data; triples (float +6–7·10⁻⁴ in periodic models) die unconditionally by sparse
  deep off-line zeros (slope 8/7 < 2).
- ceilaudit B4-PropO (N8): CGL kernels (K̂ ≤ 0 beyond |α| = 1, support A > 1) need a density slope > 2A; lemma-191 freedoms tested in CGLUnc.
Not a NOGO line (B3 is OPEN): the short-window barrier θ = 1/2 of Theorem E1 is SL's Taylor step (11), not the method's; suggested
caveat for discoveries/simple-zeros-short-windows.md (orchestrator): "below 1/2 only the expansion stops; the exact profile keeps the
dual bound in floats; the method's limit is the detector side, θ_c ≈ 0.4947 (frozen detector), ≈ 0.4940 (one-piece optimum)".

RH relevance: N8

## Files
Riemann/Deep/CeilAudit/: REPORT.md; code/ b3_profile.py (exact-profile dual term, FFT in u = H(e^{2y/H} − 1)), b3_onset.py (frozen
detector onset), b3_ceilnu.py (one-piece optimum via GL2Ceil/ceil.py, read-only import), b2_density.py (DPP form factor, density mix),
b1_complex.py (joint complex Ritz), b4_lambda.py (bandwidth λ); logs/ b3_profile_real.log (first run, killed by memcap at 432 MB in its
τ = 64 row; rows above it valid), b3_profile_real2.log, b3_profile_cplx.log, b3_onset.log, b3_ceilnu.log, b3_h.log, b2_density.log,
b1_complex.log, b4_lambda.log, novelty_b3.log, reviewer_b3.txt, *.err (memcap peaks, /usr/bin/time -l). Scratch: .scratch/ceilaudit/STATE.md.
Run: cd Riemann/Deep/CeilAudit/code && bash ../../../../.scratch/memcap.sh 290 900 uv run --no-project --with numpy --with scipy python3
b3_profile.py 0.01 0.02 1e6:500,1e6:250 (≈ 1 min, ≤ 134 MB); b3_onset.py needs python-flint, mpmath and ../../SantiAudit/params_exp010.json.

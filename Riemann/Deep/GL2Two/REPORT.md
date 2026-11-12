NO RECORD: the jointly optimised two-piece (Kühn–Robles–Zeindler) mollifier adds 1.1·10⁻¹⁴ to κ_f at its published lengths ν₁ < 5/27, ν₂ < 25/149 (optimum 0.06642206929; the record 0.0700562 needs ν₁ up to 25/128) and 2.1·10⁻¹¹ under GRC (0.08971017895, below the one-piece ceiling 0.08971026): its cross term sees Q only in boundary layers.

# gl2two — does a second mollifier piece escape the GL(2) Levinson ceiling? (nodes N8, N12)

## Claim (numerical, float; not a certified barrier)
1. Functional. KRZ's two-piece mean value (arXiv 1609.03738, Thm unsmoothed l.331–343; c_{ℓ,ℓ} l.392–401, c_{ℓ,ℓ+1} l.415–421) is, term by term,
   Bui–Conrey–Young's ζ functional (arXiv 1002.4127, eq:c1 l.344, eq:c12 l.360–366, eq:c2 l.376–389) at R′ = 2R, θ_i = ν_i/2 (KRZ l.1749):
   c = c₁ + 2c₁₂ + c₂, κ_f ≥ 1 − log c/R′. Our implementation reproduces BCY's 0.4105, BCY's printed optimal P₂ to all 6 digits, BCY's κ* 0.4058 (after a misprint fix),
   and KRZ's 0.0297607 and 0.0693872 exactly (DERIVATION §2, §5).
2. Lengths. Published theorem statements justify (ν₁, ν₂) < (5/27, 25/149) unconditionally and (1/4, 1/4) under GRC, for holomorphic
   newforms of even weight, squarefree level, trivial character (KRZ's class). AT21's ν < 25/128 covers the ψ₁-term I₁₁ verbatim (its off-diagonal
   bound is per pair (a, b), coefficient-free), and KRZ's cross-term proof needs only ν₁ + ν₂ < 1, ν₂ ≤ ν₁; so (25/128, 25/149) follows from published
   estimates but is stated nowhere. Nothing justifies ν₂ > 25/149 unconditionally (DERIVATION §4).
3. Gain. At every one of these lengths the best second piece (any P₂ with P₂(0) = P₂′(0) = P₂″(0) = 0) at the one-piece optimum raises κ_f by
   ≤ 2.1·10⁻¹¹: 1.1·10⁻¹⁴ (5/27, 25/149), 2.1·10⁻¹¹ (1/4, 1/4; GRC), 2.7·10⁻¹⁴ (25/128, 25/149), 9.1·10⁻¹⁴ (25/128, 25/128, unpublished).
   Joint re-optimisation of (P₁, P₂, Q) from 5 starts at the two published lengths and an R scan do not change this (joint = first order within
   the 3·10⁻¹³ Ritz noise; at R′ = 3, where the gain is 5·10⁻⁶, they agree to 0.2%).
   At equal lengths ν₂ = ν₁ (the clean negative): 2.5·10⁻¹⁴ (5/27), 9.1·10⁻¹⁴ (25/128), 2.1·10⁻¹¹ (1/4). For comparison ζ at θ = 4/7: +1.7·10⁻³ (BCY).
4. Mechanism (DERIVATION §3). With w = e^{Ry}Q the exponentials in c₁₂ cancel: the cross term ⟨Vψ₁, Vψ₂⟩ sees w only on [0, θ₂] ∪ [1 − θ₂, 1]
   (w(0) = 1, w(1) = e^R Q(1) = O(1)), while ‖Vψ₂‖² sees the bulk ∫w² ∝ c₁ = e^{R′(1−κ)}. Short mollifiers (θ ≈ 0.1) force R′ ≈ 8–10,
   c₁ ≈ 10³–10⁴. Along the R scan the relative gain LᵀM⁻¹L/c₁ falls roughly like c₁⁻³ (local log-log slopes −2.7 … −3.5 for c₁ from 6 to 3·10⁴).
5. Side results. (a) KRZ's improvement over Bernard (2.97% → 2.97607%, 6.93% → 6.93872%) comes from their better P₁, Q, R: their ψ₂ (P₂
   coefficients ~10⁻⁵) contributes 1.6·10⁻⁸ and 7.8·10⁻⁹. (b) BCY l.294 misprint: the x⁵ coefficient of the κ* polynomial P₂ is .00769594, not .0769594
   (printed value gives 0.3935; corrected 0.4058028 = BCY's .4058). (c) KRZ l.320, l.1761 say θ = 0 is proved (Deligne), but the θ in
   (1 − 2θ)/(4 + 2θ) is the Maass-spectrum exponent of Bernard's shifted-convolution bound; KRZ's theorem itself says "under the Ramanujan conjecture".

## Functional (KRZ; ζ-normalised R = 2R_KRZ, θ_i = ν_i/2; DERIVATION §1–2)
ψ = ψ₁ + ψ₂, ψ₁ = Σ_{h≤T^{ν₁}} μ_f(h)h^{σ₀−1/2−s}P₁[h], ψ₂ = χ_f(s+1/2−σ₀)Σ_{hk≤T^{ν₂}} (μ_f∗μ_f)(h)λ_f(k)h^{σ₀−1/2−s}k^{1/2−σ₀−(1−s)}P₂[hk] (l.189–208);
Q(0) = 1, Q(x) + Q(1−x) const (l.172); P₁(0) = 0, P₁(1) = 1; P₂(0) = P₂′(0) = P₂″(0) = 0 (l.208, l.230); 0 < ν₂ ≤ ν₁ (l.208); κ_f ≥ 1 − log(mean)/(2R_KRZ) (l.260).
c(P₁,P₂,Q,R;θ₁,θ₂) = c₁ + 2c₁₂ + c₂ (l.331–343; l.1749):
- c₁ = 1 + (1/θ₁)∫₀¹∫₀¹e^{2Rv}(QP₁′ + θ₁Q′P₁ + θ₁RQP₁)² du dv (Conrey/Bernard; KRZ l.405);
- c₁₂ = 4(θ₂/θ₁)²e^R ∂²_{xy}[∫∫_{a+b≤1}∫₀¹u²(1−u)e^{R[θ₁(y−x)+uθ₂(a−b)]}Q(−xθ₁+auθ₂)Q(1+yθ₁−buθ₂)P₁(x+y+1−(1−u)θ₂/θ₁)P₂″((1−a−b)u)]₀ (l.415–421, ℓ = 1);
- c₂ = (2/3)∂⁴_{x²y²}[∫_{[0,1]⁴}(1−r)⁴(1/θ₂+S)e^{−θ₂RS}Q(θ₂(−y+u(x+r))+t(1+θ₂S))e^{2Rt(1+θ₂S)}Q(θ₂(−x+v(y+r))+t(1+θ₂S))(x+r)(y+r)
  P₂″((1−u)(x+r))P₂″((1−v)(y+r)) dt dr du dv]₀, S = x+y−v(y+r)−u(x+r) (l.392–401, ℓ = 2).
These are BCY's eq:c1, eq:c12, eq:c2 verbatim with θ_i = ν_i/2. KRZ l.405 writes c₁,₁ with ν₁ for θ₁; the reading above is forced by Bernard's
one-piece mean value c(P,Q,2R,ν/2) and by KRZ's own printed numbers.

## Admissible lengths (DERIVATION §4)
| tier | (ν₁, ν₂) uncond. | GRC | source |
|---|---|---|---|
| stated theorem | (5/27, 25/149) | (1/4, 1/4) | KRZ Thm unsmoothed l.331–336 (θ = 7/64 / 0): I₁₁ l.822–832, I₂₂ l.838–860 |
| consequence of published estimates | (25/128, 25/149) | (1/4, 1/4) | AT21 Thm thm:nu l.227–231 for I₁₁ (Prop. prop:off-diag l.317–326 is per pair, any a, b); KRZ's I₁₂ proof uses only ν₁+ν₂ < 1, ν₁+ν₂ < 2, ν₂ ≤ ν₁ (l.347, l.1278, l.1392, l.1418) |
| plausible, unpublished | (25/128, 25/128) | — | AT21's per-pair bound inserted into KRZ's I₂₂ |
All for KRZ's class (holomorphic, even weight, squarefree level, trivial character). The longest the published theorems justify: tier 1.

## Optimisation (float; ξ = ν/2 − 10⁻⁷; one-piece values certified in Degree2Mol/logs/cert.log)
| (ν₁, ν₂) | justified by | R′* | c₁ | one piece κ₁ | min_g: LᵀM⁻¹L | two-piece κ₂ − κ₁ (first order) | joint alternating | multistart |
|---|---|---|---|---|---|---|---|---|
| (5/27, 25/149) | KRZ theorem, uncond. | 10.4093 | 16612.41 | 0.0664220692903 | 1.97·10⁻⁹ | 1.1·10⁻¹⁴ | < 4·10⁻¹² (Ritz noise floor) | 4 starts: ≤ 0.0664220692898 |
| (1/4, 1/4) | KRZ theorem, GRC | 7.6881 | 1094.85 | 0.0897101789336 | 1.78·10⁻⁷ | 2.11·10⁻¹¹ | 2.07·10⁻¹¹ (κ₂ = 0.0897101789549) | 4 starts: ≤ 0.0897101789550 |
| (25/128, 25/149) | published estimates, not stated | 9.8678 | 9667.22 | 0.0700562390188 | 2.58·10⁻⁹ | 2.7·10⁻¹⁴ | — | — |
| (25/128, 25/128) | unpublished | 9.8678 | 9667.22 | 0.0700562390188 | 8.66·10⁻⁹ | 9.1·10⁻¹⁴ | — | — |
| (5/27, 5/27) | equal lengths, hypothetical | 10.4093 | 16612.41 | 0.0664220692903 | 4.37·10⁻⁹ | 2.5·10⁻¹⁴ | — | — |
| KRZ's own P₁, Q, R | KRZ §5 | 6.42 / 5.65 | 507.2 / 192.1 | 0.0297607 / 0.0693872 | — | 1.6·10⁻⁸ / 7.8·10⁻⁹ | — | — |

R scan (one-piece optimum at each R′ + best P₂; rscan.py), lengths (5/27, 25/149) — κ₁ / gain:
R′ = 1.5: −0.324 / 7.6·10⁻⁴; 3: −0.0892 / 5.4·10⁻⁶; 5: 0.0226 / 2.6·10⁻⁸; 7: 0.0560 / 1.6·10⁻¹⁰; 9: 0.0652 / 7.9·10⁻¹³; 10.41: 0.06642207 / 1.1·10⁻¹⁴;
12: 0.0655 / < 10⁻¹⁵. GRC (1/4, 1/4): 1.5: −0.196 / 9.2·10⁻⁴; 3: −0.0024 / 7.3·10⁻⁶; 5: 0.0737 / 3.0·10⁻⁸; 7.69: 0.08971018 / 2.1·10⁻¹¹; 10: 0.0858 / 6·10⁻¹⁴.
The second piece only matters where κ₁ < 0; κ₁ + gain peaks at the one-piece R′*. Joint re-optimisation at R′ = 3 (gain 5.4·10⁻⁶): 5.416·10⁻⁶ vs
first order 5.406·10⁻⁶, so first order is accurate to 0.2% even there.

Against the thresholds of the brief: unconditional, no two-piece value is near 0.0700562390 at published lengths (0.06642207 < 0.0700562); at the
tier-2 lengths (25/128, 25/149) the two-piece value exceeds the one-piece value at the same ξ by 2.7·10⁻¹⁴. Under GRC the two-piece value at
ξ = 1/8 − 10⁻⁷ passes the literal threshold 0.0897101789 (record 0.0897101789336) by 2.1·10⁻¹¹ — 3,400 times less than the 7.2·10⁻⁸ that one
piece alone gains by moving ξ from 1/8 − 10⁻⁷ towards 1/8 (float ceiling 0.0897102510, certified ceiling ≤ 0.08971026), and only for KRZ's
narrower class of forms. The two-piece supremum over the open range (≈ 0.0897102510 + 2·10⁻¹¹) stays below the certified one-piece ceiling.
So the second piece does not escape the ceiling; a 2·10⁻¹¹ "record" would be dominated by the one-piece method itself, so no Arb certificate
was built (it would need rigorous 4D integration of c₂).

## Certificate
None: no record. The compared one-piece numbers are certified (Degree2Mol cert.py: 0.0664220692902798, 0.0897101789336033, 0.0700562390188180;
GL2Ceil ceilings ≤ 0.07005633, ≤ 0.08971026). The two-piece gains are float numbers validated as below.

## Checks
- Literature reproduction (val_bcy.py, val_bcy2.py, val_bcy3.py → logs/val_bcy*.log): BCY 0.4105121 (paper .4105; printed Q has Q(0) = 1.002,
  rescaled to 1), BCY's printed P₂ = the computed minimiser to all printed digits; κ* 0.4058028 (paper .4058, with the misprint corrected);
  KRZ 0.0297607, 0.0693872 (paper .0297607, .0693872).
- Derivative algebra (xcheck.py → logs/xcheck.log): BCY's integrands coded literally from the TeX, derivatives by Cauchy's formula in complex x, y:
  c₁₂, c₂ agree with the jet code to ≤ 7·10⁻¹⁴ relative.
- Convergence (gain.py, resol.py → logs/gain.log, logs/resol.log): gains identical to all printed digits at three quadrature resolutions and for
  K = 4, 8, 12 basis terms of P₂″. LᵀM⁻¹L is computed directly (no cancellation against c₁), so its float error is relative, not 10⁻¹⁶·c₁.
- Joint optimisation (joint.py → logs/joint_*.log): exact alternating solves in (P₁, P₂) and Q at R′*, Q-form consistency vᵀH₁₂v = gᵀL, vᵀH₂v = gᵀMg
  to 10⁻¹⁷; control run with ψ₂ ≡ 0 in the same code (logs/joint_*_control.log) fixes the noise floor (GRC 3·10⁻¹³, uncond. 4·10⁻¹²).
  Multistart (ms.sh → logs/multistart.log; starts: KRZ's Q, two random Q, the R′ = 3 optimum; first step = exact (P₁, P₂) solve):
  all 8 runs climb to the same optimum (best 0.0897101789550 GRC, 0.0664220692898 uncond.; slow starts still 10⁻¹¹ below at iteration 4).
- Reproducibility: val_bcy*.py, xcheck.py, gain.py rerun from this directory (.scratch/gl2two/rerun.sh): every value identical; only the
  rounding-level asymmetry diagnostic of M changes (t-chunked summation order). gain.py uses θ₂ = 25/298 itself, resol.py 25/298 − 10⁻⁷: same digits.
- Cross-term hypotheses audited by rg over KRZ l.1216–1640 for every constraint on ν: only ν₁+ν₂ < 1 (l.1278), ν₁+ν₂ < 2 (l.1392), ν₂ ≤ ν₁ (l.1418).
- Complex Q never helps with real P₁, P₂ (DERIVATION §6). Resources: one thread, every kept run < 300 MB except one: the first resol.py run
  peaked at 389 MB at its finest resolution (cap breach; c2_mat now chunks t and self-aborts above 280 MiB; the rerun, identical output,
  peaked at 267 MiB: logs/resol.time; first run kept as logs/resol_first_run_389MB.*).

## Review
Reviewer (decisive chain pasted inline, .scratch/gl2two/reviewer_brief.txt): "GAP (global optimality and maximal admissible lengths
unproved) … ROUTINE (useful numerical audit, but neither a certified barrier nor an improved record)." Supported in its words: "No record found in the
tested families; negligible improvement near the one-piece optimum." Response: wording now says numerical and local (descent + R scan +
multistart), not "cannot"; lengths split into explicit statements vs consequences of published estimates (the reviewer's distinction; the cross-term
hypotheses were audited: ν₁+ν₂ < 1, ν₁+ν₂ < 2, ν₂ ≤ ν₁ only, KRZ l.1278, l.1392, l.1418); multistart added; complex Q handled by the Hermitian split
(DERIVATION §6). Not done: a certified two-piece ceiling (a dual for the bilinear (P, Q) problem), so the negative stays numerical.

## Novelty
- python3 .scratch/lst.py 'mollifier|critical line|cusp form|Levinson' --zenodo '"two-piece mollifier"' (math.NT, CA, CV, SP, PR, FA, CO;
  2026-10, 2026-11; logs/novelty_lst.log): 13 title hits, none on GL(2) multi-piece mollifiers (2610.08965, 2610.23467 are ζ); Zenodo 24101490
  ("xiprime-two-piece-mollifier", a Lean repository claiming 88.45% for ζ under axioms) is unrelated to GL(2).
- OpenAlex citers of KRZ (Math. Z., W2726356247; 7; logs/novelty_openalex.log): Pratt–Robles, PRZZ, Sono (generalised mollifiers for ζ), AT21,
  two Dirichlet-L papers, one extreme-value paper: none quantifies the GL(2) two-piece gain; AT21 uses one piece.
- KRZ itself: P₂ ~ 10⁻⁵ in both parameter sets (l.1758, l.1771) without comment; it attributes its gain to ℒ = 2 (l.249, l.291). Repo: R446
  (TREE.md) noted KRZ as "unoptimised"; no earlier computation.
- Verdict: the size of the two-piece gain for GL(2) (≤ 2·10⁻¹¹) and its boundary-layer cause appear not to be in the literature; this is a
  negative (NOGO), not a discovery.

## Caveat
- Numerical: global optimality over (P₁, P₂, Q, R′) is not proved (biconvex problem; descent from several starts + R scan); float64 throughout.
- Scope: BCY/KRZ's ψ₂ = χ_f·(μ_f∗μ_f)·λ_f shape only; KRZ's ψ_ℓ, ℓ ≥ 3 (l.213–237) carry χ_f^{ℓ−1} and are weaker still (KRZ l.280–281 quote
  Sono: ψ₃ < 10⁻⁴ for ζ). Feng-type pieces (μ_f∗Λ_f^{∗k} at the same length), Hecke-twisted or derivative-based mollifiers are outside this
  computation; for ζ Feng's pieces gave +0.004 at θ = 4/7, and we know of no GL(2) mean value theorem for them.
- KRZ's class is narrower than AT21's (holomorphic, squarefree level); θ = 0 means GRC for Maass forms.

RH relevance: N8, N12 — closes the multi-piece (BCY-type) escape route from the GL(2) one-piece ceiling: at short mollifier lengths the second
piece is worth ≤ 2·10⁻¹¹, so N12 still needs a longer admissible ν or a different mollifier family, not more pieces of this shape.
Suggested NOGO line: "gl2two: BCY/KRZ two-piece mollifier for GL(2) adds ≤ 2·10⁻¹¹ to κ_f at all admissible lengths (cross term sees only the
boundary layers of e^{Ry}Q; relative gain ≈ c₁⁻³)."

## Files
- tp.py — two-piece functional (c₁, c₁₂, c₂; jets), the trusted-by-validation evaluator. onepiece.py — Ritz one-piece optimiser (copied from Degree2Mol/opt2.py).
- val_bcy.py, val_bcy2.py, val_bcy3.py — BCY/KRZ reproduction; xcheck.py — Cauchy-formula check of the derivative algebra.
- gain.py, resol.py — first-order gains at the one-piece optima, resolution check; rscan.py — R scan; joint.py, ms.sh — joint alternating optimisation, multistart.
- logs/ — all outputs and /usr/bin/time -l files. DERIVATION.md — formulas with TeX lines, lengths audit, validation, complex Q.
- Run (each < 300 MB, one thread): cd Riemann/Deep/GL2Two && uv run --no-project --with numpy python3 <script>.py [args]; params from Degree2Mol/params.json (read-only).

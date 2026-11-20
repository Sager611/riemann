BARRIER: finite free probability gives no new N6 or N2 handle: Gaussian-centred cumulant-box (CLT) certificates need |ρ₃(n)|² ≤ 1/(9(d−2)), i.e. n ≳ 4d (≈ 17d^{4/3} with the computed basin), far behind our Kim–Lee wedge; the de Bruijn–Newman flow is not a ⊞-Hermite flow on J^{d,n}; and the recorded shift-0 "4×" (Spark, 6.398·10²⁵) is a change of Jensen family, not a gain on GORZ's J^{d,0} (node N6).

(Details of the first line: the Hankel obstruction is exact; the rate ρ₃ ≈ (2/3)/√n is computed to n = 1150 and heuristic beyond.)

Short answer (pyramid):
- C3, audit, actionable. TREE N6, COROLLARIES l. 60 and R483 compare a t-family degree (Pólya's Jensen polynomials of Ξ in t)
  with our x-family degree (GORZ). The families are ⊠_d-intertwined: P_{2d} = J^{d,0} ⊙_d M_d (exact), so x ⇒ t and not
  conversely. The exact single-pair counterexample is w = 5 − i: J^{d,0} real-rooted iff d ≤ 26, P_m iff m ≤ 103.
  The x-family record stays 1.6144·10²⁵. Spark's result is ≈ 2× (1.98–2.0×; not 4×) the t-range implied by ours, 3.229·10²⁵.
- C2, brief target (b). Finite free cumulants of Jensen polynomials decay at the CLT rate: ρ₃ = (2/3 + o(1))n^{−1/2}, the
  same for ξ and DH; the limit 2/3 is heuristic (GORZ asymptotics, not proved here). The Hermite cubic basin is r₃(d) ≈ 0.16d^{−2/3}
  (exact for d ≤ 60, extrapolated beyond); rigorously r₃(d) ≤ 1/(3√(d−2)) (Hankel). Any cumulant-box (CLT-type) certificate centred
  at He_d therefore needs |ρ₃(n)| ≤ r₃(d). This is a necessary obstruction for that class, certified for n ≤ 1150, and asymptotically
  it means n ≳ 4d (rigorous form) or ≈ 17d^{4/3}. It is dominated by the existing explicit wedges.
- C1, brief target (a). FALSE PREMISE: the d = 4 heat/affine invariant ℓ₃⁴/ℓ₄³ drifts along the de Bruijn–Newman time.
  Only the d → ∞ limit holds: e^{−(τ/d²)∂²}P_d(z/d) → Ξ_τ. The backward time s(d,n) is a Γ-factor quantity, the same for DH and
  flat in τ (≤ 10⁻³): no Λ information.
- C4. explore15's screened sector lemma (neighbours heal an off-line pair) survives exact adversarial model tests:
  every S(T)-admissible configuration tried has J^{d,0} threshold ≥ 1.596× the isolated-pair (sector-lemma) value.
  A record lead (≈ 2.6·10²⁵), not RH, and not a proof.

# explore65 (Orca exploration): finite free probability → N6 (Jensen polynomials), with the brief's N2 premise tested

Tool transferred: the finite free structure of Jensen polynomials (MSS ⊞_d / ⊠_d; Arizmendi–Perales cumulants; Hermite flow).
- Reversed Jensen polynomial = symbol acting on x^d: J*^{d,n}(x) := x^d J^{d,n}(1/x) = F^{(n)}(D)x^d, F(x) = Σ γ(k)x^k/k!.
  So J* = ⊞_d over the zeros of F^{(n)} of rank-one/pair factors, heat flow = ⊞_d Hermite, log F^{(n)} (mod y^{d+1}) = additive
  finite free cumulants ℓ_k(n) (independent of d for k ≤ d).
- Jensen families are ⊠_d-images of each other (Schur–Szegő composition = MSS ⊠_d).
Data: γ(k) of ξ and of Davenport–Heilbronn (DH), k ≤ 1200, 60 digits (gamma.py; two Gauss–Legendre orders agree to 1e-59).
Normalisation: our F(x) = 2ξ(1/2 + √x/2), i.e. γ(k) = γ_GORZ(k)/4^{k+1}; real-rootedness is unchanged by γ(k) → cλ^kγ(k).
DH: F_DH(x) = Λ_f(1/2 + √x), Λ_f completed DH; all its γ_DH(k) > 0 (k ≤ 1200), like ξ.

## C3 (audit of a recorded N6 input): the "4× at shift 0" from Spark is a change of Jensen family, not a gain

Two RH-equivalent families: x-family J^{d,n} (GORZ; the N6 node; record d ≤ 1.6144·10²⁵ for all n) and t-family P_m
(Pólya/Spark: degree-m Jensen polynomial of the even function Ξ in t; Spark's Q_{n,t} control P_{2n}).
- Exact identity (c3_family.py): P_{2d}(X = −t²) = J^{d,0} ⊙_d M_d, M_d(X) = Σ_j C(d,j) 4^j (d−1/2)_j↓ X^j (⊙_d = coefficientwise
  product in the binomial basis = MSS ⊠_d). M_d(X) = d!(4X)^d L_d^{(−1/2)}(−1/(4X)) has only negative roots for every d (exact
  check d ≤ 60; reviewer). Schur–Szegő gives J^{d,0} real-rooted ⇒ P_{2d} real-rooted in X. Real roots in t need nonpositive
  X-roots, which positive coefficients (γ > 0 for ξ) supply: without them the implication fails (F = 1 − x, reviewer).
  The inverse multiplier is not real-rooted (every d = 2..60): no converse.
- Single pair F(x) = (1 + x/w)(1 + x/w̄), w = 5 − i (sin² arg w = 1/26), exact root counts:
  x-family J^{d,0} real-rooted ⟺ d ≤ 26 (GORTTW Thm 3.6, the sector lemma is sharp); t-family P_m real-rooted ⟺ m ≤ 103,
  exactly Spark's quartet criterion sin² ≤ (4m−6)/(m(m−1)). At d = 27: J^{27,0} has 0 real roots, P_54 has 4 of 4.
- So at a fixed zero angle θ the t-degree threshold is ≈ 4× the x-degree threshold (4/sin²θ vs 1/sin²θ). Spark's
  6.398·10²⁵ (and (U/a)² − 1 ≈ 6.4568·10²⁵) is a t-degree. Its x-content is d ≤ ≈ 1/sin²θ_max ≈ 1.6144·10²⁵: no change.
- Within the t-family the fair comparison is our x-range carried over by the identity: t-degree ≤ 2·1.6144·10²⁵ ≈ 3.229·10²⁵.
  Spark's t-range is ≈ 2× that (6.398/3.229 = 1.98), not 4×. Rolle (P^{(n)} of a hyperbolic P_{m+n} is hyperbolic) also gives every t-shift n with
  m + n ≤ 6.4568·10²⁵, so "shift 0 only" is not the real limit of that result.
- Source check: Spark's TeX (Zenodo 24242312, DBN_Simplicity_Transfer.tex l. 806) defines J_{d,t}(z) = Σ_k C(d,k) H_t^{(k)}(0) z^k,
  the degree-d Jensen polynomial of the even function H_t in its own variable z. That is the t-family.
- Consequence for TREE N6 (and COROLLARIES line 60, R483): read "at shift 0, J^{d,0} hyperbolic for d ≤ 6.398·10²⁵" as
  "t-family P_m (Pólya's Jensen polynomials of Ξ_t in t), m ≤ 6.398·10²⁵". The GORZ J^{d,0} record stays 1.6144·10²⁵.
  Spark's mathematics is not in question. Only our comparison mixed the two families.

## C1 (brief target a, N2): the de Bruijn–Newman flow is not a finite free heat flow on J^{d,n}

Claim tested: J^{d,n}_t = (affine map) ∘ (J^{d,n}_0 ⊞_d Hermite_{s(t)}), in X or in the reversed variable.
- Test (ffcum.py C1, d = 4): ⊞_d Hermite changes only the second cumulant ℓ₂, affine maps scale ℓ_k by λ^k, so
  I = ℓ₃⁴/ℓ₄³ must be constant in t. It is not (60-digit data, Ξ_t = ∫Φ e^{tu²} cos(zu), Polymath time):

  | n | I_rev, t = 0 | t = 0.2 | I_orig, t = 0 | t = 0.2 |
  |---|---|---|---|---|
  | 0 | −1.1712216 | −1.1687871 | −1.8766884 | −1.9940559 |
  | 10 | −2.5694947 | −2.5650118 | −782.23915 | −789.59762 |
  | 1000 | −121.41616 | −121.40305 | −3242.221 | −3241.7343 |

  The changes are 10⁻⁴ to 6·10⁻², against a numerical error below 10⁻⁵⁰. So the premise of target (a) is false (FALSE STEP).
  In GORZ coordinates, γ_t(n) = Σ_m (4t)^m (n+1/2)_m γ(n+m)/m! (our normalisation; t^m in place of (4t)^m in GORZ's). This is a Laguerre-type
  operator that reaches far up in n. It is not a shift-convolution.
- What is true (a limit, not exact): for the t-family, e^{−(τ/d²)∂²}P_d(z/d) → Ξ_τ(z) locally uniformly as d → ∞, for every
  real τ. Proof: the coefficient of z^i is Σ_m (−τ)^m Ξ^{(i+2m)}(0)/(m! i!)·(d)_{i+2m}↓/d^{i+2m}, and that ratio lies in [0, 1] and tends to 1, then dominated convergence (Ξ_τ is entire).
  Hurwitz then gives Λ ≤ liminf_d Λ_d, Λ_d := inf{τ : e^{−(τ/d²)∂²}P_d real-rooted}. RH ⟺ Λ_d ≤ 0 for every d (Pólya), so
  this is an RH-equivalent tautology at the level of signs. With Rodgers–Tao (Λ ≥ 0): RH ⇒ the backward finite free time of
  P_d is o(1/d²) (ROUTINE). A finite d gives neither an upper nor a lower bound on Λ: the lower direction would need
  limsup Λ_d ≤ Λ, and edge (unresolved) roots block that.
- The brief's s(d,n) (largest backward Hermite time keeping J* real-rooted, in units of its own variance) is a function of
  (ρ₃, …, ρ_d)(n) only (next section). It tends to the full variance as n → ∞ (the Hermite limit). It is a Taylor-coefficient
  (Turán-type) quantity and the same for DH (ρ₃√n: ξ 0.598, DH 0.613 at n = 1150). On the "zero mesh margin" hint: Hermite
  has positive margin (its backward time equals its variance). Zero margin belongs to the Jensen kernel (1+x)^d (C4).
  Exact values (c1_backward.py, bisection with flint): σ*(d, n) for ξ at t = 0, ξ at t = 0.2, and DH:

  | d | n = 0 | n = 10 | n = 100 | n = 1000 |
|---|---|---|---|---|
  | 3 | 0.2603 / 0.2594 / 0.1265 | 0.4608 / 0.4604 / 0.4124 | 0.6999 / 0.6998 / 0.6887 | 0.85267 / 0.85267 / 0.85015 |
  | 8 | 0.0234 / 0.0232 / 0.0050 | 0.0795 / 0.0792 / 0.0540 | 0.2890 / 0.2889 / 0.2718 | 0.57533 / 0.57532 / 0.56914 |

  The de Bruijn–Newman time moves σ* by at most 10⁻³ at n = 0 and by 10⁻⁵ at n = 1000, smoothly, with no signature of Λ.
  ξ and DH tend to the same curve. s(d,n) is therefore a Γ-factor quantity, not a Λ probe: NO IDEA for target (a).

## C2 (brief target b, N6): finite free CLT certificates are capped at n ≳ 17 d^{4/3}

Coordinates: after an affine change of x, J*^{d,n} = exp(−D²/2 + Σ_{k=3}^d ρ_k(n) D^k) x^d, where ρ_k = ℓ_k/(−2ℓ₂)^{k/2} and
ℓ_k = [y^k] log F^{(n)}(y)/γ(n). The ρ_k(n) do not depend on d. Hermite He_d is the point ρ = 0.
- Rate (ffcum.py C2, logs/c2_rho.log). The pattern is the finite free CLT pattern ρ_k ≍ n^{−(k−2)/2}:
  ρ₃(0) = 0.212; ρ₃√n = 0.417, 0.548, 0.598 at n = 10, 100, 1150; ρ₄n → −1.02; ρ₅n^{3/2} → 2.08.
  Asymptotic: GORZ's log γ(n+j) − log γ(n) = Aj − δ²j² + g₃j³ + …, with δ² ~ 1/(2n) and g₃ ~ 1/(6n²). A Poisson-tilt
  expansion of the log-symbol gives ℓ₃/(−2ℓ₂)^{3/2} ≈ (g₃ + 2δ⁴)/(2√2 δ³) → (1/6 + 1/2)/√n = 2/(3√n) (heuristic: no remainder control).
  DH: 0.613 at n = 1150 with the same limit. The rate is a Γ-factor property, blind to zeros.
- Cubic basin of He_d (c2_basin.py, exact rationals + flint root isolation): r₃(d) = sup{ε : exp(−y²/2 + εy³)(D)x^d real-rooted}
  = 1/3, 0.1517, 0.04156, 0.02249, 0.01174, 0.01044 at d = 3, 4, 10, 20, 50, 60, i.e. r₃(d) ≈ 0.16 d^{−2/3}
  (r₃d^{2/3} = 0.166, 0.159, 0.160 at d = 20, 40, 60: Airy-edge scaling).
  Rigorous for every d, from the centred root moments m₂ = d−1, m₃ = −3ε(d−1)(d−2), m₄ = (d−1)(2d−3) and the Hankel
  inequality m₂m₄ − m₃² − m₂³ ≥ 0: r₃(d) ≤ 1/(3√(d−2)).
- Obstruction (reviewer's wording: a necessary obstruction for independent, Gaussian-centred cumulant-box certificates).
  A CLT certificate bounds each normalised cumulant separately, |ρ_k| ≤ ε_k (a box around He_d). A box that contains J* also
  contains (ρ₃(n), 0, …, 0), so it is valid only if |ρ₃(n)| ≤ r₃(d) ≤ 1/(3√(d−2)). Certified, from the computed ρ₃: at n = 1150
  no box works for d ≥ 360 (Hankel), and none for d ≥ 28 with the exact basin (r₃(25) = 0.01897 > ρ₃ = 0.01763 > r₃(30) = 0.01661).
  Asymptotically, if ρ₃√n → 2/3 is proved: n ≥ (4 − o(1))(d−2), and n ≈ 17d^{4/3} with the extrapolated basin.
  Correlated regions, non-Gaussian centres (e.g. a Laguerre/Marchenko–Pastur centre, Schieber's d ≍ n regime) or mixed
  coordinates are outside this class; the two-thirds and Kim–Lee wedges are such methods. The class allows every affine
  normalisation and every reference Hermite variance: ρ₃ is invariant under both.
- Comparison at the N6 frontier (open region: d > D_all = 1.6144·10²⁵, n < N_KL(d)):
  | method | needs n ≥ at d = 1.6·10²⁵ |
  |---|---|
  | finite free CLT box (extrapolated basin) | ≈ 17 d^{4/3} ≈ 7·10³⁴ |
  | finite free CLT box (Hankel bound; ρ₃ rate heuristic) | ≈ 4d ≈ 6·10²⁵ |
  | two-thirds wedge n² log(n+2) ≥ 24d³ | ≈ 3.4·10³⁷ |
  | explicit Kim–Lee wedge (ours) | 4.70·10¹³ |
  The CLT ceiling n ≍ d^{4/3} beats the two-thirds exponent 3/2 by only d^{1/6}. Kim–Lee beats the two CLT rows by 12 and 21 orders.
  Below D_all every (d, n) is already proved. So finite free CLT, with any explicit constants (the Lévy-distance CLT of
  2310.15485 has C_d growing with d), adds nothing to N6: BARRIER for the box class. The obstruction itself is the exact Hankel
  bound. The ρ₃ rate is computed for n ≤ 1150 and heuristic beyond. At the frontier the margin is ≈ 10⁶ with the
  Hankel bound and ~10¹⁰ with the extrapolated basin, given that rate. Unreversed J^{d,n}: ρ₃ ≈ (1/3)(n+d)^{−1/2} asymptotically
  (0.65/(3√(n+d)) at n = 1000, d = 4, 20, 60; c2_checks.py), so its ceiling is n ≳ 2–4 d^{4/3}. The basin is an interval
  (spot checks at ±k r₃/8 and at 1.125–4 r₃, d = 10, 40).

## C4 (falsification test of explore15's screened sector lemma, N6): neighbouring zeros heal an off-line pair; the lead survives

The sector lemma (Chasse; GORTTW Thm 3.6) is sharp for sector data alone. Its extremal configuration is one isolated pair:
d_max = ⌊1/sin²θ⌋ with θ = arg w (c3: w = 5 − i gives 26). Finite free view: J* = ⊞_d over pair and rank-one factors,
and the zero-mesh kernel (1+x)^d = J_d(e^x) carries no Hermite divisor. Newton's inequality fails by e^{−2τ} for
C(d,j)e^{τj(j−1)}. So Jensen's damping (d)_j↓/d^j is not (log-cylinder heat flow) ∘ (multiplier sequence) for any τ > 0, and
de Bruijn's healing cannot be read off the operator. explore16 found this first ([ExploreQDeform](../ExploreQDeform/REPORT.md) T3'). The question is whether the zeros themselves heal.
That is explore15's screened sector lemma ([ExploreLatticeGas](../ExploreLatticeGas/REPORT.md)), heuristic ×1.6 in the S(t)-worst
case. Their gap and cluster cases above T were not tested. This section tests them exactly.
Model (c4_heal.py, c4_extra.py): x-family J^{d,0} of F = Π(1 + x/ζ²) for finite zero sets at height γ₀ = 10⁴, spacing
Δ = 0.85 = 2π/log(γ₀/2π), offset β = 1.6Δ (ξ above T = 3·10¹² has a/Δ = 0.3733/0.234 = 1.6). Off-line pairs are
(γ ± iβ), one pair = two zeros of ζ, so pair heights are 2Δ apart. Exact rationals, flint root isolation, bisection on d.
| configuration (ζ-like counts) | d_max / d_isolated |
|---|---|
| isolated pair (= sector lemma) | 1.000 |
| double pair, isolated | 1.500 |
| pair among real zeros (spacing Δ, 80 zeros) | 4.335 |
| lattice of off-line pairs (all zeros off-line, 82) | 1.886 |
| real zeros below, off-line pairs above (ξ just above T) | 1.912 |
| one-sided gap 14.84Δ | 2.746 |
| deficit spread over ±25.7Δ (half density) | 2.285 |
| symmetric S(T)-worst gap (half-width 7.42Δ: |S| ≤ 6.42 at 3·10¹², Platt–Trudgian), real beyond | 1.596 |
| same gap, off-line pairs beyond | 1.607 |
| double pair in the worst gap | 2.020 |
| DH-like offset β = 0.21Δ, real neighbours | 1.132 |
| worst gap, zeros beyond doubled (same count) | 1.611 |
| worst gap, pair 3Δ off-centre | 1.686 |
| 8 random: Wigner spacings, random gap ≤ 7.42Δ, half the neighbours off-line (denser than ζ) | 1.749 – 4.866 |
Convergence in the zero count: real neighbours 4.056 → 4.294 → 4.335 (26, 62, 82 zeros); worst gap 1.52 → 1.66 → 1.596.
Reading. All 18 configurations with neighbours and ξ's maximal offset 1.6Δ heal by ≥ 1.596× (10 are count-calibrated to ζ at T).
This matches explore15's heuristic ×1.6. A small offset (0.21Δ) heals by only 1.13×, but its own isolated threshold is
(1.6/0.21)² = 58× higher, so the maximal offset decides. Healing grows with β/Δ: it is weak for DH-like small
offsets and strong at ξ's maximal offset 1.6Δ, the case that sets the record. So "sector + explicit zero counting" would
lift the x-family range from 1.6144·10²⁵ to ≈ 2.58·10²⁵ (×1.596). All shifts follow from Rolle (J^{d,n} ∝ D^n J^{d+n,0})
plus the Kim–Lee wedge (n ≥ 1.4·10¹⁴). explore15's record-type IDEA survives this falsification attempt. It is still not a proof. Missing:
(1) a rigorous healing lemma, e.g. an MSS barrier (Stieltjes transform of the real factors) bounding the imaginary part of the
    bad pair under ⊞ with rank-one factors, with infinitely many zeros and the far cluster at x = 0;
(2) the worst case over all S(T)-admissible configurations, not 18 samples (the minimum here is a sample, not a ceiling);
(3) the dependence of the angle on height (2a/γ) and the zeros below T (Platt's list).

## Adversarial (DH and positive-definite models)
- C2: ρ₃√n is 0.598 (ξ) and 0.613 (DH) at n = 1150, with the same limit 2/3. The CLT rate is a Γ-factor property, so the
  ceiling is blind to zeros, as a barrier should be: no CLT certificate can tell ζ from DH.
- C1: σ*(d,n) differs between ξ and DH only through low Taylor coefficients (conductor 5, Γ((s+1)/2)); at n = 1000 they agree to
  3·10⁻³ (d = 3) and 6·10⁻³ (d = 8). DH's off-line zero 0.8085 + 85.70i sits far beyond the degrees where s(d,n) is computable
  (its x-family J^{d,0} fails at d of order 2·10⁴: the isolated-pair value is 1/sin²(arg ρ²) = 1.93·10⁴, and the model heals 1.13× at DH's β/Δ = 0.21).
- C3: an algebraic identity valid for every even entire function; the counterexample is a single pair, the extremal model.
- C4: healing needs neighbours within ~β. It is 1.13× at DH-like offsets (β = 0.21Δ) and ≥ 1.596× at ξ's maximal offset
  (1.6Δ). It is a range (record) statement, so it need not separate ζ from DH.

## Novelty (arXiv API and monthly listings 2026-10/11; Zenodo; our reports)
- Finite free CLT, LLN and Poisson laws for Jensen polynomials: known (2410.06403: Hermite universality via finite free
  probability; Schieber, Zenodo 22936342, Sept 2026: Marchenko–Pastur for d/(n+d) → c > 0, semicircle for c = 0, conditional on
  hyperbolicity; 2605.31356: Appell polynomials f_n(∂)z^n ↔ free infinitely divisible laws). Farmer 2008.07206: large-n
  hyperbolicity is generic.
- The box-certificate ceiling (C2) is a negative and was not found stated; it is not claimed as a discovery.
- C3 uses Schur–Malo–Szegő (classical); the finding is an internal correction. C4's mechanism is explore15's, and the
  factorisation obstruction is explore16's.
- No entry for DISCOVERIES.md.

## Suggested NOGO line
- explore65 (Orca, exploration: finite free probability — ⊞_d/⊠_d structure, Arizmendi–Perales cumulants, finite free CLT and Hermite flow → N6, N2) BARRIER: Gaussian-centred cumulant-box (CLT) certificates for Jensen hyperbolicity need |ρ₃(n)| ≤ r₃(d) ≤ 1/(3√(d−2)) (exact Hankel obstruction; r₃ ≈ 0.16d^{−2/3} computed to d = 60); with ρ₃√n → 2/3 (computed to n = 1150, heuristic limit, ξ and DH alike) this needs n ≳ 4d (≈ 17d^{4/3} with the basin), far behind the Kim–Lee wedge; the de Bruijn–Newman flow is not a ⊞-Hermite flow on J^{d,n} (d = 4 invariant ℓ₃⁴/ℓ₄³ drifts), only its d → ∞ limit (τ/d²); s(d,n) is a Γ-factor quantity; audit: the shift-0 6.398·10²⁵ (Spark) is a t-family degree (P_{2d} = J^{d,0} ⊙_d M_d; w = 5 − i: x ≤ 26, t ≤ 103), so the x-family record is unchanged and Spark's t-gain is ≈ 2× (1.98), not 4×; explore15's screened lemma survives exact adversarial models (≥ 1.596×) [ExploreFreeProb](Riemann/Deep/ExploreFreeProb/REPORT.md)

## Requested TREE fix (orchestrator)
N6 line: "at shift 0, J^{d,0} of Ξ_t hyperbolic for d ≤ 6.398·10²⁵" → "at shift 0, Pólya's t-variable Jensen polynomials P_m of Ξ_t
(degree m) hyperbolic for m ≤ 6.398·10²⁵ (≈ 2× the t-range 3.229·10²⁵ implied by our x-range via P_{2d} = J^{d,0} ⊙_d M_d and γ > 0;
GORZ's J^{d,0} unchanged at 1.6144·10²⁵, where the sector lemma is sharp)". The same applies to COROLLARIES l. 60 (6.4568·10²⁵
is a t-degree) and to R483's "4×". Within the t-family, Rolle extends Spark's range to every shift n with m + n ≤ 6.4568·10²⁵.

RH relevance: N6 (barrier for finite-free CLT certificates; audit of the shift-0 record; explore15's screened lemma survives, record only); N2 (premise of the finite-free heat-flow route false; no Λ information).

## Files
Scripts in this directory, outputs in logs/:
- gamma.py → logs/gamma.json, gamma.log (γ(k) of ξ at t = 0, 0.05, 0.1, 0.2 and of DH, k ≤ 1200, 60 digits);
- ffcum.py (C1 invariant, C2 ρ_k) → c1_invariant.log, c2_rho.log; c1_backward.py → c1_backward.log;
- c2_basin.py → c2_basin.log; c2_checks.py → c2_checks.log;
- c3_family.py → c3_family.log, footprint_c3.log;
- c4_heal.py, c4_extra.py, c4_random.py → c4_heal_small.log, c4_heal_large.log, c4_heal_counted.log, c4_extra.log, c4_random.log.
Footprint: memcap peaks ≤ 31 MB (gamma.py), all others ≤ 16 MB; CPU ≈ 45 min in total, single-threaded.

NO DISCOVERY: the claimed zero-free half-plane Re s > 11/12 for Hecke L-functions of Q(√−3) with nonzero angular type, which gives L(E,s) ≠ 0 for Re s > 17/12 for the CM elliptic curves y² = x³ + D, was already on the Hexagon preprint server on 13 Oct 2026. Rubinstein–Restrepo Borrero (2610.00220) prove the elliptic-curve case by the same route and also 11/8. Hu (2610.00216) proves 7/8 for every unitary Hecke character of every number field. Our draft also had sign and normalization errors; these are now fixed, and the fixed formulas agree with Rubinstein–Restrepo Borrero. The record B♭ = 0.874886647 is not beaten.

# A known zero-free statement, repaired, and why the record stands (2026-10-21)

## What was claimed, and what is true
- **Claim.** For every unitary Hecke character ψ of K = Q(√−3) whose angular type m (the power of n/|n| in ψ) is not 0,
  L_K(s, ψ) ≠ 0 for Re s > 11/12. A Hecke character here means a multiplicative function on ideals that may include
  a factor (n/|n|)^m.
- **Consequences.** L(E, s) ≠ 0 for Re s > 17/12 when E/Q has CM by Q(√−3) (j = 0, 54000, −12288000).
  A weight-w CM newform has no zeros for Re s > (w−1)/2 + 11/12.
- **Already known, and stronger.** Both papers below are on Hexagon, 13 Oct 2026; Hexagon is an unreviewed preprint
  server like arXiv.
  - Rubinstein–Restrepo Borrero, 2610.00220: Part I gets 17/12 for y² = x³ + d by the same adaptation of OpenAI's 11/12
    proof, with a second derivative of the cubic theta function. Part II gets 7/8 for every such character with m = ±1,
    so 11/8 for every elliptic curve with CM by Q(√−3).
  - Hu, 2610.00216: 7/8 for every unitary Hecke character of every number field, any order. It also gives every CM
    motive, after the shift by w/2.
- **The draft overclaimed.** It said "first for any GL(2)". Finite-order characters already give GL(2) examples: their
  L-functions are those of weight-one newforms, already in OpenAI's paper and Liu's. The twists by N^{it} are trivial
  infinite-order examples. The only new case beyond OpenAI was m ≠ 0, and that is the case Hu and
  Rubinstein–Restrepo Borrero did.

## What was wrong in the draft (all fixed in NOTE.md)
- **The power.** The angular power is k = 1 + m, not 1 − m. Taking the complex conjugate when m < 0 makes k ≥ 2 always.
  This avoids the derivative-free case k = 0 and the negative-k formulas.
- **The multiplier.** In the theta transformation the cubic multiplier must be conjugated. Otherwise the local factors
  change and the quadratic large sieve no longer applies.
- **The dual factor.** The dual Mellin factor is (2π)^{−2s}, not (27/(2π)²)^s.
- **The principal case.** For ψ = N^{it}, L(s, ψ) = ζ_K(s − it) has its pole at 1 + it, not 1 − it. It is handled
  as in OpenAI's paper.
- **Checked.**
  - The repaired formulas match Rubinstein–Restrepo Borrero at k = 2: the chain-rule factor (cv)⁻⁴, the reflected Gamma
    kernel, the factor (2π)^{−2s} and the conjugated multiplier.
  - `angular_check.py` (16 s, 22 MB) checks numerically, to 1e−16, the chain rule for k = 1–4, the Mellin constants
    and the one-mode functional equation.
- **What is left.** A corrected write-up of the 11/12 route for every angular type. Its result is weaker than Hu's.

## Record hunt: B♭ not beaten
- **The two OpenAI methods are already combined.** The 7/8 paper contains the 5 Oct estimate: its (17.87), and its sixth-power amplification (Lemma 17.6). The 5 Oct method alone stops at 11/12; even a perfect mean square would give only 6/7.
- **Largest lever: the sixth-power term P_a^{1/6} in the low-side Gram bound** (7/8 paper, Prop. 15.2).
  - Value: with P_a^θ instead, B = 0.8724 at θ = 1/8 and 0.8637 at θ = 0. That is below the 13/15 floor, which itself assumes θ = 1/6 (`gram_theta.py`, 30 s, 30 MB).
  - Why it stays: the loss needs ξ⁶ = 1 and an unweighted Gauss modulus. The principal signal needs both, through the ζ(6z) pole and the ζ(w) pole at w = 1. Every variant that removes the loss also removes the signal (order-5 ξ, Möbius or character weights on the modulus).
  - Open lead: apply Cauchy–Schwarz only outside the resonant part, and bound the projection of the theta rows onto it.
- **Other families.** The sextic family is forced: Möbius absorption needs the cubic theta. Quartic Gauss sums are not theta coefficients, and the quadratic family is self-dual, so the argument is circular. Other base fields give the same exponents (Liu Part III). The method has no positivity step to change.
- **ζ′.** ζ′(s) ≠ 0 for Re s < 1 − B♭ at large height. This is routine (Levinson–Montgomery).
- **What Xu–Sang refine.** The 7/8 boundary, to 0.8748867, for the same families. Three steps do it: keeping the decay of rescaled prime tuples (prime length up to 1/5), averaging the numerator over each class with the amplified fourth moment, and choosing κ in two stages.

## Sources
- OpenAI 7/8 (30 Sep 2026) and 11/12 (5 Oct 2026; TeX: build/paper2.tex), both at github.com/openai/math/tree/main/preprints.
- Rubinstein–Restrepo Borrero, Hexagon 2610.00220v1; Hu, Hexagon 2610.00216v1 (hexagonmath.org, both 13 Oct 2026).
- Liu arXiv 2610.12234 and 2610.26835; Xu–Sang arXiv 2610.15291; Cureton–Thorner Hexagon 2610.00200v2; Frei–Sofos–Widmer arXiv 2610.18491.
- Deuring (L-functions of CM elliptic curves).

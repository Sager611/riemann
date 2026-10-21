DISCOVERY: every Hecke L-function of Q(√−3) whose character (Grössencharacter) has any order and any infinity type is zero-free for Re s > 11/12, so L(E,s) ≠ 0 for Re s > 17/12 for every elliptic curve over Q with CM by Q(√−3) (y² = x³ + D, x³ + y³ = 1). This is the first fixed zero-free half-plane for infinite-order Hecke L-functions and for any GL(2) L-functions. It is proved by adapting OpenAI's 5 Oct 11/12 proof. The record B♭ = 0.874886647 is not beaten.

# A new kind of zero-free statement, and why the record stands (2026-10-21)

## What is new
- **Theorem.** For every unitary Hecke character ψ of K = Q(√−3), of finite or infinite order, L_K(s, ψ) ≠ 0 for Re s > 11/12.
- **Consequences.**
  - L(E, s) ≠ 0 for Re s > 17/12 for all E/Q with j(E) ∈ {0, 54000, −12288000}.
  - L(f, s) ≠ 0 for Re s > (k−1)/2 + 11/12 for every weight-k newform with CM by K.
  - Eisenstein-prime angles equidistribute with a power saving (routine from the theorem).
- **Novelty.** A scout searched the arXiv API and the web: no prior fixed zero-free half-plane for these characters, for CM elliptic curves or for GL(2). Cureton–Thorner, Liu and Xu–Sang treat finite order only.

## Why the proof works (full argument in NOTE.md)
- **Setup.** Write ψ = ν0·α^m·N^{it}, with α(n) = n/|n| and ν0 of finite order. In OpenAI's 5 Oct mean-square proof the target only multiplies the Möbius coefficients. Poisson summation turns this into cubic-theta coefficients times ᾱ^k, with k = 1 − m.
- **The iteration keeps k fixed.** Both Poisson steps return the same k, so the whole proof runs with ᾱ^k in place of ᾱ.
- **The one new input.** The theta transformation is needed with a k-th horizontal derivative instead of the first. A chain-rule lemma handles it: at the cusp map only the pure ∂_w^k term survives, with factor (−(cv)^{−2})^k. The Gamma factors then shift by k/2, and the dual weights keep the same bounds, with constants depending on k.
- **The case k = 0 (m = 1).** It would leave the radial theta. Its constant term and Kubota–Patterson main term are not covered by the proof. Using the complex-conjugate character instead gives k = 2.
- **Checked.**
  - §§3–7 and the appendices of the preprint were read in full, with signs taken from its TeX source.
  - `angular_check.py` (3 s) checks numerically the Bessel–Mellin identity (to 1e−16) and the chain-rule lemma for k = 2, 3 (to 1e−31).
- **Caveats.**
  - The result rests on the unrefereed 5 Oct preprint, as every current record does.
  - It was checked by hand, not by machine.
  - The 7/8 and B♭ versions for these characters are not done.

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
- Liu arXiv 2610.12234 and 2610.26835; Xu–Sang arXiv 2610.15291; Cureton–Thorner 2610.00200.
- Deuring (L-functions of CM elliptic curves).

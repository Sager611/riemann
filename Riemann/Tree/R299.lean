/-!
# Round 299 — node N9 (moments, Lindelof), second round: records unchanged; the 7/8 half-plane affects only sigma > 7/8 (log zeta << (log t)^{8(1-sigma)+eps}, Titchmarsh 14.5)

* Node: N9 moments, Lindelof, subconvexity (second round).
* State of the art: TREE.md N9 records stand as far as checked (the listings scout timed out after locating
  2610.01035's sixth-moment corollary and 2610.12234); the convexity interpolation mu(sigma) <= (13/84)(7/8 - sigma)/(3/8)
  on [1/2, 7/8] is already logged (R280, COROLLARIES.md).
* Wild: a zero-free half-plane Re s > theta gives log zeta(sigma + it) << (log t)^{(1-sigma)/(1-theta)+eps} (Titchmarsh
  14.2, 14.5), so exponent 8(1 - sigma) for 7/8 < sigma < 1: classical, KNOWN; no explicit-constant version found
  (Fiori-Jaskari 2609.23222 gives explicit psi constants only).
* Wild 2: the 7/8 half-plane enters moments/large values only through N(sigma, T) = 0 for sigma > 7/8; 2610.01035's
  density consequence targets sigma = 3/4 and its large-value ranges are amplitude exponents (<= 13/84), so the sixth
  and twelfth moments and Halasz-Montgomery bounds do not improve: HOLDS.
* Adversarial: the claim-audit scout returned no answer.
* Exploration: Igusa/motivic zeta of plane curves - the maximal-order leading coefficient of a pole is a nonnegative
  incidence sum forcing Bernstein-Sato/monodromy zeros (2608.26087): a positivity/purity mechanism, not a location
  theorem (N12, weak). Stieltjes electrostatics: every simple-root polynomial has a tailored critical field, but not a
  stable minimum (2203.01419), so no real-rootedness certificate (N2/N6, weak). Inverse spectral geometry: a geometric
  Hilbert-Polya candidate must realise wave-trace singularities at lengths log p^k with non-cancelling amplitudes
  (2511.10398 rigidity), a constraint to reverse-engineer (N13, idea).
-/

namespace Fayed.R299
end Fayed.R299

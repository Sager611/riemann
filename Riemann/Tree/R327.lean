/-!
# Round 327 — node N9 (moments, Lindelof, subconvexity), round 4: records hold (13/84; de Faveri-Pandey moments); Tao-Trudgian-Yang exponent pairs noted; half-plane interpolation of mu(sigma) routine (R280)

* Node: N9 moments, Lindelof, subconvexity (fourth round; R314 before).
* [check] Subconvexity: mu(1/2) <= 13/84 (Bourgain 2017) stands; Tao-Trudgian-Yang (arXiv 2501.16779, Thm 'New exponent
  pairs', TeX 487-498) add four exponent pairs, e.g. (89/1282, 997/1282), none improving 13/84. Noted in TREE N9.
* [check] Moments: no Oct-Nov 2026 improvement of the sixth / twelfth moment exponents of 2610.01035 (2610.27591 is about
  negative moments): HOLDS.
* Computation: convexity between mu(1/2) <= 13/84 and mu(B) = 0 (B = 0.87338196661320752) beats the classical line through
  mu(1) = 0 by 0.005-0.037 on sigma in [0.55, 0.85]; already classified routine (R280, NOGO): KNOWN. Exponent-pair bounds not
  compared (scout omitted them).
* Adversarial: a zero-free half-plane Re s > B gives log zeta(sigma + it) << (log t)^{(1-sigma)/(1-B)+eps}, hence mu(sigma) = 0
  for sigma > B (Titchmarsh 14.2, 14.5; R299): HOLDS.
* Harvest: 2610.01035 states no open question.
* Exploration: finite free convolution (2410.06403): only generic Jensen universality, no xi factorization (no idea, R278);
  super-resolution dual certificates (2609.31299): separation alone does not tell zeta from DH (no idea); Cohn-Elkies /
  Delsarte LP (1908.04876, 2310.01913): off-line cross terms are indefinite, cannot certify Q(f) >= 0 (no idea, R40).
-/

namespace Fayed.R327
end Fayed.R327

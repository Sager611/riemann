/-!
# Round 346 — node N3 (Li/Keiper), round 5: records stand; Li ranges for Dedekind zeta of abelian fields by additivity (routine); Zenodo 23977400 flaw still open

* Node: N3 Li/Keiper coefficients (fifth round; R333 before).
* [check] Zenodo and arXiv since R333: 77.03 T^2 (Fabbian) and our 106.25 T^2 stand (newest arXiv match 2610.39481, 29 Oct,
  not applicable): HOLDS. Fabbian's own limits section asks for positivity beyond 77.03 T^2 from RH to T plus the half-plane:
  that is our 106.25 T^2 (KNOWN).
* Computation (routine, COROLLARIES): Li coefficients are additive over the character factorization of Dedekind zeta, so the
  Dirichlet ranges give lambda_n(zeta_K) > 0 for every abelian K of conductor f <= 125,000 and n <= 26.12 (10^8/f)^2
  (f <= 1000: 41.30e16/f^2; Q(sqrt(-3)) 6.05e16, Q(i) 2.58e16, Q(zeta_5) 1.65e16, Q(zeta_7) 8.4e15); no Dedekind Li range on arXiv.
* Wild: the half-plane alone gives no range for lambda_n or its arithmetic part: zeros with |gamma| >> sqrt(n) have
  |1 - 1/rho|^n near 1, so the high-zero tail needs RH to a height or phase input (KNOWN, R296, R333).
* Adversarial: Zenodo 23488592 says it does not claim a proof (NOT TESTABLE). Zenodo 23977400 (Chung, N6, claimed RH): two more
  scouts; Section 6's step (critical-value sign alternation from B = Lambda A + R, Lambda > 0, |R/A| < 1) is valid in form, and
  Theorem 5.2's use of the predecessor's real roots is the induction hypothesis, not a flaw. First flaw still not located; next
  test: compute the claimed uniform floor omega >= 12.2757 (d >= 9) and the bound |R/A| < 1 for small (d, k) from exact Taylor
  coefficients of Xi.
* Exploration: anti-concentration (2609.17875 treats random or partial-sum Dirichlet polynomials, not zeta), KAM and small
  divisors (Baker bounds on sum k_p log p only make universality effective; no zero location), Fourier restriction (additive
  restriction cannot see the multiplicative phases p^{-it}): no idea.
-/

namespace Fayed.R346
end Fayed.R346

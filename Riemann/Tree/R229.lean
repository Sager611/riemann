/-!
# Round 229 — Gale's conjecture now fails in every dimension n ≥ 3 (discovery logged); hunts clean

* Discovery logged (discoveries/gale-r7-r8.md): first Gale counterexamples in R⁷ and R⁸, so the conjecture fails in every
  n ≥ 3; explicit family with a hand proof for 5 ≤ n ≤ 1000. Novelty check: existence for n ≥ 9 was already known
  (OpenAI's Borsuk counterexamples, Cor. 7.1) and the two-simplex template is a search seed in Voronov's repository;
  REPORT corrected accordingly.
* Verified (others' claims, correct): 2611.05046 (4-qubit separable counterexample, numpy residuals ≤ 1e-15),
  2611.04387 (spt congruences mod 5, 7, 13 on small n, two methods).
* No counterexample: 2611.04297 (Lemma 5.3, n ≤ 4), 2610.39967 (n = 5 in ℓ∞⁵), 2611.05271 (C_n, n ≤ 24). Unverified
  scout report: 2611.05034 Theorem 2.1 lacks an "N large" hypothesis (K₂ edge case; minor). Not testable:
  2610.41781 (asymptotic variance), 2512.24343 (unspecified Monte Carlo).
-/

namespace Fayed.R229
end Fayed.R229

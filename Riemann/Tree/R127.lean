/-!
# Round 127 — quick audits of 8 math.GM / math.CA / math.CV papers claiming proofs of open problems (concurrency 4)

* 2609.33018 Maxell–Wyer, almost-bounded trajectories of Collatz-type ("Hydra") maps on number rings: FALSE STEP
  (l. 3039–3049).
  - It bounds P(vector of run counts ≤ C₀n) ≤ P(sum of the vector ≤ C₀n). For non-negative entries the
    inclusion goes the other way: the sum bound implies the componentwise bound.
  - The next step, from the joint event to max_j P(#_j ≤ C₀n) ≪ x^{−c}, is also reversed.
  - The "head bound" is therefore unproved. Rechecked in the TeX.
* 2610.18895 Zweiman, "proof of Karp's conjecture" (GPT/Claude-assisted): overclaims. The literal conjecture
  fails for PF sequences with endpoint zeros, as the paper's own example (0,0,1,2,1) shows: the root is 0, not
  negative. Its Proposition 2 states the corrected version. Minor.
* 2609.22336 Deng, "proof" of Stieltjes' irreducibility conjecture for Legendre polynomials (open since ~1890):
  - the key step (resultant upper/lower bounds) is unverified;
  - the TeX does not compile (duplicate lemma header);
  - P_n and P_n/x are irreducible for n ≤ 40 (sympy).
  → candidate for a deep audit.
* 2610.04972 Liu–Mao, infinite log-concavity of the Taylor coefficients of ξ: every Turán iterate is positive
  for all n, via the McNamara–Sagan-type closure r ≥ (3+√5)/2, an analytic part for n ≥ 12000 and an interval
  computation below. Iterates m ≤ 7 rechecked with mpmath. No flaw at triage depth; not RH.
* 2609.09214 Tiebekabe (p^x + 6^y = z^{2n}): brute force finds no missing solution. 2609.00039 Laniewski (abc):
  nothing checkable. 2609.28529: withdrawn (404). 2609.22330 (common roots of Legendre polynomials): the scout
  timed out.
-/

namespace Fayed.R127
end Fayed.R127

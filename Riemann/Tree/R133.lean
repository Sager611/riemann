/-!
# Round 133 — quick audits of 8 new arXiv papers (2610.28685–2610.29963)

Scouts (8/8 ran). Findings below are scout claims, NOT yet verified (verification: R134).
* 2610.28688 (counterexample to the Eventual Stability Conjecture): the counterexample (z^6, −1) checks
  out (z^(6^n) + 1 = ∏ Φ_(2^(n+1)·3^j), SymPy for n ≤ 3). Possible slip: classification (i) says
  a_n(K,d,0) = d^n, but z^(d^n) has one irreducible factor (unless a_n counts multiplicity).
* 2610.28768 (Schäffer's conjecture for even k ≤ 42): the theorem as stated (x, y, n ≥ 2) is contradicted
  by 1² + … + 24² = 4900 = 70², an exception the paper itself lists (l. 447–449); completeness rests on
  Magma runs with no tables or code.
* 2610.28839 (Velani's conjecture): l. 1539–1541 state H^(−2c·log log Q) where the preceding bounds give
  only H^(−2c/log log Q): typo or gap? Also occ_q p = u·p^(−t) for u·q^(−t) (l. 930).
* 2610.29018 (no real zeros of real quadratic L(s, χ), q ≤ 10^9 even, 10^10 odd): none found; the
  Turing/Büthe constants recomputed (max f₁t² ≈ 1.321 < 3.2, max f₂t² ≈ 0.340 < 0.35); no certificates.
* 2610.29470 (exceptional Goldbach numbers among squares, E₂(N) ≪ N^(19/24+ε)): (3.11) asserts
  R₁(2m²) > |R₂(2m)|² without proof (positivity needs R₁ > |R₂|); the Markov step is omitted;
  R₂(2m) vs R₂(2m²).
* 2610.29636 (Gauss sphere problem, error O(R^(55/42+ε)), below Heath-Brown's 21/16): none found in the
  fourth-moment lemma or the exponent arithmetic (17/55, 55/42 recomputed). Record claim: second audit.
* 2610.29963 (sharp bound for additive complements): none found; no conflict with Radomskii 2610.14127.
* 2610.28685 (asymptotic for non-Wieferich numbers, conditional on W(x) = O(x/log^α x), α > 2):
  l. 343–347 replace the inner range q ≤ x by q ≤ z without justification.
-/

namespace Fayed.R133
end Fayed.R133

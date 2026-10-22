/-!
# Round 134 — verification of R133's candidate flaws, second audit of the sphere record, listing scans

Scouts (8/8 ran):
* 2610.28768 (Schäffer's conjecture for even k ≤ 42): real, but only a statement slip. The theorem
  (l. 481–485: even 2 ≤ k ≤ 42, x, y, n ≥ 2, no solutions) is false as written: 1² + … + 24² = 4900 = 70²,
  an exception the paper itself lists (l. 447–449); nothing excludes it. Fix: exclude (k,x,y,n) = (2,24,70,2).
* 2610.28685 (non-Wieferich numbers): typo, harmless. Keeping q ≤ x, Brun–Titchmarsh and partial summation
  give Σ_(q ≤ x, q ≡ 1 mod p) 1/q ≪ log₂x/p, the bound used at l. 347.
* 2610.29470 (Goldbach numbers among squares): typos, harmless. (3.10) bounds Σ|R₂(2m²)|, and
  R₁(2m²) > N² = X/2 (l. 244–248), so Markov gives E₂(N) ≪ N^(19/24+ε).
* 2610.28839 (Velani's conjecture): typos, harmless (l. 1540 should read H^(−2c/log log Q), the form used
  at l. 1416; l. 930 should read u·q^(−t)).
* 2610.28688 (eventual stability): not a flaw; a_n counts irreducible factors with multiplicity (l. 92–96).
* 2610.29636 (Gauss sphere problem, R^(55/42+ε)): one display is wrong (S_b(θ)² = Σ e(θh)U_(b,h) must be
  S_b(θ)⁴ = |Σ e(θh)U_(b,h)|², l. 348–365) and the Cauchy step after it works with the fix; the reduction
  (l. 108–245, 984–996) and the exponents check out; the key bound A_δ ≪ T^ε is imported from
  Cilleruelo–Garaev (l. 390–407); the main theorem is unconditional (only Theorem qrh assumes quasi-GRH).
  Next: check the import.
* Listing scans (math.NT; math.GM, CV, CA), next targets: 2610.29790 (Gauss–Fourier conjecture),
  2610.21472 (Smale's mean value constant ≥ 1.10186), 2610.22927 (consecutive gaps between zeta zeros),
  2610.28782 (relative Schinzel–Zassenhaus), 2610.28822 (short Möbius sums over function fields),
  2610.29211 (gradient conjecture at infinity), 2610.16605 (nil-Bohr recurrence).
* Reviewer: the Sato–Tate-in-classes entry is the general form of the sign entry. Merged into one
  DISCOVERIES entry; hypotheses (k ≥ 2, trivial nebentypus) and the discrepancy wording fixed.
-/

namespace Fayed.R134
end Fayed.R134

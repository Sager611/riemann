/-!
# Round 274 — Sun's OEIS conjectures (permanents, determinants, Hamiltonian cycles, residues): hold; one misprint

* Misprint (scout found p = 5; orchestrator checked p < 120): OEIS A372409 Conjecture 2 (Sun, 2024) states det S(p) = -det T(p)
  = 2^((p-5)/2) s(p). The first equality fails for every prime 5 <= p < 120, while det S(p) = 2^((p-5)/2) s(p) and
  det T(p) = -2 det S(p) hold throughout: a factor 2 is missing in the printed statement (Sun reports checks for p < 2000).
* Known: OEIS A321611 (the n = 12 exception is recorded there).
* Hold in tested cases: A227609, A228626, A228917, A228886, A228860 (n <= 8), A262813 (n <= 20,000), A229038, A229141,
  A319311, A319882, A319903, A309012, A320060 (primes below 97-300), A261583(iii) (n <= 10^4); MathOverflow 498447,
  479729, 507228.
* One scout timed out.
-/

namespace Fayed.R274
end Fayed.R274

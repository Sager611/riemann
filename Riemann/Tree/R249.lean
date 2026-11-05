/-!
# Round 249 — Sun's 2013 conjecture on circular permutations with primitive values fails for q = 11, a₀ = 2 (candidate); hunts clean

* Counterexample candidate (scout + orchestrator check): Z.-W. Sun's conjecture (arXiv 1309.1679, Conj. 3.10 of 2013-09-17;
  journal version Conj. 3.8, 2019; restated in 2610.06289, which proves it for q > 18,888,871) says that for every q > 7 and
  every a₀ ∈ F_q some circular arrangement of F_q^* has all a₀ + a_i a_{i+1} primitive. For q = 11, a₀ = 2 the primitive
  roots are {2, 6, 7, 8}, so neighbours satisfy xy ∈ {4, 5, 6}; vertices 2 and 9 both have exactly the neighbours {3, 8},
  forcing the closed 4-cycle 2-3-9-8: no arrangement exists (exhaustive search agrees). All a₀ for primes 13–47 work except
  a few cases the quick search left open (q = 31, 37, 43). Orca "sun38" is settling all prime powers to a bound and checking
  novelty; then reviewer and DISCOVERIES.
* Verified (others' results, correct): 2610.22622 (seven-square example, L₇* = 21), 2610.03730 (both-crossing probability
  ≤ 1/4, n ≤ 6), 2610.24502 (PSL(2,7) orbit data). Not testable: 2610.05568, 2610.33322, 2610.02009 (timed out),
  Zenodo 23765626.
-/

namespace Fayed.R249
end Fayed.R249

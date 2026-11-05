/-!
# Round 266 — Sun's OEIS congruence conjectures at the smallest primes, Guo's q-congruences: all hold; two scout misreadings

* Hold at the smallest allowed primes (exact): Sun's OEIS conjectures A154257, A304081, A154285, A155114, A232443(i-iv),
  A154290, A155860, A154940, A275027, A244973, A277640, A319894, A260667(iii), A236966, A236308, A242425, A229899, A260928,
  A260929; Guo's q-supercongruence conjectures in 1909.10294 (Conj. 1-2) and 2310.15207 (Conj. 5-6).
* Known: A154364 (Sun-Hou conjecture fails at n = 169421772576, stated there); 2610.36620 (the paper's own refutation).
* Scout errors: 2207.13039 Conj. 4.3 at p = 3 (both sides of the equivalence are false, so it holds); Zenodo 23922416
  Prop. 9.1 (pdftotext dropped the conjugation bar in zeta_C(1/2 - it) = conj zeta_C(1/2 + it); the step is right).
* Fasi-Higham diagonal case (Orca fhdiag2): claimed THEOREM for all odd m (computer-assisted; Riemann/Deep/FasiHigham2).
  Orchestrator spot check: m = 205, all c_n > 0 for n <= 2500 (mpmath at 1500 and 2500 digits agree). Adversarial review
  and novelty check running (Orca fhrev).
-/

namespace Fayed.R266
end Fayed.R266

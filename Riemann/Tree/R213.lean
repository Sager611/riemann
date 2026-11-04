/-!
# Round 213 — two scout flaw claims rejected; openai/math spot checks hold; hunts clean

* Rejected scout claims (orchestrator check): (1) Zenodo 23879511 (Collatz): the "(3x+1)/4 clone" maps x ≡ 2 (mod 4)
  to 2x and then divides by 4, i.e. x → x/2 in two steps (2 → 4 → 1), so it is a valid reformulation, not the flaw;
  (2) openai/math "Bounded-Step Walks on Gaussian Primes": its paper.pdf is the right paper (title page checked; the
  scout had read the Liénard PDF).
* Verified (others' claims, correct): 2611.05300 (happy-edge counterexample: 16-point table, 12 vs ≥ 13 flips).
* Spot checks hold: openai/math quintic Liénard (two limit cycles near amplitudes 1 and 2, ε = 0.01, RK4), openai/math
  totient fibers (finite analogue, x ≤ 10⁶), 2610.41569 Lean repository (no sorry/axiom; external results explicit).
* No counterexample: 2611.04345 (x ≤ 2·10⁶), 2602.15905 (Lagarias σ(n) ≤ H_n + e^{H_n} log H_n, n ≤ 10⁶),
  2603.29973 Conj. 2.2(i) (n ≤ 100). Not testable: openai/math planar circle-packing certificate (PDF stream errors).
-/

namespace Fayed.R213
end Fayed.R213

/-!
# Round 334 — node N4 (Robin / Lagarias), round 4: the new 10^(10^38.6) record stands (Polak has no newer version); Lagarias' form follows on the same range; Nicolas rechecked independently to 10^6

* Node: N4 Robin / Lagarias (fourth round; R316 before; robinh's record processed this stretch).
* [check] No 2026 arXiv result adds an unconditional family, a range or a proof (2610.39120, 2609.33794, 2609.26787 SEEN);
  Polak's Zenodo 21808589 has only v1.0.0 (5 Aug 2026, n <= 10^(7.1*10^22)): our 10^(10^38.60) stands.
* Computation (routine): Lagarias' sigma(n) <= H_n + e^{H_n} log H_n holds for 1 <= n <= 5040 (exact sigma; equality only at
  n = 1), and for n >= 5041 it follows from Robin's inequality (H_n >= log n + gamma gives e^{H_n} log H_n >= e^gamma n log
  log n); hence Lagarias' criterion is verified to 10^(10^38.60). TREE N4 updated.
* Adversarial (independent recheck of robinh part 1): Nicolas' inequality for every primorial with p_k <= 10^6 (mpmath,
  dps 30): minimum margin 1.487e-4 at p_k = 999983, agreeing with a direct evaluation to 18 digits: HOLDS.
* Harvest: 2610.39120 (CA2 numbers) leaves open q:bound (TeX 278-285), q:CA2 (300-303) and q:higher (382-390); its Q1
  was tested in R206 (holds for a + 1 <= 33, smallest margin 1.2e-10).
* Exploration: graph Riemann-Roch L-functions (2608.29981) have a functional equation but no RH-type positivity (R290,
  R328); random multiplicative functions (2608.21354) need independent prime values, no transfer to mu; Selberg-CLT rates
  (2308.09679, 2212.01411, 2507.04150) give no N7/N8 control: no idea.
-/

namespace Fayed.R334
end Fayed.R334

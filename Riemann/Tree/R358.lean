/-!
# Round 358 — node N11 (GRH families and Siegel zeros), round 6: the explicit Littlewood record stands; its Dusart inputs and Lemma 6 rechecked; observed minima about 6.4 times the proved constant

* Node: N11 GRH families and Siegel zeros (sixth round; R342 before).
* [check] state of the art: our explicit Littlewood bound |L(1, chi)| >= 0.11/max(3, log log q) (littlewood78) stands; the Oct-Nov
  listings and Zenodo give only 2611.01451 (GRH-conditional first moment), 2611.03356 (Beurling primes) and 2610.15483 (gap 3.0e-3,
  known): HOLDS.
* Adversarial (Dusart inputs): arXiv:1002.0442 TeX lines 660-670 (Thm 6.11: sum_{p <= x} 1/p <= log log x + B + 1/(10 log^2 x)
  + 4/(15 log^3 x) for x >= 10372) and 692-702 (Thm 6.12: sum_{p <= x} (log p)/p <= log x + E + 0.2/log x + 0.2/log^2 x for
  x >= 2974) match PROOF.md: HOLDS.
* Adversarial (Lemma 6, pairs): the worst symmetric configuration is the real pair (B, 1 - B), ratio 0.99384 at (L, h) = (8, 1),
  as PROOF.md states; a scout's FLAW used a zero set without the symmetry rho -> 1 - conj(rho) that every Dirichlet L-function
  has (false positive): HOLDS.
* Computation (truth gap): min |L(1, chi)| max(3, log log q) is 0.70908 over the 30390 fundamental discriminants -d >= -10^5
  (d = 85507, h = 22) and 0.73821 over all non-principal chi mod primes q <= 1000 (chi_{-163}), about 6.4 times the proved 0.11:
  HOLDS.
* Harvest: 2610.31402 asks (l. 982) for effective Brauer-Siegel consequences uniform in the base field and (l. 1751-1755) for a
  numerical c(k) for fixed k (candidate since R330).
* Exploration: monomer-dimer models (Heilmann-Lieb matching polynomials have integer coefficients of one parity; no weighted Xi or
  Jensen representation found), CMV matrices and OPUC (the Hermite-Biehler normalization of Xi gives a Schur function and
  Verblunsky coefficients, but no canonical Caratheodory positivity: zeta'/zeta has poles), lattice gauge theory and Witten /
  Mordell-Tornheim zeta functions (self-inversive Bernoulli polynomials have unimodular zeros, 2609.07583; Witten zeta zeros at
  s = -2 only conjectured, 1304.3543): no new route.
-/

namespace Fayed.R358
end Fayed.R358

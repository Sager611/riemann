/-!
# Round 386 — node N9 (moments, Lindelof, subconvexity), round 9: records stand (Zenodo, HAL nothing after 4 Nov); 2609.00537's explicit B-process (Cor. 'A van der Corput B-process') has a smaller constant than Patel-Yang's lemma, but in our sub-Weyl certificate it touches only negligible E-terms, so 54.6 t^{27/164} stands; de Faveri-Pandey's moment exponents recomputed exactly; explicit Lindelof beyond the half-plane gives nothing quoted; Selberg integrals, Mirzakhani recursions and flat Littlewood polynomials give no route

* Node: N9 moments, Lindelof, subconvexity (ninth round; R376 before).
* [check] Zenodo ("Lindelof", "subconvexity", explicit critical-line bounds; newest relevant item 1 Nov) and HAL (newest
  2025-10-09): no new explicit bound on |zeta(1/2+it)|, no mu(1/2) < 13/84, no new moment exponent: HOLDS.
* Computation/harvest: arXiv 2609.00537 (TeX l. 297-314, corollary 'A van der Corput B-process') replaces Patel-Yang's B-process
  constant 3 pi^-1 log(beta - alpha + 2) + 4.0 by 2 pi^-1 log(f'(a) - f'(b)) + 2.686/sqrt(lambda_2) + ... + D_f, a saving of about
  pi^-1 log D + 2.6 in that coefficient; the paper says it could refine Patel-Yang but computes no constant. Our sub-Weyl
  certificate (Riemann/Deep/SubWeyl/REPORT.md, 'Further gains') already found that this B-process only touches E-terms,
  which are negligible at the binding range log t in [648, 650]; the ABA^3B blocks (Weyl differencing, Patel-Yang's B-process
  inside them, Yang's A_5, B_5) carry 94% of the bound: no gain from the substitution (the scout's candidate was ours, seen).
* Adversarial: de Faveri-Pandey 2610.01035 (v1 only, no comment): 37/30 = 5/4 - 1/60 (TeX l. 1050-1052) and the twelfth-moment
  exponent 2 - 520/445121 (l. 1011, 1027-1032) follow exactly from the paper's own endpoint parameters: HOLDS.
* Harvest: 2606.12376 (TeX l. 162-178, 258-263: A_0, eta, the implied constants and t_0 unspecified) and 2609.15862 (l. 1216:
  a Selberg zero-density constant C) need explicit inputs that do not exist yet: low.
* Wild: an explicit Lindelof-type bound beyond the half-plane (Borel-Caratheodory from Re s > B) is routine, already implicit
  in mu(sigma) = 0 for sigma > B (COROLLARIES, R280, R343), and improves nothing quoted in N9 or N10: DEAD.
* Exploration: Selberg integrals and Jack polynomials (Keating-Snaith moments are exact Pochhammer ratios but finite-N
  positivity averages away the zero locations), Mirzakhani's volume recursion (nonnegative coefficients, no identity linking it
  to xi or the Li coefficients), flat Littlewood and Rudin-Shapiro polynomials (+-1 automatic coefficients, no Euler product):
  no route.
-/

namespace Fayed.R386
end Fayed.R386

/-!
# Round 341 — node N12 (Selberg class, function-field analogues), round 5: classification and density records stand; conductors 2, 3, 4 open (needs a linear-twist argument); the finite Rodriguez-Villegas map is not Laguerre-Polya preserving

* Node: N12 function-field and Selberg-class analogues (fifth round; R329 before).
* [check] arXiv (Oct-Nov 2026): no classification beyond degree 2 conductor 1 (2503.02492); 2610.36885 gives Omega(+-) bounds for
  degree-2 Hardy functions Z_F (any level, Maass forms unconditionally): oscillation, not a classification or zero-location record.
  Density: Cossaboom 2609.27715 (seen R129) N_pi(sigma, T) << T^{(5/2)(1 - sigma) + o(1)} for GL_2/Q on 1/2 <= sigma <= 1; the
  density hypothesis for degree <= 2 still only for sigma >= 0.9 (2408.00617). Both now in TREE N12.
* Harvest: 2503.02492 (Introduction l.165) leaves the degree-2 classification open for conductors q = 2, 3, 4; the paper names no
  finite computation, only a missing linear-twist argument (not computation-settled).
* Function field: Goss zeta (2312.01264, Thm t:Riemann-hypothesis): the Newton-polygon zero pattern holds for ordinary curves;
  non-ordinary curves and Goss L-functions open; the mechanism is non-archimedean (unique minimal valuations), no transfer to N6.
* Adversarial: the newest class-wide GRH claim is again Aksman, Zenodo 23839754 (refuted R338, eqs. (78)-(79)); no new claim.
* Computation (R338 follow-up): the kernel -> Z_E map of 2609.25564 does not preserve the Laguerre-Polya class: a Gaussian kernel
  (transform e^{-t^2/4}, no zeros) gives non-real zeros of Z_E from E = 16, a Gaussian-cos kernel (real zeros) from E = 40, while
  Riemann's kernel stays real-rooted to E = 90 (lambda = 0 and 0.5); so RH does not give 'Z_E real-rooted for all E' by
  LP-preservation ([RVApprox](Riemann/Deep/RVApprox/README.md)).
* Exploration: Krein strings / canonical systems (1204.1827 Thm 4: explicit positive Hamiltonian only for xi(1/2 + omega - iz),
  omega > 1; omega -> 0 is RH-equivalent), quantum-chaos semiclassics (2406.01828, 2406.12852, 2505.21192: heuristic or
  non-normalizable), Helton-Vinnikov determinantal representations (need strict hyperbolicity; finitely many real-rooted slices do
  not propagate, example X^2 + YZ - 10Z^2): no idea. explore7 (Orca, monotonicity formulas for the heat flow): no idea (NOGO).
-/

namespace Fayed.R341
end Fayed.R341

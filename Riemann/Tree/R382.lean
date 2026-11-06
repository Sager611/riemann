/-!
# Round 382 — node N2 (de Bruijn-Newman), round 9: Lambda <= 0.107559925318 stands (2611.06245 claims only 0.156821); the uniform Dirichlet bound 0.2788 is already recorded; Raymond's tube check is conservative (tighter arcs could sharpen Lambda, candidate, low); a Zenodo 'collision shield' RH proof is refuted; polynomial inequalities, Mosco convergence and spectral pollution give no route

* Node: N2 de Bruijn-Newman constant (ninth round; R369 before).
* [check] Zenodo, HAL, listings: Raymond's Zenodo 23752103 is still v1.0 (29 Oct); Zenodo 23842668 gives no bound on Lambda;
  arXiv 2611.06245 (4 Nov, seen) claims only Lambda < 0.156821; HAL: none: HOLDS.
* Computation: the uniform bound Lambda_chi <= 2(B - 1/2)^2 = 0.278828185983893 for every primitive Dirichlet chi (B the
  Lean-checked edge) is already in COROLLARIES (R324); 2610.27343 defines H_t^L(z) = int e^{tu^2 + izu} Phi_L(u) du, so the strip
  lemma applies with this scaling, and its uniform 953/3200 is superseded: KNOWN.
* Harvest: Raymond (Zenodo 23752103, sec. 5.3, p. 34) calls the retained tube check conservative: recomputing its 2,560 closed
  arcs with tighter bounds could lower Lambda (candidate dbntube, low); sec. 1 (p. 7) says the gains come from the source and
  geometry inequalities, so a better half-plane or source bound also moves it.
* Adversarial: Zenodo 23863923 (Franquero, Part IV, 1 Nov) reduces RH to the absence of double roots at lambda = 0; non-real
  zeros need no collision at lambda = 0 (z^4 + 1 has no repeated root and no real root), so the reduction fails: FLAW
  (REFUTATIONS).
* Wild: certify Raymond's barrier at t = 0.10 from the zeros up to 3e12 plus an explicit zero-density majorant for the omitted
  terms near sigma = 0.8734; a density bound cannot keep a single zero away from the barrier region, so the needed uniform tail
  inequality has no source: low.
* Exploration: polynomial inequalities (Bernstein, Turan, Erdos-Lax bound norms under zero-location hypotheses; no converse
  real-rootedness criterion for Jensen polynomials), Gamma/Mosco convergence (on the expanding windows [1/lambda, lambda] it
  does not prevent dilation or translation escape, so it cannot give CCM's ground-state or Fourier convergence), spectral
  pollution (second-order relative spectra enclose eigenvalues at finite lambda, where CCM's spectra are already real; no
  lambda -> infinity statement): no route.
-/

namespace Fayed.R382
end Fayed.R382

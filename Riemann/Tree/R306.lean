/-!
# Round 306 — node N2 (de Bruijn-Newman), third round: HB gates rechecked; model estimates 0.13925 / 0.13897 with Liu's newer half-planes; Dirichlet analogue Lambda_chi <= 953/3200 found

* Node: N2 de Bruijn-Newman constant (third round).
* State of the art: Lambda <= 0.13963 (Hakimoglu-Brown 2610.12803) stands as far as checked (the listings scout timed
  out). L-function analogue: 2610.27343 (Theorem, TeX lines 214-231) proves Lambda_chi <= 953/3200 = 0.2978125 for
  primitive real chi with 3 <= q <= 400, unconditionally; added to TREE.md N2. Liu's Hecke half-planes plus Platt's
  GRH heights cannot give HB-type bounds uniformly in q (half-plane alone gives about D^2/2 = 0.281; conductor-uniform
  barrier constants and verified zeros for all q are missing).
* Computation (model, not certificates): under R293's re-optimisation model Lambda(0.874058504) = 0.13925 (Liu
  Hexagon v3, audited sound given inputs) and Lambda(B_2 = 0.8733819666) = 0.13897 (Liu 2611.00082, audit running);
  slope check dLambda/dtheta = 3 dLambda/dq_0 = 0.406 at 7/8.
* Adversarial: Bernstein gates of 2610.12803 row 0 (l = 0, a = b = 0: beta = 62783/37210; l = 7, a = b = 4: 3/10)
  exceed 13/10000 exactly: HOLDS. The Rodgers-Tao dynamics scout timed out.
* Exploration: capacity (Gurvits/Brändén-Huh 2007.08390) certifies coefficient inequalities only after stability is
  known (N6, weak); Fourier interpolation with zeta zeros (Bondarenko-Radchenko-Seip 2005.02996) is unconditional and
  gives no zero location (N1/N13, weak); AAK/Hankel theory would need a Mellin intertwining that turns the
  Nyman-Beurling span into a model-space problem (N5, idea).
-/

namespace Fayed.R306
end Fayed.R306

/-!
# Round 364 — node N6 (Jensen-polynomial hyperbolicity), round 7: records stand; RH height 3e12 is still the rigorous record; Chung's Rev. 3.1 keeps the false Lemma 4.0(2); the Laguerre range of 2611.02944 is tied to the height

* Node: N6 Jensen-polynomial hyperbolicity (seventh round; R345 before).
* [check] RH height (cross-node: the N2-N6 records use it): Platt-Trudgian 2004.09765 Thm 1 (T = 3,000,175,332,800, refereed)
  is still the rigorous record; 2408.00187 only develops a verification method: HOLDS.
* [check] listings and Zenodo: only Chung's Rev. 3.1 (Zenodo 23987495, 5 Nov 2026, claims RH via all J^{d,k} hyperbolic) is new;
  no new hyperbolicity range: HOLDS.
* Adversarial: Chung's Rev. 3.1 keeps Lemma 4.0(2) (turnover index r*(z_m) = m, p. 20) and uses it at m = d - 1 and d - 2 in
  Thm 4.1; our ChungAudit check still fails at (d, k) = (9, 0): z_8 = 10.6458328564 has y*_turn = 5.4267842 < y_7 = 5.7444601,
  so r* = 7, not 8: FLAW persists (REFUTATIONS updated).
* Computation/harvest: 2611.02944 (generalized Laguerre inequalities L_n[Xi](x) >= 0 for every n on |x| <= T - 1/2, TeX l. 416-428)
  rests on RH to T = 3e12; the half-plane does not extend it; the authors call the global Laguerre inequality L_1[Xi] >= 0
  (Csordas-Varga, Problem 2) open and much stronger (TeX l. 343-346).
* Wild (de Bruijn-Newman flow): GORTTW 1910.01227 Thm 1.1 gives J^{d,n} hyperbolic only for n >= c e^{d/2} with ineffective c;
  no monotonicity in t or thresholds t_d -> 0 for the Jensen polynomials of H_t are known: NO IDEA.
* Exploration: the Minkowski question-mark function and Stern-Brocot statistics (their Fourier decay concerns continued-fraction
  cylinders, not the additive Farey discrepancy of Franel-Landau; R048, R058): no route; Witt vectors / lambda-rings and quantum
  walks: not tested (arXiv API throttled, HTTP 429).
-/

namespace Fayed.R364
end Fayed.R364

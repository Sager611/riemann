/-!
# Round 418 — node N8 (critical-line proportion), round 10: inputs of today's GL(2) discovery authenticated (Andersen-Thorner's length nu < 1/4 - theta/2 and kappa bound, Bernard's mean value and Q class, 2607.00282's sup definition and corollaries, 2508.11108's Question 5 and table); an independent recomputation agrees to 1.8e-14; one scout FLAW rejected; MRI pulse design, Bose-Einstein condensation and foams give no route

* Node: N8 critical-line proportion (tenth round; R395 before). Round type: authenticate the quoted inputs of today's discovery
  (degree2mol: at least 7.005% of the zeros of every GL(2) newform L-function on the line; discoveries/gl2-critical-line.md).
* Andersen-Thorner 2004.03581: Thm thm:nu (0 < nu < 1/4 - theta/2, every holomorphic or Hecke-Maass newform on Gamma_0(N), l.138,
  l.168-173, l.229), kappa_f >= 1 - (1/R) log c with real P, P(0) = 0, P(1) = 1 (l.194, l.214-218), M = 14 search with no optimality
  claim (l.1002-1012): HOLDS. The scout's FLAW ('7.005% absent from the source') is a misread: 6.96% / 8.96% are the 'before'.
* Bernard 1404.6432: Thm th:mm2 and the proposition (l.770-796; ranges 5/27 at theta = 7/64 and 1/4 at theta = 0); our odd-Chebyshev
  Q lies in his class (formepolynome, l.785-789; i^(k+1) real for odd k): HOLDS.
* 2607.00282: real P, Q with Conrey's conditions (l.545-547), kappa(theta) defined as a sup (l.565), Cor. fullresultpgl3per and
  simplpergl2 (l.587-623; RatioGL3.tex l.184-185), the 2 theta/3 step cites 2508.11108 (l.604-608): HOLDS.
* 2508.11108: Question 5 (l.964) verbatim; the paper proves kappa(theta) > 2 theta/3 only for theta < min(theta_0, theta_1), thresholds
  it defines but does not evaluate (l.307, l.774-790; our page said 'non-explicit', reworded); table rows kappa(2/3) = 0.364...,
  kappa(1/2) = 0.334... (l.823-830) and l.849 (8.54% under Ramanujan; AT21 0.0696 / 0.0896) as quoted: HOLDS.
* Computation: independent mpmath evaluation of c = (1 + e^{2R}Q(1)^2)/2 + (C/xi)W0 + xi B W1 from params.json gives
  kappa = 0.07005623901882 (certificate 0.0700562390188, difference 1.8e-14) and 0.08971017893360 under Ramanujan: HOLDS.
* Exploration: MRI pulse design (Shinnar-Le Roux polynomial pairs and Fejer-Riesz factorisation, already NOGO.md l.554: minimum
  phase does not exclude Davenport-Heilbronn), Bose-Einstein condensation (polylog criticality and Lee-Yang fugacity zeros live in the
  fugacity variable, already near NOGO l.67, 171, 461) and foams (Plateau, von Neumann-Mullins: local topological balance laws with
  no arithmetic invariant): no route.
-/

namespace Fayed.R418
end Fayed.R418

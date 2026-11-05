/-!
# Round 312 — node N5 (Nyman-Beurling), fourth round: open problems next to our L^p theorem - Dimitrov-Oliveira's P1 (Dirichlet L-functions) looks like a routine extension; Noor's H^p converse and weighted alpha are hard

* Node: N5 (fourth round; harvesting open questions next to the R311 discovery).
* Dirichlet L-functions: Dimitrov-Oliveira (1608.07887, Thm GBNew/GBLp; their open problem P1, TeX lines 241-251) ask
  for a true L^p natural-number endpoint theorem; the natural statement is 1_(0,1) in the L^p-closure of
  {kappa_{chi,p}(1/(kx))} (principal chi: sum b_k k^{-alpha} = 0) iff L(s, chi) != 0 on Re s > 1/p. nblp's proof should
  adapt (mu(k) -> mu(k) chi(k), completed xi_chi with gamma factor Gamma((s + alpha_chi)/2), principal pole,
  q-dependent constants), giving it unconditionally for p <= 1/0.87338196661320752 = 1.14497 by the Lean-checked
  half-plane for all Dirichlet L-functions: next Orca target (IDEA, routine extension).
* Noor's Hardy-space model (Ghosh-Kremnizer-Noor-Santos 2206.00434, Theorem 'Half-plane 1 < p < 2'): the open converse
  (zeta != 0 on Re s > 1/p implies 1 in the H^p-closure of span{h_k}, h_k(z) = log((1 + ... + z^{k-1})/k)/(k(1 - z)))
  does not follow from our theorem: their L^2 identification H^2 = l^2_omega = M in L^2(0,1) has no L^p analogue
  (HARD). Gallardo-Gutierrez-Seco (2206.02654): other weights alpha need weighted L^p(x^{-alpha-2} dx) (HARD).
  2606.16097: the 1/2 <= alpha <= 1 question is an H^2 (p = 2) problem out of reach (HARD).
* p > 2: necessity holds for all 1 < p < infinity (PROOF.md Prop. N; GGS Thm A), so 1 is NOT in the L^p-closure for
  every p > 2 because of the zero at 1/2 + 14.1347i: KNOWN-HERE.
* Exploration: interpolating the unconditional L^{1.14497} approximants with RH-conditional L^2 decay gives only
  logarithmic rates in between (N5, idea); zero-free alone does not give cyclicity in Dirichlet-series Hardy spaces
  (1908.03113: explicit zero-free non-cyclic F) (N5, barrier); open problems weaker than RH (2509.20335, 2602.04022:
  subconvexity below 13/84, pair correlation) remain out of reach of the half-plane (N9, none).
-/

namespace Fayed.R312
end Fayed.R312

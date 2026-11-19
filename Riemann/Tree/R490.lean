/-!
# Round 490 — node N8 (critical-line proportions and simple zeros), round 31, [check] of the load-bearing input of explore52's NPA lead: Hughes–Rudnick's linear statistics of ζ's zeros are unconditional and their diagonal argument gives an unconditional prime side for mixed Rudnick–Sarnak-range words, Alpöge–Furman's zero side holds in the complex frame, but the published triple correlations (Hejhal, Rudnick–Sarnak) are RH-conditional: the input is an open obligation; explore52's Claims A/B hold in their smallest case N = 3

* Node: N8 (31st round; R481 before). Round type: [check] of the load-bearing input of explore52 (ExploreNPA: unconditional evaluation of
  Rudnick–Sarnak-range triple and quartic words in the complex-zero frame) and a smallest-case falsification of its Claims A/B. Texts
  pre-extracted in .scratch/r490/txt (2608.13637, math/0208220, math/0206289, math/0208230). Load ≈ 7: all eight scouts answered.
* [check] Hughes–Rudnick, arXiv math/0208220 (linear statistics of the zeros of ζ): unconditional (l.129–131, 'will not assume the Riemann
  Hypothesis'; f̂ smooth with compact support, so f(γ) is defined for complex γ); the m-th moment of N_f is mock-Gaussian for
  supp f̂ ⊂ (−2a/m, 2a/m), H = T^a (theorem l.151–152). HOLDS.
* [check] the diagonal step for mixed words: HR's off-diagonal argument (l.437–443: H = T^a ≪ N ≤ √(MN) ≤ T^{(Σα_j − mε)/2}) goes through
  for Dirichlet polynomials of unequal lengths T^{α_j} whenever Σ_j α_j < 2a: the prime side of mixed words is unconditional exactly in the
  Rudnick–Sarnak simplex. HOLDS.
* [check] Alpöge–Furman 2608.13637: item (e) (l.591, stated 'informally') gives the prime side of tr G̃^k unconditionally for X^k ≤ T^{2−ε};
  G̃ is built from all zeros with off-line pairs ρ, 1 − ρ̄ (l.226–234, l.312–318), so the zero side holds in the complex frame for every k;
  but the Extensions (l.603) take the third trace 'under the Riemann hypothesis [Hej94, RS96]': the published triple-correlation
  evaluations are RH-conditional as cited. Net: explore52's input follows from HR's diagonal argument plus AF's complex frame, but no
  published theorem states it for mixed words up to the RS boundary Σ|ξ| = 2: an open obligation, not a citation.
* [check] explore52's 'AF's bandwidth-one ceiling 0.6818287 is this analogue's P problem at N = 256': AF's p₀ = 1 − a_N is the optimum of a
  finite 256-periodic programme whose form factor matches CUE to |256 S(j) − j| ≤ 3·10⁻⁴⁰ (l.593–600), plus the edge term 0.824|r(1)| and
  the sampled-to-continuous step (EnclOK, not kernel-checked); a scout's FLAW on the wording concerns these transfer terms only: the
  identification holds for the finite programme.
* [falsify] explore52's Claims A/B at N = 3 (M = 96, its 12-point depth grid; orchestrator, npa_cg.py on a copy, ≤ 85 MB): P 0.848059;
  P3 off-line without rays 0.868199 (= the RH world), with rays 0.848374 (1.6 % of the triple gain kept); L4 off-line with rays 0.868199
  (all kept). Both claims hold in the smallest case (the report's N = 3 line was without rays). Reruns of verify_n4.py identical (primal
  0.818792412597, residual 1.16·10⁻⁵³; dual bound 0.8429878883).
* Zenodo 19 Nov leftovers: Hackbarth 24359794 (signed areas of termwise curves of normalized η terms; no statement about zeros, off-tree);
  Qin 24364735 (claimed proof of binary Goldbach by the circle method, ~10,000 lines; off-tree, not refereed here).
* Exploration (tool transfers, tested): SOS pseudo-expectations (Grigoriev's knapsack bound; no polynomial system ties zero ordinates to
  Boolean variables), phase retrieval / homometric sets (Rosenblatt–Seymour: autocorrelation is invariant under zero flips, phase-blind),
  Marty's normal-families criterion → N9 (finite samples cannot decide normality; rescaling R12.2 and Bloch approximants R376 tried):
  no transfer.
-/

namespace Fayed.R490
end Fayed.R490

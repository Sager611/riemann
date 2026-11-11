/-!
# Round 431 — node N2 (de Bruijn–Newman), round 11, falsification of claims: Polson–Sokolov's 'Computing the Riemann Hypothesis' (Zenodo 24086101/2, v2.2.0, 10 Nov 2026) restates RH as GGC (Thorin-measure) positivity of 1/G, Prop. 2.1.1 with both directions proved: the known Laguerre–Pólya/Bondesson criterion of arXiv 1804.10043; its computations claim no record (Λ ≤ 0.2 and 0.22 repeated, fits heuristic); a moved zero pair gives complex Thorin atoms; Zenodo 24053448 reproves Λ ≥ 0; Polymath15 states its method's limit; records stand; impedance tomography, state estimation and weaving transfer no tool

* Node: N2 de Bruijn–Newman (eleventh round; R398 before). Round type: falsify new claims (AGENTS.md, 11 Nov).
* Polson–Sokolov, Computing the Riemann Hypothesis I (Zenodo 24086101, v2.2.0, 10 Nov 2026): Proposition 2.1.1 (printed p.18):
  1/G is a generalized gamma convolution (positive Thorin measure) iff RH, both directions proved (pole/Stieltjes argument and a
  gamma-product converse): KNOWN, the Laguerre–Pólya/Bondesson GGC criterion already in arXiv 1804.10043 (Thm 1, Thm 4.2).
  Test: with the first 20 zeros, λ₁ = γ₁² = 199.790454832386859 is a positive atom; moving γ₁ → γ₁ + 0.1i gives the complex
  atoms 199.780454832386859 ± 2.826945028346939i, so positivity fails as it must (RH-sensitive, equivalent to real zeros).
* Volume II (Zenodo 24086102): §4.1 p.48, Table 4.8 p.69 and §5.7.2 pp.120–121 repeat Λ ≤ 0.22 and 0.2; the fitted 0.05 and
  0.01 are explicitly heuristic; the text says finite computation 'cannot take the step'; §21.4 gives only a conditional
  Λ ≤ 0.28118553: no record, no proof claim: HOLDS.
* State of the art: Zenodo 24053448 (8 Nov 2026) Theorem 1 proves Λ ≥ 0 (Rodgers–Tao's theorem, certified); no newer upper
  bound on arXiv, Zenodo or HAL: Λ ≤ 0.107559925318 stands.
* Key lemma of the method: Polymath15 (arXiv 1904.12438, §1 l.165): 't₀ = 0.2, y₀ = 0.2 are close to the limit of our ability
  to numerically verify hypothesis (ii)'; its §6: going beyond Λ = O(1/log T) needs a major advance on RH: KNOWN barrier.
* Exploration (tool transfers, tested): Alessandrini's identity needs a Dirichlet-to-Neumann realisation (no transfer;
  impedance tomography already met in R392); observability Gramians of e^{iγt} modes are positive definite for distinct zeros
  (det = 0.9925 for γ₁, γ₂, T = 1: no bound); Seifert/Alexander invariants of braids and weaves (trefoil t⁴ − t² + 1) carry no
  ζ structure: no transfer.
-/

namespace Fayed.R431
end Fayed.R431

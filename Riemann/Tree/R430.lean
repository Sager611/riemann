/-!
# Round 430 — node N3 (Li/Keiper), round 10, falsification and rate questions: Fabbian's range authenticated (Zenodo 22839423, latest version 19 Sep 2026, Theorem 1: λ_n > 0 for n ≤ 693,318,936,098,417,738,514,687,162 = 77.03 T²); under RH the error of the trend is O(√n log n) (arXiv 1506.01755 Thm 4.1) and no oscillatory main term is proved; no new explicit Li claim to falsify; the Polson–Sokolov Thorin-measure restatement of RH (Zenodo 24086101/2, unseen) noted as a falsification candidate; archaeological seriation, paleoclimatology and cryo-EM transfer no tool

* Node: N3 Li/Keiper (tenth round; R400 before). Round type: falsification and rate questions (AGENTS.md, 11 Nov).
* [check] Fabbian, Zenodo 22839423: latest version 19 Sep 2026 (no newer one), Theorem 1: λ_n > 0 for every
  n ≤ 693,318,936,098,417,738,514,687,162 ≈ 6.93·10^26 ≈ 77.03 T² from RH to height T ≈ 3·10^12: HOLDS ([check] removed).
* Rate: under RH λ_n = (n/2)(log n − 1 − log 2π + γ) + R(n) with R(n) = O(√n log n) (arXiv 1506.01755, Thm 4.1);
  no √n-scale oscillatory main term is proved: OPEN; Johansson's certified λ_100000 (arXiv 1309.2877 §5) agrees with
  Keiper's trend; Maślanka's and Coffey's oscillation formulas were never tested against the n ≤ 10^5 data (a computable
  falsification target).
* Falsify: no new explicit claim on Li/Keiper coefficients in the October–November listings (0 title hits) or on Zenodo
  beyond known items (the orchestrator ran lst.py after the scout returned nothing); Dirichlet: 1507.03431 Prop. 1 is
  conditional on RH to height T (already R296); λ_1(χ_{−3}) = 0.0566149849 and λ_1(χ_{−4}) = 0.0777839900 recomputed
  from L′/L(1, χ): positive.
* Harvest: 2204.01036 l.581–582 (better-conditioned variants of its Li-type determinant) and 2206.03682 l.522–527 (extension
  to automorphic and Selberg-class L-functions): weak.
* New unseen Zenodo items: Polson–Sokolov, "Computing the Riemann Hypothesis" I–II (24086101, 24086102, v2.2.0, 10 Nov
  2026): RH restated as positivity of the Thorin measure of a ratio of Ξ-values; ball-arithmetic computations of the measure,
  the de Bruijn–Newman flow and a Thorin-form barrier method: candidate for an N2/N13 falsification round (does the
  criterion distinguish Davenport–Heilbronn-type counterexamples?).
* Exploration (tool transfers, tested): Fiedler-vector seriation of threshold matrices of zero ordinates (monotone, as for
  any ordered points: no bound), isotonic age–depth inversion (zero gaps are not monotone; paleoclimatology also tried in
  R410) and cryo-EM angular synchronization (needs SO(3) projection data): no transfer.
-/

namespace Fayed.R430
end Fayed.R430

/-!
# Round 436 — node N3 (Li/Keiper), round 11, falsification of explicit Li-coefficient formulas against λ_n for n ≤ 500 recomputed by FFT: Maślanka's and Arias de Reyna's trend formulas leave residuals 0.6–3.4, the oscillating part (a scout 'FLAW' compared λ_n with Maślanka's trend and had no TeX); Coffey's S_γ not independently computable; no stated monotonicity or convexity conjecture, and λ_n is not convex from n = 118; Fabbian's range stands; reaction networks, ecology and phylogenetics transfer no tool

* Node: N3 Li/Keiper (eleventh round; R430 before). Round type: falsify explicit formulas (the open liosc item of R430).
* Data (orchestrator): λ_n = n·[z^n] log ξ(1/(1 − z)) by FFT on |z| = 0.9 and 0.85 (2048 points, 70 digits; the images
  of both disks have |Im s| < 4.8, so ξ has no zeros there unconditionally); both radii agree to 5·10⁻³⁵ for n ≤ 500;
  λ_1 = 1 + γ/2 − (1/2) log 4π exactly, λ_100 = 118.603775376791, λ_500 = 991.900092992245.
* [falsify] Maślanka (math/0402168, eqs. 2.7–2.8 as quoted by the scout from the abstract; the e-print downloads
  returned error pages): his trend λ̄_n = (1/2)(1 + n log n) + (1/2)(γ − 1 − log 2π)n + 1/4 − 1/(24n) + … leaves
  λ_n − λ̄_n = 1.50, 0.63, 2.14, 2.48, 1.41, 2.66 at n = 50, 100, …, 500; these are the oscillating part, not a flaw
  (Arias de Reyna's trend residuals below minus the constant 3/4 agree): NOT TESTABLE (no flaw).
* [falsify] Arias de Reyna (2204.01036, eqs. rf–rt): λ_n = (n/2)(log n + γ − log 2π − 1) + o(n) under RH, the
  off-line sum empty for the first 60 zeros; residuals 1.378, 2.890, 3.234, 3.413 at n = 100, 200, 300, 500: HOLDS.
* [falsify] Coffey (0706.0343, Prop. 2 / Cor. 16: S₂ = S_γ + S_Λ, |S_γ(n) + n| = O(n^{1/4})): S_γ needs the Stieltjes
  constants separately, math-ph/0505052's e-print unavailable: NOT TESTABLE.
* [falsify] stated shape conjectures: 1506.01755 Thm 4.1 and math/0404394 Rem. 4.2 give only O(√n log n) with unspecified
  constants; no monotonicity or convexity conjecture found. Orchestrator recount on n ≤ 500: λ_n increases (smallest step
  λ_2 − λ_1 = 0.0693), λ_n/n increases, λ_n is not convex (first negative second difference at n = 118; 181 of 498
  negative): NOT TESTABLE (nothing stated).
* [check] state of the art since R430: Fabbian (Zenodo 22839423, Theorem 1, 77.03 T²) is still the latest version; Zenodo
  24086101/2 restate the Thorin framework, 23370435 disclaims a proof; HAL nothing after 2023: HOLDS.
* Exploration (tool transfers, tested): Horn–Jackson complex balancing (acts on positive concentrations, Weil positivity on
  signed test functions), May's stability criterion (sees only β > 0), the four-point condition on ordinate distances
  (phylogenetics tried in R373): no transfer.
-/

namespace Fayed.R436
end Fayed.R436

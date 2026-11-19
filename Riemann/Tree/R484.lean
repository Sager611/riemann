/-!
# Round 484 — node N7 (zero-free regions and zero density), round 19, [check]/[falsify] of Fodge's classical-shape records (Zenodo 24325329, 24354828): the exact Stechkin weight (sharp at √5), the ten-shift pair kernel, the Dirichlet chain and the ζ certificate hold: R₀ = 4.6892 for ζ, 4.9774 for Dirichlet L (q ≤ 400 000), Dedekind 11.0416/8.6662/−18.7991; they add to the Lean half-planes only at t < 5.39, q ∈ {3, 4} and tiny fields; Nevanlinna–Pick, Fredholm determinants and Turán–Nazarov transfer no tool

* Node: N7 (19th round; R478 before). Round type: [check]/[falsify] of fresh certified Zenodo record claims next to the node line (lessons 147, 201):
  Fodge, 'Exact Stechkin weights and improved zero-free regions for the Riemann zeta-function and Dirichlet L-functions' (Zenodo 24325329,
  17 Nov 2026): ζ(σ + it) ≠ 0 for σ ≥ 1 − 1/(4.6892 log t), t ≥ 3 (Bellotti–Trudgian–Yang 4.896); Kadiri's region for Dirichlet L with
  3 ≤ q ≤ 400 000 with 4.9774 for 5.60; and Fodge's Dedekind companion (Zenodo 24354828, 18 Nov): 11.0416 log d_K + 8.6662 n_K log|t|
  − 18.7991 n_K + 0.8548 (floor 20.55, a table to 3.09) for Das–Gaba–Lee–Savalia–Wong's 12.21124, 9.54177, −11.59548 (arXiv 2506.19319).
* [falsify] Proposition 2.2 / Corollary 3.2 (exact Stechkin weight: Π_a(d, y²) ≥ κΠ_b(d, y²) for all 0 ≤ d < a and y with κ = a/b, b = a + δ,
  iff b ≥ √5 a; Π_x(d, Y) = (x − d)/((x − d)² + Y) + (x + d)/((x + d)² + Y)): min ratio 1 for b/a = √5, 3, 10 and 0.960 at b/a = 2.2,
  0.750 at 2 (grid of 801 d and 4000 log-spaced y at four a): HOLDS and sharp (orchestrator).
* [falsify] Propositions 4.4/4.6 (ten-shift weight with shifts (2σ − 1)m): pair kernel ≥ 2·10⁻¹⁶ on a 2000 × 4000 grid at four a,
  S = 0.5000052421819419, the y → ∞ coefficient 2(1 − Σκ_m(1 + 2m)) > 0; Prop. 4.6's S > 1/2 by t → 1: HOLDS (scout w10).
* [check] Theorem 1 input: Kadiri, Mathematika 64 (2018), arXiv math/0510570: 5.60 for exactly 3 ≤ q ≤ 400 000; no earlier improvement;
  Corollary 1.2 replaces BMOR's hypothesis Z(10⁸/q, 5.6) by Z(10⁸/q, 4.9774) (scout kadiri). The orchestrator reran the two final chain
  steps (rig_dirichlet_moving_v2.py small 4.9775 4.9774 0.001329 and large … 0.001376): B ≤ 4.97693632438… and 4.97735415316…, identical.
* [check] Theorem 2: every change to BTY's argument is proved in the paper (Props. 4.4, 8.7, Lemmas 8.3, 8.5, 8.8, 8.9, Appendix C);
  imported Ford, Hiary–Patel–Yang 0.618 (below log t = 27.4), Revers, Littlewood; cert_v25_46892.log: margin > 0 on all 806 subintervals,
  smallest 2.459·10⁻⁶ on [46.9186, 46.9300] (scout zthm2).
* [check] Dedekind: DGLSW's constants as quoted (from Fodge's TeX; their abstract omits them); the floor 3.09 is attained at
  (d_K, n_K, |t|) = (3, 2, 1), (23, 3, 1), (117, 4, 1) (σ ≥ 0.676375); bundled logs ALL VERIFICATIONS PASSED, not rerun (scout dedek).
* [check] where the new regions beat our Lean half-planes (θ₁₄₁ for every Dirichlet L, 31/32 for every ζ_K): ζ only for 3 ≤ t < 5.3858;
  Dirichlet only q ∈ {3, 4} with q·max(1, |t|) < 4.8855; Dedekind only where the bracket is < 32 (|t| = 1: d_K < 505.7 quadratic,
  < 2775 cubic) (scout regime). The records are classical-shape records with little reach beyond the half-planes.
* Exploration: Nevanlinna–Pick (R466; F(0) = 1 saturates the Schur class: no transfer), Fredholm determinants of integrable kernels
  (interlacing needs the missing self-adjoint realization), Turán–Nazarov lemmas (no coefficient cancellation, no density gain): no transfer.
-/

namespace Fayed.R484
end Fayed.R484
